!------------------------------------------------------------------------------
subroutine QEq(atype, pos, q)
use atoms; use parameters
! Two vector electronegativity equilization routine
!
! The linkedlist cell size is determined by the cutoff length of bonding 
! interaction <rc> = 3A. Since the non-bonding interaction cutoff <Rcut> = 10A,
! need to take enough layers to calculate non-bonding interactoins.
!
!<Gnew>, <Gold> :: NEW and OLD squre norm of Gradient vector.
!<Est> :: ElectroSTatic energy
!-------------------------------------------------------------------------------
implicit none

real(8),intent(in) :: atype(NBUFFER), pos(NBUFFER,3)
real(8),intent(out) :: q(NBUFFER)
real(8) :: vdummy(1,1), fdummy(1,1)

integer :: i,j,l2g
integer :: i1,j1,k1, nmax
real(8) :: Gnew(2), Gold(2) 
real(8) :: Est, GEst1, GEst2, g_h(2), h_hsh(2)
real(4) :: lmin(2)
real(8) :: buf(4), Gbuf(4)
real(8) :: ssum, tsum, mu
real(8) :: qsum, gqsum
real(8) :: QCopyDr(3)
real(8) :: hshs_sum,hsht_sum

call system_clock(i1,k1)

QCopyDr(1:3)=rctap/(/lata,latb,latc/)

!--- Initialize <s> vector with current charge and <t> vector with zero.
!--- isQEq==1 Normal QEq, isQEq==2 Extended Lagrangian method, DEFAULT skip QEq 
select case(isQEq)

!=== original QEq ===!
  case (1) 
!--- In the original QEq, fictitious charges are initialized with real charges
!--- and set zero.
    qsfp(1:NATOMS)=q(1:NATOMS)
    qsfv(1:NATOMS)=0.d0
!--- Initialization of the two vector QEq 
    qs(:)=0.d0
    qt(:)=0.d0
    qs(1:NATOMS)=q(1:NATOMS)
    nmax=NMAXQEq

!=== Extended Lagrangian method ===!
  case(2)
!--- charge mixing.
    qs(1:NATOMS)=Lex_fqs*qsfp(1:NATOMS)+(1.d0-Lex_fqs)*q(1:NATOMS)
!--- the same as the original QEq, set t vector zero
    qt(1:NATOMS)=0.d0
!--- just run one step
    nmax=1

!=== else, just return ===!
  case default
     return

end select

#ifdef QEQDUMP 
open(91,file="qeqdump"//trim(rankToString(myid))//".txt")
#endif

!--- copy atomic coords and types from neighbors, used in qeq_initialize()
call COPYATOMS(MODE_COPY, QCopyDr, atype, pos, vdummy, fdummy, q)
call LINKEDLIST(atype, pos, nblcsize, nbheader, nbllist, nbnacell, nbcc, MAXLAYERS_NB)

call qeq_initialize()

#ifdef QEQDUMP 
do i=1, NATOMS
   do j1=1,nbplist(0,i)
      j = nbplist(j1,i)
      write(91,'(4i6,4es25.15)') -1, l2g(atype(i)),nint(atype(i)),l2g(atype(j)),hessian(j1,i)
   enddo
enddo
#endif

!--- after the initialization, only the normalized coords are necessary for COPYATOMS()

GEst2=1.d99
nstep_qeq=0
do 
   call COPYATOMS(MODE_QCOPY1,QCopyDr, atype, pos, vdummy, fdummy, q)
   call get_gradient(Gnew)
   if (nstep_qeq==0) then
      hs(1:NATOMS) = gs(1:NATOMS)
      ht(1:NATOMS) = gt(1:NATOMS)
  end if
  if  (nstep_qeq > 0) then
      hs(1:NATOMS) = gs(1:NATOMS) + (Gnew(1)/Gold(1))*hs(1:NATOMS)
      ht(1:NATOMS) = gt(1:NATOMS) + (Gnew(2)/Gold(2))*ht(1:NATOMS)
  end if
  call COPYATOMS(MODE_QCOPY2,QCopyDr, atype, pos, vdummy, fdummy, q)

#ifdef QEQDUMP 
  qsum = sum(q(1:NATOMS))
  call MPI_ALLREDUCE(qsum, gqsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
#endif

  call get_hsh(Est,hshs_sum,hsht_sum)

  call MPI_ALLREDUCE(Est, GEst1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

#ifdef QEQDUMP 
  if(myid==0) print'(i5,5es25.15)', nstep_qeq, 0.5d0*log(Gnew(1:2)/GNATOMS), GEst1, GEst2, gqsum
#endif

  if( ( 0.5d0*( abs(GEst2) + abs(GEst1) ) < QEq_tol) ) exit 
  if( abs(GEst2) > 0.d0 .and. (abs(GEst1/GEst2-1.d0) < QEq_tol) ) exit
  if (nstep_qeq==nmax-1) exit

  GEst2 = GEst1

!--- line minimization factor of <s> vector
  g_h(1) = dot_product(gs(1:NATOMS), hs(1:NATOMS))
  h_hsh(1) = hshs_sum

!--- line minimization factor of <t> vector
  g_h(2) = dot_product(gt(1:NATOMS), ht(1:NATOMS))
  h_hsh(2) = hsht_sum

  buf(1)=g_h(1);   buf(2)=g_h(2)
  buf(3)=h_hsh(1); buf(4)=h_hsh(2)
  Gbuf(:)=0.d0
  call MPI_ALLREDUCE(buf, Gbuf, 4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  g_h(1) = Gbuf(1);   g_h(2) = Gbuf(2)
  h_hsh(1) = Gbuf(3); h_hsh(2) = Gbuf(4)

  lmin(1:2) = g_h(1:2)/h_hsh(1:2)

!--- line minimization for each vector
  qs(1:NATOMS) = qs(1:NATOMS) + lmin(1)*hs(1:NATOMS)
  qt(1:NATOMS) = qt(1:NATOMS) + lmin(2)*ht(1:NATOMS)

!--- get a current electronegativity <mu>
  ssum = sum(qs(1:NATOMS))
  tsum = sum(qt(1:NATOMS))
  buf(1) = ssum; buf(2) = tsum

  Gbuf(:)=0.d0
  call MPI_ALLREDUCE(buf, Gbuf, 2, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  ssum=Gbuf(1); tsum=Gbuf(2)

  mu = ssum/tsum

!--- update atom charges
  q(1:NATOMS) = qs(1:NATOMS) - mu*qt(1:NATOMS)

!--- save old residues.  
  Gold(:) = Gnew(:)
  nstep_qeq = nstep_qeq + 1
enddo


call system_clock(j1,k1)
it_timer(1)=it_timer(1)+(j1-i1)

! save # of QEq iteration 
it_timer(24)=it_timer(24)+nstep_qeq

#ifdef QEQDUMP 
close(91)
#endif

return 

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------
subroutine qeq_initialize()
use atoms; use parameters; use MemoryAllocator
! This subroutine create a neighbor list with cutoff length = 10[A] and save the hessian into <hessian>.  
! <nbrlist> and <hessian> will be used for different purpose later.
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,j, ity, jty, n, m, mn, nn
integer :: c1,c2,c3, c4,c5,c6
real(4) :: dr2
real(8) :: dr(3), drtb
integer :: itb, inxn

integer :: ti,tj,tk

call system_clock(ti,tk)


!$omp parallel do schedule(runtime), default(shared), &
!$omp private(i,j,ity,jty,n,m,mn,nn,c1,c2,c3,c4,c5,c6,dr,dr2,drtb,itb,inxn)
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)

   ity=nint(atype(i))
   nbplist(0,i) = 0

   do mn = 1, nbnmesh
      c4 = c1 + nbmesh(1,mn)
      c5 = c2 + nbmesh(2,mn)
      c6 = c3 + nbmesh(3,mn)

      j = nbheader(c4,c5,c6)
      do n=1, nbnacell(c4,c5,c6)

         if(i/=j) then
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            dr2 =  sum(dr(1:3)*dr(1:3))

            if(dr2 < rctap2) then

               jty = nint(atype(j))

!--- make a neighbor list with cutoff length = 10[A]
!$omp atomic
               nbplist(0,i) = nbplist(0,i) + 1
               nbplist(nbplist(0,i),i) = j

!--- get table index and residual value
               itb = int(dr2*UDRi)
               drtb = dr2 - itb*UDR
               drtb = drtb*UDRi

               inxn = inxn2(ity, jty)

               hessian(nbplist(0,i),i) = (1.d0-drtb)*TBL_Eclmb_QEq(itb,inxn) + drtb*TBL_Eclmb_QEq(itb+1,inxn)
            endif
         endif

         j=nbllist(j)
      enddo
   enddo !   do mn = 1, nbnmesh

   if (nbplist(0,i) > MAXNEIGHBS10) then
      write(6,*) "ERROR: In qeq.F90 inside qeq_initialize, nbplist greater then MAXNEIGHBS10 on mpirank",myid,"with value",nbplist(0,i)
      call MPI_FINALIZE(ierr)
   end if
   i=nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do

!--- for array size stat
if(mod(nstep,pstep)==0) then
  nn=maxval(nbplist(0,1:NATOMS))
  i=nstep/pstep+1
  maxas(i,3)=nn
endif

call system_clock(tj,tk)
it_timer(16)=it_timer(16)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_hsh(Est,hshs_sum,hsht_sum)
use atoms; use parameters
! This subroutine updates hessian*cg array <hsh> and the electrostatic energy <Est>.  
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Est
integer :: i,j,j1, ity
real(8) :: eta_ity, Est1
real(8) :: t_hshs,t_hsht
real(8) :: hshs_sum,hsht_sum

integer :: ti,tj,tk
call system_clock(ti,tk)

Est = 0.d0
hshs_sum = 0.d0
hsht_sum = 0.d0

!$omp parallel do default(shared), schedule(runtime), private(i,j,j1,ity,eta_ity,Est1,t_hshs,t_hsht),reduction(+:Est,hshs_sum,hsht_sum)
do i=1, NATOMS
   ity = nint(atype(i))
   eta_ity = eta(ity)

   t_hshs = eta_ity*hs(i)
   t_hsht = eta_ity*ht(i)

   Est = Est + chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i)

   do j1 = 1, nbplist(0,i)
      j = nbplist(j1,i)
      t_hshs = t_hshs + hessian(j1,i)*hs(j)
      t_hsht = t_hsht + hessian(j1,i)*ht(j)
!--- get half of potential energy, then sum it up if atoms are resident.
      Est1 = 0.5d0*hessian(j1,i)*q(i)*q(j)
      Est = Est + Est1
      if(j<=NATOMS) Est = Est + Est1
   enddo

   hshs_sum = hshs_sum + t_hshs*hs(i)
   hsht_sum = hsht_sum + t_hsht*ht(i)

enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(18)=it_timer(18)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_gradient(Gnew)
! Update gradient vector <g> and new residue <Gnew>
!-----------------------------------------------------------------------------------------------------------------------
use atoms; use parameters
implicit none
real(8),intent(OUT) :: Gnew(2)
real(8) :: eta_ity, ggnew(2)
integer :: i,j,j1, ity

real(8) :: gssum, gtsum

integer :: ti,tj,tk
call system_clock(ti,tk)

!$omp parallel do default(shared), schedule(runtime), private(gssum, gtsum, eta_ity,i,j,j1,ity)
do i=1,NATOMS

   gssum=0.d0
   gtsum=0.d0
   do j1=1, nbplist(0,i) 
      j = nbplist(j1,i)
      gssum = gssum + hessian(j1,i)*qs(j)
      gtsum = gtsum + hessian(j1,i)*qt(j)
   enddo

   ity = nint(atype(i))
   eta_ity = eta(ity)

   gs(i) = - chi(ity) - eta_ity*qs(i) - gssum
   gt(i) = - 1.d0     - eta_ity*qt(i) - gtsum

enddo 
!$omp end parallel do

ggnew(1) = dot_product(gs(1:NATOMS), gs(1:NATOMS))
ggnew(2) = dot_product(gt(1:NATOMS), gt(1:NATOMS))
call MPI_ALLREDUCE(ggnew, Gnew, size(ggnew), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

call system_clock(tj,tk)
it_timer(19)=it_timer(19)+(tj-ti)

end subroutine

end subroutine QEq

!------------------------------------------------------------------------------
program rxmd
use base; use init; use pqeq_vars; use parameters
use CG; use cmdline_args
!------------------------------------------------------------------------------
implicit none
integer :: i,ity,it1,it2,irt,provided
real(8) :: ctmp, dr(3)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

if(myid==0)  print'(a30)', 'rxmd has started'

!--- process command line arguments
call get_cmdline_args(myid, eFieldDir, eFieldStrength)

!--- read ffield file
CALL GETPARAMS(FFPath,FFDescript)

!--- initialize the MD system
CALL INITSYSTEM(atype, pos, v, f, q)

if(mdmode==10) call ConjugateGradient(atype,pos)

if(isPQEq) then
   call PQEq(atype, pos, q)
else
   call QEq(atype, pos, q)
endif
call FORCE(atype, pos, f, q)

!--- Enter Main MD loop 
call system_clock(it1,irt)

do nstep=0, ntime_step-1

   if(mod(nstep,pstep)==0) then
       call PRINTE(atype, v, q)
   endif
   if(mod(nstep,fstep)==0) &
        call OUTPUT(atype, pos, v, q, GetFileNameBase(DataDir,current_step+nstep))

   if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:NATOMS,1:3)=vsfact*v(1:NATOMS,1:3)

   if(mod(nstep,sstep)==0.and.mdmode==5) then
      ctmp = (treq*UTEMP0)/( GKE*UTEMP )
      v(1:NATOMS,1:3)=sqrt(ctmp)*v(1:NATOMS,1:3)
   endif

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call INITVELOCITY(atype, v)

!--- element-wise velocity scaling
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call ScaleTemperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==8) &
      call AdjustTemperature(atype, v)

!--- update velocity
   call vkick(1.d0, atype, v, f) 

!--- update coordinates
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))
   qsfp(1:NATOMS)=qsfp(1:NATOMS)+dt*qsfv(1:NATOMS)

!--- always correct the linear momentum when electric field is applied. 
   if(isEfield) call LinearMomentum(atype, v)
   pos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+dt*v(1:NATOMS,1:3)

!--- migrate atoms after positions are updated
   call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0],atype, pos, v, f, q)
   
   if(mod(nstep,qstep)==0) then
      if(isPQEq) then
         call PQEq(atype, pos, q)
      else
         call QEq(atype, pos, q)
      endif
   endif
   call FORCE(atype, pos, f, q)

   do i=1, NATOMS
      ity = nint(atype(i))
      astr(1)=astr(1)+v(i,1)*v(i,1)*mass(ity)
      astr(2)=astr(2)+v(i,2)*v(i,2)*mass(ity)
      astr(3)=astr(3)+v(i,3)*v(i,3)*mass(ity)
      astr(4)=astr(4)+v(i,2)*v(i,3)*mass(ity)
      astr(5)=astr(5)+v(i,3)*v(i,1)*mass(ity)
      astr(6)=astr(6)+v(i,1)*v(i,2)*mass(ity)
   end do

!--- update velocity
   call vkick(1.d0, atype, v, f) 
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))

enddo

!--- save the final configurations
call OUTPUT(atype, pos, v, q,  GetFileNameBase(DataDir,current_step+nstep))

!--- update rxff.bin in working directory for continuation run
call WriteBIN(atype, pos, v, q, GetFileNameBase(DataDir,-1))

call system_clock(it2,irt)
it_timer(Ntimer)=(it2-it1)

call FinalizeMD(irt)

call MPI_FINALIZE(ierr)
end PROGRAM

!------------------------------------------------------------------------------
subroutine FinalizeMD(irt)
use atoms; use MemoryAllocator
!------------------------------------------------------------------------------
implicit none
integer :: i
integer,intent(in) :: irt ! time resolution
integer,allocatable :: ibuf(:),ibuf1(:)

!--- close summary file
if(saveRunProfile) close(RunProfileFd)

allocate(ibuf(nmaxas),ibuf1(nmaxas))
ibuf(:)=0
do i=1,nmaxas
   ibuf(i)=maxval(maxas(:,i))
enddo
call MPI_ALLREDUCE(ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

call MPI_ALLREDUCE(it_timer, it_timer_max, Ntimer, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(it_timer, it_timer_min, Ntimer, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

if(myid==0) then
   print'(a)','----------------------------------------------'
   print'(a20,i12)', 'MAXNEIGHBS: ', ibuf1(2)
   print'(a20,i12)', 'MAXNEIGHBS10: ', ibuf1(3)
   print'(a20,i12)', 'MAXNBUFFER(MOVE): ', ibuf1(1)+ibuf1(4)
   print'(a20,i12)', 'MAXNBUFFER(COPY): ', ibuf1(1)+ibuf1(5)
   print'(a20,i12)', 'QEq Iterations: ', it_timer_max(24)
   print'(a20,f12.2)','Memory (MB): ', GetTotalMemory()*1d-6
   print*

   print'(a20,f12.4,3x,f12.4)','QEq: ',  dble(it_timer_max(1))/irt, dble(it_timer_min(1))/irt
   print'(a20,f12.4,3x,f12.4)','qeq_initialize: ',  dble(it_timer_max(16))/irt, dble(it_timer_min(16))/irt
   print'(a20,f12.4,3x,f12.4)','get_hsh: ',  dble(it_timer_max(18))/irt, dble(it_timer_min(18))/irt
   print'(a20,f12.4,3x,f12.4)','get_gradient: ',  dble(it_timer_max(19))/irt, dble(it_timer_min(19))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','LINKEDLIST: ',  dble(it_timer_max(3))/irt, dble(it_timer_min(3))/irt
   print'(a20,f12.4,3x,f12.4)','COPYATOMS: ',    dble(it_timer_max(4))/irt, dble(it_timer_min(4))/irt
   print'(a20,f12.4,3x,f12.4)','send_rec: ', dble(it_timer_max(25))/irt, dble(it_timer_min(25))/irt
   print'(a20,f12.4,3x,f12.4)','store_atoms: ', dble(it_timer_max(26))/irt, dble(it_timer_min(26))/irt
   print'(a20,f12.4,3x,f12.4)','append_atoms: ', dble(it_timer_max(27))/irt, dble(it_timer_min(27))/irt

   print'(a20,f12.4,3x,f12.4)','NEIGHBORLIST: ', dble(it_timer_max(5))/irt, dble(it_timer_min(5))/irt
   print'(a20,f12.4,3x,f12.4)','GetNBPairList: ', dble(it_timer_max(15))/irt, dble(it_timer_min(15))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','BOCALC: ', dble(it_timer_max(6))/irt, dble(it_timer_min(6))/irt
   print'(a20,f12.4,3x,f12.4)','ENbond: ', dble(it_timer_max(7))/irt, dble(it_timer_min(7))/irt
   print'(a20,f12.4,3x,f12.4)','Ebond: ', dble(it_timer_max(8))/irt, dble(it_timer_min(8))/irt
   print'(a20,f12.4,3x,f12.4)','Elnpr: ', dble(it_timer_max(9))/irt, dble(it_timer_min(9))/irt
   print'(a20,f12.4,3x,f12.4)','Ehb: ', dble(it_timer_max(10))/irt, dble(it_timer_min(10))/irt
   print'(a20,f12.4,3x,f12.4)','E3b: ', dble(it_timer_max(11))/irt, dble(it_timer_min(11))/irt
   print'(a20,f12.4,3x,f12.4)','E4b: ', dble(it_timer_max(12))/irt, dble(it_timer_min(12))/irt
   print'(a20,f12.4,3x,f12.4)','ForceBondedTerms: ', dble(it_timer_max(13))/irt, dble(it_timer_min(13))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','WriteBND: ', dble(it_timer_max(20))/irt, dble(it_timer_min(20))/irt
   print'(a20,f12.4,3x,f12.4)','WritePDB: ', dble(it_timer_max(21))/irt, dble(it_timer_min(21))/irt
   print'(a20,f12.4,3x,f12.4)','ReadBIN: ', dble(it_timer_max(22))/irt, dble(it_timer_min(22))/irt
   print'(a20,f12.4,3x,f12.4)','WriteBIN: ', dble(it_timer_max(23))/irt, dble(it_timer_min(23))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','total (sec): ',dble(it_timer_max(Ntimer))/irt, dble(it_timer_min(Ntimer))/irt

   print'(a)','----------------------------------------------'

   print'(a30)', 'rxmd successfully finished'
endif

deallocate(ibuf,ibuf1)

end subroutine

!------------------------------------------------------------------------------
subroutine vkick(dtf, atype, v, f)
use atoms
!------------------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER),v(NBUFFER,3),f(NBUFFER,3)

integer :: i, ity
real(8) :: dtf

do i=1,NATOMS
   ity = nint(atype(i))
   v(i,1:3) = v(i,1:3) + dtf*dthm(ity)*f(i,1:3)
enddo

end subroutine

!----------------------------------------------------------------------------------------
subroutine PRINTE(atype, v, q)
use atoms; use parameters; use MemoryAllocator
! calculate the kinetic energy and sum up all of potential energies, then print them.
!----------------------------------------------------------------------------------------
implicit none

real(8),intent(in) :: atype(NBUFFER), q(NBUFFER)
real(8),intent(in) :: v(NBUFFER,3)

integer :: i,ity,cstep
real(8) :: qq=0.d0,tt=0.d0,ss=0.d0,buf(0:23)

i=nstep/pstep+1
maxas(i,1)=NATOMS

KE=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   KE = KE + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo
qq=sum(q(1:NATOMS))

!--- pressure 
ss=sum(astr(1:3))/3.d0

!--- potential energy 
PE(0)=sum(PE(1:13))

!--- copy data into buffer
buf(0:13) = PE(0:13)
buf(14) = KE; buf(15) = ss; buf(16) = qq
buf(17:22)=astr(1:6)
call MPI_ALLREDUCE (MPI_IN_PLACE, buf, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!--- copy data from buffer
GPE(0:13) = buf(0:13)
GKE = buf(14); ss = buf(15); qq = buf(16)
astr(1:6)=buf(17:22)

!--- compute properties
GPE(:)=GPE(:)/GNATOMS
GKE=GKE/GNATOMS
tt=GKE*UTEMP
astr(1:6)=astr(1:6)/MDBOX*USTRS/pstep
ss=ss/MDBOX*USTRS/pstep

!--- total energy
GTE = GKE + GPE(0)
if(myid==0) then
   
   cstep = nstep + current_step 

   write(6,'(a,i9,3es13.5,6es11.3,1x,3f8.2,i4,f8.2,f8.2)') 'MDstep: ',cstep,GTE,GPE(0),GKE, &
   GPE(1),sum(GPE(2:4)),sum(GPE(5:7)),sum(GPE(8:9)),GPE(10),sum(GPE(11:13)), &
   tt, ss, qq, nstep_qeq, GetTotalMemory()*1e-9, MPI_WTIME()-wt0

   !write(6,'(a,i9,6f12.6)') 'stress : ',cstep,astr(1:6)

endif

!--- reset stress tensor accumulator
astr(1:6)=0.d0

!--- save current time
wt0 = MPI_WTIME()
end subroutine

!----------------------------------------------------------------------------------------
subroutine LINKEDLIST(atype, rreal, cellDims, headAtom, atomList, NatomPerCell, Ncells, NLAYERS)
use atoms
! partitions the volume into linked-list cells <lcsize>
!----------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: atype(NBUFFER), rreal(3,NBUFFER), cellDims(3)

integer,intent(in) :: Ncells(3), NLAYERS
integer,intent(out) :: atomList(NBUFFER)
integer,intent(out) :: NatomPerCell(-NLAYERS:Ncells(1)-1+NLAYERS, &
                                    -NLAYERS:Ncells(2)-1+NLAYERS, &
                                    -NLAYERS:Ncells(3)-1+NLAYERS) 
integer,intent(out) :: headAtom(-NLAYERS:Ncells(1)-1+NLAYERS, & 
                                -NLAYERS:Ncells(2)-1+NLAYERS, &
                                -NLAYERS:Ncells(3)-1+NLAYERS) 

real(8) :: rnorm(NBUFFER,3)
integer :: n, l(3), j

integer :: ti,tj,tk
call system_clock(ti,tk)

call xu2xs(copyptr(6),rreal,rnorm)

headAtom(:,:,:) = -1; atomList(:) = 0; NatomPerCell(:,:,:)=0

!--- copyptr(6) stores the last atom index copied in COPYATOMS.
do n=1, copyptr(6) 

   if(nint(atype(n))==0) cycle

   l(1:3) = floor(rnorm(n,1:3)/cellDims(1:3))

   atomList(n) = headAtom(l(1), l(2), l(3))
   headAtom(l(1), l(2), l(3)) = n
   NatomPerCell(l(1), l(2), l(3)) = NatomPerCell(l(1), l(2), l(3)) + 1
enddo

call system_clock(tj,tk)
it_timer(3)=it_timer(3)+(tj-ti)

end subroutine 

!----------------------------------------------------------------------
subroutine NEIGHBORLIST(nlayer, atype, pos)
use atoms; use parameters
! calculate neighbor list for atoms witin cc(1:3, -nlayer:nlayer) cells.
!----------------------------------------------------------------------
implicit none
integer,intent(in) :: nlayer
real(8),intent(in) :: atype(NBUFFER), pos(NBUFFER,3)

integer :: c1,c2,c3, ic(3), c4, c5, c6
integer :: n, n1, m, m1, nty, mty, inxn
real(8) :: dr(3), dr2

integer :: i,j,i1,j1
logical :: isFound

integer :: ti,tj,tk
call system_clock(ti,tk)

nbrlist(:,0) = 0

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,inxn,dr,dr2) 
DO c1=-nlayer, cc(1)-1+nlayer
DO c2=-nlayer, cc(2)-1+nlayer
DO c3=-nlayer, cc(3)-1+nlayer

  m = header(c1, c2, c3)
  do m1=1, nacell(c1, c2, c3)
     mty = nint(atype(m))

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1, c2, c3] + [c4, c5, c6]

        n = header(ic(1),ic(2),ic(3))
        do n1=1, nacell(ic(1), ic(2), ic(3))

           if(n/=m) then
             nty = nint(atype(n))
             inxn = inxn2(mty, nty)

             dr(1:3) = pos(n,1:3) - pos(m,1:3) 
             dr2 = sum(dr(1:3)*dr(1:3))

             if(dr2<rc2(inxn)) then 
                nbrlist(m, 0) = nbrlist(m, 0) + 1
                nbrlist(m, nbrlist(m, 0)) = n
             endif 
           endif

           n=llist(n) 
        enddo
     enddo; enddo; enddo

     m = llist(m)
  enddo
enddo; enddo; enddo
!$omp end parallel do 

!--- to get the reverse information (i.e. from i,j1&j to i1), store <i1> into <nbrindx>.

!$omp parallel do default(shared) private(i,i1,j,j1,isFound)
do i=1, copyptr(6)
   do i1 = 1, nbrlist(i,0)
      j = nbrlist(i,i1)
      isFound=.false.
      do j1 = 1, nbrlist(j,0)
         if(i == nbrlist(j,j1)) then
            nbrindx(i,i1)=j1
            isFound=.true.
         endif
      enddo
      if(.not.isFound) &
      print'(a,i6,30i4)','ERROR: inconsistency between nbrlist and nbrindx found', &
           myid, i,nbrlist(i,0:nbrlist(i,0)), j, nbrlist(j,0:nbrlist(j,0))
   enddo
enddo
!$omp end parallel do

!--- error trap
n=maxval(nbrlist(1:NATOMS,0))
if(n > MAXNEIGHBS) then
   write(6,'(a45,2i5)') "ERROR: overflow of max # in neighbor list, ", myid, n
   call MPI_FINALIZE(ierr)
   stop
endif

!--- for array size stat
if(mod(nstep,pstep)==0) then
  maxas(nstep/pstep+1,2)=maxval(nbrlist(1:NATOMS,0))
endif

call system_clock(tj,tk)
it_timer(5)=it_timer(5)+(tj-ti)

end subroutine

!----------------------------------------------------------------------
subroutine GetNonbondingPairList(pos)
use atoms; use parameters
!----------------------------------------------------------------------
implicit none

real(8),intent(in) :: pos(NBUFFER,3)

integer :: c1,c2,c3,c4,c5,c6,i,j,m,n,mn
integer :: l2g
real(8) :: dr(3), dr2

integer :: ti,tj,tk

integer :: m_size=0                ! keep track of the size of the packed_indices and packed_coordinates
integer, parameter :: max_pack = 256  ! maximum packing size for packing neighborlist
integer :: packed_indices(1:max_pack)   ! contains the indicies of the neighbor for each atom
real(8) :: packed_coordinates(1:max_pack,3)  ! contains the atomic coordinates og the packed neighbor

call system_clock(ti,tk)

! reset non-bonding pair list
!nbplist(0,:) = 0

!$omp parallel do default(shared),private(c1,c2,c3,c4,c5,c6,i,j,m,n,mn,m_size,packed_indices,packed_coordinates)
do c1 = 0, nbcc(1)-1
do c2 = 0, nbcc(2)-1
do c3 = 0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)
      
      m_size = 0
      nbplist(0,i) = 0
      do mn = 1, nbnmesh
         c4 = c1 + nbmesh(1,mn)
         c5 = c2 + nbmesh(2,mn)
         c6 = c3 + nbmesh(3,mn)

         j = nbheader(c4,c5,c6)
         do n = 1, nbnacell(c4,c5,c6)
            if (i /= j) then
               m_size = m_size + 1
               packed_indices(m_size) = j
               packed_coordinates(m_size,:) = pos(j, 1:3)
            end if
            if (m_size == max_pack) then
               call calc_packed_neighbor(m_size, max_pack, pos(i,1:3), packed_coordinates, packed_indices, nbplist(:,i))
               m_size = 0
            end if
            j = nbllist(j)
         enddo
       enddo
       if (m_size > 0) then
          call calc_packed_neighbor(m_size, max_pack, pos(i,1:3), packed_coordinates, packed_indices, nbplist(:,i))
          m_size = 0
       end if
      i = nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(15)=it_timer(15)+(tj-ti)

end subroutine

!----------------------------------------------------------------------
subroutine  calc_packed_neighbor(m_size, max_pack, posi, packed_coordinates, packed_indices, nbplist_i)
use atoms; use parameters
!----------------------------------------------------------------------

implicit None
real(8), intent(in) :: posi(3)
integer, intent(in) :: m_size, max_pack
real(8), intent(in) :: packed_coordinates(1:max_pack,3)  ! contains the atomic coordinates of the packed neighbor
integer, intent(in) :: packed_indices(1:max_pack)
integer, intent(inout) :: nbplist_i(0:MAXNEIGHBS10)
real(8) :: dr2(1:max_pack), dr(3)  ! contanins distance square for the entire batch of packed atoms to the reference atom
integer :: i_pack, i_counter

! compute distance square
do i_pack = 1, m_size
   dr(1:3) = posi(1:3) - packed_coordinates(i_pack,1:3)
   dr2(i_pack) = sum(dr(1:3)*dr(1:3))
end do

! construct neighbour list
i_counter = nbplist_i(0)
do i_pack = 1, m_size
   if (dr2(i_pack) <= rctap2) then
       i_counter = i_counter + 1
       nbplist_i(i_counter) = packed_indices(i_pack)
   end if
end do
nbplist_i(0) = i_counter

end subroutine


!----------------------------------------------------------------------
subroutine angular_momentum(atype, pos, v)
use atoms; use parameters
!----------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER), pos(NBUFFER,3),v(NBUFFER,3)

integer :: i,ity
real(8) :: com(3), Gcom(3), intsr(3,3), Gintsr(3,3), intsr_i(3,3), angm(3), Gangm(3), angv(3), mm, Gmm
real(8) :: dr(3), dv(3)

!--- get center of mass
com(:)=0.d0;     Gcom(:)=0.d0
mm=0.d0; Gmm=0.d0

do i=1, NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
   com(1:3) = mass(ity)*pos(i,1:3)
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, mm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gmm = mm
call MPI_ALLREDUCE(MPI_IN_PLACE, com, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gcom = com
Gcom(1:3) = Gcom(1:3)/Gmm

!--- get the angular momentum and inertia tensor from the com
angm(:)=0.d0;    Gangm(:)=0.d0
intsr(:,:)=0.d0; Gintsr(:,:)=0.d0

do i=1, NATOMS
   dr(1:3) = pos(i,1:3) - Gcom(1:3)
   
   angm(1) = mass(ity)*( dr(2)*v(i,3)-dr(3)*v(i,2) )
   angm(2) = mass(ity)*( dr(3)*v(i,1)-dr(1)*v(i,3) )
   angm(3) = mass(ity)*( dr(1)*v(i,2)-dr(2)*v(i,1) )

   intsr(1,1) = mass(ity)*( dr(2)**2+dr(3)**2 )
   intsr(2,2) = mass(ity)*( dr(3)**2+dr(1)**2 )
   intsr(3,3) = mass(ity)*( dr(1)**2+dr(2)**2 )

   intsr(1,2) =-mass(ity)*( dr(1)*dr(2) )
   intsr(1,3) =-mass(ity)*( dr(1)*dr(3) )
   intsr(2,3) =-mass(ity)*( dr(2)*dr(3) )

   intsr(2,1) = intsr(1,2)
   intsr(3,1) = intsr(1,3)
   intsr(3,2) = intsr(2,3)
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, angm, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gangm = angm
call MPI_ALLREDUCE(MPI_IN_PLACE, intsr, 9, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gintsr = intsr

!--- get angular velocity
call matinv(Gintsr, intsr_i)

angv(1) = sum(intsr_i(1,1:3)*angm(1:3))
angv(2) = sum(intsr_i(2,1:3)*angm(1:3))
angv(3) = sum(intsr_i(3,1:3)*angm(1:3))


!--- correct rotational motion wrt CoM.
do i=1,NATOMS
   dr(1:3) = pos(i,1:3) - Gcom(1:3)
   dv(1) = angv(2)*dr(3) - angv(3)*dr(2)
   dv(2) = angv(3)*dr(1) - angv(1)*dr(3)
   dv(3) = angv(1)*dr(2) - angv(2)*dr(1)

   v(i,1:3) = v(i,1:3) - dv(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine matinv(m1,m2)
! get inverse of m1 and save to m2
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8) :: m1(3,3), m2(3,3), detm

m2(1,1) = m1(2,2)*m1(3,3)-m1(2,3)*m1(3,2)
m2(1,2) = m1(1,3)*m1(3,2)-m1(1,2)*m1(3,3)
m2(1,3) = m1(1,2)*m1(2,3)-m1(1,3)*m1(2,2)
m2(2,1) = m1(2,3)*m1(3,1)-m1(2,1)*m1(3,3)
m2(2,2) = m1(1,1)*m1(3,3)-m1(1,3)*m1(3,1)
m2(2,3) = m1(1,3)*m1(2,1)-m1(1,1)*m1(2,3)
m2(3,1) = m1(2,1)*m1(3,2)-m1(2,2)*m1(3,1)
m2(3,2) = m1(1,2)*m1(3,1)-m1(1,1)*m1(3,2)
m2(3,3) = m1(1,1)*m1(2,2)-m1(1,2)*m1(2,1)

detm = m1(1,1)*m1(2,2)*m1(3,3) + m1(1,2)*m1(2,3)*m1(3,1) &
     + m1(1,3)*m1(2,1)*m1(3,2) - m1(1,3)*m1(2,2)*m1(3,1) &
     - m1(1,2)*m1(2,1)*m1(3,3) - m1(1,1)*m1(2,3)*m1(3,2) 

m2(:,:) = m2(:,:)/detm

end subroutine

!--------------------------------------------------------------------------------------------------------------
function l2g(atype)
implicit none
!convert Local ID to Global ID 
!--------------------------------------------------------------------------------------------------------------
real(8),intent(IN) :: atype
integer :: l2g,ity

ity = nint(atype)
l2g = nint((atype-ity)*1d13)

return
end function

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs(nmax, rreal, rnorm)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: rreal(NBUFFER,3)
real(8),intent(out) :: rnorm(NBUFFER,3)
integer,intent(in) :: nmax

real(8) :: rr(3)
integer :: i

do i=1,nmax
   rr(1:3) = rreal(i,1:3)
   rnorm(i,1)=sum(HHi(1,1:3)*rr(1:3))
   rnorm(i,2)=sum(HHi(2,1:3)*rr(1:3))
   rnorm(i,3)=sum(HHi(3,1:3)*rr(1:3))
   rnorm(i,1:3) = rnorm(i,1:3) - OBOX(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs_inplace(nmax, rreal)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(inout) :: rreal(NBUFFER,3)
integer,intent(in) :: nmax
real(8) :: rr(3)
integer :: i

do i=1,nmax
   rr(1:3) = rreal(i,1:3)
   rreal(i,1)=sum(HHi(1,1:3)*rr(1:3))
   rreal(i,2)=sum(HHi(2,1:3)*rr(1:3))
   rreal(i,3)=sum(HHi(3,1:3)*rr(1:3))
   rreal(i,1:3) = rreal(i,1:3) - OBOX(1:3)
enddo

end subroutine


!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu(nmax,rnorm,rreal)
! update real coordinate from normalized coordinate
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: rnorm(NBUFFER,3)
real(8),intent(out) :: rreal(NBUFFER,3)
integer,intent(in) :: nmax

real(8) :: rr(3)
integer :: i

do i=1,nmax 
   rr(1:3) = rnorm(i,1:3) + OBOX(1:3)
   rreal(i,1)=sum(HH(1,1:3,0)*rr(1:3))
   rreal(i,2)=sum(HH(2,1:3,0)*rr(1:3))
   rreal(i,3)=sum(HH(3,1:3,0)*rr(1:3))
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu_inplace(nmax,rnorm)
! update real coordinate from normalized coordinate
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(inout) :: rnorm(NBUFFER,3)
integer,intent(in) :: nmax

real(8) :: rr(3)
integer :: i

do i=1,nmax 
   rr(1:3) = rnorm(i,1:3) + OBOX(1:3)
   rnorm(i,1)=sum(HH(1,1:3,0)*rr(1:3))
   rnorm(i,2)=sum(HH(2,1:3,0)*rr(1:3))
   rnorm(i,3)=sum(HH(3,1:3,0)*rr(1:3))
enddo

end subroutine

!-----------------------------------------------------------------------
subroutine AdjustTemperature(atype, v)
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER), v(NBUFFER,3)

integer :: i,ity
real(8) :: Ekinetic, ctmp

Ekinetic=0.d0

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic=Ekinetic+0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, Ekinetic, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
Ekinetic=Ekinetic/GNATOMS

! adjust if current temperature deviates from treq by 5%
ctmp = sqrt( (treq*UTEMP0)/(Ekinetic*UTEMP) )
if( abs(ctmp-1.d0) > 0.05d0) then

   do i=1, NATOMS
      ity=nint(atype(i))
      v(i,1:3)=ctmp*v(i,1:3)
   enddo

   call LinearMomentum(atype, v)

endif


return
end


!-----------------------------------------------------------------------
subroutine ScaleTemperature(atype, v)
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER), v(NBUFFER,3)

integer :: i,ity
integer,parameter :: MAX_ELEMENT=20
real(8) :: Ekinetic(2,MAX_ELEMENT), ctmp(MAX_ELEMENT)

Ekinetic(:,:)=0.d0

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic(1,ity)=Ekinetic(1,ity)+1.d0 
   Ekinetic(2,ity)=Ekinetic(2,ity)+0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, Ekinetic, size(Ekinetic), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

do ity=1, MAX_ELEMENT
   if(Ekinetic(1,ity)>1.d0) then
      ctmp(ity) = Ekinetic(2,ity)/Ekinetic(1,ity)
      ctmp(ity) = sqrt( (treq*UTEMP0)/(ctmp(ity)*UTEMP) )
      if(myid==0) print'(a,i4,f8.3,i9,es15.5)', &
          'atom_type, scaling_factor, num_samples, previous_Ekin: ', &
          ity, ctmp(ity), nint(Ekinetic(1,ity)), Ekinetic(2,ity)
   else
      ctmp(ity) = 0.d0
   endif
enddo


do i=1, NATOMS
   ity=nint(atype(i))
   v(i,1:3)=ctmp(ity)*v(i,1:3)
enddo

call LinearMomentum(atype, v)

return
end

!-----------------------------------------------------------------------
subroutine LinearMomentum(atype, v)
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER), v(NBUFFER,3)

integer :: i,ity
real(8) :: mm,vCM(3),sbuf(4)

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   vCM(1:3)=vCM(1:3) + mass(ity)*v(i,1:3)
   mm = mm + mass(ity)
enddo

sbuf(1)=mm; sbuf(2:4)=vCM(1:3)
call MPI_ALLREDUCE (MPI_IN_PLACE, sbuf, size(sbuf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
mm=sbuf(1); vCM(1:3)=sbuf(2:4)

!--- get the global momentum
vCM(:)=vCM(:)/mm

!--- set the total momentum to be zero 
do i=1, NATOMS
   v(i,1:3) = v(i,1:3) - vCM(1:3)
enddo

return
end

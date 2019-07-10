!------------------------------------------------------------------------------------------
module init
use fileio
use ensembles
!------------------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------------------
SUBROUTINE INITSYSTEM(atype, pos, v, f, q)
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
use parameters
use pqeq_vars
use MemoryAllocator 

implicit none

real(8),allocatable,dimension(:) :: atype, q
real(8),allocatable,dimension(:,:) :: pos,v,f

integer :: i,j,k, ity, l(3), ist=0
real(8) :: mm, gmm, dns, mat(3,3)
integer(8) :: i8
real(8) :: rcsize(3), maxrcell

wt0 = MPI_WTIME()


!--- for PQEq. cutoff length
if(isPQEq) then
  rctap = rctap0_pqeq
else
  rctap = rctap0
endif

rctap2 = rctap**2

CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(rctap)**4, &
          84.d0/(rctap)**5, -70.d0/(rctap)**6, &
          20.d0/(rctap)**7 /)

if(isPQEq) then
   call initialize_pqeq(chi,eta)
   if(isEfield) call initialize_eField(myid)
endif

astr(:)=0.d0

!--- an error trap
if(vprocs(1)*vprocs(2)*vprocs(3) /= nprocs ) then
  if(myid==0) write(6,'(a60,3i3,i5)')  &
  "ERROR: requested/allocated # of procs are not consistent: ", vprocs(1:3), nprocs
  call MPI_FINALIZE(ierr)
  stop
endif

!--- initialize charge with QEq
if(mdmode==0) then
  if(myid==0) then
    print'(a,f12.3,a,i6,a)', &
         'INFO: mdmode==0, setting isQEQ is 1. Atomic velocities are scaled to ', &
          treq, ' [K] every ', sstep, ' steps.'
  endif
  isQEq=1
endif

!--- time unit conversion from [fs] -> time unit
dt = dt/UTIME

!--- square the spring const in the extended Lagrangian method 
Lex_w2=2.d0*Lex_k/dt/dt

!--- get reduced temperature from [K]
treq=treq/UTEMP0

!--- setup the vector ID and parity for processes, in x, y and z order.
vID(1)=mod(myid,vprocs(1))
vID(2)=mod(myid/vprocs(1),vprocs(2))
vID(3)=myid/(vprocs(1)*vprocs(2))
myparity(1)=mod(vID(1),2)
myparity(2)=mod(vID(2),2)
myparity(3)=mod(vID(3),2)

k=0
do i=1,3
!--- Add (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1) to vector ID applying PBC. 
   do j=1,-1,-2

!--- save vector ID.
      l(1:3) = vID(1:3)

!--- apply PBC.
      l(i)=mod( (vID(i) + j + vprocs(i)), vprocs(i) )
     
!--- get direction index.
      k=k+1 ![123456] 

!--- convert vector ID to sequential ID.
      target_node(k) = l(1)+l(2)*vprocs(1)+l(3)*vprocs(1)*vprocs(2) 
   enddo         
                   
enddo    

!--- dt/2*mass, mass/2
call allocatord1d(dthm, 1, nso)
call allocatord1d(hmas, 1, nso)
do ity=1, nso
   dthm(ity) = dt*0.5d0/mass(ity)
   hmas(ity) = 0.5d0*mass(ity)
enddo

call allocatord1d(atype,1,NBUFFER)
call allocatord1d(q,1,NBUFFER)
call allocatord2d(pos,1,NBUFFER,1,3)
call allocatord2d(v,1,NBUFFER,1,3)
call allocatord2d(f,1,NBUFFER,1,3)

!--- For PQEq
if(isPQEq) then
   call allocatord2d(spos,1,NBUFFER,1,3)
   spos(:,:)=0.d0
endif

call allocatord1d(deltalp,1,NBUFFER)

if(isRunFromXYZ) then
  call ReadXYZ(atype, pos, v, q, f, RunFromXYZPath)
else
  call ReadBIN(atype, pos, v, q, f, trim(DataDir)//"/rxff.bin")
endif

!--- Varaiable for extended Lagrangian method
call allocatord1d(qtfp,1,NBUFFER)
call allocatord1d(qtfv,1,NBUFFER)
qtfp(:)=0.d0; qtfv(:)=0.d0

!call GetBoxParams(mat,lata,latb,latc,lalpha,lbeta,lgamma)
!do i=1, 3
!do j=1, 3
!   HH(i,j,0)=mat(i,j)
!enddo; enddo

!--- get total number of atoms per type. This will be used to determine
!--- subroutine cutofflength() 
allocate(natoms_per_type(nso)) ! NOTE 8byte int is not supported in MemoryAllocator
natoms_per_type(:)=0
do i=1, NATOMS
   ity=nint(atype(i))
   natoms_per_type(ity)=natoms_per_type(ity)+1
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, natoms_per_type, nso, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- determine cutoff distances only for exsiting atom pairs
call CUTOFFLENGTH()

!--- update box-related variables
call UpdateBoxParams()

!--- get global number of atoms
GNATOMS=NATOMS ! Convert 4 byte to 8 byte
call MPI_ALLREDUCE(MPI_IN_PLACE, GNATOMS, 1, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- Linked List & Near Neighb Parameters
call allocatori2d(nbrlist,1,NBUFFER,0,MAXNEIGHBS)
call allocatori2d(nbrindx,1,NBUFFER,1,MAXNEIGHBS)
call allocatori2d(nbplist,0,MAXNEIGHBS10,1,NBUFFER)
call allocatori1d(llist,1,NBUFFER)
call allocatori3d(header,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)
call allocatori3d(nacell,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)

!--- Bond Order Prime and deriv terms:
call allocatord3d(dln_BOp,1,3,1,NBUFFER,1,MAXNEIGHBS)
call allocatord2d(dBOp,1,NBUFFER,1,MAXNEIGHBS)
call allocatord2d(deltap,1,NBUFFER,1,3)

!--- Bond Order terms
call allocatord3d(BO,0,3,1,NBUFFER,1,MAXNEIGHBS)
call allocatord1d(delta,1,NBUFFER)
call allocatord2d(A0,1,NBUFFER,1,MAXNEIGHBS)
call allocatord2d(A1,1,NBUFFER,1,MAXNEIGHBS)
call allocatord2d(A2,1,NBUFFER,1,MAXNEIGHBS)
call allocatord2d(A3,1,NBUFFER,1,MAXNEIGHBS)
call allocatord1d(nlp,1,NBUFFER)
call allocatord1d(dDlp,1,NBUFFER)
call allocatord1d(ccbnd,1,NBUFFER)
call allocatord1d(cdbnd,1,NBUFFER)
ccbnd(:)=0.d0; cdbnd(:)=0.d0

!--- 2 vector QEq varialbes
call allocatord1d(qs,1,NBUFFER)
call allocatord1d(gs,1,NBUFFER)
call allocatord1d(qt,1,NBUFFER)
call allocatord1d(gt,1,NBUFFER)
call allocatord1d(hs,1,NBUFFER)
call allocatord1d(hshs,1,NBUFFER)
call allocatord1d(ht,1,NBUFFER)
call allocatord1d(hsht,1,NBUFFER)
call allocatord2d(hessian,1,MAXNEIGHBS10,1,NBUFFER)
qs(:)=0.d0; qt(:)=0.d0; gs(:)=0.d0; gt(:)=0.d0; hs(:)=0.d0; ht(:)=0.d0; hshs(:)=0.d0; hsht(:)=0.d0

!--- returning force index array 
call allocatord1d(frcindx,1,NBUFFER)

!--- setup potential table
call POTENTIALTABLE()

!--- get real size of linked list cell
rcsize(1) = lata/vprocs(1)/cc(1)
rcsize(2) = latb/vprocs(2)/cc(2)
rcsize(3) = latc/vprocs(3)/cc(3)
maxrcell = maxval(rcsize(1:3))

!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
call GetNonbondingMesh()

!--- get density 
mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, mm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
gmm = mm
dns=gmm/MDBOX*UDENS

!--- allocate & initialize Array size Stat variables
i=ntime_step/pstep+1
allocate(maxas(i,nmaxas))
maxas(:,:)=0

!--- for spring force
call allocatord2d(ipos,1,NBUFFER,1,3)
ipos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)

!--- Setup ensemble. treq has to be set beforehand.
call initialize_ensemble()

!--- print out parameters and open data file
if(myid==0) then
   write(6,'(a)') "----------------------------------------------------------------"
   write(6,'(a30,i9,a3,i9)') "req/alloc # of procs:", vprocs(1)*vprocs(2)*vprocs(3), "  /",nprocs
   write(6,'(a30,3i9)')      "req proc arrengement:", vprocs(1),vprocs(2),vprocs(3)
   write(6,'(a30,a70)')      "parameter set:", FFDescript
   write(6,'(a30,es12.2)')   "time step[fs]:",dt*UTIME
   write(6,'(a30,i3, i10, i10)') "MDMODE CURRENTSTEP NTIMESTPE:", &
                                  mdmode, current_step, ntime_step
   write(6,'(a30,i6,es10.1,i6,i6)') "isQEq,QEq_tol,NMAXQEq,qstep:", &
                                     isQEq,QEq_tol,NMAXQEq,qstep
   write(6,'(a30,f8.3,f8.3)') 'Lex_fqs,Lex_k:',Lex_fqs,Lex_k
   write(6,'(a30,f12.3,f8.3,i9)') 'treq,vsfact,sstep:',treq*UTEMP0, vsfact, sstep
   write(6,'(a30,2i6)') 'fstep,pstep:', fstep,pstep
   write(6,'(a30,i24,i24)') "NATOMS GNATOMS:", NATOMS, GNATOMS
   write(6,'(a30,3f12.3)') "LBOX:",LBOX(1:3)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,1,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,2,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,3,0)
   print'(a30,3f12.3)', 'lata,latb,latc:', lata,latb,latc
   print'(a30,3f12.3)', 'lalpha,lbeta,lgamma:', lalpha,lbeta,lgamma
   write(6,'(a30,3f10.4)') "density [g/cc]:",dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", &
   maxrc,lata/cc(1)/vprocs(1),latb/cc(2)/vprocs(2),latc/cc(3)/vprocs(3)
   write(6,'(a30,3i6)')  '# of linkedlist cell (NB):', nbcc(1:3)
   write(6,'(a30,3f10.2)') "lcsize [A] (NB):", &
   lata/nbcc(1)/vprocs(1),latb/nbcc(2)/vprocs(2),latc/nbcc(3)/vprocs(3)
   write(6,'(a30,2i6)') "MAXNEIGHBS, MAXNEIGHBS10:", MAXNEIGHBS,MAXNEIGHBS10
   write(6,'(a30,i6,i9)') "NMINCELL, NBUFFER:", NMINCELL, NBUFFER
   write(6,'(a30,3(a12,1x))') "FFPath, DataDir, ParmPath:", &
                          trim(FFPath), trim(DataDir), trim(ParmPath)

   print'(a30 $)','# of atoms per type:'
   do ity=1, nso
      if(natoms_per_type(ity)>0) print'(i12,a2,i2 $)',natoms_per_type(ity),' -',ity
   enddo
   print*

   if(isSpring) then
     print'(a,f8.2 $)','springConst [kcal/mol] ', springConst 
     do ity=1, size(hasSpringForce) 
        if(hasSpringForce(ity)) print'(i3 $)',ity
     enddo
     print*
   endif

   print'(a)', "----------------------------------------------------------------"
   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

END SUBROUTINE
end module

!------------------------------------------------------------------------------------------
SUBROUTINE INITVELOCITY(atype, v)
use parameters; use atoms
! Generate gaussian distributed velocity as an initial value  using Box-Muller algorithm
!------------------------------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER)
real(8) :: v(NBUFFER,3)

integer :: i, k, ity
real(8) :: vv(2), vsqr, vsl, rndm(2)
real(8) :: vCM(3), GvCM(3), mm, Gmm
real(8) :: vfactor

!--- assign velocity to two atoms together with BM algoritm. 
!--- If <NATOMS> is odd, the <NATOMS> + 1 element will be the ignored in later calculations.

do i=1, NATOMS, 2

  do k=1,3 ! three directions
     !--- generate gaussian distributed velocity
     vsqr=0.d0
     do while ( (vsqr >= 1.d0) .or. (vsqr==0.d0) ) 
        call random_number(rndm)
        vv(1) = 2.d0 * rndm(1) - 1.d0
        vv(2) = 2.d0 * rndm(2) - 1.d0
        vsqr = vv(1)**2 + vv(2)**2
     enddo

     vsl = sqrt(-2.d0 * log(vsqr)/vsqr)
     v(i,k)   = vv(1)*vsl
     v(i+1,k) = vv(2)*vsl
  enddo
  
enddo

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   vCM(1:3)=vCM(1:3) + mass(ity)*v(i,1:3)
   mm = mm + mass(ity)
enddo
 
call MPI_ALLREDUCE(MPI_IN_PLACE, vCM, size(vCM), MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, mm,1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
GvCM = vCM
Gmm = mm

!--- get the global momentum
GvCM(:)=GvCM(:)/Gmm

!--- set the total momentum to be zero and get the current kinetic energy. 
KE = 0.d0
do i=1, NATOMS
   v(i,1:3) = v(i,1:3) - GvCM(1:3)

   ity = nint(atype(i))
   KE = KE + hmas(ity)*sum( v(i,1:3)*v(i,1:3) )
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, KE, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
GKE = KE/GNATOMS

!--- scale the obtained velocity to get the initial temperature.
vfactor = sqrt(1.5d0*treq/GKE)
v(:,:) = vfactor * v(:,:)

end subroutine

!------------------------------------------------------------------------------------------
subroutine CUTOFFLENGTH()
use atoms; use parameters
!------------------------------------------------------------------------------------------
implicit none
integer :: ity,jty,inxn
real(8) :: dr,BOsig

!--- cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM()
cutoff_vpar30 = cutof2_bo*vpar30

!--- get the cutoff length based on sigma bonding interaction.

! --- Remark --- 
! sigma bond before correction is the longest, namely longer than than pi and double-pi bonds
! thus check only sigma bond convergence.
allocate(rc(nboty), rc2(nboty), stat=ast)

rc(:)=0.d0; rc2(:)=0.d0
do ity=1,nso
do jty=ity, nso

   inxn=inxn2(ity,jty)
   if(inxn==0) cycle

   dr = 1.0d0
   BOsig=1.d0
   do while (BOsig > MINBOSIG) 
      dr = dr + 0.01d0
      BOsig = exp( pbo1(inxn)*(dr/r0s(ity,jty))**pbo2(inxn) ) !<- sigma bond prime
   enddo

   rc(inxn)  = dr
   rc2(inxn) = dr*dr
enddo
enddo

!----------------------------------------------------------------------------
! In some cases, an atom that do not exist in simulation gives
! the longest bond-order cutoff length. the check below is to ignore
! such atoms to keep the linkedlist cell dimensions as small as possible.
!----------------------------------------------------------------------------
do ity=1, nso
   if(natoms_per_type(ity)==0) then
      do jty=1, nso
         inxn=inxn2(ity,jty)
         if(inxn/=0) rc(inxn)=0.d0
         inxn=inxn2(jty,ity)
         if(inxn/=0) rc(inxn)=0.d0
      enddo
   endif
enddo

!--- get the max cutoff length 
maxrc=maxval(rc(:))

end subroutine

!------------------------------------------------------------------------------------------
subroutine POTENTIALTABLE()
use atoms; use parameters; use MemoryAllocator
!------------------------------------------------------------------------------------------
implicit none
integer :: i, ity,jty,inxn
real(8) :: dr1, dr2, dr3, dr4, dr5, dr6, dr7

real(8) :: exp1, exp2
real(8) :: dr3gamij, gamwinvp, gamWij, alphaij, Dij0, rvdW0
real(8) :: Tap, dTap, fn13, dfn13, rij_vd1

real(8) :: dr_lg, dr6_lg, Elg, E_core, dE_core, dElg

real(8) :: clmb,dclmb,erf_alphaijr,derf_alphaijr
real(8) :: alpha_i, alpha_j


!--- first element in table 0: potential
!---                        1: derivative of potential
call allocatord3d(TBL_EClmb,0,1,1,NTABLE,1,nboty)
call allocatord3d(TBL_Evdw,0,1,1,NTABLE,1,nboty)
call allocatord2d(TBL_EClmb_QEq,1,NTABLE,1,nboty)

!--- unit distance in r^2 scale
UDR = rctap2/NTABLE
UDRi = 1.d0/UDR

do ity=1, nso
do jty=ity, nso

   inxn = inxn2(ity,jty)
   if(inxn/=0) then
      do i=1, NTABLE

         dr2 = UDR*i
         dr1 = sqrt(dr2)

!--- Interaction Parameters:
         gamWij = gamW(ity,jty)
         alphaij = alpij(ity,jty)
         Dij0 = Dij(ity,jty)
         rvdW0 = rvdW(ity,jty) 
         gamwinvp = (1.d0/gamWij)**pvdW1

         dr3 = dr1*dr2
         dr4 = dr2*dr2
         dr5 = dr1*dr2*dr2
         dr6 = dr2*dr2*dr2
         dr7 = dr1*dr2*dr2*dr2 

         rij_vd1 = dr2**pvdW1h
         Tap = CTap(7)*dr7 + CTap(6)*dr6 + &
               CTap(5)*dr5 + CTap(4)*dr4 + CTap(0)

         fn13 = (rij_vd1 + gamwinvp)**pvdW1inv
         exp1 = exp( alphaij*(1.d0 - fn13 / rvdW0) )
         exp2 = sqrt(exp1)

         dr3gamij = ( dr3 + gamij(ity,jty) )**( -1.d0/3.d0 )

         TBL_Evdw(0,i,inxn) = Tap*Dij0*(exp1 - 2d0*exp2)
         TBL_Eclmb(0,i,inxn) = Tap*Cclmb*dr3gamij
         TBL_Eclmb_QEq(i,inxn) = Tap*Cclmb0_qeq*dr3gamij

!--- Force Calculation:
         dTap = 7d0*CTap(7)*dr5 + 6d0*CTap(6)*dr4 + &
                5d0*CTap(5)*dr3 + 4d0*CTap(4)*dr2

         dfn13 = ((rij_vd1 + gamwinvp)**(pvdW1inv-1.d0)) * (dr2**(pvdW1h-1.d0)) 


         TBL_Evdw(1,i,inxn) = Dij0*( dTap*(exp1 - 2.d0*exp2)  &
                            - Tap*(alphaij/rvdW0)*(exp1 - exp2)*dfn13 )
         TBL_Eclmb(1,i,inxn) = Cclmb*dr3gamij*( dTap - (dr3gamij**3)*Tap*dr1 )

         if(isLG) then

!FIXME LG extension supports only C,H,O,N at this moment. We should use element name instead of fixed numbers.  
            if (ity > 4 .or. jty > 4) cycle   

            dr_lg = 2*sqrt(Re_lg(ity)*Re_lg(jty))
            dr6_lg = dr_lg**6


            Elg = -C_lg(ity,jty)/(dr6 + dr6_lg)
            E_core = ecore(ity,jty)*exp(acore(ity,jty)*(1.d0-(dr1/rcore(ity,jty))))

            dElg = C_lg(ity,jty)*(6.d0*dr5)/(dr6 + dr6_lg)**2/dr1
            dE_core = -acore(ity,jty)*E_core/rcore(ity,jty)/dr1

            TBL_Evdw(0,i,inxn) = TBL_Evdw(0,i,inxn) + Tap*(Elg+E_core)
            TBL_Evdw(1,i,inxn) = TBL_Evdw(1,i,inxn) + dTap*Elg+Tap*dElg + dTap*E_core+Tap*dE_core

         endif

      enddo
   endif

enddo
enddo

end subroutine

!----------------------------------------------------------------
subroutine GetNonbondingMesh()
use atoms; use parameters; use MemoryAllocator
! setup 10[A] radius mesh to avoid visiting unecessary cells 
!----------------------------------------------------------------
implicit none

integer :: i,j,k

real(8) :: latticePerNode(3), rr(3), dr2
real(8) :: maxrcell
integer :: imesh(3), maximesh, ii(3), i1

!--- initial estimate of LL cell dims
nblcsize(1:3)=3d0

!--- get mesh resolution which is close to the initial value of rlc.
latticePerNode(1)=lata/vprocs(1)
latticePerNode(2)=latb/vprocs(2)
latticePerNode(3)=latc/vprocs(3)
nbcc(1:3)=int(latticePerNode(1:3)/nblcsize(1:3))
nblcsize(1:3)=latticePerNode(1:3)/nbcc(1:3)
maxrcell = maxval(nblcsize(1:3))

!--- get # of linked list cell to cover up the non-bonding cutoff length
imesh(1:3)  = int(rctap/nblcsize(1:3)) + 1
maximesh = maxval(imesh(1:3))

!--- List up only cell indices within the cutoff range.
!--- pre-compute nmesh to get exact array size.
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= rctap**2) nbnmesh = nbnmesh + 1
enddo; enddo; enddo

call allocatori2d(nbmesh,1,3,1,nbnmesh)

nbmesh(:,:)=0
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= rctap**2) then
      nbnmesh = nbnmesh + 1
      nbmesh(1:3,nbnmesh) = (/i, j, k/)
   endif
enddo; enddo; enddo

call allocatori1d(nbllist,1,NBUFFER)
call allocatori3d(nbheader, &
                -MAXLAYERS_NB,nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(3)-1+MAXLAYERS_NB)
call allocatori3d(nbnacell, &
                -MAXLAYERS_NB,nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(3)-1+MAXLAYERS_NB)

!--- normalize nblcsize, like lcsize.
nblcsize(1:3)=nblcsize(1:3)/(/lata,latb,latc/)

end subroutine

!----------------------------------------------------------------
subroutine GetBoxParams(H,la,lb,lc,angle1,angle2,angle3)
!----------------------------------------------------------------
implicit none
real(8),intent(out) :: H(3,3)
real(8),intent(in) :: la,lb,lc, angle1,angle2,angle3
real(8) :: hh1, hh2 , lal, lbe, lga
real(8) :: pi=atan(1.d0)*4.d0

!--- convet unit for angles
lal=angle1*pi/180.d0
lbe=angle2*pi/180.d0
lga=angle3*pi/180.d0

!--- construct H-matrix
hh1=lc*(cos(lal)-cos(lbe)*cos(lga))/sin(lga)
hh2=lc*sqrt( 1.d0-cos(lal)**2-cos(lbe)**2-cos(lga)**2 + &
             2*cos(lal)*cos(lbe)*cos(lga) )/sin(lga)

H(1,1)=la;          H(2,1)=0.d0;        H(3,1)=0.d0
H(1,2)=lb*cos(lga); H(2,2)=lb*sin(lga); H(3,2)=0.d0
H(1,3)=lc*cos(lbe); H(2,3)=hh1;         H(3,3)=hh2

return
end subroutine

!----------------------------------------------------------------
subroutine UpdateBoxParams()
use atoms; use parameters
!----------------------------------------------------------------
implicit none

!--- get volume 
MDBOX = &
HH(1,1,0)*(HH(2,2,0)*HH(3,3,0) - HH(3,2,0)*HH(2,3,0)) + &
HH(2,1,0)*(HH(3,2,0)*HH(1,3,0) - HH(1,2,0)*HH(3,3,0)) + &
HH(3,1,0)*(HH(1,2,0)*HH(2,3,0) - HH(2,2,0)*HH(1,3,0))

!--- get inverse of H-matrix
call matinv(HH,HHi)

!--- local box dimensions (a temporary use of lbox)
LBOX(1)=lata/vprocs(1)
LBOX(2)=latb/vprocs(2)
LBOX(3)=latc/vprocs(3)

!--- get the number of linkedlist cell per domain
cc(1:3)=int(LBOX(1:3)/maxrc)

!--- local system size in the unscaled coordinate.
LBOX(1:3) = 1.d0/vprocs(1:3)

!--- get the linkedlist cell dimensions (normalized)
lcsize(1:3) = LBOX(1:3)/cc(1:3)

!--- get origin of local MD box in the scaled coordiate.
OBOX(1:3) = LBOX(1:3)*vID(1:3)

return
end

!-------------------------------------------------------------------------------------------
module base
!-------------------------------------------------------------------------------------------
! position, atom type, velocity, force & charge
real(8),allocatable,dimension(:),target :: atype, q
real(8),allocatable,dimension(:,:),target :: pos, v, f
end module

!-------------------------------------------------------------------------------------------
module atoms
!-------------------------------------------------------------------------------------------
#ifdef NOMPI
  use nompi
#else
  include 'mpif.h'
#endif

!--- For array size statistics
!  1-NATOMS, 2-nbrlist, 3-nbrlist for qeq, 4-NBUFFER for move, 5-NBUFFER for copy
!  6-NBUFFER for qeq
integer,parameter :: nmaxas=5
integer,allocatable :: maxas(:,:)

!--- lattice parameters 
real(8) :: lata,latb,latc,lalpha,lbeta,lgamma

integer :: myid, nprocs, ierr, myparity(3), vID(3)

!<sbuffer> send buffer, <rbuffer> receive buffer
real(8),allocatable :: sbuffer(:), rbuffer(:)
   
!<ns> # of atoms to be sent, <nr> # of atoms to be received, <na> # of all of transfered atoms.
!<ne> # of elements for one atom.
!     Example) In case atom type, position and velocity to be sent,  ne = 1+3+3 = 7
integer :: ns, nr, na, ne

!<NE_COPY>,<NE_MOVE>,<NE_CPBK> :: Number of Elements to COPY, MOVE atoms and CoPy BacK force. 
integer,parameter :: MODE_COPY = 1, MODE_MOVE = 2, MODE_CPBK = 3
integer,parameter :: MODE_QCOPY1 = 4, MODE_QCOPY2 = 5

integer,parameter :: NE_CPBK = 4

!<MAXLAYERS> MAXimum # of linkedlist cell LAYERS.
integer,parameter :: MAXLAYERS=5
integer,parameter :: MAXLAYERS_NB=10
        
! <target_node> stores partner node ID in the 6-communications. 
! if targe_node(i)==-1, the node doesn't have a partner in i-direction.
integer :: target_node(6)

! For benchmarking, <vprocs> and <mc> will be read from vprocs.in
integer :: vprocs(3)

!<mc(3)> # of unit cells in each directions
integer :: mc(3)

real(8),allocatable:: rc(:), rc2(:)   !<RCUT>: cutoff length for sigma-bonding.
real(8),allocatable:: rcpi(:), rcpp(:)!      : cutoff length for other bonding.

real(8),parameter :: MINBOSIG = 1d-3      !<minBOsig>: criterion to decide <rc> 
real(8),parameter :: MINBO0 = 1d-4       !<minBO0>: cutoff bond order 
real(8),parameter :: cutof2_esub = 1d-4
!real(8),parameter :: cutof2_bo = 1.d-2
real(8),parameter :: cutof2_bo = 1.d-3
integer,parameter :: is_idEh = 1

!real(8),parameter :: MINBOSIG = 1d-4      !<minBOsig>: criterion to decide <rc> 
!real(8),parameter :: MINBO0 = 0.d0       !<minBO0>: cutoff bond order 
!real(8),parameter :: cutof2_esub = 0.d0
!real(8),parameter :: cutof2_bo = 1d-4
!integer,parameter :: is_idEh = 0

! cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM()
real(8) :: cutoff_vpar30

!integer :: NBUFFER=5000
!integer,parameter :: MAXNEIGHBS=50  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
!integer,parameter :: MAXNEIGHBS10=200 !<MAXNEIGHBS>: Max # of Ngbs within the taper function cutoff. 

integer :: NBUFFER=30000
integer,parameter :: MAXNEIGHBS=30  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
integer,parameter :: MAXNEIGHBS10=1500 !<MAXNEIGHBS>: Max # of Ngbs within 10[A]. 

integer,parameter :: NMINCELL=4  !<NMINCELL>: Nr of minimum linkedlist cell <-> minimum grain size.
real(8),parameter :: MAXANGLE= 0.999999999999d0 
real(8),parameter :: MINANGLE=-0.999999999999d0
real(8),parameter :: NSMALL = 1.d-10
real(8) :: maxrc                        !<maxRCUT>: Max cutoff length. used to decide lcsize.

real(8),parameter :: pi=3.14159265358979d0
real(8),parameter :: sqrtpi_inv=1.d0/sqrt(pi)

! stress tensor
real(8) :: astr(6) 

!--- coefficient of bonding energy derivative 
real(8),allocatable :: ccbnd(:)
real(8),allocatable :: cdbnd(:)

real(8) :: HH(3,3,0:1), HHi(3,3), MDBOX, LBOX(0:3), OBOX(1:3) !MD box, local MD box, origin of box.
integer :: NATOMS         !local # of atoms
integer(8) :: GNATOMS     !global # of atoms
integer :: ALLATOMS

!<llist> Linked List
!<header> header atom of linkedlist cell.
!<nacell> Nr of atoms in a likedlist cell.
integer,allocatable :: llist(:), header(:,:,:), nacell(:,:,:)

!<nbllist> Linked List for non-bonding interaction
!<nbheader> header atom of linkedlist cell for non-bonding interaction
!<nbnacell> Nr of atoms in a likedlist cell for non-bonding interaction
integer,allocatable :: nbllist(:), nbheader(:,:,:), nbnacell(:,:,:)

!<nbrlist> neighbor list, <nbrindx> neighbor index
integer,allocatable :: nbrlist(:,:), nbrindx(:,:)

!<nbplist> neighbor list of nonbonding interaction, non-bonding pair list
integer,allocatable :: nbplist(:,:)

!<BO> Bond Order of atoms i-j (nearest neighb only) - (Eq 3a-3d)
real(8),allocatable :: BO(:,:,:) 

real(8),allocatable :: delta(:)

!--- Output variables from the BOp_CALC() subroutine:
real(8),allocatable :: deltap(:,:)

real(8),allocatable :: dln_BOp(:,:,:)

real(8),allocatable :: dBOp(:,:)

!--- For NEW DBO calc:
real(8),allocatable :: exp_delt1(:,:), exp_delt2(:,:)  ! exp( -pboc#(inxn) * deltap(i,1) ) - {inxn, i}   

!--- A[0123] coefficients for force calculation 
real(8),allocatable :: A0(:,:),A1(:,:), A2(:,:), A3(:,:) 

!--- Passed between Elnpr and E3body
real(8),allocatable :: nlp(:), dDlp(:) !Number of Lone Pairs, its derivatives.
real(8),allocatable :: deltalp(:)

! TE: Total Energy,  KE: Kinetic Energy,  PE :: Potential Energies
!  0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
!  7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
real(8) :: TE, KE, PE(0:13)
real(8) :: GTE, GKE, GPE(0:13)

!--- output file format 
logical :: isBinary=.false., isBondFile=.false., isPDB=.false., isXYZ=.false.

!--- one vector charge equlibration
! g: gradient,  h: conjugate direction, hsh: hessian * h
!real(8),allocatable :: g(:), h(:), hsh(:)

! Two vectors electrostatic energy minimization 
real(8),allocatable,target :: qs(:),qt(:),gs(:), gt(:), hs(:), ht(:), hshs(:), hsht(:)

!--- variables for extended Lagrangian method ---
!<Lex_fqs> fraction between two QEq vectors
!<Lex_w> spring constant
real(8),allocatable,target :: qsfp(:),qsfv(:),qtfp(:),qtfv(:) 
real(8),allocatable :: hessian(:,:)
real(8) :: Lex_fqs=1.0, Lex_w=1.d0, Lex_w2=1.d0, Lex_k=2.d0
 
integer :: ast ! Allocation STatus of allocatable variables.

! <lcsize> Linked list Cell SIZE. <cc> Nr of likedlist cell in local node.
real(8) :: lcsize(3), nblcsize(3)
integer :: cc(3), nbcc(3)

integer :: nmesh, nbnmesh
integer,allocatable :: mesh(:,:), nbmesh(:,:)


!--- Unit convetors. In original ReaxFF program, the units of length is [A]
!--- energy is [kcal/mol] and mass is [amu] respectively.
!--- Most of numbers written here are obtained below site.
!--- http://physics.nist.gov/cuu/Constants/Table/allascii.txt
!--- Bohr length
real(8),parameter :: Lbohr_a  = 0.5291772108d0   ! [A]
real(8),parameter :: Lbohr_m  = 0.5291772108d-10 ! [m]

!--- Electron rest mass
real(8),parameter :: Merest_amu = 5.48580152d-4  ! [amu]
real(8),parameter :: Merest_kg  = 9.10938d-31    ! [kg]
!--- Energy in Hartree
real(8),parameter :: Ehrtre_km = 6.2751d2        ! [kcal/mol]
real(8),parameter :: Ehrtre_ev  = 27.2113845d0   ! [eV]
real(8),parameter :: Ehrtre_j = 4.35974417d-18   ! [J] 
real(8),parameter :: Eev_kcal = 23.060538d0      ! [kcal/mol]

real(8),parameter :: Ekcal_j = 6.95016611d-21  ! [J]

!--- Boltzmann Constant
real(8),parameter :: BLTZMN = 1.3806503d-23  ! [m^2 kg s^-2 K^-1 ] 

real(8),parameter :: UTEMP0 = 503.398008d0    ! Ekcal_j/BLZMN [K]
real(8),parameter :: UTEMP = UTEMP0*2.d0/3.d0 ! [K]
real(8),parameter :: USTRS = 6.94728103d0     ! [GPa]
real(8),parameter :: UDENS = 1.66053886d0     ! [g/cc]
real(8),parameter :: UTIME = 1.d3/20.455d0    ! 1 = 1/20.445[ps] = 48.88780[fs]

!--- QEq variables. 
!<isQEq> flag to run QEq routine: 0-No QEq, 1-CG, 2-Extended Lagrangian
integer :: isQEq
!<NMAXQEq> Number of MAXimum iteration in QEq routine
integer :: NMAXQEq
!<QEq_thrsld> energy criterion in QEq routine
real(8) :: QEq_tol
!<nstep_qeq> counter of iteration
integer :: nstep_qeq, qstep

!-- variables for timing
integer,parameter :: Ntimer=30
integer :: it_timer(Ntimer)=0, it_timer_max(Ntimer)=0, it_timer_min(Ntimer)=0
real(8) :: wt0=0.d0

!---
! <mdmode> determines MD mode
integer :: mdmode
! <nstep> current MD step, <ntime_step> Total # of time steps in one MD run.
! <current_step> will be used for subsequent runs.
integer :: nstep=0, ntime_step, current_step
!<vsfact> velocity scaling factor, <dt> one time step
real(8) :: treq, vsfact, dt, dmt
integer :: sstep

!--- output format flags, explained in 'rxmdopt.in'
integer :: fstep, pstep

!--- <frcindx> FoRCe INDeX. Index to return calculated force to original atoms.
real(8),allocatable,target :: frcindx(:)
integer :: copyptr(0:6)

!--- stress components
real(8) :: xx,yy,zz,yz,zx,xy
!--- atomic stress index
integer :: ia,ja

!--- conjugate_gradient
REAL(8) :: ftol   ! tolerance of energy convergence 

!--- cutoff range calculation. 
integer(8),allocatable :: natoms_per_type(:)

!--- dthm=dt/(2*mass), hmas=mass/2
real(8),allocatable :: dthm(:), hmas(:)

!--- potential teble
integer,parameter :: NTABLE=5000
real(8),allocatable :: TBL_Eclmb(:,:,:), TBL_Evdw(:,:,:), TBL_Eclmb_QEq(:,:)
real(8) :: UDR, UDRi

integer(8),allocatable :: ibuf8(:)

integer,parameter :: MAXSTRLENGTH=256

!--- FF parameter description
character(MAXSTRLENGTH) :: FFDescript

!--- from cmdline_args
logical :: isSpring=.false.
real(8) :: springConst=0.d0
logical :: hasSpringForce(16)
real(8),allocatable :: ipos(:,:)

logical :: isFF=.false., isData=.false., isMDparm=.false., isRunFromXYZ=.false.
character(MAXSTRLENGTH) :: FFPath="ffield", DataDir="DAT", ParmPath="rxmd.in", RunFromXYZPath=""

logical :: saveRunProfile=.false.
character(MAXSTRLENGTH) :: RunProfilePath="profile.dat"
integer,parameter :: RunProfileFD=30 ! file descriptor for summary file

logical :: isLG=.false.

logical :: isEfield=.false.


!--- Taper function 
real(8),parameter :: rctap0 = 10.d0 ![A]
real(8),parameter :: rctap0_pqeq = 12.5d0 ![A]
real(8) :: rctap, rctap2, CTap(0:7)

!--- PQEQ variables
real(8),allocatable,target :: spos(:,:)
integer,allocatable :: inxnpqeq(:,:)

real(8),allocatable :: TBL_Eclmb_pcc(:,:,:),TBL_Eclmb_psc(:,:,:),TBL_Eclmb_pss(:,:,:)

integer :: ntype_pqeq, ntype_pqeq2
logical,allocatable :: isPolarizable(:)
character(len=:),allocatable :: Elempqeq(:)
real(8),allocatable :: X0pqeq(:),J0pqeq(:)

real(8),allocatable :: Zpqeq(:), Rcpqeq(:), Rspqeq(:), Kspqeq(:)
real(8),allocatable :: alphacc(:,:), alphasc(:,:), alphass(:,:)
real(8) :: lambda_pqeq = 0.462770d0

real(8) :: eFieldStrength=0.0
integer :: eFieldDir=1

logical :: isPQEq = .false.
character(MAXSTRLENGTH) :: PQEqParmPath

contains

!-----------------------------------------------------------------------------------------------------------------------
character(len=256) function rankToString(irank)
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: irank
write(rankToString,*) irank
rankToString = adjustl(rankToString)
end function rankToString

!-------------------------------------------------------------------------------------------
function GetFileNameBase(DataDir,nstep) result(fileNameBase)
!-------------------------------------------------------------------------------------------
integer,intent(in) :: nstep
character(MAXSTRLENGTH) :: DataDir,fileNameBase
character(9) :: a9

if(nstep>=0) then
  write(a9,'(i9.9)') nstep
  fileNameBase=trim(adjustl(DataDir))//"/"//a9
else
  fileNameBase=trim(adjustl(DataDir))//"/rxff"
endif

end function

end module atoms

!-------------------------------------------------------------------------------------------
module pqeq_vars
use atoms
!-------------------------------------------------------------------------------------------
implicit none

contains

!-------------------------------------------------------------------------------------------
subroutine initialize_eField(myid)
implicit none
!-------------------------------------------------------------------------------------------
integer,intent(in) :: myid

if(myid==0) then
   print'(a)','-----------------------------------------------------------'
   print'(a,f12.6,i6)','eField [V/A], eFiled direction : ', eFieldStrength, eFieldDir
   print'(a)','-----------------------------------------------------------'
endif

return
end subroutine

!-------------------------------------------------------------------------------------------
subroutine EEfield(Etotal,NATOMS,pos,q,f,atype,Eev_kcal)
implicit none 
!-------------------------------------------------------------------------------------------
integer,intent(in) :: NATOMS
real(8),intent(in) :: pos(NATOMS,3),q(NATOMS),atype(NATOMS),Eev_kcal

real(8) :: f(NATOMS,3), Etotal

integer :: i, ity
real(8) :: Eenergy, Eforce, qic

! NOTE: the energy function from an eField is to be determined and not computed here. 
do i=1, NATOMS

   ity = nint(atype(i))
   qic = q(i) + Zpqeq(ity)

   Eforce = -qic*eFieldStrength*Eev_kcal

   f(i,eFieldDir) = f(i,eFieldDir) + Eforce

enddo

return
end subroutine

!-------------------------------------------------------------------------------------------
subroutine get_coulomb_and_dcoulomb_pqeq(rr,alpha,Eclmb,inxn,TBL_Eclmb_PQEq,ff)
implicit none
!-------------------------------------------------------------------------------------------
real(8),intent(in) :: rr(3),alpha,TBL_Eclmb_PQEq(ntype_pqeq2,NTABLE,0:1)
integer,intent(in) :: inxn
real(8) :: Ecc,Esc,Ess,dEcc,dEsc,dEss,Eclmb,dEclmb,ff(3)

integer :: i,itb, itb1
real(8) :: drtb, drtb1, dr1, dr2, dr1i 
real(8) :: Tap, dTap, clmb, dclmb, screen, dscreen 

real(8) :: dr3,dr4,dr5,dr6,dr7

dr2 = sum(rr(1:3)*rr(1:3))

! note that shell-shell distance could be greater than cutoff though core-core is within the cutoff.
if(dr2>rctap2) return 

dr1 = sqrt(dr2)

itb = int(dr2*UDRi)
itb1 = itb+1
drtb = dr2 - itb*UDR
drtb = drtb*UDRi
drtb1= 1.d0-drtb

Eclmb = drtb1*TBL_Eclmb_PQEq(inxn,itb,0)  + drtb*TBL_Eclmb_PQEq(inxn,itb1,0)
dEclmb = drtb1*TBL_Eclmb_PQEq(inxn,itb,1)  + drtb*TBL_Eclmb_PQEq(inxn,itb1,1)

ff(1:3)=dEclmb*rr(1:3)

return

dr1i = 1.d0/dr1
clmb = dr1i
!TODO 1. Tabulate the coulomb term, 2. change the derivatives to our convention, i.e. u(r)'/r^2*dr(1:3)
!dclmb = -dr1i*dr1i*dr1i 
dclmb = -dr1i*dr1i

screen = erf(alpha*dr1)
!dscreen = 2.d0*alpha*sqrtpi_inv*exp(-alpha*alpha*dr2)*dr1i
dscreen = 2.d0*alpha*sqrtpi_inv*exp(-alpha*alpha*dr2)

!--- core-core distance is withing the taper cutoff, but core-shell & shell-shell distance can be beyond the cutoff.
!--- Directly computing the taper function for now. 
dr3 = dr1*dr2
dr4 = dr2*dr2
dr5 = dr1*dr2*dr2
dr6 = dr2*dr2*dr2
dr7 = dr1*dr2*dr2*dr2
Tap = CTap(7)*dr7 + CTap(6)*dr6 + CTap(5)*dr5 + CTap(4)*dr4 + CTap(0)
dTap = 7d0*CTap(7)*dr6 + 6d0*CTap(6)*dr5 + 5d0*CTap(5)*dr4 + 4d0*CTap(4)*dr3

Eclmb = clmb*screen*Tap
dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

ff(1:3)=dEclmb*rr(1:3)/dr1

return
end subroutine

!-------------------------------------------------------------------------------------------
subroutine set_alphaij_pqeq()
implicit none
!-------------------------------------------------------------------------------------------
integer :: ity,jty

real(8) :: alpha_ci,alpha_cj,alpha_si,alpha_sj 

alphacc(:,:)=0.d0
alphasc(:,:)=0.d0
alphass(:,:)=0.d0

do ity=1,ntype_pqeq

   alpha_ci=0.5d0*lambda_pqeq/Rcpqeq(ity)**2
   alpha_si=0.5d0*lambda_pqeq/Rspqeq(ity)**2

   do jty=1,ntype_pqeq

      alpha_cj=0.5d0*lambda_pqeq/Rcpqeq(jty)**2
      alpha_sj=0.5d0*lambda_pqeq/Rspqeq(jty)**2

      ! core(i)-core(j) term.
      alphacc(ity,jty)=sqrt( (alpha_ci*alpha_cj)/(alpha_ci + alpha_cj) )

      ! shell(i)-shell(j) term.
      if( isPolarizable(ity) .and. isPolarizable(jty) ) &
          alphass(ity,jty)=sqrt( (alpha_si*alpha_sj)/(alpha_si + alpha_sj) )

      ! shell(i)-core(j) term. the first index must be polarizable atom type. 
      ! C(r_si_cj) is fine but use alphasc(jty,ity) when refer to alphasc for C(r_ci_sj).
      if( isPolarizable(ity) ) &
          alphasc(ity,jty)=sqrt( (alpha_si*alpha_cj)/(alpha_si + alpha_cj) )

   enddo

enddo

end subroutine

!-------------------------------------------------------------------------------------------
subroutine initialize_pqeq(chi,eta)
implicit none
!-------------------------------------------------------------------------------------------
integer :: ity,jty,icounter
real(8) :: chi(ntype_pqeq),eta(ntype_pqeq)

integer :: i,inxn
real(8) :: Acc,Asc,Ass
real(8) :: dr1,dr2,dr3,dr4,dr5,dr6,dr7,dr1i,UDR,UDRi
real(8) :: clmb,dclmb,screen,dscreen,Tap,dTap,Eclmb,dEclmb

call set_alphaij_pqeq()

!--- for PQEq
do ity = 1, ntype_pqeq
  if( .not. isPolarizable(ity) ) then
    if(myid==0) &
       print'(a,i3,a)','atom type ', ity, ' is not polarizable. Setting Z & K to zero.'
    Zpqeq(ity)=0.d0
    Kspqeq(ity)=0.d0
  else
    if(myid==0) then
       print'(a $)','updating chi and eta '

       print'(a,i3,1x,a2,1x,4f12.6)', &
         'name,chi,X0pqeq,eta,J0pqeq :', &
          ity,Elempqeq(ity), chi(ity),X0pqeq(ity), eta(ity),J0pqeq(ity)
    endif

    chi(ity)=X0pqeq(ity)
    eta(ity)=J0pqeq(ity)
  endif
enddo

!--- 2x factor in param.F90
eta(:)=2.d0*eta(:)

allocate(inxnpqeq(ntype_pqeq,ntype_pqeq))

icounter=0
do ity=1, ntype_pqeq
do jty=ity, ntype_pqeq

   icounter = icounter + 1

   inxnpqeq(ity,jty)=icounter
   inxnpqeq(jty,ity)=inxnpqeq(ity,jty)
enddo; enddo

allocate(TBL_Eclmb_pcc(ntype_pqeq2,NTABLE,0:1))
allocate(TBL_Eclmb_psc(ntype_pqeq2,NTABLE,0:1))
allocate(TBL_Eclmb_pss(ntype_pqeq2,NTABLE,0:1))

!--- unit distance in r^2 scale
UDR = rctap2/NTABLE
UDRi = 1.d0/UDR

do ity=1, ntype_pqeq
do jty=ity, ntype_pqeq

   Acc = alphacc(ity,jty)
   Asc = alphasc(ity,jty)
   Ass = alphass(ity,jty)

   inxn = inxnpqeq(ity,jty)

   do i=1,NTABLE

         dr2 = UDR*i
         dr1 = sqrt(dr2)

!--- Interaction Parameters:
         dr3 = dr1*dr2
         dr4 = dr2*dr2
         dr5 = dr1*dr2*dr2
         dr6 = dr2*dr2*dr2
         dr7 = dr1*dr2*dr2*dr2 

         Tap = CTap(7)*dr7 + CTap(6)*dr6 + &
               CTap(5)*dr5 + CTap(4)*dr4 + CTap(0)

!--- Force Calculation:
         dTap = 7d0*CTap(7)*dr5 + 6d0*CTap(6)*dr4 + &
                5d0*CTap(5)*dr3 + 4d0*CTap(4)*dr2

         dr1i = 1.d0/dr1

         clmb = dr1i
         dclmb = -dr1i*dr1i*dr1i 

         !--- core-core interaction
         screen = erf(Acc*dr1)
         dscreen = 2.d0*Acc*sqrtpi_inv*exp(-Acc*Acc*dr2)*dr1i

         Eclmb = clmb*screen*Tap
         dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

         TBL_Eclmb_pcc(inxn,i,0) = Eclmb
         TBL_Eclmb_pcc(inxn,i,1) = dEclmb

         !--- core-shell interaction
         screen = erf(Asc*dr1)
         dscreen = 2.d0*Asc*sqrtpi_inv*exp(-Asc*Asc*dr2)*dr1i

         Eclmb = clmb*screen*Tap
         dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

         TBL_Eclmb_psc(inxn,i,0) = Eclmb
         TBL_Eclmb_psc(inxn,i,1) = dEclmb

         !--- shell-shell interaction
         screen = erf(Ass*dr1)
         dscreen = 2.d0*Ass*sqrtpi_inv*exp(-Ass*Ass*dr2)*dr1i

         Eclmb = clmb*screen*Tap
         dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

         TBL_Eclmb_pss(inxn,i,0) = Eclmb
         TBL_Eclmb_pss(inxn,i,1) = dEclmb

   enddo

enddo
enddo

end subroutine

end module
!-------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------
module parameters
!-------------------------------------------------------------------------------------------
!  This module stores the parameters input from the rxmda.in file.  This file is 
!    intentionally very similar to the input setup used by Adri van Duin at Caltech (to aid
!    in cross-checking and updates).  However, please note that there are a few changes made
!    aside from the obvious f77 -> f90 switch.  I have tried to clearly note these. 
!-------------------------------------------------------------------------------------------
!use cmdline_args

!Independant Parameters

integer :: nso    !Number of different types of atoms
integer :: nboty  !Number of different bonds given

! Atom Dependant (ie where they appear in input file - not implementation in code)
character(2),allocatable :: atmname(:)      !holds the Chemical Abbrev for each atomtype
real(8),allocatable :: Val(:),Valboc(:)  !Valency of atomtype (norm, boc) 
real(8),allocatable :: mass(:)           !mass of atomtype

real(8),allocatable :: pbo1(:), pbo2(:), pbo3(:)   !Bond Order terms
real(8),allocatable :: pbo4(:), pbo5(:), pbo6(:)   !Bond Order terms
real(8),allocatable :: pboc1(:), pboc2(:), pboc3(:)  !Bond Order correction terms (f1-5)
real(8),allocatable :: pboc4(:), pboc5(:)            !Bond Order correction terms (f1-5)  
real(8),allocatable :: v13cor(:) !<kn>

real(8),allocatable :: rat(:),rapt(:),vnq(:)   !r0s/r0p/r0pp for like bonds 
real(8),allocatable :: ovc(:)  !a flag to apply fn4 and fn5 !<kn>

real(8),allocatable :: Desig(:),Depi(:),Depipi(:)  !Bond Energy parameters (eq. 6)
real(8),allocatable :: pbe1(:),pbe2(:)             !Bond Energy parameters (eq. 6) 

real(8),allocatable :: Vale(:)                      !Lone Pair Energy parameters (eq. 7)
real(8),allocatable :: plp1(:), nlpopt(:), plp2(:)  !Lone Pair Energy parameters (eq.8-10)  

real(8),allocatable :: povun1(:), povun2(:), povun3(:), povun4(:)   !Overcoordination Energy (eq. 11)
real(8),allocatable :: povun5(:), povun6(:), povun7(:), povun8(:)   !Undercoordination Energy (eq. 12)

real(8),allocatable :: pval1(:), pval2(:), pval3(:), pval4(:), pval5(:)   !Valency Angle Energy (eq. 13a-g)
real(8),allocatable :: pval6(:), pval7(:), pval8(:), pval9(:), pval10(:)   
real(8),allocatable :: Valangle(:), theta00(:)

real(8),allocatable :: ppen1(:), ppen2(:), ppen3(:), ppen4(:)   !Penalty Energy (eq. 14ab)

real(8),allocatable :: pcoa1(:), pcoa2(:), pcoa3(:), pcoa4(:)   !Conjugation (3 body) Energy (eq.15)
real(8),allocatable :: Valval(:)

real(8),allocatable :: ptor1(:), ptor2(:), ptor3(:), ptor4(:)   !Torsional Energy Terms (eq.16abc)
real(8),allocatable :: V1(:), V2(:), V3(:)

real(8),allocatable :: pcot1(:), pcot2(:)  !Conjugation (4body) Energy (eq. 17ab)

real(8),allocatable :: phb1(:), phb2(:), phb3(:), r0hb(:)   !Hydrogren Bond Energy (eq. 18)

real(8),allocatable :: Dij(:,:), alpij(:,:), rvdW(:,:), gamW(:,:)  !Van der Waals Energy (eq. 21ab)
real(8) :: pvdW1, pvdW1h, pvdW1inv

! hydrogen bonding interaction cutoff
real(8),parameter :: rchb = 10.d0 ![A]
real(8),parameter :: rchb2 = rchb*rchb

!Coulomb Energy (eq. 22)  
real(8),parameter:: Cclmb0 = 332.0638d0 ! [kcal/mol/A] line 2481 in poten.f
real(8),parameter:: Cclmb0_qeq = 14.4d0 ! [ev]
real(8),parameter:: CEchrge = 23.02d0   ! [ev]
real(8) :: Cclmb = Cclmb0

real(8),allocatable :: gam(:), gamij(:,:)  

!Charge Equilibration part, <chi>  electronegativity   <eta> stiffness
real(8),allocatable :: chi(:), eta(:)

!End Lost Parameters Listing

!Not Understood Parameters: 
real(8),allocatable :: bom(:)

! 2-atom combo dependant: 
real(8),allocatable :: r0s(:,:)       !Bond Order terms
real(8),allocatable :: r0p(:,:)       !  "" 
real(8),allocatable :: r0pp(:,:)      !Bond Order terms 

! <inxn2> is type of 2-body interaction 1=C-C, 2=H-C, 3=H-H
integer,allocatable :: inxn2(:,:)

integer,allocatable :: inxn3(:,:,:)
integer,allocatable :: inxn3hb(:,:,:)
integer,allocatable :: inxn4(:,:,:,:)

!Saved calculations to prevent having to recalc lots of times.
real(8),allocatable :: cBOp1(:), cBOp3(:), cBOp5(:)
real(8),allocatable :: pbo2h(:), pbo4h(:), pbo6h(:)  

!NOTE: for debugging purpose variables
real(8)  :: vpar30,vpar1,vpar2

!--- <switch> flag to omit pi and double pi bond in bond-order prime calculation.
real(8),allocatable :: switch(:,:) 

!--- LG params 
real(8), allocatable :: C_lg(:,:), Re_lg(:)
real(8), allocatable :: rcore2(:),ecore2(:),acore2(:)
real(8), allocatable :: rcore(:,:),ecore(:,:),acore(:,:)

end module parameters 

!-------------------------------------------------------------------------------------------
module MemoryAllocator
!-------------------------------------------------------------------------------------------
implicit none
integer(8) :: totalMemory=0

contains 

subroutine AllocatorD1D(array, imin, imax)
  implicit none
  integer,intent(in) :: imin, imax
  real(8),allocatable,dimension(:) :: array
  integer :: status
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*8

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorD1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorD1D(array)
  implicit none
  real(8),allocatable,dimension(:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorD1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorD2D(array, imin1, imax1, imin2, imax2) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2
  real(8),allocatable,dimension(:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*8

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorD2D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorD2D(array) 
  implicit none
  real(8),allocatable,dimension(:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorD2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorD3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  real(8),allocatable,dimension(:,:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*8

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorD3D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorD3D(array) 
  implicit none
  real(8),allocatable,dimension(:,:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorD3D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI1D(array, imin, imax) 
  implicit none
  integer,intent(in) :: imin, imax
  integer,allocatable,dimension(:) :: array
  integer :: status
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*4

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorI1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI1D(array) 
  implicit none
  integer,allocatable,dimension(:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorI1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI2D(array, imin1, imax1, imin2, imax2) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2
  integer,allocatable,dimension(:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*4

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorI2D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI2D(array) 
  implicit none
  integer,allocatable,dimension(:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorI2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  implicit none
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  integer,allocatable,dimension(:,:,:) :: array
  integer :: status
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*4

  if(status/=0) print'(a30,i9,i3)', 'ERROR in AllocatorI3D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI3D(array) 
  implicit none
  integer,allocatable,dimension(:,:) :: array
  integer :: status
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a30,i9,i3)', 'ERROR in DeallocatorI3D: totalMemory = ', totalMemory, status

  return
end subroutine 

integer(8) function GetTotalMemory() 
  GetTotalMemory = totalMemory
  return
end function

end module MemoryAllocator
!-------------------------------------------------------------------------------------------


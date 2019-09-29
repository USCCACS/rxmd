module mpi_mod

#ifdef NOMPI
  use nompi
#else
  include 'mpif.h'
#endif

end module

!-------------------------------------------------------------------------------------------
module base
use utils, only : MAXSTRLENGTH
!-------------------------------------------------------------------------------------------

type :: force_field_class
end type

type :: mdbase_class
   class(force_field_class),pointer :: ff
end type

!integer :: NBUFFER=5000
!integer,parameter :: MAXNEIGHBS=50  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
!integer,parameter :: MAXNEIGHBS10=200 !<MAXNEIGHBS>: Max # of Ngbs within the taper function cutoff. 

integer :: NBUFFER=10000
integer,parameter :: MAXNEIGHBS=50  ! Max # of Ngbs one atom may have. 
integer,parameter :: MAXNEIGHBS10=100 ! Max # of Ngbs for nonbonding pair. 

!<NE_COPY>,<NE_MOVE>,<NE_CPBK> :: Number of Elements to COPY, MOVE atoms and CoPy BacK force. 
integer,parameter :: MODE_COPY = 1, MODE_MOVE = 2, MODE_CPBK = 3
integer,parameter :: MODE_QCOPY1 = 4, MODE_QCOPY2 = 5

integer,parameter :: MODE_COPY_FNN=11, MODE_MOVE_FNN=12

integer,parameter :: NE_CPBK = 4

!<MAXLAYERS> MAXimum # of linkedlist cell LAYERS.
integer,parameter :: MAXLAYERS=5


! position, atom type, velocity, force & charge
real(8),allocatable,dimension(:),target :: atype, q
real(8),allocatable,dimension(:,:),target :: pos, v, f

character(2),allocatable :: atmname(:)  
real(8),allocatable :: mass(:)          

!<nbrlist> neighbor list
integer,allocatable :: nbrlist(:,:)

integer :: myid, nprocs, vprocs(3), ierr, myparity(3), vID(3)
! <target_node> stores partner node ID in the 6-communications. 
! if targe_node(i)==-1, the node doesn't have a partner in i-direction.
integer :: target_node(6)

!<rc>: cutoff length for the primary cutoff. <maxrc> max cutoff length used to decide lcsize.
real(8),allocatable:: rc(:), rc2(:) 
real(8) :: maxrc                    

!<llist> Linked List
!<header> header atom of linkedlist cell.
!<nacell> Nr of atoms in a likedlist cell.
integer,allocatable :: llist(:), header(:,:,:), nacell(:,:,:)
 
! <lcsize> Linked list Cell SIZE. <cc> Nr of likedlist cell in local node.
real(8) :: lcsize(3)
integer :: cc(3)

!--- lattice parameters 
real(8) :: lata,latb,latc,lalpha,lbeta,lgamma
real(8) :: hh(3,3,0:1), hhi(3,3), mdbox, lbox(3), obox(1:3) !MD box, local MD box, origin of box.

integer :: NATOMS         !local # of atoms
integer(8) :: GNATOMS     !global # of atoms

!--- <frcindx> FoRCe INDeX. Index to return calculated force to original atoms.
real(8),allocatable,target :: frcindx(:)
integer :: copyptr(0:6)

!--- total, kinetic and potential energies (local and global)
real(8) :: TE, KE, PE0, GTE, GKE, GPE0

!--- RXMD parameters
! <mdmode> determines MD mode
integer :: mdmode
! <nstep> current MD step, <ntime_step> Total # of time steps in one MD run.
! <current_step> will be used for subsequent runs.
integer :: nstep=0, ntime_step, current_step
!<vsfact> velocity scaling factor, <dt> one time step
real(8) :: treq, vsfact, dt, dmt
integer :: sstep
!--- output file format 
logical :: isBinary=.false., isBondFile=.false., isPDB=.false., isXYZ=.false.
!--- output format flags, explained in 'rxmdopt.in'
integer :: fstep, pstep
!--- conjugate_gradient tolerance of energy convergence 
REAL(8) :: ftol   


character(len=:),allocatable :: forcefield_type
logical :: is_reaxff=.false., is_fnn=.false.

!--- natoms_per_type() is used to clear unused cutoff distance for ReaxFF
!--- see get_cutoff_bondorder()
integer,parameter :: MAX_ATOMTYPE=16
integer(8) :: natoms_per_type(MAX_ATOMTYPE)=0

!--- dthm=dt/(2*mass), hmas=mass/2
real(8),allocatable :: dthm(:), hmas(:)

!--- stress tensor
real(8) :: astr(6) = 0.d0

!--- from cmdline_args
logical :: has_initial_pos=.false.
real(8),allocatable,target :: ipos(:,:)

logical :: isSpring=.false.
real(8) :: springConst=0.d0
logical :: hasSpringForce(16)

logical :: isFF=.false., isData=.false., isMDparm=.false., isRunFromXYZ=.false.
character(MAXSTRLENGTH) :: FFPath0="ffield", DataDir0="DAT", ParmPath0="rxmd.in"
character(len=:),allocatable :: FFPath, DataDir, ParmPath, RunFromXYZPath

!<ns> # of atoms to be sent, <nr> # of atoms to be received, <na> # of all of transfered atoms.
!<ne> # of elements for one atom.
!     Example) In case atom type, position and velocity to be sent,  ne = 1+3+3 = 7
integer :: ns, nr, na, ne


!-- variables for timing
integer,parameter :: MAXTIMERS=30
integer :: it_timer(MAXTIMERS)=0
real(8) :: wt0=0.d0

interface forcefield_param_interface
  subroutine forcefield_param(path)
    character(len=:),allocatable,intent(in) :: path
  end subroutine
end interface

interface charge_model_interface
  subroutine charge_model(atype, pos, q)
    real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:)
  end subroutine
end interface

interface force_model_interface
  subroutine force_model(natoms, atype, pos, f, q)
    integer,intent(in out) :: natoms
    real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)
  end subroutine
end interface

interface velocity_scaling_interface
  subroutine velocity_scaling(atype, v, factor)
    real(8),intent(in),allocatable :: atype(:)
    real(8),intent(in out),allocatable :: v(:,:)
    real(8),intent(in),optional :: factor
  end subroutine
end interface

interface mddriver_interface 
  subroutine mddriver(mdbase, num_mdsteps)
    import :: mdbase_class
    type(mdbase_class),intent(in out) :: mdbase
    integer,intent(in) :: num_mdsteps
  end subroutine
end interface

interface mdcontext_interface 
 subroutine get_mdcontext(atype, pos, v, f, q)
   real(8),intent(in out),allocatable,dimension(:) :: atype, q
   real(8),intent(in out),allocatable,dimension(:,:) :: pos, v, f
  end subroutine
end interface

procedure(charge_model),pointer :: charge_model_func => null()
procedure(force_model),pointer :: force_model_func => null()
procedure(forcefield_param),pointer :: get_forcefield_param => null()
procedure(velocity_scaling),pointer :: velocity_scaling_func => null()
procedure(mddriver),pointer :: mddriver_func => null()
procedure(get_mdcontext),pointer :: get_mdcontext_func => null()

end module

!-------------------------------------------------------------------------------------------
module atoms
!-------------------------------------------------------------------------------------------
use utils

integer,parameter :: NMINCELL=4  !<NMINCELL>: Nr of minimum linkedlist cell <-> minimum grain size.
real(8),parameter :: MAXANGLE= 0.999999999999d0 
real(8),parameter :: MINANGLE=-0.999999999999d0
real(8),parameter :: NSMALL = 1.d-10


!--- For array size statistics
!  1-NATOMS, 2-nbrlist, 3-nbrlist for qeq, 4-NBUFFER for move, 5-NBUFFER for copy
!  6-NBUFFER for qeq
integer,parameter :: nmaxas=5
integer,allocatable :: maxas(:,:)

! TE: Total Energy,  KE: Kinetic Energy,  PE :: Potential Energies
!  0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
!  7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
real(8) :: PE(13), GPE(13)

!integer,parameter :: MAXLAYERS_NB=10
!real(8) :: nblcsize(3)
!integer :: nbcc(3), nbnmesh
!integer,allocatable :: nbmesh(:,:)

integer :: nmesh 
integer,allocatable :: mesh(:,:)

!--- QEq variables. 
!<isQEq> flag to run QEq routine: 0-No QEq, 1-CG, 2-Extended Lagrangian
integer :: isQEq
!<NMAXQEq> Number of MAXimum iteration in QEq routine
integer :: NMAXQEq
!<QEq_thrsld> energy criterion in QEq routine
real(8) :: QEq_tol
!<nstep_qeq> counter of iteration
integer :: nstep_qeq, qstep


!--- potential teble
integer,parameter :: NTABLE=5000
real(8),allocatable :: TBL_Eclmb(:,:,:), TBL_Evdw(:,:,:), TBL_Eclmb_QEq(:,:)
real(8) :: UDR, UDRi

logical :: isLG=.false.
logical :: isEfield=.false.

!!<nbllist> Linked List for non-bonding interaction
!!<nbheader> header atom of linkedlist cell for non-bonding interaction
!!<nbnacell> Nr of atoms in a likedlist cell for non-bonding interaction
!integer,allocatable :: nbllist(:), nbheader(:,:,:), nbnacell(:,:,:)
!!<nbplist> neighbor list of nonbonding interaction, non-bonding pair list
!integer,allocatable :: nbplist(:,:)

!<nbrindx> neighbor index
integer,allocatable :: nbrindx(:,:)

!--- coefficient of bonding energy derivative 
real(8),allocatable :: ccbnd(:), cdbnd(:)

!<BO> Bond Order of atoms i-j (nearest neighb only) - (Eq 3a-3d)
real(8),allocatable :: BO(:,:,:) 
real(8),allocatable :: delta(:)

!--- Output variables from the BOp_CALC() subroutine:
real(8),allocatable :: deltap(:,:)
real(8),allocatable :: dln_BOp(:,:,:)
real(8),allocatable :: dBOp(:,:)

!--- A[0123] coefficients for force calculation 
real(8),allocatable :: A0(:,:),A1(:,:), A2(:,:), A3(:,:) 

!--- Passed between Elnpr and E3body
real(8),allocatable :: nlp(:), dDlp(:) !Number of Lone Pairs, its derivatives.
real(8),allocatable :: deltalp(:)

! Two vectors electrostatic energy minimization 
real(8),allocatable,target :: qs(:),qt(:),gs(:), gt(:), hs(:), ht(:), hshs(:), hsht(:)

!--- variables for extended Lagrangian method ---
!<Lex_fqs> fraction between two QEq vectors
!<Lex_w> spring constant
real(8),allocatable,target :: qsfp(:),qsfv(:),qtfp(:),qtfv(:) 
real(8),allocatable :: hessian(:,:)
real(8) :: Lex_fqs=1.0, Lex_w=1.d0, Lex_w2=1.d0, Lex_k=2.d0

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

!------------------------------------------------------------------------------------------
subroutine mdvariables_allocator_reaxff()
use base
use memory_allocator_mod
!------------------------------------------------------------------------------------------
implicit none

!--- extra lists for ReaxFF 
call allocator(nbrindx,1,NBUFFER,1,MAXNEIGHBS)
!call allocator(nbplist,0,MAXNEIGHBS10,1,NBUFFER)

!--- Bond Order Prime and deriv terms:
call allocator(dln_BOp,1,3,1,NBUFFER,1,MAXNEIGHBS)
call allocator(dBOp,1,NBUFFER,1,MAXNEIGHBS)
call allocator(deltap,1,NBUFFER,1,3)
call allocator(deltalp,1,NBUFFER)

!--- Bond Order terms
call allocator(BO,0,3,1,NBUFFER,1,MAXNEIGHBS)
call allocator(delta,1,NBUFFER)
call allocator(A0,1,NBUFFER,1,MAXNEIGHBS)
call allocator(A1,1,NBUFFER,1,MAXNEIGHBS)
call allocator(A2,1,NBUFFER,1,MAXNEIGHBS)
call allocator(A3,1,NBUFFER,1,MAXNEIGHBS)
call allocator(nlp,1,NBUFFER)
call allocator(dDlp,1,NBUFFER)
call allocator(ccbnd,1,NBUFFER)
call allocator(cdbnd,1,NBUFFER)

!--- 2 vector QEq varialbes
call allocator(qs,1,NBUFFER)
call allocator(gs,1,NBUFFER)
call allocator(qt,1,NBUFFER)
call allocator(gt,1,NBUFFER)
call allocator(hs,1,NBUFFER)
call allocator(hshs,1,NBUFFER)
call allocator(ht,1,NBUFFER)
call allocator(hsht,1,NBUFFER)
call allocator(hessian,1,MAXNEIGHBS10,1,NBUFFER)

!--- Varaiable for extended Lagrangian method
call allocator(qtfp,1,NBUFFER)
call allocator(qtfv,1,NBUFFER)

return
end subroutine

end module atoms


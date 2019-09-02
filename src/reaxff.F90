!-------------------------------------------------------------------------------------------
module reaxff_param_mod

use base, only : atmname, mass, NBUFFER, MAXNEIGHBS, it_timer, &
                 current_step, gte, gke, gpe0, ke, pe0, nstep, pstep, astr

use memory_allocator_mod
use utils
use msd_mod, only : msd_data
use stats_mod, only : stats

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

!--- hydrogen bonding interaction cutoff
real(8),parameter :: rchb = 10.d0 ![A]
real(8),parameter :: rchb2 = rchb*rchb

!Coulomb Energy (eq. 22)  
real(8),parameter:: Cclmb0 = 332.0638d0 ! [kcal/mol/A] line 2481 in poten.f
real(8),parameter:: Cclmb0_qeq = 14.4d0 ! [ev]
real(8),parameter:: CEchrge = 23.02d0   ! [ev]
real(8) :: Cclmb = Cclmb0

! cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM()
real(8) :: cutoff_vpar30

character(MAXSTRLENGTH) :: ffFileHeader

!-------------------------------------------------------------------------------------------
!  This module stores the parameters input from the rxmda.in file.  This file is 
!    intentionally very similar to the input setup used by Adri van Duin at Caltech (to aid
!    in cross-checking and updates).  However, please note that there are a few changes made
!    aside from the obvious f77 -> f90 switch.  I have tried to clearly note these. 
!-------------------------------------------------------------------------------------------

!Independant Parameters

integer :: nso    !Number of different types of atoms
integer :: nboty  !Number of different bonds given

! Atom Dependant (ie where they appear in input file - not implementation in code)
real(8),allocatable :: Val(:),Valboc(:)  !Valency of atomtype (norm, boc) 

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

contains

!-------------------------------------------------------------------------------------------
subroutine mddriver_reaxff(mdbase, num_mdsteps)
use base
use atoms
use velocity_modifiers_mod
use communication_mod
use fileio
!-------------------------------------------------------------------------------------------
implicit none
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
integer :: i,ity
real(8) :: ctmp

!if(mdmode==10) call ConjugateGradient(atype,pos)

call charge_model_func(atype, pos, q)
call force_model_func(natoms, atype, pos, f, q)

if(stats%run_stats) call stats%measure(natoms,pos,atype)
stop 'foo'

!--- Enter Main MD loop 
do nstep=0, num_mdsteps-1

   if(mod(nstep,pstep)==0) call print_e_reaxff(atype, v, q)

   if(mod(nstep,fstep)==0) &
        call OUTPUT(atype, pos, v, q, GetFileNameBase(DataDir,current_step+nstep))

   if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:NATOMS,1:3)=vsfact*v(1:NATOMS,1:3)

   if(mod(nstep,sstep)==0.and.mdmode==5) then
      ctmp = (treq*UTEMP0)/( GKE*UTEMP )
      v(1:NATOMS,1:3)=sqrt(ctmp)*v(1:NATOMS,1:3)
   endif

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call gaussian_dist_velocity(atype, v)

!--- element-wise velocity scaling
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call scale_temperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==8) &
      call adjust_temperature(atype, v)

   if(msd_data%is_msd .and. mod(nstep,pstep)==0) &
      call msd_data%measure(NATOMS, atype, pos, ipos)

!--- update velocity
   call vkick(1.d0, atype, v, f) 

!--- update coordinates
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))
   qsfp(1:NATOMS)=qsfp(1:NATOMS)+dt*qsfv(1:NATOMS)

!--- always correct the linear momentum when electric field is applied. 
   if(isEfield) call linear_momentum(atype, v)
   pos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)+dt*v(1:NATOMS,1:3)

!--- migrate atoms after positions are updated
   call COPYATOMS(imode=MODE_MOVE,dr=[0.d0, 0.d0, 0.d0],atype=atype,pos=pos, &
                  v=v,f=f,q=q,ipos=ipos)
   
   if(mod(nstep,qstep)==0) call charge_model_func(atype, pos, q)
   call force_model_func(natoms, atype, pos, f, q)

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

if(msd_data%is_msd) call msd_data%save()

call finalize_reaxff()

end subroutine

!------------------------------------------------------------------------------------------
subroutine get_cutoff_bondorder(rcut, rcut2, maxrcut, natoms_per_type)
!------------------------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in out) :: rcut(:), rcut2(:)
real(8),intent(in out) :: maxrcut
integer(8),intent(in out) :: natoms_per_type(:)

integer :: ity,jty,inxn
real(8) :: dr, BOsig

!--- cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM()
cutoff_vpar30 = cutof2_bo*vpar30

!--- get the cutoff length based on sigma bonding interaction.

! --- Remark --- 
! sigma bond before correction is the longest, namely longer than than pi and double-pi bonds
! thus check only sigma bond convergence.
call allocator(rcut, 1, nboty)
call allocator(rcut2, 1, nboty)

do ity=1, nso
do jty=ity, nso

   inxn=inxn2(ity,jty)
   if(inxn==0) cycle

   dr = 1.0d0
   BOsig = 1.d0
   do while (BOsig > MINBOSIG) 
      dr = dr + 0.01d0
      BOsig = exp( pbo1(inxn)*(dr/r0s(ity,jty))**pbo2(inxn) ) !<- sigma bond prime
   enddo

   rcut(inxn)  = dr
   rcut2(inxn) = dr*dr
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
         if(inxn/=0) rcut(inxn)=0.d0
         inxn=inxn2(jty,ity)
         if(inxn/=0) rcut(inxn)=0.d0
      enddo
   endif
enddo

!--- get the max cutoff length 
maxrcut = maxval(rcut(:))

end subroutine

!------------------------------------------------------------------------------------------
subroutine set_potentialtables_reaxff()
use atoms 
use memory_allocator_mod
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

!--- first element in table 0: potential 1: derivative of potential
call allocator(TBL_EClmb,0,1,1,NTABLE,1,nboty)
call allocator(TBL_Evdw,0,1,1,NTABLE,1,nboty)
call allocator(TBL_EClmb_QEq,1,NTABLE,1,nboty)

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

!-------------------------------------------------------------------------------------------
subroutine get_forcefield_params_reaxff(ffFileName)
use cmdline_args
!-------------------------------------------------------------------------------------------
!  This subroutine is designed solely to obtain the parameters used in the Ecalc.f90 
!    program from the rxmda.in input file.  It is similar to the input routine used
!    by Adri in "reac.f::ffinput"
!-------------------------------------------------------------------------------------------
implicit none

character(len=:),allocatable,intent(in) ::  ffFileName  ! force field parm file

integer :: i,j,inxn   !counters for initialization
integer :: i0,i1,i2,i3,i4,ih  !Counters: # corresp to #-atom depend 

!--- Parameters that count number of entries in each field:
!integer :: nso     !* of atom types (in mod) 
!integer :: nboty  !# of bond terms given (in mod)
integer :: nodmty  !# of off-diag terms given   
integer :: npar, nvaty, ntoty, nhbty

!--- Readin Fields converted to Parameters:
integer :: nodm1, nodm2  !ID bonds to alter in off-diag terms
real(8) :: deodmh,rodmh,godmh ! off-diag term of Evdw parameters
real(8) :: rsig,rpi,rpi2 !temp storage for r0s, etc terms in off-diag part
integer :: typea,typeb   !Temp storage for filling inxn2(:,:) table

!--- NULL Transfer Fields not needed in program (used to calc other values): 
real(8) :: dnull
real(8),allocatable :: vpar(:), bo131(:), bo132(:), bo133(:)
real(8),allocatable :: rvdw1(:), eps(:), alf(:), vop(:)

!--- for LG extension
real(8) :: offdiag_C_lg

dnull = 0.d0
!--- Start Getting Parameters
open(4,file=trim(adjustl(ffFileName)),status="old")
read(4,'(a100)') ffFileHeader

read(4,*) npar  !num of parameters (independ of atom choice)

allocate(vpar(npar))

do i0=1, npar 
   read(4,1300) vpar(i0)  !temp variable...some apparently depend on atype
enddo
  
!--- Constant parameters reset to actual use:
pvdW1 = vpar(29)
pvdW1h = 0.5d0*pvdW1 
pvdW1inv = 1.d0/pvdW1

!--- a small modification in sigma-bond prime <kn>
vpar30 = vpar(30)  

!--- Parameters with 1-atom Depend,  Nr of types of atoms
read(4,'(i3)') nso    

!--- Allocation of variables:
allocate(rat(nso),rapt(nso),vnq(nso))
allocate(r0s(nso,nso),r0p(nso,nso),r0pp(nso,nso))
allocate(Val(nso),Valboc(nso))
allocate(mass(nso))
allocate(bo131(nso),bo132(nso),bo133(nso))
allocate(inxn2(nso,nso),inxn3(nso,nso,nso),inxn3hb(nso,nso,nso),inxn4(nso,nso,nso,nso))
allocate(atmname(nso))
allocate(Vale(nso), plp1(nso), nlpopt(nso), plp2(nso))
allocate(povun2(nso),povun3(nso),povun4(nso))
allocate(povun5(nso),povun6(nso),povun7(nso),povun8(nso))
!--- Valency Terms (j-dependancy only):
allocate(pval3(nso),pval5(nso), Valangle(nso),Valval(nso))
!--- Van der Waals Terms:
allocate(rvdw1(nso), rvdW(nso, nso), eps(nso), Dij(nso,nso))
allocate(alf(nso), alpij(nso,nso))
allocate(vop(nso), gamW(nso,nso))

!--- Coulomb & Charge equilibration:
allocate(chi(nso), eta(nso), gam(nso), gamij(nso,nso))

!--- LG term
if(isLG) then
   allocate(C_lg(nso, nso), Re_lg(nso))
   allocate(rcore2(nso),ecore2(nso),acore2(nso))
   allocate(rcore(nso,nso),ecore(nso,nso),acore(nso,nso))
endif
 
!--- Parameters that still don't depend on atom type yet
plp1(1:nso) = vpar(16)
povun3(1:nso) = vpar(33)
povun4(1:nso) = vpar(32)
povun6(1:nso) = vpar(7)
povun7(1:nso) = vpar(9)
povun8(1:nso) = vpar(10)

!--- skip 3 lines
read(4,*)
read(4,*)
read(4,*)

do i1=1, nso  !collect info on each type of atom
   read(4,1200) atmname(i1),rat(i1),Val(i1),mass(i1),rvdw1(i1),eps(i1),gam(i1),rapt(i1),Vale(i1)
   read(4,1250) alf(i1),vop(i1),Valboc(i1),povun5(i1),dnull,chi(i1),eta(i1),dnull
   read(4,1250) vnq(i1),plp2(i1),dnull,bo131(i1),bo132(i1),bo133(i1),dnull,dnull   

   if (isLG) then 
      read(4,1250) povun2(i1),pval3(i1),dnull,Valval(i1),pval5(i1),rcore2(i1),ecore2(i1),acore2(i1)
      read(4,1250) C_lg(i1,i1), Re_lg(i1)
   else
      read(4,1250) povun2(i1), pval3(i1),dnull,Valval(i1),pval5(i1)
   endif

enddo

!--- update for Mo
do i1=1,nso
   if(mass(i1)<21.d0 .and. Valboc(i1)/=Valval(i1) ) Valboc(i1)=Valval(i1)
enddo

nlpopt(1:nso) = 0.5d0*(Vale(1:nso) - Val(1:nso))
!--- duplicate values
Valangle(1:nso) = Valboc(1:nso)

!--- Calc default r0s, r0p, r0pp:
do i1=1,nso
   do i2=1,nso
!--- Terms for the Bond Order Calculation:
      r0s(i1,i2) = 0.5d0*(rat(i1)+rat(i2))
      r0p(i1,i2) = 0.5d0*(rapt(i1)+rapt(i2))
      r0pp(i1,i2) = 0.5d0*(vnq(i1)+vnq(i2))   
    
!--- Terms used in van der Waals calc: 
      rvdW(i1,i2) = sqrt( 4.d0*rvdw1(i1)*rvdw1(i2) )
      Dij(i1,i2) = sqrt(eps(i1)*eps(i2))
      alpij(i1,i2) = sqrt( alf(i1)*alf(i2) )
      gamW(i1,i2) = sqrt( vop(i1)*vop(i2) )  
      gamij(i1,i2) = ( gam(i1)*gam(i2) )**(-1.5d0) !<- gamcco in reac.f

      if (isLG) then
!--- for LG
        rcore(i1,i2) = sqrt( rcore2(i1)*rcore2(i2) )
        ecore(i1,i2) = sqrt( ecore2(i1)*ecore2(i2) )
        acore(i1,i2) = sqrt( acore2(i1)*acore2(i2) )
      endif

   enddo
enddo  

!--- Parameters with 2-atom Depend:
read(4,1100) nboty  !# of bonds' params given 

!--- Allocation of variables:
allocate(pbo1(nboty),pbo2(nboty),pbo3(nboty),pbo4(nboty),pbo5(nboty),pbo6(nboty),bom(nboty))
allocate(pboc1(nboty),pboc2(nboty),pboc3(nboty),pboc4(nboty),pboc5(nboty))
allocate(desig(nboty), depi(nboty),depipi(nboty),pbe1(nboty),pbe2(nboty)) 
allocate(povun1(nboty),ovc(nboty), v13cor(nboty))

!--- skip one line
read(4,*)

inxn2(1:nso,1:nso) = 0   !allows later flag to tell when a combination is not allowed
ih=0
do i2=1,nboty
    ih=ih+1
    read(4,1400) typea,typeb,Desig(ih),Depi(ih),Depipi(ih),pbe1(ih),pbo5(ih),v13cor(ih),pbo6(ih),povun1(ih)
    read(4,1450) pbe2(ih),pbo3(ih),pbo4(ih),bom(ih),pbo1(ih),pbo2(ih),ovc(ih),dnull  
    inxn2(typea,typeb) = ih
    inxn2(typeb,typea) = ih
enddo


!--- TEMP (required by input file backwards setup)
pboc1(1:nboty) = vpar(1)
pboc2(1:nboty) = vpar(2)

!--- for debugging <kn>
vpar1 = vpar(1)
vpar2 = vpar(2)
 
do i1=1,nso
  do i2=1,nso 
    inxn = inxn2(i1,i2)
    if(inxn/=0) then 
       pboc3(inxn) = sqrt(bo132(i1)*bo132(i2)) ! be careful about variable name
       pboc4(inxn) = sqrt(bo131(i1)*bo131(i2)) ! bo132 -> pboc2, bo131->pbo4  
       pboc5(inxn) = sqrt(bo133(i1)*bo133(i2))
      endif
  enddo
enddo


!--- Changes to off-diagonal terms:
read(4,1100) nodmty  !# of off-diag terms 
do i2=1, nodmty

   if (isLG) then 
     read(4,1400) nodm1,nodm2,deodmh,rodmh,godmh,rsig,rpi,rpi2, offdiag_C_lg
     C_lg(nodm1, nodm2)=offdiag_C_lg
     C_lg(nodm2, nodm1)=offdiag_C_lg
   else 
     read(4,1400) nodm1,nodm2,deodmh,rodmh,godmh,rsig,rpi,rpi2
   endif

   if(rsig.GT.0.d0) r0s(nodm1,nodm2)=rsig
   if(rsig.GT.0.d0) r0s(nodm2,nodm1)=rsig 
   if(rpi.GT.0.d0)  r0p(nodm1,nodm2)=rpi
   if(rpi.GT.0.d0)  r0p(nodm2,nodm1)=rpi
   if(rpi2.GT.0.d0) r0pp(nodm1,nodm2)=rpi2
   if(rpi2.GT.0.d0) r0pp(nodm2,nodm1)=rpi2
   if (rodmh.GT.0.d0) rvdW(nodm1,nodm2)=2.0*rodmh 
   if (rodmh.GT.0.d0) rvdW(nodm2,nodm1)=2.0*rodmh
   if (deodmh.GT.0.d0) Dij(nodm1,nodm2)=deodmh
   if (deodmh.GT.0.d0) Dij(nodm2,nodm1)=deodmh
   if (godmh.GT.0.d0) alpij(nodm1,nodm2)=godmh
   if (godmh.GT.0.d0) alpij(nodm2,nodm1)=godmh
enddo

!!--- Derived 2body parameters 
allocate(cBOp1(nboty), cBOp3(nboty), cBOp5(nboty))
allocate(pbo2h(nboty), pbo4h(nboty), pbo6h(nboty))

!--- <switch> flag to omit pi and double pi bond.
allocate(switch(1:3,nboty))

switch(:,:)=0
do i=1,nso
   do j=1,nso
   inxn = inxn2(i,j)

   if(inxn/=0) then

!!--- In BOp calculation, <switch> will be multiplied to <BOp> to remove
!!--- BOpi and BOpipi for bonding interaction of atoms with a hydrogen.
       if((rat(i)>0.d0)  .and. rat(j)>0.d0 )  switch(1,inxn)=1
       if((rapt(i)>0.d0) .and. rapt(j)>0.d0 ) switch(2,inxn)=1
       if((vnq(i)>0.d0)  .and. vnq(j)>0.d0 )  switch(3,inxn)=1

      if(r0s(i,j)<=0.d0) then
         cBOp1(inxn) = 0.d0
      else
         cBOp1(inxn) = pbo1(inxn)/(r0s(i,j)**pbo2(inxn))
      endif
      if(r0p(i,j)<=0.d0) then
         cBOp3(inxn) = 0.d0
      else
         cBOp3(inxn) = pbo3(inxn)/(r0p(i,j)**pbo4(inxn))
      endif

      if(r0pp(i,j)<=0.d0) then
         cBOp5(inxn) = 0.d0
      else
         cBOp5(inxn) = pbo5(inxn)/(r0pp(i,j)**pbo6(inxn))
      endif

      pbo2h(inxn) = 0.5d0*pbo2(inxn)
      pbo4h(inxn) = 0.5d0*pbo4(inxn)
      pbo6h(inxn) = 0.5d0*pbo6(inxn)
   endif
   enddo
enddo

!--- Input Valency Terms from Input File
inxn3(1:nso,1:nso,1:nso) = 0
read(4,1100) nvaty

allocate(pval1(nvaty),pval2(nvaty),pval4(nvaty),pval6(nvaty))
allocate(pval7(nvaty),pval8(nvaty),pval9(nvaty),pval10(nvaty))
allocate(theta00(nvaty))
allocate(ppen1(nvaty),ppen2(nvaty),ppen3(nvaty),ppen4(nvaty))
allocate(pcoa1(nvaty),pcoa2(nvaty),pcoa3(nvaty),pcoa4(nvaty))

do i=1, nvaty
   read(4,1500) i1,i2,i3,theta00(i),pval1(i),pval2(i),pcoa1(i),pval7(i),ppen1(i),pval4(i)
   inxn3(i1,i2,i3) = i
   inxn3(i3,i2,i1) = i
enddo

!--- Valency Terms which do not depend on inxn type:
pval6(1:nvaty) = vpar(15)
pval8(1:nvaty) = vpar(34)
pval9(1:nvaty) = vpar(17)
pval10(1:nvaty) = vpar(18)
!--- Penalty Terms which do not depend on inxn type:
ppen2(1:nvaty) = vpar(20)
ppen3(1:nvaty) = vpar(21)
ppen4(1:nvaty) = vpar(22)
!--- 3body Conjugation Terms which do not depend on type:
pcoa2(1:nvaty) = vpar(3)
pcoa3(1:nvaty) = vpar(39)
pcoa4(1:nvaty) = vpar(31)
!--- theta00 given in degrees, but used in radians. Convert by:
theta00(1:nvaty) = (pi/180.d0)*theta00(1:nvaty) 


read(4,1100) ntoty
allocate(ptor1(ntoty),ptor2(ntoty),ptor3(ntoty),ptor4(ntoty),V1(ntoty), V2(ntoty),V3(ntoty))
allocate(pcot1(ntoty),pcot2(ntoty))

inxn4(1:nso,1:nso,1:nso,1:nso) = 0  
do i=1,ntoty
   read(4,1600)i1,i2,i3,i4,V1(i),V2(i),V3(i),ptor1(i),pcot1(i),dnull,dnull 
!--- Set up inxn4 lookup reference array
      if(i1==0) then   !condensed input, means that all i1,i4 for this arrangement of i2,i3 are the same
         do i1=1,nso
         do i4=1,nso
         if(inxn4(i1,i2,i3,i4)==0.and.inxn4(i1,i3,i2,i4)==0) then
            inxn4(i1,i2,i3,i4) = i
            inxn4(i4,i2,i3,i1) = i
            inxn4(i1,i3,i2,i4) = i
            inxn4(i4,i3,i2,i1) = i
          endif
          enddo
          enddo
      else
          inxn4(i1,i2,i3,i4) = i
          inxn4(i4,i2,i3,i1) = i
          inxn4(i1,i3,i2,i4) = i
          inxn4(i4,i3,i2,i1) = i 
      endif
enddo

!and a few which don't depend on type 
ptor2(1:ntoty) = vpar(24)
ptor3(1:ntoty) = vpar(25)
ptor4(1:ntoty) = vpar(26)  
pcot2(1:ntoty) = vpar(28)

!--- Input Hydrogen Bond Terms
inxn3hb(1:nso,1:nso,1:nso) = 0
read(4,1100) nhbty
allocate(phb1(nhbty),phb2(nhbty),phb3(nhbty),r0hb(nhbty))

do i=1,nhbty
   read(4,1500) i1,i2,i3,r0hb(i),phb1(i),phb2(i),phb3(i)
   inxn3hb(i1,i2,i3) = i    !Note: inxn3hb(i,j,k) /= inxn3hb(k,j,i)
enddo


!--- close parameter file "ffield"
close(4)

!--- Formats:
1100 format (i3,2x,a2,3x,3d22.15)
1200 format (1x,a2,10f9.4)
1250 format (3x,10f9.4)
1300 format (f10.4)
1400 format (2i3,8f9.4)
1450 format (6x,8f9.4)
1500 format (3i3,7f9.4)
1600 format (4i3,7f9.4)


!--- coefficient of coulomb energy
Cclmb = Cclmb0  !Eclmb

!--- Parameters for charge variable routine. 
!--- In original parameter, chiEEM and etaEEM are given in [ev], not [kcal/mol]
!--- Definition of the stiffness parameter <eta> is different from 
!--- the original code and our code. It's need to be multiplied by 2.
eta(:) = eta(:)*2.d0

if(myid==0) then
   write(6,'(a)') repeat('-',60)
   write(6,'(a20,a30)') 'ffield file: ', trim(adjustl(ffFileName))
   write(6,'(a20,a60)') 'ffield header: ',  trim(adjustl(ffFileHeader))
   do i1=1, nso
      write(6,'(a3,a2,i2,a2, $)') trim(adjustl(atmname(i1))), ' -', i1, ', '
   enddo
   write(6,*)
   write(6,'(a)') repeat('-',60)
  
endif

end subroutine

!-------------------------------------------------------------------------------------------
subroutine print_e_reaxff(atype, v, q)
use mpi_mod
use base, only : hh, hhi, natoms, gnatoms, mdbox, myid, ierr, hmas, wt0
use atoms
use memory_allocator_mod
! calculate the kinetic energy and sum up all of potential energies, then print them.
!-------------------------------------------------------------------------------------------
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
PE0=sum(PE(1:13))

!--- copy data into buffer
buf(0) = PE0; buf(1:13) = PE(1:13); buf(14) = KE
buf(15) = ss; buf(16) = qq;  buf(17:22)=astr(1:6)
call MPI_ALLREDUCE (MPI_IN_PLACE, buf, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!--- copy data from buffer
GPE(1:13) = buf(1:13)/GNATOMS
GPE0 = buf(0)/GNATOMS
GKE = buf(14)/GNATOMS
ss = buf(15); qq = buf(16); astr(1:6)=buf(17:22)

!--- compute properties
tt=GKE*UTEMP
astr(1:6)=astr(1:6)/MDBOX*USTRS/pstep
ss=ss/MDBOX*USTRS/pstep

!--- total energy
GTE = GKE + GPE0
if(myid==0) then
   
   cstep = nstep + current_step 

   write(6,'(a,i9,3es13.5,6es11.3,1x,3f8.2,i4,f8.2,f8.2)') 'MDstep: ',cstep,GTE,GPE0,GKE, &
   GPE(1),sum(GPE(2:4)),sum(GPE(5:7)),sum(GPE(8:9)),GPE(10),sum(GPE(11:13)), &
   tt, ss, qq, nstep_qeq, GetTotalMemory()*1e-9, MPI_WTIME()-wt0

   !write(6,'(a,i9,6f12.6)') 'stress : ',cstep,astr(1:6)

endif

!--- reset stress tensor accumulator
astr(1:6)=0.d0

!--- save current time
wt0 = MPI_WTIME()
end subroutine

!-------------------------------------------------------------------------------------------
subroutine finalize_reaxff()
use mpi_mod
use base, only : myid, ierr, it_timer, MAXTIMERS
use atoms
use memory_allocator_mod
!-------------------------------------------------------------------------------------------
implicit none
integer :: i,it,irt
integer :: ibuf(size(maxas)),ibuf1(size(maxas))
integer :: it_timer_max(size(it_timer)), it_timer_min(size(it_timer))

it_timer_max(:) = 0; it_timer_min(:) = 0

ibuf(:)=0
do i=1,nmaxas
   ibuf(i)=maxval(maxas(:,i))
enddo
call MPI_ALLREDUCE(ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

call MPI_ALLREDUCE(it_timer, it_timer_max, MAXTIMERS, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(it_timer, it_timer_min, MAXTIMERS, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

if(myid==0) then
   call system_clock(it,irt)
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

   print'(a20,f12.4,3x,f12.4)','total (sec): ',dble(it_timer_max(MAXTIMERS))/irt, &
                                               dble(it_timer_min(MAXTIMERS))/irt

   print'(a)','----------------------------------------------'

   print'(a30)', 'rxmd successfully finished'
endif

end subroutine

end module reaxff_param_mod

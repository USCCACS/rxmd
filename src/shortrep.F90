module mod_short_repulsion

use base
use fileio, only : output
use utils, only : l2g, find_cmdline_argc, get_command_argument_str, Ekcal_j, assert, int_to_str

implicit none

include 'mpif.h'

type bond_stats

  real(8) :: dHH_min = 1.35d0, dHH_max= 1.80d0
  real(8) :: dOH_min = 0.92d0, dOH_max= 1.08d0
  real(8) :: dOO_min = 2.20d0, dOO_max= 1d99

  integer :: num_OH_max=0, num_HH_max=0, num_OO_max=0
  integer :: num_OH_min=0, num_HH_min=0, num_OO_min=0

  contains 
    procedure :: set_param => set_param_bond_stats
    procedure :: reset => reset_bond_stats
    procedure :: add => add_bond_stats
    procedure :: print => print_bond_stats

end type

type short_repulsion_type1_params

  !real(8),parameter :: alpha=0.1d0, GeGe_x=10d0, SeSe_x=1.11d0
  !real(8),parameter :: A_GeSe = 3.5d0-20, Ge_sig=0.73d0, Se_sig=2.d0 ! [J] & [A]
  !integer,parameter :: GeGe_eta = 11, GeSe_eta = 9, SeGe_eta=9, SeSe_eta = 7

  integer :: num_types
  real(8) :: alpha
  real(8),allocatable :: sigma(:)
  integer,allocatable :: eta(:,:)
  real(8),allocatable :: multiplier(:,:)
  real(8),allocatable :: Aij(:,:)
  real(8),allocatable :: coef(:,:)

  logical :: is_initialized = .false.

  contains
    procedure :: read => read_type1_params
    procedure :: print => print_type1_params
    procedure :: apply => apply_type1_short_repulsion
end type

type short_repulsion_type2_params

  real(8) :: alpha_bond, alpha_angle, vscale
  real(8) :: rc_nbr

  real(8) :: Kr, Kq, Kr_O
  real(8) :: r0, q0, r0_O

  real(8) :: rc_freeze

  real(8) :: beta_1, beta_2
  real(8) :: beta_s1, beta_l1, beta_s2, beta_l2

  integer :: htype=-1, otype=-1, getype=-1, setype=-1

  real(8) :: fcut_o=75d0, fcut_h=50d0, ffactor=0.7d0

  real(8) :: rc_inner, rc_outer, rc_hh_min, rc_oo_min

  real(8) :: stop_OH_min, stop_OH_max

  real(8),allocatable :: f_spring(:,:)

  character(2),allocatable :: atom_name(:)

  logical :: is_initialized = .false.

  logical :: does_freeze_atoms = .false.
  logical :: does_flip_velocity = .false.

  logical,allocatable :: frozen_atom_flag(:)
  real(8),allocatable :: r_oh(:,:)

  type(bond_stats) :: bstat
  real(8) :: dHH_min, dHH_max, dOH_min, dOH_max, dOO_min, dOO_max ! bond_stats parameters

  contains 
    procedure :: read => read_type2_params
    procedure :: print => print_type2_params
    procedure :: apply2 => apply_type2_short_repulsion
    procedure :: apply4 => apply_type4_short_repulsion
    procedure :: apply5 => apply_type5_short_repulsion
    procedure :: flipv => flip_frozenh_velocity
    procedure :: freezex => freeze_h_in_overstreched_bond
    procedure :: oo_spring => spring_potential_for_OObond
end type

type short_repulsion_type

  logical :: has_short_repulsion=.false.

  integer :: potential_type = 0 

  type(short_repulsion_type1_params) :: p1
  type(short_repulsion_type2_params) :: p2

end type 

type lj_potential_type
  real(8) :: sigma=2.5d0, epsiron=1.d0
  real(8) :: rcutoff=2.6d0

  contains 
    procedure :: init => set_cutoff_lj_potential
    procedure :: calc => calc_force_lj_potential
    procedure :: save_table => save_table_lj_potential
end type


type(short_repulsion_type) short_rep

type(bond_stats) bstat

type(lj_potential_type) lj_pot

contains

!-----------------------------------------------------------------------------
subroutine set_param_bond_stats(this, dHH_min, dHH_max, dOH_min, dOH_max, dOO_min, dOO_max)
!-----------------------------------------------------------------------------
   class(bond_stats),intent(in out) :: this
   real(8),intent(in) :: dHH_min, dHH_max, dOH_min, dOH_max, dOO_min, dOO_max

   this%dHH_min = dHH_min 
   this%dHH_max = dHH_max 
   this%dOH_min = dOH_min 
   this%dOH_max = dOH_max 
   this%dOO_min = dOO_min
   this%dOO_max = dOO_max
    
end subroutine

!-----------------------------------------------------------------------------
subroutine spring_potential_for_OObond(this, dr, U, Up) 
!-----------------------------------------------------------------------------
   class(short_repulsion_type2_params),intent(in) :: this
   real(8),intent(in) :: dr

   !  force,enegry,force at cutoff, energy at cutoff
   !real(8),parameter :: Kr=35.37d0, r0=2.7d0 ! changed October 7, 2:42pm
   !real(8),parameter :: Kr=17.685d0, r0=2.7d0  ! changed October 7, 8:00pm
   !real(8),parameter :: Kr=35.37d0, r0=2.6d0
   real(8) :: Kr, r0
   real(8),intent(in out) ::  U, Up
   real(8) :: U_c, Up_c 

   real(8) :: sigr, sigr12, sigr6

   Kr = this%Kr_O
   r0 = this%r0_O

   U = 0.5d0*Kr*(dr-r0)**2
   Up = Kr*(dr-r0)

   !if(dr>this%rcutoff) then
   if(dr>r0) then
     Up=0.d0; U=0.d0
     return 
   endif

   return
end subroutine

!-----------------------------------------------------------------------------
subroutine set_cutoff_lj_potential(this, rcutoff) 
!-----------------------------------------------------------------------------
   class(lj_potential_type),intent(in out) :: this
   real(8),intent(in)  :: rcutoff

   ! TODO the cutoff distance should be set from outside
   this%rcutoff = rcutoff  

end subroutine

!-----------------------------------------------------------------------------
subroutine calc_force_lj_potential(this, dr, U, Up) 
!-----------------------------------------------------------------------------
   class(lj_potential_type),intent(in) :: this
   real(8),intent(in) :: dr

   !  force,enegry,force at cutoff, energy at cutoff
   real(8),intent(in out) ::  U, Up
   real(8) :: U_c, Up_c 

   real(8) :: sigr, sigr12, sigr6


   if(dr>this%rcutoff) then
     Up=0.d0; U=0.d0
     return 
   endif


   sigr = this%sigma/this%rcutoff
   sigr6 = sigr**6
   sigr12 = sigr6*sigr6

   U_c = 4.d0*this%epsiron*(sigr12-sigr6) ! energy at cutoff
   Up_c = -48d0*this%epsiron*(sigr12-0.5d0*sigr6)/this%rcutoff ! force at cutoff


   sigr = this%sigma/dr
   sigr6 = sigr**6
   sigr12 = sigr6*sigr6

   U = 4.d0*this%epsiron*(sigr12-sigr6)-U_c - (dr-this%rcutoff)*Up_c
   Up = -48d0*this%epsiron*(sigr12-0.5d0*sigr6)/dr - Up_c

   return
end subroutine

!-----------------------------------------------------------------------------
subroutine save_table_lj_potential(this, myrank)
!-----------------------------------------------------------------------------
   class(lj_potential_type),intent(in) :: this
   integer,intent(in) :: myrank

   integer,parameter :: ntables = 100

   integer :: i,iunit
   real(8) :: dr, U, Up

   if(myrank==0) then
      open(newunit=iunit,file='lj_pot.dat')
      do i=1, ntables
         dr = dble(i)*this%rcutoff/ntables
         call this%calc(dr, U, Up)
         write(iunit,fmt='(3es15.5)') dr, U, Up
      enddo
      close(iunit)
   endif

end subroutine
   

!-----------------------------------------------------------------------------
subroutine print_bond_stats(this, myrank)
!-----------------------------------------------------------------------------
  class(bond_stats),intent(in out) :: this
  integer,intent(in) :: myrank

  call this%add()

  if(myrank==0) then
     print'(a,3i18)','bstat: num_OH_min,num_HH_min,num_OO_min :', &
        this%num_OH_min, this%num_HH_min, this%num_OO_min
     print'(a,3i18)','bstat: num_OH_max,num_HH_max,num_OO_max :', &
        this%num_OH_max, this%num_HH_max, this%num_OO_max
  endif

  call this%reset()

end subroutine

!-----------------------------------------------------------------------------
subroutine reset_bond_stats(this)
!-----------------------------------------------------------------------------
  class(bond_stats),intent(in out) :: this

  this%num_oh_max=0; this%num_hh_max=0; this%num_oo_max=0
  this%num_oh_min=0; this%num_hh_min=0; this%num_oo_min=0

end subroutine

!-----------------------------------------------------------------------------
subroutine add_bond_stats(this)
!-----------------------------------------------------------------------------
  class(bond_stats),intent(in out) :: this
  integer :: array(6)

  array(1)=this%num_hh_max
  array(2)=this%num_oh_max
  array(3)=this%num_oo_max

  array(4)=this%num_hh_min
  array(5)=this%num_oh_min
  array(6)=this%num_oo_min

  call MPI_ALLREDUCE(MPI_IN_PLACE, array, size(array), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  this%num_hh_max=array(1)
  this%num_oh_max=array(2)
  this%num_oo_max=array(3)

  this%num_hh_min=array(4)
  this%num_oh_min=array(5)
  this%num_oo_min=array(6)

end subroutine

!-----------------------------------------------------------------------------
subroutine initialize_short_repulsion(sr, atom_name)
!-----------------------------------------------------------------------------
type(short_repulsion_type),intent(in out) :: sr
character(2),allocatable,intent(in) :: atom_name(:)

integer :: idx, ity, jty, funit
real(8) :: sigma_sum

character(len=:),allocatable :: filename
logical :: fexists

if(find_cmdline_argc('--short_rep',idx)) then
   filename = 'shortrep.in'
   inquire(file=filename, exist=sr%has_short_repulsion)

   if(.not. sr%has_short_repulsion) then

     filename = get_command_argument_str(idx+1)
     inquire(file=filename, exist=sr%has_short_repulsion)

     if(.not. sr%has_short_repulsion) then
       if(myid==0) then
          print'(a)', 'ERROR: short_rep flag is given but could not find shortrep.in file. Exiting rxmd.'
          call MPI_FINALIZE(ierr)
          stop
       endif
     endif

   endif
endif

! the short_repulsion is not applied if input file is not found by here.
if(.not. sr%has_short_repulsion) return

open(newunit=funit, file=filename, form='formatted')

read(funit,*) sr%potential_type

if(myid==0) then
  print'(a)', repeat('-',60)
  print'(a)', 'parameters for short repulsion: '
  print'(a,i3)', 'potential type: ', sr%potential_type
  print'(a)', repeat('-',60)
endif

if(sr%potential_type==1) then
   call sr%p1%read(funit)
   call sr%p1%print()
else if(sr%potential_type==2.or.sr%potential_type==3.or.sr%potential_type==4.or.sr%potential_type==5) then
   call sr%p2%read(funit, atom_name, sr%potential_type) 
   call sr%p2%print(sr%potential_type)
else
   print'(a)', 'ERROR: unsupported potential type. Exiting the code', sr%potential_type
   stop
endif

if(myid==0) print'(a)', repeat('-',60)

close(funit)

end subroutine

!-----------------------------------------------------------------------------
subroutine read_type2_params(this, funit, atom_name, potential_type)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type2_params) :: this
  character(2),allocatable,intent(in) :: atom_name(:)
  integer,intent(in) :: potential_type

  integer,intent(in) :: funit 
  integer :: ity, jty
  real(8) :: sigma_sum

  call assert(2<=potential_type .and. potential_type<=5, &
     'unsupported potential type in read_type2_params(). Exiting rxmd.', myid)

  allocate(this%f_spring(NBUFFER,3))

  do ity=1, size(atom_name)
     if(index(atom_name(ity),"H") /=0) this%htype=ity
     if(index(atom_name(ity),"O") /=0) this%otype=ity
     if(index(atom_name(ity),"Ge") /=0) this%getype=ity
     if(index(atom_name(ity),"Se") /=0) this%setype=ity
  enddo

  read(funit,*) this%alpha_bond, this%alpha_angle, this%rc_nbr
  read(funit,*) this%Kr, this%Kq
  read(funit,*) this%r0, this%q0
  read(funit,*) this%Kr_O, this%r0_O

  read(funit,*) this%beta_1, this%beta_s1, this%beta_l1
  read(funit,*) this%beta_2, this%beta_s2, this%beta_l2
  read(funit,*) this%fcut_o, this%fcut_h, this%ffactor
  read(funit,*) this%dHH_min, this%dOH_max, this%dOO_min
  read(funit,*) this%dHH_max, this%dOH_min, this%dOO_max 
  read(funit,*) this%rc_inner, this%rc_outer
  read(funit,*) this%rc_hh_min, this%rc_oo_min
  read(funit,*) this%stop_OH_min, this%stop_OH_max

  this%is_initialized = .true.

end subroutine

!-----------------------------------------------------------------------------
subroutine print_type2_params(this, potential_type)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type2_params) :: this
  integer,intent(in) :: potential_type
  integer :: iunit

  if(myid==0) then
    if(this%htype>0) print'(a30,a7,i3)',  "atom type: ", "H -", this%htype
    if(this%otype>0) print'(a30,a7,i3)',  "atom type: ", "O -", this%otype
    if(this%getype>0) print'(a30,a7,i3)', "atom type: ", "Ge -", this%getype
    if(this%setype>0) print'(a30,a7,i3)', "atom type: ", "Se -", this%setype
    print'(a35,3f10.3)', 'alpha_bond, alpha_angle, rc_nbr: ', this%alpha_bond, this%alpha_angle, this%rc_nbr
    print'(a30,2f10.3)', 'Kr, Kq: ', this%Kr, this%Kq
    print'(a30,2f10.3)', 'r0, q0: ', this%r0, this%q0
    print'(a30,2f10.3)','Kr_O, r0_O: ', this%Kr_O, this%r0_O

    print'(a30,3f10.3)','beta_1, beta_s1, beta_l1: ', this%beta_1, this%beta_s1, this%beta_l1
    print'(a30,3f10.3)','beta_2, beta_s2, beta_l2: ', this%beta_2, this%beta_s2, this%beta_l2

    print'(a30,3f10.3)','fcut_o, fcut_h, ffactor: ', this%fcut_o, this%fcut_h, this%ffactor
    print'(a30,3f10.3)','dHH_min, dOH_min, dOO_min: ',this%dHH_min, this%dOH_max, this%dOO_min
    print'(a30,2f10.3,es10.1)','dHH_max, dOH_max, dOO_max: ',this%dHH_max, this%dOH_min, this%dOO_max
    print'(a30,4f10.3)','rc_inner, rc_outer: ', this%rc_inner, this%rc_outer
    print'(a30,4f10.3)','rc_hh_min, rc_oo_min: ', this%rc_hh_min, this%rc_oo_min
    print'(a30,2f10.3)','stop_OH_min, stop_OH_max: ',this%stop_OH_min, this%stop_OH_max

    open(newunit=iunit, file='shortrep.in.current',form='formatted')
    write(iunit,*) potential_type
    write(iunit,*) this%alpha_bond, this%alpha_angle, this%rc_nbr
    write(iunit,*) this%Kr, this%Kq
    write(iunit,*) this%r0, this%q0
    write(iunit,*) this%Kr_O, this%r0_O
    
    write(iunit,*) this%beta_1, this%beta_s1, this%beta_l1
    write(iunit,*) this%beta_2, this%beta_s2, this%beta_l2
    write(iunit,*) this%fcut_o, this%fcut_h, this%ffactor
    write(iunit,*) this%dHH_min, this%dOH_min, this%dOO_min
    write(iunit,*) this%dHH_max, this%dOH_max, this%dOO_max 
    write(iunit,*) this%rc_inner, this%rc_outer
    write(iunit,*) this%rc_hh_min, this%rc_oo_min
    write(iunit,*) this%stop_OH_min, this%stop_OH_max
    close(iunit)

  endif

end subroutine

!-----------------------------------------------------------------------------
subroutine read_type1_params(this, funit)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type1_params) :: this
  integer,intent(in) :: funit 
  integer :: ity, jty
  real(8) :: sigma_sum

  read(funit,*) this%num_types
  
  allocate(this%Aij(this%num_types,this%num_types), &
           this%multiplier(this%num_types,this%num_types), &
           this%coef(this%num_types,this%num_types))
  allocate(this%eta(this%num_types,this%num_types), this%sigma(this%num_types))
  
  read(funit,*) this%alpha
  read(funit,*) (this%sigma(ity), ity=1, this%num_types)
  do ity = 1, this%num_types
     read(funit,*) (this%eta(ity,jty), jty=1, this%num_types)
  enddo
  do ity = 1, this%num_types
     read(funit,*) (this%Aij(ity,jty), jty=1, this%num_types)
  enddo
  do ity = 1, this%num_types
     read(funit,*) (this%multiplier(ity,jty), jty=1, this%num_types)
  enddo
  
  do ity=1,this%num_types
  do jty=1,this%num_types
     sigma_sum = this%sigma(ity)+this%sigma(jty)
     this%coef(ity,jty) = this%alpha*(this%Aij(ity,jty)/Ekcal_j)*sigma_sum**this%eta(ity,jty)
  
     this%coef(ity,jty) = this%multiplier(ity,jty)*this%coef(ity,jty)
  enddo; enddo

  this%is_initialized = .true.
end subroutine

!-----------------------------------------------------------------------------
subroutine print_type1_params(this)
!-----------------------------------------------------------------------------
implicit none
class (short_repulsion_type1_params),intent(in) :: this
integer :: ity,jty

if(myid==0) then
   print'(a,i9)', 'num_types : ', this%num_types

   print*

   print'(9(a4,i1,a1,1x,f10.3,3x))', ('sig(',ity,')',this%sigma(ity), ity=1, this%num_types)
   do ity=1, this%num_types
      print'(9(a4,i1,a1,i1,a1,1x,i6,3x))', ('eta(',ity,',',jty,')', this%eta(ity,jty), jty=1, this%num_types)
   enddo
   do ity=1, this%num_types
      print'(9(a4,i1,a1,i1,a1,1x,es13.2,3x))', ('Aij(',ity,',',jty,')', this%Aij(ity,jty), jty=1, this%num_types)
   enddo

   print*

   print'(a,f10.3)', 'alpha : ', this%alpha
   do ity=1, this%num_types
      print'(9(a4,i1,a1,i1,a1,1x,f8.3,2x))', ('mul(',ity,',',jty,')', this%multiplier(ity,jty), jty=1, this%num_types)
   enddo

   print*

   do ity=1, this%num_types
      print'(9(a5,i1,a1,i1,a1,1x,f10.3,3x))', ('coef(',ity,',',jty,')',this%coef(ity,jty), jty=1, this%num_types)
   enddo
   print'(a)', repeat('-',60)
endif

end subroutine

!-----------------------------------------------------------------------
subroutine ForceA3(coeff,i,j,k,da0, da1, f)
! derivative of <cos_ijk>
!-----------------------------------------------------------------------
implicit none
! Caa(0,0)=Caa, Caa(0,-1)=Caa-1, Caa(-1,0)=Ca-1a, Caa(-1,-1)=Ca-1a-1
real(8),intent(in out) :: f(NBUFFER,3)
real(8),INTENT(IN) :: coeff, da0(0:3), da1(0:3)
integer,INTENT(IN) :: i,j,k
real(8) :: Caa(-2:0,-2:0), Ci(3), Ck(3)
real(8) :: fij(3), fjk(3), fijjk(3), rij(3), rjk(3)
real(8) :: CCisqr, coCC

Caa( 0,0) = da0(0)**2; Caa( 0,-1) = sum(da0(1:3)*da1(1:3))
Caa(-1,0) = Caa(0,-1); Caa(-1,-1) = da1(0)**2

CCisqr = 1.d0/( da0(0)*da1(0) )
coCC = coeff*CCisqr

!--- Some of calculations are unnecessary due to the action-reaction relation.
Ci(1) = -( Caa(0,-1)/Caa(0,0) )
Ci(2) =  1.d0
!Cj(1) =  Caa( 0,-1)/Caa( 0, 0) + 1.d0 
!Cj(2) = -( Caa(0,-1)/Caa(-1,-1) + 1.d0 ) 
Ck(1) = -1.d0
Ck(2) =  Caa(0,-1)/Caa(-1,-1)

rij(1:3) = da0(1:3)
rjk(1:3) = da1(1:3)

fij(1:3) = coCC*(Ci(1)*rij(1:3) + Ci(2)*rjk(1:3))
fjk(1:3) =-coCC*(Ck(1)*rij(1:3) + Ck(2)*rjk(1:3))
fijjk(1:3) =  -fij(1:3) + fjk(1:3) 

!$omp atomic
f(i,1) = f(i,1) + fij(1)
!$omp atomic
f(i,2) = f(i,2) + fij(2)
!$omp atomic
f(i,3) = f(i,3) + fij(3)

!$omp atomic
f(j,1) = f(j,1) + fijjk(1)
!$omp atomic
f(j,2) = f(j,2) + fijjk(2)
!$omp atomic
f(j,3) = f(j,3) + fijjk(3)

!$omp atomic
f(k,1) = f(k,1) - fjk(1)
!$omp atomic
f(k,2) = f(k,2) - fjk(2)
!$omp atomic
f(k,3) = f(k,3) - fjk(3)

!--- Check N3rd ---
!print'(a,6f20.13)','N3rd: ', &
!Ci(1)-Ci(2), -Cj(1), Ci(2), -Ck(1), Cj(2), Ck(1)-Ck(2)

end subroutine

!-----------------------------------------------------------------------------
subroutine apply_type5_short_repulsion(this, potential_type, bstat)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this
integer,intent(in) :: potential_type
type(bond_stats),intent(in out) :: bstat

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, j1, ity, jty, k, k1, k2, k3, k4, igid, jgid, kgid

integer,parameter :: MAXNEIGHBS = 100
integer :: nbr_count, nbr_index(MAXNEIGHBS)
integer :: se_list(4)
real(8) :: nbr_dist(MAXNEIGHBS), dist_1st, dist_2nd
real(8) :: rsese(0:3), rij(0:3), rik(0:3), rij0(0:3)
real(8) :: rhh(0:3), roo(0:3), roo_min
real(8) :: sine, cosine, theta, Krcoef, Kqcoef, ff(3)

real(8),parameter :: MAX_DIST = 1e9, pi=3.14159265358979d0

!call set_force_zero(this%rc_inner, this%rc_outer, this%rc_hh_min)
this%f_spring=0.d0

call bstat%reset()

do i=1, NATOMS

   ity = nint(atype(i))

   if(ity/=this%getype) cycle ! only Ge

   igid = l2g(atype(i))

   roo_min = 1d99
   nbr_count = 0; nbr_dist = MAX_DIST
   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jty = nint(atype(j))

      ! Ge-Ge force cut
      if(jty==this%getype) then

        rsese(1:3) = pos(i,1:3) - pos(j,1:3)
        rsese(0) = sqrt(sum(rsese(1:3)*rsese(1:3)))
        if(roo_min>rsese(0)) roo_min = rsese(0)

        if(rsese(0)<this%rc_oo_min) then  
           f(i,1:3)=0.d0 ! clear ML force for Ge-Ge
           f(j,1:3)=0.d0 ! clear ML force for Ge-Ge
        endif
      endif

      if(jty/=this%setype) cycle  ! only Se

      dr(1:3) = pos(i,1:3) - pos(j,1:3)
      dr2 = sum(dr(1:3)*dr(1:3))
      dr1 = sqrt(dr2)

      if(dr1>this%rc_nbr) cycle ! Ge-Se bond cutoff for neighbor list

      nbr_count = nbr_count + 1
      nbr_index(nbr_count) = j
      nbr_dist(nbr_count) = dr1

      !if(myid==0) print'(a,i6,2i6,f8.3)','jty, this%htype:',i,jty,this%htype, dr1

      !f(i,1:3) = f(i,1:3) + fcoef*dr(1:3)
   enddo

   call assert(nbr_count>=4, 'ERROR: found less than 4 Se neighbors from Ge ', myid)
   !print*,'nbr_index: ', i,nbr_count,nbr_index(1:nbr_count)
   !print*,'nbr_dist: ', i,nbr_count,nbr_dist(1:nbr_count)

   if(roo_min<bstat%dOO_min) bstat%num_oo_min = bstat%num_oo_min + 1

   ! pickup 4 Se atoms
   se_list = -1
   do j1 = 1, 4
      j = minloc(nbr_dist,dim=1)
      se_list(j1) = nbr_index(j)
      nbr_dist(j) = MAX_DIST
   enddo
   !print'(a,6i6)','myid,se_list: ', myid,i,se_list(1:4)

   do j1 = 1, size(se_list)

      j = se_list(j1)
      jgid = l2g(atype(j))

      rij(1:3) = pos(i,1:3)-pos(j,1:3)
      rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
      rij(1:3) = rij(1:3)/rij(0)

      !*** clear ML force for Ge-Se ***
      if(rij(0)<this%rc_inner.or.this%rc_outer<rij(0)) then
        !print'(a,f6.3,7e11.3)','fzero: myid,i,j,ity,jty,rij,fi,fj: '//info_ij,rij(0),f(i,1:3),f(j,1:3)
        f(j,1:3)=0.d0  
      endif

      Krcoef = this%alpha_bond * this%Kr*(rij(0)-this%r0) ! ij

      ff(1:3) = Krcoef*rij(1:3)
      this%f_spring(i,1:3) = this%f_spring(i,1:3) - ff(1:3)
      this%f_spring(j,1:3) = this%f_spring(j,1:3) + ff(1:3)

      if(rij(0)>bstat%dOH_max) print'(a,2i9,f8.3)','OUTLIER: j-atom is too far',igid,jgid,rij(0)
      if(rij(0)<bstat%dOH_min) print'(a,2i9,f8.3)','OUTLIER: j-atom is too close ',igid,jgid,rij(0)
      if(rij(0)>bstat%dOH_max) bstat%num_oh_max = bstat%num_oh_max + 1
      if(rij(0)<bstat%dOH_min) bstat%num_oh_min = bstat%num_oh_min + 1

      do k1 = j1 + 1, size(se_list)

         k = se_list(k1)
         kgid = l2g(atype(k))

         rik(1:3) = pos(i,1:3)-pos(k,1:3)
         rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
         rik(1:3) = rik(1:3)/rik(0)

         !*** clear ML force for Ge-Se ***
         if(rik(0)<this%rc_inner.or.this%rc_outer<rik(0)) then
           !print'(a,f6.3,7e11.3)','fzero: myid,i,k,ity,kty,rik,fi,fk: '//info_ik,rik(0),f(i,1:3),f(k,1:3)
           f(k,1:3)=0.d0  
         endif

         ff(1:3) = Krcoef*rik(1:3)
         this%f_spring(i,1:3) = this%f_spring(i,1:3) - ff(1:3)
         this%f_spring(k,1:3) = this%f_spring(k,1:3) + ff(1:3)

         cosine = sum(rij(1:3)*rik(1:3))
         theta = acos(cosine)
         sine = sin(theta)

         call assert(abs(sine)>1d-9, 'ERROR: too small Se-Ge-Se angle', myid)

         Kqcoef = this%alpha_angle * this%Kq*(theta*180d0/pi-this%q0)*(-1.d0/sine)

         rij0(0)=rij(0); rij0(1:3)=-rij(1:3)
         call ForceA3(Kqcoef,j,i,k,rij0,rik, this%f_spring)

         ! Se-Se distance 
         rhh(1:3) = pos(j,1:3) - pos(k,1:3)
         rhh(0) = sqrt(sum(rhh(1:3)*rhh(1:3)))

         !*** clear ML force for Se-Se ***
         if(rhh(0) < this%rc_hh_min) then
           f(j,1:3)=0d0 
           f(k,1:3)=0d0 
         endif

         if(rik(0)>bstat%dOH_max) print'(a,3i9,3f8.3)','OUTLIER: k-atom is too far ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
         if(rik(0)<bstat%dOH_min) print'(a,3i9,3f8.3)','OUTLIER: k-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
         if(rik(0)>bstat%dOH_max) bstat%num_oh_max = bstat%num_oh_max + 1
         if(rik(0)<bstat%dOH_min) bstat%num_oh_min = bstat%num_oh_min + 1

         if(rhh(0)>bstat%dHH_max) bstat%num_hh_max = bstat%num_hh_max + 1
         if(rhh(0)<bstat%dHH_min) bstat%num_hh_min = bstat%num_hh_min + 1

      enddo

      !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)
   enddo

enddo

do i=1, size(f,dim=1)
do j=1, size(f,dim=2)
   f(i,j)=f(i,j)+this%f_spring(i,j)
enddo; enddo

call bstat%print(myid)

end subroutine

!-----------------------------------------------------------------------------
subroutine apply_type4_short_repulsion(this, potential_type, bstat, lj_pot)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in) :: this
integer,intent(in) :: potential_type
type(bond_stats) :: bstat
type(lj_potential_type) :: lj_pot

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, k, l, l1, ity, jty, kty, lty, igid, jgid, kgid

real(8) :: rij(0:3), rik(0:3)
real(8) :: sine, cosine, theta, Krcoef, Kqcoef, ff(3)

real(8),parameter :: MAX_DIST = 1e9, pi=3.14159265358979d0

logical :: is_found_k, is_found_j
real(8) :: dr_k, dr_j

real(8) :: bond_coeff

character(len=:),allocatable :: info_ij, info_ik
character(len=:),allocatable :: filebase

real(8) :: rhh(0:3), roo(0:3), roo_min
real(8) :: U, Up

integer,save :: step=0
character(9) :: a9

!real(8),parameter :: rc_inner=0.89d0, rc_outer=1.11d0  ! before October 06, 9:29pm
!real(8),parameter :: rc_inner=0.93d0, rc_outer=1.06d0   ! changed October 06, 9:29pm
real(8),parameter :: rc_inner=0.92d0, rc_outer=1.07d0   ! changed October 08, 9:04pm
real(8),parameter :: rc_hh_min = 1.44d0, rc_oo_min=2.5d0
real(8),parameter :: stop_OH_min = 0.5d0, stop_OH_max = 1.5d0

call set_force_zero(this%rc_inner, this%rc_outer, this%rc_hh_min)

call bstat%reset()

!O-O NN force screening
do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   do l1 = 1, nbrlist(i,0)

      l = nbrlist(i,l1)
      lty = nint(atype(l))

      if(lty==this%otype) then

        roo(1:3) = pos(i,1:3) - pos(l,1:3)
        roo(0) = sqrt(sum(roo(1:3)*roo(1:3)))

        !if(roo(0)<2.5d0) then
        if(roo(0)<this%rc_oo_min) then
           f(i,1:3)=0.d0
           f(l,1:3)=0.d0
        endif
      endif
   enddo
enddo


do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   !if(myid==0) print'(a,2i6)','ity, this%otype:',ity,this%otype

   igid = l2g(atype(i))

   is_found_j=.false.; jgid=-1; dr_j=0d0; jty=-1
   is_found_k=.false.; kgid=-1; dr_k=0d0; kty=-1

   roo_min = 1d99
   do l1 = 1, nbrlist(i,0)

      l = nbrlist(i,l1)
      lty = nint(atype(l))

      if(lty==this%otype) then
        ! find minimum O-O distance
        roo(1:3) = pos(i,1:3) - pos(l,1:3)
        roo(0) = sqrt(sum(roo(1:3)*roo(1:3)))
        if(roo_min>roo(0)) roo_min = roo(0)

        ! apply LJ potential 
        call this%oo_spring(roo(0), U, Up)
        !if(myid==0) print'(a,i2,2es15.5)','lj_pot: ', i,l,roo(0),Up
        f(i,1:3) = f(i,1:3) - Up*roo(1:3)/roo(0)
        f(l,1:3) = f(l,1:3) + Up*roo(1:3)/roo(0)
      endif

      if(lty/=this%htype) cycle

      if(l2g(atype(l))==igid+1) then
        j = l   !H1
        is_found_j = .true.
        jgid = l2g(atype(j))
        jty = nint(atype(j))
      endif
      if(l2g(atype(l))==igid+2) then
        k = l   !H2
        is_found_k = .true.
        kgid = l2g(atype(k))
        kty = nint(atype(k))
      endif

   enddo

   info_ij = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(jgid)//' '//int_to_str(ity)//' '//int_to_str(jty)
   info_ik = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(kgid)//' '//int_to_str(ity)//' '//int_to_str(kty)

   call assert(is_found_j, 'ERROR: wrong global ID for j '//info_ij)
   call assert(is_found_k, 'ERROR: wrong global ID for k '//info_ik)

   call assert(jty==this%htype, 'ERROR: j-atom type is not H '//info_ij)
   call assert(kty==this%htype, 'ERROR: k-atom type is not H '//info_ik)

   rij(1:3) = pos(i,1:3)-pos(j,1:3)
   rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
   rij(1:3) = rij(1:3)/rij(0)

   rik(1:3) = pos(i,1:3)-pos(k,1:3)
   rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
   rik(1:3) = rik(1:3)/rik(0)

   cosine = sum(rij(1:3)*rik(1:3))
   theta = acos(cosine)
   sine = sin(theta)

   call assert(abs(sine)>1d-9, 'ERROR: too small H-O-H angle ', myid)

   bond_coeff = this%alpha_bond
   if(potential_type==4) then
     if(rij(0)<this%beta_s1.or.this%beta_l1<rij(0)) bond_coeff = bond_coeff + this%beta_1
     if(rij(0)<this%beta_s2.or.this%beta_l2<rij(0)) bond_coeff = bond_coeff + this%beta_2
   endif

   ! H-H distance 
   rhh(1:3) = pos(j,1:3) - pos(k,1:3)
   rhh(0) = sqrt(sum(rhh(1:3)*rhh(1:3)))

   !Krcoef = bond_coeff * this%Kr*(rij(0)-this%r0) ! ij
   Krcoef = bond_coeff * get_smooth_spring_coeff(rij(0),this%Kr,this%r0) ! ij
   !if(potential_type==4) then
   !  if(rij(0)<this%beta_s1.or.this%beta_l1<rij(0)) &
   !    print'(a,f8.3,f8.3,2es15.5)', 'beta1: myid,i,j,ity,jty: '//info_ij,bond_coeff,rij(0),Krcoef,rhh(0)
   !  if(rij(0)<this%beta_s2.or.this%beta_l2<rij(0)) &
   !    print'(a,f8.3,f8.3,2es15.5)', 'beta2: myid,i,j,ity,jty: '//info_ij,bond_coeff,rij(0),Krcoef,rhh(0)
   !endif

   ff(1:3) = Krcoef*rij(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(j,1:3) = f(j,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   bond_coeff = this%alpha_bond
   if(potential_type==4) then
     if(rik(0)<this%beta_s1.or.this%beta_l1<rik(0)) bond_coeff = bond_coeff + this%beta_1
     if(rik(0)<this%beta_s2.or.this%beta_l2<rik(0)) bond_coeff = bond_coeff + this%beta_2
   endif

   !Krcoef = bond_coeff * this%Kr*(rik(0)-this%r0) ! ik
   Krcoef = bond_coeff * get_smooth_spring_coeff(rik(0),this%Kr,this%r0) ! ik
   !if(potential_type==4) then
   !  if(rik(0)<this%beta_s1.or.this%beta_l1<rik(0)) &
   !     print'(a,f8.3,f8.3,2es15.5)','beta1: myid,i,k,ity,kty: '//info_ij,bond_coeff,rik(0),Krcoef,rhh(0)
   !  if(rik(0)<this%beta_s2.or.this%beta_l2<rik(0)) &
   !     print'(a,f8.3,f8.3,2es15.5)','beta2: myid,i,k,ity,kty: '//info_ik,bond_coeff,rik(0),Krcoef,rhh(0)
   !endif

   !if(rij(0)<0.85d0) print'(a,3i9,3f8.3)','OUTLIER: j-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rij(0)>1.15d0) print'(a,3i9,3f8.3)','OUTLIER: j-atom is too far',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rik(0)<0.85d0) print'(a,3i9,3f8.3)','OUTLIER: k-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rik(0)>1.15d0) print'(a,3i9,3f8.3)','OUTLIER: k-atom is too far ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rij(0)<bstat%dOH_min) print'(a,3i9,3f8.3)','OUTLIER: j-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rij(0)>bstat%dOH_max) print'(a,3i9,3f8.3)','OUTLIER: j-atom is too far',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rik(0)<bstat%dOH_min) print'(a,3i9,3f8.3)','OUTLIER: k-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   !if(rik(0)>bstat%dOH_max) print'(a,3i9,3f8.3)','OUTLIER: k-atom is too far ',igid,jgid,kgid,rij(0),rik(0),rhh(0)

   if(rij(0)>bstat%dOH_max) bstat%num_oh_max = bstat%num_oh_max + 1
   if(rik(0)>bstat%dOH_max) bstat%num_oh_max = bstat%num_oh_max + 1

   if(rij(0)<bstat%dOH_min) bstat%num_oh_min = bstat%num_oh_min + 1
   if(rik(0)<bstat%dOH_min) bstat%num_oh_min = bstat%num_oh_min + 1

   if(rhh(0)>bstat%dHH_max) bstat%num_hh_max = bstat%num_hh_max + 1
   if(rhh(0)<bstat%dHH_min) bstat%num_hh_min = bstat%num_hh_min + 1

   if(roo_min<bstat%dOO_min) bstat%num_oo_min = bstat%num_oo_min + 1

   if(rij(0)<this%stop_OH_min) print'(a,3i9,3f8.3)','ERROR: j-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   if(rij(0)>this%stop_OH_max) print'(a,3i9,3f8.3)','ERROR: j-atom is too far',igid,jgid,kgid,rij(0),rik(0),rhh(0)

   if(rik(0)<this%stop_OH_min) print'(a,3i9,3f8.3)','ERROR: k-atom is too close ',igid,jgid,kgid,rij(0),rik(0),rhh(0)
   if(rik(0)>this%stop_OH_max) print'(a,3i9,3f8.3)','ERROR: k-atom is too far ',igid,jgid,kgid,rij(0),rik(0),rhh(0)

   call assert(rij(0)>this%stop_OH_min, 'ERROR: j-atom is too close '//info_ij, val=rij(0))
   call assert(rij(0)<this%stop_OH_max, 'ERROR: j-atom is too far '//info_ij, val=rij(0))

   call assert(rik(0)>this%stop_OH_min, 'ERROR: k-atom is too close '//info_ik, val=rik(0))
   call assert(rik(0)<this%stop_OH_max, 'ERROR: k-atom is too far '//info_ik, val=rik(0))

   ff(1:3) = Krcoef*rik(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(k,1:3) = f(k,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   Kqcoef = this%alpha_angle * this%Kq*(theta*180d0/pi-this%q0)*(-1.d0/sine)

   rij(1:3)=-rij(1:3)
   !print*,Kqcoef, theta*180d0/pi, this%q0
   call ForceA3(Kqcoef,j,i,k,rij,rik, f)
enddo

call bstat%print(myid)

contains 

!-----------------------------------------------------------------------------
subroutine set_force_zero(rc_inner, rc_outer, rc_hh_min)
!-----------------------------------------------------------------------------
real(8),intent(in):: rc_inner, rc_outer, rc_hh_min
character(len=:),allocatable :: info_ij, info_ik

real(8) :: rhh(0:3)

do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   !if(myid==0) print'(a,2i6)','ity, this%otype:',ity,this%otype

   igid = l2g(atype(i))

   is_found_j=.false.; jgid=-1; dr_j=0d0; jty=-1
   is_found_k=.false.; kgid=-1; dr_k=0d0; kty=-1

   do l1 = 1, nbrlist(i,0)

      l = nbrlist(i,l1)
      lty = nint(atype(l))

      if(lty/=this%htype) cycle

      if(l2g(atype(l))==igid+1) then
        j = l   !H1
        is_found_j = .true.
        jgid = l2g(atype(j))
        jty = nint(atype(j))
      endif
      if(l2g(atype(l))==igid+2) then
        k = l   !H2
        is_found_k = .true.
        kgid = l2g(atype(k))
        kty = nint(atype(k))
      endif
   enddo

   info_ij = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(jgid)//' '//int_to_str(ity)//' '//int_to_str(jty)
   info_ik = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(kgid)//' '//int_to_str(ity)//' '//int_to_str(kty)

   call assert(is_found_j, 'ERROR(sfz): sfz wrong global ID for j '//info_ij)
   call assert(is_found_k, 'ERROR(sfz): wrong global ID for k '//info_ij)

   call assert(jty==this%htype, 'ERROR(sfz): j-atom type is not H '//info_ij)
   call assert(kty==this%htype, 'ERROR(sfz): k-atom type is not H '//info_ij)

   rij(1:3) = pos(i,1:3)-pos(j,1:3)
   rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))

   ! setting ML force zero beyond ranges rc_inner and rc_outer for O-H1
   if(rij(0)<rc_inner.or.rc_outer<rij(0)) then
     !print'(a,f6.3,7e11.3)','fzero: myid,i,j,ity,jty,rij,fi,fj: '//info_ij,rij(0),f(i,1:3),f(j,1:3)
     !f(i,1:3)=0.d0
     f(j,1:3)=0.d0
   endif

   rik(1:3) = pos(i,1:3)-pos(k,1:3)
   rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))

   ! setting ML force zero beyond ranges rc_inner and rc_outer for O-H2
   if(rik(0)<rc_inner.or.rc_outer<rik(0)) then
     !print'(a,f6.3,7e11.3)','fzero: myid,i,k,ity,kty,rik,fi,fk: '//info_ik,rik(0),f(i,1:3),f(k,1:3)
     !f(i,1:3)=0.d0
     f(k,1:3)=0.d0
   endif

   ! setting ML force zero if H-H bond is below cutoff
   rhh(1:3)=pos(j,1:3)-pos(k,1:3)
   rhh(0) = sqrt(sum(rhh(1:3)*rhh(1:3)))
   !if(rhh(0) < 1.44d0) then
   if(rhh(0) < rc_hh_min) then
     f(j,1:3)=0d0
     f(k,1:3)=0d0
   endif

enddo


end subroutine

!!-----------------------------------------------------------------------------
!subroutine force_cutoff(num_atoms, atype, f, myrank, sr)
!!-----------------------------------------------------------------------------
!   implicit none
!   integer,intent(in) :: num_atoms, myrank
!   real(8),allocatable,intent(in) ::  atype(:)
!   real(8),allocatable,intent(in out) ::  f(:,:)
!   type(short_repulsion_type),intent(in) :: sr
!
!   integer :: i,ia,ity
!   real(8) :: fcut, fset
!
!   do i = 1, num_atoms
!      ity = nint(atype(i))
!
!      if(ity==sr%p2%htype) then
!          fcut = sr%p2%fcut_h
!          fset = fcut*sr%p2%ffactor
!      else if(ity==sr%p2%otype) then
!          fcut = sr%p2%fcut_o
!          fset = fcut*sr%p2%ffactor
!      else
!          print*,'ERROR: atomtype was not found.', myid, ity, i, l2g(atype(i))
!          stop
!      endif
!
!      do ia=1,3
!         if(f(i,ia) > fcut) then
!           print'(a,5i9,es15.5,2f8.2)','INFO: max force cutoff applied.', myid, ity, i, l2g(atype(i)), ia, f(i,ia), fcut, fset
!           f(i,ia) = fset
!         endif
!         if(f(i,ia) < -fcut) then
!           print'(a,5i9,es15.5,2f8.2)','INFO: min force cutoff applied.', myid, ity, i, l2g(atype(i)), ia, f(i,ia), fcut, fset
!           f(i,ia) = -fset
!         endif
!      enddo
!
!   enddo
!end subroutine


!-----------------------------------------------------------------------------
function get_smooth_spring_coeff(r, Kr, r0_oh) result(NewSp)
!-----------------------------------------------------------------------------
!real(8),parameter :: A=9580.16d0, B=0.0286d0 ! r0_oh = 0.971d0
!real(8),parameter :: A=24520.6d0, B=0.049d0  ! r0_oh = 0.971d0
!real(8),parameter :: A=92294.2d0, B=-0.02d0    ! r0_oh = 0.970d0
!real(8),parameter :: A=92294.2d0, B=0.05d0    ! r0_oh = 0.970d0 !<- (9/11)
!real(8),parameter :: A=92294.2d0, B=0.02d0    ! r0_oh = 0.970d0 !<- (9/12)
!real(8),parameter :: A=92294.2d0, B=0.05d0    ! r0_oh = 0.970d0 !<- (9/13)
!real(8),parameter :: A=92294.2d0, B=0.03d0    ! r0_oh = 0.970d0 !<- (9/21)
real(8),parameter :: A=92294.2d0, B=0.05d0    ! r0_oh = 0.970d0 !<- (9/22)
!real(8),parameter :: A=277.78d0, B=0.0d0    ! r0_oh = 0.970d0
real(8),intent(in) :: r, Kr, r0_oh

real(8) :: S, Sp, F, Fp, NewSp

!Energy = S*F

F = A*(r - r0_oh - B)**4
Fp = 4d0*A*(r - r0_oh - B)**3
!F = A*(r - r0_oh)**2
!Fp = 2d0*A*(r - r0_oh)
F = 1.d0
Fp = 0.d0

S = 0.5d0*Kr*(r - r0_oh)**2
Sp = Kr*(r - r0_oh)

NewSp = S*Fp + Sp*F

return

end function

end subroutine

!-----------------------------------------------------------------------------
subroutine flip_frozenh_velocity(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this

integer :: i, ity
real(8) :: vv

if(.not. this%does_flip_velocity) return

call update_frozen_atom_flag(this)

do i=1, NATOMS
   if(this%frozen_atom_flag(i)) then
     ity = nint(atype(i))
     call assert(ity==this%htype, 'wrong atom type to flip velocity '//int_to_str(ity)//' '//int_to_str(i), myid)
     vv = sqrt(sum(v(i,1:3)*v(i,1:3)))
     !v(i,1:3) = vv*this%r_oh(i,1:3)
     v(i,1:3) = this%vscale*vv*this%r_oh(i,1:3)
     print'(a,3i6,4f8.3,1x,3(3f8.3,1x))','V: myid,i,ity,this%r_oh(0:3): ', &
        myid,l2g(atype(i)),ity,this%r_oh(i,0:3), pos(i,1:3),v(i,1:3),f(i,1:3)
   endif
enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine freeze_h_in_overstreched_bond(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this

integer :: i, ity

if(.not. this%does_freeze_atoms) return

call update_frozen_atom_flag(this)

do i=1, NATOMS
   if(this%frozen_atom_flag(i)) then

     ity = nint(atype(i))
     call assert(ity==this%htype, 'wrong atom type to freeze position '//int_to_str(ity)//' '//int_to_str(i), myid)

     pos(i,1:3) = pos(i,1:3)-dt*v(i,1:3)

     print'(a,3i6,4f8.3,1x,3(3f8.3,1x))','X: myid,i,ity,this%r_oh(0:3): ', &
        myid,l2g(atype(i)),ity,this%r_oh(i,0:3),pos(i,1:3),v(i,1:3),f(i,1:3)
   endif
enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine update_frozen_atom_flag(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this

real(8) :: ff0(0:3), fr, rij(0:3), rik(0:3)
integer :: i, j, k, l, ity, lty, l1, igid, jgid

this%frozen_atom_flag = .false.

do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   igid = l2g(atype(i))

   if(ity==this%otype) then
   
      do l1 = 1, nbrlist(i,0)
   
         l = nbrlist(i,l1)
         lty = nint(atype(l))
   
         if(lty/=this%htype) cycle
   
         if(l2g(atype(l))==igid+1) j = l   !H1
         if(l2g(atype(l))==igid+2) k = l   !H2
      enddo
   
      rij(1:3) = pos(i,1:3)-pos(j,1:3)
      rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
      rij(1:3) = rij(1:3)/rij(0)
      if(rij(0) > this%rc_freeze) then
         this%frozen_atom_flag(j) = .true.
         this%r_oh(j,0:3) = rij(0:3)
         print'(a,3i6,4f8.3,2(3f8.3,1x))','o->h1:myid,i,j,rij(0:3),pos(i,1:3),pos(j,1:3) ', &
            myid,l2g(atype(i)),l2g(atype(j)),rij(0:3),pos(i,1:3),pos(j,1:3)
      endif
   
      rik(1:3) = pos(i,1:3)-pos(k,1:3)
      rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
      rik(1:3) = rik(1:3)/rik(0)
      if(rik(0) > this%rc_freeze) then
         this%frozen_atom_flag(k) = .true.
         this%r_oh(k,0:3) = rik(0:3)
         print'(a,3i6,4f8.3,1x,2(3f8.3,1x))','o->h2:myid,i,k,rik(0:3),pos(i,1:3),pos(k,1:3) ', &
            myid,l2g(atype(i)),l2g(atype(k)),rik(0:3),pos(i,1:3),pos(k,1:3)
      endif
   
   else if(ity==this%htype) then
   
      do l1 = 1, nbrlist(i,0)
         l = nbrlist(i,l1)
         lty = nint(atype(l))
   
         if(lty/=this%otype) cycle
         jgid = l2g(atype(l))
   
         if(mod(igid,3)==2 .and. jgid == igid - 1) j = l ! O for H1
         if(mod(igid,3)==0 .and. jgid == igid - 2) j = l ! O for H2
      enddo
   
      rij(1:3) = pos(j,1:3)-pos(i,1:3)        ! r_O - r_H
      rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
      rij(1:3) = rij(1:3)/rij(0)

      if(rij(0) > this%rc_freeze) then
         this%frozen_atom_flag(i) = .true.    ! freeze H
         this%r_oh(i,0:3) = rij(0:3)          ! r_OH
         print'(a,3i6,4f8.3,2(3f8.3,1x))','h->o:myid,i,j,rij(0:3),pos(i,1:3),pos(j,1:3) ', &
            myid,l2g(atype(i)),l2g(atype(j)),rij(0:3),pos(j,1:3),pos(i,1:3)
      endif
   
   endif

enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine apply_type2_short_repulsion(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in) :: this

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, k, j1, k1, ity, jty

integer :: nbr_count, nbr_index(10)
real(8) :: nbr_dist(10), dist_1st, dist_2nd, rij(0:3), rik(0:3)
real(8) :: sine, cosine, theta, Krcoef, Kqcoef, ff(3)

real(8),parameter :: MAX_DIST = 1e9, pi=3.14159265358979d0

do i=1, NATOMS

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   !if(myid==0) print'(a,2i6)','ity, this%otype:',ity,this%otype

   nbr_count = 0; nbr_dist = MAX_DIST
   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jty = nint(atype(j))

      if(jty/=this%htype) cycle

      dr(1:3) = pos(i,1:3) - pos(j,1:3)
      dr2 = sum(dr(1:3)*dr(1:3))
      dr1 = sqrt(dr2)
      dri = 1d0/dr1

      if(dr1>this%rc_nbr) cycle

      nbr_count = nbr_count + 1
      nbr_index(nbr_count) = j
      nbr_dist(nbr_count) = dr1

      !if(myid==0) print'(a,i6,2i6,f8.3)','jty, this%htype:',i,jty,this%htype, dr1

      !f(i,1:3) = f(i,1:3) + fcoef*dr(1:3)
   enddo
   !print*,nbr_index(1:nbr_count)

   j1 = minloc(nbr_dist,dim=1)
   dist_1st = nbr_dist(j1)
   nbr_dist(j1) = MAX_DIST
   j = nbr_index(j1) 

   k1 = minloc(nbr_dist,dim=1)
   dist_2nd = nbr_dist(k1)
   nbr_dist(k1) = MAX_DIST
   k = nbr_index(k1) 

   if(nbr_count > 2) then
      print'(a)',repeat('-=',30)
      print'(a,2i6)', 'ERROR: more than two neighbor hydrogen found : ', myid, nbr_count
      print'(a)',repeat('-=',30)
   endif
   if(nbr_count < 2) then
      print'(a)',repeat('-=',30)
      print'(a,2i6)', 'ERROR: less than two neighbor hydrogen found : ', myid, nbr_count
      print'(a)',repeat('-=',30)
   endif

   rij(1:3) = pos(i,1:3)-pos(j,1:3)
   rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
   rij(1:3) = rij(1:3)/rij(0)

   rik(1:3) = pos(i,1:3)-pos(k,1:3)
   rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
   rik(1:3) = rik(1:3)/rik(0)

   cosine = sum(rij(1:3)*rik(1:3))
   theta = acos(cosine)
   sine = sin(theta)

   Krcoef = this%alpha_bond * this%Kr*(rij(0)-this%r0) ! ij
   ff(1:3) = Krcoef*rij(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(j,1:3) = f(j,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   Krcoef = this%alpha_bond * this%Kr*(rik(0)-this%r0) ! ik
   ff(1:3) = Krcoef*rik(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(k,1:3) = f(k,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   call assert(abs(sine)>1d-9, 'ERROR: too small H-O-H angle', myid)

   Kqcoef = this%alpha_angle * this%Kq*(theta*180d0/pi-this%q0)*(-1.d0/sine)

   rij(1:3)=-rij(1:3)
   !print*,Kqcoef, theta*180d0/pi, this%q0
   call ForceA3(Kqcoef,j,i,k,rij,rik, f)

   !if(myid==0) print'(a,3i6,3f8.3)','i,j,k,dist_1st,dist_2nd, theta: ', &
   !    i, j, k, dist_1st, dist_2nd, theta


enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine apply_type1_short_repulsion(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type1_params),intent(in) :: this

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, j1, ity, jty

do i=1, NATOMS

   ity = nint(atype(i))

   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jty = nint(atype(j))

      dr(1:3) = pos(i,1:3) - pos(j,1:3)
      dr2 = sum(dr(1:3)*dr(1:3))
      dr1 = sqrt(dr2)
      dri = 1d0/dr1
      dri_eta = dri**(this%eta(ity,jty)+1)
      fcoef = this%eta(ity,jty)*this%coef(ity,jty)*dri_eta*dri
      !if(l2g(atype(i))==1.and.dr1<3.d0)
      !print'(2i6,2f12.5)',l2g(atype(i)),l2g(atype(j)),dr1,fcoef

      f(i,1:3) = f(i,1:3) + fcoef*dr(1:3)

   enddo
enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine short_repulsion(sr)
!-----------------------------------------------------------------------------
!type(short_repulsion_type),intent(in) :: sr
type(short_repulsion_type),intent(in out) :: sr ! intent(out) because of f_spring. move it out from sr?
character(len=:),allocatable :: filebase
real(8),allocatable :: ftmp(:,:)

if(.not. sr%has_short_repulsion) return

if (sr%potential_type==1) then
  call sr%p1%apply()
else if (sr%potential_type==2) then
  call sr%p2%apply2()
else if (sr%potential_type==3.or.sr%potential_type==4) then
  call sr%p2%apply4(sr%potential_type, bstat, lj_pot)
else if (sr%potential_type==5) then
  call sr%p2%apply5(sr%potential_type, bstat)
else
  print*,'To be implemented'

endif

end subroutine

end module

module fnnin_parser

  use utils, only : assert, pi, getstr, getfilenamebase, l2g, Eev_kcal, Ekcal_j
  use base, only : force_field_class
  use fileio, only : output, writebin
  use velocity_modifiers_mod, only : vkick, linear_momentum, &
                                     scale_to_target_temperature, &
                                     maximally_preserving_bd

  use iso_fortran_env, only: int32, int64, real32, real64 

  implicit none

  integer,parameter :: rk = real32 ! rk = real64
  integer,parameter :: ik = int32 ! ik = int64

  type :: model_params
    character(len=:),allocatable :: filename
    character(len=:),allocatable :: element
    real(rk) :: mass

    ! maximum angular cutoff within a model for neighborlist construction
    real(rk) :: rad_rc
  end type

  type, extends(force_field_class) :: fnn_param
    type(model_params), allocatable :: models(:) 
    contains 
       procedure :: print => fnn_param_print
  end type

  type(fnn_param),target :: fnn_param_obj

  character(len=:),allocatable,private :: token

contains

!------------------------------------------------------------------------------
  function get_max_cutoff(fp) result(max_rc)
!------------------------------------------------------------------------------
    type(fnn_param),intent(in) :: fp 
    real(rk) :: max_rc
    integer :: ia, ib

    max_rc = -1.d0

    do ia=1,size(fp%models)
       if(max_rc<fp%models(ia)%rad_rc) max_rc = fp%models(ia)%rad_rc
    enddo

    return 
  end function

!------------------------------------------------------------------------------
  subroutine get_tokens_and_append_model(linein, models)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    type(model_params),allocatable,intent(in out) :: models(:)
    type(model_params) :: mbuf

    if (getstr(linein, token) < 0) stop 'error while reading element name'
    mbuf%element = trim(adjustl(token))
    if (getstr(linein, token) < 0) stop 'error while reading element mass'
    read(token, *) mbuf%mass
    if (getstr(linein, token) < 0) stop 'error while reading model filename'
    mbuf%filename = trim(adjustl(token))
    read(token, *) mbuf%filename
    if (getstr(linein, token) < 0) stop 'error while reading radial cutoff'
    read(token, *) mbuf%rad_rc

    ! allocate zero-sized array
    if(.not.allocated(models)) allocate(models(0)) 
    models = [models, mbuf]

    return
  end subroutine

  function get_index_of_model(element, models) result(idx)
    type(model_params),allocatable,intent(in) :: models(:)
    character(len=:),allocatable,intent(in) :: element
    integer :: idx
    do idx=1, size(models)
       if(models(idx)%element == element) return
    enddo
    idx = -1
    return
  end function

!------------------------------------------------------------------------------
  function fnn_param_ctor(path) result(c)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in) :: path
    character(256) :: linein0
    character(len=:),allocatable :: linein 

    integer :: iunit
    type(fnn_param) :: c 

    open(newunit=iunit, file=path, status='old', form='formatted')

    ! find how many models exist. 
    do while (.true.)
      read(iunit,'(a)',end=10) linein0
      linein = trim(adjustl(linein0))

      if(getstr(linein, token) > 0) then
         if(token=='model') call get_tokens_and_append_model(linein, c%models)
      endif
    end do
    10 rewind(iunit)

    if (size(c%models)<=0) stop 'ERROR: at least one model must be defined.'

  end function

!------------------------------------------------------------------------------
  subroutine fnn_param_print(this)
!------------------------------------------------------------------------------
    class(fnn_param), intent(in) :: this
    integer(ik) :: ia,ib,ic,i1,i2,i3

    write(*,fmt='(a)') repeat('=',80)
    do ia=1, size(this%models)

       associate(m=>this%models(ia))
          print'(a,i3,2a,1x,f10.3,3x,a,1x,f8.3)', &
             'element,mass,filename,rad_cut: ', & 
             get_index_of_model(m%element,this%models),'-',m%element, &
             m%mass, m%filename, m%rad_rc
       end associate
      
    enddo
    write(*,fmt='(a)') repeat('=',80)

  end subroutine

end module

!------------------------------------------------------------------------------
module fnn
!------------------------------------------------------------------------------

  use fnnin_parser

  use utils, only : pi, int_to_str
  use memory_allocator_mod
  use mpi_mod

  use base
  use lists_mod, only : getnonbondingmesh, linkedlist
  use communication_mod, only : copyatoms
  use msd_mod, only : msd_data, msd_add_initial_pos, msd_measure, msd_save

  use aenet

  use mod_short_repulsion, only : short_repulsion, short_rep 
 
  implicit none

  integer(ik) :: num_types = 0, num_pairs = 0

  integer,parameter :: NMINCELL_FNN=1

  integer,allocatable :: pair_types(:,:)

  ! timing for 1-feature calc & 2-force inference
  real(rk),save,private :: tstart(0:3)=0.0, tfinish(0:3)=0.0
  real(8) :: Epot

contains

!------------------------------------------------------------------------------
subroutine get_force_fnn(ff, num_atoms, atype, pos, f, q)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff
type(fnn_param),pointer :: fp => null()
integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

real(8) :: coo_i(3),  coo_j(3, MAXNEIGHBS), E_i, f3r(3, NBUFFER), F_i(3, num_atoms)
integer :: type_i, index_i, type_j(MAXNEIGHBS), index_j(MAXNEIGHBS)
integer :: i,j,j1,ity,n_j,stat


! not sure if this is the best way, but binding force_field_class to fnn_parm
select type(ff); type is (fnn_param) 
   fp => ff
end select

call COPYATOMS(imode = MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos, ipos=ipos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

call get_nbrlist_fnn(num_atoms, atype, pos, fp) 

call cpu_time(tstart(1))
f = 0.d0; F_i = 0.d0; f3r = 0.d0;
Epot = 0.d0
do i = 1, num_atoms
   coo_i(1:3) = pos(i,1:3)
   type_i = nint(atype(i))
   index_i = i
   n_j = nbrlist(i,0)

   !print'(a,3f8.3,3i6)',atmname(type_i),coo_i, type_i, index_i, n_j
   do j1 = 1, nbrlist(i,0)
      j = nbrlist(i,j1)
      coo_j(1:3,j1) = pos(j,1:3)
      type_j(j1) = nint(atype(j))
      index_j(j1) = j

      !print'(a3,3f8.3,2i6,f8.3)','j: ',coo_j(1:3,j1), type_j(j1), index_j(j1)
   enddo

   call aenet_atomic_energy_and_forces( &
        coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, num_atoms, &
        E_i, F_i, f3r, NBUFFER, stat) 

   Epot = Epot + E_i

   !print'(a,i5,4es15.5)','i,E_i, F_i: ', i, E_i*Eev_kcal, F_i(1:3,i)*Eev_kcal
enddo
call cpu_time(tfinish(1))

Epot = Epot*Eev_kcal

!omp simd
do i = 1, num_atoms
   f(i,1:3) = f (i,1:3) + F_i(1:3,i)*Eev_kcal
enddo
!omp simd
do i = 1, NBUFFER
   f(i,1:3) = f(i,1:3) + f3r(1:3,i)*Eev_kcal
enddo

!call force_cutoff(fcut_o=70d0, fcut_h=50d0, ffactor=0.7d0) ! before October 06, 9:45pm
!call force_cutoff(fcut_o=70d0, fcut_h=50d0, ffactor=0.0d0) ! before October 09, 12:04pm
!call force_cutoff(fcut_o=70d0, fcut_h=50d0, ffactor=0.7d0)
call force_cutoff(fcut_o=short_rep%p2%fcut_o, fcut_h=short_rep%p2%fcut_h, ffactor=short_rep%p2%ffactor)
call short_repulsion(short_rep)

CALL COPYATOMS(imode=MODE_CPBK, dr=dr_zero, atype=atype, pos=pos, f=f, q=q)

contains

subroutine force_cutoff(fcut_o, fcut_h, ffactor)
   real(8),intent(in) :: fcut_o, fcut_h, ffactor
   integer :: i,ia,ity
   real(8) :: fcut, fset
   
   do i = 1, NATOMS
      ity = nint(atype(i))

      if(ity==1) then
          fcut = fcut_h
          fset = fcut*ffactor
      else if(ity==2) then
          fcut = fcut_o
          fset = fcut*ffactor
      else
          print*,'ERROR: atomtype was not found.', myid, ity, i, l2g(atype(i))
          stop 
      endif

      do ia=1,3
         if(f(i,ia) > fcut) then
           !print'(a,5i9,es15.5,2f8.2)','max force cutoff applied.', myid, ity, i, l2g(atype(i)), ia, f(i,ia), fcut, fset
           f(i,ia) = fset
         endif
         if(f(i,ia) < -fcut) then
           !print'(a,5i9,es15.5,2f8.2)','min force cutoff applied.', myid, ity, i, l2g(atype(i)), ia, f(i,ia), fcut, fset
           f(i,ia) = -fset
         endif
      enddo

   enddo
end subroutine

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_fnn(mdbase, num_mdsteps) 
use utils, only : UTEMP, UTEMP0
use velocity_modifiers_mod, only : gaussian_dist_velocity, adjust_temperature, scale_temperature
!------------------------------------------------------------------------------
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
real(8) :: ctmp,cpu0,cpu1,cpu2,comp=0.d0

integer :: i,ity

if(reset_velocity_random) call gaussian_dist_velocity(atype, v)

if(mdmode==0) then
  call gaussian_dist_velocity(atype, v)
  call WriteBIN(atype, pos, v, q, GetFileNameBase(DataDir,-1))
  return
endif

call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)

call cpu_time(cpu0)

!--- set force model
do nstep=0, num_mdsteps-1

  if(mod(nstep,pstep)==0) call print_e_fnn(atype, v, q)

  call cpu_time(tstart(0))

  if(mod(nstep,fstep)==0) &
        call OUTPUT(GetFileNameBase(DataDir,current_step+nstep), atype, pos, v, q, v)
        !call OUTPUT(GetFileNameBase(DataDir,current_step+nstep), atype, pos, v, q, f)

  if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:NATOMS,1:3)=vsfact*v(1:NATOMS,1:3)

   if(mod(nstep,sstep)==0.and.mdmode==5) &
      call scale_to_target_temperature(atype, v, treq)

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call gaussian_dist_velocity(atype, v)

!--- element-wise velocity scaling
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call scale_temperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==8) &
      call adjust_temperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==9) &
      call maximally_preserving_BD(atype, v, vsfact) 

!--- MSD measurements
   call msd_add_initial_pos(msd_data, nstep, NATOMS, pos, ipos)
   call msd_measure(msd_data, nstep, NATOMS, atype, pos, ipos)

!--- total force may not be zero with FNN. fix linear momentum every pstep.
   if(mod(nstep,pstep)==0) call linear_momentum(atype, v)

!--- update velocity & position
   call vkick(1.d0, atype, v, f)

   pos(1:natoms,1:3)=pos(1:natoms,1:3)+dt*v(1:natoms,1:3)
   if(short_rep%has_short_repulsion .and. short_rep%potential_type == 3) call short_rep%p2%freezex()

!--- migrate atoms after positions are updated
   call COPYATOMS(imode=MODE_MOVE_FNN,dr=dr_zero,atype=atype,pos=pos, &
                  v=v,f=f,q=q,ipos=ipos)

   call cpu_time(cpu1)
   call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)
   call cpu_time(cpu2)
   comp = comp + (cpu2-cpu1)

!--- update velocity
   call vkick(1.d0, atype, v, f)

   if(short_rep%has_short_repulsion .and. short_rep%potential_type == 3) call short_rep%p2%flipv()

   call cpu_time(tfinish(0))

enddo

!--- save the final configurations
!call OUTPUT(GetFileNameBase(DataDir,current_step+nstep), atype, pos, v, q, f)
call OUTPUT(GetFileNameBase(DataDir,current_step+nstep), atype, pos, v, q, v)

!--- update rxff.bin in working directory for continuation run
call WriteBIN(atype, pos, v, q, GetFileNameBase(DataDir,-1))

!--- save result if msd_data%is_msd == true
call msd_save(msd_data)

call cpu_time(cpu2)
if(myid==0) print'(a,2f12.5)','comp, total (sec): ', comp, cpu2-cpu0

return
end subroutine

!------------------------------------------------------------------------------
subroutine get_nbrlist_fnn(num_atoms, atype, pos, fp)
!------------------------------------------------------------------------------
integer,intent(in) :: num_atoms
real(8),intent(in),allocatable :: atype(:), pos(:,:)
type(fnn_param),intent(in out) :: fp

real(rk) :: rr(3), rr2, rij, dsum 
integer(ik) :: i, j, k, i1, j1, k1, l1, l2, l3, l4, ii
integer(ik) :: c1,c2,c3,ic(3),c4,c5,c6,ity,jty,kty,inxn

nbrlist(:,0) = 0

call cpu_time(tstart(2))

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,rr,rr2,rij,fr_ij,rij_mu,eta_ij,idx) 
do c1=0, cc(1)-1
do c2=0, cc(2)-1
do c3=0, cc(3)-1

  i = header(c1, c2, c3)
  do i1=1, nacell(c1, c2, c3)
     ity = nint(atype(i))

     !print'(3i6,i6,3f10.5)',c1,c2,c3,m,pos(m,1:3)

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1+c4, c2+c5, c3+c6]

        j = header(ic(1),ic(2),ic(3))
        do j1=1, nacell(ic(1), ic(2), ic(3))

           if(i/=j) then
             jty = nint(atype(j))

             rr(1:3) = pos(i,1:3) - pos(j,1:3)
             rr2 = sum(rr(1:3)*rr(1:3))
             rij = sqrt(rr2)

             if(rij<fp%models(ity)%rad_rc) then !!! FIXME !!!
               nbrlist(i, 0) = nbrlist(i, 0) + 1
               nbrlist(i, nbrlist(i, 0)) = j
             endif

           endif

           j=llist(j)
        enddo
     enddo; enddo; enddo

     i = llist(i)
  enddo
enddo; enddo; enddo
!$omp end parallel do 
call cpu_time(tfinish(2))

return

end subroutine

!------------------------------------------------------------------------------------------
subroutine get_cutoff_fnn(rcut, rcut2, maxrcut, radial_cutoff)
!------------------------------------------------------------------------------------------
real(8),allocatable,intent(in out) :: rcut(:), rcut2(:)
real(8),intent(in out) :: maxrcut
real(rk),intent(in) :: radial_cutoff

integer :: ity,jty,inxn

!--- get the cutoff length 
call allocator(rcut, 1, num_pairs)
call allocator(rcut2, 1, num_pairs)
call allocator(pair_types, 1, num_types, 1, num_types)

inxn=0
do ity=1, num_types
do jty=ity, num_types
   inxn = inxn + 1
   pair_types(ity,jty) = inxn

   rcut(inxn)  = radial_cutoff
   rcut2(inxn) = radial_cutoff*radial_cutoff
   !print'(a,3i6,2f10.5)','ity, jty, inxn: ', ity, jty, inxn, rcut(inxn), rcut2(inxn)

   pair_types(jty,ity) = pair_types(ity,jty) 
enddo
enddo

maxrcut = maxval(rcut)

end subroutine

!-------------------------------------------------------------------------------------------
subroutine print_e_fnn(atype, v, q)
use mpi_mod
use base, only : hh, hhi, natoms, gnatoms, mdbox, myid, ierr, hmas, wt0
use atoms
use memory_allocator_mod
! calculate the kinetic energy and sum up all of potential energies, then print them.
!-------------------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: v(:,:)

integer :: i,ity,cstep
real(8) :: tt=0.d0, Etotal

ke=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE (MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
ke = ke/GNATOMS
tt = ke*UTEMP
GKE = ke ! FIXME for ctmp = (treq*UTEMP0)/( GKE*UTEMP )

call MPI_ALLREDUCE (MPI_IN_PLACE, Epot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
Epot = Epot/GNATOMS
Etotal = Epot + ke

if(myid==0) then
   
   cstep = nstep + current_step 

   write(6,'(a,i9,es13.5,f10.3,3x,4f10.5)') &
        'MDstep,Etotal,T(K),onestep(sec),force_calc,nbr_calc: ', cstep, Epot + ke, tt, tfinish(0:3)-tstart(0:3)
endif

end subroutine

end module

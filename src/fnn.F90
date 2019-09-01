module fnnin_parser

  use utils, only : assert, pi, getstr, getfilenamebase, l2g, Eev_kcal
  use base, only : force_field_class
  use fileio, only : output, writebin
  use velocity_modifiers_mod, only : vkick, linear_momentum, &
                                     scale_to_target_temperature, &
                                     maximally_preserving_bd

  use iso_fortran_env, only: int32, int64, real32, real64 

  implicit none

  integer,parameter :: rk = real32 ! rk = real64
  integer,parameter :: ik = int32 ! ik = int64

  integer,parameter :: NUM_ANG_MODES=2 ! number of modes in a triplet

  integer,parameter :: SIZE_FEATURE_TABLE=20000

  type :: layer
    real(rk),allocatable :: b(:)
    real(rk),allocatable :: w(:,:)
  end type

  type :: network
    integer(ik),allocatable :: dims(:)
    type(layer),allocatable :: layers(:)
  end type

  type :: feature_stat
    real(rk),allocatable :: mean(:)
    real(rk),allocatable :: stddev(:)
  end type

  type :: rad_feature_type
    real(rk),allocatable,dimension(:) :: mu, eta
    real(rk) :: rc, rdamp

    integer(ik),allocatable :: indices(:)
    integer(ik) :: total

    real(rk),allocatable :: ftab(:,:) ! feature calc table
    real(rk) :: ftab_dr 
  end type

  type :: ang_feature_type
    real(rk),allocatable,dimension(:) :: mu, eta, zeta
    integer(ik),allocatable,dimension(:) :: lambda
    real(rk) :: rc, rdamp

    integer(ik),allocatable :: indices(:)
    integer(ik) :: total
  end type

  type :: model_params
    character(len=:),allocatable :: element
    real(rk) :: mass
    real(rk) :: scaling_factor

    integer(ik),allocatable :: layersize(:)
    type(network),allocatable :: networks(:)
    type(feature_stat),allocatable :: fstat(:)

    type(rad_feature_type),allocatable :: rad(:) ! radial feature params 
    type(ang_feature_type),allocatable :: ang(:) ! angular feature params

    integer(ik),allocatable :: map_rad(:) ! map from (jty) to feature array index
    integer(ik),allocatable :: map_ang(:,:) ! map from (ity,kty) to feature array index

    real(rk),allocatable :: features(:,:,:) 
    integer(ik),allocatable :: feature_ptr_rad(:), feature_ptr_ang(:), num_features(:)

    ! maximum angular cutoff within a model for neighborlist construction
    real(rk) :: max_ang_rc
  end type

  type, extends(force_field_class) :: fnn_param
    type(model_params), allocatable :: models(:) 
    contains 
       procedure :: print => fnn_param_print
  end type

  type(fnn_param),target :: fnn_param_obj

  character(len=:),allocatable,private :: sbuf
  integer,private :: ibuf
  real(rk),private :: rbuf
  character(len=:),allocatable,private :: token

contains

!------------------------------------------------------------------------------
subroutine set_feature_tables_fnn(fp)
!------------------------------------------------------------------------------
type(fnn_param),intent(in out) :: fp

real(rk) :: rij, fc_ij, rij_mu, eta_ij
integer(ik) :: t, ity, jty, ii, l1, l2, idx, num_mu, num_eta
 

do ity=1, size(fp%models)

   do jty=1, size(fp%models)

      ! lookup table. should be ii == jty
     ii = fp%models(ity)%map_rad(jty)

     associate(rad=>fp%models(ity)%rad(ii)) ! rad for (ity,jty) pair
  
     num_mu = size(rad%mu)
     num_eta = size(rad%eta)
  
     if(allocated(rad%ftab)) deallocate(rad%ftab)
     allocate(rad%ftab(SIZE_FEATURE_TABLE,num_mu*num_eta))
  
     rad%ftab_dr = rad%rc/SIZE_FEATURE_TABLE
  
     do t=1, SIZE_FEATURE_TABLE

       rij = rad%ftab_dr*(t-0.5)
       fc_ij = 0.5*(1.0+cos(pi*rij/rad%rdamp))
  
       do l1=1, size(rad%mu)
         rij_mu = rij - rad%mu(l1)
  
         do l2=1, size(rad%eta) 
           eta_ij = exp( -rad%eta(l2) * rij_mu * rij_mu )
           idx = (l1-1)*size(rad%eta) + (l2-1) + 1
           rad%ftab(t,idx) = eta_ij*fc_ij

         enddo
  
       enddo
     enddo
  
     end associate

  enddo
enddo

end subroutine


!------------------------------------------------------------------------------
  function get_max_cutoff(fp) result(max_rc)
!------------------------------------------------------------------------------
    type(fnn_param),intent(in) :: fp 
    real(rk) :: max_rc
    integer :: ia, ib

    max_rc = -1.d0

    do ia=1,size(fp%models)
       do ib=1,size(fp%models(ia)%rad)
          if(max_rc<fp%models(ia)%rad(ib)%rc) max_rc = fp%models(ia)%rad(ib)%rc
       enddo
       !do ib=1,size(fp%models(ia)%ang)
       !   if(max_rc<fp%models(ia)%ang(ib)%rc) max_rc = fp%models(ia)%ang(ib)%rc
       !enddo
    enddo

    return 
  end function

!------------------------------------------------------------------------------
  function rad_ctor_from_line(linein,indices) result(rad)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    integer(ik),allocatable,intent(in) :: indices(:)
    type(rad_feature_type) :: rad

    if(.not.allocated(rad%mu)) allocate(rad%mu(0))
    if(.not.allocated(rad%eta)) allocate(rad%eta(0))
    if(.not.allocated(rad%indices)) allocate(rad%indices(0))
    rad%indices = indices;  rad%total = 0

    do while (getstr(linein, token) > 0)

99     select case(token)

         case('rc')
           call assert(getstr(linein, token)>0, &
                        'ERROR: rc in rad_ctor_from_line(): '//token)
           read(token,*) rad%rc

         case('rdamp')
           call assert(getstr(linein, token)>0, &
                        'ERROR: rdamp in rad_ctor_from_line(): '//token)
           read(token,*) rad%rdamp

         case('mu')
           do while (getstr(linein, token) > 0)
             read(token,*,err=99) rbuf
             rad%mu = [rad%mu, rbuf]
           enddo

         case('eta')
           do while (getstr(linein, token) > 0)
             read(token,*,err=99) rbuf
             rad%eta = [rad%eta, rbuf]
           enddo

         case default
           print'(a)', 'WARNING: unknown field in rad_ctor_from_line(): '//linein

       end select
    enddo

    rad%total = 1
    if(size(rad%eta)>0) rad%total = rad%total*size(rad%eta)
    if(size(rad%mu)>0) rad%total = rad%total*size(rad%mu)

    return
  end function

!------------------------------------------------------------------------------
  function ang_ctor_from_line(linein,indices) result(ang)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein 
    integer(ik),allocatable,intent(in) :: indices(:)
    type(ang_feature_type) :: ang

    if(.not.allocated(ang%lambda)) allocate(ang%lambda(0))
    if(.not.allocated(ang%zeta)) allocate(ang%zeta(0))
    if(.not.allocated(ang%eta)) allocate(ang%eta(0))
    if(.not.allocated(ang%mu)) allocate(ang%mu(0))
    if(.not.allocated(ang%indices)) allocate(ang%indices(0))
    ang%indices = indices;  ang%total = 0

    do while (getstr(linein, token) > 0)

99    select case(token)

         case('rc')
           call assert(getstr(linein, token)>0, &
                       'ERROR: rc in ang_ctor_from_line(): '//token)
           read(token,*) ang%rc

         case('rdamp')
           call assert(getstr(linein, token)>0, &
                       'ERROR: rdamp in ang_ctor_from_line(): '//token)
           read(token,*) ang%rdamp

         case('lambda')
           do while (getstr(linein, token) > 0)
             read(token,*,err=99) ibuf
             ang%lambda = [ang%lambda, ibuf]
           enddo

         case('zeta')
           do while (getstr(linein, token) > 0)
             read(token,*,err=99) rbuf
             ang%zeta = [ang%zeta, rbuf]
           enddo

         case('eta')
           do while (getstr(linein, token) > 0)
             read(token,*,err=99) rbuf
             ang%eta = [ang%eta, rbuf]
           enddo

         case('mu')
           do while (getstr(linein, token) > 0)
             read(token,*,err=99) rbuf
             ang%mu = [ang%mu, rbuf]
           enddo

         case default
           print'(a)', 'WARNING: unknown field in ang_ctor_from_line(): '//linein

      end select
    enddo

    ang%total = NUM_ANG_MODES  ! vibration & stretch
    if(size(ang%lambda)>0) ang%total = ang%total*size(ang%lambda)
    if(size(ang%zeta)>0) ang%total = ang%total*size(ang%zeta)
    if(size(ang%eta)>0) ang%total = ang%total*size(ang%eta)
    if(size(ang%mu)>0) ang%total = ang%total*size(ang%mu)

    return
  end function

!------------------------------------------------------------------------------
  subroutine get_tokens_and_append_model(linein, models)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    type(model_params),allocatable,intent(in out) :: models(:)
    type(model_params) :: mbuf

    if (getstr(linein, token) < 0) stop 'erro while reading element name'
    mbuf%element = trim(adjustl(token))
    if (getstr(linein, token) < 0) stop 'erro while reading element mass'
    read(token, *) mbuf%mass
    if (getstr(linein, token) < 0) stop 'erro while reading scaling factor'
    read(token, *) mbuf%scaling_factor

    ! allocate zero-sized array
    if(.not.allocated(mbuf%layersize)) allocate(mbuf%layersize(0))
    do while( getstr(linein, token) > 0 )
       read(token, *) ibuf
       mbuf%layersize= [mbuf%layersize, ibuf]
    enddo

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

    type(fnn_param) :: c 
    integer(ik) :: iunit, num_models, feature_ptr, ii, ia,ib,ic, ang_index
    integer(ik),allocatable :: indices(:)

    open(newunit=iunit, file=path, status='old', form='formatted')

    ! find how many models exist first. 
    do while (.true.)
      read(iunit,'(a)',end=10) linein0
      linein = trim(adjustl(linein0))

      if(getstr(linein, token) > 0) then
         if(token=='model') call get_tokens_and_append_model(linein, c%models)
      endif
    end do
    10 rewind(iunit)

    if (size(c%models)<=0) stop 'ERROR: at least one model must be defined.'

    ! radial and angular parameters are defined for atomic pairs & triplets
    num_models = size(c%models)

    do ia=1,size(c%models)
       allocate(c%models(ia)%rad(0), c%models(ia)%ang(0))
       allocate(c%models(ia)%map_rad(num_models))
       allocate(c%models(ia)%map_ang(num_models,num_models))
       c%models(ia)%map_rad=0; c%models(ia)%map_ang=0
       c%models(ia)%max_ang_rc=-1.d0
    enddo

    do while (.true.)
      read(iunit,'(a)',end=11) linein0
      linein = trim(adjustl(linein0))

      if(getstr(linein, token) > 0) then
        select case(token)

           case('rad')
             if(getstr(linein,token)>0) ia = get_index_of_model(token,c%models) !<- center atom
             if(getstr(linein,token)>0) ib = get_index_of_model(token,c%models)
             do while(.true.)
               read(iunit,'(a)') linein0
               if(linein0(1:1)=='#') cycle ! skip comment
               if(index(linein0,'end')>0) exit 
               linein = linein//' '//trim(adjustl(linein0))
             enddo
             if(c%models(ia)%map_rad(ib)==0) then
               indices = [ia,ib]
               c%models(ia)%rad = [c%models(ia)%rad, rad_ctor_from_line(linein,indices)]
               c%models(ia)%map_rad(ib) = size(c%models(ia)%rad)
             endif

           case('ang')
             if(getstr(linein,token)>0) ia = get_index_of_model(token,c%models)
             if(getstr(linein,token)>0) ib = get_index_of_model(token,c%models) !<- center atom
             if(getstr(linein,token)>0) ic = get_index_of_model(token,c%models)

             if(c%models(ib)%map_ang(ia,ic)==0) then
               indices = [ia,ib,ic]
               do while(.true.)
                 read(iunit,'(a)') linein0
                 if(linein0(1:1)=='#') cycle ! skip comment 
                 if(index(linein0,'end')>0) exit
                 linein = linein//' '//trim(adjustl(linein0))
               enddo
               c%models(ib)%ang = [c%models(ib)%ang, ang_ctor_from_line(linein,indices)]
               ang_index = size(c%models(ib)%ang)
               c%models(ib)%map_ang(ia,ic) = ang_index
               c%models(ib)%map_ang(ic,ia) = ang_index

               ! keep the max angular cutoff for neighborlist construction
               ! during feature calculation
               if(c%models(ib)%max_ang_rc < c%models(ib)%ang(ang_index)%rc) &
                   c%models(ib)%max_ang_rc = c%models(ib)%ang(ang_index)%rc 
             endif

           case default
        end select
      endif

    end do
    11 close(iunit)


! store the first array index after packing feature values:
!   radial: the first array index for an ity-jty pair is models(ity)%feature_ptr_rad(jty). 
!   angular: the first array index for an ity-jty-kty triplet is models(jty)%feature_ptr_ang(ity or kty). 
    do ia=1,size(c%models)

       associate(m=>c%models(ia)) 

         allocate(m%feature_ptr_rad(0),m%feature_ptr_ang(0),m%num_features(2))
  
         feature_ptr = 1

         do ib=1,size(m%rad)
            m%feature_ptr_rad = [m%feature_ptr_rad, feature_ptr]
            feature_ptr = feature_ptr + m%rad(ib)%total
         enddo
         m%num_features(1) = feature_ptr - 1

         do ib=1,size(m%ang)
            m%feature_ptr_ang = [m%feature_ptr_ang, feature_ptr]
            feature_ptr = feature_ptr + m%ang(ib)%total
         enddo
         m%num_features(2) = feature_ptr - 1 - m%num_features(1)

         !print'(a,5i6)','ia,size(rad),size(ang),m%feature_ptr(1:2): ', &
         !        ia,size(m%feature_ptr_rad),size(m%feature_ptr_ang),m%num_features(1:2)
       end associate

    enddo

  end function

!------------------------------------------------------------------------------
  subroutine fnn_param_print(this)
!------------------------------------------------------------------------------
    class(fnn_param), intent(in) :: this
    integer(ik) :: ia,ib,ic,i1,i2,i3

    do ia=1, size(this%models)

       associate(m=>this%models(ia))

       write(*,*)
       write(*,fmt='(a)') repeat('=',80)
       print'(a,i3,2a,2f10.3)', 'element,mass,scaling_factor : ', & 
            get_index_of_model(m%element,this%models),'-',m%element, m%mass, m%scaling_factor
       print'(a,10i5)','layer size: ', m%layersize
       print'(a,2i6)', 'rad&ang feature types: ', size(m%rad), size(m%ang)
       write(*,fmt='(a,10i6)',advance='no') 'feature_ptr_(rad,ang): ', m%feature_ptr_rad
       write(*,fmt='(a,10i6)') ', ', m%feature_ptr_ang
       print'(a,3i6)', 'num_features(rad,ang),total: ', m%num_features, sum(m%num_features)

       print'(a)',repeat('-',80)
       do ib=1,size(m%rad)
          write(*,fmt='(a)',advance='no') m%element//' '//this%models(ib)%element
          write(*,fmt='(a,2i3,i6)') '    size(mu,eta,total): ', &
             size(m%rad(ib)%mu), size(m%rad(ib)%eta), m%rad(ib)%total
       enddo

       do ib=1,size(m%rad)
          i1 = m%rad(ib)%indices(1)
          i2 = m%rad(ib)%indices(2)
          write(*,fmt='(a,2i2,1x,a,1x)',advance='no') 'rad_params: ', &
             i1,i2,this%models(i1)%element//'-'//this%models(i2)%element

          write(*,fmt='(a,2(a,f8.3))', advance='no') &
               repeat(' ',3), ' rc ', m%rad(ib)%rc, ' rdamp ', m%rad(ib)%rdamp 
          write(*,fmt='(2a)',advance='no') repeat(' ',3), ' mu '
          do ic=1,size(m%rad(ib)%mu);  write(*,fmt='(f6.2)',advance='no') m%rad(ib)%mu(ic);  enddo
          write(*,fmt='(2a)',advance='no') repeat(' ',3), ' eta '
          do ic=1,size(m%rad(ib)%eta);  write(*,fmt='(f6.2)',advance='no') m%rad(ib)%eta(ic);  enddo
          write(*,*)
       enddo

       print'(a)',repeat('-',60)
       do ib=1,size(m%ang)
          i1 = m%ang(ib)%indices(1);  i2 = m%ang(ib)%indices(2);  i3 = m%ang(ib)%indices(3)
          write(*,fmt='(a,a)', advance='no') this%models(i1)%element//' '//m%element//' '//this%models(i3)%element
          write(*,fmt='(a,i3)', advance='no') '  ang_modes: ', NUM_ANG_MODES 
          write(*,fmt='(a,4i3,i6)') '   size(lambda,zeta,mu,eta,total): ', &
                size(m%ang(ib)%lambda),size(m%ang(ib)%zeta),size(m%ang(ib)%mu),size(m%ang(ib)%eta), m%ang(ib)%total
       enddo
       do ib=1,size(m%ang)
          i1 = m%ang(ib)%indices(1);  i2 = m%ang(ib)%indices(2);  i3 = m%ang(ib)%indices(3)
          write(*,fmt='(a,3i2,1x,a,1x)', advance='no') 'ang_params: ',i1,i2,i3,&
             this%models(i1)%element//'-'//m%element//'-'//this%models(i3)%element

          write(*,fmt='(a,2(a,f8.3))',advance='no') &
               repeat(' ',3), ' rc ', m%ang(ib)%rc, ' rdamp ', m%ang(ib)%rdamp 
          write(*,fmt='(2a)',advance='no') repeat(' ',3), ' lambda '
          do ic=1,size(m%ang(ib)%lambda);  write(*,fmt='(i3)',advance='no') m%ang(ib)%lambda(ic);  enddo

          write(*,fmt='(2a)',advance='no') repeat(' ',3), ' zeta '
          do ic=1,size(m%ang(ib)%zeta);  write(*,fmt='(f6.2)',advance='no') m%ang(ib)%zeta(ic);  enddo

          write(*,fmt='(2a)',advance='no') repeat(' ',3), ' mu '
          do ic=1,size(m%ang(ib)%mu);  write(*,fmt='(f6.2)',advance='no') m%ang(ib)%mu(ic);  enddo

          write(*,fmt='(2a)',advance='no') repeat(' ',3), ' eta '
          do ic=1,size(m%ang(ib)%eta);  write(*,fmt='(f6.2)',advance='no') m%ang(ib)%eta(ic);  enddo

          write(*,*)
       enddo

       if(m%layersize(1) /= sum(m%num_features)) then
          write(*,fmt='(a)') repeat('!',80)
          write(*,fmt='(a,i5,a,i5)') 'ERROR: the size of input layer in fnn.in should be ', &
                                     sum(m%num_features),'   but ', m%layersize(1)
          write(*,fmt='(a)') repeat('!',80)
          stop
       endif

       write(*,fmt='(a)') repeat('=',80)

       end associate
      
    enddo
    write(*,*) 

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
  use lists_mod, only : get_mesh_for_nonbonding_list, linkedlist
  use communication_mod, only : copyatoms
  use msd_mod, only : msd_data

  implicit none


  integer(ik),parameter :: num_forcecomps = 1
  integer(ik),parameter :: num_networks_per_atom = 3

  integer(ik) :: num_types = 0, num_pairs = 0

  integer,parameter :: NMINCELL_FNN=1

  integer,allocatable :: pair_types(:,:)

  ! timing for 1-feature calc & 2-force inference
  real(rk),save,private :: tstart(0:3)=0.0, tfinish(0:3)=0.0

contains

!------------------------------------------------------------------------------
subroutine save_features_and_terminate(num_atoms, atype, pos, f, fp, filename) 
!------------------------------------------------------------------------------
implicit none
integer,intent(in) :: num_atoms 
real(8),intent(in),allocatable :: atype(:), pos(:,:), f(:,:)
type(fnn_param),intent(in) :: fp
character(len=:),allocatable,intent(in) :: filename

integer(ik) :: i, ity, gid, iunit, num_models
integer(ik),allocatable :: num_features(:)
integer(ik) :: rad_types, ang_types
integer(ik) :: num_features_rad, num_features_ang, num_features_total

open(newunit=iunit,file=filename//"_feature.xyz",form='formatted')
write(iunit,'(i6)') num_atoms
write(iunit,'(6f12.5)') lata,latb,latc,lalpha,lbeta,lgamma
do i=1, num_atoms
   ity = nint(atype(i))
   write(iunit,fmt='(a,6f12.5,3x)') fp%models(ity)%element,pos(i,1:3),f(i,1:3)
enddo
close(iunit)

num_models = size(fp%models)
allocate(num_features(num_models))

open(newunit=iunit,file=filename//"_feature.bin",access='stream',form='unformatted')
write(iunit) num_atoms, num_models

do ity = 1, num_models

   num_features(ity) = size(fp%models(ity)%features,dim=2)

   rad_types = size(fp%models(ity)%rad)
   ang_types = size(fp%models(ity)%ang)
   num_features_rad = fp%models(ity)%num_features(1)
   num_features_ang = fp%models(ity)%num_features(2)
   num_features_total = sum(fp%models(ity)%num_features)

   if(num_features(ity) /= num_features_total) stop

   write(iunit) rad_types, ang_types
   write(iunit) num_features_rad, num_features_ang, num_features_total

   do i=1,rad_types; write(iunit) fp%models(ity)%feature_ptr_rad(i); enddo
   do i=1,ang_types; write(iunit) fp%models(ity)%feature_ptr_ang(i); enddo

enddo

do i=1, num_atoms
   ity = nint(atype(i))
   gid = l2g(atype(i))
   write(iunit) gid, ity
   write(iunit) fp%models(ity)%features(1,1:num_features(ity),i)
   write(iunit) fp%models(ity)%features(2,1:num_features(ity),i)
   write(iunit) fp%models(ity)%features(3,1:num_features(ity),i)
enddo
close(iunit)

if(myid==0) then
  print'(a)',repeat('=',80)
  print'(a)','INFO: saved feature vectors in '//filename//'_feature.bin'
  print'(a)',repeat('=',80)
endif
call MPI_FINALIZE(ierr)
stop 

end subroutine


!------------------------------------------------------------------------------
subroutine get_force_fnn(ff, num_atoms, atype, pos, f, q)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff
type(fnn_param),pointer :: fp => null()
integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

integer(ik) :: i, ity, nn, nl, ncol, nrow

real(rk),allocatable :: x(:),y(:)

! not sure if this is the best way, but binding force_field_class to fnn_parm
select type(ff); type is (fnn_param) 
   fp => ff
end select

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos, ipos=ipos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

call get_features_fnn(num_atoms, atype, pos, fp) 

if(isRunFromXYZ) &
  call save_features_and_terminate(num_atoms,atype,pos,f,fp,RunFromXYZPath) 

call cpu_time(tstart(1))
do i=1, num_atoms 
  
   ity = atype(i)

   do nn=1, num_networks_per_atom  ! fx,fy,fz loop
 
     y = fp%models(ity)%features(nn,1:fp%models(ity)%layersize(1),i)
 
     associate(m=>fp%models(ity), n=>fp%models(ity)%networks(nn)) 
        do nl=1, size(n%dims)-1
           x = matmul(n%layers(nl)%w,y) + n%layers(nl)%b

           !print'(a,3i4,4i6)','i,ity,nn,shape(w),shape(b),shape(x): ', &
           !        i,ity,nn,shape(n%layers(nl)%w),shape(n%layers(nl)%b),shape(x)

           y = max(x,0.0) ! relu
        enddo 
     end associate
 
     f(i,nn) = x(1) ! update force
   enddo

! TODO: obtained force is scaled by the scaling_factor assuming that the trained
! weight matrix & biased are also scaled.  better to have a check mechanism on their consistency. 
   f(i,1:3)=f(i,1:3)/fp%models(ity)%scaling_factor*Eev_kcal

enddo
call cpu_time(tfinish(1))

#ifdef FORCEDUMP
do i=1, num_atoms
   ity=nint(atype(i))
   print'(2i6,6f12.5)',i,ity,pos(i,1:3),f(i,1:3)*fp%models(ity)%scaling_factor/Eev_kcal
enddo
stop 
#endif

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

call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)

call cpu_time(cpu0)

!--- set force model
do nstep=0, num_mdsteps-1

  if(mod(nstep,pstep)==0) call print_e_fnn(atype, v, q)

  call cpu_time(tstart(0))

  if(mod(nstep,fstep)==0) &
        call OUTPUT(atype, pos, v, q, GetFileNameBase(DataDir,current_step+nstep))

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

   if(msd_data%is_msd .and. mod(nstep,pstep)==0) & 
      call msd_data%measure(NATOMS, atype, pos, ipos)

!--- total force may not be zero with FNN. fix linear momentum every pstep.
   if(mod(nstep,pstep)==0) call linear_momentum(atype, v)

!--- update velocity & position
   call vkick(1.d0, atype, v, f)

   pos(1:natoms,1:3)=pos(1:natoms,1:3)+dt*v(1:natoms,1:3)

!--- migrate atoms after positions are updated
   call COPYATOMS(imode=MODE_MOVE_FNN,dr=[0.d0, 0.d0, 0.d0],atype=atype,pos=pos, &
                  v=v,f=f,q=q,ipos=ipos)

   call cpu_time(cpu1)
   call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)
   call cpu_time(cpu2)
   comp = comp + (cpu2-cpu1)

!--- update velocity
   call vkick(1.d0, atype, v, f)

   call cpu_time(tfinish(0))

enddo

!--- save the final configurations
call OUTPUT(atype, pos, v, q,  GetFileNameBase(DataDir,current_step+nstep))

!--- update rxff.bin in working directory for continuation run
call WriteBIN(atype, pos, v, q, GetFileNameBase(DataDir,-1))

if(msd_data%is_msd) call msd_data%save()


call cpu_time(cpu2)
if(myid==0) print'(a,2f12.5)','comp, total (sec): ', comp, cpu2-cpu0

return
end subroutine

!------------------------------------------------------------------------------
subroutine get_features_fnn(num_atoms, atype, pos, fp)
!------------------------------------------------------------------------------
integer,intent(in) :: num_atoms
real(8),intent(in),allocatable :: atype(:), pos(:,:)
type(fnn_param),intent(in out) :: fp

real(rk) :: rr(3), rr2, rij, dsum 
integer(ik) :: i, j, k, i1, j1, k1, l1, l2, l3, l4, ii
integer(ik) :: c1,c2,c3,ic(3),c4,c5,c6,ity,jty,kty,inxn

integer(ik) :: nnbr, lnbr(MAXNEIGHBS)
integer(ik) :: idx, idx_stride, l1_stride, l2_stride, l3_stride, l4_stride
real(rk) :: r_ij(0:3), r_kj(0:3), r_ij_norm(3), r_kj_norm(3), eta_ij, eta_kj, fc_ij, fc_kj, rij_mu, rkj_mu
real(rk) :: cos_ijk, lambda_ijk, rijk_inv, zeta_G3a, zeta_G3b, zeta_G3b_0, zeta_const

real(rk) :: G3_mu_eta, G3a_xyz(3), G3a_c1, G3a_c2, G3a, G3b_xyz(3), G3b

real(rk) :: fdr, frac
integer :: size_etamu, ia

nbrlist(:,0) = 0

call cpu_time(tstart(2))
do ity = 1, size(fp%models)
   fp%models(ity)%features(:,:,1:num_atoms) = 0.0
enddo

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

             ! lookup table. should be ii == jty
             ii = fp%models(ity)%map_rad(jty) 

             associate( rad => fp%models(ity)%rad(ii), &
                        ptr => fp%models(ity)%feature_ptr_rad(ii), &
                        max_ang_rc => fp%models(ity)%max_ang_rc, &
                        feats => fp%models(ity)%features ) 

                if(rij<rad%rc) then
  
                  if(rij<max_ang_rc) then
                    nbrlist(i, 0) = nbrlist(i, 0) + 1
                    nbrlist(i, nbrlist(i, 0)) = j
                  endif

                  rr(1:3) = rr(1:3)/rij

                  size_etamu = size(rad%mu)*size(rad%eta)

#ifdef TABLE
                  fdr = rij/rad%ftab_dr
                  idx = int(fdr)
                  frac = fdr - idx

                  do ia=1,3
                    feats(ia,ptr:ptr+size_etamu-1,i) = &
                    feats(ia,ptr:ptr+size_etamu-1,i) + rad%ftab(idx,1:size_etamu)*rr(ia)
                  enddo
#else
                  fc_ij = 0.5*( 1.0 + cos(pi*rij/rad%rdamp) ) 

                  do l1 = 1, size(rad%mu)
                     rij_mu = rij - rad%mu(l1)

                     do l2 = 1, size(rad%eta)
                        eta_ij = exp( -rad%eta(l2) * rij_mu * rij_mu )
                        idx = ptr + (l1-1)*size(rad%eta) + (l2-1)
                        feats(1:3,idx,i) = feats(1:3,idx,i) + eta_ij*fc_ij*rr(1:3)
                     enddo
  
                  enddo
#endif

                endif

             end associate

           endif

           j=llist(j)
        enddo
     enddo; enddo; enddo

     i = llist(i)
  enddo
enddo; enddo; enddo
!$omp end parallel do 
call cpu_time(tfinish(2))

call cpu_time(tstart(3))
do j=1, num_atoms

   jty = nint(atype(j))

   do i1=1, nbrlist(j,0)-1

      i = nbrlist(j,i1)
      ity = nint(atype(i))

      r_ij(1:3) = pos(j,1:3) - pos(i,1:3)
      r_ij(0) = sqrt( sum(r_ij(1:3)*r_ij(1:3)) )
      r_ij_norm(1:3) = r_ij(1:3)/r_ij(0)


      do k1=i1+1, nbrlist(j,0)

         k = nbrlist(j,k1)
         kty = nint(atype(k))

         if(i==k) cycle

         r_kj(1:3) = pos(j,1:3) - pos(k,1:3)
         r_kj(0) = sqrt( sum(r_kj(1:3)*r_kj(1:3)) )
         r_kj_norm(1:3) = r_kj(1:3)/r_kj(0)

         rijk_inv = 1.0/(r_ij(0) * r_kj(0))

         ii = fp%models(jty)%map_ang(ity,kty) 
         if(ii<1) cycle ! no triplet param exists if map_ang < 1. go to next iteration.

         associate(ang => fp%models(jty)%ang(ii), &
                   ptr_ang => fp%models(jty)%feature_ptr_ang(ii), &
                   feats => fp%models(jty)%features )

           ! the neighborlist is constructed with the max_ang_rc. 
           ! need to check the angular cutoff again.
           if(r_ij(0)>ang%rc .or. r_kj(0)>ang%rc) cycle

           l1_stride = NUM_ANG_MODES*size(ang%eta)*size(ang%zeta)*size(ang%lambda)
           l2_stride = NUM_ANG_MODES*size(ang%eta)*size(ang%zeta)
           l3_stride = NUM_ANG_MODES*size(ang%eta)
           l4_stride = NUM_ANG_MODES
  
           fc_ij = 0.5*( 1.0 + cos(pi*r_ij(0)/ang%rdamp) ) 
           fc_kj = 0.5*( 1.0 + cos(pi*r_kj(0)/ang%rdamp) ) 
  
           cos_ijk = sum( r_ij(1:3)*r_kj(1:3) ) * rijk_inv
  
           G3a_c1 = r_ij(0)-r_kj(0)*cos_ijk
           G3a_c2 = r_kj(0)-r_ij(0)*cos_ijk
  
           G3a_xyz(1:3) = (r_ij_norm(1:3)*G3a_c1 + r_kj_norm(1:3)*G3a_c2)*rijk_inv*fc_ij*fc_kj
           G3b_xyz(1:3) = (r_ij_norm(1:3) + r_kj_norm(1:3))*fc_ij*fc_kj
  
  ! l1: mu, l2:lambda, l3:zeta, l4:eta
           do l1=1, size(ang%mu)
  
              rij_mu = r_ij(0) - ang%mu(l1)
              rkj_mu = r_kj(0) - ang%mu(l1)
  
              do l2=1, size(ang%lambda)
  
                 lambda_ijk = 1.0 + ang%lambda(l2)*cos_ijk 
  
                 do l3=1, size(ang%zeta)
  
                    zeta_const = 2.0**(1.0-ang%zeta(l3))
                    zeta_G3a = ang%zeta(l3) * ang%lambda(l2) * zeta_const * (lambda_ijk**(ang%zeta(l3)-1))
                    zeta_G3b_0 = zeta_const * (lambda_ijk**ang%zeta(l3))
  
                    do l4=1, size(ang%eta)
  
                       zeta_G3b = - 2.0*ang%eta(l4) * zeta_G3b_0
  
                       eta_ij = exp( -ang%eta(l4) * rij_mu * rij_mu )
                       eta_kj = exp( -ang%eta(l4) * rkj_mu * rkj_mu )
  
                       G3_mu_eta = eta_ij*eta_kj

                       idx = ptr_ang + &
                           (l1-1)*l1_stride + (l2-1)*l2_stride + (l3-1)*l3_stride + (l4-1)*l4_stride 
  
                       !feats(1:3,idx,j) = feats(1:3,idx,j) + &
                       !    G3a_xyz(1:3)*zeta_G3a*G3_mu_eta + G3b_xyz(1:3)*zeta_G3b*G3_mu_eta
                       feats(1:3,idx,j) = feats(1:3,idx,j) + G3a_xyz(1:3)*zeta_G3a*G3_mu_eta
                       feats(1:3,idx+1,j) = feats(1:3,idx+1,j) + G3b_xyz(1:3)*zeta_G3b*G3_mu_eta
  
           enddo; enddo; enddo; enddo

         end associate
          
      enddo
   enddo

enddo

do i=1, num_atoms
   ity = nint(atype(i))
   do j = 1, 3 ! xyz-loop
      fp%models(ity)%features(j,:,i) = &
     (fp%models(ity)%features(j,:,i) - fp%models(ity)%fstat(j)%mean(:))/fp%models(ity)%fstat(j)%stddev(:)
   enddo 
enddo
call cpu_time(tfinish(3))

return
end subroutine

!------------------------------------------------------------------------------
subroutine mean_stddev_loader(mean, stddev, feature_size, path, suffix, verbose) 
!------------------------------------------------------------------------------
real(rk),allocatable,intent(in out) :: mean(:), stddev(:)
integer(ik),intent(in) :: feature_size
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix !<- ['x','y','z']
logical,optional :: verbose

integer(ik) :: funit_m, funit_s
logical :: exist_m, exist_s
character(len=:),allocatable :: filename_m, filename_s

filename_m = path//'feature_mean_'//int_to_str(feature_size)//'.'//suffix
filename_s = path//'feature_stddev_'//int_to_str(feature_size)//'.'//suffix

if(.not.allocated(mean)) allocate(mean(feature_size))
if(.not.allocated(stddev)) allocate(stddev(feature_size))

inquire(file=filename_m, exist=exist_m)
inquire(file=filename_s, exist=exist_s)

! mean and stddev must exist, otherwise no feature vector standardization.
if(exist_m .and. exist_s) then

  open(newunit=funit_m, file=filename_m, access='stream', form='formatted', status='old')
  if(myid==0) read(funit_m,*) mean
  call MPI_BCAST(mean, size(mean), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
  close(funit_m)

  open(newunit=funit_s, file=filename_s, access='stream', form='formatted', status='old')
  if(myid==0) read(funit_s,*) stddev
  call MPI_BCAST(stddev, size(stddev), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
  close(funit_s)

  if(present(verbose) .and. verbose) &
     write(*,fmt='(a15,a30,a15,a30)') 'mean: ', filename_m, ', stddev: ', filename_s

else

  if(present(verbose) .and. verbose) then
     print'(a)', repeat('- ',40)
     if(.not. exist_m) print'(a)', 'WARNING: missing '//filename_m
     if(.not. exist_s) print'(a)', 'WARNING: missing '//filename_s
     print'(a)', 'WARNING: incomplete mean & stddev data. continue with mean=0.0 & stddev=1.0'
     print'(a)', repeat('- ',40)
  endif

  mean=0.0
  stddev=1.0
endif

end subroutine

!------------------------------------------------------------------------------
type(network) function network_ctor(dims, path, suffix, verbose) result(net)
!------------------------------------------------------------------------------

integer(ik),intent(in) :: dims(:)
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix !<- ['x','y','z']
logical,optional :: verbose

character(len=:),allocatable :: filename_b, filename_w
character(len=:),allocatable :: arow, acol, alayer

integer(ik) :: i, nrow, ncol, fileunit
integer(ik) :: num_layers
logical :: has_model

net%dims = dims
num_layers = size(net%dims)

allocate(net%layers(num_layers))

do i=1, num_layers-1
  nrow = dims(i)
  ncol = dims(i+1)

  allocate(net%layers(i)%b(ncol))
  allocate(net%layers(i)%w(ncol,nrow))
  !print*,'i,nrow,ncol: ', i,nrow,ncol 

  alayer = int_to_str(i)
  arow = int_to_str(nrow)
  acol = int_to_str(ncol)

  filename_b = trim(path)//'b_'//alayer//'_'//acol//'.'//suffix
  inquire(file=filename_b, exist=has_model)
  if(has_model) then
    open(newunit=fileunit, file=filename_b, access='stream', form='formatted', status='old')
    if(myid==0) read(fileunit,*) net%layers(i)%b
    call MPI_BCAST(net%layers(i)%b, size(net%layers(i)%b), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
    close(fileunit)
  else
    if(myid==0) then
      print'(a)',repeat('-',80)
      print'(a)', 'ERROR: '//filename_b//' does not exist. continue with zero-valued bias.' 
      print'(a)',repeat('-',80)
    endif
    net%layers(i)%b=0.0
  endif

  filename_w = trim(path)//'w_'//alayer//'_'//arow//'_'//acol//'.'//suffix
  inquire(file=filename_w, exist=has_model)
  if(has_model) then
     open(newunit=fileunit, file=filename_w, access='stream', form='formatted', status='old')
     if(myid==0) read(fileunit,*) net%layers(i)%w
     call MPI_BCAST(net%layers(i)%w, size(net%layers(i)%w), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
     close(fileunit)
  else
     
    if(myid==0) then
      print'(a)',repeat('-',80)
      print'(a)', 'ERROR: '//filename_w//' does not exist. continue with zero-valued weight.' 
      print'(a)',repeat('-',80)
    endif
    net%layers(i)%w=0.0
  endif

  if(present(verbose) .and. verbose) &
     write(*, fmt='(a30,2i6,a30,i6)') &
        'w: '//filename_w, shape(net%layers(i)%w), ' b: '//filename_b, size(net%layers(i)%b)

enddo

end function

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
real(8) :: tt=0.d0

ke=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE (MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
ke = ke/GNATOMS
tt = ke*UTEMP
GKE = ke ! FIXME for ctmp = (treq*UTEMP0)/( GKE*UTEMP )

if(myid==0) then
   
   cstep = nstep + current_step 

   write(6,'(a,i9,es13.5,f10.3,3x,4f10.5)') &
        'MDstep,KE,T(K)   onestep(sec),finf,feat2b,feat3b: ', cstep, ke, tt, tfinish(0:3)-tstart(0:3)
endif

end subroutine

end module

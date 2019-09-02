module stats_mod

use utils
use lists_mod, only : nbplist, linkedlist, &
    get_mesh_for_nonbonding_list, get_nonbonding_list, &
    nblcsize, nbheader, nbllist, nbnacell

implicit none

type gr_type
  real(8),allocatable :: dat(:)
  integer,allocatable :: sample(:)

  contains
    procedure :: measure => measure_gr 
    procedure :: save => save_gr
end type

type ba_type
  real(8),allocatable :: dat(:,:)
  integer,allocatable :: sample(:,:)

  contains 
    procedure :: measure => measure_ba
    procedure :: save => save_ba
end type

type stats_type

  real(8) :: rcut
  logical :: run_stats 

  type(gr_type) :: gr
  type(ba_type),allocatable :: ba(:)

  contains 
    procedure :: measure => measure_stats
    procedure :: save => save_stats
end type

type(stats_type) :: stats

contains

function stats_type_ctor(rcut,atom_names) result(c)
  real(8),intent(in) :: rcut
  character(2),allocatable,intent(in) :: atom_names(:)

  type(stats_type) :: c
  integer :: i,na,idx

  if(find_cmdline_argc('--stats',idx).or.find_cmdline_argc('-stats',idx)) then
    c%run_stats = .true.
  else
    c%run_stats = .false.
    return
  endif

  c%rcut=rcut

  na = size(atom_names)

  allocate(c%gr%dat(na), c%gr%sample(na), c%ba(na))
  c%gr%dat=0.d0; c%gr%sample=0 
  do i=1, na
    allocate(c%ba(i)%dat(na,na), c%ba(i)%sample(na,na))
    c%ba(i)%dat=0.d0; c%ba(i)%sample=0
  enddo

  call get_mesh_for_nonbonding_list(rcut)

  return
end function

subroutine measure_gr(this)
  class(gr_type) :: this
end subroutine

subroutine save_gr(this, filename)
  class(gr_type) :: this
  character(len=:),allocatable :: filename
end subroutine

subroutine measure_ba(this)
  class(ba_type) :: this
end subroutine

subroutine save_ba(this)
  class(ba_type) :: this
end subroutine

subroutine measure_stats(this,num_atoms,pos,atype)
  class(stats_type) :: this
  integer,intent(in) :: num_atoms
  real(8),intent(in),allocatable :: atype(:), pos(:,:)

  integer :: i

  call linkedlist(atype,pos,nblcsize,nbheader,nbllist,nbnacell)
  call get_nonbonding_list(pos,this%rcut) 

  do i=1, num_atoms
     print*,i,pos(i,1:3),nbplist(0,i)
  enddo

  stop 'foo'
end subroutine

subroutine save_stats(this)
  class(stats_type) :: this
end subroutine

end module

module stats_mod

use utils
use base, only : myid
use lists_mod, only : nbplist, linkedlist, &
    get_mesh_for_nonbonding_list, get_nonbonding_list, &
    nblcsize, nbheader, nbllist, nbnacell

implicit none

type gr_type
  real(8),allocatable :: dat(:,:)
  integer,allocatable :: samples(:,:)

  contains
    procedure :: measure => measure_gr 
    procedure :: save => save_gr
end type

type ba_type
  real(8),allocatable :: dat(:,:,:)
  integer,allocatable :: samples(:,:,:)

  contains 
    procedure :: measure => measure_ba
    procedure :: save => save_ba
end type

real(8),parameter :: GR_RESOLUTION = 0.05d0 ! [A]
real(8),parameter :: BA_RESOLUTION = 1.d0  ! degree

type stats_type

  character(2),allocatable :: atom_names(:)
  real(8),allocatable :: concentration(:)
  integer,allocatable :: num_atoms_per_type(:)

  real(8) :: rcut, rho, dgr,dgri, dba,dbai
  integer :: gr_size, ba_size
  logical :: run_stats 

  type(gr_type),allocatable :: gr(:)
  type(ba_type),allocatable :: ba(:)

  contains 
    procedure :: measure => measure_stats
    procedure :: save => save_stats
end type

type(stats_type) :: stats

contains

function stats_type_ctor(rcut,atom_names,mdbox,natoms_per_type) result(c)
  real(8),intent(in) :: rcut, mdbox
  character(2),allocatable,intent(in) :: atom_names(:)
  integer(8),intent(in) :: natoms_per_type(*)

  type(stats_type) :: c
  integer :: i,na,idx
  integer(8) :: total_natoms

  if(find_cmdline_argc('--stats',idx).or.find_cmdline_argc('-stats',idx)) then
    c%run_stats = .true.
  else
    c%run_stats = .false.
    return
  endif

  c%rcut = rcut

  c%dgr = GR_RESOLUTION
  c%dgri = 1.d0/c%dgr
  c%gr_size= int(c%rcut*c%dgri)

  c%dba = BA_RESOLUTION
  c%dbai = 1.d0/c%dba
  c%ba_size= int(180*c%dbai)

  c%atom_names = atom_names
  na = size(c%atom_names)

  allocate(c%gr(na), c%ba(na), c%concentration(na), c%num_atoms_per_type(na))
  do i = 1, na
    allocate(c%gr(i)%dat(na,0:c%gr_size),    c%gr(i)%samples(na,0:c%gr_size))
    allocate(c%ba(i)%dat(na,na,0:c%ba_size), c%ba(i)%samples(na,na,0:c%ba_size))

    c%gr(i)%dat=0.d0; c%gr(i)%samples=0 
    c%ba(i)%dat=0.d0; c%ba(i)%samples=0
  enddo

  
  total_natoms = sum(natoms_per_type(1:na))
  do i = 1, na
    c%num_atoms_per_type(i) = natoms_per_type(i)
    c%concentration(i) = dble(natoms_per_type(i))/total_natoms
  enddo
  c%rho = total_natoms/mdbox

  call get_mesh_for_nonbonding_list(rcut)

  if(myid==0) then
    write(6, fmt='(a)') repeat('-',60)
    write(6, fmt='(a,f8.3,1x,f8.3, i12)') 'rcut, rho, total: ', c%rcut, c%rho, total_natoms
    write(6, fmt='(a)',advance='no') 'concentration: ' 
    do i = 1, na
      write(6, fmt='(f8.3,1x)',advance='no') c%concentration(i)
    enddo
    write(6,*)
    write(6, fmt='(a,2f8.3,i6)') 'dgr, dgri, gr_size: ', c%dgr,c%dgri,c%gr_size
    write(6, fmt='(a,2f8.3,i6)') 'dba, dbai, ba_size: ', c%dba,c%dbai,c%ba_size
    write(6, fmt='(a,2f8.3,i6)') 'dba, dbai, ba_size: ', c%dba,c%dbai,c%ba_size
    write(6, fmt='(a)') repeat('-',60)
  endif

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

  integer :: i, j, j1, ity, jty, idx
  real(8) :: dr(0:3), rr, prefactor, prefactor2

  call linkedlist(atype,pos,nblcsize,nbheader,nbllist,nbnacell)
  call get_nonbonding_list(pos,this%rcut) 

  do i=1, num_atoms
     ity = nint(atype(i))
     do j1 = 1, nbplist(0,i)
        j = nbplist(j1,i)
        jty = nint(atype(j))
        dr(1:3) = pos(j,1:3) - pos(i,1:3)
        dr(0) = sqrt(sum(dr(1:3)*dr(1:3)))
        idx = dr(0)*this%dgri
        this%gr(ity)%samples(jty,idx) = this%gr(ity)%samples(jty,idx) + 1 
     enddo
  enddo

  write(6,fmt='(a,1x)', advance='no') ' distance(A)'
  do ity = 1, size(this%atom_names)
    do jty = 1, size(this%atom_names)
      write(6,fmt='(a8,1x)', advance='no') &
        trim(this%atom_names(ity))//'-'//trim(this%atom_names(jty))
    enddo
  enddo
  write(6,*)

  do idx = 1, this%gr_size-1

    rr = idx*this%dgr
    write(6,fmt='(f8.3,1x)', advance='no') rr
    prefactor = 4.d0*pi*rr*rr*this%rho/this%dgri

    do ity = 1, size(this%atom_names)
      do jty = 1, size(this%atom_names)

        prefactor2 = this%concentration(jty)*this%num_atoms_per_type(ity)

        this%gr(ity)%dat(jty,idx) = this%gr(ity)%samples(jty,idx)/(prefactor*prefactor2)
        write(6,fmt='(f12.5,1x)', advance='no') this%gr(ity)%dat(jty,idx)

      enddo
    enddo
    write(6,*)

  enddo

  stop 'foo'
end subroutine

subroutine save_stats(this)
  class(stats_type) :: this
end subroutine

end module

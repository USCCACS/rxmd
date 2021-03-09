



module io

  implicit none
  private

  public :: io_adjustl, io_adjustl_sub,  &
            io_center,      &
            io_lower!,       &

integer, public :: nfile(2), myid, nodes

  !--------------------------------------------------------------------!
  !                       io_adjustl (interface)                       !
  !                                                                    !
  ! Use io_adjustl to adjust not only strings but also different types !
  ! for convenient output.                                             !
  !                                                                    !
  ! Available implementations:                                         !
  !   adjustl_i(int)         -- integer                                !
  !   adjustl_d(dp, digits)  -- double precision (digits is optional)  !
  !--------------------------------------------------------------------!

  interface io_adjustl
     module procedure adjustl_i, adjustl_d
  end interface

  interface io_adjustl_sub
     module procedure adjustl_sub_i, adjustl_sub_d
  end interface

contains

  !--------------------------------------------------------------------!
  !                      adjustl (implementation)                      !
  !                                                                    !
  !          Please use the module interface `io_adjustl()'.           !
  !--------------------------------------------------------------------!

  function adjustl_i(int) result(str)

    implicit none

    integer, intent(in) :: int
    character(len=50)   :: str

    write(str, *) int
    str = trim(adjustl(str))

  end function adjustl_i

  !--------------------------------------------------------------------!

  subroutine adjustl_sub_i(int,str)

    implicit none

    integer,           intent(in)  :: int
    character(len=50), intent(out) :: str

    write(str, *) int
    str = trim(adjustl(str))

  end subroutine adjustl_sub_i

  !--------------------------------------------------------------------!

  function adjustl_d(dp, digits) result(str)

    implicit none

    double precision,            intent(in) :: dp
    integer,           optional, intent(in) :: digits

    character(len=50)            :: frmt, str

    if (present(digits)) then
       write(frmt, *) digits
       frmt = '(F50.' // trim(adjustl(frmt)) // ')'
    else
       frmt = '(F50.2)'
    end if

    write(str, frmt) dp
    str = trim(adjustl(str))

  end function adjustl_d

  !--------------------------------------------------------------------!

  subroutine adjustl_sub_d(dp, str, digits)

    implicit none

    double precision,            intent(in)  :: dp
    character(len=50),           intent(out) :: str
    integer,           optional, intent(in)  :: digits

    character(len=50)            :: frmt

    if (present(digits)) then
       write(frmt, *) digits
       frmt = '(F50.' // trim(adjustl(frmt)) // ')'
    else
       frmt = '(F50.2)'
    end if

    write(str, frmt) dp
    str = trim(adjustl(str))

  end subroutine adjustl_sub_d

  !--------------------------------------------------------------------!
  !                  Lower and upper case conversion                   !
  !--------------------------------------------------------------------!

  function io_lower(str_in) result(str_out)

    implicit none

    character(len=*), intent(in) :: str_in
    character(len=len(str_in))   :: str_out

    integer, parameter :: ilowerA = ichar('a')
    integer, parameter :: iupperA = ichar('A')
    integer, parameter :: iupperZ = ichar('Z')

    integer :: i, ichr, nchr, iconv

    iconv = ilowerA - iupperA

    nchr = len(str_in)
    do i = 1, nchr
       ichr = ichar(str_in(i:i))
       if ((ichr >= iupperA) .and. (ichr <= iupperZ)) then
          str_out(i:i) = char(ichr + iconv)
       else
          str_out(i:i) = str_in(i:i)
       end if
    end do

  end function io_lower

  !--------------------------------------------------------------------!
  !                    io_center - center a string                     !
  !--------------------------------------------------------------------!

  function io_center(str, n, str2) result(l1)  !) result(str2)

    implicit none

    character(len=*), intent(in) :: str
    integer,          intent(in) :: n
    character(len=n)             :: str2

    integer :: l1

    str2 = ""
    l1 = len_trim(adjustl(str))
    if (l1 > n) return

    l1 = (n - l1)/2

    str2 = repeat(' ',l1) // trim(adjustl(str))

  end function io_center

  !--------------------------------------------------------------------!

end module io




module aeio

  use io, only: & !io_lower,    &
                io_center,  &
&               nfile, myid, nodes

  implicit none
  private
  save

  public :: & !aeio_readline,             &
            aeio_header,               &
            aeio_assert_file_exists!,   &


  !---------------------------- constants -----------------------------!
  ! LINELEN    max. length of a line read from an input file           !
  ! PATHLEN    max. number of characters available for a file path     !
  ! TYPELEN    max. number of characters for atmoic species names      !
  ! STDIN      file unit of standard in                                !
  ! STDOUT     file unit of standard out                               !
  ! STDERR     file unit of standard error                             !
  !--------------------------------------------------------------------!

  integer, parameter, public :: LINELEN = 1024
  integer, parameter, public :: PATHLEN = 80
  integer, parameter, public :: TYPELEN = 2
  integer, parameter, public :: STDIN   = 5
  integer, parameter, public :: STDOUT  = 6
  integer, parameter, public :: STDERR  = 0

contains

  !--------------------------------------------------------------------!
  !           write a centered header (for formatted output)           !
  !--------------------------------------------------------------------!

  subroutine aeio_header(str, char, unit)

    character(len=*),    intent(in) :: str
    character, optional, intent(in) :: char
    integer,   optional, intent(in) :: unit

    character :: c
    character(len=70) :: str2
    integer :: ii

    if (present(char)) then
       c = char
    else
       c = '-'
    end if

    if (present(unit)) then
       write(unit,*) repeat(c,70)
!       write(unit,*) io_center(trim(str),70)
       ii = io_center(trim(str),70,str2)
       write(unit,*) str2
       write(unit,*) repeat(c,70)
    else
       write(nfile(1),*) repeat(c,70)
       write(nfile(2),*) repeat(c,70)
!       write(*,*) io_center(trim(str),70)
       ii = io_center(trim(str),70,str2)
       write(nfile(1),*) str2
       write(nfile(2),*) str2
       write(nfile(1),*) repeat(c,70)
       write(nfile(2),*) repeat(c,70)
    end if

  end subroutine aeio_header

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine aeio_assert_file_exists(file)
    implicit none
    character(len=*), intent(in) :: file
    logical :: fexists
    inquire(file=trim(adjustl(file)), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: file not found: ", trim(adjustl(file))
       stop
    end if
  end subroutine aeio_assert_file_exists

end module aeio




module geometry

  implicit none
  private
  save

  public  :: geo_itype_of_name!, &

contains

  !====================================================================!
  !                                                                    !
  !                        auxiliary procedures                        !
  !                                                                    !
  !====================================================================!

  !--------------------------------------------------------------------!
  !                retrieve type number from type name                 !
  !--------------------------------------------------------------------!

  function geo_itype_of_name(name, typeName) result(itype)

    implicit none

    character(len=*),               intent(in) :: name
    character(len=*), dimension(:), intent(in) :: typeName
    integer                                    :: itype

    integer :: itype1, ntypes

    ntypes = size(typeName)

    itype = 0
    do itype1 = 1, nTypes
       if (trim(typeName(itype1)) == trim(name)) then
          itype = itype1
          exit
       end if
    end do

  end function geo_itype_of_name

end module geometry




module feedforward

  use io, only: nfile, myid, nodes

  implicit none
  private
  save

  public  :: load_Network,           &
             load_Network_ASCII,     &
             ff_memsize,             &
             ff_print_info,          &
             ff_eval,                &
             ff_deriv!,               &

  !--------------------------------------------------------------------!
  !                            Network type                            !
  !--------------------------------------------------------------------!

  type, public :: Network

     !-----------------------------------------------------------------!
     ! init         : .true., if the instance was initialized          !
     ! memsize      : dynamically allocated words for this instance    !
     ! nlayers      : number of layers in the network, inluding input  !
     !                and output layer                                 !
     ! nnodes(i)    : number of nodes (without bias) in the i-th layer !
     ! nnodes_max   : max. number of nodes in a layer of this network  !
     ! f_a(i)       : activation function type for the i-th layer      !
     ! W(i)         : i-th weight of the graph edges (including bias)  !
     ! iw(i)        : index of the last weight of the i-th layer       !
     !                e.g.: W(iw(2)) --> last weight in the 2nd layer  !
     ! Wsize        : total number of weights; Wsize = size(W)         !
     ! value(i)     : if (evaluated) --> value of the i-th neuron      !
     !                if not         --> undefined                     !
     ! deriv(i)     : if (evaluated) --> derivative of the i-th neuron !
     !                if not         --> undefined                     !
     ! iv(i)        : index of the last node in the i-th layer         !
     ! nvalues      : total number of nodes/neurons in the network     !
     ! evaluated    : .true., if the node values and derivatives have  !
     !                been evaluated                                   !
     ! D(i)         : matrices D_ij of partial derivatives wrt. nodes  !
     ! derivatives  : .true., if the D_ij matrices have been evaluated !
     ! work...      : scratch memory for the computations              !
     !-----------------------------------------------------------------!

     logical                                       :: init = .false.
     integer                                       :: memsize

     integer                                       :: nlayers
     integer,          dimension(:),   allocatable :: nnodes
     integer                                       :: nnodes_max

     integer,          dimension(:),   allocatable :: f_a

     double precision, dimension(:),   allocatable :: W
     integer,          dimension(:),   allocatable :: iw
     integer                                       :: Wsize

     double precision, dimension(:),   allocatable :: value
     double precision, dimension(:),   allocatable :: deriv
     integer,          dimension(:),   allocatable :: iv
     integer                                       :: nvalues
     logical                                       :: evaluated

     double precision, dimension(:),   allocatable :: D
     logical                                       :: derivatives

     double precision, dimension(:),   allocatable :: work
     double precision, dimension(:,:), allocatable :: work2, work3, work4

  end type Network

contains !-------------------------------------------------------------!

  function load_Network(file, unit) result(net)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(Network)                          :: net

    integer                      :: u
    integer,           parameter :: u_sav = 41

    logical :: fexists

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = u_sav
       inquire(file=trim(file), exist=fexists)
       if (.not. fexists) then
          write(0,*) "Error: file not found in `load_Network': ", &
                     trim(file)
          stop
       end if
       open(u, file=trim(file), status='old', &
            form='unformatted', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `save_Network'."
       stop
    end if

    read(u) net%nlayers
    read(u) net%nnodes_max
    read(u) net%Wsize
    read(u) net%nvalues

    allocate(net%nnodes(net%nlayers),                          &
             net%f_a(2:net%nlayers),                           &
             net%iw(1:net%nlayers),                            &
             net%iv(1:net%nlayers),                            &
             net%W(net%Wsize),                                 &
             net%D(net%Wsize),                                 &
             net%value(net%nvalues),                           &
             net%deriv(net%nvalues),                           &
             net%work(net%nnodes_max+1),                       &
             net%work2(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work3(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work4(net%nnodes_max+1,(net%nnodes_max+1)**2) )

    read(u) net%nnodes(:)
    read(u) net%f_a(:)
    read(u) net%iw(:)
    read(u) net%iv(:)
    read(u) net%W(:)

    if (.not. present(unit)) close(u)

    net%init        = .true.
    net%evaluated   = .false.
    net%derivatives = .false.

    ! number of words allocated for this object are (at least):
    net%memsize = ff_memsize(net)

  end function load_Network

  !--------------------------------------------------------------------!

  function load_Network_ASCII(file, unit) result(net)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(Network)                          :: net

    integer                      :: u
    integer,           parameter :: u_sav = 41

    integer                     :: i
    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'
    character(len=*), parameter :: ifrmt = '(4(1x,I17))'
    logical                     :: fexists

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = u_sav
       inquire(file=trim(file), exist=fexists)
       if (.not. fexists) then
          write(0,*) "Error: file not found in `load_Network': ", &
                     trim(file)
          stop
       end if
       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `save_Network'."
       stop
    end if

    read(u,*) net%nlayers
    read(u,*) net%nnodes_max
    read(u,*) net%Wsize
    read(u,*) net%nvalues

    allocate(net%nnodes(net%nlayers),                          &
             net%f_a(2:net%nlayers),                           &
             net%iw(1:net%nlayers),                            &
             net%iv(1:net%nlayers),                            &
             net%W(net%Wsize),                                 &
             net%D(net%Wsize),                                 &
             net%value(net%nvalues),                           &
             net%deriv(net%nvalues),                           &
             net%work(net%nnodes_max+1),                       &
             net%work2(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work3(net%nnodes_max+1,net%nnodes_max+1),     &
             net%work4(net%nnodes_max+1,(net%nnodes_max+1)**2) )

    read(u,ifrmt) (net%nnodes(i), i=1,net%nlayers)
    read(u,ifrmt) (net%f_a(i), i=2,net%nlayers)
    read(u,ifrmt) (net%iw(i), i=1,net%nlayers)
    read(u,ifrmt) (net%iv(i), i=1,net%nlayers)
    read(u,dfrmt) (net%W(i), i=1,net%Wsize)

    if (.not. present(unit)) close(u)

    net%init        = .true.
    net%evaluated   = .false.
    net%derivatives = .false.

    ! number of words allocated for this object are (at least):
    net%memsize = ff_memsize(net)

  end function load_Network_ASCII

  !--------------------------------------------------------------------!
  !         size of the allocated memory of a Network instance         !
  !--------------------------------------------------------------------!

  function ff_memsize(net) result(memsize)

    implicit none

    type(Network), intent(in) :: net
    integer                   :: memsize

    if (.not. net%init) then
       memsize = 0
    else
       memsize = (4*net%nlayers - 1)           * 1 &
               + (2*net%Wsize + 2*net%nvalues) * 2 &
               + (net%nnodes_max + 1)          * 2 &
               + 2*(net%nnodes_max+1)**2       * 2 &
               + (net%nnodes_max+1)**3         * 2
    end if

  end function ff_memsize

  !--------------------------------------------------------------------!
  !             print information about a Network instance             !
  !--------------------------------------------------------------------!

  subroutine ff_print_info(net, verbose)

    implicit none

    type(Network),     intent(in) :: net
    logical, optional, intent(in) :: verbose

    character(len=100) :: fname
    integer            :: ilayer
    integer            :: iw1, iw2
    integer            :: nf

    if (.not. net%init) then
       write(*,*) 'The network is not initialized.'
       write(*,*)
       return
    end if

 do nf = 1, 2

    write(nfile(nf),'(1x,"Number of layers : ",I3)') net%nlayers
    write(nfile(nf),*)
    write(nfile(nf),'(1x,"Number of nodes (without bias) ")')
    write(nfile(nf),'(1x,"and activation type per layer :")')
    write(nfile(nf),*)
    write(nfile(nf),'(5x,I3," : ",I5)') 1, net%nnodes(1)
    do ilayer = 2, net%nlayers
       select case(net%f_a(ilayer))
       case(0)
          fname = 'linear function (linear)'
       case(1)
          fname = 'hyperbolic tangent (tanh)'
       case(2)
          fname = 'logistic function (sigmoid)'
       case(3)
          fname = 'scaled hyperbolic tangent (mtanh)'
       case(4)
          fname = 'scaled hyperbolic tangent + linear twisting (twist)'
       end select
       write(nfile(nf),'(5x,I3," : ",I5,2x,A)') ilayer, net%nnodes(ilayer), trim(fname)
    end do
    write(nfile(nf),*)

    write(nfile(nf),'(1x,"Dynamically allocated words : ",I10," (",F10.2," KB)")') &
         net%memsize, dble(net%memsize*8)/1024.0d0
    write(nfile(nf),*)

    write(nfile(nf),'(1x,"Total number of weights (incl. bias) : ",I8)') net%Wsize
    write(nfile(nf),*)

    verbose1 : if (present(verbose)) then
     if (verbose) then
       write(nfile(nf),'(1x,"Weight matrices:")')
       write(nfile(nf),*)

       iw1 = 1
       do ilayer = 1, net%nlayers-1
          iw2 = net%iw(ilayer+1)
          call print_mat(reshape(net%W(iw1:iw2), &
               (/net%nnodes(ilayer+1),net%nnodes(ilayer)+1/) ), nfile(nf))
          write(nfile(nf),*)
       end do
     end if
    end if verbose1

 end do

  end subroutine ff_print_info

  !--------------------------------------------------------------------!
  !                                                                    !
  !                 evaluation of the network function                 !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine ff_eval(net, nx, x, ny, y)

    implicit none

    type(Network),                   intent(inout) :: net
    integer,                         intent(in)    :: nx
    double precision, dimension(nx), intent(in)    :: x
    integer,                         intent(in)    :: ny
    double precision, dimension(ny), intent(out)   :: y

    integer, dimension(2) :: Wshape
    integer               :: iw1, iw2
    integer               :: iv1, iv2
    integer               :: nnodes1, nnodes2
    integer               :: ilayer

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `eval'."
       stop
    end if

    if (nx /= net%nnodes(1)) then
       write(0,*) "Error: wrong number of input nodes in `eval': ", nx
       stop
    end if

    if (ny /= net%nnodes(net%nlayers)) then
       write(0,*) "Error: wrong number of output nodes in `eval': ", ny
       stop
    end if

    net%value(1:nx) = x(1:nx)

    iw1 = 1
    iv1 = 1
    nnodes2 = 1 ! just to avoid the gfortran warning
    nnodes1 = net%nnodes(1)
    do ilayer = 1, net%nlayers-1
       iw2 = net%iw(ilayer+1)
       iv2 = net%iv(ilayer+1)

       net%value(iv2) = 1.0d0

       nnodes2 = net%nnodes(ilayer+1)
       Wshape(1:2) = (/ nnodes2, nnodes1+1 /)

       net%work(1:nnodes2) = matmul(         &
            reshape(net%W(iw1:iw2), Wshape), &
            net%value(iv1:iv2)               )

       iv1 = iv2 + 1
       iw1 = iw2 + 1

       call ff_activate(net%f_a(ilayer+1),             &
                        net%work(1:nnodes2),           &
                        net%value(iv1:iv1+nnodes2-1),  &
                        net%deriv(iv1:iv1+nnodes2-1)   )

       nnodes1 = nnodes2
    end do

    if (nnodes1 /= ny) stop 999 ! DEBUG only

    y(1:ny) = net%value(iv1:iv1+nnodes2-1)

    net%evaluated   = .true.
    net%derivatives = .false.

  end subroutine ff_eval

  !--------------------------------------------------------------------!
  !            derivative with respect to the input vector             !
  !--------------------------------------------------------------------!

  subroutine ff_deriv(net, nx, ny, dy)

    implicit none

    type(Network),                      intent(inout) :: net
    integer,                            intent(in)    :: nx, ny
    double precision, dimension(ny,nx), intent(out)   :: dy

    integer, dimension(2) :: Wshape
    integer               :: nnodes0, nnodes1, nnodes2
    integer               :: ilayer, i
    integer               :: iv1, iv2, iw1, iw2

    if (.not. net%init) then
       write(0,*) "Error: network not initialized in `deriv'."
       stop
    end if

    if (.not. net%evaluated) then
       write(0,*) "Error: network not evaluated in `deriv'."
       stop
    end if

    if (nx /= net%nnodes(1)) then
       write(0,*) "Error: wrong number of input nodes in `deriv': ", nx
       stop
    end if

    if (ny /= net%nnodes(net%nlayers)) then
       write(0,*) "Error: wrong number of output nodes in `deriv': ", ny
       stop
    end if

    ! reset all derivatives:
    ! (important to get the zeros for the bias weights)
    net%D(:) = 0.0d0

    ! start with the unity matrix:
    net%work4(:,:) = 0.0d0
    do i = 1, net%nnodes_max+1
       net%work4(i,i) = 1.0d0
    end do

    iw1 = 1
    iv1 = 1
    nnodes2 = 1 ! just to avoid the gfortran warning
    nnodes1 = net%nnodes(1)
    nnodes0 = nnodes1
    layers : do ilayer = 1, net%nlayers-1
       iw2 = net%iw(ilayer+1)
       iv2 = net%iv(ilayer+1)
       nnodes2 = net%nnodes(ilayer+1)

       ! construct diagonal matrix from stored derivatives:
       !
       !       / f'_i1    0      ...   0    \
       !      |   0      f'_i2    0   ...    |
       ! F' = |  ...     ...          ...    |
       !       \  0      ...          f'_iN /
       !
       net%work2(1:nnodes2,1:nnodes2) = 0.0d0
       do i = 1, nnodes2
          net%work2(i,i) = net%deriv(iv2+i)
       end do

       ! smaller matrices --> we don't need the bias weights here:
       Wshape(1:2) = (/ nnodes2, nnodes1 /)

       ! matrix of derivatives:
       !
       !             / d n_i1/d n_(i-1)1  ... d n_i1/d n_(i-1)M \
       ! D = F'*W = |      ...                      ...          |
       !             \ d n_iN/d n_(i-1)1  ... d n_iN/d n_(i-1)M /
       !
       net%work3(1:nnodes2,1:nnodes1) = matmul(      &
            net%work2(1:nnodes2,1:nnodes2),          &
            reshape(net%W(iw1:iw2-nnodes2), Wshape)  )

       ! store the result:
       net%D(iw1:iw2-nnodes2) = reshape(net%work3(1:nnodes2,1:nnodes1), &
                                        (/nnodes2*nnodes1/)             )

       ! combine with the matrix of the previous layer:
       net%work4(1:nnodes2,1:nnodes0) = matmul( &
            net%work3(1:nnodes2,1:nnodes1),     &
            net%work4(1:nnodes1,1:nnodes0)      )

       iv1 = iv2 + 1
       iw1 = iw2 + 1

       nnodes1 = nnodes2
    end do layers

    dy(1:ny,1:nx) = net%work4(1:nnodes2,1:nnodes0)

    net%derivatives = .true.

  end subroutine ff_deriv

  !--------------------------------------------------------------------!
  !                                                                    !
  !                        activation functions                        !
  !                                                                    !
  ! Function types:                                                    !
  !                                                                    !
  !   0 : linear function f(x) = x                                     !
  !   1 : hyperbolic tangent, y in [-1:1]                              !
  !   2 : sigmoid,            y in [ 0:1]                              !
  !   3 : modified tanh,      y in [-1.7159:1.7159]  f(+/-1) = +/-1    !
  !   4 : tanh & linear twisting term                                  !
  ! [3 & 4 Montavon, Orr, Muller Neural Networks: Tricks of the Trade] !
  !                                                                    !
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !                 return both: value and derivative                  !
  !--------------------------------------------------------------------!

  elemental subroutine ff_activate(t, x, y, dy)

    implicit none

    integer,          intent(in)  :: t
    double precision, intent(in)  :: x
    double precision, intent(out) :: y
    double precision, intent(out) :: dy

    double precision :: tanhbx
    double precision, parameter :: a = 1.7159d0
    double precision, parameter :: b = 0.666666666666667d0
    double precision, parameter :: c = 0.1d0

    select case(t)
    case(0)
       y  = x
       dy = 1.0
    case(1)
       y  = tanh(x)
       dy = 1.0d0 - y*y
    case(2)
       y  = 1.0d0/(1.0d0 + exp(-x))
       dy = y*(1.0d0 - y)
    case(3)
       tanhbx = tanh(b*x)
       y  = a*tanhbx
       dy = a*(1.0d0 - tanhbx*tanhbx)*b
    case(4)
       tanhbx = tanh(b*x)
       y  = a*tanhbx + c*x
       dy = a*(1.0d0 - tanhbx*tanhbx)*b + c
    case default
       y  = 0.0d0
       dy = 0.0d0
    end select

  end subroutine ff_activate

  !--------------------------------------------------------------------!
  !                                                                    !
  !                some utilities, mostly for debugging                !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine print_mat(A, unit)

    implicit none

    double precision, dimension(:,:), intent(in) :: A
    integer,                          intent(in) :: unit

    integer                     :: n1, n2
    integer                     :: i1, i2

    double precision            :: val
    double precision, parameter :: EPS = 1.0d-12

    n1 = size(A(:,1))
    n2 = size(A(1,:))

    do i1 = 1, n1
       do i2 = 1, n2
          val = A(i1,i2)
          if (abs(val) < EPS) val = 0.0d0
          write(unit, '(1x,ES15.8,2x)', advance='no') val
       end do
       write(unit,*)
    end do

  end subroutine print_mat

end module feedforward





!shimamura
!-----------------------------------------------------------------------
!     chebyshev.f90 - Chebyshev polynomials and their derivatives
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!-----------------------------------------------------------------------
! 2015-12-27 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module chebyshev

  implicit none
  private

  public :: chebyshev_polynomial, &
            chebyshev_polynomial_d1

contains

  !--------------------------------------------------------------------!
  !                   Evaluate Chebyshev polynomials                   !
  !--------------------------------------------------------------------!

  function chebyshev_polynomial(r, r0, r1, n) result(T)
    ! Arguments:
    !    r        function argument
    !    r0, r1   the Chebyshev polynomials will be rescaled from [-1,1]
    !             to the interval [r0,r1]
    !    n        maximum polynomial order
    ! Returns:
    !    T(i)  with i=1,n+1  where T(i) is the Chebyshev polynomial of
    !    order (i-1)
    !
    ! The Chebyshev polynomials obey the following recurrence relation:
    !    T[0](x) = 1
    !    T[1](x) = x
    !    T[n+1](x) = 2x T[n](x) - T[n-1](x)

    implicit none

    double precision, intent(in)     :: r, r0, r1
    integer,          intent(in)     :: n
    double precision, dimension(n+1) :: T

    integer          :: i
    double precision :: x

    x = (2.0d0*r - r0 - r1)/(r1 - r0)

    T(1) = 1.0d0
    if (n > 0) then
       T(2) = x
       do i = 3, n+1
          T(i) = 2.0d0*x*T(i-1) - T(i-2)
       end do
    end if

  end function chebyshev_polynomial

  !--------------------------------------------------------------------!
  !             First derivative of Chebyshev polynomials              !
  !--------------------------------------------------------------------!

  function chebyshev_polynomial_d1(r, r0, r1, n) result(dT)
    ! Arguments:
    !    r        function argument
    !    r0, r1   the Chebyshev polynomials will be rescaled from [-1,1]
    !             to the interval [r0,r1]
    !    n        maximum polynomial order
    ! Returns:
    !    dT(i)  with i=1,n+1  where dT(i) is the first derivative of the
    !    Chebyshev polynomial of order (i-1)
    !
    ! The derivatives of the Chebyshev polynomials obey the following
    ! recurrence relation:
    !    dT[n](x)/dx = n U[n-1](x) with n = 1,...
    ! where
    !    U[0](x) = 1
    !    U[1](x) = 2x
    !    U[n+1](x) = 2x U[n](x) - U[n-1](x)
    ! are the Chebyshev polynomials of the second kind.

    implicit none

    double precision, intent(in)     :: r, r0, r1
    integer,          intent(in)     :: n
    double precision, dimension(n+1) :: dT

    integer          :: i
    double precision :: x
    double precision :: U1, U2, U3

    x = (2.0d0*r - r0 - r1)/(r1 - r0)

    dT(1) = 0.0d0
    if (n > 0) then
       U1 = 1.0d0
       dT(2) = U1
       U2 = 2.0d0*x
       do i = 3, n+1
          dT(i) = U2*dble(i-1)
          U3 = 2.0d0*x*U2 - U1
          U1 = U2
          U2 = U3
       end do
    end if

    ! inner derivative (from rescaling)
    dT = dT*2.0d0/(r1 - r0)

  end function chebyshev_polynomial_d1

end module chebyshev

!----------------------------------------------------------------------!
!                              Unit test                               !
!----------------------------------------------------------------------!

!!$ program chebyshev_test
!!$
!!$   use chebyshev, only: chebyshev_polynomial, chebyshev_polynomial_d1
!!$
!!$   implicit none
!!$
!!$   integer, parameter          :: N = 100
!!$   integer, parameter          :: O = 10
!!$   double precision, parameter :: R0 = -2.0d0
!!$   double precision, parameter :: R1 =  5.0d0
!!$
!!$   double precision, dimension(O+1) :: T, dT, T_prev
!!$
!!$   double precision    :: r, dr
!!$   integer             :: i
!!$   character(len=1024) :: frmt
!!$
!!$   write(frmt, *) O + 2
!!$   frmt = "(1x," // trim(adjustl(frmt)) // "(ES18.8,1x))"
!!$
!!$   open(20, file="chebytest-values.dat", status="replace", action="write")
!!$   open(21, file="chebytest-derivs.dat", status="replace", action="write")
!!$   open(22, file="chebytest-nderivs.dat", status="replace", action="write")
!!$
!!$   dr = (R1 - R0)/dble(N - 1)
!!$   r = r0
!!$   T_prev = 0.0d0
!!$   do i = 1, N
!!$      T = chebyshev_polynomial(r, R0, R1, O)
!!$      dT = chebyshev_polynomial_d1(r, R0, R1, O)
!!$      write(20, frmt) r, T(1:O+1)
!!$      write(21, frmt) r, dT(1:O+1)
!!$      if (i > 1) then
!!$         write(22, frmt) r - 0.5d0*dr, (T - T_prev)/dr
!!$      end if
!!$      T_prev = T
!!$      r = r + dr
!!$   end do
!!$
!!$   close(20)
!!$   close(21)
!!$   close(22)
!!$
!!$ end program chebyshev_test
!shimamura




!shimamura
!-----------------------------------------------------------------------
! sfbasis.f90 - Basis for structural fingerprints of atomic environments
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!-----------------------------------------------------------------------
! 2015-12-27 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module sfbasis

  !use io, only: io_unit
  use io,          only: nfile, myid, nodes

  use chebyshev, only: chebyshev_polynomial, &
                       chebyshev_polynomial_d1

  implicit none
  private
  save

  public :: new_SFBasis,            &
            del_SFBasis,            &
            save_SFBasis,           &
            load_SFBasis,           &
            save_SFBasis_ASCII,     &
            load_SFBasis_ASCII,     &
            sfb_print_info,         &
            sfb_set_typeid,         &
            sfb_set_typespin,       &
            sfb_eval,               &
            sfb_reconstruct_radial, &
            sfb_reconstruct_angular

  type, public :: FingerprintBasis
     logical                                       :: initialized = .false.
     integer                                       :: r_order
     double precision                              :: r_Rc
     integer                                       :: r_N
     integer                                       :: a_order
     double precision                              :: a_Rc
     integer                                       :: a_N
     integer                                       :: r_i1, r_f1
     integer                                       :: r_i2, r_f2
     integer                                       :: a_i1, a_f1
     integer                                       :: a_i2, a_f2
     integer                                       :: N
     integer                                       :: num_types
     logical                                       :: multi
     character(len=2), dimension(:),   allocatable :: atom_types
     integer,          dimension(:),   allocatable :: typeid
     double precision, dimension(:),   allocatable :: typespin
     integer                                       :: num_values
  end type FingerprintBasis

  double precision, parameter, private :: PI     = 3.14159265358979d0
  double precision, parameter, private :: PI_INV = 1.0d0/PI
  double precision, parameter, private :: PI2    = 2.0d0*PI
  double precision, parameter, private :: EPS    = 1.0d-12

contains

  !--------------------------------------------------------------------!
  !            create a new Structural Fingerprint basis               !
  !--------------------------------------------------------------------!

  function new_SFBasis(num_types, atom_types, radial_order, &
                       angular_order, radial_Rc, angular_Rc) result(sfb)
    ! Arguments:
    !   num_types       number of atomic species
    !   atom_types(i)   i-th atomic species (2 characters)
    !   radial_order    expansion order for the radial basis
    !   angular_order   expansion order for the angular basis
    !   radial_Rc       cutoff radius for radial basis
    !   angular_Rc      cutoff radius for angular basis
    !
    ! Returns:
    !   sfb             allocated instance of FingerprintBasis

    implicit none

    integer,                                intent(in) :: num_types
    character(len=*), dimension(num_types), intent(in) :: atom_types
    integer,                                intent(in) :: radial_order
    integer,                                intent(in) :: angular_order
    double precision,                       intent(in) :: radial_Rc
    double precision,                       intent(in) :: angular_Rc
    type(FingerprintBasis)                             :: sfb

    integer :: i, s

    sfb%num_types = num_types
    sfb%r_order = radial_order
    sfb%a_order = angular_order
    sfb%r_Rc = radial_Rc
    sfb%a_Rc = angular_Rc

    sfb%r_N = sfb%r_order + 1
    sfb%a_N = sfb%a_order + 1
    sfb%num_values = max(sfb%r_N, sfb%a_N)
    sfb%N = sfb%r_N + sfb%a_N
    sfb%r_i1 = 1
    sfb%r_f1 = sfb%r_i1 + sfb%r_N - 1
    sfb%a_i1 = sfb%r_f1 + 1
    sfb%a_f1 = sfb%a_i1 + sfb%a_N - 1
    sfb%r_i2 = sfb%a_f1 + 1
    sfb%r_f2 = sfb%r_i2 + sfb%r_N - 1
    sfb%a_i2 = sfb%r_f2 + 1
    sfb%a_f2 = sfb%a_i2 + sfb%a_N - 1
    if (sfb%num_types > 1) then
       sfb%multi = .true.
       sfb%N = 2*sfb%N
    else
       sfb%multi = .false.
    end if

    allocate(sfb%atom_types(num_types),      &
             sfb%typeid(num_types),          &
             sfb%typespin(num_types))
    sfb%atom_types = atom_types

    do i = 1, num_types
       sfb%typeid(i) = i
    end do

    s = -num_types/2
    do i = 1, num_types
       if ((s == 0) .and. (mod(num_types, 2) == 0)) s = s + 1
       sfb%typespin(i) = dble(s)
       s = s + 1
    end do

    
    !!shimamura for s-ANN
    !write(*,*)"num_types=",num_types
    !write(*,*)"sfb%typespin=",sfb%typespin


    sfb%initialized = .true.

  end function new_SFBasis

  !--------------------------------------------------------------------!
  !         delete (deallocate) a Structural Fingerprint basis         !
  !--------------------------------------------------------------------!

  subroutine del_SFBasis(sfb)

    implicit none

    type(FingerprintBasis), intent(inout) :: sfb

    call sfb_assert_init(sfb)

    deallocate(sfb%atom_types, sfb%typeid, sfb%typespin)
    sfb%initialized = .false.

  end subroutine del_SFBasis

  !--------------------------------------------------------------------!
  !                    print information to stdout                     !
  !--------------------------------------------------------------------!

  subroutine sfb_print_info(sfb)

    implicit none

    type(FingerprintBasis), intent(in) :: sfb
    character(len=1024) :: frmt

    write(*,'(" Radial cutoff : ",F7.3)') sfb%r_Rc
    write(*,'(" Angular cutoff: ",F7.3)') sfb%a_Rc
    write(*,'(" Radial order  : ",I3)') sfb%r_order
    write(*,'(" Angular order : ",I3)') sfb%a_order
    write(*,'(" Atom types    : ")', advance='no')
    write(frmt, *) sfb%num_types
    frmt = '(' // trim(adjustl(frmt))  // '(A2,1x))'
    write(*,frmt) sfb%atom_types
    write(*,'(" Total number of basis functions: ", I3)') sfb%N

  end subroutine sfb_print_info

  !============================ properties ============================!


  subroutine sfb_set_typeid(sfb, typeid)

    implicit none

    type(FingerprintBasis), intent(inout) :: sfb
    integer, dimension(:),  intent(in)    :: typeid

    call sfb_assert_init(sfb)

    if (size(typeid) < sfb%num_types) then
       write(0,*) "Error: incompatible type ID list in `sfb_set_typeid'."
       stop
    end if

    sfb%typeid(1:sfb%num_types) = typeid(1:sfb%num_types)

  end subroutine sfb_set_typeid

  !--------------------------------------------------------------------!

  subroutine sfb_set_typespin(sfb, typespin)

    implicit none

    type(FingerprintBasis),         intent(inout) :: sfb
    double precision, dimension(:), intent(in)    :: typespin

    call sfb_assert_init(sfb)

    if (size(typespin) < sfb%num_types) then
       write(0,*) "Error: incompatible type ID list in `sfb_set_typespin'."
       stop
    end if

    sfb%typespin(1:sfb%num_types) = typespin(1:sfb%num_types)

  end subroutine sfb_set_typespin


  !=============================== I/O ================================!


  !--------------------------------------------------------------------!
  !    Read/write Structural Fingerprint Basis from/to file or unit    !
  !                                                                    !
  !    The *_ASCII procedures read/write plain text files; the other   !
  !    procedures just use binary I/O.                                 !
  !--------------------------------------------------------------------!

  subroutine save_SFBasis(sfb, file, unit)

    implicit none

    type(FingerprintBasis),     intent(in) :: sfb
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: u

    call sfb_assert_init(sfb)

    if (present(unit)) then
       u = unit

    !shimamura
    !else if (present(file)) then
    !   u = io_unit()
    !   open(u, file=trim(file), status='replace', action='write', &
    !        form='unformatted')

    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `save_SFBasis'."
       return
    end if

    write(u) sfb%r_order
    write(u) sfb%a_order
    write(u) sfb%r_Rc
    write(u) sfb%a_Rc
    write(u) sfb%r_N, sfb%a_N, sfb%N
    write(u) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    write(u) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    write(u) sfb%num_values
    write(u) sfb%num_types
    write(u) sfb%atom_types(:)
    write(u) sfb%typeid(:)
    write(u) sfb%typespin(:)

    if (.not. present(unit)) close(u)

  end subroutine save_SFBasis

  !--------------------------------------------------------------------!

  function load_SFBasis(file, unit) result(sfb)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(FingerprintBasis)                 :: sfb

    integer :: u

    if (present(unit)) then
       u = unit

    !shimamura
    !else if (present(file)) then
    !   u = io_unit()
    !   open(u, file=trim(file), action='read', form='unformatted')

    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `load_SFBasis'."
       return
    end if

    read(u) sfb%r_order
    read(u) sfb%a_order
    read(u) sfb%r_Rc
    read(u) sfb%a_Rc
    read(u) sfb%r_N, sfb%a_N, sfb%N
    read(u) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    read(u) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    read(u) sfb%num_values
    read(u) sfb%num_types
    allocate(sfb%atom_types(sfb%num_types),  &
             sfb%typeid(sfb%num_types),      &
             sfb%typespin(sfb%num_types))
    read(u) sfb%atom_types(:)
    read(u) sfb%typeid(:)
    read(u) sfb%typespin(:)

    if (sfb%num_types > 1) then
       sfb%multi = .true.
    else
       sfb%multi = .false.
    end if
    sfb%initialized = .true.

    if (.not. present(unit)) close(u)

  end function load_SFBasis

  !--------------------------------------------------------------------!

  subroutine save_SFBasis_ASCII(sfb, file, unit)

    implicit none

    type(FingerprintBasis),     intent(in) :: sfb
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    character(len=*), parameter :: DFRMT = '(4(1x,ES24.17))'
    character(len=*), parameter :: IFRMT = '(4(1x,I17))'
    character(len=*), parameter :: AFRMT = '(4(1x,A))'

    integer :: u, i

    call sfb_assert_init(sfb)

    if (present(unit)) then
       u = unit

    !shimamura
    !else if (present(file)) then
    !   u = io_unit()
    !   open(u, file=trim(file), status='replace', action='write')

    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `save_SFBasis_ASCII'."
       return
    end if

    write(u,*) sfb%r_order
    write(u,*) sfb%a_order
    write(u,*) sfb%r_Rc
    write(u,*) sfb%a_Rc
    write(u,*) sfb%r_N, sfb%a_N, sfb%N
    write(u,*) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    write(u,*) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    write(u,*) sfb%num_values
    write(u,*) sfb%num_types
    write(u,AFRMT) (sfb%atom_types(i), i=1,sfb%num_types)
    write(u,IFRMT) (sfb%typeid(i), i=1,sfb%num_types)
    write(u,DFRMT) (sfb%typespin(i), i=1,sfb%num_types)

    if (.not. present(unit)) close(u)

  end subroutine save_SFBasis_ASCII

  !--------------------------------------------------------------------!

  function load_SFBasis_ASCII(file, unit) result(sfb)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(FingerprintBasis)                 :: sfb

    character(len=*), parameter :: DFRMT = '(4(1x,ES24.17))'
    character(len=*), parameter :: IFRMT = '(4(1x,I17))'
    character(len=*), parameter :: AFRMT = '(4(1x,A))'

    integer :: u

    if (present(unit)) then
       u = unit
   
    !shimamura
    !else if (present(file)) then
    !   u = io_unit()
    !   open(u, file=trim(file), action='read')

    else
       write(0,*) "Error: neither unit number nor file name given " // &
                  "in `load_SFBasis_ASCII'."
       return
    end if

    read(u,*) sfb%r_order
    read(u,*) sfb%a_order
    read(u,*) sfb%r_Rc
    read(u,*) sfb%a_Rc
    read(u,*) sfb%r_N, sfb%a_N, sfb%N
    read(u,*) sfb%r_i1, sfb%r_f1, sfb%r_i2, sfb%r_f2
    read(u,*) sfb%a_i1, sfb%a_f1, sfb%a_i2, sfb%a_f2
    read(u,*) sfb%num_values
    read(u,*) sfb%num_types
    allocate(sfb%atom_types(sfb%num_types),  &
             sfb%typeid(sfb%num_types),      &
             sfb%typespin(sfb%num_types))
    read(u,AFRMT) sfb%atom_types(:)
    read(u,IFRMT) sfb%typeid(:)
    read(u,DFRMT) sfb%typespin(:)

    if (sfb%num_types > 1) then
       sfb%multi = .true.
    else
       sfb%multi = .false.
    end if
    sfb%initialized = .true.

    if (.not. present(unit)) close(u)

  end function load_SFBasis_ASCII


  !========================= basis evaluation =========================!

  !shimamura
  subroutine sfb_eval(sfb, itype0, coo0, nat, itype1, coo1, nv, &
                      values, deriv0, deriv1, derivh)
  !subroutine sfb_eval(sfb, itype0, coo0, nat, itype1, coo1, nv, &
  !                    values, deriv0, deriv1)

    implicit none

    type(FingerprintBasis),                          intent(inout) :: sfb
    integer,                                         intent(in)    :: itype0
    double precision, dimension(3),                  intent(in)    :: coo0
    integer,                                         intent(in)    :: nat
    integer,          dimension(nat),                intent(in)    :: itype1
    double precision, dimension(3,nat),              intent(in)    :: coo1
    integer,                                         intent(in)    :: nv
    double precision, dimension(nv),                 intent(out)   :: values
    double precision, dimension(3,nv),     optional, intent(out)   :: deriv0
    double precision, dimension(3,nv,nat), optional, intent(out)   :: deriv1

    !shimamura
    double precision, dimension(6,nv),      optional, intent(out)   :: derivh
    double precision, dimension(:,:), allocatable :: sfb_dstress

    double precision, dimension(sfb%num_values)   :: sfb_values
    double precision, dimension(:,:), allocatable :: sfb_deriv_i
    double precision, dimension(:,:), allocatable :: sfb_deriv_j
    double precision, dimension(:,:), allocatable :: sfb_deriv_k

    logical                        :: do_deriv
    double precision, dimension(3) :: R_ij, R_ik
    double precision               :: d_ij, d_ik
    double precision               :: cos_ijk
    double precision               :: s_j, s_k
    integer                        :: j, k, i1, i2, N

    !shimamura
    logical                        :: do_strs

    call sfb_assert_init(sfb)

    if (nv /= sfb%N) then
       write(0,*) "Error: wrong number of basis functions in `sfb_eval'."
       stop
    end if

    !shimamura
    do_strs = .false.

    if (present(deriv0) .and. present(deriv1)) then
       do_deriv = .true.
       deriv0(:,:) = 0.0d0
       deriv1(:,:,:) = 0.0d0
       allocate(sfb_deriv_i(3, sfb%num_values), &
                sfb_deriv_j(3, sfb%num_values), &
                sfb_deriv_k(3, sfb%num_values))

       !shimamura
       if (present(derivh)) then
           do_strs = .true.
           derivh(:,:) = 0.0d0
           allocate(sfb_dstress(6, sfb%num_values))
       end if

    else
       do_deriv = .false.
    end if


!    !shimamura
!    write(*,*)"nv            =",nv
!    write(*,*)"sfb%num_values=",sfb%num_values
!    write(*,*)"sfb%r_N       =",sfb%r_N
!    write(*,*)"sfb%r_i1      =",sfb%r_i1
!    write(*,*)"sfb%r_f1      =",sfb%r_f1
!    write(*,*)"sfb%r_i2      =",sfb%r_i2
!    write(*,*)"sfb%r_f2      =",sfb%r_f2
!    write(*,*)"sfb%a_N       =",sfb%a_N
!    write(*,*)"sfb%a_i1      =",sfb%a_i1
!    write(*,*)"sfb%a_f1      =",sfb%a_f1
!    write(*,*)"sfb%a_i2      =",sfb%a_i2
!    write(*,*)"sfb%a_f2      =",sfb%a_f2


    values(1:sfb%N) = 0.0d0
    s_j = 1.0d0

    for_j : do j = 1, nat
       R_ij = coo1(1:3, j) - coo0(1:3)
       d_ij = sqrt(dot_product(R_ij, R_ij))
       if ((d_ij <= sfb%r_Rc) .and. (d_ij > EPS)) then

          ! evaluate radial basis functions
          i1 = sfb%r_i1
          i2 = sfb%r_f1
          N = sfb%r_N
          if (do_deriv) then

             !shimamura
             if(do_strs)then
             call sfb_radial_stress(sfb, R_ij, d_ij, sfb_values, &
                             deriv_i=sfb_deriv_i, deriv_j=sfb_deriv_j, dstress=sfb_dstress)
             derivh(1:6, i1:i2) = derivh(1:6, i1:i2) + sfb_dstress(1:6, 1:N)
             else
             call sfb_radial(sfb, R_ij, d_ij, sfb_values, &
                             deriv_i=sfb_deriv_i, deriv_j=sfb_deriv_j)
             end if

             values(i1:i2) = values(i1:i2) + sfb_values(1:N)
             deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) + sfb_deriv_i(1:3, 1:N)
             deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) + sfb_deriv_j(1:3, 1:N)

          else
             call sfb_radial(sfb, R_ij, d_ij, sfb_values)
             values(i1:i2) = values(i1:i2) + sfb_values(1:N)
          end if

          ! redundant radial basis in case of multi-component systems
          i1 = sfb%r_i2
          i2 = sfb%r_f2
          N = sfb%r_N
          if (sfb%multi) then
             s_j = sfb%typespin(sfb%typeid(itype1(j)))
             values(i1:i2) = values(i1:i2) + s_j*sfb_values(1:N)
             if (do_deriv) then
                deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) &
                                   + s_j*sfb_deriv_i(1:3, 1:N)
                deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) &
                                      + s_j*sfb_deriv_j(1:3, 1:N)
             end if

             !shimamura
             if(do_strs)then
                derivh(1:6, i1:i2) = derivh(1:6, i1:i2) + s_j*sfb_dstress(1:6, 1:N)
             end if

          end if

       end if  ! within radial cutoff

       if (d_ij > sfb%a_Rc) cycle for_j
       for_k : do k = j+1, nat
          R_ik = coo1(1:3, k) - coo0(1:3)
          d_ik = sqrt(dot_product(R_ik, R_ik))
          if ((d_ik > sfb%a_Rc) .or. (d_ik < EPS)) cycle for_k
          cos_ijk = dot_product(R_ij, R_ik)/(d_ij*d_ik)

          ! evaluate angular basis functions
          i1 = sfb%a_i1
          i2 = sfb%a_f1
          N = sfb%a_N

          !shimamura
          if (do_deriv) then

             if(do_strs)then
                call sfb_angular_stress(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, &
                                 sfb_values, deriv_i=sfb_deriv_i,      &
                                 deriv_j=sfb_deriv_j, deriv_k=sfb_deriv_k, dstress=sfb_dstress)
                derivh(1:6, i1:i2) = derivh(1:6, i1:i2) + sfb_dstress(1:6, 1:N)
             else
                call sfb_angular(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, &
                                 sfb_values, deriv_i=sfb_deriv_i,      &
                                 deriv_j=sfb_deriv_j, deriv_k=sfb_deriv_k)
             end if

             values(i1:i2) = values(i1:i2) + sfb_values(1:N)
             deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) + sfb_deriv_i(1:3, 1:N)
             deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) + sfb_deriv_j(1:3, 1:N)
             deriv1(1:3, i1:i2, k) = deriv1(1:3, i1:i2, k) + sfb_deriv_k(1:3, 1:N)
          else
             call sfb_angular(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, sfb_values)
             values(i1:i2) = values(i1:i2) + sfb_values(1:N)
          end if

          ! redundant angular basis in case of multi-component systems
          i1 = sfb%a_i2
          i2 = sfb%a_f2
          N = sfb%a_N
          if (sfb%multi) then
             s_k = sfb%typespin(sfb%typeid(itype1(k)))
             values(i1:i2) = values(i1:i2) + s_j*s_k*sfb_values(1:N)
             if (do_deriv) then
                deriv0(1:3, i1:i2) = deriv0(1:3, i1:i2) &
                                   + s_j*s_k*sfb_deriv_i(1:3, 1:N)
                deriv1(1:3, i1:i2, j) = deriv1(1:3, i1:i2, j) &
                                      + s_j*s_k*sfb_deriv_j(1:3, 1:N)
                deriv1(1:3, i1:i2, k) = deriv1(1:3, i1:i2, k) &
                                      + s_j*s_k*sfb_deriv_k(1:3, 1:N)
             end if

             !shimamura
             if(do_strs)then
                derivh(1:6, i1:i2) = derivh(1:6, i1:i2) + s_j*s_k*sfb_dstress(1:6, 1:N)
             end if

          end if
       end do for_k
    end do for_j

    if (do_deriv) deallocate(sfb_deriv_i, sfb_deriv_j, sfb_deriv_k)

    !shimamura
    if (do_strs) deallocate(sfb_dstress)


!    !shimamura
!    do i1=1,nv
!       write(*,'(I4,6(1x,ES14.6))')i1, derivh(1:6,i1)
!    end do
!    write(*,'(I4,6(1x,ES14.6))')9999,sum(derivh(1,1:22)),sum(derivh(2,1:22))&
!                                   &,sum(derivh(3,1:22)),sum(derivh(4,1:22))&
!                                   &,sum(derivh(5,1:22)),sum(derivh(6,1:22)) 
    !write(*,'(I4,6(1x,ES14.6))')9999,sum(derivh(1,1:nv)),sum(derivh(2,1:nv))&
    !                               &,sum(derivh(3,1:nv)),sum(derivh(4,1:nv))&
    !                               &,sum(derivh(5,1:nv)),sum(derivh(6,1:nv)) 

  end subroutine sfb_eval


  !======================= Basis Set Expansion ========================!


  !--------------------------------------------------------------------!
  ! reconstruct radial distribution function from basis set expansion  !
  !--------------------------------------------------------------------!

  subroutine sfb_reconstruct_radial(sfb, coeff, nx, x, y)

    implicit none

    !------------------------------------------------------------------!
    ! sfb         Instance of FingerprintBasis                         !
    ! coeff(i)    Coefficient of the i-th basis function               !
    ! nx          Grid points for function evaluation                  !
    ! x(i)        x value of the i-th grid point (output)              !
    ! y(i)        y (function) value of the i-th grid point (output)   !
    !------------------------------------------------------------------!

    type(FingerprintBasis),               intent(in)  :: sfb
    double precision, dimension(sfb%r_N), intent(in)  :: coeff
    integer,                              intent(in)  :: nx
    double precision, dimension(nx),      intent(out) :: x
    double precision, dimension(nx),      intent(out) :: y

    double precision, dimension(sfb%r_N) :: f

    double precision :: dx, r_over_Rc, w
    integer :: ix, ic

    call sfb_assert_init(sfb)

    dx = sfb%r_Rc/dble(nx - 1)
    do ix = 1, nx
       x(ix) = dble(ix - 1)*dx
    end do

    ! The Chebyshev polynomials are orthogonal with respect to the
    ! weight w(x) = 1/sqrt(1 - x**2), i.e.,
    !
    ! \int_{-1}^{1} T_n(x) T_m(x) w(x) dx = \delta_{mn} * k
    ! k = \pi/2 for m = n = 0, for all other m, n k = \pi
    !
    ! The values returned by the structural fingerprint basis are for
    ! the Chebyshev polynomials without any weights, so that we have to
    ! weight the reconstruction of the RDF here.  The expressen below
    ! looks a little bit different, because the interval has been
    ! rescaled from [-1:1] to [0:Rc].  Note that the weight for r = Rc
    ! is undefined but nevertheless evaluated here.  This means that the
    ! RDF becomes unreliable for r --> Rc
    do ix = 1, nx-1
       f = chebyshev_polynomial(x(ix), 0.0d0, sfb%r_Rc, sfb%r_order)
       r_over_Rc = x(ix)/sfb%r_Rc
       w = 0.5d0/sqrt(r_over_Rc - r_over_Rc*r_over_Rc)
       f = f*w*PI_INV
       f(1) = 0.5*f(1)
       y(ix) = 0.0d0
       do ic = sfb%r_N, 1, -1
          y(ix) = y(ix) + coeff(ic)*f(ic)
       end do
    end do
    y(nx) = 0.0d0

  end subroutine sfb_reconstruct_radial

  !--------------------------------------------------------------------!
  ! reconstruct angular distribution function from basis set expansion !
  !--------------------------------------------------------------------!

  subroutine sfb_reconstruct_angular(sfb, coeff, nx, x, y)

    implicit none

    !------------------------------------------------------------------!
    ! sfb         Instance of FingerprintBasis                         !
    ! coeff(i)    Coefficient of the i-th basis function               !
    ! nx          Grid points for function evaluation                  !
    ! x(i)        x value of the i-th grid point (output)              !
    ! y(i)        y (function) value of the i-th grid point (output)   !
    !------------------------------------------------------------------!

    type(FingerprintBasis),               intent(in)  :: sfb
    double precision, dimension(sfb%r_N), intent(in)  :: coeff
    integer,                              intent(in)  :: nx
    double precision, dimension(nx),      intent(out) :: x
    double precision, dimension(nx),      intent(out) :: y

    double precision, dimension(sfb%r_N) :: f

    double precision :: dx, r_over_PI, w
    integer :: ix, ic

    call sfb_assert_init(sfb)

    dx = PI/dble(nx - 1)
    do ix = 1, nx
       x(ix) = dble(ix - 1)*dx
    end do

    ! The Chebyshev polynomials are orthogonal with respect to a weight.
    ! See the 'radial' subroutine above for further comments.
    do ix = 1, nx-1
       f = chebyshev_polynomial(x(ix), 0.0d0, PI, sfb%a_order)
       r_over_PI = x(ix)/PI
       w = 0.5d0/sqrt(r_over_PI - r_over_PI*r_over_PI)
       f = f*w*PI_INV
       f(1) = 0.5*f(1)
       y(ix) = 0.0d0
       do ic = sfb%r_N, 1, -1
          y(ix) = y(ix) + coeff(ic)*f(ic)
       end do
    end do
    y(nx) = 0.0d0

  end subroutine sfb_reconstruct_angular


  !======================== private/auxiliary =========================!


  !--------------------------------------------------------------------!
  !        assert that a FingerprintBasis has been initialized         !
  !--------------------------------------------------------------------!

  subroutine sfb_assert_init(sfb)

    implicit none

    type(FingerprintBasis), intent(in) :: sfb

    if (.not. sfb%initialized) then
       write(0, *) "Error: FingerprintBasis not initialized."
       stop
    end if

  end subroutine sfb_assert_init


  !====================================================================!
  !                                                                    !
  !                          cutoff function                           !
  !                                                                    !
  !====================================================================!


  function sfb_fc(Rij, Rc) result(fc)

    implicit none

    double precision, intent(in) :: Rij, Rc
    double precision             :: fc

    if (Rij >= Rc) then
       fc  = 0.0d0
    else
       fc  =  0.5d0*(cos(PI/Rc*Rij) + 1.0d0)
    end if

  end function sfb_fc

  !--------------------------------------------------------------------!

  function sfb_fc_d1(Rij, Rc) result(dfc)

    implicit none

    double precision, intent(in) :: Rij, Rc
    double precision             :: dfc

    double precision :: a

    if (Rij >= Rc) then
       dfc = 0.0d0
    else
       a = PI/Rc
       dfc = -0.5d0*a*sin(a*Rij)
    end if

  end function sfb_fc_d1


  !====================================================================!
  !                                                                    !
  !                      generic basis functions                       !
  !                                                                    !
  !====================================================================!

  !shimamura
  subroutine sfb_radial_stress(sfb, R_ij, d_ij, values, deriv_i, deriv_j, dstress)

    implicit none

    type(FingerprintBasis),                     intent(inout) :: sfb
    double precision, dimension(3),             intent(in)    :: R_ij
    double precision,                           intent(in)    :: d_ij
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j

    double precision                     :: w_ij, dw_ij
    double precision, dimension(sfb%r_N) :: f, df
    integer                              :: i

    !shimamura
    double precision, dimension(:,:), optional, intent(out)   :: dstress
    !double precision                     :: DGi_Rij = 0.d0

    call sfb_assert_init(sfb)

    w_ij = sfb_fc(d_ij, sfb%r_Rc)
    f = chebyshev_polynomial(d_ij, 0.0d0, sfb%r_Rc, sfb%r_order)

    values(1:sfb%r_N) = w_ij*f(1:sfb%r_N)

    if (present(deriv_i) .and. present(deriv_j)) then
       dw_ij = sfb_fc_d1(d_ij, sfb%r_Rc)
       df = chebyshev_polynomial_d1(d_ij, 0.0d0, sfb%r_Rc, sfb%r_order)
       forall (i=1:sfb%r_N)
          deriv_i(:,i) = -R_ij/d_ij*(dw_ij*f(i) + w_ij*df(i))
       end forall

       if (present(dstress)) then
          dstress = 0.d0
          forall (i=1:sfb%r_N)
                 !DGi_Rij = (dw_ij*f(i) + w_ij*df(i))/d_ij
                 !dstress(1,i) = DGi_Rij*R_ij(1)*R_ij(1)
                 !dstress(2,i) = DGi_Rij*R_ij(2)*R_ij(2)
                 !dstress(3,i) = DGi_Rij*R_ij(3)*R_ij(3)
                 !dstress(4,i) = DGi_Rij*R_ij(2)*R_ij(3)
                 !dstress(5,i) = DGi_Rij*R_ij(3)*R_ij(1)
                 !dstress(6,i) = DGi_Rij*R_ij(1)*R_ij(2)

                 !dstress(1,i) = -(dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(1)*R_ij(1)
                 !dstress(2,i) = -(dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(2)*R_ij(2)
                 !dstress(3,i) = -(dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(3)*R_ij(3)
                 !dstress(4,i) = -(dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(2)*R_ij(3)
                 !dstress(5,i) = -(dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(3)*R_ij(1)
                 !dstress(6,i) = -(dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(1)*R_ij(2)
                 dstress(1,i) = (dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(1)*R_ij(1)
                 dstress(2,i) = (dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(2)*R_ij(2)
                 dstress(3,i) = (dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(3)*R_ij(3)
                 dstress(4,i) = (dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(2)*R_ij(3)
                 dstress(5,i) = (dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(3)*R_ij(1)
                 dstress(6,i) = (dw_ij*f(i) + w_ij*df(i))/d_ij*R_ij(1)*R_ij(2)
          end forall
       end if

       deriv_j(1:3,1:sfb%r_N) = -deriv_i(1:3,1:sfb%r_N)
    end if

  end subroutine sfb_radial_stress

  subroutine sfb_radial(sfb, R_ij, d_ij, values, deriv_i, deriv_j)

    implicit none

    type(FingerprintBasis),                     intent(inout) :: sfb
    double precision, dimension(3),             intent(in)    :: R_ij
    double precision,                           intent(in)    :: d_ij
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j

    double precision                     :: w_ij, dw_ij
    double precision, dimension(sfb%r_N) :: f, df
    integer                              :: i

    call sfb_assert_init(sfb)

    w_ij = sfb_fc(d_ij, sfb%r_Rc)
    f = chebyshev_polynomial(d_ij, 0.0d0, sfb%r_Rc, sfb%r_order)

    values(1:sfb%r_N) = w_ij*f(1:sfb%r_N)

    if (present(deriv_i) .and. present(deriv_j)) then
       dw_ij = sfb_fc_d1(d_ij, sfb%r_Rc)
       df = chebyshev_polynomial_d1(d_ij, 0.0d0, sfb%r_Rc, sfb%r_order)
       forall (i=1:sfb%r_N)
          deriv_i(:,i) = -R_ij/d_ij*(dw_ij*f(i) + w_ij*df(i))
       end forall
       deriv_j(1:3,1:sfb%r_N) = -deriv_i(1:3,1:sfb%r_N)
    end if

  end subroutine sfb_radial

  !--------------------------------------------------------------------!
  !shimamura
  subroutine sfb_angular_stress(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, values, &
                         deriv_i, deriv_j, deriv_k, dstress)

    implicit none

    type(FingerprintBasis),                     intent(inout) :: sfb
    double precision, dimension(3),             intent(in)    :: R_ij, R_ik
    double precision,                           intent(in)    :: d_ij, d_ik
    double precision,                           intent(in)    :: cos_ijk
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j
    double precision, dimension(:,:), optional, intent(out)   :: deriv_k

    !shimamura
    double precision, dimension(:,:), optional, intent(out)   :: dstress
    double precision, dimension(3)       :: djpart1=0.d0,djpart2=0.d0
    double precision, dimension(3)       :: dkpart1=0.d0,dkpart2=0.d0

    double precision                     :: w_ijk
    double precision                     :: fc_j, dfc_j, fc_k, dfc_k
    double precision, dimension(sfb%a_N) :: f, df
    double precision                     :: id_ij2, id_ik2, id_ij_ik
    double precision, dimension(3)       :: di_cos_ikj, dj_cos_ikj, dk_cos_ikj
    double precision, dimension(3)       :: di_w_ijk, dj_w_ijk, dk_w_ijk
    integer                              :: i

    call sfb_assert_init(sfb)

    fc_j = sfb_fc(d_ij, sfb%a_Rc)
    fc_k = sfb_fc(d_ik, sfb%a_Rc)
    w_ijk = fc_j*fc_k
    f = chebyshev_polynomial(cos_ijk, 0.0d0, PI, sfb%a_order)

    values(1:sfb%a_N) = w_ijk*f

    if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
       dfc_j = sfb_fc_d1(d_ij, sfb%a_Rc)
       dfc_k = sfb_fc_d1(d_ik, sfb%a_Rc)
       df = chebyshev_polynomial_d1(cos_ijk, 0.0d0, PI, sfb%a_order)
       id_ij2 = 1.0d0/(d_ij*d_ij)
       id_ik2 = 1.0d0/(d_ik*d_ik)
       id_ij_ik = 1.0d0/(d_ij*d_ik)
       ! d/dR_i (cos_ijk)
       di_cos_ikj = cos_ijk*(R_ij*id_ij2 + R_ik*id_ik2) - (R_ij+R_ik)*id_ij_ik
       ! d/dR_j (cos_ijk)
       dj_cos_ikj = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
       ! d/dR_k (cos_ijk)
       dk_cos_ikj = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik
       ! d/dR_i (w_ijk)
       di_w_ijk = -(dfc_j*fc_k*R_ij/d_ij + fc_j*dfc_k*R_ik/d_ik)
       ! d/dR_j (w_ijk)
       dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
       ! d/dR_k (w_ijk)
       dk_w_ijk = fc_j*dfc_k*R_ik/d_ik

       forall (i=1:sfb%a_N)
          ! d/dR_i (w_ijk*f)
          deriv_i(:,i) = di_w_ijk(:)*f(i) + w_ijk*df(i)*di_cos_ikj(:)
          ! d/dR_j (w_ijk*f)
          deriv_j(:,i) = dj_w_ijk(:)*f(i) + w_ijk*df(i)*dj_cos_ikj(:)
          ! d/dR_k (w_ijk*f)
          deriv_k(:,i) = dk_w_ijk(:)*f(i) + w_ijk*df(i)*dk_cos_ikj(:)
       end forall

       if(present(dstress))then
          dstress = 0.d0
          djpart1 = dfc_j*fc_k*R_ij/d_ij
          dkpart1 = fc_j*dfc_k*R_ik/d_ik
          djpart2 = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
          !djpart2 = -cos_ijk*R_ij*id_ij2 + R_ij*id_ij_ik
          dkpart2 = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik
          !dkpart2 = -cos_ijk*R_ik*id_ik2 + R_ik*id_ij_ik
          forall (i=1:sfb%a_N)
             !dstress(1,i) =-(djpart1(1)*f(i) + w_ijk*df(i)*djpart2(1))*R_ij(1) &
             !              -(dkpart1(1)*f(i) + w_ijk*df(i)*dkpart2(1))*R_ik(1)
             !dstress(2,i) =-(djpart1(2)*f(i) + w_ijk*df(i)*djpart2(2))*R_ij(2) &
             !              -(dkpart1(2)*f(i) + w_ijk*df(i)*dkpart2(2))*R_ik(2)
             !dstress(3,i) =-(djpart1(3)*f(i) + w_ijk*df(i)*djpart2(3))*R_ij(3) &
             !              -(dkpart1(3)*f(i) + w_ijk*df(i)*dkpart2(3))*R_ik(3)
             !dstress(4,i) =-(djpart1(2)*f(i) + w_ijk*df(i)*djpart2(2))*R_ij(3) &
             !              -(dkpart1(2)*f(i) + w_ijk*df(i)*dkpart2(2))*R_ik(3)
             !dstress(5,i) =-(djpart1(3)*f(i) + w_ijk*df(i)*djpart2(3))*R_ij(1) &
             !              -(dkpart1(3)*f(i) + w_ijk*df(i)*dkpart2(3))*R_ik(1)
             !dstress(6,i) =-(djpart1(1)*f(i) + w_ijk*df(i)*djpart2(1))*R_ij(2) &
             !              -(dkpart1(1)*f(i) + w_ijk*df(i)*dkpart2(1))*R_ik(2)
             dstress(1,i) = (djpart1(1)*f(i) + w_ijk*df(i)*djpart2(1))*R_ij(1) &
                           +(dkpart1(1)*f(i) + w_ijk*df(i)*dkpart2(1))*R_ik(1)
             dstress(2,i) = (djpart1(2)*f(i) + w_ijk*df(i)*djpart2(2))*R_ij(2) &
                           +(dkpart1(2)*f(i) + w_ijk*df(i)*dkpart2(2))*R_ik(2)
             dstress(3,i) = (djpart1(3)*f(i) + w_ijk*df(i)*djpart2(3))*R_ij(3) &
                           +(dkpart1(3)*f(i) + w_ijk*df(i)*dkpart2(3))*R_ik(3)
             dstress(4,i) = (djpart1(2)*f(i) + w_ijk*df(i)*djpart2(2))*R_ij(3) &
                           +(dkpart1(2)*f(i) + w_ijk*df(i)*dkpart2(2))*R_ik(3)
             dstress(5,i) = (djpart1(3)*f(i) + w_ijk*df(i)*djpart2(3))*R_ij(1) &
                           +(dkpart1(3)*f(i) + w_ijk*df(i)*dkpart2(3))*R_ik(1)
             dstress(6,i) = (djpart1(1)*f(i) + w_ijk*df(i)*djpart2(1))*R_ij(2) &
                           +(dkpart1(1)*f(i) + w_ijk*df(i)*dkpart2(1))*R_ik(2)
          end forall 

          !dstress = 0.d0

       end if

    end if

  end subroutine sfb_angular_stress

  subroutine sfb_angular(sfb, R_ij, R_ik, d_ij, d_ik, cos_ijk, values, &
                         deriv_i, deriv_j, deriv_k)

    implicit none

    type(FingerprintBasis),                     intent(inout) :: sfb
    double precision, dimension(3),             intent(in)    :: R_ij, R_ik
    double precision,                           intent(in)    :: d_ij, d_ik
    double precision,                           intent(in)    :: cos_ijk
    double precision, dimension(:),             intent(out)   :: values
    double precision, dimension(:,:), optional, intent(out)   :: deriv_i
    double precision, dimension(:,:), optional, intent(out)   :: deriv_j
    double precision, dimension(:,:), optional, intent(out)   :: deriv_k

    double precision                     :: w_ijk
    double precision                     :: fc_j, dfc_j, fc_k, dfc_k
    double precision, dimension(sfb%a_N) :: f, df
    double precision                     :: id_ij2, id_ik2, id_ij_ik
    double precision, dimension(3)       :: di_cos_ikj, dj_cos_ikj, dk_cos_ikj
    double precision, dimension(3)       :: di_w_ijk, dj_w_ijk, dk_w_ijk
    integer                              :: i

    call sfb_assert_init(sfb)

    fc_j = sfb_fc(d_ij, sfb%a_Rc)
    fc_k = sfb_fc(d_ik, sfb%a_Rc)
    w_ijk = fc_j*fc_k
    f = chebyshev_polynomial(cos_ijk, 0.0d0, PI, sfb%a_order)

    values(1:sfb%a_N) = w_ijk*f

    if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
       dfc_j = sfb_fc_d1(d_ij, sfb%a_Rc)
       dfc_k = sfb_fc_d1(d_ik, sfb%a_Rc)
       df = chebyshev_polynomial_d1(cos_ijk, 0.0d0, PI, sfb%a_order)
       id_ij2 = 1.0d0/(d_ij*d_ij)
       id_ik2 = 1.0d0/(d_ik*d_ik)
       id_ij_ik = 1.0d0/(d_ij*d_ik)
       ! d/dR_i (cos_ijk)
       di_cos_ikj = cos_ijk*(R_ij*id_ij2 + R_ik*id_ik2) - (R_ij+R_ik)*id_ij_ik
       ! d/dR_j (cos_ijk)
       dj_cos_ikj = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
       ! d/dR_k (cos_ijk)
       dk_cos_ikj = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik
       ! d/dR_i (w_ijk)
       di_w_ijk = -(dfc_j*fc_k*R_ij/d_ij + fc_j*dfc_k*R_ik/d_ik)
       ! d/dR_j (w_ijk)
       dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
       ! d/dR_k (w_ijk)
       dk_w_ijk = fc_j*dfc_k*R_ik/d_ik
       forall (i=1:sfb%a_N)
          ! d/dR_i (w_ijk*f)
          deriv_i(:,i) = di_w_ijk(:)*f(i) + w_ijk*df(i)*di_cos_ikj(:)
          ! d/dR_j (w_ijk*f)
          deriv_j(:,i) = dj_w_ijk(:)*f(i) + w_ijk*df(i)*dj_cos_ikj(:)
          ! d/dR_k (w_ijk*f)
          deriv_k(:,i) = dk_w_ijk(:)*f(i) + w_ijk*df(i)*dk_cos_ikj(:)
       end forall
    end if

  end subroutine sfb_angular

end module sfbasis
!shimamura




module symmfunc

  implicit none

  double precision, parameter, private :: PI    = 3.14159265358979d0
  double precision, parameter, private :: PI2   = 2.0d0*PI
  double precision, parameter, private :: EPS   = 1.0d-12
  integer,          parameter, private :: NG    = 5

  !--------------------------------------------------------------------!
  ! sf_nG_type(i)    : number of symmetry functions for species i      !
  !                                                                    !
  ! sf_nTypes        : number of different atom/point types            !
  ! sf_nPairs        : number type pairs (= nTypes!)                   !
  ! sf_nGmax         : max. num. symmetry functs. per type pair/triple !
  ! sf_idx2(i,j)     : pair index for species i-j                      !
  ! sf_idx3(i,j,k)   : triple index for species i-j-k                  !
  ! sf_nG(i,pij)     : num. symmetry function type i for type pair pij !
  ! sf_pG?(i,iG,pij) : parameter i for the iG-th symm.-func. of type ? !
  !                    for type pair pij                               !
  !                    Parameters:  G1:  Rc                            !
  !                                 G2:  Rc, Rs, eta                   !
  !                                 G3:  Rc, kappa                     !
  !                                 G4:  Rc, lambda, zeta, eta         !
  !                                 G5:  Rc, lambda, zeta, eta         !
  !--------------------------------------------------------------------!

  integer,          dimension(:),     allocatable, public  :: sf_nG_type

  integer,                                         private :: sf_nTypes
  integer,                                         private :: sf_nPairs
  integer,                                         private :: sf_nTriples
  integer,                                         private :: sf_nGMax
  integer,          dimension(:,:),   allocatable, private :: sf_idx2
  integer,          dimension(:,:,:), allocatable, private :: sf_idx3
  integer,          dimension(:),     allocatable, private :: sf_iG02
  integer,          dimension(:),     allocatable, private :: sf_iG03
  integer,          dimension(:,:),   allocatable, private :: sf_nG2
  integer,          dimension(:,:),   allocatable, private :: sf_nG3
  double precision, dimension(:,:),   allocatable, private :: sf_pG1
  double precision, dimension(:,:,:), allocatable, private :: sf_pG2
  double precision, dimension(:,:,:), allocatable, private :: sf_pG3
  double precision, dimension(:,:,:), allocatable, private :: sf_pG4
  double precision, dimension(:,:,:), allocatable, private :: sf_pG5

  !--------------------------------------------------------------------!

  integer, public :: stdout = 6
  integer, public :: stderr = 0

  !--------------------------------------------------------------------!

  logical, private :: isInit  = .false.
  integer, private :: memSize = 0

  !--------------------------------------------------------------------!

  integer, parameter :: norder = 1
  integer :: nbin = 100000 ! size of table
!  integer, parameter :: nbin = 10000 ! for the second-order interpolation
  real(8), allocatable, private :: cutv(:), cutva(:), cutd(:), cutda(:)
  real(8) :: dx, dxr

  real(8), allocatable, private :: exv(:), exva(:)
  real(8) :: exx, dex, dexr
  save

contains !-------------------------------------------------------------!

  subroutine sf_init(ntypes, nGmax)

    implicit none

    integer, intent(in) :: ntypes
    integer, intent(in) :: nGmax

    integer :: i, j, k, idx2, idx3

    sf_nTypes = ntypes
    sf_nGMax  = nGmax

    ! initialize type pair/triple index:
    ! note: for angular symmetry functions A-B-C = A-C-B
    allocate(sf_idx2(ntypes,ntypes), sf_idx3(ntypes,ntypes,ntypes))
    idx2 = 0
    idx3 = 0
    do i = 1, ntypes
       do j = 1, ntypes
          idx2 = idx2 + 1
          sf_idx2(i,j) = idx2
          idx3 = idx3 + 1
          sf_idx3(i,j,j) = idx3
          do k = j+1, ntypes
             idx3 = idx3 + 1
             sf_idx3(i,j,k) = idx3
             sf_idx3(i,k,j) = idx3
          end do
       end do
    end do
    sf_nPairs = idx2
    sf_nTriples = idx3

    allocate(sf_nG_type(ntypes),            &
             sf_nG2(NG, sf_nPairs),         &
             sf_nG3(NG, sf_nTriples),       &
             sf_iG02(sf_nPairs),            &
             sf_iG03(sf_nTriples),          &
             sf_pG1(nGmax, sf_nPairs),      &
             sf_pG2(3, nGmax, sf_nPairs),   &
             sf_pG3(2, nGmax, sf_nPairs),   &
             sf_pG4(4, nGmax, sf_nTriples), &
             sf_pG5(4, nGmax, sf_nTriples)  )
    sf_nG_type(:) = 0
    sf_nG2(:,:) = 0
    sf_nG3(:,:) = 0

    isInit = .true.

  end subroutine sf_init

  !--------------------------------------------------------------------!
  !                      parameter initialization                      !
  !--------------------------------------------------------------------!

  subroutine sf_add_rad(funct, type1, type2, Rc, Rs, eta, kappa)

    implicit none

    integer,                    intent(in) :: funct
    integer,                    intent(in) :: type1, type2
    double precision, optional, intent(in) :: Rc, Rs, eta, kappa

    integer :: iG, pij

    pij = sf_idx2(type1,type2)

    if (sf_nG2(funct,pij) >= sf_nGMax) then
       write(stderr,*) "Error: Memory exceeded in `sf_add()'."
       write(stderr,*) sf_nG2(funct,pij), sf_nGMax
       stop
    end if

    if (.not. present(Rc)) then
       write(stderr,*) "Error: Rc undefined in `sf_add()'"
       stop
    end if

    iG = sf_nG2(funct,pij) + 1
    select case(funct)
    case(1)
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG1(iG,pij)   = Rc
    case(2)
       if (.not. (present(Rs) .and. present(eta))) then
          write(stderr,*) "Error: needed for G2: Rs, eta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG2(1,iG,pij) = Rc
       sf_pG2(2,iG,pij) = Rs
       sf_pG2(3,iG,pij) = eta
    case(3)
       if (.not. present(kappa)) then
          write(stderr,*) "Error: needed for G3: kappa"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG3(1,iG,pij) = Rc
       sf_pG3(2,iG,pij) = kappa
    case default
       write(stderr,*) "Error: Not a radial symmetry function type: ", funct
       stop
    end select

    call sf_build_index()

  end subroutine sf_add_rad

  !--------------------------------------------------------------------!

  subroutine sf_add_ang(funct, type1, type2, type3, Rc, eta, lambda, zeta)

    implicit none

    integer,                    intent(in) :: funct
    integer,                    intent(in) :: type1, type2, type3
    double precision, optional, intent(in) :: Rc, eta, lambda, zeta

    integer :: iG, tijk

    tijk = sf_idx3(type1, type2, type3)

    if (sf_nG3(funct,tijk) >= sf_nGMax) then
       write(stderr,*) "Error: Memory exceeded in `sf_add()'."
       write(stderr,*) sf_nG3(funct,tijk), sf_nGMax
       stop
    end if

    if (.not. present(Rc)) then
       write(stderr,*) "Error: Rc undefined in `sf_add()'"
       stop
    end if

    iG = sf_nG3(funct,tijk) + 1
    select case(funct)
    case(4)
       if (.not. (present(lambda) .and. present(zeta) &
           .and. present(eta))) then
          write(stderr,*) "Error: needed for G4: eta, lambda, zeta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG3(funct,tijk) = iG
       sf_pG4(1,iG,tijk) = Rc
       sf_pG4(2,iG,tijk) = lambda
       sf_pG4(3,iG,tijk) = zeta
       sf_pG4(4,iG,tijk) = eta
    case(5)
       if (.not. (present(lambda) .and. present(zeta) &
           .and. present(eta))) then
          write(stderr,*) "Error: needed for G5: eta, lambda, zeta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG3(funct,tijk) = iG
       sf_pG5(1,iG,tijk) = Rc
       sf_pG5(2,iG,tijk) = lambda
       sf_pG5(3,iG,tijk) = zeta
       sf_pG5(4,iG,tijk) = eta
    case default
       write(stderr,*) "Error: Not an angular symmetry function type: ", funct
       stop
    end select

    call sf_build_index()

  end subroutine sf_add_ang

  !--------------------------------------------------------------------!
  ! The functions are organized in memory not in the same order as     !
  ! they are added, but rather depending on their type and species.    !
  ! The procedure `sf_build_index()' generates the indices for radial  !
  ! and angular functions.                                             !
  !--------------------------------------------------------------------!

  subroutine sf_build_index()

    implicit none

    integer :: pij, tijk, i
    integer :: itype1, itype2, itype3

    do itype1 = 1, sf_nTypes
       i = 0
       do itype2 = 1, sf_nTypes
          pij = sf_idx2(itype1, itype2)
          sf_iG02(pij) = i
          i = i + sum(sf_nG2(:,pij))
          do itype3 = itype2, sf_nTypes
             tijk = sf_idx3(itype1, itype2, itype3)
             sf_iG03(tijk) = i
             i = i + sum(sf_nG3(:,tijk))
          end do
       end do
    end do

  end subroutine sf_build_index

  !--------------------------------------------------------------------!
  !               evaluate finger print for single atom                !
  !--------------------------------------------------------------------!

  subroutine sf_fingerprint(type1, Xi, nx, X, typej, n, G, dGi, dGj, dGh)

    implicit none

    !------------------------------------------------------------------!
    ! type1     : type index of the central particle                   !
    ! Xi(1:3)   : cartesian coordinates of the central particle        !
    ! nx        : number of other particles                            !
    ! X(1:3,j)  : cartesian coordinates of the j-th particle           !
    ! typej(j)  : atom type index of j-th atom                         !
    ! n         : number of symmetry functions                         !
    ! G(j)      : (out) value of the j-th symmetry function            !
    ! dGi(1:3,j): (out, optional) derivative of the j-th symmetry      !
    !             function wrt. the central particle's coordinates     !
    ! dGj(:,j,k): (out, optional) derivative of the j-th symmetry      !
    !             function wrt. the coordinates of the k-th particle   !
    ! Note: dGi and dGj have to be present/absent simultaneously       !
    !------------------------------------------------------------------!

    integer,                                       intent(in)  :: type1
    double precision, dimension(3),                intent(in)  :: Xi
    integer,                                       intent(in)  :: nx
    double precision, dimension(3,nx),             intent(in)  :: X
    integer,          dimension(nx),               intent(in)  :: typej
    integer,                                       intent(in)  :: n
    double precision, dimension(n),                intent(out) :: G
    double precision, dimension(3,n),    optional, intent(out) :: dGi
    double precision, dimension(3,n,nx), optional, intent(out) :: dGj
    double precision, dimension(6,n),    optional, intent(out) :: dGh

    integer          :: j, k, pij, tijk
    integer          :: iG_ang0, iG
    double precision :: Rij, Rik, Rjk, cost
    logical          :: deriv, lstrs

    double precision, dimension(3) :: R1, R2, R3

    lstrs = .false.
    if (present(dGi) .and. present(dGj)) then
       deriv = .true.
       if (present(dGh)) then
           lstrs = .true.
       end if
    else
       deriv = .false.
    end if

    if (.not. isInit) return

    if (n < sf_nG_type(type1)) then
       write(stderr,*) "Error: number of symmetry functions exceeded."
       stop
    end if

    G(1:n)  = 0.0d0
    iG_ang0 = 0
    if (deriv) then
       dGi(1:3,1:n) = 0.0d0
       dGj(1:3,1:n,1:nx) = 0.0d0
       if (lstrs) then
          dGh(1:6,1:n) = 0.0d0
       end if
    end if
    jloop : do j = 1, nx
       pij = sf_idx2(type1, typej(j))

       R1  = X(:,j) - Xi(:)
       Rij = sqrt(sum(R1*R1))
       R1  = R1/Rij

       iG = sf_iG02(pij)

       ! radial symmetry functions:
       if (deriv) then
          if (lstrs) then
             call sf_G1_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j), dGh)
             call sf_G2_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j), dGh)
             call sf_G3_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j), dGh)
          else
             call sf_G1_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
             call sf_G2_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
             call sf_G3_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
          end if
       else
          call sf_G1_update(pij, R1(1:3), Rij, n, G, iG)
          call sf_G2_update(pij, R1(1:3), Rij, n, G, iG)
          call sf_G3_update(pij, R1(1:3), Rij, n, G, iG)
       end if

       ! if (iG_ang0 == 0) iG_ang0 = sum(sf_nG2(:,:))
       kloop : do k = j+1, nx
          tijk = sf_idx3(type1, typej(j), typej(k))

          R2  = X(:,k) - Xi(:)
          Rik = sqrt(sum(R2*R2))
          if (Rik <= 1.0d-8) then
             write(stderr, *) "Warning: an atom in the neighbor list is " &
                  // "identical to the central atom (sf_fingerprint)"
          else
             R2  = R2/Rik
          end if

          R3  = X(:,k) - X(:,j)
          Rjk = sqrt(sum(R3*R3))
          if (Rjk <= 1.0d-8) then
             write(stderr, *) "Warning: redundant atoms in neighbor " &
                              // "list (sf_fingerprint)"
          else
             R3  = R3/Rjk
          end if

          ! cos(theta_ijk)
          cost = sum(R1*R2)
          cost = max(cost,-1.0d0)
          cost = min(cost, 1.0d0)

          iG = sf_iG03(tijk)

          ! angular symmetry functions:
          if (deriv) then
             if (lstrs) then
                call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k), dGh)
                call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k), dGh)
             else
                call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k))
                call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                                  n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k))
             end if
          else
             call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, &
                               cost, n, G, iG)
             call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                               n, G, iG)
          end if

       end do kloop
    end do jloop

  end subroutine sf_fingerprint



  !==================== radial symmetry functions =====================!



  !--------------------------------------------------------------------!
  !            radial symm. function 1 (eq. 5 in Ref. [1])             !
  !                                                                    !
  ! here: vecRij = vec{R}_ij/||vec{R})_ij||                            !
  !--------------------------------------------------------------------!

  subroutine sf_G1_ij(Rij, Rc, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc, Rcr

    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

    Gij  = fc
    dGij = dfc

  end subroutine sf_G1_ij

  !--------------------------------------------------------------------!

  subroutine sf_G1_update(pij, vecRij, Rij, n, G, iG, dGi, dGj, dGh)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer          :: iG1
    double precision :: G1, dG1, dG1r
    double precision :: Rc

    do iG1 = 1, sf_nG2(1,pij)
       iG    = iG + 1
       Rc    = sf_pG1(iG1,pij)
       call sf_G1_ij(Rij, Rc, G1, dG1)
       G(iG) = G(iG) + G1
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG1*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG1*vecRij(1:3)
       if (present(dGh)) then
           dG1r = dG1*Rij
           dGh(1,iG) = dGh(1,iG) + dG1r*vecRij(1)*vecRij(1)
           dGh(2,iG) = dGh(2,iG) + dG1r*vecRij(2)*vecRij(2)
           dGh(3,iG) = dGh(3,iG) + dG1r*vecRij(3)*vecRij(3)
           dGh(4,iG) = dGh(4,iG) + dG1r*vecRij(2)*vecRij(3)
           dGh(5,iG) = dGh(5,iG) + dG1r*vecRij(3)*vecRij(1)
           dGh(6,iG) = dGh(6,iG) + dG1r*vecRij(1)*vecRij(2)
       end if
    end do

  end subroutine sf_G1_update

  !--------------------------------------------------------------------!
  !            radial symm. function 2 (eq. 6 in Ref. [1])             !
  !--------------------------------------------------------------------!

  subroutine sf_G2_ij(Rij, Rc, Rs, eta, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: Rs
    double precision, intent(in)  :: eta
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc, Rcr
    double precision :: fexp, arg, x, d, etag
    integer :: m

    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

    arg  = Rij - Rs
!    fexp = exp(-eta*arg*arg)
    etag = eta*arg
    x = etag*arg
    if( x > exx ) then
        fexp = 0.d0
    else
    if( norder == 1 ) then
        !---first-order interpolation
        d = x*dexr
        m = d
        d = d - m
        fexp = (1d0-d)*exv(m)+d*exv(m+1)
    else
        !---second-order interpolation
        m = x/(2.d0*dex)
        m = 2*m
        d = 0.5d0*( x*dexr - dble(m) )
        fexp = d*( (d-1.d0)*exva(m) + exv(m+2) - exv(m) ) + exv(m)
    end if
    end if

    Gij  = fexp*fc
    dGij = fexp*( dfc - 2.0d0*etag*fc  )

  end subroutine sf_G2_ij

  !--------------------------------------------------------------------!

  subroutine sf_G2_update(pij, vecRij, Rij, n, G, iG, dGi, dGj, dGh)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer          :: iG2
    double precision :: Rc, Rs, eta
    double precision :: G2, dG2, vecRij2(3)

    do iG2 = 1, sf_nG2(2,pij)
       iG    = iG + 1
       Rc    = sf_pG2(1,iG2,pij)
       Rs    = sf_pG2(2,iG2,pij)
       eta   = sf_pG2(3,iG2,pij)
       call sf_G2_ij(Rij, Rc, Rs, eta, G2, dG2)
       G(iG) = G(iG) + G2
       if (present(dGi)) then
           vecRij2(1:3) = dG2*vecRij(1:3)
           dGi(1:3,iG) = dGi(1:3,iG) - vecRij2(1:3)
!       end if
!       if (present(dGj)) then
           dGj(1:3,iG) = dGj(1:3,iG) + vecRij2(1:3)
!       end if
         if (present(dGh)) then
             vecRij2(1:3) = vecRij2(1:3)*Rij
             dGh(1,iG) = dGh(1,iG) + vecRij2(1)*vecRij(1)
             dGh(2,iG) = dGh(2,iG) + vecRij2(2)*vecRij(2)
             dGh(3,iG) = dGh(3,iG) + vecRij2(3)*vecRij(3)
             dGh(4,iG) = dGh(4,iG) + vecRij2(2)*vecRij(3)
             dGh(5,iG) = dGh(5,iG) + vecRij2(3)*vecRij(1)
             dGh(6,iG) = dGh(6,iG) + vecRij2(1)*vecRij(2)
         end if
       end if
    end do

  end subroutine sf_G2_update

  !--------------------------------------------------------------------!
  !            radial symm. function 3 (eq. 6 in Ref. [1])             !
  !--------------------------------------------------------------------!

  subroutine sf_G3_ij(Rij, Rc, kappa, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: kappa
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc, Rcr
    double precision :: fcos, fsin, arg

    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

    arg  = kappa*Rij
    fcos = cos(arg)
    fsin = sin(arg)

    Gij  = fcos*fc
    dGij = fcos*dfc - kappa*fsin*fc

  end subroutine sf_G3_ij

  !--------------------------------------------------------------------!

  subroutine sf_G3_update(pij, vecRij, Rij, n, G, iG, dGi, dGj, dGh)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer          :: iG3
    double precision :: Rc, kappa
    double precision :: G3, dG3, dG3r

    do iG3 = 1, sf_nG2(3,pij)
       iG    = iG + 1
       Rc    = sf_pG3(1,iG3,pij)
       kappa = sf_pG3(2,iG3,pij)
       call sf_G3_ij(Rij, Rc, kappa, G3, dG3)
       G(iG) = G(iG) + G3
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG3*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG3*vecRij(1:3)
       if (present(dGh)) then
           dG3r = dG3*Rij
           dGh(1,iG) = dGh(1,iG) + dG3r*vecRij(1)*vecRij(1)
           dGh(2,iG) = dGh(2,iG) + dG3r*vecRij(2)*vecRij(2)
           dGh(3,iG) = dGh(3,iG) + dG3r*vecRij(3)*vecRij(3)
           dGh(4,iG) = dGh(4,iG) + dG3r*vecRij(2)*vecRij(3)
           dGh(5,iG) = dGh(5,iG) + dG3r*vecRij(3)*vecRij(1)
           dGh(6,iG) = dGh(6,iG) + dG3r*vecRij(1)*vecRij(2)
       end if
    end do

  end subroutine sf_G3_update



  !==================== angular symmetry functions ====================!

  !--------------------------------------------------------------------!
  !          first angular symm. function (eq. 8 in Ref. [1])          !
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !                        angle dependend part                        !
  ! This is function F_1(R_ij,R_ik) in the documentation.              !
  !--------------------------------------------------------------------!

  subroutine sf_F1_ijk(vecRij, vecRik, Rijr, Rikr, cost, lambda, &
                       zeta, Fijk, dFijk_dRj, dFijk_dRk)

    implicit none

    double precision, dimension(3), intent(in)  :: vecRij, vecRik
    double precision,               intent(in)  :: Rijr, Rikr
    double precision,               intent(in)  :: cost
    double precision,               intent(in)  :: lambda
    double precision,               intent(in)  :: zeta
    double precision,               intent(out) :: Fijk
    double precision, dimension(3), intent(out) :: dFijk_dRj
    double precision, dimension(3), intent(out) :: dFijk_dRk

    double precision :: arg, prefactor, argz, prefactor1
    real(8), parameter :: epsi = 1.d-12

    if (abs(cost) > 1.0d0) write(*,*) "cos(theta) = ", cost

!    if( abs(lambda-1.0d0) < epsi ) then
!        arg = 0.5d0*(1.0d0 + cost)
!    else if( abs(lambda+1.0d0) < epsi ) then
!        arg = 0.5d0*(1.0d0 - cost)
!    else
        arg = 0.5d0*(1.0d0 + lambda*cost)
!    end if

!    prefactor = 0.5d0*zeta*lambda*arg**(zeta-1.0d0)
    if( abs(zeta-1.0d0) < epsi ) then
        argz = 1.d0
    else if( abs(zeta-2.0d0) < epsi ) then
        argz = arg
    else if( abs(zeta-3.0d0) < epsi ) then
        argz = arg*arg
    else if( abs(zeta-4.0d0) < epsi ) then
        argz = arg*arg*arg
    else
        argz = arg**(zeta-1.0d0)
    end if
!    if( abs(lambda-1.0d0) < epsi ) then
!        prefactor = 0.5d0*zeta*argz
!    else if( abs(lambda+1.0d0) < epsi ) then
!        prefactor = -0.5d0*zeta*argz
!    else
        prefactor = 0.5d0*zeta*lambda*argz
!    end if

!    Fijk  = arg**zeta
    Fijk  = argz*arg

    prefactor1 = prefactor*Rijr
!    dFijk_dRj(1:3) = prefactor1*( -cost*vecRij(1:3) + vecRik(1:3) )
    dFijk_dRj(1:3) = prefactor1*vecRij(1:3)

    prefactor1 = prefactor*Rikr
!    dFijk_dRk(1:3) = prefactor1*( -cost*vecRik(1:3) + vecRij(1:3) )
    dFijk_dRk(1:3) = prefactor1*vecRik(1:3)

  end subroutine sf_F1_ijk

  !--------------------------------------------------------------------!
  !                      distance dependend part                       !
  ! This is function F_2(R) in the documentation.                      !
  ! --> very similar to sf_G2_ij(), may be combined in future.         !
  !--------------------------------------------------------------------!

  subroutine sf_F2_ij(Rij, Rc, Rcr, eta, Fij, dFij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc, Rcr
    double precision, intent(in)  :: eta
    double precision, intent(out) :: Fij
    double precision, intent(out) :: dFij

    double precision :: fc, dfc
    double precision :: fexp, x, d, etij
    integer :: m

!    Rcr = 1.d0/Rc
    call sf_cut(Rij, Rc, Rcr, fc, dfc)

!    fexp = exp(-eta*Rij*Rij)
    etij = eta*Rij
    x = etij*Rij
    if( x > exx ) then
        fexp = 0.d0
    else
    if( norder == 1 ) then
        !---first-order interpolation
        d = x*dexr
        m = d
        d = d - m
        fexp = (1d0-d)*exv(m)+d*exv(m+1)
    else
        !---second-order interpolation
        m = x/(2.d0*dex)
        m = 2*m
        d = 0.5d0*( x*dexr - dble(m) )
        fexp = d*( (d-1.d0)*exva(m) + exv(m+2) - exv(m) ) + exv(m)
    end if
    end if

    Fij  = fexp*fc
    dFij = fexp*( dfc - 2.0d0*etij*fc  )

  end subroutine sf_F2_ij

  !--------------------------------------------------------------------!

  subroutine sf_G4_update(tijk, vecRij, vecRik, vecRjk, Rij, Rik, Rjk, &
                          cost, n, G, iG, dGi, dGj, dGk, dGh)

    implicit none

    integer,                                    intent(in)    :: tijk
    double precision, dimension(3),             intent(in)    :: vecRij, vecRik, vecRjk
    double precision,                           intent(in)    :: Rij, Rik, Rjk
    double precision,                           intent(in)    :: cost
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(3,n), optional, intent(inout) :: dGk
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer                        :: iG4
    double precision               :: Rc, Rcr, lambda, zeta, eta
    double precision               :: G4, G0
    double precision, dimension(3) :: dG4, dF1j, dF1k

    double precision :: F1, F2ij, dF2ij, F2ik, dF2ik, F2jk, dF2jk
    double precision :: Rijr, Rikr, vecRijik(3), vecRikij(3)
    double precision :: dGhh(6), dGhij(3), dGhik(3), dGhjk(3)
    double precision :: dG401, dG402, dG403, vecRij0(3), vecRik0(3), vecRjk0(3)

    Rijr = 1.d0/Rij
    Rikr = 1.d0/Rik
    vecRijik(1:3) = -cost*vecRij(1:3) + vecRik(1:3)
    vecRikij(1:3) = -cost*vecRik(1:3) + vecRij(1:3)
    do iG4 = 1, sf_nG3(4,tijk)
       iG     = iG + 1
       Rc     = sf_pG4(1,iG4,tijk)
       lambda = sf_pG4(2,iG4,tijk)
       zeta   = sf_pG4(3,iG4,tijk)
       eta    = sf_pG4(4,iG4,tijk)

       call sf_F1_ijk(vecRijik, vecRikij, Rijr, Rikr, cost, lambda, &
                      zeta, F1, dF1j, dF1k)
       Rcr = 1.d0/Rc
       call sf_F2_ij(Rij, Rc, Rcr, eta, F2ij, dF2ij)
       call sf_F2_ij(Rik, Rc, Rcr, eta, F2ik, dF2ik)
       call sf_F2_ij(Rjk, Rc, Rcr, eta, F2jk, dF2jk)

       G0    = 2.d0*F2ij*F2ik*F2jk
       G4    = F1*G0
       ! factor of 2 for k<j in sf_fingerprint()
!       G(iG) = G(iG) + 2.0d0*G4
       G(iG) = G(iG) + G4

       if (present(dGi)) then
           F1   = 2.d0*F1
          dG401 = F1*dF2ij*F2ik*F2jk
          dG402 = F1*F2ij*dF2ik*F2jk
          dG403 = F1*F2ij*F2ik*dF2jk
          vecRij0(1:3) = vecRij(1:3)*dG401
          vecRik0(1:3) = vecRik(1:3)*dG402
          vecRjk0(1:3) = vecRjk(1:3)*dG403
          dF1j(1:3) = dF1j(1:3)*G0
          dF1k(1:3) = dF1k(1:3)*G0
          dG4(1:3) = -(dF1j(1:3) + dF1k(1:3)) - vecRij0(1:3) - vecRik0(1:3)
!          dGi(1:3,iG) = dGi(1:3,iG) + 2.0d0*dG4(1:3)
          dGi(1:3,iG) = dGi(1:3,iG) + dG4(1:3)
!       end if
!       if (present(dGj)) then
          dG4(1:3) =  dF1j(1:3) + vecRij0(1:3) - vecRjk0(1:3)
!          dGj(1:3,iG) = dGj(1:3,iG) + 2.0d0*dG4(1:3)
          dGj(1:3,iG) = dGj(1:3,iG) + dG4(1:3)
!       end if
!       if (present(dGk)) then
          dG4(1:3) =  dF1k(1:3) + vecRik0(1:3) + vecRjk0(1:3)
!          dGk(1:3,iG) = dGk(1:3,iG) + 2.0d0*dG4(1:3)
          dGk(1:3,iG) = dGk(1:3,iG) + dG4(1:3)
!       end if
       if (present(dGh)) then
           dGhij(1:3) = ( vecRij0(1:3) + dF1j(1:3) )*Rij
           dGhik(1:3) = ( vecRik0(1:3) + dF1k(1:3) )*Rik
           dGhjk(1:3) =   vecRjk0(1:3)*Rjk
           dGhh(1) = dGhij(1)*vecRij(1)  &
&                  + dGhik(1)*vecRik(1)  &
&                  + dGhjk(1)*vecRjk(1)
           dGhh(2) = dGhij(2)*vecRij(2)  &
&                  + dGhik(2)*vecRik(2)  &
&                  + dGhjk(2)*vecRjk(2)
           dGhh(3) = dGhij(3)*vecRij(3)  &
&                  + dGhik(3)*vecRik(3)  &
&                  + dGhjk(3)*vecRjk(3)
           dGhh(4) = dGhij(2)*vecRij(3)  &
&                  + dGhik(2)*vecRik(3)  &
&                  + dGhjk(2)*vecRjk(3)
           dGhh(5) = dGhij(3)*vecRij(1)  &
&                  + dGhik(3)*vecRik(1)  &
&                  + dGhjk(3)*vecRjk(1)
           dGhh(6) = dGhij(1)*vecRij(2)  &
&                  + dGhik(1)*vecRik(2)  &
&                  + dGhjk(1)*vecRjk(2)
!           dGh(1:6,iG) = dGh(1:6,iG) + 2.0d0*dGhh(1:6)
           dGh(1:6,iG) = dGh(1:6,iG) + dGhh(1:6)
       end if
       end if

    end do

  end subroutine sf_G4_update

  !--------------------------------------------------------------------!
  !         second angular symm. function (eq. 9 in Ref. [1])          !
  !--------------------------------------------------------------------!

  subroutine sf_G5_update(tijk, vecRij, vecRik, Rij, Rik, cost, n, &
                          G, iG, dGi, dGj, dGk, dGh)

    implicit none

    integer,                                    intent(in)    :: tijk
    double precision, dimension(3),             intent(in)    :: vecRij, vecRik
    double precision,                           intent(in)    :: Rij, Rik
    double precision,                           intent(in)    :: cost
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(3,n), optional, intent(inout) :: dGk
    double precision, dimension(6,n), optional, intent(inout) :: dGh

    integer                        :: iG5
    double precision               :: Rc, Rcr, lambda, zeta, eta
    double precision               :: G5
    double precision, dimension(3) :: dG5, dF1j, dF1k

    double precision :: F1, F2ij, dF2ij, F2ik, dF2ik
    double precision :: Rijr, Rikr, vecRijik(3), vecRikij(3)
    double precision :: dGhh(6), dGhij, dGhik, dF1ij(3), dF1ik(3)

    Rijr = 1.d0/Rij
    Rikr = 1.d0/Rik
    vecRijik(1:3) = -cost*vecRij(1:3) + vecRik(1:3)
    vecRikij(1:3) = -cost*vecRik(1:3) + vecRij(1:3)
    do iG5 = 1, sf_nG3(5,tijk)
       iG     = iG + 1
       Rc     = sf_pG5(1,iG5,tijk)
       lambda = sf_pG5(2,iG5,tijk)
       zeta   = sf_pG5(3,iG5,tijk)
       eta    = sf_pG5(4,iG5,tijk)

       call sf_F1_ijk(vecRijik, vecRikij, Rijr, Rikr, cost, lambda, &
                      zeta, F1, dF1j, dF1k)
       Rcr = 1.d0/Rc
       call sf_F2_ij(Rij, Rc, Rcr, eta, F2ij, dF2ij)
       call sf_F2_ij(Rik, Rc, Rcr, eta, F2ik, dF2ik)

       G5    = F1*F2ij*F2ik
       ! factor of two for k<j in sf_fingerprint()
       G(iG) = G(iG) + 2.0d0*G5

       if (present(dGi)) then
          dG5(1:3) = -(dF1j(1:3)+dF1k(1:3))*F2ij*F2ik        &
                   - vecRij(1:3)*F1*dF2ij*F2ik &
                   - vecRik(1:3)*F1*F2ij*dF2ik
          dGi(1:3,iG) = dGi(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGj)) then
          dG5(1:3) = dF1j(1:3)*F2ij*F2ik        &
                   + vecRij(1:3)*F1*dF2ij*F2ik
          dGj(1:3,iG) = dGj(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGk)) then
          dG5(1:3) = dF1k(1:3)*F2ij*F2ik        &
                   + vecRik(1:3)*F1*F2ij*dF2ik
          dGk(1:3,iG) = dGk(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGh)) then
           dGhij = F1*dF2ij* F2ik*Rij
           dGhik = F1* F2ij*dF2ik*Rik
           dF1ij(1:3) = F2ij*F2ik*dF1j(1:3)*Rij
           dF1ik(1:3) = F2ij*F2ik*dF1k(1:3)*Rik
           dGhh(1) = dGhij*vecRij(1)*vecRij(1)  &
&                  + dGhik*vecRik(1)*vecRik(1)  &
&                  +        dF1ij(1)*vecRij(1)  &
&                  +        dF1ik(1)*vecRik(1)
           dGhh(2) = dGhij*vecRij(2)*vecRij(2)  &
&                  + dGhik*vecRik(2)*vecRik(2)  &
&                  +        dF1ij(2)*vecRij(2)  &
&                  +        dF1ik(2)*vecRik(2)
           dGhh(3) = dGhij*vecRij(3)*vecRij(3)  &
&                  + dGhik*vecRik(3)*vecRik(3)  &
&                  +        dF1ij(3)*vecRij(3)  &
&                  +        dF1ik(3)*vecRik(3)
           dGhh(4) = dGhij*vecRij(2)*vecRij(3)  &
&                  + dGhik*vecRik(2)*vecRik(3)  &
&                  +        dF1ij(2)*vecRij(3)  &
&                  +        dF1ik(2)*vecRik(3)
           dGhh(5) = dGhij*vecRij(3)*vecRij(1)  &
&                  + dGhik*vecRik(3)*vecRik(1)  &
&                  +        dF1ij(3)*vecRij(1)  &
&                  +        dF1ik(3)*vecRik(1)
           dGhh(6) = dGhij*vecRij(1)*vecRij(2)  &
&                  + dGhik*vecRik(1)*vecRik(2)  &
&                  +        dF1ij(1)*vecRij(2)  &
&                  +        dF1ik(1)*vecRik(2)
           dGh(1:6,iG) = dGh(1:6,iG) + 2.0d0*dGhh(1:6)
       end if

    end do

  end subroutine sf_G5_update


  !======================= auxiliary functions ========================!



  !--------------------------------------------------------------------!
  !                      cosine cut-off function                       !
  !                        (eq. 4 in Ref. [1])                         !
  !--------------------------------------------------------------------!

  subroutine sf_cut(Rij, Rc, Rcr, fc, dfc)

    implicit none

    double precision, intent(in)  :: Rij, Rc, Rcr
    double precision, intent(out) :: fc, dfc
    real(8) :: x, d
    integer :: m

    if (Rij >= Rc) then
       fc  = 0.0d0
       dfc = 0.0d0
       return
    end if

    fc  =  0.5d0*(cos(PI*Rij/Rc) + 1.0d0)
    dfc = -0.5d0*PI/Rc*sin(PI*Rij/Rc)

    return  ! kn


    x = PI*Rij*Rcr
    if( norder == 1 ) then
        !---first-order interpolation
        d = x*dxr
        m = d
        d = d - m
         fc =  (1d0-d)*cutv(m)+d*cutv(m+1)
        dfc = ((1d0-d)*cutd(m)+d*cutd(m+1))*Rcr
    else
        !---second-order interpolation
        m = x/(2.d0*dx)
        m = 2*m
        d = 0.5d0*( x*dxr - dble(m) )
         fc =  d*( (d-1.d0)*cutva(m) + cutv(m+2) - cutv(m) ) + cutv(m)
        dfc = (d*( (d-1.d0)*cutda(m) + cutd(m+2) - cutd(m) ) + cutd(m))*Rcr
    end if

  end subroutine sf_cut


  subroutine sf_set_table

    implicit none
    integer :: k
    real(8) :: x

    if( norder == 1 ) then
        nbin = 1000000
    else
        nbin = 10000
    end if

    allocate( cutv(0:nbin), cutd(0:nbin), exv(0:nbin) )

!    pi = acos(-1.d0)

    !---cos, sin for sf_cut
    dx  = pi/(nbin-1)
    dxr = 1.d0/dx

     cutv(:) = 0.d0
     cutd(:) = 0.d0
    do k = 0, nbin
       x = dble(k)*dx
       cutv(k) = 0.5d0*(cos(x) + 1.d0)
       cutd(k) = -0.5d0*pi*sin(x)
    end do

    !---exp
    exx = 40.d0
    dex  = exx/(nbin-1)
    dexr = 1.d0/dex

    exv(:) = 0.d0
    do k = 0, nbin
       x = dble(k)*dex
       exv(k) = exp(-x)
    end do


    if( norder /= 1 ) then
        !---for the second-order interpolation
        allocate( cutva(0:nbin), cutda(0:nbin), exva(0:nbin) )

        cutva(:) = 0.d0
        cutda(:) = 0.d0
         exva(:) = 0.d0
        do k = 0, nbin - 2
           cutva(k) = 2.d0*( cutv(k) + cutv(k+2) - 2.d0*cutv(k+1) )
           cutda(k) = 2.d0*( cutd(k) + cutd(k+2) - 2.d0*cutd(k+1) )
            exva(k) = 2.d0*(  exv(k) +  exv(k+2) - 2.d0* exv(k+1) )
        end do
    end if


  end subroutine sf_set_table

end module symmfunc




module lclist

  implicit none
  save

  public  :: lcl_nmax_nbdist!,        &

contains

  !--------------------------------------------------------------------!
  !          max. number of (real) neighbours within cut-off           !
  !--------------------------------------------------------------------!

  function lcl_nmax_nbdist(rmin, rmax) result(nmax)

    implicit none

    double precision, intent(in) :: rmin, rmax
    integer                      :: nmax

    double precision :: V_atom, V_cut

    V_atom = (0.5d0*rmin)**3
    V_cut  = (rmax+0.5d0*rmin)**3

    ! max number of atoms assuming close packing
    ! pi/(s*sqrt(2)) ~ 0.7405
    nmax = ceiling(V_cut/V_atom*0.7405d0)

  end function lcl_nmax_nbdist

end module lclist




module sfsetup

  use aeio,     only: &!aeio_readline, &
                      TYPELEN

  use io,       only: io_lower,      &
                      io_adjustl, io_adjustl_sub,  &
&                     nfile, myid, nodes

  !shimamura
  use sfbasis,  only: FingerprintBasis, &
                      new_SFBasis,      &
                      del_SFBasis,      &
                      sfb_eval

  use symmfunc, only: sf_init,       &
                      sf_add_rad,    &
                      sf_add_ang,    &
                      sf_fingerprint

  implicit none
  private
  save

  public  :: load_Setup,            &
             load_Setup_ASCII,      &
             stp_init,              &
             stp_nsf_max,           &
             stp_eval,              &
             stp_print_info!,        &

  !shimamura
  !private :: setup_symmfunc_Behler2011!,  &
  private :: setup_symmfunc_Behler2011,  &
             setup_basis_chebyshev

  !--------------------------------------------------------------------!
  !                    structural fingerprint setup                    !
  !--------------------------------------------------------------------!

  type, public :: Setup

     !-----------------------------------------------------------------!
     ! init          .true., if the setup has been initialized         !
     ! neval         number of evaluations                             !
     ! description   an optional description from the setup file       !
     ! atomtype      species of the central atom                       !
     ! nenv          number of different surrounding species           !
     ! envtypes      species of the surrounding atoms                  !
     !                                                                 !
     ! The global atom type index is determined by the order of type   !
     ! names in the input file.  The *local* type index is given by    !
     ! the order the basis functions were first set up.  For the       !
     ! evaluation of the basis functions, the types from the input     !
     ! structure (global index) need to be assigned to local type IDs. !
     !                                                                 !
     ! gtype(i)      global atom type ID for local type i              !
     ! ltype(i)      local atom type ID for global type i              !
     !                                                                 !
     ! Rc_min/max    minimal interaction radius and max. cutoff        !
     ! sftype        basis function type (e.g. Behler2011)             !
     ! nsf           number of structural fingerprint basis functions  !
     ! sf(i)         function kind of the particular basis type        !
     ! nsfparam      the max. number of parameters of a basis function !
     ! sfparam(i,j)  i-th parameter of the j-th basis function         !
     !               i <= nsfparam                                     !
     ! sfenv(i,j)    i-th environment species for j-th basis function  !
     !                                                                 !
     ! sfval_min(i)  lowest so far encountered value of the i-th SF    !
     ! sfval_max(i)  largest so far encountered value of the i-th SF   !
     ! sfval_avg(i)  current average value of the i-th symm. function  !
     ! sfval_cov(i)  current covariance of the i-th symm. function     !
     ! --> min/max/avg/cov are updated whenever UNSCALED evaluation    !
     !     is requested.  This is usually during the screening of the  !
     !     training set.  During the training and prediciton a scaled  !
     !     value of the SFs is useful, that lies within [-1:1].        !
     !                                                                 !
     ! The scaling is implemented as                                   !
     !                                                                 !
     !  f(i) = max(sval_avg(i)-sfval_min(i),sfval_max(i)-sval_avg(i))  !
     !  s(i) = sfval_min(i)                                            !
     !  sfval(i) = (sfval(i)-s(i))/f(i)                                !
     !                                                                 !
     !-----------------------------------------------------------------!

     logical                                             :: init
     integer                                             :: neval

     character(len=1024)                                 :: description

     character(len=TYPELEN)                              :: atomtype
     integer                                             :: nenv
     character(len=TYPELEN), dimension(:),   allocatable :: envtypes

     integer                                             :: ntypes_global
     integer,                dimension(:),   allocatable :: gtype
     integer,                dimension(:),   allocatable :: ltype

     double precision                                    :: Rc_min
     double precision                                    :: Rc_max

     character(len=100)                                  :: sftype
     integer                                             :: nsf
     integer,                dimension(:),   allocatable :: sf
     integer                                             :: nsfparam
     double precision,       dimension(:,:), allocatable :: sfparam
     integer,                dimension(:,:), allocatable :: sfenv

     double precision,       dimension(:),   allocatable :: sfval_min
     double precision,       dimension(:),   allocatable :: sfval_max
     double precision,       dimension(:),   allocatable :: sfval_avg
     double precision,       dimension(:),   allocatable :: sfval_cov

  end type Setup

  !--------------- memory for basis function operations ---------------!
  ! sfval(i)         value of the i-th basis function                  !
  ! sfderiv_i(i,j)   i-th component of the derivative of the j-th SF   !
  !                  with respect to the central atom                  !
  !                  sfderiv_i(3,nsf_max)                              !
  ! sfderiv_j(i,j,k) i-th component of the derivative of the j-th SF   !
  !                  with respect to the coordinates of atom k         !
  !                  sfderiv_j(3,nsf_max,nnb_max)                      !
  !--------------------------------------------------------------------!

  double precision, dimension(:),     allocatable, public :: sfval
  double precision, dimension(:,:),   allocatable, public :: sfderiv_i
  double precision, dimension(:,:,:), allocatable, public :: sfderiv_j
  double precision, dimension(:,:),   allocatable, public :: sfstrs
  integer,                                         public :: nsf_max
  integer,                                         public :: nnb_max

  !shimamura
  !------------------------- Chebyshev basis --------------------------!
  ! sfb(i)      structural fingerprint basis of atom type i            !
  !--------------------------------------------------------------------!

  type(FingerprintBasis), dimension(:), allocatable, private :: sfb

  !--------------------------------------------------------------------!

  logical, private :: isInit = .false.

  !---------------------------- constants -----------------------------!
  ! NSFPARAM    maximum  number of SF parameters                       !
  ! NENV_MAX    maximum number of types involved in single function    !
  !             (e.g., 1 = distance, 2 = angle, 3 = dihedral)          !
  !--------------------------------------------------------------------!

  integer, parameter :: NSFPARAM = 4
  integer, parameter :: NENV_MAX = 2

  integer                            :: nx = 0
  integer, dimension(:), allocatable :: type1_loc

  character(len=50) :: str

contains

  !--------------------------------------------------------------------!
  !           print out information about a specific set-up            !
  !--------------------------------------------------------------------!

  subroutine stp_print_info(stp)

    implicit none

    type(Setup), intent(in) :: stp

    integer :: i, nf

    call stp_assert_init(stp)

 do nf = 1, 2

    write(nfile(nf),*) 'Structural fingerprint (SF) set-up for ', &
         trim(adjustl(stp%atomtype))
    write(nfile(nf),*)

    call stp_print_descr(stp%description, nfile(nf))

    if (stp%nenv>0) then
       write(nfile(nf),'(1x,"environment types: ")', advance='no')
       do i = 1, stp%nenv
          write(nfile(nf),'(A2,1x)', advance='no') stp%envtypes(i)
       end do
    end if
    write(nfile(nf),*)
!    write(*,*) 'minimal distance : ', trim(io_adjustl(stp%Rc_min,2)), ' Angstrom'
!    write(*,*) 'maximal cut-off  : ', trim(io_adjustl(stp%Rc_max,2)), ' Angstrom'
!    write(*,*) 'size of basis    : ', trim(io_adjustl(stp%nsf))
!    write(*,*) 'evaluations      : ', trim(io_adjustl(stp%neval))
    call io_adjustl_sub(stp%Rc_min,str,2)
    write(nfile(nf),*) 'minimal distance : ', trim(str), ' Angstrom'
    call io_adjustl_sub(stp%Rc_max,str,2)
    write(nfile(nf),*) 'maximal cut-off  : ', trim(str), ' Angstrom'
    call io_adjustl_sub(stp%nsf,str)
    write(nfile(nf),*) 'size of basis    : ', trim(str)
    call io_adjustl_sub(stp%neval,str)
    write(nfile(nf),*) 'evaluations      : ', trim(str)
    write(nfile(nf),*)

    select case(trim(io_lower(stp%sftype)))

    !shimamura
    case('chebyshev')
       call print_info_chebyshev(stp, nfile(nf))

    case('behler2011')
       call print_info_Behler2011(stp, nfile(nf))
    end select

 end do


  end subroutine stp_print_info

  !--------------------------------------------------------------------!

  subroutine stp_print_descr(descr, unit)

    implicit none

    character(len=*), intent(in) :: descr
    integer,          intent(in) :: unit

    integer :: i1, i2, l

    l = len_trim(descr)

    i1 = 1
    i2 = index(descr,'$$')
    do while(i2 >= i1)
       write(unit,*) descr(i1:i2-1)
       i1 = i2 + 2
       i2 = i1 + index(descr(i1:l),'$$') - 1
    end do
    if (len_trim(descr(i1:l)) > 0) then
       write(unit,*) trim(descr(i1:l))
    end if
    write(unit,*)

  end subroutine stp_print_descr

  !--------------------------------------------------------------------!

  !shimamura for s-ANN
  function load_Setup(global_types, efile, file, unit) result(stp)
  !function load_Setup(global_types, file, unit) result(stp)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    integer,          optional,     intent(in) :: unit
    type(Setup)                                :: stp

    integer :: u

    !shimamura for s-ANN
    logical,                        intent(in) :: efile


    if (present(unit)) then
       u = unit
!    else if (present(file)) then
!       u = io_unit()
!       open(u, file=trim(file), status='old', action='read', &
!            form='unformatted')
    else
       write(0,*) "Error: neither unit number nor file name given"
       write(0,*) "       in `load_Setup'."
       return
    end if

    stp%ntypes_global = size(global_types(:))

    read(u) stp%description
    read(u) stp%atomtype
    read(u) stp%nenv
    allocate(stp%envtypes(stp%nenv))
    read(u) stp%envtypes
    read(u) stp%Rc_min
    read(u) stp%Rc_max
    read(u) stp%sftype
    read(u) stp%nsf
    read(u) stp%nsfparam
    allocate(stp%sf(stp%nsf),                    &
             stp%sfparam(stp%nsfparam, stp%nsf), &
             stp%sfenv(NENV_MAX,stp%nsf),        &
             stp%sfval_min(stp%nsf),             &
             stp%sfval_max(stp%nsf),             &
             stp%sfval_avg(stp%nsf),             &
             stp%sfval_cov(stp%nsf)              )
    read(u) stp%sf(:)
    read(u) stp%sfparam(:,:)
    read(u) stp%sfenv(:,:)
    read(u) stp%neval
    read(u) stp%sfval_min(:)
    read(u) stp%sfval_max(:)
    read(u) stp%sfval_avg(:)
    read(u) stp%sfval_cov(:)

    if (.not. present(unit)) then
       close(u)
    end if


    !shimamura for s-ANN
    !! connect local atom type IDs with global ones
!    write(*,*)"efile = in load_Setup",efile
!    write(*,*)"stp%ntypes_global=",stp%ntypes_global
    if(efile)then
        allocate(stp%gtype(stp%nenv), stp%ltype(stp%nenv))
        call stp_set_global_types_sANN(stp, stp%ntypes_global, global_types)
    else 
        allocate(stp%gtype(stp%nenv), stp%ltype(stp%ntypes_global))
        call stp_set_global_types(stp, stp%ntypes_global, global_types)
    end if

    stp%init = .true.

  end function load_Setup

  !--------------------------------------------------------------------!

  !shimamura for s-ANN
  function load_Setup_ASCII(global_types, efile, file, unit) result(stp)
  !function load_Setup_ASCII(global_types, file, unit) result(stp)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    integer,          optional,     intent(in) :: unit
    type(Setup)                                :: stp

    integer :: i, j
    integer :: u

    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'
    character(len=*), parameter :: ifrmt = '(4(1x,I17))'

    !shimamura for s-ANN
    logical,                        intent(in) :: efile


    !!shimamura
    !write(*,*)"In load_Setup_ASCII"

    if (present(unit)) then
       u = unit
!    else if (present(file)) then
!       u = io_unit()
!       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit number nor file name given"
       write(0,*) "       in `load_Setup'."
       return
    end if

    stp%ntypes_global = size(global_types(:))

    read(u,'(A)') stp%description
    read(u,'(A)') stp%atomtype
    read(u,*) stp%nenv
    allocate(stp%envtypes(stp%nenv))
    read(u,'(A)') (stp%envtypes(i), i=1,stp%nenv)
    read(u,*) stp%Rc_min
    read(u,*) stp%Rc_max
    read(u,'(A)') stp%sftype
    read(u,*) stp%nsf
    read(u,*) stp%nsfparam
    allocate(stp%sf(stp%nsf),                    &
             stp%sfparam(stp%nsfparam, stp%nsf), &
             stp%sfenv(NENV_MAX,stp%nsf),        &
             stp%sfval_min(stp%nsf),             &
             stp%sfval_max(stp%nsf),             &
             stp%sfval_avg(stp%nsf),             &
             stp%sfval_cov(stp%nsf)              )
    read(u,ifrmt) (stp%sf(i), i=1,stp%nsf)
    read(u,dfrmt) ((stp%sfparam(i,j), i=1,stp%nsfparam), j=1,stp%nsf)
!    read(u,ifrmt) ((stp%sfenv(i,j), i=1,stp%nenv), j=1,stp%nsf)
    read(u,ifrmt) ((stp%sfenv(i,j), i=1,NENV_MAX), j=1,stp%nsf)
    read(u,*) stp%neval
    read(u,dfrmt) (stp%sfval_min(i), i=1,stp%nsf)
    read(u,dfrmt) (stp%sfval_max(i), i=1,stp%nsf)
    read(u,dfrmt) (stp%sfval_avg(i), i=1,stp%nsf)
    read(u,dfrmt) (stp%sfval_cov(i), i=1,stp%nsf)

    if (.not. present(unit)) then
       close(u)
    end if


    !!shimamura
    !write(*,*)"after read"


    !shimamura for s-ANN
    !! connect local atom type IDs with global ones
!    write(*,*)"efile = in load_Setup_ASCII",efile
!    write(*,*)"stp%ntypes_global=",stp%ntypes_global
    if(efile)then
        allocate(stp%gtype(stp%nenv), stp%ltype(stp%nenv))
        call stp_set_global_types_sANN(stp, stp%ntypes_global, global_types)
    else
        allocate(stp%gtype(stp%nenv), stp%ltype(stp%ntypes_global))
        call stp_set_global_types(stp, stp%ntypes_global, global_types)
    end if


    stp%init = .true.

  end function load_Setup_ASCII


!shimamura for single ANN
  !--------------------------------------------------------------------!

  subroutine stp_set_global_types_sANN(stp, ntypes_global, global_types)

    implicit none

    type(Setup),                                intent(inout) :: stp
    integer,                                    intent(in)    :: ntypes_global
    character(len=*), dimension(ntypes_global), intent(in)    :: global_types

    integer :: i, j


    !!shimamura for s-ANN
    !write(*,*)"check1, stp_set_global_types_sANN"
    !write(*,*)stp%nenv


    ! atom type indices that have no corresponding entry are set to 0

    do i = 1, stp%nenv
       stp%ltype(i) = 0
       env : do j = 1, stp%nenv
          if (trim(stp%envtypes(j)) == trim(stp%envtypes(i))) then
             stp%ltype(i) = j
             exit env
          end if
       end do env
    end do

    ! reverse direction, because we do not know, if the sets of types
    ! are identical
    do i = 1, stp%nenv
       stp%gtype(i) = 0
       global : do j = 1, stp%nenv
          if (trim(stp%envtypes(j)) == trim(stp%envtypes(i))) then
             stp%gtype(i) = j
             exit global
          end if
       end do global
    end do


    !!shimamura for s-ANN
    !write(*,*)"check2, stp_set_global_types_sANN"
    !write(*,*)stp%ltype
    !write(*,*)size(stp%ltype)
    !write(*,*)stp%gtype
    !write(*,*)size(stp%gtype)
    !write(*,*)stp%envtypes


  end subroutine stp_set_global_types_sANN

  !--------------------------------------------------------------------!

  subroutine stp_set_global_types(stp, ntypes_global, global_types)

    implicit none

    type(Setup),                                intent(inout) :: stp
    integer,                                    intent(in)    :: ntypes_global
    character(len=*), dimension(ntypes_global), intent(in)    :: global_types

    integer :: i, j

    ! atom type indices that have no corresponding entry are set to 0

    do i = 1, ntypes_global
       stp%ltype(i) = 0
       env : do j = 1, stp%nenv
          if (trim(io_lower(global_types(i))) == trim(io_lower(stp%envtypes(j)))) then
             stp%ltype(i) = j
             exit env
          end if
       end do env
    end do

    ! reverse direction, because we do not know, if the sets of types
    ! are identical
    do i = 1, stp%nenv
       stp%gtype(i) = 0
       global : do j = 1, ntypes_global
          if (trim(io_lower(global_types(j))) == trim(io_lower(stp%envtypes(i)))) then
             stp%gtype(i) = j
             exit global
          end if
       end do global
    end do

  end subroutine stp_set_global_types

  !====================================================================!
  !                                                                    !
  !       structural fingerprint basis function set-up module          !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !                    initialization/finalization                     !
  !--------------------------------------------------------------------!

  subroutine stp_init(ntypes, stp, N_nb_max)

    implicit none

    integer,                        intent(in)  :: ntypes
    type(Setup), dimension(ntypes), intent(in)  :: stp
    integer,                        intent(in)  :: N_nb_max

    integer            :: itype
    character(len=100) :: sftype
    integer            :: nG_max

    if (isInit) then
       write(0,*) "Error: module already initialized in `stp_init'."
       stop
    end if

    sftype = trim(stp(1)%sftype)
    do itype = 1, ntypes
       if (trim(stp(itype)%sftype) /= trim(sftype)) then
          write(0,*) "Error: Mixing of basis functions of different " &
                         // "types not yet implemented."
          write(0,*) trim(sftype), trim(stp(itype)%sftype)
          stop
       end if
    end do

    select case(trim(io_lower(sftype)))

    !shimamura
    case('chebyshev')
       allocate(sfb(ntypes))
       do itype = 1, ntypes
          call setup_basis_chebyshev(stp(itype), sfb(itype))
       end do

    case('behler2011')
       nG_max = 0
       do itype = 1, ntypes
          nG_max = max(nG_max, stp(itype)%nsf)
       end do
       call sf_init(ntypes, nG_max)
       do itype = 1, ntypes
          call setup_symmfunc_Behler2011(itype, stp(itype))
       end do
    case default
       write(0,*) "Error: Unknown basis function type : ", trim(sftype)
       stop
    end select

    ! store max number of SFs in module
    nsf_max = stp_nsf_max(stp)

    ! max number of atoms in interaction range
    nnb_max = N_nb_max

    ! knomura
    if(N_nb_max > 2000) then
       !print*,'INFO : N_nb_max is too large. setting it to 2000.', nnb_max
       nnb_max = 2000
    endif

    ! allocate workspace for basis function evaluation:
    allocate(sfval(nsf_max), sfderiv_i(3,nsf_max), sfderiv_j(3,nsf_max,nnb_max),  &
&            sfstrs(6,nsf_max))
    sfval(:) = 0.0d0
    sfderiv_i(:,:) = 0.0d0
    sfderiv_j(:,:,:) = 0.0d0
    sfstrs(:,:) = 0.0d0

    isInit = .true.

  end subroutine stp_init

  !--------------------------------------------------------------------!
  !                      basis function evaluation                     !
  !--------------------------------------------------------------------!

  subroutine stp_eval(itype0, coo0, n, coo1, type1, stp, deriv, scaled, lstrs)

    implicit none

    integer,                                            intent(in)    :: itype0
    double precision, dimension(3),                     intent(in)    :: coo0
    integer,                                            intent(in)    :: n
    double precision, dimension(3,n),                   intent(in)    :: coo1
    integer,          dimension(n),                     intent(in)    :: type1
    type(Setup),                                        intent(inout) :: stp
    logical,                                  optional, intent(in)    :: deriv
    logical,                                  optional, intent(in)    :: scaled
    logical,                                  optional, intent(in)    :: lstrs

    integer               :: type0_loc
!    integer, dimension(n) :: type1_loc

    integer          :: i, nsf
    logical          :: do_deriv, do_scale, do_strs


    !---memory allocation if necessary
    if( max( n, 1 ) > nx ) then
        !-----if already allocated, deallocate arrays
        if( allocated(type1_loc) ) then
            deallocate( type1_loc )
        end if
        nx = max( n, 1 )
        !------allocate memory
        allocate( type1_loc(nx) )
    end if

    call stp_assert_moduleinit()
    call stp_assert_init(stp)

    if (present(deriv)) then
       do_deriv = deriv
    else
       do_deriv = .false.
    end if

    if (present(scaled)) then
       do_scale = scaled
    else
       do_scale = .false.
    end if

    if (present(lstrs)) then
       do_strs = lstrs
    else
       do_strs = .false.
    end if

    ! convert global atom type IDs to setup local IDs
    type0_loc = stp%ltype(itype0)
    do i = 1, n
       type1_loc(i) = stp%ltype(type1(i))
    end do

    select case(trim(io_lower(stp%sftype)))

    !shimamura
    case('chebyshev')
       nsf = stp%nsf
       if (do_deriv) then
           if (do_strs) then
               call sfb_eval(sfb(itype0), type0_loc, coo0, n, type1_loc, coo1, &
                                 nsf, sfval(1:nsf), sfderiv_i(1:3,1:nsf), &
                                 sfderiv_j(1:3,1:nsf,1:n), sfstrs(1:6,1:nsf))
           else
               call sfb_eval(sfb(itype0), type0_loc, coo0, n, type1_loc, coo1, &
                                 nsf, sfval(1:nsf), sfderiv_i(1:3,1:nsf), &
                                 sfderiv_j(1:3,1:nsf,1:n))
           end if
       else
          call sfb_eval(sfb(itype0), type0_loc, coo0, n, type1_loc, coo1, &
                        nsf, sfval(1:nsf))
       end if

    case('behler2011')
       nsf = stp%nsf
       if (do_deriv) then
           if (do_strs) then
               call sf_fingerprint(type0_loc, coo0, n, coo1, type1_loc, nsf, &
                              sfval(1:nsf), sfderiv_i(1:3,1:nsf), &
                              sfderiv_j(1:3,1:nsf,1:n), sfstrs(1:6,1:nsf))
           else
               call sf_fingerprint(type0_loc, coo0, n, coo1, type1_loc, nsf, &
                              sfval(1:nsf), sfderiv_i(1:3,1:nsf), &
                              sfderiv_j(1:3,1:nsf,1:n))
           end if
       else
          call sf_fingerprint(type0_loc, coo0, n, coo1, type1_loc, nsf, &
                              sfval(1:nsf))
       end if
    end select

    if (do_scale) then
       call stp_normalize(stp, n, deriv=do_deriv, lstrs=do_strs)
    else
       if (stp%neval == 0) then
          stp%sfval_min(1:nsf) = sfval(1:nsf)
          stp%sfval_max(1:nsf) = sfval(1:nsf)
          stp%sfval_avg(1:nsf) = sfval(1:nsf)
          stp%sfval_cov(1:nsf) = sfval(1:nsf)*sfval(1:nsf)
       else
          do i = 1, stp%nsf
             stp%sfval_min(i) = min(stp%sfval_min(i), sfval(i))
             stp%sfval_max(i) = max(stp%sfval_max(i), sfval(i))
             stp%sfval_avg(i) = (dble(stp%neval)*stp%sfval_avg(i) &
                                + sfval(i))/(dble(stp%neval+1))
             stp%sfval_cov(i) = (dble(stp%neval)*stp%sfval_cov(i) &
                                + sfval(i)*sfval(i))/(dble(stp%neval+1))
          end do
       end if
    end if

    stp%neval = stp%neval + 1

  end subroutine stp_eval

  !--------------------------------------------------------------------!
  !          normalization of basis function values to [-1,1]          !
  !--------------------------------------------------------------------!

  subroutine stp_normalize(stp, n, deriv, lstrs)

    implicit none

    type(Setup),       intent(inout) :: stp
    integer,           intent(in)    :: n
    logical, optional, intent(in)    :: deriv
    logical, optional, intent(in)    :: lstrs

    double precision :: scale, shift, s
    integer          :: isf
    logical          :: do_deriv, do_strs

    if (present(deriv)) then
       do_deriv = deriv
    else
       do_deriv = .false.
    end if

    if (present(lstrs)) then
       do_strs = lstrs
    else
       do_strs = .false.
    end if

    do isf = 1, stp%nsf
       shift = stp%sfval_avg(isf)
       ! scale covariance to 1
       ! s = sqrt(stp%sfval_cov(isf) + shift*shift - 2.0d0*shift*stp%sfval_avg(isf))
       s = sqrt(stp%sfval_cov(isf) - shift*shift)
       if (s <= 1.0d-10) then
!kn          write(0,*) "Warning: Invalid scaling encountered in ", &
!kn                               "'stp_normalize()'."
!kn          write(0,*) "         This means at least one fingerprint ", &
!kn                               "function for ", trim(adjustl(stp%atomtype)), &
!kn                               " is always equal to zero!"
!kn          write(0,*) "         Maybe an atomic species is not present ", &
!kn                               "in the reference set?"
!kn!          write(0,*) "         type       = ", trim(adjustl(stp%sftype)), &
!kn!                               " ", trim(io_adjustl(stp%sf(isf)))
          call io_adjustl_sub(stp%sf(isf),str)
!kn          write(0,*) "         type       = ", trim(adjustl(stp%sftype)), &
!kn                               " ", trim(str)
!kn          write(0,*) "         covariance = ", stp%sfval_cov(isf)
!kn          write(0,*) "         average    = ", stp%sfval_avg(isf)
!kn          write(0,*) "         min, max   = ", stp%sfval_min(isf), &
!kn                                               stp%sfval_max(isf)
          scale = 0.0d0
       else
          scale = 1.0d0/s
       end if
       sfval(isf) = scale*(sfval(isf)-shift)
       if (do_deriv) then
          sfderiv_i(1:3,isf)   = scale*sfderiv_i(1:3,isf)
!          sfderiv_j(1:3,isf,:) = scale*sfderiv_j(1:3,isf,:)
          sfderiv_j(1:3,isf,1:n) = scale*sfderiv_j(1:3,isf,1:n)
          if (do_strs) then
              sfstrs(1:6,isf)   = scale*sfstrs(1:6,isf)
          end if
       end if
    end do

  end subroutine stp_normalize

  !--------------------------------------------------------------------!
  !            maximum number of basis functions per setup             !
  !--------------------------------------------------------------------!

  function stp_nsf_max(stp) result(N_sf_max)

    implicit none

    type(Setup), dimension(:), optional, intent(in) :: stp
    integer                                         :: N_sf_max

    integer :: itype, nTypes

    if (present(stp)) then
       nTypes = size(stp(:))
       N_sf_max = 0
       do itype = 1, nTypes
          call stp_assert_init(stp(itype))
          N_sf_max = max(N_sf_max, stp(itype)%nsf)
       end do
    else if (isInit) then
       N_sf_max = nsf_max
    else
       write(0,*) "Error: module not initialized in `stp_nsf_max()'."
       stop
    end if

  end function stp_nsf_max

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine stp_assert_init(stp)
    implicit none
    type(Setup), intent(in) :: stp
    if (.not. stp%init) then
       write(*,*) 'Error: The basis function set-up is not initialized.'
       write(*,*)
       stop
    end if
  end subroutine stp_assert_init

  !--------------------------------------------------------------------!

  subroutine stp_assert_moduleinit()
    implicit none
    if (.not. isInit) then
       write(*,*) 'Error: Structural fingerprint setup module NOT initialized.'
       write(*,*)
       stop
    end if
  end subroutine stp_assert_moduleinit

  !--------------------------------------------------------------------!
  ! set up symmetry functions of the Behler2011 type as specified by   !
  ! the provided set-up `stp'                                          !
  !--------------------------------------------------------------------!

  subroutine setup_symmfunc_Behler2011(itype, stp)

    implicit none

    integer,     intent(in) :: itype
    type(Setup), intent(in) :: stp

    integer          :: type1, type2, type3
    double precision :: Rc, Rs, eta, kappa, lambda, zeta
    integer          :: isf, funct

    type1 = itype

    SFs : do isf = 1, stp%nsf

       funct = stp%sf(isf)

       select case(funct)
       case(1)
          type2  = stp%sfenv(1,isf)
          Rc     = stp%sfparam(1,isf)
          call sf_add_rad(funct, type1, type2, Rc=Rc)
       case(2)
          type2  = stp%sfenv(1,isf)
          Rc     = stp%sfparam(1,isf)
          Rs     = stp%sfparam(2,isf)
          eta    = stp%sfparam(3,isf)
          call sf_add_rad(funct, type1, type2, Rc=Rc, Rs=Rs, eta=eta)
       case(3)
          type2  = stp%sfenv(1,isf)
          Rc     = stp%sfparam(1,isf)
          kappa  = stp%sfparam(2,isf)
          call sf_add_rad(funct, type1, type2, Rc=Rc, kappa=kappa)
       case(4)
          type2  = stp%sfenv(1,isf)
          type3  = stp%sfenv(2,isf)
          Rc     = stp%sfparam(1,isf)
          lambda = stp%sfparam(2,isf)
          zeta   = stp%sfparam(3,isf)
          eta    = stp%sfparam(4,isf)
          call sf_add_ang(funct, type1, type2, type3, Rc=Rc, lambda=lambda, &
                      zeta=zeta, eta=eta)
       case(5)
          type2  = stp%sfenv(1,isf)
          type3  = stp%sfenv(2,isf)
          Rc     = stp%sfparam(1,isf)
          lambda = stp%sfparam(2,isf)
          zeta   = stp%sfparam(3,isf)
          eta    = stp%sfparam(4,isf)
          call sf_add_ang(funct, type1, type2, type3, Rc=Rc, lambda=lambda, &
                      zeta=zeta, eta=eta)
       case default
!          write(0,*) "Error: Symmetry function type not implemented : ", &
!                     trim(io_adjustl(stp%sf(isf)))
          call io_adjustl_sub(stp%sf(isf),str)
          write(0,*) "Error: Symmetry function type not implemented : ", &
                     trim(str)
          write(0,*) "in   : `setup_symmfunc_Behler2011'"
          stop
       end select

   end do SFs

  end subroutine setup_symmfunc_Behler2011

  !--------------------------------------------------------------------!
  !                        print info to stdout                        !
  !--------------------------------------------------------------------!

  subroutine print_info_Behler2011(stp, unit)

    implicit none

    type(Setup), intent(in) :: stp
    integer,     intent(in) :: unit

    integer            :: isf, iG
    double precision   :: Rc, Rs, eta, kappa, lambda, zeta
    character(len=128) :: frmt
    integer            :: type2, type3

    write(unit,*) 'Basis function type Behler2011'
    write(unit,*)
    write(unit,*) '[see also: J. Behler, J. Chem. Phys. 134 (2011) 074106]'
    write(unit,*)
    write(unit,*)

    write(unit,'(6x,"G",2x,"parameters")')
    write(unit,*)

    do isf = 1, stp%nsf
       iG = stp%sf(isf)
       select case(iG)
       case(1)
          type2 = stp%sfenv(1,isf)
          Rc = stp%sfparam(1,isf)
          write(unit,'(1x,I3,2x,I1,2x,"type2=",A2,"  Rc = ",F7.3)') &
               isf, iG, stp%envtypes(type2), Rc
       case(2)
          type2 = stp%sfenv(1,isf)
          Rc  = stp%sfparam(1,isf)
          Rs  = stp%sfparam(2,isf)
          eta = stp%sfparam(3,isf)
          write(unit,'(1x,I3,2x,I1,2x,"type2=",A2,"  Rc = ",F7.3,"  Rs = ",F7.3,"  eta = ",F7.3)') &
               isf, iG, stp%envtypes(type2), Rc, Rs, eta
       case(3)
          type2 = stp%sfenv(1,isf)
          Rc    = stp%sfparam(1,isf)
          kappa = stp%sfparam(2,isf)
          write(unit,'(1x,I3,2x,I1,2x,"type2=",A2,"  Rc = ",F7.3,"  kappa = ",F7.3)') &
               isf, iG, stp%envtypes(type2), Rc, kappa
       case(4,5)
          type2  = stp%sfenv(1,isf)
          type3  = stp%sfenv(2,isf)
          Rc     = stp%sfparam(1,isf)
          lambda = stp%sfparam(2,isf)
          zeta   = stp%sfparam(3,isf)
          eta    = stp%sfparam(4,isf)
          frmt = '(1x,I3,2x,I1,2x,"type2=",A2,"  type3=",A2,' &
               // '"  Rc = ",F7.3,"  lambda = ",F7.3,' &
               // '"  zeta = ",F7.3,"  eta = ",F7.3)'
          write(unit,frmt) isf, iG, stp%envtypes(type2), &
               stp%envtypes(type3), Rc, lambda, zeta, eta
       end select
    end do
    write(unit,*)

  end subroutine print_info_Behler2011

  !shimamura
  !--------------------------------------------------------------------!

  subroutine print_info_chebyshev(stp, unit)

    implicit none

    type(Setup), intent(in) :: stp
    integer,     intent(in) :: unit

    double precision :: r_Rc, a_Rc
    integer          :: r_N, a_N

    r_Rc = stp%sfparam(1,1)
    r_N = nint(stp%sfparam(2,1))
    a_Rc = stp%sfparam(3,1)
    a_N = nint(stp%sfparam(4,1))

    write(unit,*) 'Basis function type Chebyshev'
    write(unit,*) '[N. Artrith and A. Urban (2016)]'
    write(unit,*)
    write(unit,*) 'Radial Rc     : ' // trim(io_adjustl(r_Rc))
    write(unit,*) 'Angular Rc    : ' // trim(io_adjustl(a_Rc))
    write(unit,*) 'Radial order  : ' // trim(io_adjustl(r_N))
    write(unit,*) 'Angular order : ' // trim(io_adjustl(a_N))
    write(unit,*)

  end subroutine print_info_chebyshev

  !shimamura
  !--------------------------------------------------------------------!

  subroutine setup_basis_chebyshev(stp, sfb)

    implicit none

    type(Setup),            intent(in)  :: stp
    type(FingerprintBasis), intent(out) :: sfb

    double precision :: r_Rc, a_Rc
    integer          :: r_N, a_N

    r_Rc = stp%sfparam(1,1)
    r_N = nint(stp%sfparam(2,1))
    a_Rc = stp%sfparam(3,1)
    a_N = nint(stp%sfparam(4,1))

    sfb = new_SFBasis(stp%nenv, stp%envtypes, r_N, a_N, r_Rc, a_Rc)

  end subroutine setup_basis_chebyshev

end module sfsetup




module trainset

  use aeio,    only: aeio_header,           &
                     TYPELEN, PATHLEN

  use io,      only: io_adjustl, io_adjustl_sub,  &
&                    nfile, myid, nodes

  implicit none
  private
  save

  public  :: load_TrnSet_info,        &
             load_TrnSet_info_ASCII,  &
             ts_print_info!,           &

  type, public :: TrnSet

     !-----------------------------------------------------------------!
     ! init        .true., if the training set has been initialized    !
     ! file        name of the corresponding training set file         !
     ! unit        unit number of that file                            !
     ! mode        current access mode; 'read', 'write', 'info'        !
     ! normalized  .true., if the input and output values have been    !
     !             normalized to the interval [-1,1] ('read' mode only)!
     !                                                                 !
     ! if (normalized == .true.)                                       !
     ! scale       energy scaling factor used for the normalization    !
     ! shift       atomic energy shift used for energy normalization   !
     !                                                                 !
     ! nTypes      number of atomic species in the training set        !
     ! nAtomsTot   total number of atoms in the training set           !
     ! typeName(i) name of i-th atomic species                         !
     ! nStrucs     total number of structures in the training set      !
     !             --> when open in 'write' mode, not necessarily all  !
     !                 files have yet been parsed                      !
     ! iStruc      current file record position (0=before first file)  !
     !-----------------------------------------------------------------!

     logical                                           :: init = .false.

     character(len=PATHLEN)                            :: file
     integer                                           :: unit
     character(len=5)                                  :: mode

     logical                                           :: normalized
     double precision                                  :: scale, shift

     integer                                           :: nTypes
     integer                                           :: nAtomsTot
     character(len=TYPELEN), dimension(:), allocatable :: typeName
     double precision,       dimension(:), allocatable :: E_atom
     integer                                           :: nStrucs
     integer                                           :: iStruc

     double precision                                  :: E_min, E_max, E_av

  end type TrnSet

contains

  !--------------------------------------------------------------------!

  function load_TrnSet_info(file, unit) result(ts)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(TrnSet)                           :: ts

    integer :: u

    if (present(unit)) then
       u = unit
!    else if (present(file)) then
!       u = io_unit()
!       open(u, file=trim(file), status='old', &
!            form='unformatted', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in ", &
                  "`load_TrnSet_info()'."
       stop
    end if

    read(u) ts%file
    read(u) ts%normalized
    read(u) ts%scale
    read(u) ts%shift
    read(u) ts%nTypes
    allocate(ts%typeName(ts%nTypes), ts%E_atom(ts%nTypes))
    read(u) ts%typeName(1:ts%nTypes)
    read(u) ts%E_atom(1:ts%nTypes)
    read(u) ts%nAtomsTot
    read(u) ts%nStrucs
    read(u) ts%E_min, ts%E_max, ts%E_av

    ts%iStruc = ts%nStrucs
    ts%init = .true.
    ts%mode = 'info'

    if (.not. present(unit)) close(u)

  end function load_TrnSet_info

  !--------------------------------------------------------------------!

  function load_TrnSet_info_ASCII(file, unit) result(ts)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(TrnSet)                           :: ts

    integer :: i
    integer :: u

    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'

    if (present(unit)) then
       u = unit
!    else if (present(file)) then
!       u = io_unit()
!       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `load_TrnSet_info()'."
       stop
    end if

    read(u,'(A)') ts%file
    read(u,*) ts%normalized
    read(u,*) ts%scale
    read(u,*) ts%shift
    read(u,*) ts%nTypes
    allocate(ts%typeName(ts%nTypes), ts%E_atom(ts%nTypes))
    read(u,'(A)') (ts%typeName(i), i=1,ts%nTypes)
    read(u,dfrmt) (ts%E_atom(i), i=1,ts%nTypes)
    read(u,*) ts%nAtomsTot
    read(u,*) ts%nStrucs
    read(u,*) ts%E_min, ts%E_max, ts%E_av

    ts%iStruc = ts%nStrucs
    ts%init = .true.
    ts%mode = 'info'

    if (.not. present(unit)) close(u)

  end function load_TrnSet_info_ASCII

  !--------------------------------------------------------------------!
  !                    info about the training set                     !
  !--------------------------------------------------------------------!

  subroutine ts_print_info(ts)

    implicit none

    type(TrnSet), intent(in) :: ts

    integer :: itype, nf
    character(len=50) :: str

    call ts_assert_init(ts)

    call aeio_header("Training set info.")

 do nf = 1, 2
 
   write(nfile(nf),*)

    write(nfile(nf),*) 'Training set file                   : ', trim(adjustl(ts%file))
!    write(*,*) 'Number of structures in the data set: ', trim(io_adjustl(ts%nStrucs))
    call io_adjustl_sub(ts%nStrucs,str)
    write(nfile(nf),*) 'Number of structures in the data set: ', trim(str)
    if (ts%iStruc /= ts%nStrucs) then
       if (trim(ts%mode) == 'write') then
!          write(*,*) '  Structures included so far        : ', trim(io_adjustl(ts%iStruc))
          call io_adjustl_sub(ts%iStruc,str)
          write(nfile(nf),*) '  Structures included so far        : ', trim(str)
       else
!          write(*,*) '  Structures read so far            : ', trim(io_adjustl(ts%iStruc))
          call io_adjustl_sub(ts%iStruc,str)
          write(nfile(nf),*) '  Structures read so far            : ', trim(str)
       end if
    end if
    write(nfile(nf),*)

!    write(*,*) 'Atomic species in training set      : ', trim(io_adjustl(ts%nTypes))
    call io_adjustl_sub(ts%nTypes,str)
    write(nfile(nf),*) 'Atomic species in training set      : ', trim(str)
    write(nfile(nf),'(1x,"  Species :")', advance='no')
    do itype = 1, ts%nTypes
       write(nfile(nf),'(1x,A)', advance='no') trim(ts%typeName(itype))
    end do
    write(nfile(nf),*)
    write(nfile(nf),*)

    if (ts%normalized .or. (ts%iStruc == ts%nStrucs)) then
!       write(*,*) 'Average energy (eV/atom) : ', trim(io_adjustl(ts%E_av,6))
!       write(*,*) 'Minimum energy (eV/atom) : ', trim(io_adjustl(ts%E_min,6))
!       write(*,*) 'Maximum energy (eV/atom) : ', trim(io_adjustl(ts%E_max,6))
       call io_adjustl_sub(ts%E_av,str,6)
       write(nfile(nf),*) 'Average energy (eV/atom) : ', trim(str)
       call io_adjustl_sub(ts%E_min,str,6)
       write(nfile(nf),*) 'Minimum energy (eV/atom) : ', trim(str)
       call io_adjustl_sub(ts%E_max,str,6)
       write(nfile(nf),*) 'Maximum energy (eV/atom) : ', trim(str)
       write(nfile(nf),*)
    end if

    if (ts%normalized) then
       write(nfile(nf),*) 'The input and output values have been normalized to [-1.0, 1.0].'
       write(nfile(nf),*) 'Structures outside of this interval will not be used for training.'
!       write(*,*) '  Energy scaling factor: ', trim(io_adjustl(ts%scale,6))
!       write(*,*) '  Atomic energy shift  : ', trim(io_adjustl(ts%shift,6))
       call io_adjustl_sub(ts%scale,str,6)
       write(nfile(nf),*) '  Energy scaling factor: ', trim(str)
       call io_adjustl_sub(ts%shift,str,6)
       write(nfile(nf),*) '  Atomic energy shift  : ', trim(str)
    else
       write(nfile(nf),*) 'The input and output values have not yet been normalized.'
    end if
    write(nfile(nf),*)

 end do

  end subroutine ts_print_info

  !--------------------------------------------------------------------!
  !                            state checks                            !
  !--------------------------------------------------------------------!

  subroutine ts_assert_init(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (.not. ts%init) then
       write(0,*) "Error: training set not initialized."
       stop
    end if
  end subroutine ts_assert_init

end module trainset




module potential

  use io,          only: nfile, myid, nodes

  use aeio,        only: PATHLEN, TYPELEN, &
                         aeio_assert_file_exists

  use geometry,    only: geo_itype_of_name

  use feedforward, only: Network,          &
                         load_Network,     &
                         load_Network_ASCII, &
                         ff_print_info

  use sfsetup,     only: Setup,            &
                         load_Setup,       &
                         load_Setup_ASCII, &
                         stp_print_info

  use trainset,    only: TrnSet,           &
                         load_TrnSet_info, &
                         load_TrnSet_info_ASCII,  &
                         ts_print_info

  implicit none
  private
  save

  public :: load_NNPot,     &
            pot_get_range,  &
            pot_print_info!, &

  !--------------------------------------------------------------------!
  !                      Neural Network potential                      !
  !--------------------------------------------------------------------!

  type, public :: NNPot

     !-----------------------------------------------------------------!
     ! init        .true. if potential has been initialized            !
     ! file        path to the potential file, if known                !
     ! unit        unit number, if not directly loaded from file 'unit'!
     !             will be set to -1, if the potential was loaded      !
     !             directly from a file                                !
     ! typeName    species of the central atom                         !
     !                                                                 !
     ! E_scale     energy scaling factor                               !
     ! E_shift     shift of the atomic energy                          !
     ! E_atom      atomic reference energy of the central atom         !
     !             This is the shifted atomic energy!  If you need the !
     !             reference atmic energy, use E_atom - E_shift .      !
     !                                                                 !
     ! stp         structural fingerprint basis setup (type: Setup)    !
     ! net         trained neural network (type: Network)              !
     ! ts          training set info (type: TrnSet)                    !
     !-----------------------------------------------------------------!

     logical                :: init = .false.

     character(len=PATHLEN) :: file
     integer                :: unit

     character(len=TYPELEN) :: typeName
     double precision       :: E_scale
     double precision       :: E_shift
     double precision       :: E_atom

     type(Setup)            :: stp
     type(Network)          :: net
     type(TrnSet)           :: ts

  end type NNPot

  !--------------------------------------------------------------------!

  logical :: isInit = .false.

contains


  !shimamura for s-ANN
  function load_NNPot(global_types, file, lascii, efile, unit) result(pot)
  !function load_NNPot(global_types, file, lascii, unit) result(pot)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    logical,          optional,     intent(in) :: lascii
    integer,          optional,     intent(in) :: unit
    type(NNPot)                                :: pot

    integer :: itype
    integer :: u
    logical :: lasc

    !shimamura for s-ANN
    logical,                        intent(in) :: efile


    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Error: neither file nor file unit specified in `load_NNPot()'."
       stop
    end if

    if (present(lascii)) then
        lasc = lascii
    else
        lasc = .false.
    end if

    if (present(unit)) then
       u = unit
       pot%unit = u
       pot%file = ''
    else
!       u = io_unit()
       call allocate_unit_number( u )
       call aeio_assert_file_exists(file)
       if( .not.lasc ) then
           open(u, file=trim(adjustl(file)), status='old', action='read', &
                form='unformatted')
       else
           open(u, file=trim(adjustl(file)), status='old', action='read', &
                form='formatted')
       end if
       pot%file = trim(adjustl(file))
       pot%unit = -1
    end if


    !shimamura for s-ANN
    if( .not.lasc ) then
        pot%net = load_Network(unit=u)
        pot%stp = load_Setup(global_types, efile, unit=u)
        pot%ts  = load_TrnSet_info(unit=u)
    else
        pot%net = load_Network_ASCII(unit=u)
        pot%stp = load_Setup_ASCII(global_types, efile, unit=u)
        pot%ts  = load_TrnSet_info_ASCII(unit=u)
    end if
    !if( .not.lasc ) then
    !    pot%net = load_Network(unit=u)
    !    pot%stp = load_Setup(global_types, unit=u)
    !    pot%ts  = load_TrnSet_info(unit=u)
    !else
    !    pot%net = load_Network_ASCII(unit=u)
    !    pot%stp = load_Setup_ASCII(global_types, unit=u)
    !    pot%ts  = load_TrnSet_info_ASCII(unit=u)
    !end if


    pot%E_scale     = 1.0d0/pot%ts%scale

    ! central atom
    pot%typeName = pot%stp%atomtype
    pot%E_shift  = pot%ts%shift
    itype = geo_itype_of_name(pot%typeName, pot%ts%typeName)
    pot%E_atom   = pot%ts%E_atom(itype)


!shimamura for s-ANN in load_NNPot
!write(*,*)"pot%typeName=",pot%typeName
!write(*,*)"pot%E_shift =",pot%E_shift
!write(*,*)"itype       =",itype
!write(*,*)"pot%E_atom  =",pot%E_atom
!write(*,*)"global_types=",global_types


    if (.not. present(unit)) then
         close(u)
         call deallocate_unit_number( u )
    end if

    pot%init = .true.

  end function load_NNPot

  !--------------------------------------------------------------------!
  !                   print info about NN potential                    !
  !--------------------------------------------------------------------!

  subroutine pot_print_info(pot)

    implicit  none

    type(NNPot), intent(in) :: pot

    call pot_assert_init(pot)

    write(nfile(1),*) 'Atomic species : ', trim(pot%typeName)
    write(nfile(2),*) 'Atomic species : ', trim(pot%typeName)
    write(nfile(1),*) 'File name      : ', trim(pot%file)
    write(nfile(2),*) 'File name      : ', trim(pot%file)
    write(nfile(1),*)
    write(nfile(2),*)

    call ts_print_info(pot%ts)
    call ff_print_info(pot%net)
    call stp_print_info(pot%stp)

  end subroutine pot_print_info

  !--------------------------------------------------------------------!
  !      determine interaction range for over several potentials       !
  !                (see also stp_get_range in sfsetup)                 !
  !--------------------------------------------------------------------!

  subroutine pot_get_range(nTypes, pot, Rc_min, Rc_max)

    implicit none

    integer,                        intent(in)  :: nTypes
    type(NNPot), dimension(nTypes), intent(in)  :: pot
    double precision,               intent(out) :: Rc_min
    double precision,               intent(out) :: Rc_max

    integer :: itype

    call pot_assert_init(pot(1))

    Rc_min = pot(1)%stp%Rc_min
    Rc_max = pot(1)%stp%Rc_max

    do itype = 1, nTypes
       call pot_assert_init(pot(itype))
       Rc_min = min(Rc_min, pot(itype)%stp%Rc_min)
       Rc_max = max(Rc_max, pot(itype)%stp%Rc_max)
    end do

  end subroutine pot_get_range

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine pot_assert_init(pot)
    implicit none
    type(NNPot), intent(in) :: pot
    if (.not. pot%init) then
       write(*,*) 'Error: NN potential is not initialized.'
       write(*,*)
       stop
    end if
  end subroutine pot_assert_init

end module potential




module aenet

  use aeio,        only: TYPELEN, PATHLEN

  use feedforward, only: ff_eval, ff_deriv

  use lclist,      only: lcl_nmax_nbdist

  use potential,   only: NNPot,               &
                         load_NNPot,          &
                         pot_print_info,      &
                         pot_get_range

  use sfsetup,     only: stp_init,            &
                         stp_nsf_max,         &
                         stp_eval,            &
                         sfval, sfderiv_i, sfderiv_j, sfstrs

  implicit none
  private
  save
  integer, parameter :: c_double = 8, c_int = 4, c_char = 1, c_bool = 1

  public :: aenet_init,                     &
            aenet_load_potential,           &
            aenet_atomic_energy_and_forces, &
            aenet_atomic_energy_and_forces_and_stresses, &
            aenet_free_atom_energy,         &
            aenet_print_info

  !---------------------------- constants -----------------------------!

  ! return status
!  integer(kind=c_int), bind(C, name='AENET_OK'),         public :: AENET_OK = 0_c_int
!  integer(kind=c_int), bind(C, name='AENET_ERR_INIT'),   public :: AENET_ERR_INIT = 1_c_int
!  integer(kind=c_int), bind(C, name='AENET_ERR_MALLOC'), public :: AENET_ERR_MALLOC = 2_c_int
!  integer(kind=c_int), bind(C, name='AENET_ERR_IO'),     public :: AENET_ERR_IO = 3_c_int
!  integer(kind=c_int), bind(C, name='AENET_ERR_TYPE'),   public :: AENET_ERR_TYPE = 4_c_int
  integer(kind=c_int),         public :: AENET_OK = 0_c_int
  integer(kind=c_int),   public :: AENET_ERR_INIT = 1_c_int
  integer(kind=c_int), public :: AENET_ERR_MALLOC = 2_c_int
  integer(kind=c_int),     public :: AENET_ERR_IO = 3_c_int
  integer(kind=c_int),   public :: AENET_ERR_TYPE = 4_c_int

!  integer(kind=c_int), bind(C, name='AENET_TYPELEN'),    public :: AENET_TYPELEN = TYPELEN
!  integer(kind=c_int), bind(C, name='AENET_PATHLEN'),    public :: AENET_PATHLEN = PATHLEN
  integer(kind=c_int),    public :: AENET_TYPELEN = TYPELEN
  integer(kind=c_int),    public :: AENET_PATHLEN = PATHLEN

  !---------------------------- variables -----------------------------!

!  integer(kind=c_int), bind(C, name='aenet_nsf_max'), public :: aenet_nsf_max
!  integer(kind=c_int), bind(C, name='aenet_nnb_max'), public :: aenet_nnb_max
!  real(kind=c_double), bind(C, name='aenet_Rc_min'),  public :: aenet_Rc_min
!  real(kind=c_double), bind(C, name='aenet_Rc_max'),  public :: aenet_Rc_max
  integer(kind=c_int), public :: aenet_nsf_max
  integer(kind=c_int), public :: aenet_nnb_max
  real(kind=c_double), public :: aenet_Rc_min
  real(kind=c_double), public :: aenet_Rc_max

  !----------------------------- private ------------------------------!

  logical, private :: aenet_is_init = .false.
  logical, private :: aenet_is_loaded = .false.

  integer,                                           private :: aenet_ntypes
  character(len=TYPELEN), dimension(:), allocatable, private :: aenet_atom_types
  type(NNPot),            dimension(:), allocatable, private :: aenet_pot
  double precision,       dimension(:), allocatable, private :: aenet_dE_dG

  integer,                                       public :: aenet_nnb_maxx = 0
  integer,          dimension(:),   allocatable, public :: nblist
  double precision, dimension(:,:), allocatable, public :: nbcoo
  double precision, dimension(:),   allocatable, public :: nbdist
  integer,          dimension(:),   allocatable, public :: nbtype

contains

  !--------------------------------------------------------------------!
  !                  Initialization and Finalization                   !
  !--------------------------------------------------------------------!

  subroutine aenet_init(atom_types, stat)

    implicit none

    character(len=*), dimension(:), intent(in)  :: atom_types
    integer,                        intent(out) :: stat

    integer :: ok

    stat = AENET_OK
    if (aenet_is_init) then
       stat = AENET_ERR_INIT
       return
    end if

    aenet_ntypes = size(atom_types)
    allocate(aenet_pot(aenet_ntypes),        &
             aenet_atom_types(aenet_ntypes), &
             stat=ok)
    if (ok /= 0) then
       aenet_ntypes = 0
       stat = AENET_ERR_MALLOC
       return
    end if
    aenet_atom_types = atom_types
    aenet_is_init = .true.

  end subroutine aenet_init

  !--------------------------------------------------------------------!
  !                               output                               !
  !--------------------------------------------------------------------!

  subroutine aenet_print_info() !bind(C)

    implicit none

    integer :: ipot

    if (.not. aenet_is_init) then
       write(*,*) "Nothing to print. AenetLib is not initialized."
    else
       do ipot = 1, aenet_ntypes
          if (aenet_pot(ipot)%init) then
             call pot_print_info(aenet_pot(ipot))
          end if
       end do
    end if

  end subroutine aenet_print_info

  !--------------------------------------------------------------------!
  !                        Load ANN potentials                         !
  !--------------------------------------------------------------------!

  !shimamura for s-ANN
  subroutine aenet_load_potential(type_id, filename, lascii, efile, stat)
  !subroutine aenet_load_potential(type_id, filename, lascii, stat)

    implicit none

    integer,          intent(in)  :: type_id
    character(len=*), intent(in)  :: filename
    logical,          intent(in)  :: lascii
    integer,          intent(out) :: stat

    integer :: ok
    logical :: fexists

    !shimamura for s-ANN
    logical,          intent(in)  :: efile


    stat = AENET_OK
    if (.not. aenet_is_init) then
       stat = AENET_ERR_INIT
       return
    end if

    if ((type_id <= 0) .or. (type_id > aenet_ntypes)) then
       stat = AENET_ERR_TYPE
       return
    end if

    inquire(file=trim(filename), exist=fexists)
    if (.not. fexists) then
       stat = AENET_ERR_IO
       return
    end if


    !nomura print*,'foo', type_id, filename 
    !shimamura for s-ANN
    aenet_pot(type_id) = load_NNPot(aenet_atom_types, filename, lascii, efile)
    !aenet_pot(type_id) = load_NNPot(aenet_atom_types, filename, lascii)


    ! when all potentials are loaded, determine array sizes
    if (aenet_all_loaded()) then
       call pot_get_range(aenet_ntypes, aenet_pot, aenet_Rc_min, aenet_Rc_max)
       if (stat /= 0) return
       aenet_nnb_max = lcl_nmax_nbdist(aenet_Rc_min, aenet_Rc_max)
       call stp_init(aenet_ntypes, aenet_pot(1:aenet_ntypes)%stp, aenet_nnb_max)
       aenet_nsf_max = stp_nsf_max()
       allocate(aenet_dE_dG(aenet_nsf_max), stat=ok)
       if (ok /= 0) then
          stat = AENET_ERR_MALLOC
          return
       end if
       aenet_is_loaded = .true.
    end if

  end subroutine aenet_load_potential

  function aenet_all_loaded() result(all_loaded) !bind(C)

    implicit none

    logical(kind=c_bool) :: all_loaded
    integer :: ipot

    all_loaded = .true.
    do ipot = 1, aenet_ntypes
       if (.not. aenet_pot(ipot)%init) then
          all_loaded = .false.
          return
       end if
    end do

  end function aenet_all_loaded

  !--------------------------------------------------------------------!
  !                            information                             !
  !--------------------------------------------------------------------!

  function aenet_free_atom_energy(type_id) result(E_atom) !bind(C)

    implicit none

    integer(kind=c_int), intent(in) :: type_id
    real(kind=c_double)             :: E_atom

    E_atom = aenet_pot(type_id)%E_atom

  end function aenet_free_atom_energy

  !--------------------------------------------------------------------!
  !                             Evaluation                             !
  !                                                                    !
  ! Attention: all routines require synchronized atom type IDs, i.e.,  !
  !            the IDs passed to the evaluation routines must be       !
  !            compatible with the ANN potential type IDs.             !
  !                                                                    !
  ! Notes:     * Coordinates are Cartesian.                            !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine aenet_atomic_energy_and_forces( &
       coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, natoms, &
       E_i, F, f3r, nmd_buffer, stat) !bind(C)

    implicit none

    real(kind=c_double), dimension(3),        intent(in)    :: coo_i
    integer(kind=c_int), value,               intent(in)    :: type_i
    integer(kind=c_int), value,               intent(in)    :: index_i
    integer(kind=c_int), value,               intent(in)    :: n_j
    real(kind=c_double), dimension(3,n_j),    intent(in)    :: coo_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: type_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: index_j
    integer(kind=c_int), value,               intent(in)    :: natoms
    real(kind=c_double),                      intent(out)   :: E_i
    real(kind=c_double), dimension(3,natoms), intent(inout) :: F
    integer(kind=c_int), value,               intent(in)    :: nmd_buffer
    real(kind=c_double), dimension(3,nmd_buffer), intent(inout) :: f3r
    integer(kind=c_int),                      intent(out)   :: stat

    double precision, dimension(1)              :: E_i_arr
    integer                                     :: nsf, j, kz

    stat = AENET_OK
    if (.not. (aenet_is_init .and. aenet_is_loaded)) then
       stat = aenet_ERR_INIT
       return
    end if

    nsf = aenet_pot(type_i)%stp%nsf
    call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
                  aenet_pot(type_i)%stp, scaled=.true., deriv=.true.)

    call ff_eval(aenet_pot(type_i)%net, nsf, sfval, 1, E_i_arr)
    call ff_deriv(aenet_pot(type_i)%net, nsf, 1, aenet_dE_dG(1:nsf))

    E_i = aenet_pot(type_i)%E_scale*E_i_arr(1) + aenet_pot(type_i)%E_shift
    E_i = E_i + aenet_pot(type_i)%E_atom

    F(1:3, index_i) = F(1:3, index_i) - aenet_pot(type_i)%E_scale &
                    * matmul(sfderiv_i(1:3,1:nsf), aenet_dE_dG(1:nsf))

    do j = 1, n_j
       if( index_j(j) <= natoms ) then
           F(1:3, index_j(j)) = F(1:3, index_j(j)) - aenet_pot(type_i)%E_scale &
                          * matmul(sfderiv_j(1:3,1:nsf,j), aenet_dE_dG(1:nsf))
       else
           !kz = index_j(j) - natoms
           kz = index_j(j)  !knomura
           f3r(1:3,kz) = f3r(1:3,kz) - aenet_pot(type_i)%E_scale &
                          * matmul(sfderiv_j(1:3,1:nsf,j), aenet_dE_dG(1:nsf))
       end if
    end do

  end subroutine aenet_atomic_energy_and_forces

  subroutine aenet_atomic_energy_and_forces_and_stresses( &
       coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, natoms, &
       E_i, F, f3r, nmd_buffer,  &
&      pit1lr, pit1sr, lcatomic, wstrs, stat)

    implicit none

    real(kind=c_double), dimension(3),        intent(in)    :: coo_i
    integer(kind=c_int), value,               intent(in)    :: type_i
    integer(kind=c_int), value,               intent(in)    :: index_i
    integer(kind=c_int), value,               intent(in)    :: n_j
    real(kind=c_double), dimension(3,n_j),    intent(in)    :: coo_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: type_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: index_j
    integer(kind=c_int), value,               intent(in)    :: natoms
    real(kind=c_double),                      intent(out)   :: E_i
    real(kind=c_double), dimension(3,natoms), intent(inout) :: F
    integer(kind=c_int), value,               intent(in)    :: nmd_buffer
    real(kind=c_double), dimension(3,nmd_buffer), intent(inout) :: f3r
    integer(kind=c_int),                      intent(out)   :: stat
real(8), intent(inout) :: pit1lr(6), pit1sr(6)
logical, intent(in)    :: lcatomic
real(8), intent(inout) :: wstrs(6)

    double precision, dimension(1)              :: E_i_arr
    integer                                     :: nsf, j, kz
    double precision :: strs(6)

    stat = AENET_OK
    if (.not. (aenet_is_init .and. aenet_is_loaded)) then
       stat = aenet_ERR_INIT
       return
    end if

    nsf = aenet_pot(type_i)%stp%nsf
    call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
                  aenet_pot(type_i)%stp, scaled=.true., deriv=.true., lstrs=.true.)

    call ff_eval(aenet_pot(type_i)%net, nsf, sfval, 1, E_i_arr)
    call ff_deriv(aenet_pot(type_i)%net, nsf, 1, aenet_dE_dG(1:nsf))

    E_i = aenet_pot(type_i)%E_scale*E_i_arr(1) + aenet_pot(type_i)%E_shift
    E_i = E_i + aenet_pot(type_i)%E_atom


!!shimamura
!write(*,*)aenet_pot(type_i)%E_scale,aenet_pot(type_i)%E_shift,aenet_pot(type_i)%E_atom


    F(1:3, index_i) = F(1:3, index_i) - aenet_pot(type_i)%E_scale &
                    * matmul(sfderiv_i(1:3,1:nsf), aenet_dE_dG(1:nsf))

    do j = 1, n_j
       if( index_j(j) <= natoms ) then
           F(1:3, index_j(j)) = F(1:3, index_j(j)) - aenet_pot(type_i)%E_scale &
                          * matmul(sfderiv_j(1:3,1:nsf,j), aenet_dE_dG(1:nsf))
       else
           kz = index_j(j) - natoms
           f3r(1:3,kz) = f3r(1:3,kz) - aenet_pot(type_i)%E_scale &
                          * matmul(sfderiv_j(1:3,1:nsf,j), aenet_dE_dG(1:nsf))
       end if
    end do

    strs(1:6) = - aenet_pot(type_i)%E_scale &
              * matmul(sfstrs(1:6,1:nsf), aenet_dE_dG(1:nsf))
    pit1lr(1:6) = pit1lr(1:6) + strs(1:6)
    if( lcatomic ) wstrs(1:6) = wstrs(1:6) + strs(1:6)

  end subroutine aenet_atomic_energy_and_forces_and_stresses

end module aenet



module md_potential_ANN
!-----------------------------------------------------------------------
! type declaration of parameters for Artificial Neural Network (ANN) potential
!-----------------------------------------------------------------------
implicit none

integer :: nc   = 1                           ! the number of components
character(80), allocatable :: ANNfiles(:)  ! network file name
character(2),  allocatable :: atom_types(:)
character(10) :: ANNlunit = '', ANNeunit = ''
real(8) :: unitlength = 1.d0
real(8) :: unitenergy = 1.d0
logical, allocatable :: lANNascii(:)

!!-----two-body potential parameters
!real*8  :: rc_2body                                  !  cutoff length in [a.u.]

!-----arrays for soft core
logical :: lANNsccal = .false.
logical, allocatable :: lANNsc(:)
real*8,  allocatable, dimension(:) :: ANNscepsilon, ANNscsigma, ANNscpower

!-----arrays for two-body potential
integer, parameter :: nbin = 4096 * 2                ! size of potential table
real*8,  allocatable, dimension(:,:) :: rcij, rcij2  ! cutoff length for i-j pair
real*8,  allocatable, dimension(:,:,:,:) :: v        ! two-body potential table
real*8  :: dr2i, rij2max

real*8,  allocatable, dimension(:,:) :: f3r          !

!-----cutoff length for pair list = rc + rc_extend
real*8 :: rc            != rc_2body                  !  cutoff length in [a.u.]
real*8,  parameter :: rc_extend = 4.d0               ! <----- should be optimized !!!


logical :: ladopted = .false.


!shimamura for s-ANN
logical :: efile = .false.

!-----for Coulomb potential
logical :: lANNcoulomb = .false.
real*8,  allocatable, dimension(:) :: zz

save

end module




subroutine set_lANNcoulomb( nfile, myid, nodes, lANNcoulomb_ )
!-----------------------------------------------------------------------
!     set lANNcoulomb
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
logical :: lANNcoulomb_

lANNcoulomb = lANNcoulomb_

return
end subroutine




subroutine adopted_ANN( nfile, myid, nodes, nc_ )
!-----------------------------------------------------------------------
!     set ladopted
!-----------------------------------------------------------------------
use io, only: nfile_loc => nfile,  &
&              myid_loc => myid,  &
&             nodes_loc => nodes
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
integer :: nc_

!-----declare local variables
integer :: status
real*8  :: the_mem


ladopted = .true.

if( nc_ <= 0 ) return

nc = nc_

!------allocate memory
allocate( ANNfiles(nc), atom_types(nc), lANNascii(nc),  &
& lANNsc(nc), ANNscepsilon(nc), ANNscsigma(nc), ANNscpower(nc),  &
& zz(nc),  &
& stat=status )

the_mem = ( 83.d0 + 25.d0 + 8d0 ) * nc

!------error trap
call check_alloc_accum( nfile, myid, nodes,  &
& status, the_mem, 'adopted_ANN', .true. )

nfile_loc(1:2) = nfile(1:2)
 myid_loc = myid
nodes_loc = nodes

lANNsc(:) = .false.
ANNscepsilon(:) = 0.d0
ANNscsigma(:) = 0.d0
ANNscpower(:) = 1.d0
zz(:) = 0.d0


return
end




subroutine set_rc_buffer_ANN( nfile, myid, nodes, rc_buffer )
!-----------------------------------------------------------------------
!     set buffer length
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
real*8  :: rc_buffer


!-----length of buffer region
rc_buffer = rc + rc_extend


return
end




subroutine set_rc_ANN( nfile, myid, nodes, rc_ )
!-----------------------------------------------------------------------
!     set cutoff length
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: myid, nodes
integer, dimension(*) :: nfile
real*8  :: rc_


!-----cutoff length
rc_ = rc


return
end




!subroutine chkdis_ANN( nfile, myid, nodes,  &
!& x, is, n, ntype, h, lfixion, xrec, imts )
!!-----------------------------------------------------------------------
!!    check atomic displacements
!!    set imts = 0, if neighbor list must be updated
!!-----------------------------------------------------------------------
!use md_potential_ANN
!implicit none
!integer :: nfile(*), myid, nodes
!integer :: n, ntype
!real*8  :: x(3,n)
!integer :: is(n)
!real*8  :: h(3,3)
!logical :: lfixion(ntype)
!real*8  :: xrec(3,n)
!integer :: imts
!
!
!call chkdis2( nfile, myid, nodes,  &
!& x, is, n, ntype, h, lfixion, xrec, imts, rc_extend )
!
!
!return
!end




subroutine md_potential_ANN_prealloc( nfile, myid, nodes, alloc_mem, nc_ )
!-----------------------------------------------------------------------
!     allocate memory for variables for MD potentials
!-----------------------------------------------------------------------
use constants
use md_potential_ANN
use aenet
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nc_

!-----declare local variables
integer :: itype, unit
integer :: stat
real*8  :: the_mem

!shimamura for s-ANN
LOGICAL,DIMENSION(nc_) :: mask
LOGICAL                :: bklANNascii
character(80)          :: bkANNfile


if( nc /= nc_ ) call fstop( nfile, myid, nodes,  &
& 'wrong # of species in md_potential_ANN_prealloc' )

!---set units of length and energy in the ANN files
if(  ANNlunit == '(ang)     ' ) then
     unitlength = audang
end if
if(  ANNeunit == '(ev)      ' ) then
     unitenergy = hrdev
else if(  ANNeunit == '(ry)      ' ) then
     unitenergy = 2.d0
end if


    !shimamura for s-ANN
!    write(*,*)"Check1, md_potential_ANN_prealloc"
!    write(*,*)"nc=",nc
!    write(*,*)"ANNlunit=",ANNlunit
!    write(*,*)"atom_types=",atom_types
!    write(*,*)"size(atom_types)=",size(atom_types)
!    write(*,*)"ANNfiles(:)=",ANNfiles(:)
!    write(*,*)"lANNascii(:)=",lANNascii(:)
    mask=(ANNfiles(1)==ANNfiles)
    efile=ALL(mask)
    !WRITE(*,*) ALL(mask)
    write(*,*)"efile = in md_potential_ANN_prealloc",efile

    if(efile)then
       nc = 1
       bkANNfile = ANNfiles(1)
       bklANNascii = lANNascii(1)
       deallocate(atom_types,ANNfiles,lANNascii)
       allocate(atom_types(1),ANNfiles(1),lANNascii(1))
       atom_types(1) = "A"
       ANNfiles(1) = bkANNfile
       lANNascii(1) = bklANNascii
    end if

    ! initialize aenet
    call aenet_init(atom_types, stat)
    if (stat /= 0) then
       write(0,*) 'Error: aenet initialization failed'
!       call finalize()
       stop
    end if


    !!shimamura for s-ANN
    !write(*,*)"Check1, md_potential_ANN_prealloc"
    !write(*,*)"nc=",nc
    !write(*,*)"ANNlunit=",ANNlunit
    !write(*,*)"atom_types=",atom_types
    !write(*,*)"ANNfiles(:)=",ANNfiles(:)
    !write(*,*)"lANNascii(:)=",lANNascii(:)
    !mask=(ANNfiles(1)==ANNfiles)
    !efile=ALL(mask)
    !!WRITE(*,*) ALL(mask)
    !write(*,*)"efile = in md_potential_ANN_prealloc",efile
        

    ! load ANN potentials
    do itype = 1, nc

       !shimamura for s-ANN
       call aenet_load_potential(itype, ANNfiles(itype), lANNascii(itype), efile,  stat)
       !call aenet_load_potential(itype, ANNfiles(itype), lANNascii(itype),  stat)

       if (stat /= 0) then
       write(0,*) 'Error: could not load ANN potentials'
!          call finalize()
          stop
       end if
    end do

  if ( myid == 0 ) call aenet_print_info()

!---set potential cutoff length
rc = aenet_Rc_max/unitlength


    !shimamura for s-ANN
!    write(*,*)"Check2, md_potential_ANN_prealloc"
    if(efile)then    
       nc = nc_
    end if

return
end




subroutine md_potential_ANN_alloc( nfile, myid, nodes,  &
& alloc_mem, nc_, nprmax, nmd_buffer, n, nmdx )
!-----------------------------------------------------------------------
!     allocate memory for variables for MD potentials
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
real*8  :: alloc_mem
integer :: nc_
integer :: nprmax
integer :: nmd_buffer
integer :: n
integer :: nmdx

!-----declare local variables
integer :: status
real*8  :: the_mem


if( nc /= nc_ ) call fstop( nfile, myid, nodes,  &
& 'wrong # of species in md_potential_ANN_alloc' )

!------allocate memory
allocate( rcij(nc,nc), rcij2(nc,nc), v(nbin,nc,nc,0:2),  &
& f3r(3,nmd_buffer),  &
& stat=status )

the_mem =  &
&  8.d0 * ( size(rcij) + size(rcij2) + size(v) + size(f3r) )

!------error trap
call check_alloc( nfile, myid, nodes,  &
& status, alloc_mem, the_mem, 'md_potential_ANN_alloc', .true. )

 rcij(:,:) = 0.d0
rcij2(:,:) = 0.d0
v(:,:,:,:) = 0.d0


return
end subroutine




subroutine potpreset_ANNatom( nfile, myid, nodes, it, zatom )
!-----------------------------------------------------------------------
!     set potential parameters
!-----------------------------------------------------------------------
use constants
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
integer :: it
real(8) :: zatom


if( .not.ladopted ) return

if( .not.allocated(atom_types) ) return

atom_types(it) = aname(nint(zatom))


return
end subroutine




subroutine potpreset_ANNfiles( nfile, myid, nodes,  &
& it, ANNfiles_, ANNlunit_, ANNeunit_, lANNascii_ )
!-----------------------------------------------------------------------
!     set potential parameters
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
integer :: it
character(80) :: ANNfiles_
character(10) :: ANNlunit_, ANNeunit_
logical :: lANNascii_


if( .not.ladopted ) return

if( .not.allocated(ANNfiles) ) return

 ANNfiles(it) =  ANNfiles_
lANNascii(it) = lANNascii_

if( ANNlunit == '' ) then
    ANNlunit = ANNlunit_
else
    !---error trap
    if( ANNlunit /= ANNlunit_ ) then
        call fstop( nfile, myid, nodes,  &
& 'wrong length unit for ANNfiles: '//trim(ANNfiles(it)) )
    end if
end if

if( ANNeunit == '' ) then
    ANNeunit = ANNeunit_
else
    !---error trap
    if( ANNeunit /= ANNeunit_ ) then
        call fstop( nfile, myid, nodes,  &
& 'wrong energy unit for ANNfiles: '//trim(ANNfiles(it)) )
    end if
end if


return
end subroutine




subroutine potpreset_ANNsc( nfile, myid, nodes,  &
& it, lANNsc_, ANNscepsilon_, ANNscsigma_, ANNscpower_ )
!-----------------------------------------------------------------------
!     set potential parameters
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
integer :: it
logical :: lANNsc_
real(8) :: ANNscepsilon_, ANNscsigma_, ANNscpower_


if( .not.ladopted ) return

if( .not.allocated(lANNsc) ) return

if( .not.lANNsc_ ) return

       lANNsccal = .true.
      lANNsc(it) =  lANNsc_
ANNscepsilon(it) = ANNscepsilon_
  ANNscsigma(it) = ANNscsigma_
  ANNscpower(it) = ANNscpower_


return
end subroutine




subroutine potpreset_ANNzv( nfile, myid, nodes, it, ANNzv )
!-----------------------------------------------------------------------
!     set potential parameters
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile(*), myid, nodes
integer :: it
real(8) :: ANNzv


if( .not.ladopted ) return

if( .not.allocated(zz) ) return

if( .not.lANNcoulomb ) return

zz(it) =  ANNzv


return
end subroutine




subroutine out_ANN_potential( nfile )
!-----------------------------------------------------------------------
!     out potential parameters
!-----------------------------------------------------------------------
use md_potential_ANN
implicit none
integer :: nfile

!-----declare local variables
integer :: it


write(nfile,*)
write(nfile,'(a)') '   units in the ANN potential files:'//  &
&                  '   length = '//trim(ANNlunit)//  &
&                  '   energy = '//trim(ANNeunit)

write(nfile,*)
write(nfile,'(a)') '       ANNscepsilon   ANNscsigma   ANNscpower'
do it = 1, nc
   write(nfile,'(i2,a,2x,es13.5,f12.6,f10.2)') it, ':',  &
& ANNscepsilon(it), ANNscsigma(it), ANNscpower(it)
end do

if( lANNcoulomb ) then
    write(nfile,*)
    write(nfile,'(a)') '       ANN charge'
    do it = 1, nc
       write(nfile,'(i2,a,2x,f12.6,f10.2)') it, ':',  zz(it)
    end do
end if


return
end subroutine




!subroutine potset_ANN( nfile, myid, nodes )
!!-----------------------------------------------------------------------
!!     prepare buffer length & potential table v
!!-----------------------------------------------------------------------
!use symmfunc, only: sf_set_table
!use md_potential_ANN
!implicit none
!integer :: nfile(*), myid, nodes
!
!!-----declare local variables
!real*8,  allocatable, dimension(:,:,:) :: vrc
!real*8,  allocatable, dimension(:,:) :: sig12, avpower
!integer :: status, ic, jc, k, l
!real*8  :: rijmax
!real*8  :: rc2, dr2, r, ri, ri2, ri3, ri4, ri12, rep, rep1, rep2, rcv
!real(8) :: facc, r2, col, col1, col2, exbd
!real*8  :: pirub, pirub2, expgr, DERFNC
!real*8  :: gamma, pi
!
!
!do ic = 1, nc
!do jc = 1, nc
!     rcij(ic,jc) = rc
!    rcij2(ic,jc) = rcij(ic,jc)* rcij(ic,jc)
!end do
!end do
!
!
!!---set function tables
!call sf_set_table
!
!
!!---set soft core potentials if necessary
!if( .not.lANNsccal .and. .not.lANNcoulomb ) return
!
!
!!------allocate memory
!allocate( vrc(nc,nc,0:2), sig12(nc,nc), avpower(nc,nc), stat=status )
!
!
!if( lANNcoulomb ) then
!    !-----Get Ewald parameters
!    call potset_rec( nfile, myid, nodes,  &
!& nc, zz, rc, gamma )
!
!    pi     = acos(-1.d0)
!    pirub  = 2.d0*gamma/sqrt(pi)
!    pirub2 = 2.d0*gamma*gamma*pirub
!end if
!
!
!!-----average potential parameters
!avpower(:,:) = 0.d0
!  sig12(:,:) = 0.d0
!do ic = 1, nc
!do jc = 1, nc
!   if( .not.lANNsc(ic) .or. .not.lANNsc(jc) ) cycle
!   avpower(ic,jc) = 0.5d0*( ANNscpower(ic) + ANNscpower(jc) )
!   sig12(ic,jc) = sqrt( ANNscepsilon(ic)*ANNscsigma(ic)**ANNscpower(ic)  &
!&                     * ANNscepsilon(jc)*ANNscsigma(jc)**ANNscpower(jc) )
!end do
!end do
!
!
!
!!-----2-body potential & derivative table, V
!
!rc2 = rc*rc
!dr2 = rc2/nbin
!dr2i= 1d0/dr2
!!-----To avoid segmentation fault of potential arrays
!rij2max=rc2-dr2
!rijmax =sqrt(rij2max)
!
!!sig12 = epsilon*sigma**power
!!-----Evaluate the potential & derivative at the truncation length, RCIJ
!do ic=1,nc
!do jc=1,nc
!   r=rcij(ic,jc)
!   ri=1d0/r
!
!   if( lANNsccal ) then
!       ri12=ri**avpower(ic,jc)
!       rep  = sig12(ic,jc)*ri12
!       rep1 = -avpower(ic,jc)*rep*ri
!   else
!       rep  = 0d0
!       rep1 = 0d0
!   end if
!
!   if( lANNcoulomb ) then
!       facc=zz(ic)*zz(jc)
!       r2 = gamma*r
!       col =facc*ri*DERFNC(r2,1.d-11)
!       col1= - (col+facc*pirub*dexp(-r2*r2))*ri
!   else
!       col  = 0d0
!       col1 = 0d0
!   end if
!
!   vrc(ic,jc,0) = rep  + col
!   vrc(ic,jc,1) = rep1 + col1
!end do
!end do
!
!
!!-----2-body potential & derivative table, V
!do k = 1, nbin
!   r = dsqrt(dble(k)*dr2)
!
!  do ic = 1, nc
!  do jc = 1, nc
!    if( r > min( rcij(ic,jc), rijmax ) ) then
!
!        do l = 0, 2
!           v(k,ic,jc,l)=0d0
!        end do
!
!    else
!
!      ri=1d0/r
!      ri2=ri*ri
!      ri3=ri2*ri
!
!      if( lANNsccal ) then
!          ri4=ri2*ri2
!          ri12=ri**avpower(ic,jc)
!          rep  = sig12(ic,jc)*ri12
!          rep1 = avpower(ic,jc)*rep*ri2
!          rep2 = -avpower(ic,jc)*(avpower(ic,jc)+2.d0)*rep*ri4
!      else
!          rep  = 0d0
!          rep1 = 0d0
!          rep2 = 0d0
!      end if
!
!      if( lANNcoulomb ) then
!          facc=zz(ic)*zz(jc)
!          r2 = gamma*r
!          exbd=DERFNC(r2,1.d-11)
!          expgr = dexp(-r2*r2)
!          col  = facc*ri*exbd
!          col1 = facc*(ri*exbd+pirub*expgr)*ri2
!          col2 =-facc*(3d0*ri2*exbd+3d0*ri*pirub*expgr+pirub2*r*expgr)*ri3
!      else
!          col  = 0d0
!          col1 = 0d0
!          col2 = 0d0
!      end if
!
!      rcv = rcij(ic,jc)
!
!      !-----2-body potential
!      v(k,ic,jc,0) = rep + col
!      !-----Truncate the potential at RCIJ(IC,JC)
!      v(k,ic,jc,0)=v(k,ic,jc,0)-vrc(ic,jc,0)-(r-rcv)*vrc(ic,jc,1)
!
!      !-----the 1st derivative / r
!      v(k,ic,jc,1) = rep1 + col1
!      !-----Truncate the 1st derivative at RCIJ(IC,JC)
!      v(k,ic,jc,1)=v(k,ic,jc,1)+vrc(ic,jc,1)*ri
!
!      !-----the 2nd derivative
!      v(k,ic,jc,2) = rep2 + col2
!      !-----------Truncate 2nd derivative at RCIJ(IC,JC)
!      v(k,ic,jc,2)=v(k,ic,jc,2)-vrc(ic,jc,1)*ri**3
!
!    endif
!
!  end do
!  end do
!
!!-----Enddo potential grids, R
!end do
!
!!-----deallocate memory
!deallocate( vrc, sig12, avpower, stat=status )
!
!
!return
!end




!subroutine accel_ANN( nfile, myid, nodes,  &
!& epot, x, is, n, npb, a, h, nmdx, lsprlr, nprmax, nmd_buffer,  &
!& ba_call, pintlr, pintsr, lcstress, lcatomic, wepot, wstrs, nmdmax__ )
!!-----------------------------------------------------------------------
!!     classical energy and forces
!!-----------------------------------------------------------------------
!use aenet
!use md_potential_ANN
!implicit none
!integer :: nfile(*), myid, nodes
!real*8  :: epot
!integer :: n, npb, nmdx
!real*8  :: x(3,nmdx)
!integer :: is(nmdx)
!real*8  :: a(3,n)
!real*8  :: h(3,3)
!integer :: nprmax
!integer :: lsprlr(0:nprmax,n)
!integer :: nmd_buffer
!logical :: ba_call
!real*8  :: pintlr(3,3), pintsr(3,3)
!logical :: lcstress, lcatomic
!integer :: nmdmax__
!real*8  :: wepot(nmdmax__), wstrs(6,nmdmax__)
!
!!-----declare local variables
!integer :: i, j, k, ic, jc, ip, ir, ia
!real*8  :: sxi, syi, szi
!real*8  :: xi, yi, zi, xj, yj, zj
!real*8  :: xij, yij, zij
!real*8  :: rij2
!integer :: type_i, iatom, nnb
!real(8) :: coo_i(3)
!real*8  :: fr
!real*8  :: v0, v1, fx, fy, fz
!real*8,  dimension(6) :: pit1lr, pit1sr, dbuf
!real*8,  dimension(6) :: rrij, raux2
!real(8) :: Ecoh, E_i, epotANN
!logical :: bintra, lANNshort
!integer :: stat
!
!
!     !!shimamura for s-ANN
!     !write(*,*)"efile= in accel_ANN",efile 
!
!
!    !---memory allocation if necessary
!    if( aenet_nnb_max > aenet_nnb_maxx ) then
!        !-----if already allocated, deallocate arrays
!        if( allocated(nblist) ) then
!            deallocate( nblist, nbcoo, nbdist, nbtype )
!        end if
!        aenet_nnb_maxx = aenet_nnb_max
!        !------allocate memory
!        allocate( nblist(aenet_nnb_maxx), nbcoo(3,aenet_nnb_maxx),  &
!&                 nbdist(aenet_nnb_maxx), nbtype(aenet_nnb_maxx) )
!    end if
!
!
!lANNshort = lANNsccal .or. lANNcoulomb
!
!!-----Reset potential energy, forces
!epot = 0.d0
!epotANN = 0.d0
!  a(1:3,1:n) = 0.d0
!f3r(1:3,1:nmd_buffer) = 0.d0
!!-----Reset pressure tensor, PIT
!pit1lr(:) = 0d0
!pit1sr(:) = 0d0
! wepot(:) = 0d0
!wstrs(:,:)= 0d0
!
!Ecoh = 0.d0
!
!  atomi: do i = 1, n
!     !-----The species & coordinates of atom i
!     ic  = is(i)
!     sxi = x(1,i)
!     syi = x(2,i)
!     szi = x(3,i)
!
!     !-----Physical coordinate XI = H*SXI
!     xi = h(1,1)*sxi+h(1,2)*syi+h(1,3)*szi
!     yi = h(2,1)*sxi+h(2,2)*syi+h(2,3)*szi
!     zi = h(3,1)*sxi+h(3,2)*syi+h(3,3)*szi
!
!!---coo_i  ... real coordinate of the central atom
!!---type_i ... type            of the central atom
!!---iatom  ... atom id         of the central atom
!!---nnb    ... # of neighbors around the central atom
!!---nbcoo  ... real coordinates of neighbors around the central atom
!!---nbtype ... type             of neighbors around the central atom
!!---nblist ... atom id          of neighbors around the central atom
!     coo_i(1) = xi
!     coo_i(2) = yi
!     coo_i(3) = zi
!     iatom = i
!
!
!     !shimamura for s-ANN
!     if(efile)then
!        type_i = 1
!     else
!        type_i = ic
!     end if
!
!
!     nnb = 0
!     !-----Scan pair atom j
!     atomj: do ip = 1, lsprlr(0,i)
!        j = lsprlr(ip,i)
!        !-----The species & address of atom j
!        jc  = is(j)
!        sxi = x(1,j)
!        syi = x(2,j)
!        szi = x(3,j)
!
!        !-----Physical coordinate XJ = H*SXI
!        xj = h(1,1)*sxi+h(1,2)*syi+h(1,3)*szi
!        yj = h(2,1)*sxi+h(2,2)*syi+h(2,3)*szi
!        zj = h(3,1)*sxi+h(3,2)*syi+h(3,3)*szi
!
!        !-----Physical pair vector XIJ = H*SXJ
!        xij = xi - xj
!        yij = yi - yj
!        zij = zi - zj
!        rij2= xij*xij + yij*yij + zij*zij
!
!        if( rij2 < rcij2(ic,jc) ) then
!            !---neighbor atom list for ANN potential
!            nnb = nnb + 1
!            nbcoo(1,nnb) = xj
!            nbcoo(2,nnb) = yj
!            nbcoo(3,nnb) = zj
!            nbtype(nnb) = jc
!            nblist(nnb) = j
!
!            !-----Calculate 2-body forces for intranode pairs i < j, and all
!            !-----the copied partners, if necessary
!
!!            if( j > i .and. rij2 < rcij2(ic,jc) ) then
!            if( lANNshort .and.  j > i ) then
!
!                !-----To avoid segmentation fault of V
!                rij2=dmin1(rij2,rij2max)
!                ir=(rij2*dr2i)
!                fr=(rij2*dr2i)-ir
!                v1=(1d0-fr)*v(ir,ic,jc,1)+fr*v(ir+1,ic,jc,1)
!                !---unit change
!                v1 = v1 /unitlength*unitenergy
!
!                fx = v1*xij
!                fy = v1*yij
!                fz = v1*zij
!
!                !-----Store 2-body forces on atom i in A
!                a(1,i) = a(1,i) + fx
!                a(2,i) = a(2,i) + fy
!                a(3,i) = a(3,i) + fz
!
!                !-----Calculate (partial) two-body potential
!                v0=(1d0-fr)*v(ir,ic,jc,0)+fr*v(ir+1,ic,jc,0)
!                !---unit change
!                v0 = v0*unitenergy
!                if( lcatomic ) wepot(i) = wepot(i) + 0.5d0*v0
!
!                bintra = j >= 1 .and. j <= n
!                if( bintra ) then
!                    !-----Store 2-body forces on atom j if a resident 
!                    a(1,j) = a(1,j) - fx
!                    a(2,j) = a(2,j) - fy
!                    a(3,j) = a(3,j) - fz
!                    if( lcatomic ) wepot(j) = wepot(j) + 0.5d0*v0
!                    epot = epot + v0
!                else
!                    epot = epot + 0.5d0*v0
!                end if
!
!                if( lcstress ) then
!                    rrij(1)=xij*xij
!                    rrij(2)=yij*yij
!                    rrij(3)=zij*zij
!                    rrij(4)=yij*zij
!                    rrij(5)=zij*xij
!                    rrij(6)=xij*yij
!
!                    !---unit change for stress
!                    v1 = v1 *unitlength
!
!                    rrij(1:6)  = rrij(1:6)*v1
!                    raux2(1:6) = 0.5d0*rrij(1:6)
!                    if( lcatomic ) wstrs(1:6,i) = wstrs(1:6,i) + raux2(1:6)
!
!                    if( bintra ) then
!                        !-----Store 2-body forces on atom j if a resident 
!                        pit1lr(1:6)=pit1lr(1:6)+rrij(1:6)
!                        if( lcatomic ) wstrs(1:6,j) = wstrs(1:6,j) + raux2(1:6)
!                    else
!                        pit1lr(1:6)=pit1lr(1:6)+raux2(1:6)
!                    end if
!                end if
!
!            !-----Endif i<j
!            end if
!
!        end if
!
!     !-----Enddo pair atom j
!     end do atomj
!
!     !---bohr -> A
!     coo_i(1:3) = coo_i(1:3)*unitlength
!     nbcoo(1:3,1:nnb) = nbcoo(1:3,1:nnb)*unitlength
!
!     if( .not.lcstress ) then
!         call aenet_atomic_energy_and_forces( &
!               coo_i, type_i, iatom, nnb, nbcoo, nbtype, nblist, &
!               n, E_i, a, f3r, nmd_buffer, stat)
!     else
!         call aenet_atomic_energy_and_forces_and_stresses( &
!               coo_i, type_i, iatom, nnb, nbcoo, nbtype, nblist, &
!               n, E_i, a, f3r, nmd_buffer,  &
!&              pit1lr, pit1sr, lcatomic, wstrs(1,min(i,nmdmax__)), stat)
!     end if
!
!     epotANN = epotANN + E_i
!        Ecoh =    Ecoh + E_i - aenet_free_atom_energy(type_i)
!
!     if( lcatomic ) wepot(i) = wepot(i) + E_i
!
!  end do atomi
!
!!---eV/A -> hartree
!epot    = epot   /unitenergy
!epotANN = epotANN/unitenergy
!a(1:3,1:n) = a(1:3,1:n)*unitlength/unitenergy
!f3r(1:3,1:nmd_buffer) = f3r(1:3,1:nmd_buffer)*unitlength/unitenergy
!if( lcstress ) then
!    !---Caution! pit1lr and pit1sr will be devided by volume in [a.u.]
!    !---         in subroutine get_force.
!    pit1lr(1:6) = pit1lr(1:6)/unitenergy !*unitlength**3  
!    pit1sr(1:6) = pit1sr(1:6)/unitenergy !*unitlength**3
!end if
!if( lcatomic ) then
!    wepot(:) = wepot(:)/unitenergy
!    if( lcstress ) then
!        wstrs(:,:) = wstrs(:,:)/unitenergy !*unitlength**3
!    end if
!end if
!
!!-----Send & receive reaction on cached atoms : f3r
!!      if( imts /= 0 ) then 
!
!if( npb > 0 .and. ba_call ) then
!    call get_reaction( nfile, myid, nodes,  &
!& a, n, npb, f3r, nmd_buffer )
!end if
!
!
!call gdsum(pit1lr,6,dbuf)
!pintlr(1,1)=pit1lr(1)
!pintlr(2,2)=pit1lr(2)
!pintlr(3,3)=pit1lr(3)
!pintlr(2,3)=pit1lr(4)
!pintlr(3,1)=pit1lr(5)
!pintlr(1,2)=pit1lr(6)
!pintlr(3,2)=pintlr(2,3)
!pintlr(1,3)=pintlr(3,1)
!pintlr(2,1)=pintlr(1,2)
!!          call matscl(pintlr,1d0/volume,pintlr)
!
!!      call gdsum(pit1sr,6,dbuf)
!pintsr(1,1)=pit1sr(1)
!pintsr(2,2)=pit1sr(2)
!pintsr(3,3)=pit1sr(3)
!pintsr(2,3)=pit1sr(4)
!pintsr(3,1)=pit1sr(5)
!pintsr(1,2)=pit1sr(6)
!pintsr(3,2)=pintsr(2,3)
!pintsr(1,3)=pintsr(3,1)
!pintsr(2,1)=pintsr(1,2)
!!      call matscl(pintsr,1d0/volume,pintsr)
!
!
!if( lANNcoulomb ) then
!    !----- direct coulomb interaction in reciprocal space
!    call accel_rec( nfile, myid, nodes,  &
!& epot, x, is, n, a, h, nmdx, pintlr,  &
!& lcstress, lcatomic, wepot, wstrs, nmdmax__ )
!end if
!
!!v1 = epot
!!call gdsum(v1,1,dbuf)
!!if( myid == 0 ) write(*,*) 'Ecoulomb=', v1
!
!epot = epot + epotANN
!
!
!return
!end subroutine

subroutine fstop( nfile, myid, nodes, message )
!-----------------------------------------------------------------------
!    error trap: display message and stop program
!-----------------------------------------------------------------------
implicit none
integer :: nfile(*), myid, nodes
character(*) :: message
integer :: ierr

return
end subroutine

subroutine deallocate_unit_number( iunit )
!-----------------------------------------------------------------------
!  deallocate unit number
!-----------------------------------------------------------------------
implicit none
integer :: iunit
return
end subroutine

subroutine allocate_unit_number( iunit )
!-----------------------------------------------------------------------
!  deallocate unit number
!-----------------------------------------------------------------------
implicit none
integer :: iunit
integer,save :: current_unit = 900

iunit = current_unit
current_unit = current_unit + 1

return
end subroutine

subroutine check_alloc_accum( nfile, myid, nodes, status, the_mem, sub_name, lstop )
implicit none
integer :: nfile(*), myid, nodes
integer :: status
real*8  :: the_mem
character(*) :: sub_name
logical :: lstop
return
end subroutine


subroutine check_alloc( nfile, myid, nodes,  status, alloc_mem, the_mem, sub_name, lstop )
implicit none
integer :: nfile(*), myid, nodes
integer :: status
real*8  :: alloc_mem, the_mem
character(*) :: sub_name
logical :: lstop
return
end subroutine

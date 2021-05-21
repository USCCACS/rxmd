module nompi

   use iso_fortran_env, only: int32, int64, real32, real64, real128

   implicit none

   integer,parameter :: NOMPI_MAXSTRLENGTH=256 

   integer,parameter :: mpi_status_size = 1
   integer,parameter :: mpi_in_place = 999999

   integer :: mpi_sum=1, mpi_max=2, mpi_min=3
 

   integer :: mpi_comm_world, mpi_info, mpi_file, mpi_datatype, mpi_status
   integer :: mpi_mode_rdonly, mpi_mode_wronly, mpi_mode_create, mpi_seek_set
   integer :: mpi_info_null, mpi_status_ignore

   integer,parameter :: mpi_double = real64
   integer,parameter :: mpi_comm = int32
   integer,parameter :: mpi_character = int32
   integer,parameter :: mpi_integer = int32
   integer,parameter :: mpi_integer8 = int64
   integer,parameter :: mpi_double_precision = real64
   integer,parameter :: mpi_offset_kind = real64

   integer :: nompi_file_handler = 0
   character(NOMPI_MAXSTRLENGTH) :: nompi_file_path

   integer :: nompi_sbuf_scalar_i4, nompi_rbuf_scalar_i4
   real(8) :: nompi_sbuf_scalar_r8, nompi_rbuf_scalar_r8
   integer,allocatable,dimension(:) :: nompi_sbuf_vector_i4, nompi_rbuf_vector_i4
   real(8),allocatable,dimension(:) :: nompi_sbuf_vector_r8, nompi_rbuf_vector_r8

   interface mpi_file_write
      module procedure mpi_file_write_r8, mpi_file_write_i4, mpi_file_write_c
   end interface

   interface mpi_file_read
      module procedure mpi_file_read_r8, mpi_file_read_i4
   end interface

   interface mpi_bcast
      module procedure mpi_bcast_scalar_i4, mpi_bcast_scalar_r8, mpi_bcast_i4, mpi_bcast_r8
   end interface

   interface mpi_send
      module procedure mpi_send_vector_i4, mpi_send_scalar_i4
      module procedure mpi_send_vector_r8, mpi_send_scalar_r8
   end interface

   interface mpi_recv
      module procedure mpi_recv_vector_i4, mpi_recv_scalar_i4
      module procedure mpi_recv_vector_r8, mpi_recv_scalar_r8
   end interface

   interface mpi_allreduce
      module procedure mpi_allreduce_scalar_i8, mpi_allreduce_1d_i8
      module procedure mpi_allreduce_scalar_i4, mpi_allreduce_1d_i4, mpi_allreduce_2d_i4, mpi_allreduce_1d_i4_minmax
      module procedure mpi_allreduce_scalar_r8, mpi_allreduce_1d_r8, mpi_allreduce_2d_r8
   end interface

contains

!all mpi objects (e.g., mpi_datatype, mpi_comm) are of type integer in fortran.
!https://www.mpich.org/static/docs/v3.2/www3/mpi_comm_rank.html

   subroutine mpi_init( ierr )
     integer :: ierr
   end subroutine

   subroutine mpi_comm_rank(comm, rank, ierror)
     integer, intent(in) :: comm
     integer, intent(out) :: rank
     integer, optional, intent(out) :: ierror
     rank = 0
   end subroutine

   subroutine mpi_comm_size(comm, size, ierror)
     integer, intent(in) :: comm
     integer, intent(out) :: size
     integer, optional, intent(out) :: ierror
     size = 1
   end subroutine

   subroutine mpi_finalize(ierror)
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_barrier(comm, ierror)
     integer, intent(in) :: comm
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_send_scalar_i4(buf, icount, datatype, dest, tag, comm, ierror)
     integer, intent(in) :: buf
     integer, intent(in) :: icount, dest, tag, comm, datatype
     integer, optional, intent(out) :: ierror

     nompi_sbuf_scalar_i4 = buf
   end subroutine

   subroutine mpi_send_scalar_r8(buf, icount, datatype, dest, tag, comm, ierror)
     real(8), intent(in) :: buf
     integer, intent(in) :: icount, dest, tag, comm, datatype
     integer, optional, intent(out) :: ierror

     nompi_sbuf_scalar_r8 = buf
   end subroutine

   subroutine mpi_send_vector_i4(buf, icount, datatype, dest, tag, comm, ierror)
     integer, dimension(:), intent(in) :: buf
     integer, intent(in) :: icount, dest, tag, comm, datatype
     integer, optional, intent(out) :: ierror
     
     nompi_sbuf_vector_i4 = buf
   end subroutine

   subroutine mpi_send_vector_r8(buf, icount, datatype, dest, tag, comm, ierror)
     real(8), dimension(:), intent(in) :: buf
     integer, intent(in) :: icount, dest, tag, comm, datatype
     integer, optional, intent(out) :: ierror
     
     nompi_sbuf_vector_r8 = buf
   end subroutine

   subroutine mpi_recv_vector_i4(buf, icount, datatype, source, tag, comm, status, ierror)
     integer, dimension(:) :: buf
     integer, intent(in) :: icount, source, tag, comm, datatype
     integer :: status(mpi_status_size)
     integer, optional, intent(out) :: ierror

     buf = nompi_sbuf_vector_i4
   end subroutine

   subroutine mpi_recv_vector_r8(buf, icount, datatype, source, tag, comm, status, ierror)
     real(8), dimension(:) :: buf
     integer, intent(in) :: icount, source, tag, comm, datatype
     integer :: status(mpi_status_size)
     integer, optional, intent(out) :: ierror

     buf = nompi_sbuf_vector_r8
   end subroutine

   subroutine mpi_recv_scalar_i4(buf, icount, datatype, source, tag, comm, status, ierror)
     integer :: buf
     integer, intent(in) :: icount, source, tag, comm, datatype
     integer :: status(mpi_status_size)
     integer, optional, intent(out) :: ierror

     buf = nompi_sbuf_scalar_i4
   end subroutine

   subroutine mpi_recv_scalar_r8(buf, icount, datatype, source, tag, comm, status, ierror)
     real(8) :: buf
     integer, intent(in) :: icount, source, tag, comm, datatype
     integer :: status(mpi_status_size)
     integer, optional, intent(out) :: ierror

     buf = nompi_sbuf_scalar_r8
   end subroutine

   subroutine mpi_allreduce_scalar_i8(inplace, buf, icount, datatype, mpiop, comm, ierror)
     integer :: inplace
     integer(8) :: buf
     integer, intent(in) :: icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine

   subroutine mpi_allreduce_1d_i8(inplace, buf, icount, datatype, mpiop, comm, ierror)
     integer(8), dimension(:) :: buf
     integer, intent(in) :: inplace, icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine


   subroutine mpi_allreduce_scalar_i4(sbuf, rbuf, icount, datatype, mpiop, comm, ierror)
     integer :: sbuf, rbuf
     integer, intent(in) :: icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     if (sbuf == mpi_in_place) return
     rbuf = sbuf
   end subroutine

   subroutine mpi_allreduce_1d_i4(inplace, buf, icount, datatype, mpiop, comm, ierror)
     integer, dimension(:) :: buf
     integer, intent(in) :: inplace, icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine

   subroutine mpi_allreduce_1d_i4_minmax(sbuf, rbuf, icount, datatype, mpiop, comm, ierror)
     integer, dimension(:) :: rbuf, sbuf
     integer, intent(in) :: icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     rbuf = sbuf
   end subroutine

   subroutine mpi_allreduce_2d_i4(inplace, buf, icount, datatype, mpiop, comm, ierror)
     integer, dimension(:,:) :: buf
     integer, intent(in) :: inplace, icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine


   subroutine mpi_allreduce_1d_r8(inplace, buf, icount, datatype, mpiop, comm, ierror)
     real(8), dimension(:) :: buf
     integer, intent(in) :: inplace, icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine

   subroutine mpi_allreduce_2d_r8(inplace, buf, icount, datatype, mpiop, comm, ierror)
     real(8), dimension(:,:) :: buf
     integer, intent(in) :: inplace, icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine

   subroutine mpi_allreduce_scalar_r8(inplace, buf, icount, datatype, mpiop, comm, ierror)
     real(8) :: buf
     integer, intent(in) :: inplace, icount, mpiop, comm, datatype
     integer, optional, intent(out) :: ierror

     return
   end subroutine


   subroutine mpi_scan(sendbuf, recvbuf, icount, datatype, op, comm, ierror)
     integer :: sendbuf, recvbuf
     integer :: icount, datatype, op, comm, ierror
   end subroutine

   subroutine mpi_probe(source, tag, comm, status, ierror)
     integer :: source, tag, comm, status(mpi_status_size), ierror
   end subroutine

   subroutine mpi_get_count(status, datatype, icount, ierror)
     integer :: status(mpi_status_size), datatype, icount, ierror
   end subroutine

   subroutine mpi_bcast_scalar_i4(buffer, icount, datatype, root, comm, ierror)
     integer, intent(in) :: buffer
     integer, intent(in) :: icount, root
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_bcast_scalar_r8(buffer, icount, datatype, root, comm, ierror)
     real(8), intent(in) :: buffer
     integer, intent(in) :: icount, root
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_bcast_i4(buffer, icount, datatype, root, comm, ierror)
     integer, intent(in) :: buffer(*)
     integer, intent(in) :: icount, root
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_bcast_r8(buffer, icount, datatype, root, comm, ierror)
     real(8), intent(in) :: buffer(*)
     integer, intent(in) :: icount, root
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_set_size(fh, size, ierror)
     integer, intent(in) :: fh
     integer(kind=mpi_offset_kind), intent(in) :: size
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_open(comm, filename, amode, info, fh, ierror)
     integer, intent(in) :: comm
     character(len=*), intent(in) :: filename
     integer, intent(in) :: amode
     integer, intent(in) :: info, fh
     integer, optional, intent(out) :: ierror

     open(newunit=nompi_file_handler, file=trim(filename), form='unformatted', access='stream')
   end subroutine

   subroutine mpi_file_write_r8(fh, buf, icount, datatype, status, ierror)
     real(8), intent(in) :: buf(icount)
     integer, intent(in) :: fh, icount, datatype, status
     integer, optional, intent(out) :: ierror

     write(nompi_file_handler) buf
   end subroutine

   subroutine mpi_file_write_i4(fh, buf, icount, datatype, status, ierror)
     integer, intent(in) :: buf(icount)
     integer, intent(in) :: fh, icount, datatype, status
     integer, optional, intent(out) :: ierror

     write(nompi_file_handler) buf
   end subroutine

   subroutine mpi_file_write_c(fh, buf, icount, datatype, status, ierror)
     integer, intent(in) :: fh, icount, datatype, status
     character(len=icount), intent(in) :: buf
     integer, optional, intent(out) :: ierror

     write(nompi_file_handler) buf
   end subroutine

   subroutine mpi_file_read_i4(fh, buf, icount, datatype, status, ierror)
     integer, dimension(icount) :: buf
     integer, intent(in) :: fh, icount, datatype, status 
     integer, optional, intent(out) :: ierror

     read(nompi_file_handler) buf(:)
   end subroutine

   subroutine mpi_file_read_r8(fh, buf, icount, datatype, status, ierror)
     real(8), dimension(icount) :: buf
     integer, intent(in) :: fh, icount, datatype, status 
     integer, optional, intent(out) :: ierror

     read(nompi_file_handler) buf
   end subroutine

   subroutine mpi_file_seek(fh, offset, whence, ierror)
     integer :: fh, whence, ierror
     integer(kind=mpi_offset_kind) :: offset
   end subroutine

   subroutine mpi_file_close(fh, ierror)
     integer,  intent(inout) :: fh
     integer, optional, intent(out) :: ierror

     close(nompi_file_handler)
   end subroutine

   real(8) function mpi_wtime() result(t)
     call cpu_time(t) 
   end function

   subroutine mpi_get_processor_name(hostname, resultlen, ierr)
      character(*) :: hostname
      integer :: resultlen, ierr
      hostname = 'localhost'
      resultlen = 8
      ierr = 0
   end subroutine

end module

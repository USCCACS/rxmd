!-----------------------------------------------------------------------------------------
program main
  use voxels_mod
!-----------------------------------------------------------------------------------------
  implicit none
  type(string_array),allocatable :: filenames(:)
  type(mdframe) :: oneframe
  type(mdframe),allocatable :: mdframes(:)
  type(analysis_context) :: ac
  type(voxel),allocatable :: voxels(:,:,:)
  real(8) :: cell_size(3)
  integer :: num_cells(3)

  integer :: n, i,j,k, i1,j1,k1, ity, jty
  character(len=256) :: argv
  character(len=:),allocatable :: dirpath
  logical :: is_there

  ! glob all files
  call getarg(1,argv)
  dirpath = trim(adjustl(argv))
  print'(a)', 'input data path: '//dirpath
  inquire(file=dirpath, exist=is_there) 
  if(.not. is_there) stop 'ERROR: input dir does not exist.'

  filenames = get_filelist_xyz(dirpath)

  ! allocate z-element array
  allocate(mdframes(0))
  do i=1, size(filenames)
     mdframes = [mdframes, get_mdframe_from_xyzfile(filenames(i)%str)]
  enddo

  ! pick one frame to setup parameters
  oneframe = mdframes(1)
  ac = get_analysis_context_from_mdframe(oneframe)
  call get_voxel_params_from_mdframe(oneframe, num_cells, cell_size)

  !voxels = get_voxels_from_mdcontext(oneframe, num_cells, ac%elems)
  !call save_voxels_into_xyz(voxels,ac%elems,'foo.xyz')

  ! receiver needs to be indexed with the same range as passer
  allocate(voxels(DIMS(1):DIMS(4),DIMS(2):DIMS(5),DIMS(3):DIMS(6)))

  do n=1, size(mdframes)

     ! increment the number of sampled mdframes
     ac%num_sample_frames = ac%num_sample_frames + 1

     voxels = get_voxels_from_mdcontext(mdframes(n), num_cells, ac%elems)
    
     do i=DIMS0(1), DIMS0(4)
     do j=DIMS0(2), DIMS0(5)
     do k=DIMS0(3), DIMS0(6)
        do i1=-1, 1
        do j1=-1, 1
        do k1=-1, 1
           !print'(6i4,3x,3i4)',i,j,k, i1,j1,k1, i+i1,j+j1,k+k1
           call compute_interatomic_distance(ac, voxels(i,j,k), voxels(i+i1,j+j1,k+k1))
        enddo; enddo; enddo
     enddo; enddo; enddo

  enddo

  call ac%save_stat()

end program

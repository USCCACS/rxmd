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

  integer :: n, i,j,k,l, i1,j1,k1, ity,jty,kty, idx, idx2
  real(8) :: rr(0:3), rr1(0:3), rr2(0:3), cosine, theta
  character(len=256) :: argv
  character(len=:),allocatable :: dirpath
  logical :: is_there

  type(nbrlist_type),allocatable :: nbrlist(:)

  ! glob all files
  call getarg(1,argv)
  dirpath = trim(adjustl(argv))
  print'(a)', 'input data path: '//dirpath
  inquire(file=dirpath, exist=is_there) 
  if(.not. is_there) stop 'ERROR: input dir does not exist.'

  filenames = get_filelist_xyz(dirpath)

  ! pick one frame to setup parameters
  oneframe = get_mdframe_from_xyzfile(filenames(1)%str)
  ac = get_analysis_context_from_mdframe(oneframe)

  ! load all mdframes
  allocate(mdframes(0))
  do i=1, size(filenames)
     mdframes = [mdframes, get_mdframe_from_xyzfile(filenames(i)%str, ac=ac)]
  enddo

  call get_voxel_params_from_mdframe(oneframe, num_cells, cell_size)

  !voxels = get_voxels_from_mdcontext(oneframe, num_cells, ac%elems)
  !call save_voxels_into_xyz(voxels,ac%elems,'foo.xyz')

  ! receiver needs to be indexed with the same range as passer
  allocate(voxels(DIMS(1):DIMS(4),DIMS(2):DIMS(5),DIMS(3):DIMS(6)))

  do n=1, size(mdframes)
     print*,'filename: ', mdframes(n)%filename

     ! increment the number of sampled mdframes
     ac%num_sample_frames = ac%num_sample_frames + 1

     nbrlist = nbrlist_ctor_from_mdframe(mdframes(n))

     voxels = get_voxels_from_mdcontext(mdframes(n), num_cells, ac%elems)
    
     ! reset neighbor list
     do i=DIMS0(1), DIMS0(4)
     do j=DIMS0(2), DIMS0(5)
     do k=DIMS0(3), DIMS0(6)
        do i1=-1, 1
        do j1=-1, 1
        do k1=-1, 1
           !print'(6i4,3x,3i4)',i,j,k, i1,j1,k1, i+i1,j+j1,k+k1
           call compute_interatomic_distance(nbrlist, voxels(i,j,k), voxels(i+i1,j+j1,k+k1))
        enddo; enddo; enddo
     enddo; enddo; enddo

     do i = 1, size(nbrlist)

        associate(ia=>nbrlist(i)) 
        ity = ia%itype
        do j = 1, ia%counter
          jty = ia%nbrs(j)%itype
          idx = ia%nbrs(j)%ir
          ac%gr(ity,jty,idx)=ac%gr(ity,jty,idx)+1.d0
        enddo

        do j = 1, ia%counter-1
          if(ia%nbrs(j)%rr>BA_CUTOFF) cycle
          jty = ia%nbrs(j)%itype

          rr1(1:3) = ia%nbrs(j)%pos(1:3)-ia%pos(1:3)
          rr1(0) = sqrt(sum(rr1(1:3)*rr1(1:3)))
          rr1(1:3) = rr1(1:3)/rr1(0)

          do k = j+1, ia%counter
             if(ia%nbrs(k)%rr>BA_CUTOFF) cycle
             kty = ia%nbrs(k)%itype

             rr2(1:3) = ia%nbrs(k)%pos(1:3)-ia%pos(1:3)
             rr2(0) = sqrt(sum(rr2(1:3)*rr2(1:3)))
             rr2(1:3) = rr2(1:3)/rr2(0)
        
             cosine = sum(rr1(1:3)*rr2(1:3))
             theta = acos(cosine)*180d0/pi
             idx = int(theta) + 1

             do l = 1, NUM_BA
               if(ia%nbrs(j)%rr<ac%rc_ba(l) .and. ia%nbrs(k)%rr<ac%rc_ba(l)) then
                 ac%ba(jty,ity,kty,l,idx) = ac%ba(jty,ity,kty,l,idx)+0.5d0
                 ac%ba(kty,ity,jty,l,idx) = ac%ba(kty,ity,jty,l,idx)+0.5d0
               endif
             enddo 

           enddo

        enddo
        end associate
     enddo

     deallocate(nbrlist)

  enddo

  call ac%save_stat()

end program

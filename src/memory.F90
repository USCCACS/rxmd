!-------------------------------------------------------------------------------------------
module memory_allocator_mod
!-------------------------------------------------------------------------------------------
   implicit none

   integer :: status
   integer(8) :: totalMemory=0

   interface allocator
      module procedure AllocatorD1D, AllocatorD2D, AllocatorD3D, AllocatorI1D, &
                       AllocatorI2D, AllocatorI3D, AllocatorI81D, AllocatorF2D, AllocatorF3D
   end interface

   interface deallocator
      module procedure DeallocatorD1D, DeallocatorD2D, DeallocatorD3D, DeallocatorI1D, &
                       DeallocatorI2D, DeallocatorI3D, DeallocatorI81D, DeallocatorF2D, DeallocatorF3D
   end interface

contains 

subroutine AllocatorD1D(array, imin, imax)
  integer,intent(in) :: imin, imax
  real(8),allocatable,dimension(:) :: array
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*8

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorD1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorD1D(array)
  real(8),allocatable,dimension(:) :: array
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorD1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorD2D(array, imin1, imax1, imin2, imax2) 
  integer,intent(in) :: imin1, imax1, imin2, imax2
  real(8),allocatable,dimension(:,:) :: array
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*8

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorD2D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorD2D(array) 
  real(8),allocatable,dimension(:,:) :: array
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)

  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorD2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorD3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  real(8),allocatable,dimension(:,:,:) :: array
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*8

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorD3D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorD3D(array) 
  real(8),allocatable,dimension(:,:,:) :: array
  
  totalMemory = totalMemory - size(array)*8
  deallocate(array, stat=status)
  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorD3D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI81D(array, imin, imax) 
  integer,intent(in) :: imin, imax
  integer(8),allocatable,dimension(:) :: array
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*4

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorI1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI81D(array) 
  integer(8),allocatable,dimension(:) :: array
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorI1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI1D(array, imin, imax) 
  integer,intent(in) :: imin, imax
  integer,allocatable,dimension(:) :: array
  
  allocate(array(imin:imax), stat=status)
  totalMemory = totalMemory + size(array)*4

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorI1D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI1D(array) 
  integer,allocatable,dimension(:) :: array
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorI1D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI2D(array, imin1, imax1, imin2, imax2) 
  integer,intent(in) :: imin1, imax1, imin2, imax2
  integer,allocatable,dimension(:,:) :: array
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*4

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorI2D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI2D(array) 
  integer,allocatable,dimension(:,:) :: array
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorI2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorI3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  integer,allocatable,dimension(:,:,:) :: array
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*4

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorI3D: totalMemory = ', totalMemory, status

  return 
end subroutine 

subroutine DeallocatorI3D(array) 
  integer,allocatable,dimension(:,:,:) :: array
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)
  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorI3D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine AllocatorF2D(array, imin1, imax1, imin2, imax2) 
  integer,intent(in) :: imin1, imax1, imin2, imax2
  real(4),allocatable,dimension(:,:) :: array
  
  allocate(array(imin1:imax1,imin2:imax2), stat=status)
  totalMemory = totalMemory + size(array)*4

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorF2D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine AllocatorF3D(array, imin1, imax1, imin2, imax2, imin3, imax3) 
  integer,intent(in) :: imin1, imax1, imin2, imax2, imin3, imax3
  real(4),allocatable,dimension(:,:,:) :: array
  
  allocate(array(imin1:imax1,imin2:imax2,imin3:imax3), stat=status)
  totalMemory = totalMemory + size(array)*4

  array = 0

  if(status/=0) print'(a,i9,i3)', 'ERROR in AllocatorF3D: totalMemory = ', totalMemory, status

  return 
end subroutine

subroutine DeallocatorF2D(array) 
  real(4),allocatable,dimension(:,:) :: array
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)

  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorF2D: totalMemory = ', totalMemory, status

  return
end subroutine 

subroutine DeallocatorF3D(array) 
  real(4),allocatable,dimension(:,:,:) :: array
  
  totalMemory = totalMemory - size(array)*4
  deallocate(array, stat=status)

  if(status/=0) print'(a,i9,i3)', 'ERROR in DeallocatorF2D: totalMemory = ', totalMemory, status

  return
end subroutine 

integer(8) function GetTotalMemory() 
  GetTotalMemory = totalMemory
  return
end function

end module memory_allocator_mod
!-------------------------------------------------------------------------------------------

module fileio

  use mpi_mod
  use utils, only : l2g, get_boxparameters
  use base
  use atoms
  use memory_allocator_mod

contains

!----------------------------------------------------------------------------------------
subroutine OUTPUT(atype, pos, v, q, fileNameBase)
!----------------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: pos(:,:),v(:,:)
character(len=:),allocatable,intent(in) :: fileNameBase

if(isBinary) then
  call WriteBIN(atype,pos,v,q,fileNameBase)
endif

if(isBondFile) call WriteBND(fileNameBase)
if(isPDB) call WritePDB(fileNameBase)
if(isXYZ) call WriteXYZ(fileNameBase)

return

Contains 

!--------------------------------------------------------------------------
subroutine WriteBND(fileNameBase)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileNameBase

integer :: i, ity, j, j1, jty, m
real(8) :: bndordr(MAXNEIGHBS)
integer :: igd,jgd,bndlist(0:MAXNEIGHBS)

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer :: BNDLineSize, baseCharPerAtom
integer,parameter :: MaxBNDLineSize=4096
character(MaxBNDLineSize) :: BNDOneLine
real(8),parameter :: BNDcutoff=0.3d0

character(len=:),allocatable :: BNDAllLines

integer :: scanbuf

integer :: ti,tj,tk
call system_clock(ti,tk)

! precompute the total # of neighbors
m=0
do i=1, NATOMS
   do j1 = 1, nbrlist(i,0)
!--- don't count if BO is less than BNDcutoff.
       if(BO(0,i,j1) > BNDcutoff) then 
           m=m+1
       endif
   enddo
enddo

200 format(i12.12,1x,3f12.3,1x,2i3,20(1x,i12.12,f6.3)) 

! get local datasize based on above format and the total # of neighbors
baseCharPerAtom=12+1+3*12+1+2*3 +1 ! last 1 for newline
localDataSize=NATOMS*(baseCharPerAtom)+m*(1+12+6)

if( (baseCharPerAtom+MAXNEIGHBS*(1+12+6)) > MaxBNDLineSize) then
    print'(a,i6,2i12)', 'ERROR: MaxBNDLineSize is too small @ WriteBND', &
                    myid, baseCharPerAtom+MAXNEIGHBS*(1+12+6), MaxBNDLineSize
endif

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".bnd", &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(character(len=localDataSize) :: BNDAllLines)
BNDALLLines=""

BNDLineSize=0
do i=1, NATOMS
   ity = nint(atype(i))
!--- get global ID for i-atom
   igd = l2g(atype(i))

!--- count the number bonds to be shown.
   bndlist(0)=0
   do j1 = 1, nbrlist(i,0)
      j = nbrlist(i,j1)
      jty = nint(atype(j))

!--- get global ID for j-atom
      jgd = l2g(atype(j))

!--- if bond order is less than 0.3, ignore the bond.
      if( BO(0,i,j1) < 0.3d0) cycle

      bndlist(0) = bndlist(0) + 1
      bndlist(bndlist(0)) = jgd
      bndordr(bndlist(0)) = BO(0,i,j1)
   enddo

   BNDOneLine=""
   write(BNDOneLine,200) igd, pos(i,1:3),nint(atype(i)),bndlist(0), &
         (bndlist(j1),bndordr(j1),j1=1,bndlist(0))

   ! remove space and add new_line
   BNDOneLine=trim(adjustl(BNDOneLine))//NEW_LINE('A')
   BNDLineSize=BNDLineSize+len(trim(BNDOneLine))

   BNDAllLines=trim(BNDAllLines)//trim(BNDOneLine)
enddo

if(localDataSize>0) then
   call MPI_File_Write(fh,BNDAllLines,localDataSize, &
        MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(BNDAllLines) 

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(20)=it_timer(20)+(tj-ti)

return
end subroutine

!--------------------------------------------------------------------------
subroutine WritePDB(fileNameBase)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileNameBase

integer :: i, ity, igd 
real(8) :: tt=0.d0, ss=0.d0

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer,parameter :: PDBLineSize=67
character(PDBLineSize) :: PDBOneLine

character(len=:),allocatable :: PDBAllLines

integer :: scanbuf

integer :: ti,tj,tk
call system_clock(ti,tk)

! get local datasize
localDataSize=NATOMS*PDBLineSize

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".pdb", &
     MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

allocate(character(len=localDataSize) :: PDBAllLines)
PDBAllLines=""

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

do i=1, NATOMS

  ity = nint(atype(i))
!--- calculate atomic temperature 
  tt = hmas(ity)*sum(v(i,1:3)*v(i,1:3))
  tt = tt*UTEMP*1d-2 !scale down to use two decimals in PDB format 

!--- sum up diagonal atomic stress components 
  ss = sum(astr(1:3))/3.d0*USTRS

  tt = q(i)

  igd = l2g(atype(i))
  write(PDBOneLine,100)'ATOM  ',0, atmname(ity), igd, pos(i,1:3), tt, ss

  PDBOneLine(PDBLineSize:PDBLineSize)=NEW_LINE('A')
  PDBAllLines=trim(PDBAllLines)//trim(PDBOneLine)

enddo

if(localDataSize>0) then
    call MPI_File_Write(fh,PDBAllLines,localDataSize, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(PDBAllLines)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

100 format(A6,I5,1x,A2,i12,4x,3f8.3,f6.2,f6.2)

call system_clock(tj,tk)
it_timer(21)=it_timer(21)+(tj-ti)


end subroutine

!--------------------------------------------------------------------------
subroutine WriteXYZ(fileNameBase)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileNameBase

integer :: i, ity, idx1, idx0 , igd

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer :: OneLineSize, MetaDataSize
character(60) :: a60
character(len=:),allocatable :: OneLine,AllLines

integer :: scanbuf

integer :: ti,tj,tk
call system_clock(ti,tk)

MetaDataSize = 9 + 60 + 2
write(a60,'(3f12.5,3f8.3)')  lata,latb,latc,lalpha,lbeta,lgamma

if(isPQEq) then
  ! name + pos(i,1:3) + q(i) + gid + spos(i,1:3) + newline
  OneLineSize = 3 + 60 + 20 + 9 + 60 + 1 
else
  ! name + pos(i,1:3) + q(i) + gid + newline
  OneLineSize = 3 + 36 + 8 + 9 + 1 
endif

! get local datasize
localDataSize=NATOMS*OneLineSize
if(myid==0) localDataSize = localDataSize + MetaDataSize

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".xyz", &
     MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(character(OneLineSize) :: OneLine)
allocate(character(localDataSize) :: AllLines)

idx0=1 ! allline index
if(myid==0) then
  write(AllLines(idx0:idx0+8),'(i9)') GNATOMS; idx0=idx0+9
  write(AllLines(idx0:idx0),'(a1)') new_line('A'); idx0=idx0+1
  write(AllLines(idx0:idx0+59),'(a60)') a60; idx0=idx0+60
  write(AllLines(idx0:idx0),'(a1)') new_line('A'); idx0=idx0+1
endif

do i=1, NATOMS
  idx1 = 1 ! oneline index

  ity = nint(atype(i))
  igd = l2g(atype(i))

! element name
  write(OneLine(idx1:idx1+2),'(a3)') atmname(ity); idx1=idx1+3

! atom positions. with high precision for dielectric calculation
  if(isPQEq) then
     write(OneLine(idx1:idx1+59),'(3es20.12)') pos(i,1:3); idx1=idx1+60
     write(OneLine(idx1:idx1+19),'(es20.12)') q(i); idx1=idx1+20
  else
     write(OneLine(idx1:idx1+35),'(3f12.5)') pos(i,1:3); idx1=idx1+36
     write(OneLine(idx1:idx1+7),'(3f8.3)') q(i); idx1=idx1+8
  endif

! global Id
  write(OneLine(idx1:idx1+8),'(i9)') igd; idx1=idx1+9

! shell charge positions with high precision for dielectric calculation
  if(isPQEq) then
     write(OneLine(idx1:idx1+59),'(3es20.12)') spos(i,1:3); idx1=idx1+60
  endif
  write(OneLine(idx1:idx1),'(a1)') new_line('A'); idx1=idx1+1

  write(AllLines(idx0:idx0+OneLineSize-1),'(a)') OneLine; idx0=idx0+OneLineSize

enddo

if(localDataSize>0) then
    call MPI_File_Write(fh,AllLines,localDataSize, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(AllLines, OneLine)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(21)=it_timer(21)+(tj-ti)


end subroutine

end subroutine OUTPUT

!--------------------------------------------------------------------------
subroutine ReadXYZ(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:),intent(inout) :: atype,q
real(8),allocatable,dimension(:,:),intent(inout) :: rreal,v,f

real(8) :: dbuf6(6), mat(3,3), mati(3,3)

integer :: ti,tj,tk

integer :: i,j,i1, ntot, iigd, num_unit, funit
real(8),allocatable :: pos0(:)
integer,allocatable :: atype0(:)
character(2) :: c2

call system_clock(ti,tk)

if(myid==0) then
  open(newunit=funit, file=trim(filename), status='old', form='formatted')

  read(funit,fmt=*) num_unit
  allocate(atype0(num_unit), pos0(3*num_unit))

  read(funit,fmt=*) dbuf6(1:6)

  do i = 1, num_unit
     read(funit,fmt=*) c2,pos0(i*3-2:i*3) 

!TODO: use atmname to get the element-to-integer mapping.
!      atmname has to be set based on the input files, e.g. ffield or ffn.in.
     do i1=1, size(atmname)
       if(c2(1:2) == atmname(i1)) exit
     enddo
     atype0(i)=i1
  enddo

  close(funit)
endif

call MPI_Bcast(num_unit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myid>0) allocate(pos0(3*num_unit), atype0(num_unit))

call MPI_Bcast(dbuf6,size(dbuf6),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
lata=dbuf6(1); latb=dbuf6(2); latc=dbuf6(3)
lalpha=dbuf6(4); lbeta=dbuf6(5); lgamma=dbuf6(6)

call MPI_Bcast(atype0,size(atype0),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(pos0,size(pos0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!--- allocate arrays
if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0

iigd = num_unit*myid ! for global ID
do i=1,num_unit
   pos(i,1:3) = pos0(i*3-2:i*3)
   atype(i) = dble(atype0(i)) + (iigd+i)*1d-13
enddo 

!--- assuming XYZ format is in the real coordinates given for one domain. 
!    To run multiple domain job, normalize the coords, update H-matrix, 
!    update the box origin, then scale back to the real coords. 
call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
call matinv(mat,mati)
call xu2xs_inplace(mati,[0.d0,0.d0,0.d0],num_unit,pos)

NATOMS=num_unit

!--- update to glocal cell parameters
lata=lata*vprocs(1)
latb=latb*vprocs(2)
latc=latc*vprocs(3)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call xs2xu_inplace(hh,obox,num_unit,pos)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end

!--------------------------------------------------------------------------
subroutine ReadBIN(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:) :: atype,q
real(8),allocatable,dimension(:,:) :: rreal,v,f

integer :: i,i1

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: idata(:)
real(8),allocatable :: dbuf(:)
real(8) :: ddata(6), d10(10)

real(8) :: rnorm(NBUFFER,3), mat(3,3)
integer :: j

integer :: ti,tj,tk
call system_clock(ti,tk)


! Meta Data: 
!  Total Number of MPI ranks and MPI ranks in xyz (4 integers)
!  Number of resident atoms per each MPI rank (nprocs integers) 
!  current step (1 integer) + lattice parameters (6 doubles)

nmeta=4+nprocs+1
allocate(idata(nmeta))
metaDataSize = 4*nmeta + 8*6

if(myid==0) write(6,'(2a)') 'INFO: Opening file in ReadBIN(): ', trim(fileName)
call MPI_File_Open(MPI_COMM_WORLD,trim(fileName),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

if(ierr > 0) then
   if(myid==0) write(6,'(2a)') 'MPI_File_Open Error in ReadBIN(): ', trim(fileName)
   call MPI_Finalize(ierr)
   stop
endif

! read metadata at the beginning of file
offsettmp=0
call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
call MPI_File_Read(fh,idata,nmeta,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

offsettmp=4*nmeta
call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
call MPI_File_Read(fh,ddata,6,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

NATOMS = idata(4+myid+1)
current_step = idata(nmeta)
deallocate(idata)
lata=ddata(1); latb=ddata(2); latc=ddata(3)
lalpha=ddata(4); lbeta=ddata(5); lgamma=ddata(6)

! Get local datasize: 10 doubles for each atoms
localDataSize = 8*NATOMS*10

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! Since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset = scanbuf + metaDataSize

! nprocs-1 rank has the total data size
fileSize = offset
!call MPI_Bcast(fileSize,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
!call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize
call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(dbuf(10*NATOMS))
call MPI_File_Read(fh,dbuf,10*NATOMS,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0

do i=1, NATOMS
    i1=10*(i-1)
    rnorm(i,1:3)=dbuf(i1+1:i1+3)
    v(i,1:3)=dbuf(i1+4:i1+6)
    q(i)=dbuf(i1+7)
    atype(i)=dbuf(i1+8)
    qsfp(i)=dbuf(i1+9)
    qsfv(i)=dbuf(i1+10)
enddo
deallocate(dbuf)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo

! Check: box params are updated first time here but cutoff distance is determined 
! based on ffield params. give a flag? 
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call xs2xu(hh,obox,NATOMS,rnorm,rreal)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end

!--------------------------------------------------------------------------
subroutine WriteBIN(atype, rreal, v, q, fileNameBase)
!--------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: rreal(:,:),v(:,:)
character(len=:),allocatable,intent(in) :: fileNameBase

integer :: i,j

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: ldata(:),gdata(:)
real(8) :: ddata(6)
real(8),allocatable :: dbuf(:)

real(8) :: rnorm(NBUFFER,3)

integer :: ti,tj,tk
call system_clock(ti,tk)

call xu2xs(hhi,obox,NATOMS,rreal,rnorm)

if(.not. isBinary) return

! Meta Data: 
!  Total Number of MPI ranks and MPI ranks in xyz (4 integers)
!  Number of resident atoms per each MPI rank (nprocs integers) 
!  current step (1 integer) + lattice parameters (6 doubles)
nmeta=4+nprocs+1
metaDataSize = 4*nmeta + 8*6

! Get local datasize: 10 doubles for each atoms
localDataSize = 8*NATOMS*10

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".bin",MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! Since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset = scanbuf + metaDataSize

! save metadata at the beginning of file
allocate(ldata(nmeta),gdata(nmeta))
ldata(:)=0
ldata(4+myid+1)=NATOMS
call MPI_ALLREDUCE(MPI_IN_PLACE,ldata,nmeta,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
gdata = ldata
gdata(1)=nprocs
gdata(2:4)=vprocs
gdata(nmeta)=nstep+current_step

ddata(1)=lata; ddata(2)=latb; ddata(3)=latc
ddata(4)=lalpha; ddata(5)=lbeta; ddata(6)=lgamma

if(myid==0) then
   offsettmp=0
   call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
   call MPI_File_Write(fh,gdata,nmeta,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

   offsettmp=4*nmeta
   call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
   call MPI_File_Write(fh,ddata,6,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
endif
deallocate(ldata,gdata)

! set offset at the beginning of the local write
offset=offset-localDataSize
call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(dbuf(10*NATOMS))
do i=1, NATOMS
   j = (i - 1)*10
   dbuf(j+1:j+3)=rnorm(i,1:3)
   dbuf(j+4:j+6)=v(i,1:3)
   dbuf(j+7)=q(i)
   dbuf(j+8)=atype(i)
   dbuf(j+9)=qsfp(i)
   dbuf(j+10)=qsfv(i)
enddo
call MPI_File_Write(fh,dbuf,10*NATOMS,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
deallocate(dbuf)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(23)=it_timer(23)+(tj-ti)

return
end

end module

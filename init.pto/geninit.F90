module params
implicit none
integer :: vprocs(3)=(/1,1,1/)
integer :: mc(3)=(/1,1,1/)

integer :: nprocs, mctot
integer,allocatable :: lnatoms(:), lnatoms1(:), lnatoms2(:)
real(8) :: L1,L2,L3,Lalpha,Lbeta,Lgamma
real(8) :: lbox(3),obox(3), dtype
real(8) :: H(3,3), Hi(3,3), rr(3),rr1(3),vv(3),qq,rmin(3), rmax(3)
real(8) :: qfsp0=0.d0, qfsv0=0.d0
integer :: natoms,  ntot, ix,iy,iz
real(8),allocatable :: pos0(:,:), pos1(:,:)
character(3),allocatable :: ctype0(:), ctype1(:)
integer,allocatable :: itype0(:)
real(8),allocatable :: itype1(:)

character(256) :: inputFileName="input.xyz"
character(256) :: ffieldFileName="../ffield"
character(256) :: outputDirName="."

character(256) :: fnote

character(len=3),allocatable :: atomNames(:)
integer :: numParams, numAtomNames

logical :: getReal=.false., getNorm=.false., noCoordinateShift=.false.

logical :: isLG, is_fnn, is_reaxff

contains

!-------------------------------------------------------------------------------------------
integer function getstr(linein,lineout)
implicit none
!-------------------------------------------------------------------------------------------

character(len=:),allocatable,intent(in out) :: linein,lineout
integer :: pos1

!--- remove whitespace 
linein = adjustl(linein)

!--- return if black line
if(len(linein)==0) then
  getstr=-2
  return
endif

!--- return if it's a comment line or entirely whitespace
if(linein(1:1)=='#' .or. linein == repeat(' ', len(linein)) ) then
   getstr=-1
   return
endif

! find position in linein to get a token. if whitespace is not found, take entire line
pos1=index(linein,' ')-1
if(pos1==-1) pos1=len(linein)

lineout=linein(1:pos1)
linein=linein(pos1+1:)
getstr=len(lineout)

end function

!----------------------------------------------------------------
subroutine convertAndDumpCoordinate(L1,L2,L3,lalpha,lbeta,lgamma,mc,natoms,ctype,pos,note)
implicit none
!----------------------------------------------------------------
real(8),intent(in) :: L1,L2,L3,lalpha,lbeta,lgamma
integer,intent(in) :: mc(3)
integer,intent(in) :: natoms
real(8),intent(inout) :: pos(3,natoms)
character(3),intent(in) :: ctype(natoms)
character(256),intent(in) :: note

integer :: i,j,k,n
real(8) :: hh(3,3), rr(3)
character(8) :: outFile

call getBox(L1,L2,L3,lalpha,lbeta,lgamma)

write(6,'(a)') '------------------------------------------------------------'

if(getReal) then
  print'(a)','converting coordinates : normalized to real'
  outFile="real.xyz"
else if(getNorm) then
  print'(a)','converting coordinates : real to normalized'
  outFile="norm.xyz"
else
  print'(a)','no coordinate conversion detected '
  return
endif
open(40,file=outFile,form="formatted")

!--- write header part
write(40,'(i12,3x,a)') natoms*mc(1)*mc(2)*mc(3),'"'//trim(adjustl(note))//'"'
write(40,'(3f10.5, 3f10.3)') L1*mc(1),L2*mc(2),L3*mc(3),lalpha,lbeta,lgamma
write(6,'(a)') "header after conversion :"
write(6,'(i12,3x,a)') natoms*mc(1)*mc(2)*mc(3),'"'//trim(adjustl(note))//'"'
write(6,'(3f10.5, 3f10.3)') L1*mc(1),L2*mc(2),L3*mc(3),lalpha,lbeta,lgamma

write(6,'(a)') '------------------------------------------------------------'

if(getReal) then !--- from normalized to real coordinates

  call getBox(L1*mc(1),L2*mc(2),L3*mc(3),lalpha,lbeta,lgamma)
  do i=0,mc(1)-1
  do j=0,mc(2)-1
  do k=0,mc(3)-1
     do n=1,natoms

        rr(1:3)=pos(1:3,n)+(/i,j,k/)
        rr(1:3)=rr(1:3)/mc(1:3)

        rr(1)=sum(h(1,1:3)*rr(1:3))
        rr(2)=sum(h(2,1:3)*rr(1:3))
        rr(3)=sum(h(3,1:3)*rr(1:3))

        write(40,'(a3,1x,3f15.9)') ctype(n),rr(1:3)
     enddo
  enddo; enddo; enddo

else if(getNorm) then !--- from real to normalized coordinates

  do i=0,mc(1)-1
  do j=0,mc(2)-1
  do k=0,mc(3)-1
     do n=1,natoms

        rr(1)=sum(hi(1,1:3)*pos(1:3,n))
        rr(2)=sum(hi(2,1:3)*pos(1:3,n))
        rr(3)=sum(hi(3,1:3)*pos(1:3,n))

        rr(1:3)=rr(1:3)+(/i,j,k/)
        rr(1:3)=rr(1:3)/mc(1:3)

        write(40,'(a3,1x,3f15.9)') ctype(n),rr(1:3)
     enddo
  enddo; enddo; enddo

else
endif

close(40)
print'(2a)','coordinates are saved in ',outFile

stop

end subroutine
!----------------------------------------------------------------
subroutine get_atomnames_for_fnn(fileName)
implicit none
!----------------------------------------------------------------
character(len=256),intent(in) :: fileName

character(len=256) :: linein0
character(len=:),allocatable :: linein, token
character(len=3) :: c3
integer :: funit, i

open(newunit=funit,file=trim(fileName),status='old',form='formatted')

do while (.true.)
  read(funit,'(a)',end=10) linein0
  linein = trim(adjustl(linein0))

  if(getstr(linein, token) > 0) then
   
    if(token=='model') then
      if(getstr(linein, token)>0) then
        c3 = token
        if(.not.allocated(atomNames)) then
           allocate(character(len=3)::atomNames(1))
           atomNames(1) = c3
        else
           atomNames = [atomNames,c3]
        endif
      else
        print*,'ERROR: while processing', trim(linein0)
      endif
    endif
  endif

end do

10 close(funit)

numAtomNames = size(atomNames)
do i=1, numAtomNames
   print'(i3,a,a2 $)',i,'-',atomNames(i)
enddo
print*

end subroutine

!----------------------------------------------------------------
subroutine getAtomNames(fileName)
implicit none
!----------------------------------------------------------------
integer :: i
character(256) :: fileName

write(6,'(a)') '------------------------------------------------------------'

write(6,'(2a)') 'element name and IDs from ', trim(fileName)

open(20,file=fileName)
read(20,*)
read(20,*) numParams
do i=1, numParams
   read(20,*)
enddo
read(20,*) numAtomNames
read(20,*)
read(20,*)
read(20,*)
allocate( character(3) :: atomNames(numAtomNames) )
do i=1, numAtomNames
   read(20,*) atomNames(i)
   read(20,*)
   read(20,*)
   read(20,*)
   if(isLG) read(20,*)
   print'(i3,a,a2 $)',i,'-',atomNames(i)
enddo
close(20)

write(6,*)
write(6,'(a)') '------------------------------------------------------------'

end subroutine

!----------------------------------------------------------------
subroutine getBox(la,lb,lc,angle1,angle2,angle3)
!----------------------------------------------------------------
implicit none
real(8),intent(in) :: la,lb,lc, angle1,angle2,angle3
real(8) :: hh1, hh2 , lal, lbe, lga
real(8) :: pi=datan(1.d0)*4.d0

!--- convet unit for angles
lal=angle1*pi/180.d0
lbe=angle2*pi/180.d0
lga=angle3*pi/180.d0

!--- construct H-matrix
hh1=lc*(cos(lal)-cos(lbe)*cos(lga))/sin(lga)
hh2=lc*sqrt( 1.d0-cos(lal)**2-cos(lbe)**2-cos(lga)**2 + &
             2*cos(lal)*cos(lbe)*cos(lga) )/sin(lga)

H(1,1)=la
H(2,1)=0.d0
H(3,1)=0.d0
H(1,2)=lb*cos(lga)
H(2,2)=lb*sin(lga)
H(3,2)=0.d0
H(1,3)=lc*cos(lbe)
H(2,3)=hh1
H(3,3)=hh2

call matinv(H,Hi)

return 
end subroutine


!----------------------------------------------------------------
subroutine matinv(m1,m2)
! get inverse of m1 and save to m2
!----------------------------------------------------------------
implicit none
real(8) :: m1(3,3), m2(3,3), detm

m2(1,1) = m1(2,2)*m1(3,3)-m1(2,3)*m1(3,2)
m2(1,2) = m1(1,3)*m1(3,2)-m1(1,2)*m1(3,3)
m2(1,3) = m1(1,2)*m1(2,3)-m1(1,3)*m1(2,2)
m2(2,1) = m1(2,3)*m1(3,1)-m1(2,1)*m1(3,3)
m2(2,2) = m1(1,1)*m1(3,3)-m1(1,3)*m1(3,1)
m2(2,3) = m1(1,3)*m1(2,1)-m1(1,1)*m1(2,3)
m2(3,1) = m1(2,1)*m1(3,2)-m1(2,2)*m1(3,1)
m2(3,2) = m1(1,2)*m1(3,1)-m1(1,1)*m1(3,2)
m2(3,3) = m1(1,1)*m1(2,2)-m1(1,2)*m1(2,1)

detm = m1(1,1)*m1(2,2)*m1(3,3) + m1(1,2)*m1(2,3)*m1(3,1) &
     + m1(1,3)*m1(2,1)*m1(3,2) - m1(1,3)*m1(2,2)*m1(3,1) &
     - m1(1,2)*m1(2,1)*m1(3,3) - m1(1,1)*m1(2,3)*m1(3,2)

m2(:,:) = m2(:,:)/detm

return
end subroutine

end module


!================================================================
program geninit
use params
implicit none
!================================================================
integer :: i, j, k, n, ia, myid, sID
integer(8) :: ii=0
character(6) :: a6
character(64) :: argv

!--- get input parameters
do i=1, command_argument_count()
   call get_command_argument(i,argv)
   select case(adjustl(argv))
     case("-help","-h")
       print'(a)', "./geninit -mc 1 1 1 -vprocs 1 1 1 -inputxyz input.xyz -ffield ffield [-r or -n]"
       stop
     case("-mc","-m")
       call get_command_argument(i+1,argv); read(argv,*) mc(1)
       call get_command_argument(i+2,argv); read(argv,*) mc(2)
       call get_command_argument(i+3,argv); read(argv,*) mc(3)
     case("-vprocs","-v")
       call get_command_argument(i+1,argv); read(argv,*) vprocs(1)
       call get_command_argument(i+2,argv); read(argv,*) vprocs(2)
       call get_command_argument(i+3,argv); read(argv,*) vprocs(3)
     case("-inputxyz","-i")
       call get_command_argument(i+1,argv)
       inputFileName=adjustl(argv)
     case("-ffield","-f")
       call get_command_argument(i+1,argv)
       ffieldFileName=adjustl(argv)
     case("-output","-o")
       call get_command_argument(i+1,argv)
       outputDirName=adjustl(argv)
     case("-fnn")
       call get_command_argument(i+1,argv)
       ffieldFileName=adjustl(argv)
       is_fnn = .true.
     case("-lg")
       isLG = .true. 
     case("-getreal","-r")
       getReal=.true.
     case("-getnorm","-n")
       getNorm=.true.
     case("-noCoordinateShift","-nocs")
       noCoordinateShift=.true.
       print*,'INFO: enabling noCoordinateShift flag'
     case default
   end select
enddo

inquire(file='../fnn.in', exist=is_fnn)
if(is_fnn) then
   print'(a)',repeat('-',60)
   ffieldFileName='../fnn.in'
   print'(a)','Found fnn.in file: '//trim(adjustl(ffieldFileName))
   print'(a)','Generate rxff.bin for FNN simulation' 
   print'(a)',repeat('-',60)
endif

inquire(file='../ffield', exist=is_reaxff)
if(is_reaxff) then
   print'(a)',repeat('-',60)
   ffieldFileName='../ffield'
   print'(a)','Found ffield file: '//trim(adjustl(ffieldFileName))
   print'(a)','Generate rxff.bin for ReaxFF simulation' 
   print'(a)',repeat('-',60)
endif


mctot=mc(1)*mc(2)*mc(3)
nprocs=vprocs(1)*vprocs(2)*vprocs(3)
allocate(lnatoms(0:nprocs-1),lnatoms1(0:nprocs-1), lnatoms2(0:nprocs-1))

open(1,file=inputFileName,form="formatted")
write(6,'(a60)') repeat('-',60)
write(6,'(a20,a)') ' input file: ', trim(inputFileName)
write(6,'(a20,a)') ' ffield file: ', trim(ffieldFileName)
write(6,'(a20,a)') ' output dir: ', trim(outputDirName)
write(6,'(a20,i9,3i6)') ' nprocs,vprocs: ',nprocs, vprocs(1:3)
write(6,'(a20,i9,3i6)') ' mctot,mc: ',mctot, mc(1:3)
write(6,'(a60)') repeat('-',60)

if(is_fnn) then
  call get_atomnames_for_fnn(ffieldFileName)
else
  call getAtomNames(ffieldFileName)
endif

!--- read # of atoms and file description
read(1,*) natoms, fnote
!--- read lattice parameters
read(1,*) L1, L2, L3, Lalpha, Lbeta, Lgamma
print'(a60)',repeat('-',60)
print'(a)', '** UNIT CELL PROFIEL **'
print'(a)', trim(fnote)
print'(a,3x,6f12.3)','1, L2, L3, Lalpha, Lbeta, Lgamma: ', & 
           L1, L2, L3, Lalpha, Lbeta, Lgamma
print'(a60)',repeat('-',60)

!--- allocate arrays for atom position and type
allocate(ctype0(natoms),pos0(3,natoms),itype0(natoms))
allocate(ctype1(natoms*mctot),pos1(3,natoms*mctot),itype1(natoms*mctot))
do i=1, natoms
   read(1,*) ctype0(i),pos0(1:3,i)
   ctype0(i)=adjustl(ctype0(i))

   do j=1, numAtomNames
      if(ctype0(i)==atomNames(j)) then
         !print*, ctype0(i), atomNames(j)
         itype0(i)=j
         exit
      endif
   enddo
enddo
close(1)

if(getNorm .or. getReal) & 
  call convertAndDumpCoordinate(L1,L2,L3,lalpha,lbeta,lgamma, &
                                mc,natoms,ctype0,pos0,fnote)

!--- repeat the unit cell
ntot=0
do ix=0, mc(1)-1
do iy=0, mc(2)-1
do iz=0, mc(3)-1
do i=1, natoms
   ntot=ntot+1
   rr(1:3)=pos0(1:3,i) 
   rr(1)=(rr(1)+ix)/mc(1)
   rr(2)=(rr(2)+iy)/mc(2)
   rr(3)=(rr(3)+iz)/mc(3)
   pos1(1:3,ntot)=rr(1:3)
   ctype1(ntot)=ctype0(i)
   itype1(ntot)=itype0(i)+ntot*1d-13+1d-14
   !print'(a3,1x,i3,3f8.3,3f)',ctype0(i),itype0(i),pos0(1:3,i), rr(1:3)
enddo 
enddo; enddo; enddo

!--- shift coordinates, then shift a bit to avoid zero coordinates
if (noCoordinateShift) then
   do i=1, 3
      rmin(i)=minval(pos1(i,:))
      pos1(i,:)=pos1(i,:)-rmin(i)
   enddo
endif

!--- wrap back atom coordinates if necessary
do i=1, ntot
   do j=1, 3
      pos1(j,i)=mod(pos1(j,i),1.d0)
   enddo
   pos1(1:3,i)=pos1(1:3,i)+1d-9 ! shift by a small value to avoid coordinates at 0
enddo

!--- check coordinates
do i=1, 3
   rmin(i)=minval(pos1(i,:))
   rmax(i)=maxval(pos1(i,:))
enddo
print'(a,3es15.5)','rmin(1:3): ', rmax(1:3)
print'(a,3es15.5)','rmax(1:3): ', rmin(1:3)

!--- update natoms & H-matrix 
natoms=natoms*mctot
L1=L1*mc(1); L2=L2*mc(2); L3=L3*mc(3)
call getBox(L1,L2,L3,lalpha,lbeta,lgamma)

!--- count how many atoms per MPI domain, get lnatoms()
lnatoms(:)=0
do n=1,natoms 
   i=pos1(1,n)*vprocs(1)
   j=pos1(2,n)*vprocs(2)
   k=pos1(3,n)*vprocs(3)
   sID = i + j*vprocs(1) + k*vprocs(1)*vprocs(2)
   lnatoms(sID)=lnatoms(sID)+1
enddo

!--- prefix sum
lnatoms1(0)=0
do i=1, nprocs-1
   lnatoms1(i)=lnatoms1(i-1)+lnatoms(i-1)
enddo

!--- To avoid opening too many files, sort atom data using disk
lbox(1:3)=1.d0/vprocs(1:3)
lnatoms2(:)=0 ! for atom counting & result check
open(1,file="all.bin",form="unformatted",access="stream")
do n=1,natoms 
   i=pos1(1,n)*vprocs(1)
   j=pos1(2,n)*vprocs(2)
   k=pos1(3,n)*vprocs(3)
   sID = i + j*vprocs(1) + k*vprocs(1)*vprocs(2)
   obox(1:3)=lbox(1:3)*(/i,j,k/)
   pos1(1:3,n)=pos1(1:3,n)-obox(1:3)

!--- total number of atoms before n-th atom, 
!--- each atom has 32 (=3*8 + 8) bytes data
   ii=lnatoms1(sID)+lnatoms2(sID) 
   write(1,pos=ii*32+1) pos1(1:3,n),itype1(n)

   lnatoms2(sID)=lnatoms2(sID)+1
enddo
close(1)

!--- error check
do i=0, nprocs-1
   if(lnatoms(i)/=lnatoms2(i)) then
      print'(a,3i9)', 'Error,i,lnatoms(i),lnatoms2(i)',i,lnatoms(i),lnatoms2(i)
      stop
   endif
enddo

qq=0.d0; vv(:)=0.d0; qfsp0=0.d0; qfsv0=0.d0
open(1,file="all.bin",form="unformatted",access="stream")
open(20,file="geninit.xyz") 
open(30,file=trim(adjustl(outputDirName))//"/rxff.bin",form="unformatted",access="stream")
write(30) nprocs, vprocs(1:3)
write(30) lnatoms(0:nprocs-1)
write(30) 0
write(30) L1, L2, L3, Lalpha, Lbeta, Lgamma

write(20,'(i12)') sum(lnatoms(:))
write(20,'(6f12.3,3x,a)') L1, L2, L3, Lalpha, Lbeta, Lgamma, trim(inputFileName)
do myid=0,nprocs-1
   write(a6(1:6),'(i6.6)') myid

   i=mod(myid,vprocs(1))
   j=mod(myid/vprocs(1),vprocs(2))
   k=myid/(vprocs(1)*vprocs(2))
   obox(1:3)=lbox(1:3)*(/i,j,k/)

   do n=1,lnatoms(myid)
      read(1) rr(1:3),dtype

      write(30)rr(1:3)
      write(30)vv(1:3)
      write(30)qq
      write(30)dtype
      write(30)qfsp0
      write(30)qfsv0

      rr1(1:3)=rr(1:3)+obox(1:3)
      rr(1)=sum(H(1,1:3)*rr1(1:3))
      rr(2)=sum(H(2,1:3)*rr1(1:3))
      rr(3)=sum(H(3,1:3)*rr1(1:3))

      write(20,'(a3, $)') atomNames(nint(dtype))
      write(20,'(3f12.5)') rr(1:3)
   enddo

enddo
close(1, status='delete')
close(20)
close(30)

print'(a60)',repeat('-',60)
print'(a)', '** SUPER CELL PROFIEL **'
print'(a,3i6)', 'number of domains: ', vprocs(1:3)
print'(a,i12)', 'total number of atoms: ', sum(lnatoms)
print*,'atoms per domain: ',lnatoms(:)
print'(a,3x,6f12.3)','L1, L2, L3, Lalpha, Lbeta, Lgamma: ', & 
           L1, L2, L3, Lalpha, Lbeta, Lgamma
print'(a60)',repeat('-',60)
print'(a)','geninit successfully finished'
end


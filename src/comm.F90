!--------------------------------------------------------------------------------------------------------------
subroutine COPYATOMS(imode, dr, atype, pos, v, f, q)
use atoms
!
! TODO: update notes here
!
! In this subroutine, boundary atoms are copied to neighbour nodes. Three subroutines, 
! <store_atoms()>, <send_recv()> and  <append_atoms()> work together, shareing variables 
! <sbuffer>, <rbuffer>, <ne>, <na>, <ns> and <nr>. 
!
!--- Variables --- 
!<dflag>:  Direction FLAG specifing the direction of communication.
!          1 +x, 2 -x, 3 +y, 4 -y, 5 +z, 6 -z. 
!<parity>:  PARITY of a node in <dflag> direction. 
!<nlayer>:  Nr of LAYERs to be copied. COPYATOMS() features two different modes, migration mode and copy mode.
!           <nlayer> will be used as a flag to distinguish them. 
!
!--- Shared Variables ---
!<ne>: # of elements one atom has.
!   Example) in copy mode, three positions[xyz] + charge + atomtype = 5 elements. 
!<na>, <ns> & <nr>: Nr of All of transfered elemetns, Number of elements to be Sent, and Recieved one respectivly.
!        
!--- Subroutines ---
!<store_atoms()>: store boundary atoms information.
!<send_recv()):   send & receive. 
!<append_atoms()>: append received infromation into array.  deallocate <sbuffer> & <rbuffer>
!
!--------------------------------------------------------------------------------------------------------------
implicit none

type pack2darray
   real(8),pointer :: ptr(:,:)
   logical :: shift=.false.
   character(8) :: name=''
end type 

type pack1darray
   real(8),pointer :: ptr(:)
   logical :: cpbk=.false.
   character(8) :: name=''
end type 

integer,intent(IN) :: imode 
real(8),intent(IN) :: dr(3)
real(8),target :: atype(NBUFFER), q(NBUFFER)
real(8),target :: pos(NBUFFER,3),v(NBUFFER,3),f(NBUFFER,3)

real(8),target :: nipos(NBUFFER,3)  

logical :: commflag(NBUFFER)
type(pack1darray),allocatable :: pack1d(:)
type(pack2darray),allocatable :: pack2d(:)
type(pack2darray),allocatable :: norm2d(:)

integer :: ixyz,tn1,tn2, dflag
integer :: ni, ity

integer :: ti,tj,tk,tti,ttj

integer,parameter :: dinv(6)=(/2,1,4,3,6,5/)
integer,parameter :: cptridx(6)=(/0,0,2,2,4,4/)
integer,parameter :: is_xyz(6)=(/1,1,2,2,3,3/)  !<- [xxyyzz]

call system_clock(tti,tk)

call initialize(imode)

do dflag=1, 6

   tn1 = target_node(dflag)
   tn2 = target_node(dinv(dflag))
   ixyz = is_xyz(dflag)
   
   if(imode==MODE_CPBK) then  ! communicate with neighbors in reversed order
      tn1 = target_node(7-dinv(dflag)) ! <-[563412] 
      tn2 = target_node(7-dflag) ! <-[654321] 
      ixyz = is_xyz(7-dflag)        ! <-[332211]
   endif

   call step_preparation(dflag, dr, commflag)

   call store_atoms(tn1, dflag, imode)
   call send_recv(tn1, tn2, myparity(ixyz))
   call append_atoms(dflag, imode)

enddo

call finalize(imode)

!write(*,*)  "MPI send stats: ",myid,ns_atoms(1:6)
!write(*,*)  "MPI receive stats: ",myid,nr_atoms(1:6)

!--- for array size stat
if(mod(nstep,pstep)==0) then
  ni=nstep/pstep+1
  if(imode==MODE_MOVE) maxas(ni,4)=na/ne
  if(imode==MODE_COPY) maxas(ni,5)=na/ne
endif

call system_clock(ttj,tk)
it_timer(4)=it_timer(4)+(ttj-tti)

return

CONTAINS 
!--------------------------------------------------------------------------------------------------------------
subroutine initialize(imode)
implicit none
!--------------------------------------------------------------------------------------------------------------
integer,intent(in) :: imode
integer :: a, np2d, np1d, nn2d

!--- clear total # of copied atoms, sent atoms, recieved atoms
na=0;ns=0;nr=0

!--- Since cached atoms are stored after the resident atoms (i.e. i > NATOMS), 
!--- initialize the cache atom pointer 0th elem with NATOMS.
copyptr(0)=NATOMS

!--- set the number of data per atom 
select case(imode)
   case(MODE_COPY)

      np2d=1; np1d=7; nn2d=1
      if(isPQEq) np2d=np2d+1 ! for spos

      allocate(pack2d(np2d),pack1d(np1d),norm2d(nn2d))
      ne = size(pack2d)*3+size(pack1d)

      a=1
      pack2d(a)%ptr=>pos; pack2d(a)%shift=.true.; a=a+1
      if(isPQEq) then
         pack2d(a)%ptr=>spos; pack2d(a)%shift=.false.; a=a+1
      endif

      a=1
      pack1d(a)%ptr=>atype; a=a+1
      pack1d(a)%ptr=>q;     a=a+1
      pack1d(a)%ptr=>qs;    a=a+1
      pack1d(a)%ptr=>qt;    a=a+1
      pack1d(a)%ptr=>hs;    a=a+1
      pack1d(a)%ptr=>ht;    a=a+1
      pack1d(a)%ptr=>frcindx; pack1d(a)%cpbk=.true.; a=a+1

      a=1
      norm2d(a)%ptr=>pos

      do a=1, NATOMS
         frcindx(a)=a
      enddo

   case(MODE_MOVE)

      np2d=2; np1d=6; nn2d=1

      if(isPQEq) np2d=np2d+1 ! for spos
      if(isSpring) then
         np2d=np2d+1 ! for nipos
         nn2d=nn2d+1 ! for nipos
      endif

      allocate(pack2d(np2d),pack1d(np1d),norm2d(nn2d))
      ne = size(pack2d)*3+size(pack1d)

      a=1
      pack2d(a)%ptr=>pos; pack2d(a)%shift=.true.; a=a+1
      pack2d(a)%ptr=>v; a=a+1 
      if(isPQEq) then
         pack2d(a)%ptr=>spos; pack2d(a)%shift=.false.; a=a+1
      endif
      if(isSpring) then
         pack2d(a)%ptr=>nipos; pack2d(a)%shift=.true.; a=a+1
      endif

      a=1
      pack1d(a)%ptr=>atype; a=a+1
      pack1d(a)%ptr=>q;     a=a+1
      pack1d(a)%ptr=>qs;    a=a+1
      pack1d(a)%ptr=>qt;    a=a+1
      pack1d(a)%ptr=>qsfp;  a=a+1
      pack1d(a)%ptr=>qsfv;  a=a+1

      a=1
      norm2d(a)%ptr=>pos;  a=a+1

      if(isSpring) then
         norm2d(a)%ptr=>nipos; a=a+1
      endif

   case(MODE_QCOPY1)

      np2d=0; np1d=2; nn2d=1
      allocate(pack1d(np1d),norm2d(nn2d))
      ne = size(pack1d)

      a=1
      pack1d(a)%ptr=>qs; a=a+1
      pack1d(a)%ptr=>qt; a=a+1

      a=1
      norm2d(a)%ptr=>pos

   case(MODE_QCOPY2)

      np2d=0; np1d=3; nn2d=1
      allocate(pack1d(np1d),norm2d(nn2d))
      ne = size(pack1d)

      a=1
      pack1d(a)%ptr=>hs; a=a+1
      pack1d(a)%ptr=>ht; a=a+1
      pack1d(a)%ptr=>q; a=a+1

      a=1
      norm2d(a)%ptr=>pos

   case(MODE_CPBK)
      ne = NE_CPBK

   case default
      print'(a,i3)', "ERROR: imode doesn't match in COPYATOMS: ", imode

end select

!--- Get the normalized local coordinate will be used through this function. 
if(allocated(norm2d)) then
   do a=1,size(norm2d)
      call xu2xs_inplace(max(copyptr(6),NATOMS),norm2d(a)%ptr)
   enddo
endif

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine finalize(imode)
implicit none
!--------------------------------------------------------------------------------------------------------------
integer,intent(in) :: imode
integer :: i,a

if(imode==MODE_MOVE) then
!--- remove atoms which are transfered to neighbor nodes.
!--- if atype is smaller than zero (this is done in store_atoms), ignore the atom.
   ni=0
   do i=1, copyptr(6)
      if(nint(atype(i))>0) then
        ni=ni+1
        do a=1,size(pack2d)
           pack2d(a)%ptr(ni,1:3) = pack2d(a)%ptr(i,1:3)
        enddo
        do a=1,size(pack1d)
           pack1d(a)%ptr(ni) = pack1d(a)%ptr(i)
        enddo
      endif
   enddo 

!--- update the number of resident atoms
   NATOMS=ni

endif

!--- by here, we got new atom positions in the normalized coordinate, need to update real coordinates.
if(allocated(norm2d)) then
   do a=1,size(norm2d)
      call xs2xu_inplace(copyptr(6),norm2d(a)%ptr)
   enddo
endif

if(allocated(pack1d)) deallocate(pack1d)
if(allocated(pack2d)) deallocate(pack2d)
if(allocated(norm2d)) deallocate(norm2d)

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine step_preparation(dflag, dr, commflag)
implicit none
!--------------------------------------------------------------------------------------------------------------
integer,intent(in) :: dflag
real(8),intent(in) :: dr(3)
logical,intent(inout) :: commflag(NBUFFER)
integer :: i

!--- start buffering data depending on modes. all copy&move modes use buffer size, dr, to select atoms.
do i=1, copyptr(cptridx(dflag))
   commflag(i) = inBuffer(dflag,dr,pos(i,is_xyz(dflag)))
enddo

return
end subroutine


!--------------------------------------------------------------------------------------------------------------
subroutine send_recv(tn1, tn2, mypar)
use atoms
! shared variables::  <ns>, <nr>, <na>, <sbuffer()>, <rbuffer()>
! This subroutine only takes care of communication part. won't be affected by wether atom migration or atom 
! copy mode. 
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) ::tn1, tn2, mypar 
integer :: recv_stat(MPI_STATUS_SIZE)
real(8) :: recv_size

!--- if myid is the same of target-node ID, don't use MPI call.
!--- Just copy <sbuffer> to <rbuffer>. Because <send_recv()> will not be used,
!--- <nr> has to be updated here for <append_atoms()>.
if(myid==tn1) then
   if(ns>0) then
      nr=ns
      call CheckSizeThenReallocate(rbuffer,nr)
      rbuffer(1:ns) = sbuffer(1:ns)
   else
      nr=0
   endif

   return
endif

call system_clock(ti,tk)

if (mypar == 0) then

     ! the number of elements per data packet has to be greater than 1, for example NE_COPY = 10.
     ! if ns == 0, send one double to tell remote rank that there will be no atom data to be sent. 
     if (ns > 0) then 
       call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn1, 10, MPI_COMM_WORLD, ierr)
     else
       call MPI_SEND(1, 1, MPI_DOUBLE_PRECISION, tn1, 10, MPI_COMM_WORLD, ierr)
     endif

     call MPI_Probe(tn2, 11, MPI_COMM_WORLD, recv_stat, ierr)
     call MPI_Get_count(recv_stat, MPI_DOUBLE_PRECISION, nr, ierr)

     call CheckSizeThenReallocate(rbuffer,nr)

     call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn2, 11, MPI_COMM_WORLD, recv_stat, ierr)

     ! the number of elements per data packet has to be greater than 1, for example NE_COPY = 10.
     ! nr == 1 means no atom data to be received. 
     if(nr==1) nr=0 

elseif (mypar == 1) then

     call MPI_Probe(tn2, 10, MPI_COMM_WORLD, recv_stat, ierr)
     call MPI_Get_count(recv_stat, MPI_DOUBLE_PRECISION, nr, ierr)

     call CheckSizeThenReallocate(rbuffer,nr)

     call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn2, 10, MPI_COMM_WORLD, recv_stat, ierr)

     ! the number of elements per data packet has to be greater than 1, for example NE_COPY = 10.
     ! nr == 1 means no atom data to be received. 
     if(nr==1) nr=0 

     if (ns > 0) then
        call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn1, 11, MPI_COMM_WORLD, ierr)
     else
        call MPI_SEND(1, 1, MPI_DOUBLE_PRECISION, tn1, 11, MPI_COMM_WORLD, ierr)
     endif

endif

call system_clock(tj,tk)
it_timer(25)=it_timer(25)+(tj-ti)

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine store_atoms(tn, dflag, imode)
use atoms
! <nlayer> will be used as a flag to change the behavior of this subroutine. 
!    <nlayer>==0 migration mode
!            > 0 copy mode
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: tn, dflag, imode 

integer :: n,ni,is,a,b,ioffset
real(8) :: sft 

call system_clock(ti,tk)

!--- reset the number of atoms to be sent
ns=0

if(imode==MODE_CPBK) then

   is = 7 - dflag !<- [654321] reversed order direction flag

   n = copyptr(is) - copyptr(is-1) + 1
   call CheckSizeThenReallocate(sbuffer,n*ne)

   do n=copyptr(is-1)+1, copyptr(is)
      sbuffer(ns+1) = dble(frcindx(n))
      sbuffer(ns+2:ns+4) = f(n,1:3)
      ns=ns+ne
   enddo

else

!--- # of elements to be sent. should be more than enough. 
   ni = copyptr(cptridx(dflag))*ne

!--- <sbuffer> will deallocated in store_atoms.
   call CheckSizeThenReallocate(sbuffer,ni)

!--- get the coordinate Index to be Shifted.
   is = int((dflag-1)/2) !<- [012] means [xyz]

!--- When atom moves to neighbor nodes, their coordinates must be shifted. 
!--- xshift() returns the edge length of one node assuming all of node size is same.
   sft=xshift(dflag)

!--- start buffering data depending on modes. all copy&move modes use buffer size, dr, to select atoms.
   do n=1, copyptr(cptridx(dflag))

      if(commflag(n)) then

        ioffset=ns

        if(allocated(pack2d)) then
           do a=1,size(pack2d)
              sbuffer(ioffset+1:ioffset+3)=pack2d(a)%ptr(n,1:3)
              if(pack2d(a)%shift) sbuffer(ioffset+1+is) = sbuffer(ioffset+1+is) + sft
              ioffset=ioffset+3
           enddo
        endif

        if(allocated(pack1d)) then
           do a=1,size(pack1d)
              if(pack1d(a)%cpbk) then
                 sbuffer(ioffset+1)=n
              else
                 sbuffer(ioffset+1)=pack1d(a)%ptr(n)
              endif
              ioffset=ioffset+1
           enddo
        endif

!--- In append_atoms subroutine, atoms with <atype>==-1 will be removed
        if(imode==MODE_MOVE) atype(n) = -1.d0 

!--- increment the number of atoms to be sent 
        ns = ns + ne
     endif

   enddo

endif

! store the number of atoms being transfered in dflag direction
ns_atoms(dflag)= ns/ne

call system_clock(tj,tk)
it_timer(26)=it_timer(26)+(tj-ti)

end subroutine store_atoms

!--------------------------------------------------------------------------------------------------------------
subroutine append_atoms(dflag, imode)
use atoms
! <append_atoms> append copied information into arrays
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, imode 
integer :: m, i, ine, a, ioffset

call system_clock(ti,tk)

if( (na+nr)/ne > NBUFFER) then
   print'(a,i4,5i8)', "ERROR: over capacity in append_atoms; myid,na,nr,ne,(na+nr)/ne,NBUFFER: ", &
        myid,na,nr,ne,(na+nr)/ne, NBUFFER
   call MPI_FINALIZE(ierr)
   stop
endif

if(imode == MODE_CPBK) then  

   do i=0, nr/ne-1
!--- get current index <ine> in <rbuffer(1:nr)>.
!--- Append the transferred forces into the original position of force array.
      ine=i*ne
      m = nint(rbuffer(ine+1))
      f(m,1:3) = f(m,1:3) + rbuffer(ine+2:ine+4)
   enddo

else
   !--- store the last transfered atom index to <dflag> direction.
   !--- if no data have been received, use the index of previous direction.
   if(nr==0) then 
      copyptr(dflag) = copyptr(dflag-1)
   else
      copyptr(dflag) = copyptr(dflag-1) + nr/ne
   endif

! number of atomd received from the neighbor
nr_atoms(dflag)= nr/ne

!--- go over the buffered atom
   do i=0, nr/ne-1

!--- get current index <ine> in <rbuffer(1:nr)>.
      ine=i*ne

!--- current atom index; resident + 1 + stored atoms so far + atoms in buffer
      m = copyptr(dflag-1) + 1 + i

      ioffset = ine

      if(allocated(pack2d)) then
         do a=1,size(pack2d)
            pack2d(a)%ptr(m,1:3)=rbuffer(ioffset+1:ioffset+3)
            ioffset=ioffset+3
         enddo
      endif

      if(allocated(pack1d)) then
         do a=1,size(pack1d)
            pack1d(a)%ptr(m)=rbuffer(ioffset+1)
            ioffset=ioffset+1
         enddo
      endif
      
   enddo

endif   

!--- update the total # of transfered elements.
na=na+nr

call system_clock(tj,tk)
it_timer(27)=it_timer(27)+(tj-ti)

end subroutine append_atoms

!--------------------------------------------------------------------------------------------------------------
real(8) function xshift(dflag)
use atoms
! This subroutine returns a correction value in coordinate for transfered atoms.
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag 
integer :: i, j
  
  i = mod(dflag,2)  !<- i=1 means positive, i=0 means negative direction 

  if(i==1) then
      xshift =-lbox(is_xyz(dflag))
  else if(i==0) then
      xshift = lbox(is_xyz(dflag))
  endif
   
!print'(a,i,3f10.3)',"shift: ",myid, xshift
end function

!--------------------------------------------------------------------------------------------------------------
function inBuffer(dflag, dr, rr) result(isInside)
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag
real(8),intent(IN) :: dr(3), rr
logical :: isInside

select case(dflag)
   case(1) 
      isInside = lbox(1) - dr(1) < rr
   case(2) 
      isInside = rr <= dr(1)
   case(3) 
      isInside = lbox(2) - dr(2) < rr
   case(4) 
      isInside = rr <= dr(2)
   case(5) 
      isInside = lbox(3) - dr(3) < rr
   case(6) 
      isInside = rr <= dr(3)
   case default
      write(6,*) "ERROR: no matching directional flag in isInside: ", dflag
end select

end function

!--------------------------------------------------------------------------------------------------------------
subroutine CheckSizeThenReallocate(buffer,nsize)
implicit none
!--------------------------------------------------------------------------------------------------------------
real(8),allocatable :: buffer(:)
integer,intent(IN) :: nsize
integer :: ast

if(allocated(buffer)) then
   if(nsize > size(buffer)) then
       deallocate(buffer)
       allocate(buffer(2*nsize), stat=ast)
   endif
else
   allocate(buffer(nsize), stat=ast)
endif

end subroutine

end subroutine COPYATOMS

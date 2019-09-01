module lists_mod

  use mpi_mod
  implicit none

contains

!----------------------------------------------------------------------------------------
subroutine LINKEDLIST(atype, rreal, cellDims, headAtom, atomList, NatomPerCell)
use utils, only : xu2xs
use base, only : hhi, obox, copyptr, it_timer
! partitions the volume into linked-list cells <lcsize>
!----------------------------------------------------------------------------------------
real(8),allocatable,intent(in) :: atype(:), rreal(:,:)
real(8),intent(in) ::  cellDims(3)

integer,allocatable,intent(in out) :: atomList(:)
integer,allocatable,intent(in out) :: NatomPerCell(:,:,:)
integer,allocatable,intent(in out) :: headAtom(:,:,:)

real(8),allocatable :: rnorm(:,:)
integer :: n, l(3), j

integer :: ti,tj,tk
call system_clock(ti,tk)

!--- allocate rnorm with same shape of rreal
if(.not.allocated(rnorm)) allocate(rnorm, mold=rreal)
if(.not.allocated(atomList)) allocate(atomList(size(atype)))

call xu2xs(hhi,obox,copyptr(6),rreal,rnorm)

headAtom(:,:,:) = -1; atomList(:) = 0; NatomPerCell(:,:,:)=0

!--- copyptr(6) stores the last atom index copied in COPYATOMS.
do n=1, copyptr(6) 

   if(nint(atype(n))==0) cycle

   l(1:3) = floor(rnorm(n,1:3)/cellDims(1:3))

   atomList(n) = headAtom(l(1), l(2), l(3))
   headAtom(l(1), l(2), l(3)) = n
   NatomPerCell(l(1), l(2), l(3)) = NatomPerCell(l(1), l(2), l(3)) + 1
enddo

call system_clock(tj,tk)
it_timer(3)=it_timer(3)+(tj-ti)

end subroutine 

!----------------------------------------------------------------------
subroutine neighborlist(nlayer, atype, pos, pair_types, skip_check)
use base, only : natoms, myid, maxneighbs, nbrlist, header, llist, nacell, &
                 cc, rc2, nstep, pstep, copyptr, it_timer
use atoms, only : nbrindx, maxas
!----------------------------------------------------------------------
! calculate neighbor list for atoms witin cc(1:3, -nlayer:nlayer) cells.
integer,intent(in) :: nlayer
real(8),allocatable,intent(in) :: atype(:), pos(:,:)
integer,allocatable,intent(in) :: pair_types(:,:)
logical,intent(in),optional :: skip_check

integer :: c1,c2,c3, ic(3), c4, c5, c6, ierr
integer :: n, n1, m, m1, nty, mty, inxn
real(8) :: dr(3), dr2

integer :: i,j,i1,j1
logical :: isFound

integer :: ti,tj,tk
call system_clock(ti,tk)

nbrlist(:,0) = 0

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,inxn,dr,dr2) 
DO c1=-nlayer, cc(1)-1+nlayer
DO c2=-nlayer, cc(2)-1+nlayer
DO c3=-nlayer, cc(3)-1+nlayer

  m = header(c1, c2, c3)
  do m1=1, nacell(c1, c2, c3)
     mty = nint(atype(m))

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1, c2, c3] + [c4, c5, c6]

        n = header(ic(1),ic(2),ic(3))
        do n1=1, nacell(ic(1), ic(2), ic(3))

           if(n/=m) then
             nty = nint(atype(n))
             inxn = pair_types(mty, nty)

             dr(1:3) = pos(n,1:3) - pos(m,1:3) 
             dr2 = sum(dr(1:3)*dr(1:3))

             if(dr2<rc2(inxn)) then 
                nbrlist(m, 0) = nbrlist(m, 0) + 1
                nbrlist(m, nbrlist(m, 0)) = n
             endif 
           endif

           n=llist(n) 
        enddo
     enddo; enddo; enddo

     m = llist(m)
  enddo
enddo; enddo; enddo
!$omp end parallel do 

!--- to get the reverse information (i.e. from i,j1&j to i1), store <i1> into <nbrindx>.

if(.not. present(skip_check)) then
!$omp parallel do default(shared) private(i,i1,j,j1,isFound)
   do i=1, copyptr(6)
      do i1 = 1, nbrlist(i,0)
         j = nbrlist(i,i1)
         isFound=.false.
         do j1 = 1, nbrlist(j,0)
            if(i == nbrlist(j,j1)) then
               nbrindx(i,i1)=j1
               isFound=.true.
            endif
         enddo
         if(.not.isFound) &
         print'(a,i6,30i4)','ERROR: inconsistency between nbrlist and nbrindx found', &
              myid, i,nbrlist(i,0:nbrlist(i,0)), j, nbrlist(j,0:nbrlist(j,0))
      enddo
   enddo
!$omp end parallel do
   
   !--- error trap
   n=maxval(nbrlist(1:NATOMS,0))
   if(n > MAXNEIGHBS) then
      write(6,'(a45,2i5)') "ERROR: overflow of max # in neighbor list, ", myid, n
      call MPI_FINALIZE(ierr)
      stop
   endif

endif

!--- for array size stat
if(mod(nstep,pstep)==0) then
  maxas(nstep/pstep+1,2)=maxval(nbrlist(1:NATOMS,0))
endif

call system_clock(tj,tk)
it_timer(5)=it_timer(5)+(tj-ti)

end subroutine

!----------------------------------------------------------------------
subroutine get_nonbonding_list(pos, rcut)
use base, only : it_timer
use atoms, only : nbnmesh, nbheader, nbmesh, nbplist, &
                  nbllist, nbnacell, nbcc
!----------------------------------------------------------------------
real(8),allocatable,intent(in) :: pos(:,:)
real(8),intent(in) :: rcut

integer :: c1,c2,c3,c4,c5,c6,i,j,m,n,mn,iid,jid
real(8) :: dr(3), dr2, rcut2

integer :: ti,tj,tk
call system_clock(ti,tk)

rcut2 = rcut*rcut

! reset non-bonding pair list
nbplist(0,:)=0

!$omp parallel do default(shared),private(c1,c2,c3,c4,c5,c6,i,j,m,n,mn,iid,jid,dr,dr2)
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)

      do mn = 1, nbnmesh
         c4 = c1 + nbmesh(1,mn)
         c5 = c2 + nbmesh(2,mn)
         c6 = c3 + nbmesh(3,mn)

         j = nbheader(c4,c5,c6)
         do n=1, nbnacell(c4,c5,c6)

            !if(i<j .or. NATOMS<j) then
            if(i/=j) then
               dr(1:3) = pos(i,1:3) - pos(j,1:3)
               dr2 = sum(dr(1:3)*dr(1:3))

               if(dr2<=rcut2) then
                 nbplist(0,i)=nbplist(0,i)+1
                 nbplist(nbplist(0,i),i)=j
               endif

            endif

            j=nbllist(j)
         enddo
       enddo

      i=nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(15)=it_timer(15)+(tj-ti)

end subroutine

!----------------------------------------------------------------
subroutine get_mesh_for_nonbonding_list(rcut)
! setup 10[A] radius mesh to avoid visiting unecessary cells 
!----------------------------------------------------------------
use base, only : lata, latb, latc, vprocs, nbuffer
use atoms, only : maxlayers_nb, nblcsize, &
                  nbheader, nbllist, nbmesh, nbnacell, nbnmesh, nbcc
use memory_allocator_mod
implicit none

real(8),intent(in) :: rcut

real(8) :: maxrcell
real(8) :: latticePerNode(3), rr(3), dr2
integer :: i,j,k, imesh(3), maximesh, ii(3), i1


!--- initial estimate of LL cell dims
nblcsize(1:3)=3d0

!--- get mesh resolution which is close to the initial value of rlc.
latticePerNode(1)=lata/vprocs(1)
latticePerNode(2)=latb/vprocs(2)
latticePerNode(3)=latc/vprocs(3)
nbcc(1:3)=int(latticePerNode(1:3)/nblcsize(1:3))
nblcsize(1:3)=latticePerNode(1:3)/nbcc(1:3)
maxrcell = maxval(nblcsize(1:3))

!--- get # of linked list cell to cover up the non-bonding cutoff length
imesh(1:3)  = int(rcut/nblcsize(1:3)) + 1
maximesh = maxval(imesh(1:3))

!--- List up only cell indices within the cutoff range.
!--- pre-compute nmesh to get exact array size.
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= rcut**2) nbnmesh = nbnmesh + 1
enddo; enddo; enddo

call allocator(nbmesh,1,3,1,nbnmesh)

nbmesh(:,:)=0
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= rcut**2) then
      nbnmesh = nbnmesh + 1
      nbmesh(1:3,nbnmesh) = (/i, j, k/)
   endif
enddo; enddo; enddo

call allocator(nbllist,1,NBUFFER)
call allocator(nbheader, &
                -MAXLAYERS_NB,nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(3)-1+MAXLAYERS_NB)
call allocator(nbnacell, &
                -MAXLAYERS_NB,nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(3)-1+MAXLAYERS_NB)

!--- normalize nblcsize, like lcsize.
nblcsize(1:3)=nblcsize(1:3)/(/lata,latb,latc/)

end subroutine

end module

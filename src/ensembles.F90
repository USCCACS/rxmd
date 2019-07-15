module ensembles

use base
use atoms
use parameters, only : mass

implicit none

integer,parameter :: num_chain_particle = 8 
integer,parameter :: num_chain_barostat = 8 
real(8),parameter :: chain_mass_particle = 1.d0
real(8),parameter :: chain_mass_barostat = 5.d0
real(8),parameter :: barostat_mass = 200d0

!--- external pressure
real(8),parameter :: Pext = 1.d0*USTRS ! [Gpa]
real(8) :: dt2, dt4, dt8

logical :: is_nvt=.false., is_npt=.false.

type nhc_type
   integer :: num_chains
   real(8),allocatable :: force(:)
   real(8),allocatable :: momentum(:)
   real(8),allocatable :: pos(:)
   real(8),allocatable :: mass(:)
end type

type barostat_type
   real(8) :: force(3,3), momentum(3,3), pos(3,3), velocity(3,3), mass
   real(8) :: eigv(3,3), eige(3)
end type

real(8),parameter :: identity_matrix(3,3) = &
                     [[1d0,0d0,0d0], [0d0,1d0,0d0], [0d0,0d0,1d0]]

real(8) :: KbT
type(nhc_type) :: nhc_particle
type(nhc_type) :: nhc_barostat

type(barostat_type) :: barostat

contains

!-----------------------------------------------------------------------
subroutine nhcp_L_g2(deltaT)
! calculate the barostat force and update its momentum. 
! this routine also update associated eigen value&vectors. 
!-----------------------------------------------------------------------
   real(8),intent(in) :: deltaT
   integer :: ia,ib,num_rotation
   real(8) :: volume, kene

   do ia=1,3; do ib=1,3
      barostat%pos(ia,ib)=HH(ia,ib,0)
   enddo; enddo

   !if(myid==0) then
   !   print'(a,3f8.3)','barostat%pos(:,1)',barostat%pos(:,1)
   !   print'(a,3f8.3)','barostat%pos(:,2)',barostat%pos(:,2)
   !   print'(a,3f8.3)','barostat%pos(:,3)',barostat%pos(:,3)
   !endif

   volume = get_determinant(barostat%pos)
   kene = get_kinetic_energy()*2.d0
   !if(myid==0) print'(a,3es15.5)','volume,kene: ', volume, kene
  
   barostat%force = volume*(get_pint() - Pext*identity_matrix) + kene/(3*GNATOMS)*identity_matrix

   !if(myid==0) then
   !   print'(a)',repeat('-',60)
   !   print'(a,3f15.5)','barostat%force(:,1)',barostat%force(:,1)
   !   print'(a,3f15.5)','barostat%force(:,2)',barostat%force(:,2)
   !   print'(a,3f15.5)','barostat%force(:,3)',barostat%force(:,3)
   !endif

   barostat%momentum = barostat%momentum + deltaT*barostat%force
   !if(myid==0) then
   !   print'(a)',repeat('-',60)
   !   print'(a,3f15.5)','barostat%momentum(:,1)',barostat%momentum(:,1)
   !   print'(a,3f15.5)','barostat%momentum(:,2)',barostat%momentum(:,2)
   !   print'(a,3f15.5)','barostat%momentum(:,3)',barostat%momentum(:,3)
   !endif

   barostat%velocity = barostat%momentum / barostat%mass
   call jacobi(barostat%velocity,3,3,barostat%eige,barostat%eigv,num_rotation)
   !if(myid==0) then
   !   print'(a,i6,3f10.5)','num_rotation,barostat%eige(1:3)',num_rotation,barostat%eige(:)
   !   print'(a,3es15.5)','barostat%eigv(:,1)',barostat%eigv(:,1)
   !   print'(a,3es15.5)','barostat%eigv(:,2)',barostat%eigv(:,2)
   !   print'(a,3es15.5)','barostat%eigv(:,3)',barostat%eigv(:,3)
   !endif

end subroutine

!-----------------------------------------------------------------------
subroutine nhcp_L_2(deltaT)
! update the particle velocities using the eigen value&vectors of
! the barostat momentum. 
!-----------------------------------------------------------------------
   real(8),intent(in) :: deltaT
   integer :: i,ity,ia
   real(8) :: arg, vg_trace, deltaTh
   real(8) :: dmat1(3,3),dmat2(3,3), odo1(3,3), odo2(3,3)

   deltaTh = 0.5d0*deltaT

   vg_trace = get_trace(barostat%momentum)/barostat%mass

   dmat1=0.d0; dmat2=0.d0
   do ia=1,3
      arg = barostat%eige(ia) + vg_trace/(3*GNATOMS)
      dmat1(ia,ia) = exp(deltaT*arg)
      dmat2(ia,ia) = exp(deltaTh*arg)*get_sinhx(deltaTh*arg)
   enddo

   odo1 = matmul(barostat%eigv,matmul(dmat1,transpose(barostat%eigv)))
   odo2 = matmul(barostat%eigv,matmul(dmat2,transpose(barostat%eigv)))

   do i = 1, NATOMS
      ity = nint(atype(i))
      do ia = 1,3 
         v(i,ia) = sum(odo1(ia,1:3)*v(i,1:3)) + deltaTh*sum(odo2(ia,1:3)*f(i,1:3))/mass(ity)
      enddo
   enddo

end subroutine

!-----------------------------------------------------------------------
subroutine nhcp_L_1(deltaT)
! update the particle positions using the eigen value&vectors of
! the barostat momentum. 
!-----------------------------------------------------------------------
   real(8),intent(in) :: deltaT
   integer :: i,ity,ia
   real(8) :: deltaTh, dmat1(3,3), dmat2(3,3), odo1(3,3), odo2(3,3)

   deltaTh = 0.5d0*deltaT

   dmat1=0.d0; dmat2=0.d0
   do ia=1,3
      dmat1(ia,ia) = exp(deltaT*barostat%eige(ia))
      dmat2(ia,ia) = exp(deltaTh*barostat%eige(ia))*get_sinhx(deltaTh*barostat%eige(ia))
   enddo

   odo1 = matmul(barostat%eigv,matmul(dmat1,transpose(barostat%eigv)))
   odo2 = matmul(barostat%eigv,matmul(dmat2,transpose(barostat%eigv)))
   do i = 1, NATOMS
      do ia = 1,3 
         pos(i,ia) = sum(odo1(ia,1:3)*pos(i,1:3)) + deltaT*sum(odo2(ia,1:3)*v(i,1:3))
      enddo
   enddo

end subroutine

!-----------------------------------------------------------------------
subroutine nhcp_L_g1(deltaT)
! update the box size using the eigen value&vectors of the barostat momentum. 
!-----------------------------------------------------------------------
   real(8),intent(in) :: deltaT
   integer :: i,ity,ia,ib
   real(8) :: deltaTh, dmat1(3,3), odo1(3,3)
   integer :: ic
   real(8) :: odo2(3,3), vec(3)

   deltaTh = 0.5d0*deltaT

   !dmat1=0.d0
   !do ia=1,3
   !   dmat1(ia,ia) = exp(deltaT*barostat%eige(ia))
   !enddo
   !odo1 = matmul(barostat%eigv,matmul(dmat1,transpose(barostat%eigv)))
   !barostat%pos = matmul(transpose(odo1), barostat%pos)

   barostat%pos = barostat%pos + deltaT*barostat%velocity

   !if(myid==0) then
   !   print'(a)',repeat('-',60)
   !   print'(a,3(3f15.5,1x))','barostat%pos(:,1),v,f: ',barostat%pos(:,1), barostat%velocity(:,1), barostat%force(:,1)
   !   print'(a,3(3f15.5,1x))','barostat%pos(:,2),v,f: ',barostat%pos(:,2), barostat%velocity(:,2), barostat%force(:,2)
   !   print'(a,3(3f15.5,1x))','barostat%pos(:,3),v,f: ',barostat%pos(:,3), barostat%velocity(:,3), barostat%force(:,3)
   !endif

   do ia=1,3; do ib=1,3
      HH(ia,ib,0) = barostat%pos(ia,ib)
   enddo; enddo

   call xu2xs_inplace(NATOMS, pos)
   call UpdateBoxParams()
   call xs2xu_inplace(NATOMS, pos)

end subroutine


!-----------------------------------------------------------------------
function nhc_ctor(num,ctime) result(c)
!-----------------------------------------------------------------------
   integer,intent(in) :: num
   real(8),intent(in) :: ctime ! characteristic time
   type(nhc_type) :: c

   c%num_chains = num
   allocate(c%force(c%num_chains),c%pos(c%num_chains),c%mass(c%num_chains))
   allocate(c%momentum(c%num_chains+1)) ! +1 for convenience, see nhc_update_chains()

   c%pos=0.d0; c%momentum=0.d0; c%force=0.d0

   c%mass=KbT*ctime**2

   return
end function

!-----------------------------------------------------------------------
subroutine update_thermostats() 
!-----------------------------------------------------------------------

if(is_npt) then

   call nhc_update_chains(nhc_particle)
   call nhc_update_chains(nhc_barostat, barostat)

   !call print_nhc_chains(nhc_barostat) 

else if(is_nvt) then

   call nhc_update_chains(nhc_particle)

endif

end subroutine

!-----------------------------------------------------------------------
subroutine initialize_ensemble() 
!-----------------------------------------------------------------------
integer :: ia,ib

!is_nvt=.true.
is_npt=.true.

KbT = treq
if(is_nvt .or. is_npt) then
   if(myid==0) print'(a,2l,f8.3)','is_nvt,is_npt,KbT: ',is_nvt,is_npt,KbT
endif

!--- Nose-Hoover chain
if(is_npt) then

   nhc_particle = nhc_ctor(num=num_chain_particle, ctime=100d0*dt)
   nhc_barostat = nhc_ctor(num=num_chain_barostat, ctime=100d0*dt)

   barostat%force=0d0;  barostat%momentum=0d0; barostat%velocity=0d0
   barostat%mass = KbT*(100d0*dt)**2*(3*GNATOMS+3)
   do ia=1,3
   do ib=1,3
      barostat%pos(ia,ib) = HH(ia,ib,0)
   enddo; enddo

   if(myid==0) print'(a,es15.5)', 'barostat%mass: ', barostat%mass

else if(is_nvt) then

   nhc_particle = nhc_ctor(num=num_chain_particle, ctime=100d0*dt)

endif

dt2 = dt*0.5d0
dt4 = dt*0.25d0
dt8 = dt*0.125d0

end subroutine

!-----------------------------------------------------------------------
subroutine print_nhc_chains(nhc) 
!-----------------------------------------------------------------------
   type(nhc_type),intent(in out) :: nhc
   integer :: i, m
   
   if(myid==0) then
      do i = 1, size(nhc%mass)
         print'(a,i3,3es15.5)','i,momentum,force,mass: ', &
                    i,nhc%momentum(i),nhc%force(i),nhc%mass(i)
      enddo
   endif

end subroutine

!-----------------------------------------------------------------------
subroutine nhc_update_chains(nhc, barostat) 
!-----------------------------------------------------------------------
   type(nhc_type),intent(in out) :: nhc
   type(barostat_type),optional :: barostat

   real(8) :: kene, NKbT, scaling_factor, exp1, exp2
   real(8) :: mm(3,3)
   integer :: i, m
   
   if(present(barostat)) then
      mm = matmul(transpose(barostat%momentum),barostat%momentum)
      kene = get_trace(mm)/barostat%mass
      NKbT = 6d0*KbT !see the paragraph below Eq 46 in 2010, Chemi Phys 307, 294-305, (2010)
      !if(myid==0) then
      !  print'(a,l,3es15.5)','momentum(:,1): ',present(barostat),barostat%momentum(:,1)
      !  print'(a,l,3es15.5)','momentum(:,2): ',present(barostat),barostat%momentum(:,2)
      !  print'(a,l,3es15.5)','momentum(:,3): ',present(barostat),barostat%momentum(:,3)
      !  print'(a,l,2es15.5)','kene,NKbT: ', present(barostat),kene, NKbT
      !endif
   else
      kene = get_kinetic_energy()*2d0
      NKbT = 3d0*GNATOMS*KbT
      !if(myid==0) print'(a,l,2es15.5)','kene,NKbT: ', present(barostat),kene, NKbT
   endif

   !NOTE: momentum(num_chains+1) is always zero so that exp1=exp2=1.0
   do i = nhc%num_chains, 1, -1
      exp1 = exp(-dt8*nhc%momentum(i+1))
      exp2 = exp1*exp1

      if(i==1) then
         nhc%force(i) = kene - NKbT
      else
         nhc%force(i) = nhc%momentum(i-1)**2/nhc%mass(i-1) - KbT
      endif

      nhc%momentum(i) = nhc%momentum(i)*exp2 + dt4*nhc%force(i)*exp1
   enddo

   do i = 1, nhc%num_chains
      nhc%pos(i) = nhc%pos(i) + dt2*nhc%momentum(i)/nhc%mass(i)
   enddo

   scaling_factor = exp(-dt2*nhc%momentum(1)/nhc%mass(1))

   !if(myid==0) print'(a,l,9es15.5)','scaling_factor,momentum,force: ', &
   !present(barostat), scaling_factor, nhc%momentum(1), nhc%force(1), KbT

   if(present(barostat)) then
!if(myid==0) then
!print'(a,3f15.5)', 'barostat%momentum: ', barostat%momentum(1:3,1)
!print'(a,3f15.5)', 'barostat%momentum: ', barostat%momentum(1:3,2)
!print'(a,3f15.5)', 'barostat%momentum: ', barostat%momentum(1:3,3)
!print'(a,3f15.5)', 'scaling_factor: ', scaling_factor
!endif
      barostat%momentum = barostat%momentum*scaling_factor
   else
      v(1:NATOMS,1:3) = v(1:NATOMS,1:3)*scaling_factor
   endif
   kene = kene*scaling_factor**2

   !NOTE: momentum(num_chains+1) is always zero so that exp1=exp2=1.0
   do i = 1, nhc%num_chains
      exp1 = exp(-dt8*nhc%momentum(i+1))
      exp2 = exp1*exp1

      if(i==1) then
         nhc%force(i) = kene - NKbT
      else
         nhc%force(i) = nhc%momentum(i-1)**2/nhc%mass(i-1) - KbT
      endif

      nhc%momentum(i) = nhc%momentum(i)*exp2 + dt4*nhc%force(i)*exp1
   enddo

end subroutine

!-----------------------------------------------------------------------
function get_pint() result(pint)
!-----------------------------------------------------------------------
integer :: i,ity, ia, ib
real(8) :: pint(3,3), pint0(3,3), volume

pint0=0.d0
do i=1, NATOMS
   ity = atype(i)
   do ia=1,3
   do ib=1,3
      pint0(ia,ib) = pint0(ia,ib) + mass(ity)*v(i,ia)*v(i,ib)
      pint0(ia,ib) = pint0(ia,ib) + f(i,ia)*pos(i,ib)
   enddo; enddo
enddo


!--- symmetrize internal pressure
do ia=1,3
do ib=1,3
   pint(ia,ib) = 0.5d0*(pint0(ia,ib) + pint0(ib,ia))
enddo; enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, pint, size(pint), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

volume = get_determinant(barostat%pos)
pint = pint/volume

return
end function

!-----------------------------------------------------------------------
subroutine npt_test()
!-----------------------------------------------------------------------
integer :: i,ity, ia, ib
integer,save :: counter = 0
real(8),save :: pint1_sum(3,3) = 0.d0, pint2_sum(3,3) = 0.d0
real(8) :: pint1(3,3), pint2(3,3)
real(8) :: volume, temperature, uni

volume = &
HH(1,1,0)*(HH(2,2,0)*HH(3,3,0) - HH(3,2,0)*HH(2,3,0)) + &
HH(2,1,0)*(HH(3,2,0)*HH(1,3,0) - HH(1,2,0)*HH(3,3,0)) + &
HH(3,1,0)*(HH(1,2,0)*HH(2,3,0) - HH(2,2,0)*HH(1,3,0))

pint1=0.d0; pint2=0.d0
do i=1, NATOMS
   ity = atype(i)
   do ia=1,3
   do ib=1,3
      pint1(ia,ib) = pint1(ia,ib) + mass(ity)*v(i,ia)*v(i,ib)
      pint2(ia,ib) = pint2(ia,ib) + f(i,ia)*pos(i,ib)
   enddo; enddo
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, pint1, size(pint1), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, pint2, size(pint2), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

temperature = 0.5d0*(pint1(1,1)+pint1(2,2)+pint1(3,3))/GNATOMS*UTEMP

pint1_sum = pint1_sum + pint1/volume 
pint2_sum = pint2_sum + pint2/volume 

counter = counter + 1

uni = USTRS/counter

if(myid==0) then
   print'(a)',repeat('-',60)
   print'(a,2es15.5)','volume,temperature: ', volume, temperature
   print'(a,3es15.5,2x,3es15.5)','pint[12](:,1):  ', pint1_sum(:,1)*uni, pint2_sum(:,1)*uni
   print'(a,3es15.5,2x,3es15.5)','pint[12](:,2):  ', pint1_sum(:,2)*uni, pint2_sum(:,2)*uni
   print'(a,3es15.5,2x,3es15.5)','pint[12](:,3):  ', pint1_sum(:,3)*uni, pint2_sum(:,3)*uni
   print'(a,3es15.5,2x)','pint[1+2](:,1): ', (pint1_sum(:,1)+pint2_sum(:,1))*uni
   print'(a,3es15.5,2x)','pint[1+2](:,2): ', (pint1_sum(:,2)+pint2_sum(:,2))*uni
   print'(a,3es15.5,2x)','pint[1+2](:,3): ', (pint1_sum(:,3)+pint2_sum(:,3))*uni
   print'(a)',repeat('-',60)
endif

end subroutine

!-----------------------------------------------------------------------
function get_kinetic_energy() result(kene)
!-----------------------------------------------------------------------
   real(8) :: kene
   integer :: i, ity

   kene = 0.d0
   do i=1, NATOMS
      ity = nint(atype(i))
      kene = kene + hmas(ity)*sum(v(i,1:3)*v(i,1:3)) ! 0.5*mass*v**2
   enddo
 
   call MPI_ALLREDUCE(MPI_IN_PLACE, kene, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

   return
end function


!----------------------------------------------------------------------!
subroutine jacobi(a,n,np,d,v,nrot)
!----------------------------------------------------------------------!
!  Diagonalizes a real symmetric matrix by Jacobi transformation.
!  From "Numerical Recipes".
!----------------------------------------------------------------------!
implicit none
integer,parameter :: NMAX=500

integer,intent(in) :: n,np
integer,intent(in out) :: nrot
real(8),intent(in out) :: d(np),v(np,np)
real(8),value :: a(np,np)

integer :: i,ip,iq,j
real(8) :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

do ip=1,n
  do iq=1,n
    v(ip,iq)=0d0
  enddo 
  v(ip,ip)=1d0
enddo 

do ip=1,n
  b(ip)=a(ip,ip)
  d(ip)=b(ip)
  z(ip)=0d0
enddo

nrot=0
do i=1, 50

  sm=0.d0
  do ip=1,n-1
    do iq=ip+1,n
      sm=sm+abs(a(ip,iq))
    enddo
  enddo

  if(sm.eq.0d0) return

  if(i.lt.4)then
    tresh=0.2d0*sm/n**2
  else
    tresh=0d0
  endif

  do ip=1,n-1
    do iq=ip+1,n

      g=100d0*abs(a(ip,iq))

      if( (i.gt.4).and.(abs(d(ip))+ g.eq.abs(d(ip))) .and. &
          (abs(d(iq))+g.eq.abs(d(iq)))) then

        a(ip,iq)=0d0

      else if(abs(a(ip,iq)).gt.tresh)then

        h=d(iq)-d(ip)
        if(abs(h)+g.eq.abs(h))then
          t=a(ip,iq)/h
        else
          theta=0.5d0*h/a(ip,iq)
          t=1d0/(abs(theta)+sqrt(1d0+theta**2))
          if(theta.lt.0d0)t=-t
        endif
        c=1d0/sqrt(1d0+t**2)
        s=t*c
        tau=s/(1d0+c)
        h=t*a(ip,iq)
        z(ip)=z(ip)-h
        z(iq)=z(iq)+h
        d(ip)=d(ip)-h
        d(iq)=d(iq)+h
        a(ip,iq)=0d0

        do j=1,ip-1
          g=a(j,ip)
          h=a(j,iq)
          a(j,ip)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        enddo

        do j=ip+1,iq-1
          g=a(ip,j)
          h=a(j,iq)
          a(ip,j)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        enddo

        do j=iq+1,n
          g=a(ip,j)
          h=a(iq,j)
          a(ip,j)=g-s*(h+g*tau)
          a(iq,j)=h+s*(g-h*tau)
        enddo

        do j=1,n
          g=v(j,ip)
          h=v(j,iq)
          v(j,ip)=g-s*(h+g*tau)
          v(j,iq)=h+s*(g-h*tau)
        enddo

        nrot=nrot+1

      endif

    enddo
  enddo

  do ip=1,n
    b(ip)=b(ip)+z(ip)
    d(ip)=b(ip)
    z(ip)=0d0
  enddo

enddo

print*, 'ERROR: too many iterations in jacobi: '

return
end subroutine

!----------------------------------------------------------------------!
function get_determinant(m1) result (det)
!----------------------------------------------------------------------!
implicit none
real(8),intent(in) :: m1(3,3)
real(8) :: det

det = m1(1,1)*m1(2,2)*m1(3,3) + m1(1,2)*m1(2,3)*m1(3,1) &
    + m1(1,3)*m1(2,1)*m1(3,2) - m1(1,3)*m1(2,2)*m1(3,1) &
    - m1(1,2)*m1(2,1)*m1(3,3) - m1(1,1)*m1(2,3)*m1(3,2)

return
end function

!----------------------------------------------------------------------!
function get_trace(m1) result (trace)
!----------------------------------------------------------------------!
implicit none
real(8),intent(in) :: m1(3,3)
real(8) :: trace

trace = m1(1,1)+m1(2,2)+m1(3,3)

return
end function

!----------------------------------------------------------------------!
function get_sinhx(x) result (s)
!----------------------------------------------------------------------!
implicit none
real(8),intent(in) :: x
real(8) :: s
real(8) :: x2, x4, x6, x8, x10
real(8),parameter :: a2=1d0/6d0, a4=1d0/120d0, a6=1d0/5040d0, & 
                     a8=1d0/36288d0, a10=1d0/39916800d0

x2=x*x; x4=x2*x2; x6=x2*x4; x8=x4*x4; x10=x2*x8
s = 1.d0 + a2*x2 + a4*x4 + a6*x6 + a8*x8 + a10*x10

return
end function

end module

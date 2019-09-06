module stat_mod

  implicit none
  real(8),parameter :: pi = 4.d0*datan(1.d0)
  integer,parameter :: NTABLES = 2000
  real(8),parameter :: RCUT = 40d0, DRI = NTABLES/RCUT
  !real(8),parameter :: RCUT = 15d0, DRI = NTABLES/RCUT
  real(8),parameter :: QCUT = 10d0, DQ = QCUT/NTABLES

  type NSD_type ! Neutron Scattering Data type
     character(len=2) :: name
     real(8) :: length
  end type

  ! neutron scattering length data are from 
  !  https://www.nist.gov/ncnr/neutron-scattering-lengths-list
  type(NSD_type) :: NSD0(11)=[ NSD_type(name='Ge',length=8.185d-5), & 
                              NSD_type(name='Se',length=7.970d-5), &
                              NSD_type(name='Sb',length=5.57d-5), &
                              NSD_type(name='Te',length=5.80d-5), & 
                              NSD_type(name='C',length=6.646d-5), &
                              NSD_type(name='Si',length=4.1491d-5), &
                              NSD_type(name='O',length=5.803d-5), &
                              NSD_type(name='Al',length=3.449d-5), &
                              NSD_type(name='H',length=-3.7390d-5), &
                              NSD_type(name='Na',length=3.63d-5), &
                              NSD_type(name='Cl',length=9.5770d-5) ]

  type string_array
     character(len=:),allocatable :: str
  end type

  type mdframe
     character(len=:),allocatable :: filename
     type(string_array),allocatable :: elems(:)
     real(8),allocatable :: pos(:,:), v(:,:), f(:,:), q(:)
     real(8) :: lattice(6)
     integer :: num_atoms
  end type

  type analysis_context

     type(string_array),allocatable :: elems(:)
     type(NSD_type),allocatable :: NSD(:)

     real(8),allocatable :: concentration(:), gr(:,:,:), nr(:,:,:), ba(:,:,:,:), sq(:,:,:)

     integer :: num_atoms, num_atom_types
     real(8) :: volume, lattice(6), kvector(3)
     integer,allocatable :: num_atoms_per_type(:)

     integer :: num_sample_frames 

  contains

     procedure :: print => print_analysis_context 
     procedure :: save_stat => save_analysis_results

  end type

contains

  subroutine print_analysis_context(this) 
     class(analysis_context) :: this
     character(len=:),allocatable :: name
     integer :: i

     print'(a)',repeat('-',60)
     print'(a,i6)','num_atoms: ', this%num_atoms
     print'(a,3f8.3, 3f8.2)','lattice: ', this%lattice
     print'(a,es15.5)','volume: ', this%volume
     print'(a,3f8.3)','k-vector: ', this%kvector

     print'(a, i6)', 'num_elements: ', size(this%elems)
     print'(a $)', 'elements: '
     do i = 1, size(this%elems)
        name = this%elems(i)%str
        print'(i1, a1, a, a $)',get_index(this%elems, name),'-',adjustl(name), '   '
     enddo
     print*
     print'(a,10i6)', 'num_atoms_per_type: ', this%num_atoms_per_type
     print'(a,10f8.5)', 'concentration: ', this%concentration
     print'(a $)', 'neutron scattering: '
     do i=1, size(this%NSD)
         print'(i1,a1,a2,es10.3,a $)', i, '-', this%NSD(i)%name, this%NSD(i)%length, ', '
     enddo
     print*
     print'(a)',repeat('-',60)

  end subroutine

  subroutine save_analysis_results(this) 
     class(analysis_context) :: this
     integer :: iunit,ity,jty,k,kk
     real(8) :: dr, rho, dqk, Snq, Snq_denom, Gnr, Gnr_denom, prefactor, prefactor2

     ! get the number density, rho
     rho = this%num_atoms/this%volume

     open(newunit=iunit,file='gr.dat',form='formatted')

     write(unit=iunit,fmt='(a $)') ' distance'
     do ity=1,size(this%elems)
     do jty=1,size(this%elems)
        !print'(a5,i1,a1,i1 $)',',    ',ity,'-',jty
        write(unit=iunit,fmt='(a5,a5 $)') & 
           '     ',this%elems(ity)%str//'-'//this%elems(jty)%str
     enddo; enddo
     write(unit=iunit,fmt='(a)') '  Gnr'

     Gnr_denom = sum(this%NSD(:)%length * this%concentration(:))
     Gnr_denom = Gnr_denom**2

     do k=1,size(this%gr,dim=3)

        dr = k/DRI
        prefactor = 4.d0*pi*dr*dr*rho/DRI
        write(unit=iunit,fmt='(f12.5 $)') dr

        Gnr = 0.d0
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)

          prefactor2 = this%concentration(jty) * this%num_atoms_per_type(ity) * this%num_sample_frames
          this%gr(ity,jty,k) = this%gr(ity,jty,k)/(prefactor*prefactor2)
          write(iunit, fmt='(f12.5,3x $)') this%gr(ity,jty,k)

          Gnr = Gnr + this%gr(ity,jty,k) * &
                      this%concentration(ity) * this%concentration(jty) * & 
                      this%NSD(ity)%length * this%NSD(jty)%length
        enddo; enddo

        write(iunit, fmt='(f12.5)') Gnr/Gnr_denom
     enddo

     close(iunit)


     open(newunit=iunit,file='sq.dat',form='formatted')
     write(unit=iunit,fmt='(a $)') ' wave_number '
     do ity=1,size(this%elems)
     do jty=1,size(this%elems)
        write(unit=iunit,fmt='(a5,a5 $)') & 
           '     ',this%elems(ity)%str//'-'//this%elems(jty)%str
     enddo; enddo
     write(unit=iunit,fmt='(a)') '  Snq'

     ! get the denominator of Sn(q) 
     Snq_denom = sum(this%NSD(:)%length * this%concentration(:))
     Snq_denom = Snq_denom**2

     do kk=1, size(this%sq,dim=3)

        dqk = DQ*kk
        write(unit=iunit,fmt='(f12.5 $)') dqk

        Snq = 0.d0
        do ity=1,size(this%elems)
        do jty=1,size(this%elems)

           ! initialize sq for the r-value integral 
           this%sq(ity,jty,kk) = 0.d0

           do k=1,size(this%gr,dim=3)
              dr = k/DRI
              !prefactor = (this%gr(ity,jty,k)-1.d0)*cos(0.5d0*pi*dr/RCUT) ! (gr - 1)*taper
              prefactor = (this%gr(ity,jty,k)-1.d0) ! (gr - 1)*taper
              this%sq(ity,jty,kk) = this%sq(ity,jty,kk) + &
                  dr*dr*prefactor*sin(dr*dqk)/(dr*dqk)/DRI
           enddo

           ! get prefactor 
           prefactor = 4d0*pi*rho*sqrt(this%concentration(ity)*this%concentration(jty))
           this%sq(ity,jty,kk) = this%sq(ity,jty,kk)*prefactor

           ! add delta function if ity==jty
           if(ity==jty) this%sq(ity,jty,kk) = this%sq(ity,jty,kk) + 1.d0
           write(iunit, fmt='(f12.5,3x $)') this%sq(ity,jty,kk)

           Snq = Snq + this%sq(ity,jty,kk) * &
                       sqrt(this%concentration(ity) * this%concentration(jty)) * & 
                       this%NSD(ity)%length * this%NSD(jty)%length
        enddo; enddo
        write(iunit, fmt='(f12.5)') Snq/Snq_denom

     enddo
     close(iunit)

  end subroutine

  function get_analysis_context_from_mdframe(oneframe) result(c)
     type(mdframe),intent(in) :: oneframe
     type(analysis_context) :: c
     integer :: i, j, ne, ity, jty, kty, idx

     c%num_atoms = oneframe%num_atoms
     c%lattice = oneframe%lattice
     c%volume = oneframe%lattice(1)*oneframe%lattice(2)*oneframe%lattice(3)

     allocate(c%elems(0))
     do i=1, size(oneframe%elems)
        if( get_index(c%elems,oneframe%elems(i)%str) < 0 ) &
            c%elems = [c%elems, oneframe%elems(i)]
     enddo

     ! setup neutron scattering length data for existing atom types
     allocate(c%NSD(size(c%elems)))
     do i=1, size(c%elems)
        do j=1, size(NSD0)
           if(c%elems(i)%str==NSD0(j)%name) c%NSD(i) = NSD0(j)
        enddo
     enddo

     c%kvector(1:3)=2d0*pi/c%lattice(1:3)

     ne = size(c%elems)
     c%num_atom_types = ne

     allocate(c%gr(ne,ne,NTABLES),c%nr(ne,ne,NTABLES),c%sq(ne,ne,NTABLES),c%ba(ne,ne,ne,NTABLES))
     c%gr=0.d0;  c%nr=0.d0;  c%sq=0.d0; c%ba=0.d0

     allocate(c%concentration(ne), c%num_atoms_per_type(ne))
     c%concentration=0.0d0; c%num_atoms_per_type=0

     do i=1, size(oneframe%elems)
        idx = get_index(c%elems, oneframe%elems(i)%str) 
        if(idx>0) then 
            c%num_atoms_per_type(idx) = c%num_atoms_per_type(idx)+1
        else
            print'(a,a,i6)', 'ERROR: uncategorized atom found: in ', oneframe%filename, i
        endif
     enddo

     c%concentration = dble(c%num_atoms_per_type)/c%num_atoms

     c%num_sample_frames = 0

     call c%print()
  end function

  function get_index(elems, name) result(idx)
    type(string_array),allocatable,intent(in) :: elems(:)
    character(len=:),allocatable,intent(in) :: name
    integer :: idx

    do idx = 1, size(elems)
       if (elems(idx)%str==name) return 
    enddo
    idx = -1

    return
  end

end module

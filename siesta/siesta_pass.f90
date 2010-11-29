      module siesta_pass
      !APE modules (modified)
      use global
      use mesh, only: mesh_null, mesh_type, mesh_generation, mesh_transfer, mesh_end
      !QE modules
      use radial_grids, only: radial_grid_type,ndmx
      use ld1inc, only: zed,enne,zval,nwf,nn,ll,jj,el,rcut,octsc,enltsc,&
      phis,lloc,nbeta,rhos,pawsetup,enls,vpsloc,betas,rho,lls,isws,nlcc,&
      rhoc,rcore,nns,nwfs
      use ld1_parameters, only: nwfx,nwfsx

      type ps_io_type
      integer :: file_format
  
      type(mesh_type) :: m
  
      !Generalities
      character(3) :: symbol
      real(R8) :: z_nuc
      real(R8) :: z_val
      real(R8) :: z_ion
      integer :: wave_eq
      integer :: nspin
      integer :: scheme
  
      !Pseudo potentials
      logical :: have_psp
      integer :: n_psp!the number of pseudopotentials
      integer,  pointer :: psp_l(:)
      real(R8), pointer :: psp_j(:)
      real(R8), pointer :: psp_v(:,:)
  
      !Pseudo wavefunctions
      logical :: have_wfs
      integer :: n_wfs 
      integer, pointer :: wfs_n(:)
      integer, pointer :: wfs_l(:)
      real(R8), pointer :: wfs_j(:)
      character(len=5), pointer :: wfs_label(:)
      real(R8), pointer :: wfs_rc(:)
      real(R8), pointer :: wfs_occ(:,:)
      real(R8), pointer :: wfs_ev(:)
      real(R8), pointer :: wfs(:,:)
      real(R8), pointer :: rho_val(:)
  
      !XC
      integer :: ixc!exchange correlation functional
      !ape seems to only really allow certain combinations of X and C. set this as constant for now
  
      !Non-linear core-corrections
      logical :: nlcc
      real(R8) :: nlcc_rc
      real(R8), pointer :: rho_core(:)
  
      !Kleinman-Bylander projectors
      logical :: have_kb!
      integer :: kb_l_local!
      integer :: kb_n_proj! 
      integer,  pointer :: kb_l(:)!
      real(R8), pointer :: kb_j(:)!
      real(R8), pointer :: kb_v_local(:)!
      real(R8), pointer :: kb_e(:)!
      real(R8), pointer :: kb_proj(:,:)!
      end type ps_io_type
 
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!MAIN FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine siesta_output(grid_in)
       type(radial_grid_type), intent(in):: grid_in
       type(ps_io_type) :: info
       print *, "call the ps fill up routine"
       call ps_io_type_fill(grid_in,info)
       call ps_io_save(info)
      endsubroutine siesta_output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ps_io_type_fill(grid_in,info)
       type(radial_grid_type), intent(in) :: grid_in
       type(ps_io_type) ,intent(out) :: info
       !mesh parameters
       type(mesh_type) :: m
       integer :: no_mesh_points,i,j
       real(R8) :: rn,r1
       !generalities
       character, external :: atom_name*2

       info%file_format=M_SIESTA! these values are defined in the global.o file

       !set the mesh
       info%m%type=2!set the mesh type to LOG1 ri=b*exp(a*i)
       !set number of mesh points
       no_mesh_points=grid_in%mesh
       allocate(m%r(no_mesh_points))
       m%np=no_mesh_points
      
       !set the mesh parameters 
       r1=grid_in%r(1)
       rn=grid_in%r(no_mesh_points)
       m%a = log(rn/r1)/real(no_mesh_points - 1,R8)
       m%b = r1/exp(m%a)

       !allocate mesh points
       m%r(:)=grid_in%r(:)
       !put this back into the info type
       info%m=m

       !set the generalities
       !set nuclear charge
       info%z_nuc=zed
       !set charge of the ion
       info%z_ion=zed-enne
       !set name
       info%symbol=atom_name(nint(zed))
       !set valence charge
       info%z_val=zval
       !wave_eq to integrate one of schrod. dirac or scalar_rel (1,2,3)
       info%wave_eq=SCHRODINGER!lets assume this for now
       !number of spin channels
       info%nspin=1!assume this for now
       !pseudo potential scheme one of HAM TM RTM MRPP (1,2,3,4)
       info%scheme=2!assume just TM for now

       !set the wavefunction stuff now
       info%have_wfs=.true.
       info%n_wfs=nwf
       allocate(info%wfs_n(nwfx))
       info%wfs_n(:)=nn(:)
       allocate(info%wfs_l(nwfx))
       info%wfs_l(:)=ll(:)
       allocate(info%wfs_j(nwfx))
       info%wfs_j(:)=jj(:)
       allocate(info%wfs_label(nwfx))
       info%wfs_label(:)=el(:)
       allocate(info%wfs_rc(nwfsx))!cut off for pseudos (NO US allowed)
       info%wfs_rc(:)=rcut(:)
       allocate(info%wfs_occ(nwfsx,info%nspin))
       !here we are using the variale configuration occupation listing, just the first one though
       info%wfs_occ(:,1)=octsc(:,1)!pretend there's only a single spin channel
       allocate(info%wfs_ev(nwfsx))
       !info%wfs_ev(:)=enltsc(:,1)
       info%wfs_ev(:)=enls(:)
       allocate(info%wfs(ndmx,nwfsx))
       info%wfs(:,:)=phis(:,:)
       allocate(info%rho_val(ndmx))
       info%rho_val(:)=rhos(:,1)
       !rhos(ndmx,2), i guess this is up and down density. with nspin=1 only 1 channel here is nonzero
       !(check this again later, i might have the wrong channel ^)

       !set the projector stuff
       info%have_kb=.true.
       info%kb_l_local=lloc
       info%kb_n_proj=nbeta
       !allocate(info%kb_l(info%kb_n_proj))
       !info%kb_l(:)=pawsetup%l(:)
       allocate(info%kb_v_local(ndmax))
       info%kb_v_local(:)=vpsloc(:)
       allocate(info%kb_proj(ndmx,nwfsx))
       info%kb_proj(:,:)=betas(:,:)
       allocate(info%kb_l(nbeta))
       info%kb_l=lls
       allocate(info%kb_j(nbeta))
       !just pretend that we can calculate this from j=l+/-s
       do i=1,nbeta,1!most likly wrong all the s values are 1
        info%kb_j(i)=lls(i)+isws(i)!isws is the spin of the pseudowavefunction
       enddo
       allocate(info%kb_e(nbeta))
       !just pretend that we can calculate this from wfs_ev
       do i=1,nbeta,1
        info%kb_e(i)=info%wfs_ev(i)
       enddo

       !set xc stuff
       info%ixc=1!this should be set as a constant for now, check code later on to see how to integrate it properly
       info%nlcc=nlcc
       allocate(info%rho_core(ndmx))
       info%rho_core(:)=rhoc(:)
       info%nlcc_rc=rcore

       !set pseudo_potentials
       info%have_psp=.true.
       info%n_psp=nwfs
       allocate(info%psp_l(nwfs))
       info%psp_l=lls(:nwfs)
       allocate(info%psp_j(nwfs))
       do i=1,nwfs,1
       info%psp_j(i)=i!not needed for non spin polarised 
       enddo
       !this is the semi-local potential. we can also ignore this 
       !allocate(info%psp_v(,))
       !info%psp_v(,)
      endsubroutine ps_io_type_fill

      subroutine ps_io_save(info)!this seems only necessary to setup the filename and assign io units
      !-----------------------------------------------------------------------!
      ! Writes the pseudo-atom information to a file so it can be used by     !
      ! other codes.                                                          !
      !-----------------------------------------------------------------------!
      integer :: unit, iostat
      character(len=20) :: filename
      type(ps_io_type), intent(in) :: info
  
      !Open file
      !call io_assign(unit)
      unit=999!just make sure this is out of the range of the ld1x (seems to be ok though)
      filename = trim(info%symbol)//'.psf'
      open(unit, file=trim(filename), status='unknown', iostat=iostat)
      if (iostat > 0) then
        print *, "Error: something wrong with filename"
        !write(message(1),'("Unable to open file ''",A,"''")') trim(filename)
        !call write_fatal(1)
      end if
  
      !Write data to the file
      call siesta_save(unit, .false.,info)
      write(unit,*) "This is siesta pseudo output"

      close(unit)

      end subroutine ps_io_save
  
    subroutine siesta_save(unit, parsec, info)
    !-----------------------------------------------------------------------!
    ! Write the pseudo-potential stuff using the format of the SIESTA code. !
    !-----------------------------------------------------------------------!
    integer, intent(in) :: unit
    logical, intent(in) :: parsec
    type(ps_io_type), intent(in) :: info

    integer :: i, l, n_dn, n_up, values(8), n, ir
    real(R8) :: j, occ
    character(3)    :: irel
    character(4)    :: icore
    character(60)   :: header
    character(70)   :: title
    type(mesh_type) :: new_m
    real(R8), allocatable :: dum(:)

    !Get spin/relativistic mode and number of states
    select case (info%wave_eq)
    case (DIRAC)
      irel = "rel"
      n_up = 0
      n_dn = 0
      do i = 1, info%n_psp
        if (info%psp_l(i) == 0 .or. (info%psp_l(i) /= 0 .and. &
             info%psp_j(i) == info%psp_l(i) - M_HALF)) then
          n_dn = n_dn + 1          
        elseif (info%psp_l(i) /= 0 .and. &
                (info%psp_j(i) == info%psp_l(i) + M_HALF)) then
          n_up = n_up + 1

        end if
      end do
    case (SCHRODINGER, SCALAR_REL)
      n_dn = info%n_psp
      select case (info%nspin)
      case (1)
        irel = "nrl"
        n_up = 0
      case (2)
        irel = "isp"
        n_up = n_dn
      end select
    end select
      
    call date_and_time(values=values)
    write(header,'(A,A,3X,I4.4,"/",I2.2,"/",I2.2)') "APE Version-", PACKAGE_VERSION, values(1), values(2), values(3)
    select case (info%scheme)
    case (HAM)
      header = trim(header)//"   Hamann"
    case (TM)
      header = trim(header)//"   Troullier-Martins"
    end select
    
    !Core corrections mode
    if (info%nlcc) then
      icore = "pcec"
    else
      icore = "nc  "
    end if

    !
    title = ""
    do l = 0, 3
      n = minval(info%wfs_n, mask=info%wfs_l == l)
      j = minval(info%wfs_j, mask=(info%wfs_l == l .and. info%wfs_n == n))

      do i = 1, info%n_wfs
        if (info%wfs_n(i) == n .and. info%wfs_l(i) == l .and. info%wfs_j(i) == j) then
          select case (irel)
          case('nrl')
            write(title,'(A)') trim(title)
            write(title,'(A,A2,F5.2,"  r=",F5.2,"/")') trim(title), &
                            info%wfs_label(i), info%wfs_occ(i, 1), info%wfs_rc(i)
          case('isp')
            write(title,'(A,A2,F4.2,1X,F4.2,1X,F4.2,"/")') trim(title), &
                     info%wfs_label(i), info%wfs_occ(i, 1), info%wfs_occ(i, 2), &
                     info%wfs_rc(i)
          case ('rel')
            occ = sum(info%wfs_occ(:, 1), mask=(info%wfs_l == l .and. &
                 (info%wfs_j == l + M_HALF .or. info%wfs_j == l - M_HALF)))
            write(title,'(A,A2,F5.2,"  r=",F5.2,"/")') trim(title), &
                                           info%wfs_label(i), occ, info%wfs_rc(i)
          end select
        end if
      end do

    end do

    !Mesh
    call mesh_null(new_m)
    call mesh_generation(new_m, LOG2, info%m%r(1), info%m%r(info%m%np), info%m%np) 
    
    !General info
    write(unit,'(1X,A2,1X,A2,1X,A3,1X,A4)') info%symbol(1:2), ixc_to_icorr(info%ixc), irel, icore
    write(unit,'(1X,A60)') header
    write(unit,'(1X,A70)') title
    write(unit,'(1X,2I3,I5,3F20.10)') n_dn, n_up, new_m%np, new_m%b, new_m%a, info%z_val
      
    !Write radial grid
    write(unit,'(" Radial grid follows")')
    write(unit,'(4(g20.12))') (new_m%r(i),i = 1, new_m%np)

    allocate(dum(new_m%np))

    !Down pseudopotentials
    do l = 0, 3
      if (info%wave_eq == DIRAC) then
        j = max(l - M_HALF, M_HALF)
      else
        j = M_ZERO
      end if

      do i = 1, info%n_psp
        if (info%psp_l(i) == l .and. info%psp_j(i) == j) then
          write(unit,'(" Down Pseudopotential follows (l on next line)")')
          write(unit,'(1X,I2)') info%psp_l(i)
          
          call mesh_transfer(info%m, info%psp_v(:, i), new_m, dum, 3)
          where (abs(dum) < 1.0E-30)
            dum = M_ZERO
          end where
          write(unit,'(4(g20.12))') (dum(ir)*new_m%r(ir)*M_TWO, ir = 1, new_m%np)
        end if
      end do
    end do

    !Up pseudopotentials
    if (n_up /= 0) then
      do l = 0, 3
        do i = 1, info%n_psp
          if (l /= 0 .and. info%wave_eq == DIRAC) then
            j = l + M_HALF
          else
            j = M_ZERO
          end if
          
          if (info%psp_l(i) == l .and. info%psp_j(i) == j) then
            write(unit,'(" Up Pseudopotential follows (l on next line)")')
            write(unit,'(1X,I2)') info%psp_l(i)
          
            call mesh_transfer(info%m, info%psp_v(:, i), new_m, dum, 3)
            where (abs(dum) < 1.0E-30)
              dum = M_ZERO
            end where
            write(unit,'(4(g20.12))') (dum(ir)*new_m%r(ir)*M_TWO, ir = 1, new_m%np)
          end if
        end do
      end do
    end if

    !Write core charge
    write(unit,'(" Core charge follows")')
    if (icore == "nc  ") then
      dum = M_ZERO
    else
      call mesh_transfer(info%m, info%rho_core, new_m, dum, 3)
      where (abs(dum) < 1.0E-30)
        dum = M_ZERO
      elsewhere
        dum = new_m%r**2*M_FOUR*M_PI*dum
      end where
    end if
    write(unit,'(4(g20.12))') (dum(i), i = 1, new_m%np)

    !Write valence charge
    write(unit,'(" Valence charge follows")')
    call mesh_transfer(info%m, info%rho_val, new_m, dum, 3)
    where (abs(dum) < 1.0E-30)
      dum = M_ZERO
    end where
    write(unit,'(4(g20.12))') (dum(i)*new_m%r(i)**2*M_FOUR*M_PI, i = 1, new_m%np)

    if (parsec) then
      !Write pseudo-wave-functions
      do l = 0, 3
        n = minval(info%wfs_n, mask=info%wfs_l == l)
        j = minval(info%wfs_j, mask=(info%wfs_l == l .and. info%wfs_n == n))

        do i = 1, info%n_wfs
          if (info%wfs_n(i) == n .and. info%wfs_l(i) == l .and. info%wfs_j(i) == j) then

            call mesh_transfer(info%m, info%wfs(:, i), new_m, dum, 3)

            write(unit,'(1X,A,A2)') 'Pseudo-wave-function follows (l, zelect, rc)  ',info%wfs_label(i)
            write(unit,'(I2,F6.2,2X,F6.2)') info%wfs_l(i), &
                 info%wfs_occ(i, 1), info%wfs_rc(i)
            write(unit,'(1P4E19.11)') (dum(ir)*new_m%r(ir), ir = 1, new_m%np)

          end if
        end do

      end do
    end if
    
    deallocate(dum)
    call mesh_end(new_m)
    
  end subroutine siesta_save
  end module siesta_pass

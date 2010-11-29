      module siesta_pass
      !APE modules (modified)
      use global
      use mesh, only: mesh_null, mesh_type
      !QE modules
      use radial_grids, only: radial_grid_type,ndmx
      use ld1inc, only: zed,enne,zval,nwf,nn,ll,jj,el,rcut,octsc,enltsc,&
      phis,lloc,nbeta,rhos,pawsetup,enls,vpsloc,betas,rho,lls,isws,nlcc,&
      rhoc,rcore
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
!!  
!!      !Pseudo potentials
!!      logical :: have_psp
!!      integer :: n_psp
!!      integer,  pointer :: psp_l(:)
!!      real(R8), pointer :: psp_j(:)
!!      real(R8), pointer :: psp_v(:,:)
!!  
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

      subroutine siesta_output(grid_in)
      type(radial_grid_type), intent(in):: grid_in
      print *, "here's some siesta output"
      call ps_io_type_fill(grid_in)
      endsubroutine siesta_output

      subroutine ps_io_type_fill(grid_in)
       type(radial_grid_type), intent(in):: grid_in
       type(ps_io_type) :: info
       !mesh parameters
       type(mesh_type) :: m
       integer :: no_mesh_points,i,j
       real(R8) :: rn,r1
       !generalities
       character, external :: atom_name*2

       info%file_format=M_SIESTA! remember to include the M_SIESTA etc definitions at the top of this file

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
       !test

       print *, "xc"
       print *, info%nlcc
       print *, info%rho_core

      endsubroutine ps_io_type_fill
      
      end module siesta_pass

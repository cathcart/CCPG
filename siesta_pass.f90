      module siesta_pass
      use global
      use mesh
      type ps_io_type
      private
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
      integer :: n_psp
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
      integer :: ixc
  
      !Non-linear core-corrections
      logical :: nlcc
      real(R8) :: nlcc_rc
      real(R8), pointer :: rho_core(:)
  
      !Kleinman-Bylander projectors
      logical :: have_kb
      integer :: kb_l_local
      integer :: kb_n_proj 
      integer,  pointer :: kb_l(:)
      real(R8), pointer :: kb_j(:)
      real(R8), pointer :: kb_v_local(:)
      real(R8), pointer :: kb_e(:)
      real(R8), pointer :: kb_proj(:,:)
      end type ps_io_type
 
      contains

      subroutine siesta_output()
      print *, "here's some siesta output"
      endsubroutine siesta_output

      subroutine ps_io_type_fill()
      endsubroutine ps_io_type_fill
      
      end module siesta_pass

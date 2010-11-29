!! Copyright (C) 2004-2007 M. Oliveira, F. Nogueira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

!#include "global.h"

module mesh
  use global
  !use oct_parser
  !use messages
  use splines
  !use units
  implicit none


                    !---Interfaces---!

!!  interface assignment (=)
!!     module procedure mesh_copy
!!  end interface
!!
!!  interface operator (==)
!!     module procedure equal_mesh
!!  end interface


                    !---Derived Data Types---!

  type mesh_type
    integer           :: type !mesh type
    real(R8)          :: a    !mesh parameters
    real(R8)          :: b    !
    integer           :: np   !mesh number of points
    real(R8), pointer :: r(:) !mesh points
  end type mesh_type


                    !---Global Variables---!

  !Mesh types
  integer, parameter :: LIN  = 1, &
                        LOG1 = 2, &
                        LOG2 = 3


                    !---Public/Private Statements---!

  private
  public :: mesh_type, &
            mesh_null, &
!!            mesh_init, &
            mesh_generation, &
!!            mesh_save, &
!!            mesh_load, &
            mesh_transfer, &
!!            mesh_output_params, &
            mesh_end, &
!!            assignment(=), &
!!            operator(==), &
            LOG1, &
            LIN, &
            LOG2

contains

  subroutine mesh_null(m)
    !-----------------------------------------------------------------------!
    ! Nullifies and sets to zero all the components of mesh m.              !
    !-----------------------------------------------------------------------!
    type(mesh_type), intent(out) :: m

    !call push_sub("mesh_null")!this just looks like a timing function, comment out

    m%type = 0
    m%a    = M_ZERO
    m%b    = M_ZERO
    m%np   = 0
    nullify(m%r)

    !call pop_sub()
  end subroutine mesh_null

!!  subroutine mesh_init(z, m)
!!    !-----------------------------------------------------------------------!
!!    ! Initializes the calculation mesh by reading the mesh parameter from   !
!!    ! the input file. The default parameters depend on the nuclear charge.  !
!!    !                                                                       !
!!    !  z - nuclear charge                                                   !
!!    !  m - mesh object                                                      !
!!    !-----------------------------------------------------------------------!
!!    real(R8),        intent(in)    :: z
!!    type(mesh_type), intent(inout) :: m
!!
!!    integer  :: type, is_def_np, is_def_a, np
!!    real(R8) :: mesh_r1, mesh_rnp, a
!!
!!    !call push_sub("mesh_init")
!!
!!    !ASSERT(m%type == 0)
!!
!!    !Read input options
!!    call oct_parse_int('MeshType', LOG1, type)
!!    select case (type)
!!    case (LOG1, LOG2, LIN)
!!    case default
!!      !.(1) =  "Illegal MeshType."
!!      call write_fatal(1)
!!    end select
!!
!!    call oct_parse_double('MeshStartingPoint', sqrt(z)*1.0E-5, mesh_r1)
!!    if (mesh_r1 <= M_ZERO) then
!!      !.(1) =  "MeshStartingPoint can not take negative values."
!!      call write_fatal(1)
!!    end if
!!    !mesh_r1 = mesh_r1*units_in%length%factor!assume for now that we can simply operate in on fixed unit type => factor ==1 or constant
!!    mesh_r1 = mesh_r1*M_ONE!QE uses atomic Ry units internally, that is bohr for length
!!
!!    call oct_parse_double('MeshOutmostPoint', sqrt(z)*M_THIRTY, mesh_rnp)
!!    if (mesh_rnp <= M_ZERO .or. mesh_rnp < mesh_r1) then
!!      !.(1) = "MeshOutmostPoint can''t take negative values"
!!      !.(2) = "and must be greater than the mesh start point."
!!      call write_fatal(2)
!!    end if
!!    !mesh_rnp = mesh_rnp*units_in%length%factor
!!    mesh_rnp = mesh_rnp*M_ONE
!!    
!!    is_def_np = oct_parse_isdef('MeshNumberOfPoints')
!!    is_def_a  = oct_parse_isdef('MeshParameter')
!!    if (is_def_np == 1 .and. is_def_a == 1) then
!!      !.(1) = "MeshNumberOfPoints and MeshParameter input options"
!!      !.(2) = "can not be used at the same time."
!!      call write_fatal(2)
!!    end if
!!
!!    if (is_def_a == 0) then
!!      call oct_parse_int('MeshNumberOfPoints', int(sqrt(z))*200, np)
!!      if (np <= 3) then
!!        !.(1) = "Mesh must have at least 3 points."
!!        call write_fatal(1)
!!      end if
!!    else
!!      call oct_parse_double('MeshParameter', M_DIME, a)
!!      if (a < M_ZERO) then
!!        !.(1) = "MeshParameter can not take negative values."
!!        call write_fatal(1)
!!      end if
!!    end if
!!
!!    !Generate mesh
!!    if (is_def_a == 0) then
!!      call mesh_generation(m, type, mesh_r1, mesh_rnp, n=np)
!!    else
!!      call mesh_generation(m, type, mesh_r1, mesh_rnp, a=a)
!!    end if
!!
!!    !call pop_sub()
!!  end subroutine mesh_init
!!
  subroutine mesh_generation(m, type, r1, rn, n, a)
    !-----------------------------------------------------------------------!
    ! Initializes a mesh and generates the mesh points. Note that only one  !
    ! of the optional variables n and a should be used.                     !
    !                                                                       !
    !  m    - mesh                                                          !
    !  type - mesh type                                                     !
    !  r1   - starting point                                                !
    !  rn   - ending point                                                  !
    !  n    - number of points                                              !
    !  a    - mesh parameter a                                              !
    !-----------------------------------------------------------------------!
    type(mesh_type), intent(inout) :: m
    integer,         intent(in)    :: type
    real(R8),        intent(in)    :: r1, rn
    integer,         intent(in), optional :: n
    real(R8),        intent(in), optional :: a

    integer :: i
    real(R8) :: a1, a2, n1, n2, am, nm, f1, fm

    !call push_sub("mesh_generation")

    !Check optional arguments
    if (present(n) .and. present(a)) then
      !.(1) = "Only one optional argument can be set in mesh_generation"
      print *, "Something wrong with the mesh arguments"
      !call write_fatal(1)
    end if

    !Set the mesh type, the number of points, the first point and the parameter a
    m%type = type
    if (present(n)) then
      m%np = n
      allocate(m%r(m%np))
      m%r(1) = r1
      select case (m%type)
      case (LIN)
        m%a = (rn - r1)/real(n-1, R8)
      case (LOG1)
        m%a = log(rn/r1)/real(n - 1,R8)
      case (LOG2)
        a1 = 1.0e-8_r8
        f1 = func(r1, rn, real(n,R8), a1)
        a2 = M_ONE
        do
          am = (a2 + a1)*M_HALF
          fm = func(r1, rn, real(n,R8), am)
          if (M_HALF*abs(a1 - a2) < 1.0e-12) exit
          if (fm*f1 > M_ZERO) then
            a1 = am; f1 = fm
          else
            a2 = am
          end if
        end do
        m%a = am
      end select

    elseif (present(a)) then
      select case (m%type)
      case (LIN)
        m%np = int((rn - r1)/a)
      case (LOG1)
        m%np = int(log(rn/r1)/a) + 1
      case (LOG2)
        n1 = M_ONE
        f1 = func(r1, rn, n1, a)
        n2 = 10000.0_r8
        do
          nm = (n2 + n1)*M_HALF
          fm = func(r1, rn, nm, a)
          if (M_HALF*abs(n1 - n2) < 1.0e-12) exit
          if (fm*f1 > M_ZERO) then
            n1 = nm; f1 = fm
          else
            n2 = nm
          end if
        end do
        m%np = int(nm)
      end select
      allocate(m%r(m%np))
      m%r(1) = r1
      m%a = a

    end if

    !Set the parameter b and the remaining points
    select case (m%type)
    case (LIN)
      m%b = M_ZERO
      do i = 2, m%np
        m%r(i) = m%r(i-1) + m%a
      end do
    case (LOG1)
      m%b = r1/exp(m%a)
      do i = 2, m%np
        m%r(i) = exp(m%a)*m%r(i-1)
      end do
    case (LOG2)
      m%b = r1/(exp(m%a) - M_ONE)
      do i = 2, m%np
        m%r(i) = m%r(i-1)*exp(m%a) + r1
      end do
    end select

    !call pop_sub()!this is also a timing function
  contains

    real(R8) function func(r1, rn, n, a)
      real(R8), intent(in) :: r1, rn, a, n
      func = exp(n*a)*r1 - M_ONE*r1 - rn*exp(a) + rn*M_ONE
    end function func

  end subroutine mesh_generation
!!
!!  subroutine mesh_copy(m_a, m_b)
!!    !-----------------------------------------------------------------------!
!!    ! Copies m_a mesh to m_b.                                               !
!!    !-----------------------------------------------------------------------!
!!    type(mesh_type), intent(inout) :: m_a
!!    type(mesh_type), intent(in)    :: m_b
!!
!!    !call push_sub("mesh_copy")
!!
!!    call mesh_end(m_a)
!!    m_a%type = m_b%type
!!    m_a%a    = m_b%a
!!    m_a%b    = m_b%b
!!    m_a%np   = m_b%np
!!    allocate(m_a%r(m_a%np))
!!    m_a%r = m_b%r
!!
!!    !call pop_sub()
!!  end subroutine mesh_copy
!!
!!  function equal_mesh(m_a, m_b)
!!    !-----------------------------------------------------------------------!
!!    ! Returns true if m_a and m_b are the same mesh; false otherwise.       !
!!    !-----------------------------------------------------------------------!
!!    type(mesh_type), intent(in) :: m_a, m_b
!!    logical :: equal_mesh
!!
!!    !call push_sub("equal_mesh")
!!
!!    equal_mesh = (m_a%type == m_b%type .and. m_a%a == m_b%a .and. &
!!                  m_a%b == m_b%b .and. m_a%np == m_b%np)
!!
!!    !call pop_sub()
!!  end function equal_mesh
!!
!!  subroutine mesh_save(unit, m)
!!    !-----------------------------------------------------------------------!
!!    ! Writes the mesh information to a file.                                !
!!    !                                                                       !
!!    !  unit - file unit number                                              !
!!    !  m    - mesh to be written                                            !
!!    !-----------------------------------------------------------------------!
!!    integer,         intent(in) :: unit
!!    type(mesh_type), intent(in) :: m
!!
!!    integer :: i
!!
!!    !call push_sub("mesh_save")
!!
!!    !ASSERT(m%type /= 0)
!!
!!    write(unit) m%type, m%a, m%b, m%np
!!    write(unit) (m%r(i), i=1, m%np)
!!
!!    !call pop_sub()
!!  end subroutine mesh_save
!!
!!  subroutine mesh_load(unit, m)
!!    !-----------------------------------------------------------------------!
!!    ! Reads the mesh information from a file.                               !
!!    !                                                                       !
!!    !  unit - file unit number                                              !
!!    !  m    - mesh to be read                                               !
!!    !-----------------------------------------------------------------------!
!!    integer,         intent(in)    :: unit
!!    type(mesh_type), intent(inout) :: m
!!
!!    integer :: i
!!
!!    !call push_sub("mesh_load")
!!
!!    !ASSERT(m%type == 0)
!!
!!    read(unit) m%type, m%a, m%b, m%np
!!    allocate(m%r(m%np))
!!    read(unit) (m%r(i), i=1, m%np)
!!
!!    !call pop_sub()
!!  end subroutine mesh_load
!!
  subroutine mesh_transfer(m_a, fa, m_b, fb, interp_type)
    !-----------------------------------------------------------------------!
    ! Having a function on a mesh m_a, this routines returns the values of  !
    ! that function on a mesh m_b, by interpolating the function.           !
    !                                                                       !
    !  m_a         - mesh m_a                                               !
    !  fa          - values of the function on mesh A                       !
    !  x_b         - mesh m_b                                               !
    !  fa          - values of the function on mesh A                       !
    !  interp_type - interpolation type                                     !
    !-----------------------------------------------------------------------!
    type(mesh_type), intent(in)  :: m_a, m_b
    real(R8),        intent(in)  :: fa(m_a%np)
    real(R8),        intent(out) :: fb(m_b%np)
    integer,         intent(in)  :: interp_type

    !call push_sub("mesh_transfer")

    !ASSERT(m_a%type /= 0 .and. m_b%type /= 0)

    call spline_mesh_transfer(m_a%np, m_a%r, fa, m_b%np, m_b%r, fb, interp_type)

    !call pop_sub()
  end subroutine mesh_transfer
!!
!!  subroutine mesh_output_params(m, unit, verbose_limit)
!!    !-----------------------------------------------------------------------!
!!    ! Writes the mesh input options and the mesh paremeters in a nice       !
!!    ! readable format. If a unit number is provided, it writes it to a file !
!!    ! otherwise it writes it to the standard output.                        !
!!    !-----------------------------------------------------------------------!
!!    type(mesh_type), intent(in) :: m
!!    integer,         intent(in), optional :: unit, verbose_limit
!!
!!    !call push_sub("mesh_output_params")
!!
!!    !ASSERT(m%type /= 0)
!!
!!    !message(1) = ""
!!    !message(2) = "Mesh information:"
!!    select case (m%type)
!!    case (LIN)
!!      !message(3) = "  Type: linear"
!!    case (LOG1)
!!      !message(3) = "  Type: logarithmic [ri = b*exp(a*i)]"
!!    case (LOG2)
!!      !message(3) = "  Type: logarithmic [ri = b*(exp(a*i) - 1)]"
!!    end select
!!    !write(message(4), '("  Mesh starting point:   ",ES8.2E2,1X,A)') m%r(1)/units_out%length%factor, trim(units_out%length%abbrev)
!!    !write(message(5), '("  Mesh outmost point:    ",F7.3,1X,A)') m%r(m%np)/units_out%length%factor, trim(units_out%length%abbrev)
!!    !write(message(6), '("  Mesh parameters (a, b): ",ES12.5E2,", ",ES12.5E2)') m%a, m%b
!!
!!!!    if (present(unit)) then !not too sure this is necessary
!!!!      call write_info(6, unit=unit)
!!!!    else
!!!!      if (present(verbose_limit)) then
!!!!        call write_info(6, verbose_limit)
!!!!      else
!!!!        call write_info(6)
!!!!      end if
!!!!    end if
!!!! 
!!!!    call pop_sub() !this is a timer thing, not interested
!!  end subroutine mesh_output_params
!!
  subroutine mesh_end(m)
    !-----------------------------------------------------------------------!
    ! Frees all memory associated to the mesh object m.                     !
    !-----------------------------------------------------------------------!
    type(mesh_type), intent(inout) :: m

    !call push_sub("mesh_end")

    m%type = 0
    m%a    = M_ZERO
    m%b    = M_ZERO
    m%np   = 0
    if (associated(m%r)) then
      deallocate(m%r)
    end if

    !call pop_sub()
  end subroutine mesh_end

end module mesh

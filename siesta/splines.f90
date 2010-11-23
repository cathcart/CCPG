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

module splines
  use global
  use gsl_interface
  !use messages
  implicit none


                    !---Interfaces---!

  interface assignment (=)
     module procedure copy_splines
  end interface

  interface operator (==)
     module procedure equal_splines
  end interface


                    !---Derived Data Types---!

  type spline_type
    !This derived data type is the only one acessible to other modules.
    !It is just a pointer to the splines information
    private
    type(spline_info_type), pointer :: info
  end type spline_type

  type spline_info_type
    !acc and spl are the gsl objects
    !n_ptrs is the number of spline_type objects pointing to this spline
    integer(POINTER_SIZE) :: acc, spl
    integer               :: n_ptrs
  end type spline_info_type


                    !---Public/Private Statements---!

  private
  public :: spline_type, &
            spline_null, &
            spline_init, &
            assignment(=), &
            operator(==), &
            spline_end, &
            spline_eval, &
            spline_eval_deriv, &
            spline_eval_deriv2, &
            spline_eval_integ, &
            spline_integ, &
            spline_deriv, &
            spline_deriv2, &
            spline_mesh_transfer


contains

  subroutine spline_null(spline)
    !-----------------------------------------------------------------------!
    ! Nullifies the components of the spline.                               !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(out) :: spline

    nullify(spline%info)

  end subroutine spline_null

  subroutine spline_init(spline, n, x, y, interp_type)
    !-----------------------------------------------------------------------!
    ! Given a set of data points (x_i,y_i) it computes an interpolation     !
    ! function of a chosen type.                                            !
    !                                                                       !
    !  spline      - pointer to the spline information                      !
    !  n           - number of data-points                                  !
    !  x           - x coordinates of the data points                       !
    !  y           - y coordinates of the data points                       !
    !  interp_type - interpolation type. Possible values are:               !
    !                 1 -> linear                                           !
    !                 2 -> polynomial                                       !
    !                 3 -> cubic spline with natural boundary conditions    !
    !                 4 -> cubic spline with periodic boundary conditions   !
    !                 5 -> Akima spline with natural boundary conditions    !
    !                 6 -> Akima spline with periodic boundary conditions   !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(inout) :: spline
    integer,           intent(in) :: interp_type, n
    real(R8),          intent(in) :: x(n), y(n)

    real(R8), allocatable :: auxx(:), auxy(:)

    !ASSERT(.not.associated(spline%info))!this is an error message thing

    !Allocate memory
    allocate(spline%info)

    !Store all the information about the interpolation function
    spline%info%n_ptrs = 1
    spline%info%acc = 0
    spline%info%spl = 0
    call gsl_interp_accel_alloc(spline%info%acc)
    call gsl_spline_alloc(spline%info%spl, n, interp_type)
    allocate(auxx(n), auxy(n))
    auxy = y !The use of an auxiliary array prevents some types of errors
    auxx = x
    call gsl_spline_init(spline%info%spl, n, auxx(1), auxy(1))
    deallocate(auxx, auxy)

  end subroutine spline_init

  subroutine copy_splines(spline_a, spline_b)
    !-----------------------------------------------------------------------!
    ! Makes spline_a equal to spline_b. That is done by making both splines !
    ! component info point to the same interpolation information.           !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(inout) :: spline_a
    type(spline_type), intent(in)    :: spline_b

    !ASSERT(associated(spline_b%info))

    !Deallocate the memory previously assigned to spline_a
    call spline_end(spline_a)

    !Make spline_a point to the same interpolation info than spline_b
    spline_a%info => spline_b%info

    !Now we have an addicional pointer to the spline information
    spline_a%info%n_ptrs = spline_a%info%n_ptrs + 1

  end subroutine copy_splines

  function equal_splines(spline_a, spline_b)
    !-----------------------------------------------------------------------!
    ! Returns true if spline_a and spline_b interpolation objects are equal.!
    ! Returns false if they are not.                                        !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(in) :: spline_a, spline_b
    logical :: equal_splines

    equal_splines = associated(spline_a%info, spline_b%info)

  end function equal_splines

  subroutine spline_end(spline)
    !-----------------------------------------------------------------------!
    ! Frees the interpolation object spline.                                !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(inout) :: spline

    if (associated(spline%info)) then
      if (spline%info%n_ptrs > 1) then
        !There are more pointers pointing to this spline information
        spline%info%n_ptrs = spline%info%n_ptrs - 1
        nullify(spline%info)
      else
        !This is the only pointer pointing to this spline information
        call gsl_spline_free(spline%info%spl)
        call gsl_interp_accel_free(spline%info%acc)
        deallocate(spline%info)
      end if
    end if

  end subroutine spline_end

  function spline_eval(spline, x)
    !-----------------------------------------------------------------------!
    ! Returns the interpolated value for a given point.                     !
    !                                                                       !
    !  spline      - interpolation object                                   !
    !  x           - point where to evaluate the function                   ! 
    !-----------------------------------------------------------------------!
    type(spline_type), intent(in) :: spline
    real(R8),          intent(in) :: x
    real(R8) :: spline_eval

    spline_eval = gsl_spline_eval(x, spline%info%spl, spline%info%acc)
 
  end function spline_eval

  function spline_eval_deriv(spline, x)
    !-----------------------------------------------------------------------!
    ! Returns the derivative of an interpolated function for a given point. !
    !                                                                       !
    !  spline            - interpolation object                             !
    !  x                 - point where to evaluate the derivative           ! 
    !-----------------------------------------------------------------------!
    type(spline_type), intent(in) :: spline
    real(R8),          intent(in) :: x
    real(R8) :: spline_eval_deriv
    
    spline_eval_deriv = gsl_spline_eval_deriv(x, spline%info%spl, spline%info%acc)
 
  end function spline_eval_deriv

  function spline_eval_deriv2(spline, x)
    !-----------------------------------------------------------------------!
    ! Returns the second derivative of an interpolated function for a given !
    ! point.                                                                !
    !                                                                       !
    !  spline            - interpolation object                             !
    !  x                 - point where to evaluate the second derivativ     !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(in) :: spline
    real(R8),          intent(in) :: x
    real(R8) :: spline_eval_deriv2
    
    spline_eval_deriv2 = gsl_spline_eval_deriv2(x, spline%info%spl, spline%info%acc)
 
  end function spline_eval_deriv2

  function spline_eval_integ(spline, a, b)
    !-----------------------------------------------------------------------!
    ! Returns the numerical integral of an interpolated function over the   !
    ! range [a,b].                                                          !
    !                                                                       !
    !  spline - interpolation object                                        !
    !  a, b   - upper and lower integration bounds                          !
    !-----------------------------------------------------------------------!
    type(spline_type), intent(in) :: spline
    real(R8),          intent(in) :: a, b
    real(R8) :: spline_eval_integ

    spline_eval_integ = gsl_spline_eval_integ(a, b, spline%info%spl, spline%info%acc)
 
  end function spline_eval_integ

  function spline_integ(np, x, y, a, b, interp_type)
    !-----------------------------------------------------------------------!
    ! Returns the numerical integral of a function over the range [1,] by   !
    ! interpolating a set of data points.                                   !
    !                                                                       !
    !  np          - number of data points                                  !
    !  x           - x coordinates of the data points                       !
    !  y           - y coordinates of the data points                       !
    !  a, b        - upper and lower integration bounds                     !
    !  interp_type - interpolation type                                     !
    !-----------------------------------------------------------------------!
    integer,  intent(in) :: np, interp_type
    real(R8), intent(in) :: x(np), y(np)
    real(R8), intent(in) :: a, b
    real(R8) :: spline_integ

    type(spline_type) :: spl
    
    call spline_null(spl)
    call spline_init(spl, np, x, y, interp_type)
    spline_integ = spline_eval_integ(spl, a, b)
    call spline_end(spl)

  end function spline_integ

  function spline_deriv(np, x, y, interp_type)
    !-----------------------------------------------------------------------!
    ! Returns the numerical derivative of a function by interpolating a set !
    ! of data points.                                                       !
    !                                                                       !
    !  np          - number of data points                                  !
    !  x           - x coordinates of the data points                       !
    !  y           - y coordinates of the data points                       !
    !  interp_type - interpolation type                                     !
    !-----------------------------------------------------------------------!
    integer,  intent(in)  :: np, interp_type
    real(R8), intent(in)  :: x(np), y(np)
    real(R8) :: spline_deriv(np)

    integer :: i
    type(spline_type) :: spl
    
    call spline_null(spl)
    call spline_init(spl, np, x, y, interp_type)
    do i = 1, np
      spline_deriv(i) = spline_eval_deriv(spl, x(i))
    end do
    call spline_end(spl)

  end function spline_deriv

  function spline_deriv2(np, x, y, interp_type)
    !-----------------------------------------------------------------------!
    ! Returns the numerical second derivative of a function by              !
    ! interpolating a set of data points.                                   !
    !                                                                       !
    !  np          - number of data points                                  !
    !  x           - x coordinates of the data points                       !
    !  y           - y coordinates of the data points                       !
    !  interp_type - interpolation type                                     !
    !  d2ydx2      - second derivative of y(x)                              !
    !-----------------------------------------------------------------------!
    integer,  intent(in)  :: np, interp_type
    real(R8), intent(in)  :: x(np), y(np)
    real(R8) :: spline_deriv2(np)

    integer :: i
    type(spline_type) :: spl
    
    call spline_null(spl)
    call spline_init(spl, np, x, y, interp_type)
    do i = 1, np
      spline_deriv2(i) = spline_eval_deriv2(spl, x(i))
    end do
    call spline_end(spl)

  end function spline_deriv2

  subroutine spline_mesh_transfer(np_a, x_a, y_a, np_b, x_b, y_b, interp_type)
    !-----------------------------------------------------------------------!
    ! Having a function on a given mesh A, this routines returns the values !
    ! of that function on a mesh B, by interpolating the function.          !
    !                                                                       !
    !  np_a        - mesh A number of points                                !
    !  x_a         - mesh A points                                          !
    !  y_a         - values of the function on mesh A                       !
    !  np_a        - mesh B number of points                                !
    !  x_b         - mesh B points                                          !
    !  y_b         - values of the function on mesh B                       !
    !  interp_type - interpolation type                                     !
    !-----------------------------------------------------------------------!
    integer,  intent(in)  :: np_a, np_b, interp_type
    real(R8), intent(in)  :: x_a(np_a), y_a(np_a), x_b(np_b)
    real(R8), intent(out) :: y_b(np_b)

    integer :: i
    type(spline_type) :: spl

    !call push_sub("spline_mesh_transfer")

    call spline_null(spl)
    call spline_init(spl, np_a, x_a, y_a, interp_type)
    do i = 1, np_b
      y_b(i) = spline_eval(spl, x_b(i))
    end do
    call spline_end(spl)

    !call pop_sub()
  end subroutine spline_mesh_transfer

end module splines

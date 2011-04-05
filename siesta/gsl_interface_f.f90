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

! This module contains the explicit interfaces for the external GSL routines.
! It also contains a very short description of the routines. 
! For more information please refer to the GSL documentation.


module gsl_interface
implicit none

  ! Error Handling
  interface
    subroutine gsl_strerror(err, res)
      !Returns a string describing the error code
      use global
      implicit none
      integer, intent(in) :: err
      character(len=*), intent(out) :: res
    end subroutine gsl_strerror
  end interface

  ! Mathematical Functions
  interface
    elemental function gsl_asinh(x)
      !Computes the value of arcsinh(x)
      use global
      implicit none
      real(R8), intent(in) :: x
      real(R8) :: gsl_asinh
    end function gsl_asinh
  end interface


  ! Special Functions
  interface
    function gsl_sf_bessel_knu_scaled(xnu, x)
      !Computes the scaled irregular modified Bessel function of 
      !fractional order nu, \exp(+|x|) K_\nu(x) for x>0, \nu>0.
      use global
      implicit none
      real(R8), intent(in) :: x, xnu
      real(R8) :: gsl_sf_bessel_knu_scaled
    end function gsl_sf_bessel_knu_scaled
  end interface


  ! Vectors and Matrices
  interface
    subroutine gsl_vector_alloc(n, v)
      !Creates a vector of length n, returning a pointer
      !to a newly initialized vector struct
      use global
      implicit none
      integer, intent(in) :: n
      integer(POINTER_SIZE), intent(inout) :: v
    end subroutine gsl_vector_alloc

    subroutine gsl_vector_free(v)
      !Frees a previously allocated vector v
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: v
    end subroutine gsl_vector_free

    subroutine gsl_vector_set(v, i, x)
      !Sets the value of the i-th element of a vector v to x
      use global
      implicit none
      integer, intent(in) :: i
      real(R8), intent(in) :: x
      integer(POINTER_SIZE), intent(inout) :: v
    end subroutine gsl_vector_set

    pure function gsl_vector_get(v, i)
      !Returns the i-th element of a vector v
      use global
      implicit none
      integer, intent(in) :: i
      integer(POINTER_SIZE), intent(in) :: v
      real(R8) :: gsl_vector_get
    end function gsl_vector_get

    subroutine gsl_matrix_alloc(n1, n2, m)
      !Creates a matrix of size n1 rows by n2 columns, returning 
      !a pointer to a newly initialized matrix struct
      use global
      implicit none
      integer, intent(in) :: n1, n2
      integer(POINTER_SIZE), intent(inout) :: m
    end subroutine gsl_matrix_alloc

    subroutine gsl_matrix_calloc(n1, n2, m)
      !Allocates memory for a matrix of size n1 rows by n2 columns 
      !and initializes all the elements of the matrix to zero
      use global
      implicit none
      integer, intent(in) :: n1, n2
      integer(POINTER_SIZE), intent(inout) :: m
    end subroutine gsl_matrix_calloc

    subroutine gsl_matrix_free(m)
      !Frees a previously allocated matrix m
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: m
    end subroutine gsl_matrix_free

    subroutine gsl_matrix_set(m, i, j, x)
      !Sets the value of the (i,j)th element of a matrix m to x
      use global
      implicit none
      integer, intent(in) :: i, j
      real(R8), intent(in) :: x
      integer(POINTER_SIZE), intent(inout) :: m
    end subroutine gsl_matrix_set

    pure function gsl_matrix_get(m, i, j)
      !Returns the (i,j)th element of a matrix m
      use global
      implicit none
      integer, intent(in) :: i, j
      integer(POINTER_SIZE), intent(in) :: m
      real(R8) :: gsl_matrix_get
    end function gsl_matrix_get
  end interface


  ! Permutations
  interface
    subroutine gsl_permutation_alloc(n, p)
      !Allocates memory for a new permutation of size n
      use global
      implicit none
      integer, intent(in) :: n
      integer(POINTER_SIZE), intent(inout) :: p
    end subroutine gsl_permutation_alloc

    subroutine gsl_permutation_free(p)
      !Frees all the memory used by the permutation p
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: p
    end subroutine gsl_permutation_free
  end interface


  ! Linear Algebra
  interface
    function gsl_linalg_lu_decomp(m, p, signum)
      !Factorizes the square matrix A into the LU decomposition PA = LU
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: m, p
      integer, intent(out) :: signum
      integer :: gsl_linalg_lu_decomp
    end function gsl_linalg_lu_decomp

    function gsl_linalg_lu_invert(lu, p, inverse) 
      !Computes the inverse of a matrix A from its LU decomposition (LU,p)
      use global
      implicit none
      integer(POINTER_SIZE), intent(in) :: lu, p
      integer(POINTER_SIZE), intent(inout) :: inverse
      integer :: gsl_linalg_lu_invert
    end function gsl_linalg_lu_invert

    function gsl_linalg_lu_solve(a, p, b, x)
      !Solves the system A x = b using the LU decomposition of A into (LU, p)
      use global
      implicit none
      integer(POINTER_SIZE), intent(in) :: b
      integer(POINTER_SIZE), intent(inout) :: a, p, x
      integer :: gsl_linalg_lu_solve
    end function gsl_linalg_lu_solve
  end interface


  ! Ordinary Differential Equations
  interface
    subroutine gsl_odeiv_step_alloc(stepping_func, dim, stp)
      !Returns a pointer to a newly allocated instance of a stepping function
      use global
      implicit none
      integer, intent(in) :: stepping_func, dim
      integer(POINTER_SIZE), intent(inout) :: stp
    end subroutine gsl_odeiv_step_alloc

    subroutine gsl_odeiv_step_free(stp)
      !Frees all the memory associated with the stepping function stp
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: stp
    end subroutine gsl_odeiv_step_free

    subroutine gsl_odeiv_evolve_alloc(dim, evl)
      !Returns a pointer to a newly allocated instance of an evolution function
      use global
      implicit none
      integer, intent(in) :: dim
      integer(POINTER_SIZE), intent(inout) :: evl
    end subroutine gsl_odeiv_evolve_alloc

    function gsl_odeiv_evolve_reset(evl)
      !Resets the evolution function evl
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: evl
      integer :: gsl_odeiv_evolve_reset
    end function gsl_odeiv_evolve_reset

    subroutine gsl_odeiv_evolve_free(evl)
      !Frees all the memory associated with the evolution function evl
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: evl
    end subroutine gsl_odeiv_evolve_free

    subroutine gsl_odeiv_control_standart_new(ctrl, eps_abs, eps_rel, a_y, a_dydt)
      !Creates a new control function. The standard control object is a four 
      !parameter heuristic based on absolute and relative errors eps_abs and
      !eps_rel, and scaling factors a_y and a_dydt for the system state y(t) 
      !and derivatives dydt(t) respectively
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: ctrl
      real(R8), intent(in) :: eps_abs, eps_rel, a_y, a_dydt
    end subroutine gsl_odeiv_control_standart_new

    subroutine gsl_odeiv_control_free(ctrl)
      !Frees all the memory associated with the control function ctrl
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: ctrl
    end subroutine gsl_odeiv_control_free
  end interface


  ! Interpolation
  interface
    subroutine gsl_spline_alloc(spl, n, opt)
      !Returns a pointer to a newly allocated interpolation object
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: spl
      integer, intent(in) :: opt, n
    end subroutine gsl_spline_alloc

    subroutine gsl_spline_init(spl, n, x, f)
      !Initializes the interpolation object spl for the data (x,f)
      use global
      implicit none
      integer, intent(in) :: n
      real(R8), intent(in) :: x, f
      integer(POINTER_SIZE), intent(inout) :: spl
    end subroutine gsl_spline_init

    subroutine gsl_spline_free(spl)
      !Frees the interpolation object spl
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: spl
    end subroutine gsl_spline_free

    subroutine toms_gsl_interp_accel_alloc(acc)
      !Returns a pointer to an accelerator object, which
      !is a kind of iterator for interpolation lookups
      use global
      implicit none
      integer(POINTER_SIZE), intent(in) :: acc
      !integer(POINTER_SIZE), intent(inout) :: acc
    end subroutine toms_gsl_interp_accel_alloc

    subroutine gsl_interp_accel_free(acc)
      !Frees the accelerator object acc
      use global
      integer(POINTER_SIZE), intent(inout) :: acc
    end subroutine gsl_interp_accel_free

    elemental function gsl_spline_eval(x, spl, acc)
      !Returns the interpolated value for a given point x
      use global
      implicit none
      real(R8), intent(in) :: x
      integer(POINTER_SIZE), intent(in) :: spl, acc
      real(R8) :: gsl_spline_eval
    end function gsl_spline_eval

    elemental function gsl_spline_eval_deriv(x, spl, acc)
      !Returns the derivative of an interpolated function for a given point x
      use global
      implicit none
      real(R8), intent(in) :: x
      integer(POINTER_SIZE), intent(in) :: spl, acc
      real(R8) :: gsl_spline_eval_deriv
    end function gsl_spline_eval_deriv

    elemental function gsl_spline_eval_deriv2(x, spl, acc)
      !Returns the second derivative d2 of an interpolated function for a given 
      !point x
      use global
      implicit none
      real(R8), intent(in) :: x
      integer(POINTER_SIZE), intent(in) :: spl, acc
      real(R8) :: gsl_spline_eval_deriv2
    end function gsl_spline_eval_deriv2

    function gsl_spline_eval_integ(a, b, spl, acc)
      !Returns the numerical integral of an interpolated function over the range
      ![a, b]
      use global
      implicit none
      real(R8), intent(in) :: a, b
      integer(POINTER_SIZE), intent(in) :: spl, acc
      real(R8) :: gsl_spline_eval_integ
    end function gsl_spline_eval_integ
  end interface


  !Multidimensional Root-Finding
  interface
    subroutine gsl_multiroot_fsolver_alloc(s, t, n)
      !Returns a pointer to a newly allocated instance of a solver of type t for
      !a system of n dimensions
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: s
      integer, intent(in) :: t, n
    end subroutine gsl_multiroot_fsolver_alloc

    subroutine gsl_multiroot_fsolver_free(s)
      !Frees all the memory associated with the solver s
      use global
      implicit none
      integer(POINTER_SIZE), intent(inout) :: s
    end subroutine gsl_multiroot_fsolver_free
  end interface

end module gsl_interface

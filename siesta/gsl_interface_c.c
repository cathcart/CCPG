/*
 Copyright (C) 2004-2007 M. Oliveira, F. Nogueira

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.
*/

#include <config.h>
#include <string_f.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multiroots.h>


/* This file contains interface routines that allow Fortran 90/95
routines to call GSL routines. The names of these routines are
compiler dependent: Fortran 90/95 compilers have different
name-mangling schemes.  For a description of the routines please refer
to the GSL documentation. */


          /* Error Handling */

void FC_FUNC_(gsl_strerror, GSL_STRERROR)
     (int *err, STR_F_TYPE res STR_ARG1)
{
  char *c;

  c = gsl_strerror(*err);
  TO_F_STR1(c, res);
}


          /* Mathematical Functions */

double FC_FUNC_(gsl_asinh, GSL_ASINH)
     (const double *x) 
{
  return gsl_asinh(*x);
}


          /* Special Functions */

double FC_FUNC_(gsl_sf_bessel_knu_scaled, GSL_SF_BESSEL_KNU_SCALED)
     (double *nu, double *x) 
{
  return gsl_sf_bessel_Knu_scaled(*nu,*x);
}


          /* Vectors and Matrices */

void FC_FUNC_(gsl_vector_alloc, GSL_VECTOR_ALLOC)
     (int *n, void **v)
{
  *v = (void*) gsl_vector_alloc (*n);
}

void FC_FUNC_(gsl_vector_free, GSL_VECTOR_FREE)
     (void **v)
{
  gsl_vector_free ((gsl_vector *)(*v));
}

double FC_FUNC_(gsl_vector_get, GSL_VECTOR_GET) (void **v, int *i)
{
  return gsl_vector_get ((gsl_vector *)(*v), *i);
}

void FC_FUNC_(gsl_vector_set, GSL_VECTOR_SET)
     (void **v, int *i, double *x)
{
  gsl_vector_set ((gsl_vector *)(*v), *i, *x);
} 

void FC_FUNC_(gsl_matrix_alloc, GSL_MATRIX_ALLOC)
     (int *n1, int *n2, void **m)
{
  *m = (void*) gsl_matrix_alloc (*n1, *n2);
}

void FC_FUNC_(gsl_matrix_calloc, GSL_MATRIX_CALLOC)
     (int *n1, int *n2, void **m)
{
  *m = (void*) gsl_matrix_calloc (*n1, *n2);
}

void FC_FUNC_(gsl_matrix_free, GSL_MATRIX_FREE)
     (void **m)
{
  gsl_matrix_free ((gsl_matrix *)(*m));
}

double FC_FUNC_(gsl_matrix_get, GSL_MATRIX_GET)
     (void **m, int *i, int *j)
{
  return gsl_matrix_get ((gsl_matrix *)(*m), *i, *j);
}

void FC_FUNC_(gsl_matrix_set, GSL_MATRIX_SET)
     (void **m, int *i, int *j, double *x)
{
  gsl_matrix_set ((gsl_matrix *)(*m), *i, *j, *x);
}


          /* Permutations */


void FC_FUNC_(gsl_permutation_alloc, GSL_PERMUTATION_ALLOC)
     (int *n, void **p)
{
  *p = (void*) gsl_permutation_alloc (*n);
}

void FC_FUNC_(gsl_permutation_free, GSL_PERMUTATION_FREE)
     (void **p)
{
  gsl_permutation_free ((gsl_permutation *)(*p));
}


          /* Linear Algebra */


int FC_FUNC_(gsl_linalg_lu_decomp, GSL_LINALG_LU_DECOMP)
     (void **m, void **p, int *signum)
{
  int s,ierr;
  ierr = gsl_linalg_LU_decomp ((gsl_matrix *)(*m) , (gsl_permutation *)(*p), &s);
  *signum = s;
  return (ierr);
}

int FC_FUNC_(gsl_linalg_lu_invert, GSL_LINALG_LU_INVERT)
     (void **LU, void **p, void **inverse)
{

  return gsl_linalg_LU_invert ((gsl_matrix *)(*LU),(gsl_permutation *)(*p), (gsl_matrix *)(*inverse)); 
}

int FC_FUNC_(gsl_linalg_lu_solve, GSL_LINALG_LU_SOLVE)
     (void **LU, void **p, void **b, void **x)
{
  return gsl_linalg_LU_solve ((gsl_matrix *)(*LU), (gsl_permutation *)(*p), (gsl_vector *)(*b), (gsl_vector *)(*x));
}


          /* Ordinary Differential Equations */

void FC_FUNC_(gsl_odeiv_step_alloc, GSL_ODEIV_STEP_ALLOC)
     (const int *stepping_func, const int *dim, void **stp)
{
  /* WARNING: This routine has a small modification: the user needs to
   pass an integer in order to choose the stepping function instead of
   a variable of type gsl_odeiv_step_type. */

  switch(*stepping_func) {
  case 1 : // Embedded 2nd order Runge-Kutta with 3rd order error estimate
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rk2, *dim);
    break;
  case 2 : // 4th order (classical) Runge-Kutta
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rk4, *dim);
    break;
  case 3 : // Embedded 4th order Runge-Kutta-Fehlberg method with 5th order error estimate. This method is a good general-purpose integrator
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, *dim);
    break;
  case 4 : // Embedded 4th order Runge-Kutta Cash-Karp method with 5th order error estimate
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rkck, *dim);
    break;
  case 5 : // Embedded 8th order Runge-Kutta Prince-Dormand method with 9th order error estimate
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, *dim);
    break;
  case 6 : // Implicit 2nd order Runge-Kutta at Gaussian points 
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rk2imp, *dim);
    break;
  case 7 : // Implicit 4th order Runge-Kutta at Gaussian points
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_rk4imp, *dim);
    break;
  case 8 : // Implicit Bulirsch-Stoer method of Bader and Deuflhard. This algorithm requires the Jacobian
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_bsimp, *dim);
    break;
  case 9 : // M=1 implicit Gear method 
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_gear1, *dim);
    break;
  case 10 : // M=2 implicit Gear method
    *stp = (void *)gsl_odeiv_step_alloc(gsl_odeiv_step_gear2, *dim);
    break;
  }
}

void FC_FUNC_(gsl_odeiv_step_free, GSL_ODEIV_STEP_FREE)
     (void **stp)
{
  gsl_odeiv_step_free((gsl_odeiv_step *)(*stp));
}

void FC_FUNC_(gsl_odeiv_evolve_alloc, GSL_ODEIV_EVOLVE_ALLOC)
     (const int *dim, void **evl)
{
  *evl = (void *)gsl_odeiv_evolve_alloc (*dim);
}

int FC_FUNC_(gsl_odeiv_evolve_reset, GSL_ODEIV_EVOLVE_RESET)
     (void **evl)
{
  return gsl_odeiv_evolve_reset ((gsl_odeiv_evolve *)(*evl));
}

void FC_FUNC_(gsl_odeiv_evolve_free, GSL_ODEIV_EVOLVE_FREE)
     (void **evl)
{
  gsl_odeiv_evolve_free((gsl_odeiv_evolve *)(*evl));
}

void FC_FUNC_(gsl_odeiv_control_standart_new, GSL_ODEIV_CONTROL_STANDART_NEW)
     (void **ctrl, double *epsabs, double *epsrel, double *ay, double *adydt)
{
  *ctrl = (void*) gsl_odeiv_control_standard_new (*epsabs, *epsrel, *ay, *adydt);
}

void FC_FUNC_(gsl_odeiv_control_free, GSL_ODEIV_CONTROL_FREE)
     (void **ctrl)
{
  gsl_odeiv_control_free((gsl_odeiv_control *)(*ctrl));
}


          /* Interpolation */

void FC_FUNC_(gsl_interp_accel_alloc, GSL_INTERP_ACCEL_ALLOC)
     (void **acc)
{
     *acc = (void *)gsl_interp_accel_alloc();
}

void FC_FUNC_(gsl_spline_alloc, GSL_SPLINE_ALLOC)
     (void **spl, int *n, short int *opt)
{
  /* WARNING: This routine has a small modification: the user needs to
     pass an integer in order to choose the interpolation type instead
     of a variable of type gsl_interp_type. */

  switch (*opt) {
  case 1 : // linear interpolation
    *spl = (void *)gsl_spline_alloc(gsl_interp_linear, *n);
    break;
  case 2 : // polynomial interpolation
    *spl = (void *)gsl_spline_alloc(gsl_interp_polynomial, *n);
    break;
  case 3 : // cubic spline with natural boundary conditions
    *spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *n);
    break;
  case 4 : // cubic spline with periodic boundary conditions 
    *spl = (void *)gsl_spline_alloc(gsl_interp_cspline_periodic, *n);
    break;
  case 5 : // Akima spline with natural boundary conditions
    *spl = (void *)gsl_spline_alloc(gsl_interp_akima, *n);
    break;
  case 6 : // Akima spline with periodic boundary conditions 
    *spl = (void *)gsl_spline_alloc(gsl_interp_akima_periodic, *n);
    break;
  }
}

void FC_FUNC_(gsl_spline_init, GSL_SPLINE_INIT)
     (void **spl, int *n, double *x, double *f)
{
  gsl_spline_init((gsl_spline *)(*spl), x, f, *n);
}

double FC_FUNC_(gsl_spline_eval, GSL_SPLINE_EVAL)
     (double *x, void **spl, void **acc)
{
  return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

double FC_FUNC_(gsl_spline_eval_deriv, GSL_SPLINE_EVAL_DERIV)
     (double *x, void **spl, void **acc)
{
  return gsl_spline_eval_deriv((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

double FC_FUNC_(gsl_spline_eval_deriv2, GSL_SPLINE_EVAL_DERIV2)
     (double *x, void **spl, void **acc)
{
  return gsl_spline_eval_deriv2((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

void FC_FUNC_(gsl_spline_free, GSL_SPLINE_FREE)
     (void **spl)
{
  gsl_spline_free((gsl_spline *)(*spl));
      
}

void FC_FUNC_(gsl_interp_accel_free, GSL_INTERP_ACCEL_FREE)
     (void **acc)
{
  gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

double FC_FUNC_(gsl_spline_eval_integ, GSL_SPLINE_EVAL_INTEG)
     (double *a, double *b, void **spl, void **acc)
{
  return gsl_spline_eval_integ((gsl_spline *)(*spl), *a, *b, (gsl_interp_accel *)(*acc));
}


          /* Multidimensional Root-Finding */

 void FC_FUNC_(gsl_multiroot_fsolver_alloc, GSL_MULTIROOT_FSOLVER_ALLOC)
     (void **s, const int *solver_type, const int *n)
{
  switch(*solver_type) {
  case 1 : // hybrid algorithm with internal scaling.
    *s = (void *)gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, *n);
    break;
  case 2 : // hybrid algorithm without internal scaling.
    *s = (void *)gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, *n);
    break;
  case 3 : // discrete Newton algorithm
    *s = (void *)gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, *n);
    break;
  case 4 : // Broyden algorithm
    *s = (void *)gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_broyden, *n);
    break;
  }
}

void FC_FUNC_(gsl_multiroot_fsolver_free, GSL_MULTIROOT_FSOLVER_FREE)
     (void **s)
{
  gsl_multiroot_fsolver_free((gsl_multiroot_fsolver *)(*s));
}

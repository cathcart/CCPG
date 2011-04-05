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

!#include "global.h"!i don't think this is necessary, looks like preprocessor stuff

module global
  implicit none

                    !---Names of kind types parameters---!

  !8 bytes real
  integer, parameter :: R8 = 8

                    !---Parameters---!

  real(R8), parameter :: M_ZERO     = 0.0_r8
  real(R8), parameter :: M_ONE      = 1.0_r8
  real(R8), parameter :: M_TWO      = 2.0_r8
  real(R8), parameter :: M_THREE    = 3.0_r8
  real(R8), parameter :: M_FOUR     = 4.0_r8
  real(R8), parameter :: M_FIVE     = 5.0_r8
  real(R8), parameter :: M_SIX      = 6.0_r8
  real(R8), parameter :: M_SEVEN    = 7.0_r8
  real(R8), parameter :: M_EIGHT    = 8.0_r8
  real(R8), parameter :: M_NINE     = 9.0_r8
  real(R8), parameter :: M_TEN      = 10.0_r8
  real(R8), parameter :: M_ELEVEN   = 11.0_r8
  real(R8), parameter :: M_TWELVE   = 12.0_r8
  real(R8), parameter :: M_FOURTEEN = 14.0_r8
  real(R8), parameter :: M_SIXTEEN  = 16.0_r8
  real(R8), parameter :: M_EIGHTEEN = 18.0_r8
  real(R8), parameter :: M_TWENTY   = 20.0_r8
  real(R8), parameter :: M_THIRTY   = 30.0_r8
  real(R8), parameter :: M_FIFTY    = 50.0_r8
  real(R8), parameter :: M_SIXTY    = 60.0_r8
  real(R8), parameter :: M_HUNDRED  = 100.0_r8
  real(R8), parameter :: M_HALF     = 0.5_r8
  real(R8), parameter :: M_THIRD    = M_ONE/M_THREE
  real(R8), parameter :: M_TWOTHIRD = M_TWO/M_THREE
  real(R8), parameter :: M_DIME     = 0.1_r8
  real(R8), parameter :: M_CENT     = 0.01_r8
  real(R8), parameter :: M_PI       = 3.141592653589793238462643383279502884197_r8
  real(R8), parameter :: M_C        = 137.03599976_r8
  real(R8), parameter :: M_C2       = 18778.86524_r8
  real(R8), parameter :: M_LOGC     = 4.92024366328042038029011278105_r8
  real(R8), parameter :: M_EPSILON  = 1.0e-20_r8

  !Other stuff
  logical :: in_debug_mode = .false.
  
  !something i put in here to cover the pointer size problem
  integer ,parameter :: POINTER_SIZE = 8

  !some extra crap to keep siesta_pass nice and clean
  real, parameter :: PACKAGE_VERSION = 1.0
  integer, parameter :: HAM   = 1, &
                        TM    = 2, &
                        RTM   = 3, &
                        MRPP  = 4

  integer, parameter :: AVERAGED = 1, &
                        J_DEP    = 2

  integer, parameter :: SCHRODINGER = 1, &
                        DIRAC       = 2, &
                        SCALAR_REL  = 3

  integer, parameter :: M_SIESTA  = 1, &
                        M_FHI     = 2, &
                        M_ABINIT5 = 3, &
                        M_ABINIT6 = 4, &
                        M_UPF     = 5, &
                        M_PARSEC  = 6
  !Mesh types
  integer, parameter :: LIN  = 1, &
                        LOG1 = 2, &
                        LOG2 = 3



end module global

!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
MODULE mod_lu
  implicit none

  interface LUDCMP
     module procedure LUDCMP_d, LUDCMP_c
  end interface LUDCMP
  interface LUBKSB
     module procedure LUBKSB_d, LUBKSB_c
  end interface LUBKSB


CONTAINS

  !  ***************************************************************
  !  * Given an N x N matrix A, this routine replaces it by the LU *
  !  * decomposition of a rowwise permutation of itself. A and N   *
  !  * are input. INDX is an output vector which records the row   *
  !  * permutation effected by the partial pivoting; D is output   *
  !  * as -1 or 1, depending on whether the number of row inter-   *
  !  * changes was even or odd, respectively. This routine is used *
  !  * in combination with LUBKSB to solve linear equations or to  *
  !  * invert a matrix. Return code is 1, if matrix is singular.   *
  !  ***************************************************************
  Subroutine LUDCMP_d(A,N,INDX,D,CODE)
    integer, parameter           :: NMAX=10
    double precision, parameter  :: TINY=1.5D-16
    double precision             :: AMAX,DUM, SUM, A(N,N),VV(NMAX)
    INTEGER                      :: CODE, D, INDX(N)
    integer                      :: N
    ! .. local ..
    integer                      :: I, J, K, IMAX

    D=1; CODE=0

    DO I=1,N
       AMAX=0.d0
       DO J=1,N
          IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
          CODE = 1
          RETURN
       END IF
       VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
       DO I=1,J-1
          SUM = A(I,J)
          DO K=1,I-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
          SUM = A(I,J)
          DO K=1,J-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
          DUM = VV(I)*DABS(SUM)
          IF(DUM.GE.AMAX) THEN
             IMAX = I
             AMAX = DUM
          END IF
       END DO ! i loop  

       IF(J.NE.IMAX) THEN
          DO K=1,N
             DUM = A(IMAX,K)
             A(IMAX,K) = A(J,K)
             A(J,K) = DUM
          END DO ! k loop
          D = -D
          VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

       IF(J.NE.N) THEN
          DUM = 1.d0 / A(J,J)
          DO I=J+1,N
             A(I,J) = A(I,J)*DUM
          END DO ! i loop
       END IF
    END DO ! j loop

    RETURN
  END subroutine LUDCMP_d
  !=============================================================================
  Subroutine LUDCMP_c(A,N,INDX,D,CODE)
    DOUBLE PRECISION, PARAMETER :: TINY=1.5D-16
    DOUBLE COMPLEX AMAX,DUM, SUM, A(N,N)
    INTEGER CODE, D, INDX(N)
    DOUBLE COMPLEX, pointer ::  VV(:)     !complex auxiliary vector (n)
    integer                      :: N
    ! .. local ..
    integer                      :: I, J, K, IMAX, ialloc

    allocate(VV(N),stat=ialloc)

    D=1; CODE=0

    DO I=1,N
       AMAX=CMPLX(0.0d0,0.0d0,kind(1.0d0))
       DO J=1,N
          IF (ABS(A(I,J)).GT.ABS(AMAX))  AMAX=(A(I,J))
       END DO ! j loop
       IF(ABS(AMAX).LT.TINY) THEN
          CODE = 1
          RETURN
       END IF
       VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
       DO I=1,J-1
          SUM = A(I,J)
          DO K=1,I-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
       END DO ! i loop
       AMAX = CMPLX(0.0d0,0.0d0,kind(1.0d0))
       DO I=J,N
          SUM = A(I,J)
          DO K=1,J-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
          DUM = VV(I)*ABS(SUM)
          IF(ABS(DUM).GE.ABS(AMAX)) THEN
             IMAX = I
             AMAX = DUM
          END IF
       END DO ! i loop  

       IF(J.NE.IMAX) THEN
          DO K=1,N
             DUM = A(IMAX,K)
             A(IMAX,K) = A(J,K)
             A(J,K) = DUM
          END DO ! k loop
          D = -D
          VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(ABS(A(J,J)) < TINY) A(J,J) = CMPLX(TINY,0.0d0,kind(1.0d0))

       IF(J.NE.N) THEN
          DUM = 1.0 / A(J,J)
          DO I=J+1,N
             A(I,J) = A(I,J)*DUM
          END DO ! i loop
       END IF
    END DO ! j loop

    deallocate(VV)

    RETURN
  END subroutine LUDCMP_c
  !=============================================================================

  !  ******************************************************************
  !  * Solves the set of N linear equations A . X = B.  Here A is     *
  !  * input, not as the matrix A but rather as its LU decomposition, *
  !  * determined by the routine LUDCMP. INDX is input as the permuta-*
  !  * tion vector returned by LUDCMP. B is input as the right-hand   *
  !  * side vector B, and returns with the solution vector X. A, N and*
  !  * INDX are not modified by this routine and can be used for suc- *
  !  * cessive calls with different right-hand sides. This routine is *
  !  * also efficient for plain matrix inversion.                     *
  !  ******************************************************************
  Subroutine LUBKSB_d(A,N,INDX,B)
    double precision        :: SUM, A(N,N),B(N)
    INTEGER                 :: INDX(N), N
    ! .. local ..
    integer                 :: II, LL, I, J

    II = 0

    DO I=1,N
       LL = INDX(I)
       SUM = B(LL)
       B(LL) = B(I)
       IF(II.NE.0) THEN
          DO J=II,I-1
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       ELSE IF(SUM.NE.0.d0) THEN
          II = I
       END IF
       B(I) = SUM
    END DO ! i loop

    DO I=N,1,-1
       SUM = B(I)
       IF(I < N) THEN
          DO J=I+1,N
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       END IF
       B(I) = SUM / A(I,I)
    END DO ! i loop

    RETURN
  END subroutine LUBKSB_d
  !=============================================================================
  Subroutine LUBKSB_c(A,N,INDX,B)
    DOUBLE COMPLEX SUM, A(N,N),B(N)
    INTEGER                 :: INDX(N), N
    ! .. local ..
    integer                 :: II, LL, I, J

    II = 0

    DO I=1,N
       LL = INDX(I)
       SUM = B(LL)
       B(LL) = B(I)
       IF(II.NE.0) THEN
          DO J=II,I-1
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       ELSE IF(ABS(SUM).NE.0.0) THEN
          II = I
       END IF
       B(I) = SUM
    END DO ! i loop

    DO I=N,1,-1
       SUM = B(I)
       IF(I < N) THEN
          DO J=I+1,N
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       END IF
       B(I) = SUM / A(I,I)
    END DO ! i loop

    RETURN
  END subroutine LUBKSB_c
  !=============================================================================
END MODULE mod_lu
!=============================================================================

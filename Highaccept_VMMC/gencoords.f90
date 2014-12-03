      PROGRAM GENERATE_COORDINATES
!     Generates random configuraions for an atomic or a rigid-body system
!     Set the parameters: NCONFIG, NMOL, RADIUS and RBAAT 
!     NCONFIG: Number of independent configurations required to be generated
!     NMOL: Number of molecules 
!     RADIUS: Radius of the spherical container holding the coordinates
!     RBAAT: Set TRUE for a rigid-body system
      IMPLICIT NONE 

      INTEGER, PARAMETER :: NCONFIG = 1, NMOL = 24
      DOUBLE PRECISION, PARAMETER :: RADIUS = 2.8D0 
      LOGICAL, PARAMETER :: RBAAT = .TRUE.
      INTEGER          :: I, J, K
      DOUBLE PRECISION ::   COORDS(NCONFIG,6*NMOL), RANF, DUMMY, Q(4), P(3), V(3)
      CHARACTER (LEN = 20) :: CI, FILENAME, ORTNFILE
      LOGICAL :: OVERLAP

      DO I = 1, NCONFIG

         WRITE(CI,'(I2)') I
         CI = ADJUSTL(CI)
         FILENAME = 'initialpos.dat'//TRIM(CI)
         OPEN (UNIT = 1, FILE = FILENAME, STATUS = 'UNKNOWN')

         J=1
         K=3*J
         CALL RANDVEC(V)
         COORDS(I,K-2:K) = RADIUS*RANF(DUMMY)*V(:)
        
         DO J = 2, NMOL
           OVERLAP=.FALSE.
           DO WHILE (OVERLAP.EQV..FALSE.)
             K = J*3 !*3 IN ORDER TO GENERATE X, Y, Z COORDINATES
             CALL RANDVEC (V)
             COORDS(I,K-2:K) = RADIUS*RANF(DUMMY)*V(:) !USE RANF FUNCTION TO GENERATE RANDOM COORDINATE WITHIN RADIUS SELECTED
             CALL OVERSUB(COORDS, OVERLAP, I, K)
           END DO
         END DO

         IF (.NOT. RBAAT) RETURN
         
         ORTNFILE = 'initialortn.dat' 
         OPEN(UNIT = 2, FILE = ORTNFILE, STATUS = 'UNKNOWN')
         DO J = 1, NMOL
            CALL RANDQTN(Q)
            !CALL RQTOAA(Q,P)
            !K = 3*NMOL + 3*J 
            !COORDS(I,K-2:K) = P(1:3)
            WRITE(2,*) Q(1), Q(2), Q(3), Q(4)
         ENDDO

         DO J = 1, NMOL
            K = 3*J
            WRITE(1,*) COORDS(I,K-2), COORDS(I,K-1), COORDS(I,K)
         ENDDO

      ENDDO
         
      END PROGRAM GENERATE_COORDINATES

!     ==============================================================================================

      SUBROUTINE RANDVEC(V)
!     uniformly random unit vector
!     See Marsaglia, Ann. Math. Stat. 43, 645 (1972); Vesely, J. Comput. Phys. 47, 291-296 (1982).

      IMPLICIT NONE
      DOUBLE PRECISION :: V(3), DUMMY, RANF, XISQ, XI1, XI2, XI

      XISQ = 1.0
!     iterative loop
      DO WHILE (XISQ >= 1.D0)
         XI1  = RANF(DUMMY)*2.D0 - 1.D0
         XI2  = RANF(DUMMY)*2.D0 - 1.D0
         XISQ = XI1*XI1 + XI2*XI2
      ENDDO

      XI = SQRT(1.D0 - XISQ)
      V(1) = 2.D0 * XI1*XI
      V(2) = 2.D0 * XI2*XI
      V(3) = 1.D0 - 2.D0*XISQ

      END SUBROUTINE RANDVEC

!     ==============================================================================================

      SUBROUTINE RANDQTN(Q)
!     uniformly random unit quaternion
!     See Marsaglia, Ann. Math. Stat. 43, 645 (1972); Vesely, J. Comput. Phys. 47, 291-296 (1982).

      IMPLICIT NONE
      DOUBLE PRECISION :: Q(4), DUMMY, RANF, FCTR, XIASQ, XIBSQ, XI1, XI2, XI3, XI4 

      XIASQ = 1.0
      XIBSQ = 1.0
!     iterative loop
      DO WHILE (XIASQ >= 1.D0)
         XI1  = RANF(DUMMY)*2.D0 - 1.D0
         XI2  = RANF(DUMMY)*2.D0 - 1.D0
         XIASQ = XI1*XI1 + XI2*XI2
      ENDDO
      DO WHILE (XIBSQ >= 1.D0)
         XI3  = RANF(DUMMY)*2.D0 - 1.D0
         XI4  = RANF(DUMMY)*2.D0 - 1.D0
         XIBSQ = XI3*XI3 + XI4*XI4
      ENDDO
      FCTR = SQRT((1.D0-XIASQ)/XIBSQ)
      Q(1) = XI1
      Q(2) = XI2
      Q(3) = XI3*FCTR
      Q(4) = XI4*FCTR

      END SUBROUTINE RANDQTN

!     ==============================================================================================

      SUBROUTINE RQTOAA(Q,P)

!     transforms a unit quaternion Q to the corresponding angle-axis variables P  

      IMPLICIT NONE
      DOUBLE PRECISION :: Q(4), P(3), THETA, FCT

      THETA  = 2.D0*ACOS(Q(1))

      IF (THETA <= 1.D-12) THEN
         P (1:3) = 0.D0
      ELSE
         FCT = DSQRT(DOT_PRODUCT(Q(2:4),Q(2:4)))
         P(1:3) = THETA*Q(2:4)/FCT
      ENDIF

      END SUBROUTINE RQTOAA

!     ==============================================================================================

      DOUBLE PRECISION FUNCTION RANF(DUMMY)
!     drawing a uniform random variate between 0 and 1

      INTEGER, PARAMETER :: L = 1029, C = 221591, M = 1048576
      INTEGER ::          SEED
      DOUBLE PRECISION :: DUMMY
      SAVE             SEED
      DATA             SEED / 0 /

      SEED = MOD(SEED * L + C, M)
      RANF = DFLOAT(SEED) / DFLOAT(M)

      END FUNCTION RANF

! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE OVERSUB(COORDS, OVERLAP, I, K)

      IMPLICIT NONE

      INTEGER :: I, K, X, Y
      DOUBLE PRECISION :: COORDS(I,K), NEWXYZ(3), RIJ(3), DIS, MINDISVAL
      LOGICAL :: OVERLAP

      MINDISVAL=0.9D0
      X= (K/3)
      WRITE(*,*) X
      NEWXYZ(1)=COORDS(I,(K-2))
      NEWXYZ(2)=COORDS(I,(K-1))
      NEWXYZ(3)=COORDS(I,K)

      DO Y =1, X-1
        RIJ(1)=ABS(NEWXYZ(1)-COORDS(I,(3*Y-2)))
        RIJ(2)=ABS(NEWXYZ(2)-COORDS(I,(3*Y-1)))
        RIJ(3)=ABS(NEWXYZ(3)-COORDS(I,(3*Y)))

        DIS=SQRT(DBLE(RIJ(1)**2)+DBLE(RIJ(2)**2)+DBLE(RIJ(3)**2))
        WRITE(*,*) DIS
        IF (DIS.GE.MINDISVAL) THEN
          OVERLAP=.TRUE.
        ELSE
          OVERLAP=.FALSE.
        END IF

        IF (OVERLAP.eqv..FALSE.) EXIT

      END DO

      END SUBROUTINE OVERSUB

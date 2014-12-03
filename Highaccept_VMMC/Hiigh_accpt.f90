!     ===============================================================================================

SUBROUTINE TRANS_VMMC(PE, ACCPTC)

USE COMMONS

IMPLICIT NONE

INTEGER :: I, J, K, L, X, Z, Y, V, W, T, INDXP, CACCPT, OCOUNT, FCOUNT, BUCOUNT, BVCOUNT, ACCPTC
INTEGER, ALLOCATABLE :: FFAIL(:), OFAIL(:,:), C(:), PROPD(:,:), BU(:,:), BV(:,:), RN(:,:)
DOUBLE PRECISION :: RIJ(3), RIJSQ, DIS, UPE, UIPE, UJPE, M(3), RANF, DUMMY, PE, PET
DOUBLE PRECISION, PARAMETER :: DCUT = 1.5D0
LOGICAL :: GROWTH, OUTEST, FTEST, ALPROP, FRUSTEST, PRETEST, ACCPTCLUS, INSIDET, PERCT

ALLOCATE(PROPD(NPART,NPART), C(NPART), OFAIL(NPART,NPART), FFAIL(NPART), RN(NPART,3))

DO X = 1, NPART

  DO I = 1, NPART
    DO J = 1, NPART
      PROPD(I,J) = 0
      OFAIL(I,J) = 0 ! initialises array storing outright failed links
    END DO
  END DO

  DO I = 1, NPART
    C(I) = 0 ! intialises array storing members of the cluster
    FFAIL(I) = 0 ! initialise array storing frustarted failed links
  END DO



  CACCPT = 1 ! number of particles accepted to the group

  GROWTH = .TRUE. ! logical that determines wether the group is still growing or

  INDXP = INT(NPART*RANF(DUMMY)) + 1 ! select random root partcle

  C(1) = INDXP ! allocates the randomly selected particle as part of the group

  DO I = 1, 3
    M(I) = (2.D0*RANF(DUMMY) - 1.D0)*MAXDTR ! random translational displacement
  END DO

  OCOUNT = 0 ! a counter to determine the number of outright failed links

  FCOUNT = 0 ! a counter to determine the number of frustrated failed links

  Z = 0 ! Z is a counter to help determine when the cluster has finished growing

! recursive pairwise selection of particle in initial state u

  DO WHILE (GROWTH)

    Z = Z + 1

    I = C(Z)
!    WRITE(*,*) Z
    DO J = 1, NPART

      ALPROP = .FALSE. ! logical to fill conditions that for pair (i,j) j is not already a memeber of the group and a link between (i,j) has not already been proposed

      OUTEST = .FALSE.

      FTEST = .FALSE.

      IF (I==J) CYCLE

      DO K = 1, CACCPT ! do loop checks to make sure j is not already a member of the group
        L = C(K)
        IF (L==J) THEN
          ALPROP = .TRUE.
        END IF
      END DO

      IF (ALPROP.EQV..TRUE.) CYCLE

      IF (PROPD(J,I)==1) THEN ! if statement checks to see that a link has not already been proposed to pair (i,j)
        ALPROP = .TRUE.
      ENDIF

      IF (ALPROP.EQV..TRUE.) CYCLE

! distance criterion follows this comment

      RIJ(:) = R(J,:) - R(I,:)
      RIJSQ = DOT_PRODUCT(RIJ(:),RIJ(:))
      DIS = SQRT(RIJSQ)
      IF (DIS > DCUT) CYCLE ! if the distance is larger than the cut off cycle
      PROPD(I,J) = 1 ! if distance is less than dcut it records that a link is the proposed between (i,j)

! calculation of enrgies and link formation probabilities follows this comment

      CALL IJENERGY(I, J, UPE) ! calculates energy of (i,j) in initial state u

      CALL MIJENERGY(I, J, M, UIPE) ! calculates energy of (i,j) after applying M to particle i

      CALL LINKO(I, J, UPE, UIPE, OCOUNT,  OFAIL, OUTEST) ! link formation test 1

      IF (OUTEST.EQV..TRUE.) CYCLE

      CALL IMJENERGY(I, J, M, UJPE) ! calculates energy of (i,j) after applying M to j

      CALL LINKF(I, J, UPE, UIPE, UJPE, FCOUNT, FFAIL, FTEST) ! link formation test 2

      IF (FTEST.EQV..TRUE.) CYCLE

! accept new group member

      CACCPT = CACCPT + 1
      C(CACCPT) = J
!      WRITE(*,*) CACCPT
    END DO

! check to see if cluster selection is finished

    IF (Z<CACCPT) THEN
      GROWTH = .TRUE.
    ELSE
      GROWTH = .FALSE.
    END IF

  END DO

! recursive pairwise selection of particles in state ui

  Z = 0

  DO WHILE (GROWTH)

    Z = Z + 1

    I = C(Z)
!    WRITE(*,*) Z
    DO J = 1, NPART

      ALPROP = .FALSE. ! logical to fill conditions that for pair (i,j) j is not already a memeber of the group and a link between (i,j) has not already been proposed

      OUTEST = .FALSE.

      FTEST = .FALSE.

      IF (I==J) CYCLE

      DO K = 1, CACCPT ! do loop checks to make sure j is not already a member of the group
        L = C(K)
        IF (L==J) THEN
          ALPROP = .TRUE.
        END IF
      END DO

      IF (ALPROP.EQV..TRUE.) CYCLE

      IF (PROPD(J,I)==1) THEN ! if statement checks to see that a link has not already been proposed to pair (i,j)
        ALPROP = .TRUE.
      ENDIF

      IF (ALPROP.EQV..TRUE.) CYCLE

! distance criterion follows this comment

      DO K = I, NPART
        RN(:,:) = R(:,:)
      END DO

      RN(I,:) = RN(I,:) + M(:)

      RIJ(:) = RN(J,:) - RN(I,:)
      RIJSQ = DOT_PRODUCT(RIJ(:),RIJ(:))
      DIS = SQRT(RIJSQ)
      IF (DIS > DCUT) CYCLE ! if the distance is larger than the cut off cycle
      PROPD(I,J) = 1 ! if distance is less than dcut it records that a link is the proposed between (i,j)

! calculation of enrgies and link formation probabilities follows this comment

      CALL IJENERGY(I, J, UPE) ! calculates energy of (i,j) in initial state u

      CALL MIJENERGY(I, J, M, UIPE) ! calculates energy of (i,j) after applying M to particle i

      CALL LINKO(I, J, UPE, UIPE, OCOUNT,  OFAIL, OUTEST) ! link formation test 1

      IF (OUTEST.EQV..TRUE.) CYCLE

      CALL IMJENERGY(I, J, M, UJPE) ! calculates energy of (i,j) after applying M to j

      CALL LINKF(I, J, UPE, UIPE, UJPE, FCOUNT, FFAIL, FTEST) ! link formation test 2

      IF (FTEST.EQV..TRUE.) CYCLE

! accept new group member

      CACCPT = CACCPT + 1
      C(CACCPT) = J
!      WRITE(*,*) CACCPT
    END DO

! check to see if cluster selection is finished

    IF (Z<CACCPT) THEN
      GROWTH = .TRUE.
    ELSE
      GROWTH = .FALSE.
    END IF

  END DO

! TEST for the existence of frustrated links in the cluster boundary

  FRUSTEST = .FALSE.

  IF (FCOUNT>0) THEN
    DO T = 1, FCOUNT
      Y = FFAIL(FCOUNT)
      PRETEST = .FALSE.
      DO V = 1, CACCPT
        W = C(V)
        IF (Y==W) THEN
          PRETEST = .TRUE.
        END IF
      END DO
      IF (PRETEST.EQV..FALSE.) THEN
        FRUSTEST = .TRUE.
      END IF
    END DO
  END IF

! if frustrated links present reject the cluster move
  WRITE(*,*) FRUSTEST, 'frustest'
  IF (FRUSTEST.EQV..TRUE.) CYCLE

   
  DO I = 1, CACCPT
    J = C(I)
    R(J,:) = R(J,:) + M(:)
    RS(J,:) = RS(J,:) + M(:)
  END DO

  CALL CONTAINER(INSIDET)

  IF (INSIDET.EQV..FALSE.) THEN
    DO I = 1, CACCPT
      J = C(I)
      R(J,:) = R(J,:) - M(:)
      RS(J,:) = RS(J,:) - M(:)
    END DO
  END IF

  IF (INSIDET.EQV..FALSE.) CYCLE

  ACCPTC = ACCPTC + 1
  WRITE(*,*) ACCPTC
  CALL ENRG (PET)

  PE = PET

END DO

END SUBROUTINE

!==================================================================================================================

SUBROUTINE IJENERGY(I, J, UPE)

USE COMMONS

IMPLICIT NONE

INTEGER ::  I, J
DOUBLE PRECISION :: RIJ(3), RIJSQ, ABSRIJ, R2, R4, EXPFCT, RSS(3), RSSSQ
DOUBLE PRECISION :: UPE, ALP, BET, GAM, VR, VA, VB, VG, DPFCT, NR(3), EI(3), EJ(3)

UPE  = 0.D0

EI(:) = E(I,:)
EJ(:) = E(J,:)

RIJ(:) = R(J,:) - R(I,:)
RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))

! Isotropic interaction between the spherical cores described by the Yukawa potential

ABSRIJ = DSQRT(RIJSQ)
R2     = 1.D0/RIJSQ
EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
UPE = UPE + EXPFCT/ABSRIJ

! Dipolar contribution

EJ(:) = E(J,:)
RSS(:) = RS(J,:) - RS(I,:)
R2     = DOT_PRODUCT(RSS(:),RSS(:))
ABSRIJ = DSQRT(R2)
NR(:)  = RSS(:)/ABSRIJ
R2     = 1.D0/R2
R4     = R2*R2
ALP    = DOT_PRODUCT(NR(:),EJ(:))
BET    = DOT_PRODUCT(NR(:),EI(:))
GAM    = DOT_PRODUCT(EJ(:),EI(:))
DPFCT  = 3.D0*DPMU*DPMU
UPE = UPE + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

IF (FIELDT) THEN
  UPE = UPE - DPMU*FIELD*E(I,3)
ENDIF

!WRITE(*,*) UPE

END SUBROUTINE IJENERGY

! =======================================================================================

SUBROUTINE MIJENERGY(I, J, M, UIPE)

USE COMMONS

IMPLICIT NONE

INTEGER :: I, J, X
DOUBLE PRECISION :: RIJ(3), RIJSQ, ABSRIJ, R2, R4, EXPFCT, RSS(3), RSSSQ
DOUBLE PRECISION :: UIPE, ALP, BET, GAM, VR, VA, VB, VG, DPFCT, NR(3), EI(3), EJ(3), M(3)
DOUBLE PRECISION, ALLOCATABLE :: RN(:,:), RNS(:,:)

ALLOCATE(RN(NPART,3), RNS(NPART,3))

DO X = 1, NPART
  RN(X,:) = R(X,:)
  RNS(X,:) = RS(X,:)
END DO

RN(I,:) = RN(I,:) + M(:)
RNS(I,:) = RNS(I,:) + M(:)

UIPE  = 0.D0

EI(:) = E(I,:)
EJ(:) = E(J,:)

RIJ(:) = RN(J,:) - RN(I,:)
RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))

! Isotropic interaction between the spherical cores described by the Yukawa potential

ABSRIJ = DSQRT(RIJSQ)
R2     = 1.D0/RIJSQ
EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
UIPE = UIPE + EXPFCT/ABSRIJ

! Dipolar contribution

EJ(:) = E(J,:)
RSS(:) = RNS(J,:) - RNS(I,:)
R2     = DOT_PRODUCT(RSS(:),RSS(:))
ABSRIJ = DSQRT(R2)
NR(:)  = RSS(:)/ABSRIJ
R2     = 1.D0/R2
R4     = R2*R2
ALP    = DOT_PRODUCT(NR(:),EJ(:))
BET    = DOT_PRODUCT(NR(:),EI(:))
GAM    = DOT_PRODUCT(EJ(:),EI(:))
DPFCT  = 3.D0*DPMU*DPMU
UIPE = UIPE + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

IF (FIELDT) THEN
  UIPE = UIPE - DPMU*FIELD*E(I,3)
ENDIF

! WRITE(*,*) UIPE

END SUBROUTINE MIJENERGY

! ============================================================================================

SUBROUTINE LINKO(I, J, UPE, UIPE, OCOUNT,  OFAIL, OUTEST)

USE COMMONS, ONLY : TMP, NPART

IMPLICIT NONE

INTEGER :: I, J, OCOUNT, OFAIL(NPART,NPART)
DOUBLE PRECISION :: UPE, UIPE, IDIFF, LOPROB, RANF, DUMMY
LOGICAL :: OUTEST

IDIFF = UIPE - UPE
LOPROB = 1 - (DEXP(- IDIFF / TMP))

IF (LOPROB.LT.0.D0) THEN
  LOPROB = 0.D0
END IF

IF ((RANF(DUMMY)>LOPROB)) THEN
  OUTEST = .TRUE.
  OCOUNT = OCOUNT + 1
  OFAIL(I, J) = 1
ELSE
  OUTEST = .FALSE.
ENDIF

!WRITE(*,*) IDIFF, OUTEST

END SUBROUTINE

! ============================================================================================

SUBROUTINE IMJENERGY(I, J, M, UJPE)

USE COMMONS

IMPLICIT NONE

INTEGER :: I, J, X
DOUBLE PRECISION :: RIJ(3), RIJSQ, ABSRIJ, R2, R4, EXPFCT, RSS(3), RSSSQ
DOUBLE PRECISION :: UJPE, ALP, BET, GAM, VR, VA, VB, VG, DPFCT, NR(3), EI(3), EJ(3), M(3)
DOUBLE PRECISION, ALLOCATABLE :: RN(:,:), RNS(:,:)

ALLOCATE(RN(NPART,3), RNS(NPART,3))

DO X = 1, NPART
  RN(X,:) = R(X,:)
  RNS(X,:) = RS(X,:)
END DO

RN(J,:) = RN(J,:) + M(:)
RNS(J,:) = RNS(J,:) + M(:)

UJPE  = 0.D0

EI(:) = E(I,:)
EJ(:) = E(J,:)

RIJ(:) = RN(J,:) - RN(I,:)
RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))

! Isotropic interaction between the spherical cores described by the Yukawa potential

ABSRIJ = DSQRT(RIJSQ)
R2     = 1.D0/RIJSQ
EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
UJPE = UJPE + EXPFCT/ABSRIJ

! Dipolar contribution

EJ(:) = E(J,:)
RSS(:) = RNS(J,:) - RNS(I,:)
R2     = DOT_PRODUCT(RSS(:),RSS(:))
ABSRIJ = DSQRT(R2)
NR(:)  = RSS(:)/ABSRIJ
R2     = 1.D0/R2
R4     = R2*R2
ALP    = DOT_PRODUCT(NR(:),EJ(:))
BET    = DOT_PRODUCT(NR(:),EI(:))
GAM    = DOT_PRODUCT(EJ(:),EI(:))
DPFCT  = 3.D0*DPMU*DPMU
UJPE = UJPE + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ


IF (FIELDT) THEN
  UJPE = UJPE - DPMU*FIELD*E(I,3)
ENDIF

END SUBROUTINE IMJENERGY

! =========================================================================================================

SUBROUTINE LINKF(I, J, UPE, UIPE, UJPE, FCOUNT, FFAIL, FTEST)

USE COMMONS, ONLY : TMP, NPART

IMPLICIT NONE

INTEGER :: I, J, FCOUNT, FFAIL(NPART)
DOUBLE PRECISION :: UPE, UIPE, UJPE, IDIFF, JDIFF, FPI, FPJ, LFPROB, RANF, DUMMY
LOGICAL :: FTEST

IDIFF = UIPE - UPE
JDIFF = UJPE - UPE

FPI = 1 - (DEXP(-IDIFF/TMP))
FPJ = 1 - (DEXP(-JDIFF/TMP))

LFPROB = FPJ/FPI

IF (LFPROB.GT.1.D0) THEN
  LFPROB = 1.D0
END IF

IF ((RANF(DUMMY)>LFPROB)) THEN
  FTEST = .TRUE.
  FCOUNT = FCOUNT + 1
  FFAIL(FCOUNT) = J
ELSE IF (LFPROB==1.D0) THEN
  FTEST = .FALSE.
ELSE
  FTEST = .FALSE.
END IF

IF (FTEST.EQV..FALSE..AND.IDIFF.GT.100.D0) THEN
  FTEST = .TRUE.
  FCOUNT = FCOUNT + 1
  FFAIL(FCOUNT) = J
END IF


END SUBROUTINE

! ========================================================================================================


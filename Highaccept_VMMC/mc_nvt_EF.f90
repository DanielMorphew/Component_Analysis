      MODULE COMMONS

      IMPLICIT NONE

      INTEGER :: NPART, NRBSITE, NTSITE, MYUNIT, ISTEP, QCOUNT 
      DOUBLE PRECISION :: YKAPPA, DPMU, RSHIFT, FIELD
      DOUBLE PRECISION :: MAXDTR, MAXDRT, TMP, RADIUS
      DOUBLE PRECISION, ALLOCATABLE :: RBSITE(:,:), RBUV(:,:), R(:,:), Q(:,:), RS(:,:), E(:,:)
      LOGICAL          :: FIELDT

      END MODULE COMMONS

!     ==============================================================================================

      PROGRAM SHIFTDP_MC_NVT

      USE COMMONS

      IMPLICIT NONE
             
 
      INTEGER          :: I, K, ICNT, ACCPTC, NSTART, J1, J2, J3
      INTEGER          :: NSTEP, NEQ, IDUMP1, IDUMP2, IDUMP3
      DOUBLE PRECISION :: VLM, W, PRS, SUMPRS, AVPRS, PRSFIX, PET
      DOUBLE PRECISION :: PE, PEPP, SUMPE, AVPE, PEEND, QI(4), RM(3,3)
      DOUBLE PRECISION :: SUMVLM, SUMDNS, AVVLM, AVDNS, DNS, RACCPC

      OPEN (UNIT = 1, FILE = 'parameter_mc.inp', STATUS = 'UNKNOWN')

!     Read input parameters

      READ (1, *)
      READ (1, *) NPART, NRBSITE 
      READ (1, *) 
      READ (1, *) YKAPPA, RSHIFT, DPMU, FIELD 
      READ (1, *)
      READ (1, *) TMP
      READ (1, *)
      READ (1, *) RADIUS
      READ (1, *)
      READ (1, *) MAXDTR, MAXDRT
      READ (1, *)
      READ (1, *) NSTEP, NEQ
      READ (1, *) 
      READ (1, *) IDUMP1, IDUMP2, IDUMP3
 
      CLOSE (UNIT = 1)

      QCOUNT = 0
      MYUNIT = 21
      OPEN (MYUNIT, FILE = 'run.dat', STATUS = 'UNKNOWN')

      NTSITE = NPART*NRBSITE

      IF (FIELD /= 0.D0) THEN
         FIELDT = .TRUE.
      ELSE
         FIELDT = .FALSE.
      ENDIF

      ALLOCATE(RBSITE(NRBSITE,3), RBUV(NRBSITE,3), R(NPART,3), Q(NPART,4), RS(NTSITE,3), E(NTSITE,3))

      CALL DEFSHIFTDP()
     
      CALL EQCON()

!     Start simulation

      ISTEP  = 0  
      ICNT   = 0
      SUMPE  = 0.D0
!      SUMPRS = 0.D0
!      SUMVLM = 0.D0
!      SUMDNS = 0.D0
      ACCPTC = 0
        
!      CALL ENRG(PE,W)
      CALL ENRG(PE)
      
      DO WHILE (ISTEP < NSTEP)            

!     Sample the configuration space
           
!         CALL MOVE(PE, ACCPTC)

         CALL TRANS_VMMC(PE, ACCPTC)

         CALL CENTRE()

         DO J1 = 1, NPART
            QI(:) = Q(J1,:)

            CALL ROTMAT(QI, RM)

            DO J2 = 1, NRBSITE
               J3 = NRBSITE*(J1-1) + J2
               RS(J3,:) = R(J1,:) + MATMUL(RM,RBSITE(J2,:))
!               E(J3,:)  = MATMUL(RM(:,:),RBUV(J2,:))
            ENDDO
         ENDDO


!         CALL MVVLM (PE, W, PRSFIX, VLM, BMAX, ACCPTV)

!     Calculate average values 

         PEPP  = PE / DFLOAT(NPART)
         SUMPE = SUMPE + PEPP
!         PRS     = (NP * TMP + W / 3.D0) / VLM
!         SUMPRS  = SUMPRS + PRS
!         SUMVLM  = SUMVLM + VLM
!         DNS     = DFLOAT(NPART)/ VLM
!         SUMDNS  = SUMDNS + DNS
         ICNT    = ICNT + 1
         ISTEP   = ISTEP + 1
 
         IF (MOD(ISTEP, IDUMP1) == 0) THEN

            AVPE   = SUMPE / DFLOAT(ICNT)
!            print *, PEPP, SUMPE, ICNT, AVPE
!            AVPRS  = SUMPRS / DFLOAT(ICNT)
!            AVVLM  = SUMVLM / DFLOAT(ICNT)
!            AVDNS  = SUMDNS / DFLOAT(ICNT)

            OPEN (UNIT=3, FILE='energy.dat', STATUS='UNKNOWN',ACCESS ='APPEND')  
!            WRITE (3, *) PEPP, AVPE, PRS, AVPRS
            WRITE (3, *) PEPP, AVPE
            CLOSE (UNIT = 3, STATUS = 'KEEP')

!            OPEN (UNIT=9, FILE='dns.dat', STATUS='UNKNOWN', ACCESS ='APPEND')
!            WRITE (9, *) BOXL(1), BOXL(2), BOXL(3), AVVLM, AVDNS
!            CLOSE (UNIT = 9, STATUS = 'KEEP')

         ENDIF

         IF (ISTEP > NEQ .AND. MOD(ISTEP, IDUMP2) == 0) THEN   
            OPEN (UNIT=7, FILE='pos.dat',  STATUS='UNKNOWN', ACCESS='APPEND')
            OPEN (UNIT=8, FILE='ortn.dat', STATUS='UNKNOWN', ACCESS='APPEND')
 
            DO I = 1, NPART
               WRITE (7, *) R(I,1), R(I,2), R(I,3)
               WRITE (8, *) Q(I,1), Q(I,2), Q(I,3), Q(I,4)
            ENDDO

            CLOSE (UNIT=7, STATUS='KEEP')  
            CLOSE (UNIT=8, STATUS='KEEP')  

!            CALL QUENCH(PE)

         ENDIF
               
         IF (ISTEP == NEQ) THEN
            ICNT   = 0
            SUMPE  = 0.D0    
            SUMPRS = 0.D0
            SUMVLM = 0.D0
            SUMDNS = 0.D0
         ENDIF  

         IF (MOD(ISTEP, IDUMP3) == 0) THEN
            RACCPC = DFLOAT(ACCPTC)/DFLOAT(IDUMP3*NPART)
!            IF (RACCPC < 0.45D0) THEN
!               MAXDTR = MAXDTR * 0.975D0
!               MAXDRT = MAXDRT * 0.975D0
!            ELSE
!               MAXDTR = MAXDTR * 1.025D0
!               MAXDRT = MAXDRT * 1.025D0
!            ENDIF  
            WRITE (MYUNIT, *) RACCPC, MAXDTR, MAXDRT
            ACCPTC = 0
         ENDIF          
            
      ENDDO

      CALL ENRG (PEEND)

      WRITE (MYUNIT, *) PEEND, PE, RACCPC

      OPEN (UNIT = 31, FILE = 'finalpos.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 32, FILE = 'finalortn.dat', STATUS = 'UNKNOWN')
      
      DO I = 1, NPART
         WRITE (31, *) R(I,1), R(I,2), R(I,3)
         WRITE (32, *) Q(I,1), Q(I,2), Q(I,3), Q(I,4)
      ENDDO

      CALL VIEWCONFIG()

      CLOSE(UNIT = 31)
      CLOSE(UNIT = 32)

!      WRITE(MYUNIT, *) 'Number of successful quenches = ', QCOUNT
      CLOSE(MYUNIT)
 
      STOP  
      END PROGRAM SHIFTDP_MC_NVT

!     ============================================================================================== 

      SUBROUTINE MOVE(PE, ACCPTC)

      USE COMMONS

      IMPLICIT NONE
     

      INTEGER          :: J, J1, J2, INDXP, ACCPTC
      DOUBLE PRECISION :: RO(3), QO(4), QN(4), RSO(NRBSITE,3), RM(3,3), EO(NRBSITE,3)
      DOUBLE PRECISION :: PE, DELE, ENRGO, ENRGN 
      DOUBLE PRECISION :: RANF, DUMMY, BLTZMN, rand   
      LOGICAL          :: INSIDET, PERCT, REJECTT

      DO J = 1, NPART

         INDXP = INT(NPART*RANF(DUMMY)) + 1

         IF (INDXP < 1 .OR. INDXP > NPART) THEN
            WRITE (MYUNIT,'(A,I5)') 'Error in the index of the selected particle', INDXP
            STOP
         ENDIF
!     Calculate energy contribution of the selected particle with the current position and orientation

         CALL SNENRG(INDXP,ENRGO)

         RO(:) = R(INDXP,:)
         QO(:) = Q(INDXP,:)
         DO J1 = 1, NRBSITE
            J2 = NRBSITE*(INDXP-1) + J1
            RSO(J1,:) = RS(J2,:)
         ENDDO

!     Propose perturbation in the translational coordinates of the selected particle
         R(INDXP,:) = R(INDXP,:) + (2.D0*RANF(DUMMY) - 1.D0)*MAXDTR

!     Now propose perturbation in the rotational coordinates of the same particle

         CALL QTNSTEP(QO,MAXDRT,QN)

         CALL ROTMAT(QN,RM)

         DO J1 = 1, NRBSITE
            J2 = NRBSITE*(INDXP-1) + J1
            EO(J1,:) = E(J2,:)
            RS(J2,:) = R(INDXP,:) + MATMUL(RM,RBSITE(J1,:))
            E(J2,:) = MATMUL(RM,RBUV(J1,:))
         ENDDO
         Q(INDXP,:) = QN(:)

!     Calculate energy contribution of the selected particle with the proposed position and orientation

         CALL SNENRG(INDXP,ENRGN)

         DELE = ENRGN - ENRGO 

         REJECTT = .FALSE.

         IF (DELE > 0.D0) THEN
            BLTZMN = DEXP(- DELE / TMP)
            IF ((RANF(DUMMY) > BLTZMN)) THEN
               REJECTT = .TRUE.         ! REJECTION
            ELSE
!               CALL PERCOLATE(PERCT)
!               IF (PERCT) THEN
!                  ACCPTC = ACCPTC + 1
!               ELSE
!                  REJECTT = .TRUE.
!               ENDIF
               CALL CONTAINER(INSIDET)
               IF (INSIDET) THEN
                  ACCPTC = ACCPTC + 1
               ELSE
                  REJECTT = .TRUE.
               ENDIF
            ENDIF
         ELSE
!            CALL PERCOLATE(PERCT)
!            IF (PERCT) THEN
!               ACCPTC = ACCPTC + 1
!            ELSE
!               REJECTT = .TRUE.
!            ENDIF
             CALL CONTAINER(INSIDET)
             IF (INSIDET) THEN
                ACCPTC = ACCPTC + 1
             ELSE
                REJECTT = .TRUE.
             ENDIF
         ENDIF

         IF (REJECTT) THEN
            Q(INDXP,:) = QO(:)
            R(INDXP,:) = RO(:)
            DO J1 = 1, NRBSITE
               J2 = NRBSITE*(INDXP-1) + J1
               RS(J2,:) = RSO(J1,:)
               E(J2,:) = EO(J1,:)
            ENDDO
            DELE = 0.D0
         ENDIF

         PE = PE + DELE

      ENDDO

      END SUBROUTINE MOVE

!     ============================================================================================== 
           
      SUBROUTINE SNENRG(J2,PES)

      USE COMMONS

      IMPLICIT NONE
     
      INTEGER          :: J1, J2, J7, J8, I, J
      DOUBLE PRECISION :: RIJ(3), RIJSQ, ABSRIJ, R2, R4, EXPFCT, RSS(3), RSSSQ
      DOUBLE PRECISION :: PES, ALP, BET, GAM, VR, VA, VB, VG, DPFCT, NR(3), EI(3), EJ(3) 

      PES  = 0.D0

      DO J = 1, NRBSITE
         J8 = NRBSITE*(J2-1) + J
         EJ(:) = E(J8,:)
         DO J1 = 1, NPART 
            IF (J1 == J2) CYCLE
            RIJ(:) = R(J1,:) - R(J2,:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
!     Isotropic interaction between the spherical cores described by the Yukawa potential
            ABSRIJ = DSQRT(RIJSQ)
            R2     = 1.D0/RIJSQ
            EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
            PES = PES + EXPFCT/ABSRIJ
!     Dipolar contribution
            DO I = 1, NRBSITE
               J7    = NRBSITE*(J1-1) + I
               EI(:) = E(J7,:)
               RSS(:) = RS(J7,:) - RS(J8,:)
               R2     = DOT_PRODUCT(RSS(:),RSS(:))
               ABSRIJ = DSQRT(R2)
               NR(:)  = RSS(:)/ABSRIJ
               R2     = 1.D0/R2
               R4     = R2*R2
               ALP    = DOT_PRODUCT(NR(:),EI(:))
               BET    = DOT_PRODUCT(NR(:),EJ(:))
               GAM    = DOT_PRODUCT(EI(:),EJ(:))
               DPFCT  = 3.D0*DPMU*DPMU
               PES = PES + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ
            ENDDO
         ENDDO
      ENDDO
      IF (FIELDT) THEN
               PES = PES - DPMU*FIELD*E(J2,3)
      ENDIF
      END SUBROUTINE SNENRG

!     ============================================================================================== 

      SUBROUTINE ENRG(PE)

      USE COMMONS

      IMPLICIT NONE
     
      INTEGER          :: J1, J2, J3, J4, J5, J7, J8, I, J
      DOUBLE PRECISION :: RI(3), RIJ(3), QI(4), RIJSQ, ABSRIJ, R2, R4, EXPFCT, RSS(3), RSSSQ, RM(3,3)
      DOUBLE PRECISION :: PE, ALP, BET, GAM, VR, VA, VB, VG, DPFCT, NR(3), EI(3), EJ(3) 

      PE  = 0.D0

      DO J1 = 1, NPART
         RI(:) = R(J1,:)
         QI(:) = Q(J1,:)

         CALL ROTMAT(QI, RM)

         DO J2 = 1, NRBSITE
            J3 = NRBSITE*(J1-1) + J2
            RS(J3,:) = RI(:) + MATMUL(RM,RBSITE(J2,:))
            E(J3,:)  = MATMUL(RM(:,:),RBUV(J2,:))
         ENDDO
      ENDDO

      DO J1 = 1, NPART - 1
         J3 = 3*J1
         DO J2 = J1 + 1, NPART
            J4 = 3*J2
            RIJ(:) = R(J1,:) - R(J2,:)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))

!     Isotropic interaction between the spherical cores described by the Yukawa potential
            ABSRIJ = DSQRT(RIJSQ)
            R2     = 1.D0/RIJSQ
            EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
            PE     = PE + EXPFCT/ABSRIJ
!     Dipolar contribution
            DO I = 1, NRBSITE
               J7    = NRBSITE*(J1-1) + I
               EI(:) = E(J7,:)
               DO J = 1, NRBSITE
                  J8     = NRBSITE*(J2-1) + J
                  EJ(:)  = E(J8,:)
                  RSS(:) = RS(J7,:) - RS(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
                  R2     = 1.D0/R2
                  R4     = R2*R2
                  ALP    = DOT_PRODUCT(NR(:),EI(:))
                  BET    = DOT_PRODUCT(NR(:),EJ(:))
                  GAM    = DOT_PRODUCT(EI(:),EJ(:))
                  DPFCT  = 3.D0*DPMU*DPMU
                  PE     = PE + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF (FIELDT) THEN
         DO J1 = 1, NPART
            J3 = 3*J1
            DO I = 1, NRBSITE
               J7    = NRBSITE*(J1-1) + I
               PE = PE - DPMU*FIELD*E(J7,3)
            ENDDO
         ENDDO
      ENDIF


!      stop

      END SUBROUTINE ENRG

!     ==============================================================================================

      SUBROUTINE ROTMAT(Q,RM)

!     provides the rotation matrix RM corresponding to quaternion Q
!     right-handed rotation in the right-handed coordinate system

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: Q(4)
      DOUBLE PRECISION, INTENT(OUT) :: RM(3,3)

      RM(1,1) = Q(1)**2 + Q(2)**2 - Q(3)**2 - Q(4)**2
      RM(1,2) = 2.D0*(Q(2)*Q(3) - Q(1)*Q(4))
      RM(1,3) = 2.D0*(Q(2)*Q(4) + Q(1)*Q(3))
      RM(2,1) = 2.D0*(Q(2)*Q(3) + Q(1)*Q(4))
      RM(2,2) = Q(1)**2 + Q(3)**2 - Q(2)**2 - Q(4)**2
      RM(2,3) = 2.D0*(Q(3)*Q(4) - Q(1)*Q(2))
      RM(3,1) = 2.D0*(Q(2)*Q(4) - Q(1)*Q(3))
      RM(3,2) = 2.D0*(Q(3)*Q(4) + Q(1)*Q(2))
      RM(3,3) = Q(1)**2 + Q(4)**2 - Q(2)**2 - Q(3)**2

      END SUBROUTINE

!     ============================================================================================== 

      SUBROUTINE QTNSTEP(Q1,MAXDRT,Q2)
!     proposes a rotational step in terms of a unit quaternion such that a small step size
!     corresponds to a small perturbation of the current orientation and large step size produces
!     a completely random orientation uncorrelated with the exisiting orientation
!     See Trosset and Scheraga, J. Comput. Chem. 20, 412-427 (1999)

      IMPLICIT NONE
      DOUBLE PRECISION :: Q1(4), Q2(4), MAXDRT, DUMMY, GAUSS, FCTR, XISQ, XI1, XI2, XI3, XI4 

      XISQ = 1.0
!     iterative loop
      DO WHILE (XISQ >= 1.D0)
         XI2  = GAUSS(DUMMY)*MAXDRT
         XI3  = GAUSS(DUMMY)*MAXDRT
         XI4  = GAUSS(DUMMY)*MAXDRT
         XISQ = XI2*XI2 + XI3*XI3 + XI4*XI4
      ENDDO
      XI1 = SQRT(1.D0-XISQ)
      Q2(:) = (/XI1, XI2, XI3, XI4/)

!     Now transform the current quaternion Q1 corresponding to rotation via quaternion Q2 and return the
!     transformed quaternion as Q2

      XI1  = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
      XI2  = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
      XI3  = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
      XI4  = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

      Q2(:) = (/XI1, XI2, XI3, XI4/)

      END SUBROUTINE QTNSTEP

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

      SUBROUTINE EQCON()

      USE COMMONS

      IMPLICIT NONE

      INTEGER :: I

      OPEN (UNIT=13, FILE='initialpos.dat', STATUS='UNKNOWN')
      OPEN (UNIT=14, FILE='initialortn.dat', STATUS='UNKNOWN')

      DO I = 1, NPART  
         READ(13, *) R(I,1), R(I,2), R(I,3)
         READ(14, *) Q(I,1), Q(I,2), Q(I,3), Q(I,4)
      ENDDO

      CLOSE (13)
      CLOSE (14)

      END SUBROUTINE EQCON

!     ==============================================================================================

      DOUBLE PRECISION FUNCTION RANF(DUMMY)
!     Draws a uniform random variate between 0 and 1 

      INTEGER, PARAMETER :: L = 1029, C = 221591, M = 1048576
      INTEGER ::          SEED
      DOUBLE PRECISION :: DUMMY
      SAVE             SEED
      DATA             SEED / 0 /

      SEED = MOD(SEED * L + C, M)
      RANF = DFLOAT(SEED) / DFLOAT(M)

      END FUNCTION RANF

!     ==============================================================================================

      DOUBLE PRECISION FUNCTION GAUSS(DUMMY)
!     Gaussian distribution with zero mean and unit variance
      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: SUMRND, RV, RVSQ
      DOUBLE PRECISION :: RANF, DUMMY
      DOUBLE PRECISION, PARAMETER :: A1 = 3.949846138D0, A3 = 0.252408784D0
      DOUBLE PRECISION, PARAMETER :: A5 = 0.076542912D0, A7 = 0.008355968D0, A9 = 0.029899776D0

      SUMRND = 0.D0

      DO J1 = 1, 12
         SUMRND = SUMRND + RANF(DUMMY)
      ENDDO

      RV    = (SUMRND - 6.D0)/4.D0
      RVSQ  = RV*RV
      GAUSS = ((((A9*RVSQ + A7)*RVSQ + A5)*RVSQ + A3)*RVSQ + A1)*RV

      END FUNCTION GAUSS

!     ==============================================================================================

      SUBROUTINE DEFSHIFTDP()

      USE COMMONS, ONLY: RBSITE, RBUV, RSHIFT

      IMPLICIT NONE

      RBSITE(1,:) = (/0.D0, 0.D0, RSHIFT/)
      RBUV(1,:)   = (/ 0.D0, 0.D0, 1.D0/)

      END SUBROUTINE DEFSHIFTDP

!     ==============================================================================================

      SUBROUTINE VIEWCONFIG()

      USE COMMONS

      IMPLICIT NONE

      INTEGER :: J1, J2
      DOUBLE PRECISION :: QI(4), RM(3,3), SR(NRBSITE,3), UV(3)

      OPEN (UNIT = 31, FILE = 'finalconfig.xyz', STATUS = 'UNKNOWN')

      WRITE(31,'(I6)') NPART*(NRBSITE+1)
      WRITE(31, *)
      DO J1 = 1, NPART
         QI(:) = Q(J1,:)
         CALL ROTMAT(QI, RM)
         WRITE(31,'(A5,1X,3F20.10)') 'O ', R(J1,1), R(J1,2), R(J1,3)
         DO J2 = 1, NRBSITE
            SR(J2,:) = R(J1,:) + MATMUL(RM,RBSITE(J2,:))
            UV(:)     = MATMUL(RM(:,:),RBUV(J2,:))
            WRITE(31,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'C', SR(J2,1), SR(J2,2), SR(J2,3), &
                'atom_vector', UV(1), UV(2), UV(3)
         ENDDO
      ENDDO

      CLOSE (UNIT = 31)

      END SUBROUTINE VIEWCONFIG

!     ==============================================================================================

      SUBROUTINE QUENCH(PE)

!
      USE COMMONS, ONLY: NPART, MYUNIT, QCOUNT, R, Q

      IMPLICIT NONE

      INTEGER          :: J1, J2, NOPT, ITDONE 
      DOUBLE PRECISION :: XCOORDS(6*NPART), ENERGY, PE, QJ(4), P(3)
      LOGICAL          :: IFLAG

      NOPT = 6*NPART

      DO J1 = 1, NPART

         QJ(:) = Q(J1,:)

         CALL RQTOAA(QJ,P)         

         J2 = 3*J1
         XCOORDS(J2-2:J2) = R(J1,:)
         XCOORDS(3*NPART+J2-2:3*NPART+J2) = P(:)

      ENDDO

      CALL MYLBFGS(NOPT,XCOORDS,.FALSE.,IFLAG,ENERGY,ITDONE,PE)

      IF (IFLAG) THEN

         OPEN (UNIT=27, FILE='qconfig.dat',  STATUS='UNKNOWN', ACCESS='APPEND')
         OPEN (UNIT=28, FILE='qenergy.dat',  STATUS='UNKNOWN', ACCESS='APPEND')

         QCOUNT = QCOUNT + 1

         WRITE(28, *) QCOUNT, ENERGY

         DO J1 = 1, 2*NPART
            J2 = 3*J1
            WRITE (27, *) XCOORDS(J2-2), XCOORDS(J2-1), XCOORDS(J2)
         ENDDO

         CLOSE (UNIT=27, STATUS='KEEP')
         CLOSE (UNIT=28, STATUS='KEEP')


      ELSE
         WRITE(MYUNIT, *) ' WARNING: quench is not converged'
      ENDIF

      END SUBROUTINE QUENCH

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

      SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)

      IMPLICIT NONE

      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, I3(3,3), E(3,3), ESQ(3,3)
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      LOGICAL          :: GTEST

!     P(3)     : rotation vector 
!     RM(3,3)  : rotation matrix
!     DRMk(3,3): derivative of the rotation matrix with respect to the kth component of the rotation vector
!     GTEST    : true if derivatives are to be found  
!     PN(3)    : the unit vector parallel to P
!     THETA    : the modulus of the rotation vector P, equivalent to the angle of rotation
!     THETA2   : THETA squared
!     THETA3   : THETA**(-3)
!     CT       : cos(THETA)
!     ST       : sin(THETA)
!     I3(3,3)  : 3x3 identity matrix
!     E(3,3)   : the skew-symmetric matrix obtained from a unit vector parallel to P (equation (2) in the paper)
!     ESQ(3,3) : the square of E
!     DEk      : derivate of E with respect to the kth component of P

!     Set the values of the idenity matrix I3
      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

!     Calculate the value of THETA2 as the square modulus of P
      THETA2  = DOT_PRODUCT(P,P)

      IF (THETA2 < 1.0D-12) THEN
!     Execute if the angle of rotation is zero
!     In this case the rotation matrix is the identity matrix
         RM(:,:) = I3(:,:)

!     First order corrections to rotation matrix
         RM(1,2) =-P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) =-P(2)
         RM(2,3) =-P(1)
         RM(3,2) = P(1)

!     If derivatives do not need to found, we're finished
         IF (.NOT. GTEST) RETURN

!     This is the special case described in the paper, where DRMk is equal to E(k), which is the skew-symmetric matrix 
!     obtained from P with P(k) equal to 1 and other components equal to zero
!         PN        = (/1.D0, 0.D0, 0.D0/)
!         E(:,:)    = 0.D0
!         E(2,3)    = -PN(1)
!         E(3,2)    = -E(2,3)
!         DRM1(:,:) = E(:,:)

!         PN        = (/0.D0, 1.D0, 0.D0/)
!         E(:,:)    = 0.D0
!         E(1,3)    =  PN(2)
!         E(3,1)    = -E(1,3)
!         DRM2(:,:) = E(:,:)

!         PN        = (/0.D0, 0.D0, 1.D0/)
!         E(:,:)    = 0.D0
!         E(1,2)    = -PN(3)
!         E(2,1)    = -E(1,2)
!         DRM3(:,:) = E(:,:)

!     Now up to the linear order in theta
         E(:,:)    = 0.D0
         E(1,1)    = 0.0D0
         E(1,2)    = P(2)
         E(1,3)    = P(3)
         E(2,1)    = P(2)
         E(2,2)    =-2.0D0*P(1)
         E(2,3)    =-2.0D0
         E(3,1)    = P(3)
         E(3,2)    = 2.0D0
         E(3,3)    =-2.0D0*P(1)
         DRM1(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    =-2.0D0*P(2)
         E(1,2)    = P(1)
         E(1,3)    = 2.0D0
         E(2,1)    = P(1)
         E(2,2)    = 0.0D0
         E(2,3)    = P(3)
         E(3,1)    =-2.0D0
         E(3,2)    = P(3)
         E(3,3)    =-2.0D0*P(2)
         DRM2(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    =-2.0D0*P(3)
         E(1,2)    =-2.0D0
         E(1,3)    = P(1)
         E(2,1)    = 2.0D0
         E(2,2)    =-2.0D0*P(3)
         E(2,3)    = P(2)
         E(3,1)    = P(1)
         E(3,2)    = P(2)
         E(3,3)    = 0.0D0
         DRM3(:,:) = 0.5D0*E(:,:)

      ELSE
!     Execute for the general case, where THETA dos not equal zero
!     Find values of THETA, CT, ST and THETA3
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)

!     Set THETA to 1/THETA purely for convenience
         THETA   = 1.D0/THETA

!     Normalise P and construct the skew-symmetric matrix E
!     ESQ is calculated as the square of E
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         ESQ     = MATMUL(E,E)

!     RM is calculated from Rodrigues' rotation formula (equation (1) in the paper)
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

!     If derivatives do not need to found, we are finished
         IF (.NOT. GTEST) RETURN

!     Set up DEk using the form given in equation (4) in the paper
         DE1(:,:) = 0.D0
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
         DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

!     Use equation (3) in the paper to find DRMk
         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) + CT*PN(1)*E(:,:) + ST*DE1(:,:)
         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) + CT*PN(2)*E(:,:) + ST*DE2(:,:)
         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) + CT*PN(3)*E(:,:) + ST*DE3(:,:)

      ENDIF
      END SUBROUTINE RMDRVT

!     ==============================================================================================

      SUBROUTINE MYLBFGS(NOPT,XCOORDS,DIAGCO,IFLAG,ENERGY,ITDONE,PE)

!     NOPT    An INTEGER variable corresponding to the number of variables to be optimised
! 
!     M       An INTEGER variable set to the number of corrections used in the BFGS update. 
!             Values of M less than 3 are not recommended; large values of M will result in 
!             excessive computing time. 3<= M <=7 is recommended. Restriction: M>0.
! 
!     XCOORDS A DOUBLE PRECISION array of length NOPT holding the coordinates
!             On exit with IFLAG=TRUE, it contains the values of the variables at the stationary point.
!
!     DIAGCO  A LOGICAL variable that must be set to .TRUE. if the user  wishes to provide the 
!             diagonal matrix Hk0 at each iteration. Otherwise it should be set to .FALSE., in which
!             case  LBFGS will use a default value described below. If DIAGCO is set to .TRUE. the 
!             routine will return at each iteration of the algorithm with IFLAG=2, and the diagonal
!             matrix Hk0  must be provided in the array DIAG.
!
!     IFLAG  A LOGICAL variable set to .TRUE. if optimisation is converged 
! 
!     GMAX    A positive  DOUBLE PRECISION variable that must be set and determines the accuracy of 
!             convergence.
!
!     ITMAX   An INTEGER variable defining the maximum number of LBFGS iterations allowed
!
!     ENERGY  A DOUBLE PRECISION variable containing the value of the function at the point XCOORDS.
!
!     ITDONE  An INTEGER variable counting the number LBFGS iterations performed
! 
!     GRAD    A DOUBLE PRECISION array of length NOPT containing the components of the gradient GRAD at
!             the point XCOORDS.
!
!     DIAG    A DOUBLE PRECISION array of length NOPT. If DIAGCO=.TRUE., then on initial entry or on 
!             re-entry with IFLAG=2, DIAG it must be set by the user to contain the values of the 
!             diagonal matrix Hk0.  Restriction: all elements of DIAG must be positive.

      USE COMMONS, ONLY: MYUNIT, NPART, R, Q

      IMPLICIT NONE

      INTEGER          :: NOPT, M, ITDONE, ITR, NDECREASE, NFAIL, I, J, ITMAX
      INTEGER          :: BOUND, CP, ISPT, IYCN, IYPT, INMC, ISCN, NPT, POINT
      DOUBLE PRECISION :: XCOORDS(NOPT), GRAD(NOPT), ENERGY, DIAG(NOPT), GNORM, PE 
      DOUBLE PRECISION :: DUMMY, SQ, YR, YS, YY, BETA
      DOUBLE PRECISION :: STP, SLENGTH, DOT1, DOT2, OVERLAP, XSAVE(NOPT), GNEW(NOPT), ENEW
      DOUBLE PRECISION :: DGUESS, GMAX, MAXBFGS, MAXEFALL, MAXERISE, RMS
      DOUBLE PRECISION, ALLOCATABLE :: W(:)
      LOGICAL          :: DIAGCO, IFLAG, RESET, DEBUG

!     Assign parameter values

      DGUESS    = 0.1D0
      GMAX      = 1.D-06
      MAXBFGS   = 0.1D0
      MAXEFALL  =-HUGE(1.0D0)
      MAXERISE  = 1.d-10
      ITMAX     = 100000
      M         = 4
      DEBUG     = .FALSE.

      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(NOPT*(2*M+1)+2*M))       

!     Initialise
      ITDONE = 0
      NFAIL  = 0
      ITR    = 0

      CALL POTENTIAL(XCOORDS,GRAD,RMS,ENERGY,.TRUE.)

      IF (ABS(PE-ENERGY) > 1.D-04) THEN
         WRITE(MYUNIT,*) ' ERROR in coordinate transformation for quenching'
         WRITE(MYUNIT,*) ' Original energy = ', PE
         WRITE(MYUNIT,*) ' Transformed energy = ', ENERGY
         WRITE(MYUNIT,*) ' Original coordinates'
         DO I = 1, NPART
            WRITE (MYUNIT, *) R(I,1), R(I,2), R(I,3)
            WRITE (MYUNIT, *) Q(I,1), Q(I,2), Q(I,3), Q(I,4)
         ENDDO
         WRITE(MYUNIT,*) ' Transformed coordinates'
         DO I = 1, NPART
            J = 3*I
            WRITE (MYUNIT, *) XCOORDS(J-2), XCOORDS(J-1), XCOORDS(J)
         ENDDO
         STOP
      ENDIF

      IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force = ', ENERGY, RMS,' after ', ITDONE, ' LBFGS steps'

!     Check for convergence and hence termination

10    IFLAG = .FALSE.

      IF (RMS <= GMAX) THEN
         IFLAG = .TRUE.
         IF (IFLAG) THEN
            WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force = ', ENERGY, RMS,      &
     &                                                          ' after ', ITDONE, ' LBFGS steps'
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE == ITMAX) THEN
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ', DIAG(1)
         RETURN
      ENDIF

      IF (ITR == 0) THEN
         IF (NOPT <= 0 .OR. M <= 0) THEN
            WRITE(MYUNIT,240)
240         FORMAT(' IMPROPER INPUT PARAMETERS (NOPT OR M IS NOT POSITIVE)')
            STOP
         ENDIF
         IFLAG = .FALSE.
         IF (DIAGCO) THEN
            WRITE(MYUNIT,'(A)') 'using estimate of the inverse diagonal elements'
            DO I = 1, NOPT
               IF (DIAG(I) <= 0.0D0) THEN
                  WRITE(MYUNIT,235) I
235               FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
!     Initial guess for diagonal elements of the inverse Hessian, used whenever the LBFGS optimiser 
!     is reset. The default is DGUESS = 0.1.
            DO I = 1, NOPT
               DIAG(I) = DGUESS
            ENDDO
         ENDIF

!     The work vector W is divided as follows:
!     ---------------------------------------
!     The first NOPT elements are used to store the gradient information.
!     Elements (NOPT+1),...,(NOPT+M) store the scalars RHO.
!     Elements (NOPT+M+1),...,(NOPT+2M) store the numbers ALPHA used in the formula computing H*G.
!     Elements (NOPT+2M+1),...,(NOPT+2M+NOPT*M) store the last M search steps: s_{k} = x_{k+1} - x_{k}. 
!     Elements (NOPT+2M+NOPT*M+1),...,(NOPT+2M+2NOPT*M) store the last M gradient 
!     differences: y_{k} = g_{k+1) - g_{k}.
!
!     The search steps and gradient differences are stored in a circular order controlled by the 
!     parameter POINT.

         POINT  = 0    
         ISPT   = NOPT + 2*M            ! index for storage of search steps 
         IYPT   = ISPT + NOPT*M         ! index for storage of gradient differences

!     Initial search direction: d_{0} = -H_{0}g_{0} 
 
         DO I = 1, NOPT

            DUMMY     =-GRAD(I)*DIAG(I)
            W(ISPT+I) = DUMMY               
            W(I)      = DUMMY
         
         ENDDO

         GNORM = DSQRT(DOT_PRODUCT(GRAD(1:NOPT),GRAD(1:NOPT)))

!     First guess for the step length

         STP = MIN(1.0D0/GNORM,GNORM)

      ELSE

!     Compute -H*G using the formula given by J. Nocedal, 
!     Mathematics of Computation, Vol.35, No.151, pp. 773-782 (1980).
!     ------------------------------------------------------------------

         BOUND = ITR
         IF (ITR > M) BOUND = M
         
         YS = DOT_PRODUCT(W(IYPT+NPT+1:IYPT+NPT+NOPT),W(ISPT+NPT+1:ISPT+NPT+NOPT))
        
         IF (.NOT. DIAGCO) THEN
            YY= DOT_PRODUCT(W(IYPT+NPT+1:IYPT+NPT+NOPT),W(IYPT+NPT+1:IYPT+NPT+NOPT))
            IF (YY == 0.0D0) THEN
               WRITE(MYUNIT,'(A)') 'WARNING, resetting YY to 1 in mylbfgs'
               YY = 1.0D0
            ENDIF
            IF (YS == 0.0D0) THEN
               WRITE(MYUNIT,'(A)') 'WARNING, resetting YS to 1 in mylbfgs'
               YS = 1.0D0
            ENDIF
            DO I = 1, NOPT
               DIAG(I) = YS/YY
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)') 'using estimate of the inverse diagonal elements'
            DO I = 1, NOPT
               IF (DIAG(I) <= 0.0D0) THEN
                  WRITE(MYUNIT,235) I
                  STOP
               ENDIF
            ENDDO
         ENDIF

         CP = POINT
         IF (POINT == 0) CP = M         !???
         W(NOPT+CP) = 1.0D0/YS          !???
         
         W(1:NOPT)= -GRAD(1:NOPT)
         
         CP = POINT

         DO I = 1, BOUND

            CP      = CP - 1
            IF (CP == -1) CP = M - 1
            SQ        = DOT_PRODUCT(W(ISPT+CP*NOPT+1:ISPT+(CP+1)*NOPT),W(1:NOPT)) ! s*q
            INMC      = NOPT+M+CP+1
            IYCN      = IYPT+CP*NOPT
            W(INMC)   = W(NOPT+CP+1) * SQ ! alpha = rho*s*q
            W(1:NOPT) = W(1:NOPT) - W(INMC)*W(IYCN+1:IYCN+NOPT) ! q = q - alpha * y

         ENDDO

         DO I = 1, NOPT
            W(I) = DIAG(I)*W(I) ! r_{0}
         ENDDO

         DO I = 1,BOUND

            YR        = DOT_PRODUCT(W(IYPT+CP*NOPT+1:IYPT+(CP+1)*NOPT),W(1:NOPT)) ! y*r
            BETA      = W(NOPT+CP+1) * YR ! beta = rho*y*r 
            INMC      = NOPT+M+CP+1
            BETA      = W(INMC) - BETA ! alpha - beta
            ISCN      = ISPT+CP*NOPT
            W(1:NOPT) = W(1:NOPT) +  BETA*W(ISCN+1:ISCN+NOPT) 
            CP     = CP + 1
            IF (CP == M) CP = 0

         ENDDO

         DO I = 1, NOPT
            W(ISPT+POINT*NOPT+I)= W(I)
         ENDDO

         STP = 1.0D0
         
      ENDIF

!     Propose Steplength 

      DOT1 = DSQRT(DOT_PRODUCT(GRAD(1:NOPT),GRAD(1:NOPT)))
      DOT2 = DSQRT(DOT_PRODUCT(W(1:NOPT),W(1:NOPT)))
      OVERLAP = 0.0D0
      IF (DOT1*DOT2 /= 0.0D0) THEN
         OVERLAP = DOT_PRODUCT(GRAD(1:NOPT),W(1:NOPT))/(DOT1*DOT2)
      ENDIF
      IF (OVERLAP > 0.0D0) THEN
         IF (DEBUG) WRITE(MYUNIT,'(A)') 'Search direction has positive projection onto gradient - reversing step'
         DO I = 1, NOPT
            W(ISPT+POINT*NOPT+I)=-W(I)  !!! DJW, reverses step
         ENDDO
      ENDIF

      W(1:NOPT) = GRAD(1:NOPT)

      SLENGTH = 0.0D0
      DO I = 1, NOPT
         SLENGTH = SLENGTH + W(ISPT+POINT*NOPT+I)**2
      ENDDO
      SLENGTH = SQRT(SLENGTH)
      IF (STP*SLENGTH > MAXBFGS) STP = MAXBFGS/SLENGTH

!     We now have the proposed step.

!     Save XCOORDS here so that we can undo the step reliably.

      XSAVE(1:NOPT) = XCOORDS(1:NOPT)
   
!     Take step

!      XCOORDS(1:NOPT) = XCOORDS(1:NOPT) + STP*W(ISPT+POINT*NOPT+1:ISPT+(POINT+1)*NOPT) ! x_{k+1) = x_{k} + step_{k}*d_{k}

      DO I = 1, NOPT
         XCOORDS(I) = XCOORDS(I)+STP*W(ISPT+POINT*NOPT+I)
      ENDDO

      NDECREASE = 0

20    CONTINUE

      CALL POTENTIAL(XCOORDS,GNEW,RMS,ENEW,.TRUE.)      ! g_{k+1}

!     IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force = ', ENEW, RMS,' after ', ITDONE, ' LBFGS steps'

!     Check whether the proposed step is accepted or not - alternative to line search

      IF (((ENEW - ENERGY) <= MAXERISE) .AND. ((ENEW-ENERGY) > MAXEFALL)) THEN

         ITR = ITR + 1
         ITDONE = ITDONE + 1
         ENERGY = ENEW
         GRAD(1:NOPT) = GNEW(1:NOPT)
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A,F13.10)') ' Energy and RMS force=',ENERGY,RMS,      &
     &' after ', ITDONE,' LBFGS steps, step:', STP*SLENGTH

      ELSEIF ((ENEW - ENERGY) <= MAXEFALL) THEN

!     Energy decreased too much - try again with a smaller step size

!     Resetting to XSAVE and adding half the step should be the same as subtracting half the step. 

         XCOORDS(1:NOPT) = XSAVE(1:NOPT)
         STP = 0.5D0*STP
         DO I = 1, NOPT
            XCOORDS(I) = XCOORDS(I) + STP*W(ISPT+POINT*NOPT+I)
         ENDDO  
         NDECREASE = NDECREASE + 1
         IF (DEBUG) WRITE(MYUNIT,'(A,F19.10,A,F16.10,A,F15.8)')                 & 
     &   ' energy decreased too much from ', ENERGY,' to ', ENEW,' decreasing step to ', STP*SLENGTH
         GOTO 20

         IF (NDECREASE > 5) THEN
            NFAIL = NFAIL + 1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL = ', NFAIL

!     Resetting to XSAVE should be the same as subtracting the step. 

            XCOORDS(1:NOPT) = XSAVE(1:NOPT)
            GRAD(1:NOPT)    = GNEW(1:NOPT)        ! GRAD contains the gradient at the lowest energy point
            ITR = 0                               ! try resetting
            IF (NFAIL > 20) THEN
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               RETURN
            ENDIF
            GOTO 30
         ENDIF

      ELSE

!     Energy increased - try again with a smaller step size

         XCOORDS(1:NOPT)=XSAVE(1:NOPT)
         STP = 0.1D0*STP
         DO I = 1, NOPT
            XCOORDS(I) = XCOORDS(I) + STP*W(ISPT+POINT*NOPT+I)
         ENDDO
         NDECREASE = NDECREASE + 1
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,A,F20.10,A,F20.10)')                &
     &              ' energy increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         GOTO 20

         IF (NDECREASE > 5) THEN 
            NFAIL = NFAIL + 1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL = ', NFAIL
            XCOORDS(1:NOPT) = XSAVE(1:NOPT)
            ITR = 0                                 !  try resetting
            IF (NFAIL > 20) THEN
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               RETURN
            ENDIF
            GOTO 30
         ENDIF

      ENDIF

30    NPT = POINT*NOPT

      DO I = 1, NOPT

         W(ISPT+NPT+I) = STP*W(ISPT+NPT+I)    ! s_{k} = x_{k+1} - x_{k} = step_{k}*d_{k}; save the step taken
         W(IYPT+NPT+I) = GRAD(I) - W(I)       ! y_{k} = g_{k+1} - g_{k}; save gradient difference - W(1:NOPT) contains the old gradient

      ENDDO

      POINT = POINT + 1
      IF (POINT == M) POINT = 0

      GOTO 10

      END SUBROUTINE MYLBFGS


!     ==============================================================================================

      SUBROUTINE POTENTIAL(X, G, RMS, ENERGY, GTEST)

      USE COMMONS, ONLY: NPART, NRBSITE, RBSITE, RBUV, DPMU, YKAPPA, FIELD, FIELDT

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET 
      DOUBLE PRECISION :: X(6*NPART), G(6*NPART)
      DOUBLE PRECISION :: ENERGY, RMS, RLJN, R2LJN, R2, R4, ABSRIJ, RIJSQ, DVDR, EXPFCT, FCTR
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), RSS(3), NR(3), P(3), EI(3), EJ(3)
      DOUBLE PRECISION :: R(NPART*NRBSITE,3), E(NPART*NRBSITE,3)
      DOUBLE PRECISION :: DR1(NPART*NRBSITE,3), DR2(NPART*NRBSITE,3), DR3(NPART*NRBSITE,3) 
      DOUBLE PRECISION :: DE1(NPART*NRBSITE,3), DE2(NPART*NRBSITE,3), DE3(NPART*NRBSITE,3)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: ALP, BET, GAM, VR, VA, VB, VG, FIJN, FIJEI, FIJEJ, FIJ(3), DPFCT
      DOUBLE PRECISION :: DOTI1, DOTI2, DOTI3, DOTJ1, DOTJ2, DOTJ3
      DOUBLE PRECISION :: DADPI1, DADPI2, DADPI3, DADPJ1, DADPJ2, DADPJ3
      DOUBLE PRECISION :: DBDPI1, DBDPI2, DBDPI3, DBDPJ1, DBDPJ2, DBDPJ3
      DOUBLE PRECISION :: DGDPI1, DGDPI2, DGDPI3, DGDPJ1, DGDPJ2, DGDPJ3
      LOGICAL          :: GTEST

      ENERGY = 0.D0
 
      IF (GTEST) G(:) = 0.D0

      OFFSET = 3*NPART
  
      DO J1 = 1, NPART
         J3 = 3*J1
         J5 = OFFSET + J3
         RI(:) = X(J3-2:J3)
         P(:)  = X(J5-2:J5)
         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITE
            J4        = NRBSITE*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI(:,:),RBSITE(J2,:))
            E(J4,:)   = MATMUL(RMI(:,:),RBUV(J2,:))
            IF (GTEST) THEN
               DR1(J4,:) = MATMUL(DRMI1(:,:),RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2(:,:),RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3(:,:),RBSITE(J2,:))

               DE1(J4,:) = MATMUL(DRMI1(:,:),RBUV(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBUV(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBUV(J2,:))
            ENDIF
         ENDDO
      ENDDO

      DO J1 = 1, NPART  
         J3 = 3*J1
         J5 = OFFSET + J3
         DO J2 = J1 + 1, NPART
            J4 = 3*J2
            J6 = OFFSET + J4
            RIJ(:) = X(J3-2:J3) - X(J4-2:J4)
            RIJSQ  = DOT_PRODUCT(RIJ(:),RIJ(:))
            ABSRIJ = DSQRT(RIJSQ)
            R2     = 1.D0/RIJSQ
            EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
            ENERGY = ENERGY + EXPFCT/ABSRIJ
            DVDR   =-EXPFCT*R2*(YKAPPA + 1.D0/ABSRIJ)
            IF (GTEST) THEN
               G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJ(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJ(:)
            ENDIF
!     Dipolar contribution
            DO I = 1, NRBSITE
               J7    = NRBSITE*(J1-1) + I
               EI(:) = E(J7,:)
               DO J = 1, NRBSITE
                  J8     = NRBSITE*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  R2     = DOT_PRODUCT(RSS(:),RSS(:))
                  ABSRIJ = DSQRT(R2)
                  NR(:)  = RSS(:)/ABSRIJ
                  R2     = 1.D0/R2
                  R4     = R2*R2
                  EJ(:)  = E(J8,:)
                  ALP    = DOT_PRODUCT(NR(:),EI(:))
                  BET    = DOT_PRODUCT(NR(:),EJ(:))
                  GAM    = DOT_PRODUCT(EI(:),EJ(:))
                  DPFCT  = 3.D0*DPMU*DPMU
                  ENERGY = ENERGY + DPFCT*R2*(GAM/3.D0 - ALP*BET)/ABSRIJ

                  IF (GTEST) THEN
                     VR     = -DPFCT*R4*(GAM - 3.D0*ALP*BET)
                     VA     = -DPFCT*BET*R2/ABSRIJ
                     VB     = -DPFCT*ALP*R2/ABSRIJ
                     VG     =  DPFCT*R2/(3.D0*ABSRIJ)

                     FIJN   = VR - (VA*ALP+VB*BET)/ABSRIJ
                     FIJEI  = VA/ABSRIJ
                     FIJEJ  = VB/ABSRIJ
                     FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

                     G(J3-2:J3) = G(J3-2:J3) + FIJ(:)
                     G(J4-2:J4) = G(J4-2:J4) - FIJ(:)

                     DOTI1 = DOT_PRODUCT(RSS,DR1(J7,:))
                     DOTI2 = DOT_PRODUCT(RSS,DR2(J7,:))
                     DOTI3 = DOT_PRODUCT(RSS,DR3(J7,:))

                     DOTJ1 =-DOT_PRODUCT(RSS,DR1(J8,:))
                     DOTJ2 =-DOT_PRODUCT(RSS,DR2(J8,:))
                     DOTJ3 =-DOT_PRODUCT(RSS,DR3(J8,:))

                     DADPI1 = DOT_PRODUCT(DR1(J7,:),EI(:))/ABSRIJ - ALP*R2*DOTI1 + DOT_PRODUCT(NR(:),DE1(J7,:))
                     DADPI2 = DOT_PRODUCT(DR2(J7,:),EI(:))/ABSRIJ - ALP*R2*DOTI2 + DOT_PRODUCT(NR(:),DE2(J7,:))
                     DADPI3 = DOT_PRODUCT(DR3(J7,:),EI(:))/ABSRIJ - ALP*R2*DOTI3 + DOT_PRODUCT(NR(:),DE3(J7,:))
                                   
                     DADPJ1 =-DOT_PRODUCT(DR1(J8,:),EI(:))/ABSRIJ - ALP*R2*DOTJ1
                     DADPJ2 =-DOT_PRODUCT(DR2(J8,:),EI(:))/ABSRIJ - ALP*R2*DOTJ2
                     DADPJ3 =-DOT_PRODUCT(DR3(J8,:),EI(:))/ABSRIJ - ALP*R2*DOTJ3

                     DBDPI1 = DOT_PRODUCT(DR1(J7,:),EJ(:))/ABSRIJ - BET*R2*DOTI1
                     DBDPI2 = DOT_PRODUCT(DR2(J7,:),EJ(:))/ABSRIJ - BET*R2*DOTI2
                     DBDPI3 = DOT_PRODUCT(DR3(J7,:),EJ(:))/ABSRIJ - BET*R2*DOTI3

                     DBDPJ1 =-DOT_PRODUCT(DR1(J8,:),EJ(:))/ABSRIJ - BET*R2*DOTJ1 + DOT_PRODUCT(NR(:),DE1(J8,:))
                     DBDPJ2 =-DOT_PRODUCT(DR2(J8,:),EJ(:))/ABSRIJ - BET*R2*DOTJ2 + DOT_PRODUCT(NR(:),DE2(J8,:))
                     DBDPJ3 =-DOT_PRODUCT(DR3(J8,:),EJ(:))/ABSRIJ - BET*R2*DOTJ3 + DOT_PRODUCT(NR(:),DE3(J8,:))
                                   
                     DGDPI1 = DOT_PRODUCT(DE1(J7,:),EJ(:))
                     DGDPI2 = DOT_PRODUCT(DE2(J7,:),EJ(:))
                     DGDPI3 = DOT_PRODUCT(DE3(J7,:),EJ(:))

                     DGDPJ1 = DOT_PRODUCT(EI(:),DE1(J8,:))
                     DGDPJ2 = DOT_PRODUCT(EI(:),DE2(J8,:))
                     DGDPJ3 = DOT_PRODUCT(EI(:),DE3(J8,:))

                     G(J5-2) = G(J5-2) + VR*DOTI1/ABSRIJ + VA*DADPI1 + VB*DBDPI1 + VG*DGDPI1
                     G(J5-1) = G(J5-1) + VR*DOTI2/ABSRIJ + VA*DADPI2 + VB*DBDPI2 + VG*DGDPI2
                     G(J5)   = G(J5)   + VR*DOTI3/ABSRIJ + VA*DADPI3 + VB*DBDPI3 + VG*DGDPI3

                     G(J6-2) = G(J6-2) + VR*DOTJ1/ABSRIJ + VA*DADPJ1 + VB*DBDPJ1 + VG*DGDPJ1
                     G(J6-1) = G(J6-1) + VR*DOTJ2/ABSRIJ + VA*DADPJ2 + VB*DBDPJ2 + VG*DGDPJ2
                     G(J6)   = G(J6)   + VR*DOTJ3/ABSRIJ + VA*DADPJ3 + VB*DBDPJ3 + VG*DGDPJ3
                  ENDIF 
               ENDDO
            ENDDO
 
            IF (FIELDT) THEN
               ENERGY = ENERGY - DPMU*FIELD*EI(3)
               IF (GTEST) THEN
                  G(J5-2) = G(J5-2) - DPMU*FIELD*DE1(J7,3)
                  G(J5-1) = G(J5-1) - DPMU*FIELD*DE2(J7,3)
                  G(J5)   = G(J5)   - DPMU*FIELD*DE3(J7,3)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (GTEST) RMS = DSQRT(DOT_PRODUCT(G,G)/(6*NPART))

      END SUBROUTINE POTENTIAL

!     ==============================================================================================

      SUBROUTINE CENTRE()
!
!     Subroutine CENTRE moves the centre of geometry to the origin.
!
      USE COMMONS

      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: C(3)

      C(:) = 0.0D0

      DO J1 =1, NPART
         C(:) = C(:) + R(J1,:)
      ENDDO

      C(:) = C(:)/NPART

      DO J1 =1, NPART
         R(J1,:) = R(J1,:) - C(:)
      ENDDO
      
      END SUBROUTINE CENTRE

!     ==============================================================================================
!
      SUBROUTINE CONTAINER(INSIDET)
!     Check nothing has moved outside the container radius 
      USE COMMONS, ONLY: NPART, R, RADIUS

      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: DIST2, RADSQ, C(3)
      LOGICAL          :: INSIDET

      RADSQ = RADIUS*RADIUS
      INSIDET = .TRUE.

      C(:) = 0.0D0

      DO J1 =1, NPART
         C(:) = C(:) + R(J1,:)
      ENDDO

      C(:) = C(:)/NPART
      
      DO J1 = 1, NPART
         DIST2 = (R(J1,1) - C(1))**2 + (R(J1,2) - C(2))**2 + (R(J1,3) - C(3))**2
         IF (DIST2 > RADSQ) THEN
            INSIDET = .FALSE.
         ENDIF
      ENDDO

      END SUBROUTINE CONTAINER

!     ==============================================================================================

      SUBROUTINE PERCOLATE(PERCT)

!     Written by DC - August, 2013 

!     This subroutine checks if the cluster is a percolating structure in a graph theoretical 
!     treatment that employs the depth-first search (DFS) algorithm to determine the number of 
!     connected components. The structure is reprsented by a graph using the cut-off parameter
!     PERCCUT for the maximum sqared distance between any two neighbours to belong to the same 
!     component. 

!     Each particle represents a vertex of the graph and an edge between two vertices exits if 
!     the corresponding squared distance is not greater than PERCCUT


!     An array MARK is used to keep track of visited vertices. All vertices are initially 
!     'unvisited' or equivalently the array MARK is initialised to 'zero'.

!     After performing the following procedure, NCOMP will be the number of connected components
!     in the graph represented by the adjacency matrix G, and for each vertex V, MARK[V] will 
!     have the index of the connected component to which V belongs.

      USE COMMONS, ONLY: NPART, R

      IMPLICIT NONE

      INTEGER :: NCOMP, J1, J2
      INTEGER :: MARK(NPART)
      DOUBLE PRECISION :: DIST2, PERCCUT2
      DOUBLE PRECISION, PARAMETER :: PERCCUT = 2.5D0
      LOGICAL :: G(NPART,NPART), PERCT

      PERCCUT2 = PERCCUT*PERCCUT

!     Construct the adjacency matrix G for the graph
      G(:,:) = .FALSE.
      PERCT  = .FALSE.

      DO J1 = 1, NPART-1
         DO J2 = J1+1, NPART
            DIST2 = (R(J1,1)-R(J2,1))**2 + (R(J1,2)-R(J2,2))**2 + (R(J1,3)-R(J2,3))**2
            IF (DIST2 <= PERCCUT2) THEN
               G(J1,J2) = .TRUE.
               G(J2,J1) = .TRUE.
            ENDIF
         ENDDO
      ENDDO

      MARK(1:NPART) = 0
      NCOMP        = 0

      DO J1 = 1, NPART
         IF (MARK(J1) == 0) THEN
            NCOMP = NCOMP + 1
            CALL DFS(J1,NCOMP,G,MARK)
         ENDIF
      ENDDO

      IF (NCOMP == 1) PERCT = .TRUE.


      END SUBROUTINE PERCOLATE

!     ==============================================================================================

      RECURSIVE SUBROUTINE DFS(V,NCOMP,G,MARK)
!     The depth-first-search algorithm adapted to find the number of connected components on a graph

      USE COMMONS, ONLY: NPART

      INTEGER :: NCOMP, V, J2
      INTEGER :: MARK(NPART)
      LOGICAL :: G(NPART,NPART)

      MARK(V) = NCOMP
      DO J2 = 1, NPART
         IF (G(V,J2)) THEN
            IF (MARK(J2) == 0) THEN
               CALL DFS(J2,NCOMP,G,MARK)
            ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE DFS
!     ===============================================================================================

SUBROUTINE TRANS_VMMC(PE, ACCPTC)

USE COMMONS

IMPLICIT NONE

INTEGER :: I, J, K, L, X, Z, Y, V, W, T, INDXP, CACCPT, OCOUNT, FCOUNT, BUCOUNT, BVCOUNT, ACCPTC
INTEGER, ALLOCATABLE :: FFAIL(:), OFAIL(:,:), C(:), PROPD(:,:), BU(:,:), BV(:,:), RN(:,:)
DOUBLE PRECISION :: RIJ(3), RIJSQ, DIS, UPE, UIPE, UJPE, M(3), RANF, DUMMY, PE, PET
DOUBLE PRECISION, PARAMETER :: DCUT = 1.5D0
LOGICAL :: GROWTH, OUTEST, FTEST, ALPROP, FRUSTEST, PRETEST, ACCPTCLUS, INSIDET, PERCT

ALLOCATE(PROPD(NPART,NPART), C(NPART), OFAIL(NPART,NPART), FFAIL(NPART), BU(NPART,NPART), BV(NPART,NPART), RN(NPART,3))

DO X = 1, NPART

  DO I = 1, NPART
    DO J = 1, NPART
      PROPD(I,J) = 0
      OFAIL(I,J) = 0 ! initialises array storing outright failed links
      BU(I,J) = 0
      BV(I,J) = 0
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


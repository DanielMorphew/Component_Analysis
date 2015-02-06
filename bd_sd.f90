      MODULE COMMONS

      IMPLICIT NONE
       
      INTEGER, PARAMETER :: NPART = 24
      DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)
      INTEGER          :: G, NC, NCSQMAX
      DOUBLE PRECISION :: ALPHA, ALPSQ, BOXL, INVINR, INVMAS, INVRPI, MASS, INERT, RCUTSQ, CNSTA, CNSTB
      DOUBLE PRECISION :: R(NPART,3), E(NPART,3), V(NPART,3), U(NPART,3), FORCE(NPART,3), GP(NPART,3)
      DOUBLE PRECISION :: RSTORE(NPART,3), ESTORE(NPART,3)
      DOUBLE PRECISION :: DT, FBDTSQ, HALFDT, DTSQB2, GRFCTR, GSFCTR, GUFCTR, X(6*NPART)
      DOUBLE PRECISION :: TMPFIX, ETA, XI, Q, RADIUS
      DOUBLE PRECISION :: DPMUSQ, YKAPPA, RSHIFT
 
      END MODULE COMMONS 
      
!     ============================================================================================== 
 
      PROGRAM MD_STOCK

      USE COMMONS

      IMPLICIT NONE

      INTEGER          :: NSTEP, NBATH, NEQ, NEVAP, IDUMP, ISTEP, ICNT, J1, J2, IVRNO
      DOUBLE PRECISION :: RHO, TMP, TMPTR, TMPRT, KETR, KERT, PRS, VLM
      DOUBLE PRECISION :: EPP, PEPP, KEPP
      DOUBLE PRECISION :: HEXT, AVHEXT, SUMHEXT
      DOUBLE PRECISION :: BOXLH, DOT, DPMU, RCUT, RCUT12, RCUT13, PE, KE
      DOUBLE PRECISION :: SUME, SUMPE, SUMKE, SUMTMPT, SUMTMPR, SUMTMP, SUMPRS
      DOUBLE PRECISION :: AVE, AVPE, AVKE, AVTMPT, AVTMPR, AVTMP, AVPRS
      LOGICAL          :: DEBUG, INSIDET 
      DOUBLE PRECISION :: DELX, ERRLIM, DFA, DFN, FM, FP
 
      OPEN (UNIT = 1, FILE = 'parameter_sd.inp', STATUS = 'UNKNOWN')

!     read input parameters

!     READ (1, *) NPART
      READ (1, *) 
      READ (1, *) MASS, INERT
      READ (1, *) 
      READ (1, *) DPMUSQ, YKAPPA, RSHIFT 
      READ (1, *) 
      READ (1, *) TMPFIX
      READ (1, *) 
      READ (1, *) RADIUS
      READ (1, *)
      READ (1, *) DT
      READ (1, *) 
      READ (1, *) NSTEP, NEQ, IDUMP
 
      CLOSE(1)
      
!     calculate from input parameters
      G      = 5*NPART - 3
      HALFDT = 0.5D0*DT
      DTSQB2 = DT*HALFDT
      FBDTSQ = 2.D0/DTSQB2
      
!     Some useful quantities
      VLM    = DFLOAT(NPART)/RHO
      BOXL   = VLM**(1.D0/3.D0)
      BOXLH  = 0.5D0*BOXL
      INVMAS = 1.D0/MASS
      INVINR = 1.D0/INERT
      RCUTSQ = RCUT*RCUT      ! for soft sphere repulsion
      RCUT12 = (1.D0/RCUT)**12
      RCUT13 = RCUT12/RCUT
      CNSTA  = 12.D0*RCUT13
      CNSTB  = 13.D0*RCUT13

      INVRPI = 1.D0/SQRT(PI)
      ALPHA = ALPHA/BOXL
      ALPSQ  = ALPHA*ALPHA 
      GUFCTR = 2.D0*PI*DPMUSQ/BOXL**3       ! Cubic box 
      GRFCTR = 4.D0*PI*GUFCTR/BOXL
      GSFCTR = 2.D0*GUFCTR

!     start with a fcc lattice configuration
!      CALL FCC(BOXL,BOXLH) 
 
      CALL EQCON()     
!     previous configuration
!     CALL READCONFIG

!     set initial linear velocities
      CALL INTVEL(TMPFIX)
      
!     set initial orientation and angular velocities
      CALL INTORN(TMPFIX)

!     calculate initial temperature
      KE = 0.D0
      DO J1 = 1, NPART
         KETR = KETR + MASS*DOT_PRODUCT(V(J1,:),V(J1,:))
         KERT = KERT + INERT*DOT_PRODUCT(U(J1,:),U(J1,:))
         KE   = 0.5*(KETR + KERT)
      ENDDO        
      TMP = KE/DFLOAT(G)
      
!     write parameters
      OPEN (UNIT = 2, FILE = 'runset.dat', STATUS = 'UNKNOWN')
      WRITE (2, *) 'Total number of Stockmayer particles = ', NPART
      WRITE (2, *) 'Mass = ', MASS, 'Moment of inertia = ', INERT
      WRITE (2, *) 'Temperature = ', TMPFIX
      WRITE (2, *) 'Density = ', RHO
      WRITE (2, *) 'Timestep = ', DT
      WRITE (2, *) 'Total number of MD steps = ', NSTEP
      WRITE (2, *) 'Number of steps in contact with the thermal bath = ', NBATH 
      WRITE (2, *) 'Number of eqilibration steps = ', NEQ
      WRITE (2, *) 'Number of steps between two dumped states = ', IDUMP
      WRITE (2, *) 'Total number of degrees of freedom = ', G
      WRITE (2, *) 'Initial temperature = ', TMP
      WRITE (2, *) 'Edge length of cubic box =', BOXL
      WRITE (2, *) 'Potential parameters: dipole stength, rcut = ', DPMUSQ, RCUT
      WRITE (2, *) 'Initial bath parameters: ETA, XI, Q = ', ETA, XI, Q

      CLOSE (2)

      OPEN (UNIT = 41, FILE = 'evap.dat', STATUS = 'UNKNOWN')
!     start simulation
      ISTEP   = 0 
      ICNT    = 0
      SUME    = 0.D0
      SUMHEXT = 0.D0
      SUMKE   = 0.D0
      SUMPE   = 0.D0
      SUMTMPT = 0.D0
      SUMTMPR = 0.D0
      SUMTMP  = 0.D0
      SUMPRS  = 0.D0
      NEVAP   = 0

!      CALL CHECKDRVTS()
!      STOP

!      CALL FORQUE(PE,FORCE,GP)

      DO WHILE (ISTEP < NSTEP)            

         ICNT  = ICNT + 1
         ISTEP = ISTEP + 1
         IF (MOD(ISTEP,100) == 0) THEN
            RSTORE(:,:) = R(:,:)
            ESTORE(:,:) = E(:,:)
         ENDIF

      CALL FORQUE(PE,FORCE,GP)
!     advance positions, orientations, and their time derivatives
        
!        CALL MVNVEA()
!        CALL MVNVTA()

         CALL MOVE(DT,TMPFIX)          
110     CALL CENTRE()

        CALL CONTAINER(INSIDET)

        IF (.NOT. INSIDET) THEN

           IF (ISTEP < 100) THEN
              PRINT *, 'EVAPORATION TOO SOON'
              STOP
           ENDIF

           NEVAP = NEVAP + 1
           WRITE(41, *) NEVAP
           
           R = RSTORE; E = ESTORE

!     set initial linear velocities
!           CALL INTVEL(TMPFIX)

!     set initial orientation and angular velocities
!           CALL INTORN(TMPFIX)

           CALL FORQUE(PE,FORCE,GP)

!           CALL MVNVEA()
!           CALL MVNVTA()

           CALL MOVE(DT,TMPFIX)          
           GO TO 110

        ENDIF


!     calculate forces and gorques at time t
       
!        CALL FORQUE(PE,FORCE,GP)  

!        CALL MVNVEB(KETR,KERT)
!        CALL MVNVTB(KETR,KERT)
      
!     calculate temperatures, pressure, and energies at time t

!         TMPTR    = KETR/DFLOAT(3*NPART - 3)
!         TMPRT    = KERT/DFLOAT(2*NPART)           
!         KE       = KETR + KERT
!         TMP      = KE/DFLOAT(G)            ! Here KE is twice the actual KE 
!         KE       = 0.5D0*KE
         PEPP     = PE/DFLOAT(NPART)
!         print *, pepp
!         KEPP     = KE/DFLOAT(NPART)
!         EPP      = PEPP + KEPP
!         PRS      = (KETR + W)/(3.D0*VLM)
!         SUMKE    = SUMKE + KEPP
         SUMPE    = SUMPE + PEPP
!         SUMTMPT  = SUMTMPT + TMPTR
!         SUMTMPR  = SUMTMPR + TMPRT
!         SUMTMP   = SUMTMP  + TMP
!         SUMPRS   = SUMPRS + PRS
 
    
!         HEXT    = PE + KE + 0.5D0*Q*XI*XI + TMPFIX*DFLOAT(G)*ETA
!         HEXT    = HEXT/DFLOAT(NPART)
!         SUMHEXT = SUMHEXT + HEXT
!         AVHEXT  = SUMHEXT/DFLOAT(ICNT)

!         IF (MOD(ISTEP,IDUMP) == 0) THEN
!            OPEN(UNIT = 15, FILE = 'hext.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!            WRITE (15, *) HEXT, AVHEXT
!            CLOSE (UNIT = 15, STATUS = 'KEEP')
!         ENDIF
 
         IF (MOD(ISTEP,IDUMP) == 0) THEN 

!            AVE     = SUME/DFLOAT(ICNT)
!            AVKE    = SUMKE/DFLOAT(ICNT)
            AVPE    = SUMPE/DFLOAT(ICNT)
!            AVTMPT  = SUMTMPT/DFLOAT(ICNT)
!            AVTMPR  = SUMTMPR/DFLOAT(ICNT)
!            AVTMP   = SUMTMP/DFLOAT(ICNT)
!            AVPRS   = SUMPRS/DFLOAT(ICNT)

!            OPEN (UNIT = 3, FILE ='energy.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!            WRITE (3,*) KEPP, PEPP, EPP
!            CLOSE (UNIT = 3, STATUS = 'KEEP')
            OPEN (UNIT = 4, file = 'ave.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (4,*) PEPP, AVPE
!            print *, avpe 
            CLOSE (UNIT = 4, STATUS = 'KEEP')
!            OPEN (UNIT = 14, FILE = 'tmp.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!            WRITE (14,*) AVTMPT, AVTMPR, AVTMP
!            CLOSE (UNIT= 14, STATUS = 'KEEP')  

         ENDIF

         IF (ISTEP > NEQ .AND. MOD(ISTEP,IDUMP) ==  0) THEN

            OPEN (UNIT = 7,  FILE = 'pos.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            OPEN (UNIT = 8,  FILE = 'ortn.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!            OPEN (UNIT = 9,  FILE = 'linvel.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!            OPEN (UNIT = 10, FILE = 'angvel.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
       
            DO J1 = 1, NPART
               WRITE (7,*)  R(J1,1), R(J1,2), R(J1,3)
               WRITE (8,*)  E(J1,1), E(J1,2), E(J1,3)
!               WRITE (9,*)  V(J1,1), V(J1,2), V(J1,3)
!               WRITE (10,*) U(J1,1), U(J1,2), U(J1,3)
            ENDDO
        
            CLOSE (UNIT = 7,  STATUS = 'KEEP')  
            CLOSE (UNIT = 8,  STATUS = 'KEEP')  
!            CLOSE (UNIT = 9,  STATUS = 'KEEP')  
!            CLOSE (UNIT = 10, STATUS = 'KEEP')  
            
         ENDIF

         IF (ISTEP == NBATH .OR. ISTEP == NEQ) THEN
            ICNT    = 0
            SUME    = 0.D0
            SUMHEXT = 0.D0
            SUMKE   = 0.D0
            SUMPE   = 0.D0
            SUMTMPT = 0.D0
            SUMTMPR = 0.D0
            SUMTMP  = 0.D0
            SUMPRS  = 0.D0
         ENDIF
      
      ENDDO 

      OPEN (UNIT = 17,  FILE = 'finalpos.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
      OPEN (UNIT = 18,  FILE = 'finalortn.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!      OPEN (UNIT = 19,  FILE = 'finallinvel.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
!      OPEN (UNIT = 20, FILE = 'finalangvel.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
       
      DO J1 = 1, NPART
         WRITE (17,*)  R(J1,1), R(J1,2), R(J1,3)
         WRITE (18,*)  E(J1,1), E(J1,2), E(J1,3)
!         WRITE (19,*)  V(J1,1), V(J1,2), V(J1,3)
!         WRITE (20,*) U(J1,1), U(J1,2), U(J1,3)
      ENDDO
         
      CLOSE (UNIT = 17,  STATUS = 'KEEP')  
      CLOSE (UNIT = 18,  STATUS = 'KEEP')  
!      CLOSE (UNIT = 19,  STATUS = 'KEEP')  
!      CLOSE (UNIT = 20, STATUS = 'KEEP')  
     
      CALL VIEWCONFIG() 

      END PROGRAM MD_STOCK
      
!     ======================================SUBROUTINES=============================================

      SUBROUTINE FORQUE(PE,FORCE,GP)

      USE COMMONS, ONLY: NPART, R, E, DPMUSQ, YKAPPA, RSHIFT

      IMPLICIT NONE

      INTEGER          :: J1, J2
      DOUBLE PRECISION :: RIJ(3), NR(3), RSQ, R2, R12, R6, R3, R4, ABSRIJ, RS1(3), RS2(3), RSS(3)
      DOUBLE PRECISION :: EI(3), EJ(3), FIJ(3), GI(3), GJ(3)
      DOUBLE PRECISION :: ALP, BET, GAM, DADEI(3), DADEJ(3), DBDEI(3), DBDEJ(3)
      DOUBLE PRECISION :: FIJN, FIJEI, FIJEJ, VR, VA, VB, VG, DPFCT, EXPFCT
      DOUBLE PRECISION :: PE, FORCE(NPART,3), GP(NPART,3)

      PE = 0.D0
      FORCE(:,:) = 0.D0
      GP(:,:) = 0.D0

      DO J1 = 1, NPART - 1
         EI(:) = E(J1,:)
         DO J2 = J1 + 1, NPART
            EJ(:)  = E(J2,:)
            RIJ(:) = R(J1,:) - R(J2,:)
            RSQ    = DOT_PRODUCT(RIJ(:),RIJ(:))
            ABSRIJ = DSQRT(RSQ)
            R2     = 1.D0/RSQ
            EXPFCT = EXP(-YKAPPA*(ABSRIJ - 1.D0))
            PE     = PE + EXPFCT/ABSRIJ
            FIJ(:) =-EXPFCT*R2*(YKAPPA + 1.D0/ABSRIJ)*RIJ(:)
            FORCE(J1,:) = FORCE(J1,:) - FIJ(:)
            FORCE(J2,:) = FORCE(J2,:) + FIJ(:)
!     Dipolar part
            RS1(:) = R(J1,:) + RSHIFT*EI(:)
            RS2(:) = R(J2,:) + RSHIFT*EJ(:)
            RSS(:) = RS1(:) - RS2(:)
            R2     = DOT_PRODUCT(RSS(:),RSS(:))
            ABSRIJ = DSQRT(R2)
            NR(:)  = RSS(:)/ABSRIJ
            R2     = 1.D0/R2
            R3     = R2/ABSRIJ
            R4     = R2*R2
            ALP    = DOT_PRODUCT(NR(:),EI(:))
            BET    = DOT_PRODUCT(NR(:),EJ(:))
            GAM    = DOT_PRODUCT(EI(:),EJ(:))
            DPFCT  = 3.D0*DPMUSQ
            PE     = PE + DPFCT*R3*(GAM/3.D0 - ALP*BET)

            VR     =-DPFCT*R4*(GAM - 3.D0*ALP*BET)
            VA     =-DPFCT*BET*R3
            VB     =-DPFCT*ALP*R3
            VG     = DPFCT*R3/3.D0
            FIJN   = VR - (VA*ALP + VB*BET)/ABSRIJ
            FIJEI  = VA/ABSRIJ
            FIJEJ  = VB/ABSRIJ
            FIJ(:) = FIJN*NR(:) + FIJEI*EI(:) + FIJEJ*EJ(:)

            DADEI(:) = NR(:) + RSHIFT*EI(:)/ABSRIJ - ALP*R2*RSHIFT*RSS(:) 
            DADEJ(:) =       - RSHIFT*EI(:)/ABSRIJ + ALP*R2*RSHIFT*RSS(:) 

            DBDEI(:) =         RSHIFT*EJ(:)/ABSRIJ - BET*R2*RSHIFT*RSS(:) 
            DBDEJ(:) = NR(:) - RSHIFT*EJ(:)/ABSRIJ + BET*R2*RSHIFT*RSS(:) 

            GI(:)  = DPFCT*(RSHIFT*GAM*R4*NR(:) - R3*EJ(:)/3.D0)      &
     &              -3.D0*DPFCT*RSHIFT*R4*ALP*BET*NR(:) + DPFCT*R3*(BET*DADEI(:) + ALP*DBDEI(:))
            GJ(:)  =-DPFCT*(RSHIFT*GAM*R4*NR(:) + R3*EI(:)/3.D0)      &
     &              +3.D0*DPFCT*RSHIFT*R4*ALP*BET*NR(:) + DPFCT*R3*(BET*DADEJ(:) + ALP*DBDEJ(:))

            FORCE(J1,:) = FORCE(J1,:) - FIJ(:)
            FORCE(J2,:) = FORCE(J2,:) + FIJ(:)

            GP(J1,:) = GP(J1,:) + GI(:)
            GP(J2,:) = GP(J2,:) + GJ(:)

         ENDDO

      ENDDO

!      DO J1 = 1, NPART
!         GP(J1,:)    = GP(J1,:) - DOT_PRODUCT(GP(J1,:),E(J1,:))*E(J1,:)
!      ENDDO

      END SUBROUTINE FORQUE

!     ============================================================================================== 

      SUBROUTINE MVNVEA()

      USE COMMONS, ONLY: NPART, INVMAS, INVINR, DT, DTSQB2, FBDTSQ, HALFDT, R, E, U, V, FORCE, GP
      
      IMPLICIT NONE
      
      INTEGER          :: J1
      DOUBLE PRECISION :: UH(3), UHSQ, FLMD, DFLMD, DLAMBDA, NRMFCT, LAMBDA
      DOUBLE PRECISION, PARAMETER :: ER = 1.0D-10
      LOGICAL          :: READY
      
      DO J1 = 1, NPART

         R(J1,:) = R(J1,:) + DT*V(J1,:) + DTSQB2*FORCE(J1,:)*INVMAS
         V(J1,:) = V(J1,:) + HALFDT*FORCE(J1,:)*INVMAS

         UH(:)  = U(J1,:) + HALFDT*GP(J1,:)*INVINR
         UHSQ   = DOT_PRODUCT(UH(:),UH(:))
         READY  = .FALSE.
         LAMBDA = UHSQ 
         DO WHILE (.NOT. READY)
            FLMD    = LAMBDA*LAMBDA + FBDTSQ*(LAMBDA + UHSQ)
            DFLMD   = 2.D0*LAMBDA + FBDTSQ
            DLAMBDA = FLMD/DFLMD
            LAMBDA  = LAMBDA - DLAMBDA               
            IF (ABS(DLAMBDA) < ER) READY = .TRUE.
         ENDDO
         U(J1,:) = UH(:) + HALFDT*LAMBDA*E(J1,:)
         E(J1,:) = E(J1,:) + DT*U(J1,:)
         NRMFCT = 1.D0/DSQRT(DOT_PRODUCT(E(J1,:),E(J1,:)))
         E(J1,:) = NRMFCT*E(J1,:)

      ENDDO

      END SUBROUTINE MVNVEA

!     ==============================================================================================

      SUBROUTINE MVNVEB(KETR,KERT)
      
      USE COMMONS, ONLY: NPART, E, U, V, FORCE, GP, MASS, INERT, HALFDT, INVMAS, INVINR, HALFDT

      IMPLICIT NONE
      
      INTEGER          :: J1
      DOUBLE PRECISION :: DOT, KETR, KERT

      KETR  = 0.D0
      KERT  = 0.D0
      
      DO J1 = 1, NPART
         V(J1,:) = V(J1,:) + HALFDT*FORCE(J1,:)*INVMAS 

         DOT     = DOT_PRODUCT(E(J1,:),U(J1,:))
         U(J1,:) = U(J1,:) + HALFDT*GP(J1,:)*INVINR - DOT*E(J1,:)

         KETR = KETR + DOT_PRODUCT(V(J1,:),V(J1,:))
         KERT = KERT + DOT_PRODUCT(U(J1,:),U(J1,:))
      ENDDO
      
      KETR = MASS*KETR
      KERT = INERT*KERT

      END SUBROUTINE MVNVEB

!     ==============================================================================================

      SUBROUTINE MVNVTA()

      USE COMMONS, ONLY: NPART, MASS, INERT, INVMAS, INVINR, DT, DTSQB2, FBDTSQ, HALFDT, R, E, U, V, &
     &                   FORCE, GP, ETA, XI, Q, G, TMPFIX

      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: UH(3), UHSQ, FLMD, DFLMD, DLAMBDA, NRMFCT, LAMBDA, FCTR, KEHT, XINEW
      DOUBLE PRECISION, PARAMETER :: ER = 1.0D-10
      LOGICAL          :: READY

      FCTR = 1.D0 - HALFDT*XI
      KEHT = 0.D0

      DO J1 = 1, NPART

         R(J1,:) = R(J1,:) + DT*V(J1,:) + DTSQB2*INVMAS*FORCE(J1,:)
         V(J1,:) = FCTR*V(J1,:) + HALFDT*INVMAS*FORCE(J1,:)
         KEHT    = KEHT + MASS*DOT_PRODUCT(V(J1,:),V(J1,:))

         UH(:)  = FCTR*U(J1,:) + HALFDT*INVINR*GP(J1,:)
         UHSQ   = DOT_PRODUCT(UH(:),UH(:))
         READY  = .FALSE.
         LAMBDA = UHSQ
         DO WHILE (.NOT. READY)
            FLMD    = LAMBDA*LAMBDA + FBDTSQ*(LAMBDA + UHSQ)
            DFLMD   = 2.D0*LAMBDA + FBDTSQ
            DLAMBDA = FLMD/DFLMD
            LAMBDA  = LAMBDA - DLAMBDA
            IF (ABS(DLAMBDA) < ER) READY = .TRUE.
         ENDDO
         U(J1,:) = UH(:) + HALFDT*LAMBDA*E(J1,:)
         KEHT    = KEHT + INERT*DOT_PRODUCT(U(J1,:),U(J1,:))
         E(J1,:) = E(J1,:) + DT*U(J1,:)
         NRMFCT = 1.D0/DSQRT(DOT_PRODUCT(E(J1,:),E(J1,:)))
         E(J1,:) = NRMFCT*E(J1,:)

      ENDDO

      XINEW = XI + (KEHT - DFLOAT(G)*TMPFIX)*DT/Q
      ETA   = ETA + HALFDT*(XI + XINEW)
      XI    = XINEW

      END SUBROUTINE MVNVTA

!     ==============================================================================================

      SUBROUTINE MVNVTB(KETR,KERT)

      USE COMMONS, ONLY: NPART, E, U, V, FORCE, GP, MASS, INERT, HALFDT, INVMAS, INVINR, HALFDT, ETA, XI, Q

      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: DOT, KETR, KERT, FCTR

      KETR  = 0.D0
      KERT  = 0.D0
      FCTR  = 1.D0/(1.D0 + HALFDT*XI)

      DO J1 = 1, NPART
         V(J1,:) = FCTR*(V(J1,:) + HALFDT*INVMAS*FORCE(J1,:))

         DOT     = DOT_PRODUCT(E(J1,:),U(J1,:))
         U(J1,:) = FCTR*(U(J1,:) + HALFDT*INVINR*GP(J1,:) - DOT*E(J1,:))

         KETR = KETR + DOT_PRODUCT(V(J1,:),V(J1,:))
         KERT = KERT + DOT_PRODUCT(U(J1,:),U(J1,:))

      ENDDO

      KETR = MASS*KETR
      KERT = INERT*KERT

      END SUBROUTINE MVNVTB

!     ==============================================================================================

      SUBROUTINE FCC(BOXL,BOXLH)

      USE COMMONS, ONLY: NPART, R, E

      IMPLICIT NONE

      INTEGER :: J1, IX, IY, IZ, M, IREF, NUC
      DOUBLE PRECISION :: BOXL, BOXLH, UCL, UCLH, RROOT3
      
!     Number of unit cells

      J1 = 1
      DO WHILE (4*J1**3 < NPART)
         J1 = J1 + 1
      ENDDO

      NUC = J1
!     Unit cell length
      UCL     = BOXL/DFLOAT(NUC) 
      UCLH    = 0.5D0*UCL
      RROOT3  = 1.D0/DSQRT(3.D0)

!     Build the unit cell
!     Sublattice A
      R(1,:) = 0.D0
      E(1,:) = RROOT3

!     Sublattice B
      R(2,:) = (/UCLH, UCLH, 0.D0/)
      E(2,:) = (/RROOT3,-RROOT3,-RROOT3/)

!     Sublattice C
      R(3,:) = (/0.D0, UCLH, UCLH/)
      E(3,:) = (/-RROOT3, RROOT3,-RROOT3/)
!     Sublattice D
      R(4,:) = (/UCLH, 0.D0, UCLH/)
      E(4,:) = (/-RROOT3,-RROOT3, RROOT3/)

!     Construct the lattice from the unit cell 
      M = 0
      DO IZ = 1, NUC
         DO IY = 1, NUC
            DO IX = 1, NUC
               DO IREF = 1, 4
                  R(IREF+M,1) = R(IREF,1) + UCL*DFLOAT(IX-1)
                  R(IREF+M,2) = R(IREF,2) + UCL*DFLOAT(IY-1)
                  R(IREF+M,3) = R(IREF,3) + UCL*DFLOAT(IZ-1)
                  E(IREF+M,1:3) = E(IREF,1:3)
               ENDDO
               M = M + 4
            ENDDO
         ENDDO
      ENDDO

!     Shift the centre of the box to the origin
      R(:,:) = R(:,:) - BOXLH

      END SUBROUTINE FCC

!     ==============================================================================================

      SUBROUTINE READCONFIG

      USE COMMONS, ONLY: NPART, R, E, U, V

      IMPLICIT NONE

      INTEGER :: J1

      OPEN (UNIT = 33, FILE = 'finalpos1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 34, FILE = 'finalortn1.dat', STATUS = 'UNKNOWN')
!      OPEN (UNIT = 35, FILE = 'finallinvel.dat', STATUS = 'UNKNOWN')
!      OPEN (UNIT = 36, FILE = 'finalangvel.dat', STATUS = 'UNKNOWN')

      DO J1 = 1, NPART
         READ(33,*) R(J1,1), R(J1,2), R(J1,3)
         READ(34,*) E(J1,1), E(J1,2), E(J1,3)
!         READ(35,*) V(J1,1), V(J1,2), V(J1,3)
!         READ(36,*) U(J1,1), U(J1,2), U(J1,3)
      ENDDO

      CLOSE(UNIT = 33)
      CLOSE(UNIT = 34)
      CLOSE(UNIT = 35)
      CLOSE(UNIT = 36)

      END SUBROUTINE READCONFIG

!     ==============================================================================================

      SUBROUTINE INTVEL(TMPFIX)

      USE COMMONS, ONLY: V, INVMAS, NPART 

      IMPLICIT NONE

      INTEGER :: J1, J2
      DOUBLE PRECISION :: SDVEL, SUMV(3), AVV(3), GAUSS, DUMMY, TMPFIX
      
      SDVEL = DSQRT(TMPFIX*INVMAS)
      
      DO J1 = 1, NPART
         DO J2 = 1, 3
            V(J1,J2) = SDVEL*GAUSS(DUMMY)
         ENDDO
      ENDDO
      
      SUMV(:) = 0.D0
      DO J1 = 1, NPART
         DO J2 = 1, 3
            SUMV(J2) = SUMV(J2) + V(J1,J2)
         ENDDO
      ENDDO

      AVV(:) = SUMV(:)/DFLOAT(NPART)

      DO J1 = 1, NPART
         DO J2 = 1, 3
            V(J1,J2) = V(J1,J2) - AVV(J2)
         ENDDO
      ENDDO

      END SUBROUTINE INTVEL

!     ==============================================================================================

      SUBROUTINE INTORN(TMPFIX) 
      
      USE COMMONS, ONLY: NPART, E, U, INERT

      IMPLICIT NONE 

      INTEGER :: J1
      DOUBLE PRECISION :: MEAN, XI, XISQ, XI1, XI2, DOT, USQ, NRMFCR, RANF, DUMMY, TMPFIX

      MEAN = 2.D0*TMPFIX/INERT

      DO J1 = 1, NPART

!     Set direction of the angular velocity 
!     Choose a random vector in space

         XISQ = 1.D0

1000     IF (XISQ >= 1.D0) THEN
            XI1  = RANF(DUMMY)*2.D0 - 1.D0
            XI2  = RANF(DUMMY)*2.D0 - 1.D0
            XISQ = XI1*XI1 + XI2*XI2

            GO TO 1000

         ENDIF 

         XI      = DSQRT(1.D0 - XISQ)
         U(J1,1) = 2.D0*XI1*XI
         U(J1,2) = 2.D0*XI2*XI
         U(J1,3) = 1.D0 - 2.D0*XISQ

!     Constrain the vector to be perpendicular to the molecules
         DOT     = DOT_PRODUCT(U(J1,:),E(J1,:))
         U(J1,:) = U(J1,:) - DOT*E(J1,:)

!     Renormalize 
         USQ     = DOT_PRODUCT(U(J1,:),U(J1,:))
         NRMFCR  = 1.D0/DSQRT(USQ)
         U(J1,:) = NRMFCR*U(J1,:)

!     Choose the magnitude of the angular velocity 
         USQ   =-MEAN*DLOG(RANF(DUMMY))
         U(J1,:) = DSQRT(USQ)*U(J1,:)

      ENDDO

      END SUBROUTINE INTORN

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

      DOUBLE PRECISION FUNCTION RANF(DUMMY)
!     Drawing a uniform random variate between 0 and 1 
      IMPLICIT NONE
      INTEGER, PARAMETER :: L = 1029, C = 221591, M =1048576
      INTEGER :: SEED
      DOUBLE PRECISION :: DUMMY
      SAVE SEED
      DATA SEED / 0 /

      SEED = MOD(SEED*L + C, M)
      RANF = DFLOAT(SEED)/DFLOAT(M)

      END FUNCTION RANF

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

      SUBROUTINE EQCON()

      USE COMMONS

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: P(3), RM(3,3)

!      OPEN (UNIT=23, FILE='initialconfig.dat', STATUS='UNKNOWN')
      OPEN (UNIT=13, FILE='initialpos.dat', STATUS='UNKNOWN')
      OPEN (UNIT=14, FILE='initialortn.dat', STATUS='UNKNOWN')

      DO I = 1, NPART  
!         READ(23, *) R(I,1), R(I,2), R(I,3)
         READ(13, *) R(I,1), R(I,2), R(I,3)
         READ(14, *) E(I,1), E(I,2), E(I,3)
      ENDDO

!      DO I = 1, NPART
!         READ(23, *) P(1), P(2), P(3)
!         CALL ROTMAT(P,RM)
!         E(I,:) = RM(:,3)
!      ENDDO

!      CLOSE (23)
      CLOSE (13)
      CLOSE (14)

      END SUBROUTINE EQCON

!     ==============================================================================================

      SUBROUTINE ROTMAT(P, RM)

      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, CT, ST, I3(3,3), E(3,3), RM(3,3)

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

      THETA2 = DOT_PRODUCT(P,P)

      IF (THETA2 == 0.D0) THEN
         RM(:,:) = I3(:,:)
      ELSE
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         RM      = I3(:,:) + (1.D0-CT)*MATMUL(E(:,:),E(:,:)) + ST*E(:,:)
      ENDIF

      END SUBROUTINE ROTMAT

!     ==============================================================================================

      SUBROUTINE VIEWCONFIG()

      USE COMMONS

      IMPLICIT NONE

      INTEGER :: J1, J2
      DOUBLE PRECISION :: SR(3)

      OPEN (UNIT = 31, FILE = 'finalconfig.xyz', STATUS = 'UNKNOWN')

      WRITE(31,'(I6)') NPART*2
      WRITE(31, *)
      DO J1 = 1, NPART
         WRITE(31,'(A5,1X,3F20.10)') 'O ', R(J1,1), R(J1,2), R(J1,3)
         SR(:) = R(J1,:) + RSHIFT*E(J1,:)
         WRITE(31,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'C', SR(1), SR(2), SR(3), &
                'atom_vector', E(J1,1), E(J1,2), E(J1,3)
      ENDDO

      CLOSE (UNIT = 31)

      END SUBROUTINE VIEWCONFIG

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
            RETURN
         ENDIF
      ENDDO

      END SUBROUTINE CONTAINER

!     ==============================================================================================

      SUBROUTINE CHECKDRVTS()

      USE COMMONS, ONLY: NPART, R, E

      IMPLICIT NONE

      INTEGER          :: IVRNO1, J1, J2
      DOUBLE PRECISION :: PE, FORCE(NPART,3), GP(NPART,3), FM, FP, DFA, DFN
      DOUBLE PRECISION, PARAMETER :: ERRLIM = 1.D-06, DELX = 1.D-06

!     Checks gradients

      IVRNO1 = 0

      DO J1 = 1, NPART
   
         DO J2 = 1, 3

            IVRNO1 = IVRNO1 + 1

            WRITE(*, *) IVRNO1

            R(J1,J2) = R(J1,J2) - DELX

            CALL FORQUE(FM,FORCE,GP)
!           WRITE(*, *) 'Energy minus = ', FM

            R(J1,J2) = R(J1,J2) + 2.D0*DELX
            CALL FORQUE(FP,FORCE,GP)
!           WRITE(*, *) 'Energy plus  = ', FP

            R(J1,J2) = R(J1,J2) - DELX
            CALL FORQUE(PE,FORCE,GP)
            DFN = (FP - FM) / (2.D0*DELX)
            DFA =-FORCE(J1,J2)

            WRITE(*, *) 'Gradient numerical  = ', DFN
            WRITE(*, *) 'Gradient analytical = ', DFA

            IF (ABS(DFN - DFA) > ERRLIM) WRITE(*, *) 'Error:', IVRNO1, DFN, DFA, ABS(DFN-DFA)

         ENDDO

         DO J2 = 1, 3

            IVRNO1 = IVRNO1 + 1

            WRITE(*, *) IVRNO1

            E(J1,J2) = E(J1,J2) - DELX

            CALL FORQUE(FM,FORCE,GP)
!           WRITE(*, *) 'Energy minus = ', FM

            E(J1,J2) = E(J1,J2) + 2.D0*DELX
            CALL FORQUE(PE,FORCE,GP)
!           WRITE(*, *) 'Energy plus  = ', FP

            E(J1,J2) = E(J1,J2) - DELX
            CALL FORQUE(FP,FORCE,GP)
            DFN = (FP - FM) / (2.D0*DELX)
            DFA =-GP(J1,J2)

            WRITE(*, *) 'Gradient numerical  = ', DFN
            WRITE(*, *) 'Gradient analytical = ', DFA

            IF (ABS(DFN - DFA) > ERRLIM) WRITE(*, *) 'Error:', IVRNO1, DFN, DFA, ABS(DFN-DFA)
         
         ENDDO

      ENDDO

      END SUBROUTINE CHECKDRVTS

!     ==============================================================================================

      SUBROUTINE MOVE (DT,TMPFIX)

      USE COMMONS, ONLY: NPART, R, E, FORCE, GP

      IMPLICIT NONE
      INTEGER   ::  J1, J2
      DOUBLE PRECISION :: DT, TMPFIX
      DOUBLE PRECISION :: XI(3), GAUSS, DUMMY, EI(3), GI(3), TI(3), A(3), B(3)

      DO J1 = 1, NPART
!     Calculate uncorrelated random normal deviates with zero mean and variance 2*DT
         DO J2 = 1, 3
            XI(J2) = GAUSS (DUMMY)*SQRT(2.D0*DT)
         ENDDO

         R(J1,:) = R(J1,:) + DT/TMPFIX*FORCE(J1,:)+ XI(:)

         DO J2 = 1, 3
            XI(J2) = GAUSS (DUMMY)*SQRT(6.D0*DT)
         ENDDO
         EI(:) = E(J1,:)
         GI(:) = GP(J1,:)
         TI(:) = (/EI(2)*GI(3) - EI(3)*GI(2), EI(3)*GI(1) - EI(1)*GI(3), EI(1)*GI(2) - EI(2)*GI(1)/)
         A(:)  = (/TI(2)*EI(3) - TI(3)*EI(2), TI(3)*EI(1) - TI(1)*EI(3), TI(1)*EI(2) - TI(2)*EI(1)/)
         B(:)  = (/XI(2)*EI(3) - XI(3)*EI(2), XI(3)*EI(1) - XI(1)*EI(3), XI(1)*EI(2) - XI(2)*EI(1)/)

         E(J1,:) = EI(:) + 3.D0*DT/TMPFIX*A(:)+ B(:)

         E(J1,:) = E(J1,:)/DSQRT(DOT_PRODUCT(E(J1,:),E(J1,:)))
      ENDDO

      END SUBROUTINE MOVE

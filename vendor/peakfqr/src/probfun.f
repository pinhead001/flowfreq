C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C     SET OF SUBROUTINES TO DEAL WITH PROBABILITY FUNCTIONS
C
C     COPYRIGHT, TIMOTHY A. COHN, 1995
C     PROPERTY OF US GOVERNMENT, GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C     
C
C     AUTHOR.......TIM COHN
C           DATE.....JANUARY 5, 1995
C
C     Minor modifications and additions
C     by John England, USACE john.f.england@usace.army.mil
C
C     modifications by John England denoted with C JFE
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C    07-SEP-2017 extended H-S routine arrays to 20000
C          
C
C
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
CC
C     N.B.  IMSL SUBROUTINES ARE USED:
C                 DGAMMA
C                 DGAMDF
C                 DNORDF
C                 DNORIN
C                 DCHIIN
C
C     N.B.  CDFLIB ROUTINES ARE USED
C                 CDFBET
C
C     ROUTINES AVAILABLE
C           STANDARD NORMAL DISTRIBUTION
C                 FP_Z_PDF(X)
C                 FP_Z_CDF(X)
C                 FP_Z_ICDF(P)
C                 SP_Z(X,F,PHI,LOGPHI,A,G)
C
C
C           GAMMA DISTRIBUTION
C                 FP_G1_ICDF(P,ALPHA)
C                 FP_G2_ICDF(P,PARMS)  -- PARMS=(ALPHA,BETA)
C                 FP_G3_ICDF(P,PARMS)  -- PARMS=(TAU,ALPHA,BETA)
C                 FP_G1_MOM_TRA(ALPHA,T,K) 
C
C           BETA DISTRIBUTION
C                 FP_B_ICDF(P,ALPHA,BETA)
C
C     NAMING CONVENTIONS:
C
C     1)    TYPE OF ROUTINE
C                 FP_         =     DOUBLE PRECISION FUNCTIONS
C                 SP_         =     SUBROUTINES
C
C     2)    NEXT COMES LETTERS DEFINING THE DISTRIBUTION:
C
C                 Z_          =     PHI(X)                              STANDARD NORMAL DISTRIBUTION
C                 N_          =     N(MU,SIGMA)                   2-P NORMAL DISTRIBUTION 
C                 G1_         =     GAMMA(X,ALPHA)                1-P GAMMA DISTRIBUTION 
C                 G2_         =     GAMMA(X,ALPHA,BETA)           2-P GAMMA DISTRIBUTION 
C                 G3_         =     GAMMA(X,ALPHA,BETA,TAU) 3-P GAMMA DISTRIBUTION 
C
C     3)    NEXT COMES FUNCTION TYPE:
C                 
C                 PDF_        =     PROBABILITY DENSITY FUNCTION
C                 CDF_        =     CUMULATIVE DENSITY FUNCTION
C                 ICDF_       =     INVERSE CDF
C                 MOM_        =     MOMENTS
C                 HAZ_        =     F/(1-F) 
C                 RAT_        =     F/F
C
C     4)    NEXT COMES WHETHER LOGS ARE USED
C
C                 LN_               =     USE NATURAL LOG
C
C     5)    NEXT COMES FUNCTION TYPE
C
C                 JAC_  =     JACOBIAN (VECTOR OF FIRST DERIVATIVES WRT PARAMS)
C                 HESS_ =     HESSIAN OF PDF (MATRIX OF SECOND DERIVATIVES WRT PARAMS)
C                 
C     6)    NEXT COMES TRUNCATION
C
C                 TRA_  =     DISTRIBUTION TRUNCATED ABOVE (I.E. ALL VALUES "LESS THAN")
C                 TRB_  =     DISTRIBUTION TRUNCATED BELOW
C 
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_Z_PDF(X)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE PDF OF STANDARD NORMAL DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      
      CALL SP_Z(X,F,PHI,XLNPHI,RAT,HAZ)
      
      FP_Z_PDF    =     F
      
      RETURN
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_Z_CDF(X)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE CDF OF STANDARD NORMAL DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      
      CALL SP_Z(X,F,PHI,XLNPHI,RAT,HAZ)
      
      FP_Z_CDF    =     PHI
      
      RETURN
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_Z_ICDF(P)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE INVERSE CDF OF STANDARD NORMAL DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      
      FP_Z_ICDF   =     DNORIN(P)
      
      RETURN
      END



C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G1_PDF(X,ALPHA)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE PDF OF 1-P GAMMA DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      use iso_fortran_env, only: real64
      implicit none
  
      REAL(real64) X, ALPHA, DLNGAM, TEMP

      
      IF(X .GT. 0.D0 .AND. ALPHA .GT. 0.D0) THEN
            IF((ALPHA+X) .LT. 100.D0) THEN
                  FP_G1_PDF   =     X**(ALPHA-1.D0) * 
     1                                 EXP(-X)/gamma(ALPHA)
            ELSE
                  TEMP        =     (ALPHA-1.D0) * LOG(X) - X -  
     1                                 DLNGAM(ALPHA)
                  FP_G1_PDF   =     EXP(TEMP)
            ENDIF 
      ELSE
            FP_G1_PDF   =     0.D0
      ENDIF
      
      RETURN
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G2_PDF(X,PARMS)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE PDF OF 2-P GAMMA DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DIMENSION PARMS(2)
      
      ALPHA =     PARMS(1)
      BETA  =     PARMS(2)
      
      FP_G2_PDF   =     FP_G1_PDF(X/BETA,ALPHA)/ABS(BETA)
      
      RETURN
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G3_PDF(X,PARMS)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE PDF OF 3-P GAMMA DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DIMENSION PARMS(3)
      
      TAU         =     PARMS(1)
      ALPHA =     PARMS(2)
      BETA  =     PARMS(3)
      
      FP_G3_PDF   =     FP_G2_PDF(X-TAU,PARMS(2))
      
      RETURN
      END

!*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~

      DOUBLE PRECISION FUNCTION FP_G1_CDF(X,ALPHA)
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
! 
!     PROGRAM TO COMPUTE CDF OF 1-P GAMMA DISTRIBUTION
!
!     AUTHOR.....TIM COHN
!     DATE.......DECEMBER 28, 1994
!
!     MODIFIED SEPTEMBER 4, 2024 (SAS) - FOR 128 BIT X IN DGAMDF
!
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      use iso_fortran_env, only: real64, real128
      implicit none

      REAL(real64) X, ALPHA, DGAMDF
      REAL(real128) XQ 
      
      XQ = REAL(X, real128) !DGAMDF needs 128 bit X
      
      FP_G1_CDF   =     DGAMDF(XQ,ALPHA)
      
      RETURN
      END
      
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~

      DOUBLE PRECISION FUNCTION FP_G2_CDF(X,PARMS)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE CDF OF 2-P GAMMA DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DIMENSION PARMS(2)
      
      ALPHA =     PARMS(1)
      BETA  =     PARMS(2)
      
      IF(BETA .GT. 0.D0) THEN
            FP_G2_CDF   =     FP_G1_CDF(X/BETA,ALPHA)
      ELSE IF(BETA .LT. 0.D0) THEN
            FP_G2_CDF   =     1.D0 - FP_G1_CDF(X/BETA,ALPHA)
      ELSE
        IF(X .GT. 0.D0) THEN
              FP_G2_CDF   =     1.D0
        ELSE
              FP_G2_CDF   =     0.D0
      ENDIF
      ENDIF
      
      RETURN
      END
      
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~

      DOUBLE PRECISION FUNCTION FP_G3_CDF(X,PARMS)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE CDF OF 2-P GAMMA DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DIMENSION PARMS(3)
      
      TAU         =     PARMS(1)
C     ALPHA =     PARMS(2)
C     BETA  =     PARMS(3)
      
            FP_G3_CDF   =     FP_G2_CDF(X-TAU,PARMS(2))

      RETURN
      END
      
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G1_ICDF(P,ALPHA)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE INVERSE OF GAMMA(ALPHA)
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      use iso_fortran_env, only: real64
      implicit none

      REAL(real64) P, ALPHA, DCHIIN
      
      IF(P .GT. 0.999999999) THEN
            FP_G1_ICDF  =     9.D99
            RETURN
      ELSE IF(P .LT. 0.000000001) THEN
            FP_G1_ICDF  =     0.D0
            RETURN      
      ENDIF
      
            FP_G1_ICDF  =     0.5*DCHIIN(P,2.D0*ALPHA)
      RETURN
      END
      
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C
      DOUBLE PRECISION FUNCTION FP_G2_ICDF(P,PARMS)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE INVERSE OF GAMMA(P,ALPHA,BETA)
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      use iso_fortran_env, only: real64
      implicit none

      REAL(real64) P, PARMS(2), ALPHA, BETA, FP_G1_ICDF
      
      ALPHA =     PARMS(1)
      BETA  =     PARMS(2)
      
      IF(BETA .GT. 0.D0) THEN
            FP_G2_ICDF  =     BETA * FP_G1_ICDF(P,ALPHA)
      ELSE IF(BETA .LT. 0.D0) THEN
            FP_G2_ICDF  =     BETA * FP_G1_ICDF(1.D0-P,ALPHA)
      ELSE
            ERROR STOP ' INVALID PARAMETERS IN FP_G2_ICDF:  BETA = 0'
            RETURN
      ENDIF

      RETURN
      END
      
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G3_ICDF(P,PARMS)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE INVERSE OF GAMMA(P,TAU,ALPHA,BETA)
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      use iso_fortran_env, only: real64
      implicit none

      REAL(real64) P, PARMS(3), TAU, FP_G2_ICDF
      
      TAU         =     PARMS(1)
C     ALPHA =     PARMS(2)
C     BETA  =     PARMS(3)
      
      FP_G3_ICDF  =     TAU + FP_G2_ICDF(P,PARMS(2))
      
      RETURN
      END

!*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G1_MOM_TRC(ALPHA,TL,TU,K)
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
!     PROGRAM TO COMPUTE THE EXPECTED VALUE OF THE K-TH MOMENT OF 
!     A 1-PARAMETER GAMMA VARIATE BETWEEN TL AND TU WITH PARAMETER ALPHA
!
!     AUTHOR.....TIM COHN
!     DATE.......DECEMBER 28, 1994
!     
!     MODIFIED...DECEMBER 23, 1996 (TAC)
!     MODIFIED...JULY 25, 1998 (TAC)
!     MODIFIED SEPTEMBER 4, 2024 (SAS) - FOR 128 BIT X IN DGAMDF
!

      use iso_fortran_env, only: real64, real128
      implicit none

      INTEGER J, K
      REAL(real64) ALPHA,TL,TU,TL1,ANS,DOWN,UP,DGAMDF
      REAL(real128) TUQ, TL1Q

            TL1   =     MAX(MIN(0.D0,TU),TL)
            
            TUQ = REAL(TU,real128)
            TL1Q = REAL(TL1, real128)
      IF(TL1 .EQ. TU) THEN
          ANS = TL1**K
      ELSE
        DOWN =  DGAMDF(TUQ ,ALPHA) - DGAMDF(TL1Q,ALPHA)
      IF(DOWN .GT. 0.D0) THEN
        UP    =  DGAMDF(TUQ,ALPHA+K) - DGAMDF(TL1Q,ALPHA+K)
        ANS    =  UP/DOWN
!     THIS SECTION COMPUTES RATIO GAMMA(ALPHA+K)/GAMMA(ALPHA)
            DO 10 J=0,K-1
                  ANS   =     ANS*(ALPHA+J)
10          CONTINUE
      ELSE
        IF(TL1 .GT. ALPHA) THEN
           ANS = TL1**K
        ELSE
           ANS = TU**K
        ENDIF
      ENDIF
      ENDIF

      FP_G1_MOM_TRC     =     ANS

      RETURN
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G2_MOM_TRC(PARMS,TL,TU,K)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C     PROGRAM TO COMPUTE THE EXPECTED VALUE OF THE K-TH MOMENT OF 
C     A 2-PARAMETER GAMMA VARIATE BETWEEN TL AND TU WITH PARAMETERS PARMS

C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 28, 1994
C

      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C      SAVE
      DIMENSION PARMS(2)
      
      ALPHA =     PARMS(1)
      BETA  =     PARMS(2)
      
      IF(BETA .EQ. 0.D0) THEN
            ERROR STOP ' INVALID PARAMETERS IN FP_G2_MOM_TRC: BETA = 0'
            RETURN
      ELSE IF(BETA .GT. 0.D0) THEN
            FP_G2_MOM_TRC     =     BETA**K *  
     1            FP_G1_MOM_TRC(ALPHA,TL/BETA,TU/BETA,K)
      ELSE
            FP_G2_MOM_TRC     =     BETA**K *  
     1            FP_G1_MOM_TRC(ALPHA,TU/BETA,TL/BETA,K)
      ENDIF
            
      RETURN
      END
      
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_G3_MOM_TRC(PARMS,TL,TU,K)
C****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C     PROGRAM TO COMPUTE KTH NON-CENTRAL MOMENT OF PEARSON TYPE III 
C     RANDOM VARIABLE BETWEEN TL AND TU WITH PARAMETERS PARMS
C
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C      SAVE
      DIMENSION PARMS(3)

            TAU         =     PARMS(1)
            ALPHA =     PARMS(2)
            BETA  =     PARMS(3)

            SUM   =     FP_G2_MOM_TRC(PARMS(2),TL-TAU,TU-TAU,K)
            
      IF(TAU .NE. 0.D0) THEN
            DO 10 J=0,K-1
                  SUM   =     SUM + CHOOSE(K,J)*TAU**(K-J)*
     1                       FP_G2_MOM_TRC(PARMS(2),TL-TAU,TU-TAU,J)
10          CONTINUE
      ENDIF
      
            FP_G3_MOM_TRC     =     SUM

      RETURN
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
       SUBROUTINE SP_Z(XSI,F,PHI,LOGPHI,A,G)
C======================================================================
C
C      SUBROUTINE TO COMPUTE GAUSSIAN FUNCTIONS AND THEIR
C      DERIVATIVES, WHERE XSI IS A STANDARD NORMAL VARIATE
C
C           PHI       =  STANDARD NORMAL CDF
C           LOGPHIT   =  LOG(PHIT)
C           AT        =  F/PHIT
C           GT        =  F/(1.D0-PHIT)
C
C           THE ROUTINE WORKS ACCURATELY OVER A VERY LONG RANGE, AND PROVIDES
C           A 2ND ORDER CONTINUITY CORRECTION WHERE ESTIMATION METHODS CHANGE.
C
C
C        AUTHOR......TIM COHN
C        DATE........APRIL 1, 1987
C          ERROR CORRECTED OCTOBER 4, 1992 (TAC)
C                 IMPROVED VERSION....DECEMBER 20, 1994
C
C======================================================================
       IMPLICIT DOUBLE PRECISION (A-Z)
C      SAVE

       DATA LIMIT/6.5D0/,LIMITI/6.D0/,CONST/0.3989422804014327/
!      SAS: CONST is 1/sqrt(2*pi)
         
           Z(B)      =  1.D0/( SQRT(B)*(1.D0-1.D0*B*(1.D0-3.D0*B*
     1                  (1.D0-5.D0*B*(1.D0-7.D0*B*(1.D0-9.D0*B*
     2                  (1.D0-11.D0*B*(1.D0-13.D0*B*(1.D0-15.D0*B*
     3                  (1.D0-17.D0*B*(1.D0-19.D0*B*(1.D0-21.D0*B*
     4                  (1.D0-23.D0*B*(1.D0-12.5D0*B))))))))))))) )
              
         F         =  CONST*EXP(-0.5*XSI**2)
       IF(XSI .LT. -LIMITI) THEN
         PHI       =  F/Z(1.D0/XSI**2) 
         A         =  Z(1.D0/XSI**2)
         LOGPHI    =  LOG(CONST/A) - 0.5*XSI**2
         G         =  F/(1.D0-PHI) 
       ELSE IF (XSI .GT. LIMITI) THEN
         PHI       =  1.D0 - F/Z(1.D0/XSI**2) 
         LOGPHI    =  LOG(PHI) 
         A         =  F/PHI 
         G         =  Z(1.D0/XSI**2)
       ENDIF
         
       IF(ABS(XSI) .LT. LIMIT) THEN
           PHIT       =  DNORDF(XSI)
           LOGPHIT    =  LOG(PHIT)
           AT         =  F/PHIT
           GT         =  F/(1.D0-PHIT)
             IF(ABS(XSI) .GT. LIMITI) THEN
                  IF(ABS(XSI) .GT. LIMIT) THEN
                        ARG   =     1.D0
                  ELSE
                     ARGT      =  (ABS(XSI)-(LIMIT+LIMITI)/2.D0) 
     1                                   /(LIMIT-LIMITI)
                     ARG       =  (1.D0+SIN(3.14159*ARGT))/2.D0
                  ENDIF
             ELSE
               ARG       =  0.D0
             ENDIF
               PHI       =  ARG*PHI + (1.D0-ARG)*PHIT
               LOGPHI    =  ARG*LOGPHI + (1.D0-ARG)*LOGPHIT
               A         =  ARG*A + (1.D0-ARG)*AT
               G         =  ARG*G + (1.D0-ARG)*GT
       ENDIF
       RETURN
       END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION CHOOSE(N,K)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE COMBINATORIAL CHOOSE FUNCTION
C
C                 / N \
C                 \ K /
C
C     N.B.  N AND K MUST BE INTEGERS;  CHOOSE IS DOUBLE PRECISION
C
C     AUTHOR.....TIM COHN
C     DATE.......DECEMBER 21, 1994
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      INTEGER N,K,TAB1(0:12,0:12)
      
      DATA TAB1/
     * 1,  12*0,
     * 1,1,  11*0,
     * 1,2,1,  10*0,
     * 1,3,3,1,  9*0,
     * 1,4,6,4,1,  8*0,
     * 1,5,10,10,5,1,  7*0,
     * 1,6,15,20,15,6,1,  6*0,
     * 1,7,21,35,35,21,7,1,  5*0,
     * 1,8,28,56,70,56,28,8,1,  4*0,
     * 1,9,36,84,126,126,84,36,9,1,  3*0,
     * 1,10,45,120,210,252,210,120,45,10,1,  2*0,
     * 1,11,55,165,330,462,462,330,165,55,11,1,  1*0,
     * 1,12,66,220,495,792,924,792,495,220,66,12,1 /
      
      IF(N .LT. MAX(0,K) .OR. K .LT. 0) THEN
            ERROR STOP 'ERROR IN CHOOSE'
      ENDIF
      
      IF(N .LE. 12) THEN
            CHOOSE      =     TAB1(K,N)
      ELSE IF(N .LT. 170.D0) THEN
            CHOOSE      =     DFAC(N)/(DFAC(N-K)*DFAC(K))
      ELSE
            ARG         =     DLNGAM(N+1.D0)-DLNGAM(N-K+1.D0)- 
     1                                   DLNGAM(K+1.D0)
            IF(ARG .LT. 500.D0) THEN
                  CHOOSE      =     EXP(ARG)
            ELSE
                  CHOOSE      =     -99.D0
            ENDIF
      ENDIF
      
      RETURN
      END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      SUBROUTINE DPDM(PARMS,JACT)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE JACOBIAN DP/DM TRANSFORMATION CONVERTING
C        PEARSON TYPE 3 NON-CENTRAL MOMENTS (M1,M2,M3)
C        INTO PARAMETERS (TAU,ALPHA,BETA) 
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 10, 1999
C
C    N.B.  JAC IS COMPUTED;  JACT = TRANSPOSE(JAC) IS RETURNED
C
C            JAC(1,1) =  DT/DM1
C            JAC(1,2) =  DA/DM1
C            JAC(1,3) =  DB/DM1
C
C            JAC(2,1) =  DT/DM2
C            JAC(2,2) =  DA/DM2
C            JAC(2,3) =  DB/DM2
C
C            JAC(3,1) =  DT/DM3
C            JAC(3,2) =  DA/DM3
C            JAC(3,3) =  DB/DM3
C
C    N.B.: THIS CODE APPEARS TO HANDLE ALL CASES CORRECTLY (TAC 4/10/99)
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      
      DOUBLE PRECISION PARMS(3),M(3),M1,M2,M3,S,A,B,D,T
      DOUBLE PRECISION JAC(3,3),JACT(3,3)
            
      CALL P2MN(PARMS,M)
        M1  =  M(1)
        M2  =  M(2)
        M3  =  M(3)
      
C   V IS THE VARIANCE

      V        =  M2 - M1**2
        DVDM1  =  -2.D0*M1
        DVDM2  =  1.D0
        DVDM3  =  0.D0

C   S  IS THE SKEW       
      S        =  M3 - 3.D0*M2*M1 + 2.D0*M1**3
        DSDM1  =  -3.D0*M2 + 6.D0*M1**2
        DSDM2  =  -3.D0*M1
        DSDM3  =  1.D0

C   A  IS THE PARAMETER

      A        =  4.D0*(V**3)/S**2
        DADM1  =  4.D0*(3.D0*V**2*DVDM1 * S**2 - 2*S*DSDM1*V**3)/S**4
        DADM2  =  4.D0*(3.D0*V**2*DVDM2 * S**2 - 2*S*DSDM2*V**3)/S**4
        DADM3  =  4.D0*(3.D0*V**2*DVDM3 * S**2 - 2*S*DSDM3*V**3)/S**4
        
        
C   D IS B**2

      D        =  V/A
        DDDM1  =  (DVDM1*A-DADM1*V)/A**2
        DDDM2  =  (DVDM2*A-DADM2*V)/A**2
        DDDM3  =  (DVDM3*A-DADM3*V)/A**2
        
C   B IS SIGN(S)*SQRT(D)

      B        =  SIGN(1.D0,S)*SQRT(D)
        DBDM1  =  0.5*DDDM1/B
        DBDM2  =  0.5*DDDM2/B
        DBDM3  =  0.5*DDDM3/B

C   T

      T        =  M1 - A*B
        DTDM1  =  1.D0 - (DADM1*B + A*DBDM1)
        DTDM2  =       - (DADM2*B + A*DBDM2)
        DTDM3  =       - (DADM3*B + A*DBDM3)

      JAC(1,1) =  DTDM1
      JAC(1,2) =  DADM1
      JAC(1,3) =  DBDM1

      JAC(2,1) =  DTDM2
      JAC(2,2) =  DADM2
      JAC(2,3) =  DBDM2

      JAC(3,1) =  DTDM3
      JAC(3,2) =  DADM3
      JAC(3,3) =  DBDM3
      
      DO I=1,3
        DO J=1,3
          JACT(I,J) = JAC(J,I)
        END DO
      END DO

      RETURN
      END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C        
        SUBROUTINE P2M(P,M)
C        SAVE
        IMPLICIT NONE
        DOUBLE PRECISION P(3),M(3),T,A,B
        T=P(1)
        A=P(2)
        B=P(3)
        M(1) =  A*B + T
        M(2) =  A*B**2
        M(3) =  2.D0/SQRT(A)
        IF(B .LT. 0.D0) M(3) = -M(3)
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
        SUBROUTINE M2P(M,P)
        IMPLICIT NONE
C        SAVE
        DOUBLE PRECISION P(3),A,B,T,M(3),M1,M2,M3
        M1=M(1)
        M2=M(2)
        M3=M(3)
        A = 4.D0/M3**2
        B = SQRT(M2/A)
        IF(M3 .LT. 0.D0) B=-B
        T = M1 - A*B
        P(1)  =  T
        P(2)  =  A
        P(3)  =  B
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C        
        SUBROUTINE P2MN(P,M)
C        SAVE
        IMPLICIT NONE
        DOUBLE PRECISION P(3),M(3),T,A,B
        T=P(1)
        A=P(2)
        B=P(3)
        M(1) =  A*B + T
        M(2) =  A*(1 + A)*B**2 + 2*A*B*T + T**2
        M(3) =  A*(1 + A)*(2 + A)*B**3 + 3*A*(1 + A)*B**2*T + 
     1             3*A*B*T**2 + T**3
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
        SUBROUTINE MN2P(M,P)
        IMPLICIT NONE
C      SAVE
        DOUBLE PRECISION P(3),M(3),M1,M2,M3
        M1=M(1)
        M2=M(2)
        M3=M(3)
        P(1) = (M1**2*M2 - 2*M2**2 + M1*M3)/(2*M1**3 - 3*M1*M2 + M3)
        P(2) = -4*(M1**2 - M2)**3/(2*M1**3 - 3*M1*M2 + M3)**2
        P(3) = (2*M1**3 - 3*M1*M2 + M3)/(-2*M1**2 + 2*M2)
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C        
        SUBROUTINE M2MN(M,MN)
C      SAVE
        IMPLICIT NONE
        DOUBLE PRECISION M(3),MN(3)
        MN(1) =  M(1)
        MN(2) =  M(2)+M(1)**2
        MN(3) =  M(3)*M(2)**1.5 + 3.D0*MN(2)*M(1) - 2.D0 * M(1)**3
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C        
        SUBROUTINE MN2M(MN,M)
C      SAVE
        IMPLICIT NONE
        DOUBLE PRECISION M(3),MN(3)
        M(1) =  MN(1)
        M(2) =  MN(2)-MN(1)**2
        M(3) =  (MN(3) - 3.D0*MN(2)*MN(1) + 2.D0 * MN(1)**3)
     1              /M(2)**1.5
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
         SUBROUTINE EXPMOMDERIV(PARMS,XMIN,XMAX,JAC)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE DERIVATIVE OF EXPECTED VALUE OF 
C         MOMENTS OF GAMMA VARIATE TRUNCATED BETWEEN XMIN AND XMAX
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 10, 1999
C
C     REQUIRES FUNCTIONS FOR LOG-GAMMA AND PSI
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DOUBLE PRECISION 
     1  PARMS(3),XMIN,XMAX,JAC(3,3),
     2  JAC1(3,3),JAC2(3,3),ALPHA,BETA,TAU,M(3)

          TAU   = PARMS(1)
          ALPHA = PARMS(2)
          BETA  = PARMS(3)
        CALL DEXPECT(TAU,ALPHA,BETA,XMIN,XMAX,M,JAC1)
         
        CALL DPDM(PARMS,JAC2)
        
        CALL DMRRRR(3,3,JAC1,3,3,3,JAC2,3,3,3,JAC,3)
        
        RETURN
        END
!*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
         SUBROUTINE EXPMOMCDERIV(PARMS, TL, TU, MNE, DEDMC)
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
! 
!     PROGRAM TO COMPUTE DERIVATIVE OF EXPECTED VALUE OF CENTRAL MOMENTS
!         OF GAMMA VARIATE TRUNCATED BETWEEN XMIN AND XMAX
!
!     FOR IMPLEMENTATION OF GREG SCHWARZ'S HALLOWEEN SKEW WEIGHTING METHOD
!
!     AUTHOR.....SETH SIEFKEN
!     DATE.......NOVEMBER 4, 2022
!     MODIFIED...DECEMBER 22, 2023
!
!     INPUT VARIABLES:
!          PARMS CONTAINS P3 PARAMETER VALUES (TAU, ALPHA, BETA)
!          TL = LOWER CENSORING THRESHOLD
!          TU = UPPER CENSORING THRESHOLD
!
!     RETURN VARIABLES:
!          MNE[J]      CONTAINS E[X^J | X IS CENSORED]
!          DEDMC[J,K]   CONTAINS D[ E[MC[J] | X IS CENSORED] ]/MC[K]
!
!          DEDMC[1,1]  =   D[ E[MU | X IS CENSORED] ]/D_MU
!          DEDMC[1,2]  =   D[ E[MU | X IS CENSORED] ]/D_SIGMA2
!          DEDMC[1,3]  =   D[ E[MU | X IS CENSORED] ]/D_SKEW
!
!          DEDMC[2,1]  =   D[ E[SIGMA2 | X IS CENSORED] ]/D_MU
!          DEDMC[2,2]  =   D[ E[SIGMA2 | X IS CENSORED] ]/D_SIGMA2
!          DEDMC[2,3]  =   D[ E[SIGMA2 | X IS CENSORED] ]/D_SKEW
!
!          DEDMC[3,1]  =   D[ E[SKEW | X IS CENSORED] ]/D_MU
!          DEDMC[3,2]  =   D[ E[SKEW | X IS CENSORED] ]/D_SIGMA2
!          DEDMC[3,3]  =   D[ E[SKEW | X IS CENSORED] ]/D_SKEW
!
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION   PARMS(3), TL, TU, DEDMC(3,3), MC(3)
      DOUBLE PRECISION   db(3,3), B(3,3), tp(3,3)
      DOUBLE PRECISION   ALPHA, BETA, TAU, JACL(3,3), JACG(3,3), D2(3,3)
      DOUBLE PRECISION   D3(3,3), BM(3,3), D12(3,3), PIL, PIG
      DOUBLE PRECISION   MNEL(3), MNEG(3), MNE(3), INF
      DOUBLE PRECISION   PIJACL(3,3), PIJACG(3,3), JAC(3,3), PICK
      
      data inf/20d0/ !Make sure infinity is consistent throughout program
      
      
        CALL P2M(PARMS,MC) !Compute central moments, mc(2) = variance
        
        
        S = SQRT(MC(2))
        
        TAU   = PARMS(1)
        ALPHA = PARMS(2)
        BETA  = PARMS(3)
          


        !Compute probability of censoring
        PIL = FP_G3_CDF(TL,PARMS) 
        PIG = 1 - FP_G3_CDF(TU,PARMS) 
        PICK = PIL + PIG

        !Compute individual Jacobians for upper and lower censoring
        CALL DMMULT(PIL,3,3,JACL,3,3,3,PIJACL,3) 
        CALL DMMULT(PIG,3,3,JACG,3,3,3,PIJACG,3) 
        
        
        IF (PICK .GT. 1d-10) THEN
        
        !Compute individual Jacobians for upper and lower censoring
            CALL DEXPECT(TAU,ALPHA,BETA,-inf,tl,MNEL,JACL)
            CALL DEXPECT(TAU,ALPHA,BETA,tu, inf,MNEG,JACG)
            
            CALL DMMULT(PIL,3,3,JACL,3,3,3,PIJACL,3) 
            CALL DMMULT(PIG,3,3,JACG,3,3,3,PIJACG,3) 
            
        !Compute expected values of noncentral moments when X is censored
            MNE(1) = (PIL*MNEL(1) + PIG*MNEG(1))/PICK
            MNE(2) = (PIL*MNEL(2) + PIG*MNEG(2))/PICK
            MNE(3) = (PIL*MNEL(3) + PIG*MNEG(3))/PICK
            
        !Compute Jacobian when X is censored
            CALL DMSUM(3,3,PIJACL,3,3,3,PIJACG,3,3,3,JAC,3)
            CALL DMMULT(1/PICK,3,3,JAC,3,3,3,JAC,3) 
        ELSE
        !When X is not (effectively) censored can just use DEXPECT
            CALL DEXPECT(TAU,ALPHA,BETA,-inf,inf,MNE,JAC)
        
        ENDIF
        
        
        
        
        
        
!       Compute derivative of b vector from Schwarz (2022)
        db(1,1) = 0
        db(1,2) = 0
        db(1,3) = 0

        db(2,1) = 2*MC(1)
        db(2,2) = 0
        db(2,3) = 0

        db(3,1) = -3*mc(1)**2/S**3
        db(3,2) = 3*mc(1)**3/(2*S**5)
        db(3,3) = 0
        


!       Compute B matrix from Schwarz (2022)
        B(1,1) = 1d0
        B(1,2) = 0d0
        B(1,3) = 0d0

        B(2,1) = -2*mc(1)
        B(2,2) = 1d0
        B(2,3) = 0d0

        B(3,1) = 3*mc(1)**2/S**3
        B(3,2) = -3*mc(1)/S**3
        B(3,3) = 1/S**3

!       Compute theta prime matrix from Schwarz (2022)
        tp(1,1) = 1d0
        tp(1,2) = -1d0/(S*mc(3))
        tp(1,3) = 2*S/mc(3)**2
        
        tp(2,1) = 0d0
        tp(2,2) = 0d0
        tp(2,3) = -8/mc(3)**3

        tp(3,1) = 0d0
        tp(3,2) = mc(3)/(4*S)
        tp(3,3) = S/2

!       Compute third term of derivative
        D3(1,1) = 0d0
        D3(1,2) = 0d0
        D3(1,3) = 0d0

        D3(2,1) = -2d0*MNE(1)
        D3(2,2) = 0d0
        D3(2,3) = 0d0

      D3(3,1)=6*mc(1)/S**3*MNE(1) - 3/S**3*MNE(2)
      D3(3,2)=(-4.5*mc(1)**2*MNE(1)+4.5*mc(1)*MNE(2)-1.5*MNE(3))/S**5
      D3(3,3)=0d0

        

        
        CALL DMRRRR(3,3,B,3,3,3,JAC,3,3,3,BM,3) 
        CALL DMRRRR(3,3,BM,3,3,3,tp,3,3,3,D2,3) !2nd term in EQ 27
        
        CALL DMSUM(3,3,db,3,3,3,D2,3,3,3,D12,3)
        CALL DMSUM(3,3,D12,3,3,3,D3,3,3,3,DEDMC,3)
        
        
        RETURN
        END


        
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION DDGAM(ALPHA,X)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE DERIVATIVE OF INCOMPLETE GAMMA FUNCTION
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 9, 1999
C
C     REQUIRES FUNCTIONS FOR LOG-GAMMA AND PSI
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DOUBLE PRECISION I,TOL,LOGX,SUM
      
       DATA TOL/1.D-11/
      
       IF(ALPHA .LE. 0.D0 .OR. X .LE. 0.D0) THEN
               RESULT = 0.D0
               GOTO 100
       END IF
              
       IF(ABS(X-ALPHA)/SQRT(ALPHA+1.D0) .GT. 7.D0) THEN
               RESULT = 0.D0
               GOTO 100
       END IF
      
          A        = ALPHA
          LOGX     =  LOG(X)
C
      IF(X .LT. 5.D0) THEN
          T        =  X**A
          R        =  (A*LOGX-1.D0)/A**2
          SUM      =  T*R
        DO 10 I2=1,1000
          I        =  I2
          AI       =  A+I
          T        =  -T*X/I
          R        =  (AI*LOGX-1.D0)/AI**2
          DEL      =  R*T
          SUM      =  SUM + DEL
           IF(I .GT. 1 .AND. ABS(DEL) .LT. (1.D0+ABS(SUM))*TOL) THEN
             RESULT =  SUM/EXP(DGAMLN(A,IERR)) - DPSI(A)*FP_G1_CDF(X,A)
             GOTO 100
           ENDIF
10      CONTINUE
          GOTO 99
C
      ELSE IF(X .GT. A+30.D0) THEN
          T        =  EXP(-X + (A-1.D0)*LOGX - DGAMLN(A,IERR) )
          R        =  LOGX - DPSI(A)
          SUM      =  R*T
        DO 20 I2=1,INT(A-1.D0)
          I         =  I2
          AMI       =  A-I
          T        =  T*(AMI)/X
          R        =  LOGX - DPSI(AMI)
          DEL      =  R*T
          SUM      =  SUM + DEL
           IF(I .GT. 1 .AND. ABS(DEL) .LT. (1.D0+ABS(SUM))*TOL) THEN
             RESULT = -SUM
             GOTO 100
           ENDIF
20      CONTINUE
          GOTO 99
C
      ELSE
          T        =  EXP(-X + A*LOGX - DGAMLN(A+1.D0,IERR) )
          R        =  LOGX - DPSI(A+1.D0)
          SUM      =  R*T
        DO 30 I2=1,10000
          I        =  I2
          T        =  T*X/(A+I)
          R        =  LOGX - DPSI(A+I+1.D0)
          DEL      =  R*T
          SUM      =  SUM + DEL
           IF(I .GT. 1 .AND. ABS(DEL) .LT. (1.D0+ABS(SUM))*TOL) THEN
             RESULT = SUM
             GOTO 100
           ENDIF
30      CONTINUE
          GOTO 99
      ENDIF
C
100   CONTINUE
        DDGAM      =  RESULT
      RETURN
99    CONTINUE
!             WRITE(*,*) 'SOMETHING BAD (DDGAM).....'
!             WRITE(*,*) 'Contact Tim Cohn (703/648-5711).....'
!             WRITE(*,*) 'tacohn@usgs.gov.....'
!             WRITE(*,*) 'ALPHA: ',ALPHA
!             WRITE(*,*) 'X:     ',X
              ERROR STOP 'ERROR IN DDGAM'
      END

C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      SUBROUTINE DEXPECT(TAU,ALPHA,BETA,TL,TU,M3,DM3)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE DERIVATIVE OF INCOMPLETE GAMMA FUNCTION
C       TRUNCATED BELOW AT TL AND ABOVE AT TU, WITH PARMS TAU, ALPHA, BETA
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 9, 1999
C
C     REQUIRES FUNCTIONS FOR LOG-GAMMA AND PSI
C
C     NOTE RETURN VARIABLES:
C         M3[J]      CONTAINS E[X^J]
C         DM3[J,K]   CONTAINS D[ E[X^J] ]/P[K]
C
C          DM3[1,1]  =   D[ E[X] ]/D_TAU
C          DM3[1,2]  =   D[ E[X] ]/D_ALPHA
C          DM3[1,3]  =   D[ E[X] ]/D_BETA
C
C          DM3[2,1]  =   D[ E[X^2] ]/D_TAU
C          DM3[2,2]  =   D[ E[X^2] ]/D_ALPHA
C          DM3[2,3]  =   D[ E[X^2] ]/D_BETA
C
C          DM3[3,1]  =   D[ E[X^3] ]/D_TAU
C          DM3[3,2]  =   D[ E[X^3] ]/D_ALPHA
C          DM3[3,3]  =   D[ E[X^3] ]/D_BETA
C
C       WHERE P[1] = TAU, P[2] = ALPHA, P[3] = BETA
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE
      DOUBLE PRECISION 
     1    ALPHA,BETA,TAU,TL,TU,
     2    GU,GL,DG1(0:3),G1(0:3),
     3    ADJ(0:3),DADJ(0:3),B(0:3),DB(0:3),
     4    M1(3),DM1(0:3,3),M2(3),DM2(3,3),M3(3),DM3(3,3),
     5    LSTDL,LSTDU,FU(0:3),FL(0:3),STDL,STDU

C
C    RETURN IF THE RESULTS DO NOT DEPEND ON THE PARAMETERS
C
        IF(TU .LE. TL) GOTO 99
        
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C    STANDARDIZE THE THRESHOLDS;  COMPUTE DERIVATIVES
C
        IF(BETA .GT. 0.D0) THEN
          STDL  =  MAX(0.D0,(TL-TAU)/BETA)
          STDU  =  (TU-TAU)/BETA
        ELSE IF(BETA .LT. 0.D0) THEN
          STDL  =  MAX(0.D0,(TU-TAU)/BETA)
          STDU  =  (TL-TAU)/BETA
        ELSE IF(BETA .EQ. 0.D0) THEN
          GOTO 99
        ENDIF
          LSTDL   =  LOG(MAX(1.D-99,STDL))
          LSTDU   =  LOG(MAX(1.D-99,STDU))

        IF(LSTDL .EQ. LSTDU) GOTO 99 !Added by SAS, for Halloween weighting

          DSTDLT  =  -1.D0/BETA
          DSTDUT  =  -1.D0/BETA
          DSTDLB  =  -STDL/BETA
          DSTDUB  =  -STDU/BETA

C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C    1.  FIRST FOR THE ONE-PARAMETER GAMMA (ASSUME TAU=0; BETA=1)
C         -- NOTE THE FOLLOWING:
C              E[X^P] = GAMMA[ALPHA+P,TL,TU] / GAMMA[ALPHA,TL,TU]
C                     = ADJ[P]*B[P]
C              ADJ[P] = GAMMA[ALPHA+P]/GAMMA[ALPHA]
C              B[P]   = GCDF[ALPHA+P,TL,TU] / GCDF[ALPHA,TL,TU]
C              GCDF[ALPHA,TL,TU] = GAMMA[ALPHA,TL,TU]/GAMMA[ALPHA,0,INF]
C

C
C      1.A   FIRST COMPUTE THE ADJUSTMENT FACTOR
C
          ADJ(0)    =  1.D0
          DADJ(0)   =  0.D0
        DO 40 J=1,3
          ADJ(J)    =  (ALPHA+J-1.D0) * ADJ(J-1)
          DADJ(J)   =  ADJ(J-1) + (ALPHA+J-1.D0) * DADJ(J-1)
40      CONTINUE

C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C      1.B   SECOND COMPUTE THE INCOMPLETE G1 CDF AND ITS DERIVATIVES
C              -- NB:  B(J) IS THE NORMALIZED CDF
C                      DB(J) IS COMPUTED FROM RATIO
C                      G1 CONTAINS TAU, BETA THROUGH STDL
C
      DO 20 J=0,3
        GU      =  FP_G1_CDF(STDU,ALPHA+J)
        GL      =  FP_G1_CDF(STDL,ALPHA+J)
        G1(J)   =  GU-GL
        DG1(J)  =  DDGAM(ALPHA+J,STDU) - DDGAM(ALPHA+J,STDL)
        B(J)    =  G1(J)/G1(0)
        DB(J)   =  (DG1(J)*G1(0) - DG1(0)*G1(J))/G1(0)**2
        FU(J)   =  FP_G1_PDF(STDU,ALPHA+J)
        FL(J)   =  FP_G1_PDF(STDL,ALPHA+J)
        DM1(J,1) = (
     1     (FU(J)*DSTDUT - FL(J)*DSTDLT) * G1(0) -
     2     (FU(0)*DSTDUT - FL(0)*DSTDLT) * G1(J) )/G1(0)**2
        DM1(J,3) = (
     1     (FU(J)*DSTDUB - FL(J)*DSTDLB) * G1(0) -
     2     (FU(0)*DSTDUB - FL(0)*DSTDLB) * G1(J) )/G1(0)**2
20    CONTINUE
    
C
C      1.C   NOW COMPUTE THE DERIVATIVE OF THE PRODUCT
C
      DO 50 J=1,3
        M1(J)     =  ADJ(J)*B(J)
        DM1(J,2)  =  ADJ(J)*DB(J) + DADJ(J)*B(J)
        DM1(J,1)  =  ADJ(J)*DM1(J,1)
        DM1(J,3)  =  ADJ(J)*DM1(J,3)
50    CONTINUE

C
C
C    2.  NOW GENERALIZE TO THE 2-PARAMETER GAMMA (ASSUME TAU=0)
C        --NOTE THAT ORDER OF PARAMETERS IS TAU,ALPHA,BETA
C          (TAC -- ANOTHER GOOFY CONVENTION); THEREFORE
C              DM2(J,1) REFERS TO TAU
C              DM2(J,2) REFERS TO ALPHA
C              DM2(J,3) REFERS TO BETA
C
      DO 60 J=1,3
        M2(J)     = BETA**J * M1(J)
        DM2(J,1)  = BETA**J * DM1(J,1)
        DM2(J,2)  = BETA**J * DM1(J,2)
        DM2(J,3)  = J * BETA**(J-1) * M1(J) + BETA**J * DM1(J,3)
60    CONTINUE
C
C
C    3.  NOW GENERALIZE TO THE 3-PARAMETER GAMMA
C
C
      DO J=1,3
          M3(J)     = M2(J) + TAU**J
          DM3(J,1)  = J * POWER(TAU,(J-1)) + DM2(J,1)
          DM3(J,2)  = DM2(J,2)
          DM3(J,3)  = DM2(J,3)
        DO K=1,J-1
          CH        =  CHOOSE(J,K)
          M3(J)     =  M3(J) + CH*TAU**(J-K) * M2(K)
          DM3(J,1)  =  DM3(J,1) + CH*(
     1                   (J-K)*POWER(TAU,(J-K-1)) * M2(K) +
     2                   TAU**(J-K) * DM2(K,1) )
          DM3(J,2)  =  DM3(J,2) + CH*TAU**(J-K) * DM2(K,2)
          DM3(J,3)  =  DM3(J,3) + CH*TAU**(J-K) * DM2(K,3)
        END DO
      END DO

      RETURN
C 
C  NON-FUNCTIONAL EXIT
C
99    CONTINUE
          CALL DSET(9,0.D0,DM3,1)
          CALL DSET(3,0.D0,M3,1)
          RETURN
      END
C
      DOUBLE PRECISION FUNCTION POWER(X,I)
C      SAVE
      DOUBLE PRECISION X
      INTEGER I
      IF(X .EQ. 0.D0 .AND. I .EQ. 0) THEN
         POWER = 1.D0
      ELSE
         POWER = X**I
      ENDIF
      RETURN
      END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_TNC_CDF(TP,NU,DELTA)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE NON-CENTRAL T CDF
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 9, 1999
C
C     REQUIRES FUNCTIONS FOR Z
C
C     REFERENCE:  APPROXIMATE FORMULA ON PAGE 949, A&S
C                 TNC FUNCTION IS ALGORITHM AS 243  
C                    APPL. STATIST. (1989), VOL.38, NO. 1
C
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT NONE
C      SAVE
      
      INTEGER IFAULT
      
      DOUBLE PRECISION TP,NU,DELTA,Z,FP_Z_CDF,TNC,ANS
      
      IF( NU .GT. 20.D0) THEN
       Z  = (TP*(1.D0-1.D0/(4.D0*NU))-DELTA)/SQRT(1.D0+TP**2/(2.D0*NU))
       ANS=  FP_Z_CDF(Z) 
      ELSE
       ANS=  TNC(TP,NU,DELTA,IFAULT)
       IF(IFAULT .NE. 0) THEN
            ERROR STOP 'ERROR IN FP_TNC_CDF'
       Z  = (TP*(1.D0-1.D0/(4.D0*NU))-DELTA)/SQRT(1.D0+TP**2/(2.D0*NU))
       ANS=  FP_Z_CDF(Z) 
       ENDIF
      ENDIF
      
      FP_TNC_CDF =  ANS
      
      RETURN
      END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FP_TNC_ICDF(P_IN,NU_IN,DELTA_IN)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE NON-CENTRAL T INVERSE CDF
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 9, 1999
C
C     REQUIRES FUNCTIONS FOR Z-INVERSE, GAMMA INVERSE
C
C     NOTE:  THIS USED TO BE LABELED FP_TNC2_ICDF (TAC: 20 SEP 2010)
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      IMPLICIT NONE
C      SAVE
      
      INTEGER IFLAG
      
      DOUBLE PRECISION 
     1   P_IN,NU_IN,DELTA_IN,
     2   P,NU,DELTA,
     3   AE,RE,B,C
     
      DOUBLE PRECISION FTNCINV
      
      EXTERNAL FTNCINV
      
      COMMON /ZNCT02/P,NU,DELTA

      DATA RE/1.D-8/,AE/1.D-8/

      P     =  P_IN
      NU    =  NU_IN
      DELTA =  DELTA_IN
      
c  change made to limits 20 Sep 2010 (TAC)
c      X1    =  FP_Z_ICDF(P)+DELTA
c                PARMS(1)  =  NU/2.D0
c                PARMS(2)  =  2.D0
c      X2    =  SQRT(FP_G2_ICDF(1.D0-P,PARMS)/NU)
c      
c      IF(X2 .LE. 0.D0) GOTO 99
c      
c        B     =  MIN(0.D0,X1,1.D0/X2,X1/X2)
c        C     =  MAX(0.D0,X1,1.D0/X2,X1/X2)

         B     =  -32.D0  !  Works for 0.01<=p<=0.99,nu=1,nc=0  07/15/2013 (TAC)
         C     =   32.D0 
                      
      CALL DZEROIN(FTNCINV,B,C,RE,AE,IFLAG)
      
      IF(IFLAG .NE. 1) GOTO 99
      
      FP_TNC_ICDF =  B
      
      RETURN
      
99    CONTINUE
        ERROR STOP 'ERROR IN FP_TNC_ICDF'
        FP_TNC_ICDF =  0.D0
      
      RETURN
      END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
      DOUBLE PRECISION FUNCTION FTNCINV(TP)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE INVERSE CDF OF NON-CENTRAL T DISTRIBUTION
C
C     AUTHOR.....TIM COHN
C     DATE.......APRIL 19, 1999
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c        SAVE
        DOUBLE PRECISION TP,ANS,P,NU,DELTA,FP_TNC_CDF
      
        COMMON /ZNCT02/P,NU,DELTA

          ANS     =    FP_TNC_CDF(TP,NU,DELTA) - P
          FTNCINV =    ANS
          
        RETURN
        END
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
C
      DOUBLE PRECISION FUNCTION FP_B_ICDF(P,A,B)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C 
C     PROGRAM TO COMPUTE INVERSE CDF OF BETA(P,A,B)
C
C     AUTHOR.....TIM COHN
C     DATE.......08 MAR 2010
C
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

      DOUBLE PRECISION
     1  P,X,Y,A,B,BOUND
      
      INTEGER
     1  WHICH,STATUS
       
      DATA WHICH/2/

      CALL CDFBET ( WHICH, P, 1.D0-P, X, Y, A, B, STATUS, BOUND )
        IF(STATUS .NE. 0) ERROR STOP ' ERROR IN CDFBET/FP_B_ICDF'
      
      FP_B_ICDF = X
      
      RETURN
      END
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
      INTEGER FUNCTION MGBTP(X,N,PVALUEW)
!*** REVISION 2012.09
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
! 
!       IDENTIFIES THE NUMBER OF LOW OUTLIERS USING THE 
!       GENERALIZED GRUBBS-BECK TEST
!       NOTES:  SEE 2011 MANUSCRIPT BY COHN, STEDINGER, ENGLAND, ET AL.
! 
!       AUTHOR.........TIM COHN
!       DATE...........1 MAR 2010 (TAC)
!               MODIFIED....... 8 MAR 2010 (TAC)
!               MODIFIED.......29 APR 2010 (TAC) [OFFSET ADDED FOR SMALL N]
!               MODIFIED.......21 OCT 2010 (TAC) [OFFSET REMOVED; ESTIMATOR
!                                  FOR E[S] REPLACED WITH GAMMA FUNCTION
!               MODIFIED.......5 JAN 2011 (TAC)
!               MODIFIED.......16 FEB 2011 (TAC) RETURNS PVALUEW
!                 MODIFICATIONS FOR MGBT TESTING
!                  John England, Bureau of Reclamation, 17-SEP-2012                     
!               MODIFIED.......26 SEP 2012 (TAC) REWORKED CODE TO CLARIFY LOGIC
!                                                AND AVOID STDOUT CALLS
!
!c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
!       calling routine: gbtest
!
!       input variables (from subroutine gbtest):
!       ------------------------------------------------------------------------
!            n          i*4  number of observations (uncensored, with ql=qu)
!            x(n)       r*8  vector of floods in real space (qs from gbtest)
!
!                       note: see subroutine gbtest for input and call
!                       input is systematic 'Syst' records where ql_in=qu_in
!                       and tl_in <= gage_base from peakfqSA
!
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
!       output variables:
!       ------------------------------------------------------------------------
!           pvaluew(n)  r*8  vector of MGB p-values
!                           (function of n,k,w) returned from ggbcritp
!
!****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
!       controlling parameters:
!       ------------------------------------------------------------------------
!           Alphaout    r*8  significance level for outward sweep from median
!                            valid value range: >0
!
!           Alphain     r*8  significance level for inward sweep from smallest obs
!                            valid value range: any; controlled by tiny
!
!           Alphazeroin r*8  significance level for inward sweep from smallest obs
!                            valid value range: >0
!
!       default parameter values stored in block data
!         data Alphaout/0.005d0/,Alphain/0.d0/,Alphazeroin/0.1d0/      
!
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
      IMPLICIT NONE
      
      INTEGER DIM1,N
C JFE      PARAMETER (DIM1=10000)
      PARAMETER (DIM1=20000)
      
      DOUBLE PRECISION 
     * X(N),PVALUEW(DIM1),W(DIM1),ZT(DIM1),S1,S2,XM,XV,
     1 GGBCRITP,Alphaout,Alphain,Alphazeroin

      INTEGER 
     * I,J1,J2,J3,N2,NC

      common /MGB001/Alphaout,Alphain,Alphazeroin

      DO 10 I=1,N
         ZT(I) = LOG10(MAX(1D-88,X(I)))
         PVALUEW(I) = -99.D0    ! MISSING VALUES
 10   CONTINUE

!     sort log flows from smallest to largest
      CALL DSVRGN(N,ZT,ZT)
        
!     set starting point for MGBT search at approximate median position (1/2 N)
      N2    = N/2

      S1 = 0.D0
      S2 = 0.D0

      DO 20 I=N,N2+2,-1
         S1 = S1 + ZT(I)
         S2 = S2 + ZT(I)**2
 20   CONTINUE

      DO 30 I=N2,1,-1
         S1  = S1 + ZT(I+1)
         S2  = S2 + ZT(I+1)**2
         NC  = N-I
         XM  = S1/DBLE(NC)
         XV  = (S2-NC*XM**2)/(NC-1.D0)
         W(I)       = (ZT(I)-XM)/SQRT(XV)
         PVALUEW(I) = GGBCRITP(N,I,W(I))
 30   CONTINUE
!
!     Determine number of low outliers in 2 or 3 steps.
!     Based on TAC original code and JRS recommendations.
!
!     Step 1. Outward sweep from median (always done).
!             alpha level of test = Alphaout
!             number of outliers = J1
!
!     Step 2. Inward sweep from largest low outlier identified in Step 1.
!             alpha level of test = Alphain
!             number of outliers = J2
!
!     Step 3. Inward sweep from smallest observation
!             alpha level of test = Alphazeroin
!             number of outliers = J3
!
!     Initialize counters
      J1 = 0  ! Outward sweep number of low outliers
      J2 = 0  ! Inward sweep number of low outliers
      J3 = 0  ! Oth Inward sweep number of low outliers

!     1) Outward sweep check: Loop over low flows up to median
      DO 110 I=N2,1,-1
         IF(PVALUEW(I) .LT. Alphaout) GOTO 111
110   CONTINUE
111   CONTINUE
        J1=I

!     2) Inward sweep check with Alphain
           J2 = J1
      DO 120 I=J1+1,N2  
        IF( (PVALUEW(I) .GE. Alphain)) GOTO 121
120   CONTINUE
121   CONTINUE
        J2 = I-1

!     3) Inward sweep check with Alphazeroin
      DO 130 I=1,N2
        IF( (PVALUEW(I) .GE. Alphazeroin)) GOTO 131
130   CONTINUE
131   CONTINUE
        J3=I-1
        
!     Set  number of low outliers as max of 3 sweeps
      MGBTP = MAX(J1,J2,J3)
      RETURN
      END
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
      blockdata MGBblk
      implicit none
      double precision Alphaout,Alphain,Alphazeroin     
      common /MGB001/Alphaout,Alphain,Alphazeroin      
      data Alphaout/0.005d0/,Alphain/0.d0/,Alphazeroin/0.1d0/      
      end

        DOUBLE PRECISION FUNCTION GGBCRITP(N,R,ETA)
!*** REVISION 2011.01
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
! 
!       PROGRAM TO COMPUTE P-VALUES (GGCRITP) AND CRITICAL POINTS  (GGCRITK)
!         FOR A GENERALIZED GRUBBS-BECK TEST
!       NOTES:  SEE 2010 MANUSCRIPT BY COHN, ENGLAND, ET AL.
! 
!       AUTHOR.........TIM COHN
!       DATE...........1 MAR 2010 (TAC)
!               MODIFIED....... 8 MAR 2010 (TAC)
!               MODIFIED.......29 APR 2010 (TAC) [OFFSET ADDED FOR SMALL N]
!               MODIFIED.......21 OCT 2010 (TAC) [OFFSET REMOVED; ESTIMATOR
!                                  FOR E[S] REPLACED WITH GAMMA FUNCTION
!               MODIFIED........5 JAN 2011 (TAC) MGBT REVISED
!
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
      IMPLICIT NONE
      
      INTEGER
     1 LIMIT
     
      PARAMETER
     1 (LIMIT=40000)
     
      INTEGER
     * N,R,N_IN,R_IN,
     1 NEVAL,IER,LENW,LAST,IWORK(LIMIT),
     2 KEY
     
      DOUBLE PRECISION 
     * ETA,ETA_IN,
     1 FGGB,EPSABS,EPSREL,RESULT,ABSERR,WORK(4*LIMIT)
     

      EXTERNAL FGGB
      
      COMMON /GGB001/ETA_IN,N_IN,R_IN
 
      DATA EPSABS/1.D-5/,EPSREL/1.D-5/,KEY/3/,LENW/160000/
      
     
      IF(N .LT. 2 .OR. R .GT. N/2) THEN
        ERROR STOP ' INVALID INPUT FOR GGBCRITP'
!       GGBCRITP = 0.5D0
!       RETURN
      ELSE
          N_IN   = N
          R_IN   = R
          ETA_IN = ETA
      ENDIF
      
      CALL DQAG(FGGB,0.D0,1.D0,EPSABS,EPSREL,KEY,RESULT,ABSERR,
     *    NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

      GGBCRITP   =  RESULT
      
      RETURN
      END
            
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
      DOUBLE PRECISION FUNCTION FGGB(PZR)
!*** REVISION 2011.01
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
!****|==|====-====|====-====|====-====|====-====|====-====|====-====|==////////
!
      IMPLICIT NONE
      
      INTEGER
     * N,R
     
      DOUBLE PRECISION 
     * PZR,ETA,DF,MuM, MuS2,VarM,VarS2,CovMS2,EX1,EX2,EX3,EX4,
     2 CovMS,VarS,alpha,beta,dlngam,
     3 MuMP,EtaP,H,Lambda,MuS,ncp,q,VarMP,PR,ZR,N2,
     4 FP_B_ICDF,FP_Z_ICDF,FP_Z_PDF,
     5 FP_TNC_CDF

      COMMON /GGB001/ETA,N,R
!
!  Compute the value of the r-th smallest obs. based on its order statistic
!
      N2   = dble(N-R)
      PR   = FP_B_ICDF(PZR,dble(R),dble(N+1-R))
      ZR   = FP_Z_ICDF(PR)

!
!  Calculate the expected values of M, S2, S and their variances/covariances
!
      H        =  FP_Z_PDF(ZR)/(Max(1.d-10,1.d0-PR))
      EX1      =  H
      EX2      =  1 + H * ZR
      EX3      =  2 * EX1 + H * ZR**2
      EX4      =  3 * EX2 + H * ZR**3
      
      MuM      =  EX1
      MuS2     =  (EX2-EX1**2)
      VarM     =  MuS2/N2
      VarS2    =  ((EX4 - 4*EX3*EX1 + 6*EX2*EX1**2 - 3*EX1**4) - 
     1                MuS2**2)/N2 + 2.d0/((N2-1.d0)*n2)*MuS2**2

       alpha    =  MuS2**2/VarS2
       beta     =  MuS2/alpha

      CovMS2   =  (EX3 - 3*EX2*EX1 + 2*EX1**3)/SQRT(N2*(N2-1.D0))
      
      MuS      =  SQRT(beta) * exp(dlngam(alpha+0.5d0)-dlngam(alpha))
      CovMS    =  CovMS2/(2*MuS)
      VarS     =  MuS2 - MuS**2
      
        Lambda = CovMS/VarS
        EtaP   = Eta + Lambda
        MuMP   = MuM - Lambda * MuS
        VarMP  = VarM - CovMS**2/VarS
        df     = 2.d0*alpha
        ncp    = (MuMP -zr)/SQRT(VarMP)
        q      = -sqrt(MuS2/VarMP) * EtaP
       FGGB = 1.d0 - FP_TNC_CDF(q,df,ncp)

      return
      end
c*_*-*~*_*-*~*             new program begins here            *_*-*~*_*-*~*_*-*~
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c****|subroutine plotposHS
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c    
c    this routine provides plotting positions for censored data 
c      using the Hirsch/Stedinger plotting position formula (hirsch, 1987)
c      Purpose is operational version of peakfq (bulletin 17b
c      implementation) which would include the expected moments
c      algorithm (ema; cohn et al. 1997; 2001)
c
c    N.B:  The returned values are "exceedance probabilities" which are 
c             equal to 1-NonExceedance probabilities
c       
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c    
c    References
c    ----------
c
c    @article{hirsch1987plotting,
c      title={{Plotting positions for historical floods and their precision}},
c      author={Hirsch, R.M. and Stedinger, J.R.},
c      journal={Water Resources Research},
c      volume={23},
c      number={4},
c      pages={715--727},
c      year={1987}
c    }
c
c    @article{cohn1997algorithm,
c      title={{An algorithm for computing moments-based flood quantile estimates when historical flood information is available}},
c      author={Cohn, TA and Lane, WL and Baier, WG},
c      journal={Water Resources Research},
c      volume={33},
c      number={9},
c      pages={2089--2096},
c      year={1997}
c    }
c
c    @article{cohn2001confidence,
c      title={{Confidence intervals for Expected Moments Algorithm flood quantile estimates}},
c      author={Cohn, T.A. and Lane, W.L. and Stedinger, J.R.},
c      journal={Water resources research},
c      volume={37},
c      number={6},
c      pages={1695--1706},
c      year={2001}
c    }
c
c
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c    development history
c
c    timothy a. cohn        15 apr 2010
c       modified            16 apr 2010
c
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c    subroutine calls
c
c    ppfit                  available in probfun.f
c    arrange                available in probfun.f
c    porder                 available in probfun.f
c    ppsolve                available in probfun.f
c    plpos                  available in probfun.f
c    shellsort              available in probfun.f
c
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c       n.b.  use of plotposHS requires two distinct types of information:
c
c            a.  "observations," which are described by {ql,qu};
c                  specifically, the true value of q is assumed to lie
c                  wihin the interval ql(i) <= q(i) <= qu(i);
c                  if q(i) is known exactly, then ql(i) = q(i) = qu(i);
c
c                  {ql,qu} are required by ema to fit p3 distn to data
c
c            b.  "censoring pattern," which is described by {tl,tu}
c                  {tl,tu} define the interval on which data get observed:
c                  if q is not inside {tl,tu}, q is either left or right 
c                  censored.
c
c                  examples:
c                    {-inf,+inf}  => systematic data (exact value of q
c                      would be measured)
c                    {t,   +inf}  => historical data 
c                      q>=t would be known exactly; 
c                      q<t would be reported as censored ("less-than")
c                
c                    {tl,tu} are required by ema to estimate confidence 
c                      intervals; they are not required for frequency
c                      curve determination
c
c       n.b.  broken record is represented in the time-series data by a
c             threshold value, tl(i), equal to 1.d12 or larger.  This
c             indicates that only flows greater than 1.d12 would have 
c             been recorded.  No such floods have ever occurred on the
c             planet [TAC, 2010, personal communication].
c
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c       input variables:
c       ------------------------------------------------------------------------
c            n          i*4  number of observations (censored, uncensored, or 
c                              other)
c            ql(n)      r*8  vector of lower bounds on floods
c            qu(n)      r*8  vector of upper bounds on floods
c            tl(n)      r*8  vector of lower bounds on flood threshold
c            tu(n)      r*8  vector of upper bounds on flood threshold
c            alpha      r*8  alpha value to be used in defining plotting
c                              position pp(i) = (i-alpha)/(n+1-2*alpha)
c
c       n.b.  input data are the same as the first 5 arguments to emafit
c
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c       output variables:
c       ------------------------------------------------------------------------
c
c            q(n)       r*8  sorted discharges
c            peX(n)     r*8  estimated exceedance probabilities
c            nt         i*4  number of censoring thresholds
c            thr(nt)    r*8  vector of censoring thresholds
c            peT(nt)    r*8  estimated exceedance probabilities of thresholds
c            nb(nt)     i*4  number of observations censored at each threshold
c
c       n.b.  output data correspond to every observation that is not identified
c             as broken record.  Censored data are assigned relative plotting
c             positions according to their order in the input time series.
c             
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
      subroutine plotposHS(n,ql,qu,tl,tu,alpha,q,peX,nt,thr,peT,nb)

      implicit none
      
      integer
     1 n,i,j,n_true,nmax
     
C JFE      parameter (nmax=10000)
      parameter (nmax=20000)
      
      integer
     1 ix(nmax),nt,nb(nmax)
      
      double precision
     1 ql(n),qu(n),tl(n),tu(n),peX(*),peT(*),infty,
     2 q(nmax),t(nmax),mu,sigma,pp(nmax),xsynth(nmax),
     3 alpha,thr(nmax)
     
      data infty/1.d12/
     
        n_true = 0
      do 10 i=1,n
        if(tl(i) .lt. infty) then
          n_true = n_true + 1
c  for interval data, use average of endpoints to determine order for plotting
          q(n_true) = (ql(i)+qu(i))/2.d0
          t(n_true) = tl(i)
          ix(n_true) = i
        else
          peX(i) = -99.d0
        endif
10    continue

      
      call ppfit2(n_true,q,t,mu,sigma,pp,xsynth,
     1            alpha,nt,thr,peT,nb)
      
      do 20 j=1,n_true
        peX(ix(j)) = 1.d0-pp(j)
20    continue
      
      return
      end
C*_*-*~*_*-*~*             NEW PROGRAM BEGINS HERE            *_*-*~*_*-*~*_*-*~
         SUBROUTINE PPFIT2(N,X,XD,MU,SIGMA,PPALL,XSYNTH,
     1       ALPHA,NT,THR,PE_OFF,NBT)

C===============================================================================
C
C        SUBROUTINE TO FIT LN2 DISTRIBUTION TO MULTIPLY-CENSORED DATA
C
C        AUTHOR........TIM COHN
C        DATE..........APRIL 7, 1986
C        MODIFIED......12 FEBRUARY 2003 (TAC)
C         --REMOVED ERRORS INTRODUCED IN 1986 REVISION TO CODE
C         --IMPLICIT NONE ADDED; ALL VARIABLES DECLARED
C         --DOUBLE PRECISION
C        MODIFIED......10 JUNE 2011 (TAC)
C         --COMPUTE PLOTTING POSITIONS FOR ''GLYPHS''
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        N         I*4       INPUT NUMBER OF OBSERVATIONS (NO-DETECTS COUNT)
C        X(N)      R*8       INPUT VECTOR CONTAINING DATA;  '0' INDICATES
C                              NO-DETECT
C        XD(N)     R*8       INPUT VECTOR CONTAINING CENSORING THRESHOLD FOR 
C                              EACH PT.
C        MU        R*8       OUTPUT FITTED VALUE OF MU
C        SIGMA     R*8       OUTPUT FITTED VALUE OF SIGMA
C        PPALL(N)  R*8       OUTPUT VECTOR OF PLOTTING POSITIONS 
C                              CORRESPONDING TO X
C                              (SEE HELSEL & COHN, WRR 24(12), 1988)
C        XSYNTH(N) R*8       OUTPUT VECTOR OF REAL AND SYNTHETIC OBSERVATIONS
C                              XSYNTH(I) = X(I);  X(I)>=XD(I) 
C                              XSYNTH(I) = EXP(MU+Z(PPALL(I))*SIGMA); X(I)<XD(I)
C        NT        I*4       OUTPUT NUMBER OF DISTINCT THRESHOLD
C        THR(NT)   R*8       OUTPUT VECTOR OF THRESHOLDS
C        PE_OFF(NT)R*8       OUTPUT VECTOR NON-EXCEEDANCE PROBS CORR. TO
C                              THRESHOLDS; NOTE THAT PE(0) HAS BEEN OMITTED
C        NBT(NT)    I*4       OUTPUT NO. NO-DETECTS CENSORED AT THR(I)
C                              (OBSVS SPECIFIED AS "<THR(I)")
C
C===============================================================================
C
C        N.B.   1)  "SYNTHETIC" VALUES ARE ASSIGNED TO THE "LESS-THAN" VALUES.
C               2)  THE "SYNTHETIC VALUES" MAY EXCEED THE CENSORING THRESHOLD
C                   (THIS OCCURS FOR TWO REASONS:
C                    1)  AN LN2 MODEL IS EMPLOYED TO MODEL THE CENSORED DATA
C                    2)  THE PLOTTING POSITION AND MR FORMULAS ASSUME A
C                        TYPE-II CENSORING MODEL.  IN FACT, CENSORING OF WATER
C                        QUALITY DATA ALMOST INVARIABLY RESULTS FROM A TYPE-I
C                        CENSORING PROCESS)   
C
C===============================================================================
C
C        NT        I*4       NUMBER OF THRESHOLDS
C        NBT(NT)   I*4       NUMBER OF OBSERVATIONS BELOW EACH THRESHOLD
C        ND(NT)    I*4       NUMBER OF DETECTS ABOVE THIS THR. AND BELOW NEXT
C        NVAL      I*4       TOTAL NUMBER OF DETECTS
C
C===============================================================================

         IMPLICIT NONE

         INTEGER N_PTS,N_T
C JFE
C         PARAMETER (N_T=10000,N_PTS=10000)
         PARAMETER (N_T=20000,N_PTS=20000)

         DOUBLE PRECISION
     1     X(*),XD(*),MU,SIGMA,
     2     Y(N_PTS),PP(N_PTS),ALPHA,
     3     PPC(N_PTS),PPALL(*),XSYNTH(*),
     4     THR(N_T),PE(0:N_T),PE_OFF(N_T)
     
         INTEGER
     1     N,ND(N_T),NBT(N_T),NB(N_T),NT,NVAL,I

         CALL ARRANGE2(X,XD,N,NT,ND,NB,NVAL,Y,THR,NBT)

         CALL PPLOT2(NT,ND,NB,ALPHA,PP,PPC,PE)
           DO 10 I=1,NT
             PE_OFF(I) = PE(I)      ! ELIMINATE PE(0), WHICH EQUALS 1.0
   10      CONTINUE
         
         CALL PPSOLVE(Y,PP,NVAL,MU,SIGMA)
         
         CALL PLPOS(N,X,XD,PP,PPC,MU,SIGMA,PPALL,XSYNTH)

         RETURN
         END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
         SUBROUTINE ARRANGE2(X,XD,N,NT,ND,NB,NVAL,Y,THRESH,NBT)
C===============================================================================
C
C        SUBROUTINE DESIGNED TO ARRANGE DATA FOR PPFIT, WHICH EMPLOYS DATA TO
C        FIT LN2 DISTRIBUTION TO MULTIPLY-CENSORED DATA
C
C        AUTHOR........TIM COHN
C        DATA..........APRIL 7, 1986
C        MODIFIED......10 JUNE 2011 (TAC)
C         --COMPUTE PLOTTING POSITIONS FOR ''GLYPHS''
C        MODIFIED......12 FEBRUARY 2003 (TAC)
C         --REMOVED ERRORS INTRODUCED IN 1986
C           WHEN CODE WAS REVISED BY DRH.  
C         --IMPLICIT NONE ADDED; ALL VARIABLES DECLARED
C         --DOUBLE PRECISION
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        X(N)      R*8       INPUT VECTOR CONTAINING DATA;  '0' INDICATES 
C                              NO-DETECT
C        XD(N)     R*8       INPUT VECTOR CONTAINING CENSORING THRESHOLD FOR 
C                              EACH PT.
C        N         I*4       INPUT NUMBER OF OBSERVATIONS (NO-DETECTS COUNT)
C        NT        I*4       OUTPUT NUMBER OF THRESHOLDS
C        NB(NT)    I*4       OUTPUT TOTAL NO. OBS. BELOW EACH THRESHOLD
C        NBT(NT)   I*4       OUTPUT NO. OBS. BELOW SPECIFIC THRESHOLD
C        ND(NT)    I*4       OUTPUT NUMBER OF DETECTS ABOVE THIS THR. AND 
C                              BELOW NEXT
C        NVAL      I*4       OUTPUT TOTAL NUMBER OF DETECTS
C        Y(N)      R*8       OUTPUT VECTOR OF SORTED DETECTS
C        THRESH(NT)R*8       OUTPUT VECTOR OF SORTED THRESHOLDS
C
C===============================================================================

         IMPLICIT NONE

         DOUBLE PRECISION
     1      X(*),XD(*),Y(*),
     2      THRESH(200),PLUSINF
     
         INTEGER
C JFE
C     1      N,NT,NVAL,ND(*),NBT(*),NB(*),IP(10000),
     1      N,NT,NVAL,ND(*),NBT(*),NB(*),IP(20000),
     2      I,ICT,J

         DATA PLUSINF/1.0D30/

         CALL PORDER(XD,N,IP)

           ICT          =  0
           J            =  1
           THRESH(1)    =  XD(IP(1))
           NBT(1)       =  0

         DO 10 I=1,N

            IF(XD(IP(I)) .NE. THRESH(J)) THEN
              J    =  J+1
              THRESH(J) =  XD(IP(I))
              ND(J)     =  0
              NBT(J)    =  0
            ENDIF
              
            IF(X(IP(I)) .GE. THRESH(J)) THEN
              ICT       =  ICT+1
              Y(ICT)    =  X(IP(I))
            ELSE
              NBT(J)    =  NBT(J)+1
            ENDIF
   10    CONTINUE

         NT             =  J
         NVAL           =  ICT
         THRESH(NT+1)   =  PLUSINF

         CALL SHELLSORT(Y,ICT)
                
           J       =  1
           ND(1)   =  0
         DO 20 I=1,NVAL
   30      IF(Y(I) .GE. THRESH(J+1)) THEN
              J    =  J+1
              ND(J)=  0
              GOTO 30
           ELSE
              ND(J)=  ND(J)+1
           ENDIF
   20    CONTINUE

           NB(1)       =  NBT(1)
         DO 40 I=2,NT
           NB(I)   =  NB(I-1)+NBT(I)+ND(I-1)
   40    CONTINUE

         RETURN
         END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
         SUBROUTINE PPLOT2(NT,ND,NB,ALPHA,PP,PPC,PE)
C===============================================================================
C
C        SUBROUTINE TO FIND PLOTTING POSITIONS FOR DATA CENSORED AT MULTIPLE 
C          THRESHOLDS
C
C        AUTHOR........TIM COHN
C        DATE..........APRIL 7, 1986
C        MODIFIED......12 FEBRUARY 2003 (TAC)
C         --REMOVED ERRORS INTRODUCED IN 1986
C           WHEN CODE WAS REVISED BY DRH.  
C         --IMPLICIT NONE ADDED; ALL VARIABLES DECLARED
C         --DOUBLE PRECISION
C        MODIFIED......15 JUNE 2011 (TAC)
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        NT        I*4       INPUT NUMBER OF THRESHOLDS
C        ND(NT)    I*4       INPUT NUMBER OF DETECTS ABOVE THIS THR. AND BELOW 
C                              NEXT
C        NB(NT)    I*4       INPUT NUMBER OF OBSERVATIONS BELOW EACH THRESHOLD
C        ALPHA     R*8       INPUT PLOTTING POSITION ALPHA (I-A)/(N+1-2A)
C        PP(NA)    R*8       OUTPUT VECTOR OF PLOTTING POSITIONS
C        PPC(NB)   R*8       OUTPUT VECTOR OF PLOTTING POSITIONS CORRESONDING TO 
C                              X<XD
C        PE(0:NT)  R*8       OUTPUT VECTOR OF PLOTTING POSITIONS CORRESONDING TO 
C                              THRESHOLDS
C
C===============================================================================

         IMPLICIT NONE

         INTEGER N_T
C JFE
C
C         PARAMETER (N_T=10000)
         PARAMETER (N_T=20000)

         DOUBLE PRECISION
     1     PP(*),PE(0:N_T),ALPHA,PD,PPT,PPC(*)
     
         INTEGER
     1     ND(*),NB(*),NT,
     2     ISUM,I,J,JSUM,NBT


           PE(NT+1)  =  0.D0
         DO 10 I=NT,1,-1
           PD        =  ND(I)/FLOAT(MAX(1,ND(I))+NB(I))
           PE(I)     =  PE(I+1) + (1.D0-PE(I+1)) * PD
   10    CONTINUE

           ISUM    =  0
           JSUM    =  0
           PE(0)   =  1.00D0

         DO 20 I=1,NT
           DO 30 J=1,ND(I)
               PPT  =  (J-ALPHA)/(ND(I) + 1.0 - 2.0*ALPHA)
               PP(ISUM+J)    =  (1.D0-PE(I)) + (PE(I)-PE(I+1))*PPT
   30      CONTINUE
               ISUM =  ISUM+ND(I)

           IF(I .EQ. 1) THEN
              NBT = NB(1)
           ELSE
              NBT = NB(I)-NB(I-1)-ND(I-1)
           ENDIF
           
           DO 40 J=1,NBT
             PPT  = (J-ALPHA)/(NBT+1.0-2*ALPHA)
             PPC(JSUM+J) = (1-PE(I))*PPT
   40      CONTINUE
             JSUM=JSUM+NBT

   20    CONTINUE

         RETURN
         END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
         SUBROUTINE PPSOLVE(Y,PP,NVAL,MU,SIGMA)
C===============================================================================
C
C        SUBROUTINE TO SOLVE FOR LEAST SQUARES ESTIMATES OF MU AND SIGMA
C        ASSUMING A LINEAR MODEL
C
C        AUTHOR........TIM COHN
C        DATE..........APRIL 7, 1986
C        MODIFIED......12 FEBRUARY 2003 (TAC)
C         --REMOVED ERRORS INTRODUCED IN 1986
C           WHEN CODE WAS REVISED BY DRH.  
C         --IMPLICIT NONE ADDED; ALL VARIABLES DECLARED
C         --DOUBLE PRECISION
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        Y(N)      R*8       INPUT VECTOR CONTAINING SORTED DETECT-DATA
C        PP(NA)    R*8       INPUT VECTOR OF PLOTTING POSITIONS
C        NVAL      I*4       INPUT TOTAL NUMBER OF DETECTS
C        MU        R*8       OUTPUT FITTED VALUE OF MU
C        SIGMA     R*8       OUTPUT FITTED VALUE OF SIGMA
C
C===============================================================================

         IMPLICIT NONE

         INTEGER N_PTS
C JFE
C         PARAMETER (N_PTS=10000)
         PARAMETER (N_PTS=20000)

         DOUBLE PRECISION
     1     Y(*),PP(*),MU,SIGMA,
     2     PSIG(N_PTS),LOGY(N_PTS),
     3     FP_Z_ICDF,
     4     SUMY,SUMP,SUMPY,SUMP2,YBAR,PBAR
     
         INTEGER
     1     I,NVAL

              SUMY =  0.D0
              SUMP =  0.D0
              SUMPY=  0.D0
              SUMP2=  0.D0

         DO 10 I=1,NVAL
              PSIG(I)   =  FP_Z_ICDF(PP(I))
              LOGY(I)   =  LOG(Y(I))
              SUMY      =  SUMY    +  LOGY(I)
              SUMPY     =  SUMPY   +  LOGY(I)*PSIG(I)
              SUMP      =  SUMP    +  PSIG(I)
              SUMP2     =  SUMP2   +  PSIG(I)**2
   10    CONTINUE
              YBAR      =  SUMY/NVAL
              PBAR      =  SUMP/NVAL

              SIGMA     =  (SUMPY-NVAL*PBAR*YBAR)/(SUMP2-NVAL*PBAR**2)
              MU        =  YBAR - SIGMA*PBAR
         RETURN
         END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
         SUBROUTINE PLPOS(N,X,XD,PP,PPC,MU,SIGMA,PPALL,XSYNTH)
C===============================================================================
C
C        SUBROUTINE PLPOS GIVES THE PLOTTING POSITIONS FOR ALL OF THE
C          OBSERVATIONS, CENSORED AND UNCENSORED, IN THE ORIGINAL INPUT
C          DATA SET
C
C
C        AUTHOR....TIM COHN
C        DATE......12 FEBRUARY 2003 (TAC)
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        N         I*4       INPUT NUMBER OF OBSERVATIONS (NO-DETECTS COUNT)
C        X(N)      R*8       INPUT VECTOR CONTAINING DATA;  '0' INDICATES 
C                              NO-DETECT
C        XD(N)     R*8       INPUT VECTOR CONTAINING CENSORING THRESHOLD FOR 
C                              EACH PT.
C        PP(NA)    R*8       INPUT VECTOR OF PLOTTING POSITIONS CORRESONDING TO 
C                              X>=XD
C                              (SEE HELSEL & COHN, WRR 24(12), 1988)
C        PPC(NB)   R*8       INPUT VECTOR OF PLOTTING POSITIONS CORRESONDING TO 
C                              X<XD
C        MU        R*8       INPUT FITTED VALUE OF MU
C        SIGMA     R*8       INPUT FITTED VALUE OF SIGMA
C        PPALL(N)  R*8       OUTPUT VECTOR OF PLOTTING POSITIONS CORRESONDING TO X
C                              (SEE HELSEL & COHN, WRR 24(12), 1988)
C        XSYNTH    R*8       OUTPUT VECTOR OF REAL AND SYNTHETIC OBSERVATIONS
C                              XSYNTH(I) = X(I);                       
C                              X(I)>=XD(I) 
C                              XSYNTH(I) = EXP(MU+Z(PPALL(I))*SIGMA);  
C                              X(I)<XD(I)
C
C===============================================================================

         IMPLICIT NONE
         
         INTEGER N_PTS
C JFE
C
C         PARAMETER (N_PTS=10000)
         PARAMETER (N_PTS=20000)

         DOUBLE PRECISION 
     1     X(*),XD(*),MU,SIGMA,PPALL(*),PP(*),PPC(*),XSYNTH(*),
     2     FP_Z_ICDF
         
         INTEGER
     1     IX(N_PTS),IXD(N_PTS),N,NA,NB,I

         CALL PORDER(X,N,IX)
         CALL PORDER(XD,N,IXD)
         
           NA  =  0
           NB  =  0
         DO 10 I=1,N
           IF(X(IX(I)) .GE. XD(IX(I))) THEN
             NA            =  NA+1
             PPALL(IX(I))  =  PP(NA)
             XSYNTH(IX(I)) =  X(IX(I))
           ENDIF
           IF(X(IXD(I)) .LT. XD(IXD(I))) THEN
             NB            =  NB+1
             PPALL(IXD(I)) =  PPC(NB)
             XSYNTH(IXD(I))=  EXP(MU+FP_Z_ICDF(PPC(NB))*SIGMA)
           ENDIF
10       CONTINUE

         RETURN
         END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
         SUBROUTINE PORDER(X,N,IX)
C===============================================================================
C
C        SUBROUTINE PORDER GIVES THE PERMUTATION OF A USER-SUPPLIED
C        VECTOR.  TREE-SORT ALGORITHM IS EMPLOYED (WHY NOT?)
C        RETURN RESULTING IX ARRAY WITH THE ASCENDING PERMUTATIONS
C
C          X(IX(1)) = SMALLEST VALUE IN X
C          X(IX(N)) = LARGEST VALUE IN X
C
C
C        AUTHOR....TIM COHN
C        DATE......APRIL 1,  1986
C        REVISED...AUGUST 9, 1986      (TAC)
C        MODIFIED......12 FEBRUARY 2003 (TAC)
C         --REMOVED ERRORS INTRODUCED IN 1986
C           WHEN CODE WAS REVISED BY DRH.  
C         --IMPLICIT NONE ADDED; ALL VARIABLES DECLARED
C         --DOUBLE PRECISION
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        X(N)      R*8       INPUT VECTOR OF UNORDERED DATA
C        N         I*4       INPUT NUMBER OF OBSERVATIONS IN X
C        IX(N)     I*4       OUTPUT VECTOR CONTAINING PERMUTATION OF X
C                               X(IX(1)) IS THE SMALLEST VALUE OF X
C                                 . . . 
C                               X(IX(N)) IS THE LARGEST VALUE OF X
C
C===============================================================================

         IMPLICIT NONE
         
         INTEGER N_PTS
C JFE
C
C         PARAMETER (N_PTS=10000)
         PARAMETER (N_PTS=20000)

         DOUBLE PRECISION X(*)
         
         INTEGER
     1     L(0:N_PTS),R(0:N_PTS),P(0:N_PTS),
     2     IX(*),I2,I,INDX,N,ICT

C===============================================================================
C
C    FIRST CHECK TO SEE IF WE HAVE AN ORDERED DATA VECTOR TO BEGIN WITH
C
         DO 50 I2=2,N
              IX(I2)    =  I2
              IF(X(I2) .LT. X(I2-1)) GOTO 1
   50    CONTINUE
              IX(1)     =  1
              RETURN

    1    CONTINUE

         L(1) =  0
         R(1) =  0
         P(1) =  0

        DO 10 I=2,N
         INDX =  1
         L(I)   =  0
         R(I)   =  0

   20    CONTINUE
         IF(X(I) .GE. X(INDX)) THEN
              IF(R(INDX) .EQ. 0) THEN
                   R(INDX)   =  I
                   P(I)      =  INDX
                   GOTO 10
              ELSE
                   INDX =  R(INDX)
                   GOTO 20
              ENDIF
         ELSE
              IF(L(INDX) .EQ. 0) THEN
                   L(INDX)   =  I
                   P(I)      =  INDX
                   GOTO 10
              ELSE
                   INDX =  L(INDX)
                   GOTO 20
              ENDIF
         ENDIF
   10   CONTINUE

          INDX =  1
         DO 40 ICT=1,N

   30    CONTINUE
         IF(L(INDX) .EQ. 0) THEN
              IX(ICT)     =  INDX
              P(R(INDX))   =  P(INDX)
              L(P(INDX))   =  R(INDX)
              INDX         =  P(INDX)
         ELSE
              INDX =  L(INDX)
              GOTO 30
         ENDIF
   40    CONTINUE
         RETURN
         END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
         SUBROUTINE SHELLSORT(X,N)
C===============================================================================
C
C        SUBROUTINE TO SORT A VECTOR IN PLACE
C
C        AUTHOR........TIM COHN
C        DATE..........12 FEBRUARY 2003 (TAC)
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        X(N)      R*8       INPUT/OUTPUT VECTOR OF DATA
C        N         I*4       INPUT SAMPLE SIZE
C
C===============================================================================

      IMPLICIT NONE
      
      DOUBLE PRECISION X(*),R
      
      INTEGER N,M,I,J,S,Z
      
      M=N
  100 M=INT(M/2)
      IF (M.EQ.0) GOTO 210
      DO 190 S=1,M
  130 I=S
      J=S+M
      Z=0
  140 IF (X(I).LE.X(J)) GOTO 160
      Z=1
      R=X(I)
      X(I)=X(J)
      X(J)=R
  160 I=J
      J=J+M
      IF (J.LT.(N+1)) GOTO 140
      IF (Z.EQ.1) GOTO 130
  190 CONTINUE
      GOTO 100
  210 RETURN
      END
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c    silly little function to print asterisks for p-value significance
c
c    timothy a. cohn        9 mar 2011
c
      character*4 function signif(p)
      double precision p
      if(p .lt. 0.01) then
        signif = ' ***'
      else if(p .lt. 0.05d0) then
        signif = ' ** '
      else if(p .lt. 0.10d0) then
        signif = ' *  '
      else if(p .ge. 0.10d0) then
        signif = '    '
      endif
        return
      end
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
c
c    another silly little function to grab last 4 characters in filename
c      (usually the ".ext") and make it available.
c
c    timothy a. cohn        9 mar 2011
c
      character*4 function filext(fname)
      character*(*) fname
      integer ichar
        ichar=len_trim(fname)
      filext=fname(max(1,(ichar-3)):ichar)
      return
      end
c****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////

C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
       SUBROUTINE NORMQUAD(N,T,W)
C===============================================================================
C
C        SUBROUTINE NORMQUAD COMPUTES QUADRATURE NODES AND WEIGHTS 
C          CORRESPONDING TO THE STANDARD NORMAL DISTRIBUTION
C
C          INT(-INF,INF) g(z)f(z)dz = sum(i=1,n) w(i)*g(t(i))
C
C
C        AUTHOR....TIM COHN
C        DATE......05 OCTOBER 2012
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        N         I*4       INPUT NUMBER OF NODES/ORDER OF OTHOGONAL POLYNOMIAL
C        T(N)      R*8       OUTPUT VECTOR OF NODES
C        W(N)      R*8       OUTPUT VECTOR OF WEIGHTS
C
C===============================================================================
C           THIS SET OF ROUTINES COMPUTES THE NODES T(J) AND WEIGHTS
C        W(J) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
C
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX
C
C                              N
C        BY                   SUM W  F(T )
C                             J=1  J    J
C
C        (NOTE W(X) AND W(J) HAVE NO CONNECTION WITH EACH OTHER.)
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
C
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
C        ORTHOGONAL POLYNOMIALS.  THE NODES T(J) ARE JUST THE ZEROES
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.
C
C     INPUT PARAMETERS (ALL REAL NUMBERS ARE IN DOUBLE PRECISION)
C
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
C                 QUADRATURE RULE:
C
C        KIND = 1:  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
C        KIND = 2:  CHEBYSHEV QUADRATURE OF THE FIRST KIND
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
C        KIND = 3:  CHEBYSHEV QUADRATURE OF THE SECOND KIND
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)
C        KIND = 4:  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
C                   (-INFINITY, +INFINITY)
C        KIND = 5:  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.
C                   NOTE: KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
C        KIND = 6:  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1
C===============================================================================
C       
       IMPLICIT NONE
       
       DOUBLE PRECISION T(*),W(*),ENDPTS(2),B(1000),SQ2,SQPI
       INTEGER I,KIND,KPTS,N
       
       DATA KPTS/0/,ENDPTS/2*0.D0/
       DATA SQ2/1.41421356237/,SQPI/0.56418958354/
       
	   	   !SAS: Initialize B
	   DO I = 1,1000
	      B(I) = 0D0
       END DO

         KIND=4
       CALL GAUSSQ(KIND, N, 0.D0, 0.D0, KPTS, ENDPTS, B, T, W)
       DO 10 I=1,N
         T(I) = T(I)*SQ2
         W(I) = W(I)*SQPI
10     CONTINUE
       
       RETURN
       END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
       SUBROUTINE GAMMAQUAD(N,ALPHA,BETA,T,W)
C===============================================================================
C
C        SUBROUTINE GAMMAQUAD COMPUTES QUADRATURE NODES AND WEIGHTS 
C          CORRESPONDING TO A GAMMA DISTRIBUTION WITH PARAMETERS
C          ALPHA, BETA
C
C          INT(0,INF) g(y)f(y)dy = sum(i=1,n) w(i)*g(t(i))
C
C
C        AUTHOR....TIM COHN
C        DATE......05 OCTOBER 2012
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        N         I*4       INPUT NUMBER OF NODES/ORDER OF OTHOGONAL POLYNOMIAL
C        ALPHA     R*8       SHAPE PARAMETER OF GAMMA DISTRIBUTION
C        BETA      R*8       SCALE PARAMETER OF GAMMA DISTRIBUTION
C        T(N)      R*8       OUTPUT VECTOR OF NODES
C        W(N)      R*8       OUTPUT VECTOR OF WEIGHTS
C
C===============================================================================
C       
       IMPLICIT NONE
       
       DOUBLE PRECISION 
     1   T(*),W(*),ALPHA,BETA,ENDPTS(2),B2(1000),SQ2,SQPI,DGI,
     2   ALPHAMAX,A,B,C
       INTEGER I,KIND,KPTS,N
       DATA ALPHAMAX/160.D0/
       
       DATA KPTS/0/,ENDPTS/2*0.D0/
       DATA SQ2/1.41421356237/,SQPI/0.56418958354/
	   
	   !SAS: Initialize B2
	   DO I = 1,1000
	      B2(I) = 0D0
       END DO
       
       IF(ALPHA .LE. ALPHAMAX) THEN
         A = ALPHA
         B = BETA
         C = 0.D0
       ELSE
         A = ALPHAMAX
         B = BETA*SQRT(ALPHA/A)
         C = ALPHA*BETA-A*B
       ENDIF
           KIND=6
         CALL GAUSSQ(KIND, N, A-1.D0, 0.D0, KPTS, ENDPTS, B2, T, W)
           DGI = 1.D0/gamma(A)
         DO 10 I=1,N
           T(I) = C + T(I)*B
           W(I) = W(I)*DGI
10       CONTINUE
       
       RETURN
       END
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
       SUBROUTINE CHOL33(S,V,IFLAG)
C===============================================================================
C
C        SUBROUTINE CHOL33 COMPUTES A SORT-OF CHOLESKY DECOMPOSITION (NOT LU)
C          OF A SYMMETRIC POSITIVE-DEFINITE 3x3 MATRIX WITH THE FOLLOWING
C          PROPERTIES:
C            
C          ASSUME WE START WITH A (Nx3) VARIABLE X (THINK {M,S2,G})
C          X HAS 3x3 COVARIANCE MATRIX S AND P=CHOL33(S,P)
C            --THE COLUMNS OF XxP ARE ORTHOGONAL
c            --THE SECOND COLUMN OF P HAS TWO ZEROES (CHOLESKY ITS COL 1)
C
C        N.B. THIS ROUTINE IS USEFUL FOR ONLY ONE PURPOSE: GENERATING
C             SAMPLES USING QUADRATURE POINTS WHERE THE MIDDLE VARIABLE
C             IS GAMMA AND THE TWO OUTER ONES ARE BOTH NORMAL
C             NO CHECK IS DONE TO VERIFY THAT MATRIX HAS NON-ZERO
C             MAIN DIAGONAL OR ANYTHING ELSE...
C
C        FOR THE ALGEBRA, SEE COHN [2012] (NOTEBOOK 3)
C          
C        AUTHOR....TIM COHN
C        DATE......05 OCTOBER 2012
C
C===============================================================================
C
C     PROPERTY OF US GOVERNMENT, U.S. GEOLOGICAL SURVEY
C
C     *** DO NOT MODIFY WITHOUT AUTHOR''S CONSENT ***
C
C     AUTHOR CAN BE CONTACTED AT:  TACOHN@USGS.GOV (703/648-5711)     
C
C===============================================================================
C
C        S(3,3)    R*8       INPUT MATRIX
C        V(3,3)    R*8       OUTPUT MATRIX, MODIFIED CHOLESKY DECOMPOSITION OF S
C        IFLAG     I*4       OUTPUT FLAG; 0 IF OK; 1 IF S NOT PSD
C
C  THE OUTPUT IS NOT UPPER TRIANGULAR, BUT RATHER UT ONCE COLUMNS 1 AND 2 ARE 
C  SWAPPED AND ROWS 1 AND 2 ARE SWAPPED. 
C
C        > chol(var(x))
C                 [,1]      [,2]        [,3]
C        [1,] 0.938117 0.1361061 -0.04404994
C        [2,] 0.000000 0.9575892 -0.07947178
C        [3,] 0.000000 0.0000000  0.84128956
C        > chol33(var(x))
C                  [,1]      [,2]        [,3]
C        [1,] 0.9287822 0.0000000 -0.03242836
C        [2,] 0.1320116 0.9672135 -0.08487969
C        [3,] 0.0000000 0.0000000  0.84128956
C
C===============================================================================
C
      IMPLICIT NONE
      DOUBLE PRECISION S(3,3),V(3,3),T1
      INTEGER IFLAG
      
      IF(S(1,1).LE.0.D0 .OR. S(2,2).LE.0.D0 .OR. S(3,3).LE.0.D0) GOTO 99
      
      V(1,2)=0.D0
      V(3,2)=0.D0
      V(3,1)=0.D0
      V(2,2)=SQRT(S(2,2))
      V(2,1)=S(2,1)/V(2,2)
      V(2,3)=S(2,3)/V(2,2)
	  
	  !Initialize last entries before computation
	  V(1,1)=0.D0
      V(1,3)=0.D0
	  V(3,3)=0.D0
	  
	  !Now try to compute (success NOT guaranteed)
        T1=S(1,1)-V(2,1)**2
        IF(T1 .LT. 0.D0) GOTO 99
      V(1,1)=SQRT(T1)
      V(1,3)=(S(3,1)-V(2,3)*V(2,1))/V(1,1)
        T1=S(3,3)-V(2,3)**2-V(1,3)**2
        IF(T1 .LT. 0.D0) GOTO 99
      V(3,3)=SQRT(T1)
        IFLAG=0
      RETURN
99    CONTINUE
!      MATRIX S NOT POSITIVE SEMI-DEFINITE
!      VERY IMPORTANT THAT FLAG IS PROPERLY HANDLED (SAS)
        IFLAG=1 
        RETURN
      END
C****
      DOUBLE PRECISION FUNCTION MEANW (N,X,W)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C  COMPUTES THE WEIGHTED MEAN OF A VECTOR X
C
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION X(*),W(*),WSUM,XSUM

        WSUM = 0.D0
        XSUM = 0.D0
      DO 10 I=1,N
        WSUM = WSUM+W(I)
        XSUM = XSUM+W(I)*X(I)
10    CONTINUE
        MEANW=XSUM/WSUM
      RETURN
      END
C****
      DOUBLE PRECISION FUNCTION COVW (N,X,Y,W)
C****|===|====-====|====-====|====-====|====-====|====-====|====-====|==////////
C
C  COMPUTES THE WEIGHTED MEAN OF A VECTOR X
C
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION X(*),Y(*),W(*),WSUM,XSUM,MEANW,MX,MY

        MX  = MEANW(N,X,W)
        MY  = MEANW(N,Y,W)      
        WSUM = 0.D0
        XSUM = 0.D0
      DO 10 I=1,N
        WSUM = WSUM+W(I)
        XSUM = XSUM+W(I)*(X(I)-MX)*(Y(I)-MY)
10    CONTINUE
        COVW = XSUM/WSUM
      RETURN
      END
     
      

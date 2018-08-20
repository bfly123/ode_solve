	      SUBROUTINE RADCOV(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,NS,
     &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,SAFE,THET,QUOT1,QUOT2,NIT1,IJOB,STARTN,NIND1,
     &   NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1,NSMIN,NSMAX,NNMS,NM1NS,
     &   NMEE,IMPLCT,BANDED,LDJAC,LDE1,LDMAS,ZZ,Y0,SCAL,FF,FJAC,
     &   E1,EE2,FMAS,CONT,IP1,IP2,IPHES,VITU,VITD,HHOU,HHOD,
     &   NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,RPAR,IPAR)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR RADAU
C     PARAMETERS SAME AS IN RADAU WITH WORKSPACE ADDED
C ----------------------------------------------------------
C         DECLARATIONS
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NSDIM=7)
C --- THIS PARAMETER HAS TO BE CHANGED IF NUMBER OF STAGES IS >=7
      DIMENSION Y(N),ZZ(NNMS),FF(NNMS),Y0(N),SCAL(N)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),CONT(NNMS+N)
      DIMENSION ATOL(1),RTOL(1),RPAR(*),IPAR(*)
      DIMENSION E1(LDE1,NM1),EE2(LDE1,NMEE)
      DIMENSION ALPH(NSDIM),BETA(NSDIM),DD(NSDIM)
      DIMENSION ALPHN(NSDIM),BETAN(NSDIM)
      INTEGER IP1(NM1),IP2(NMEE/2),IPHES(NM1)
      COMMON/WEIGHT/NN,NSCON,XSOL,HSOL,C(0:NSDIM)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      COMMON /COE3/ T311,T312,T313,T321,T322,T323,T331,
     *    TI311,TI312,TI313,TI321,TI322,TI323,TI331,TI332,TI333
      COMMON /COE5/ T511,T512,T513,T514,T515,T521,T522,T523,T524,T525,
     *    T531,T532,T533,T534,T535,T541,T542,T543,T544,T545,T551,
     *    TI511,TI512,TI513,TI514,TI515,TI521,TI522,TI523,TI524,TI525,
     *    TI531,TI532,TI533,TI534,TI535,TI541,TI542,TI543,TI544,TI545,
     *    TI551,TI552,TI553,TI554,TI555
      COMMON /COE7/ T711,T712,T713,T714,T715,T716,T717,
     *              T721,T722,T723,T724,T725,T726,T727,
     *              T731,T732,T733,T734,T735,T736,T737,
     *              T741,T742,T743,T744,T745,T746,T747,
     *              T751,T752,T753,T754,T755,T756,T757,
     *              T761,T762,T763,T764,T765,T766,T767,T771,
     *              TI711,TI712,TI713,TI714,TI715,TI716,TI717,
     *              TI721,TI722,TI723,TI724,TI725,TI726,TI727,
     *              TI731,TI732,TI733,TI734,TI735,TI736,TI737,
     *              TI741,TI742,TI743,TI744,TI745,TI746,TI747,
     *              TI751,TI752,TI753,TI754,TI755,TI756,TI757,
     *              TI761,TI762,TI763,TI764,TI765,TI766,TI767,
     *              TI771,TI772,TI773,TI774,TI775,TI776,TI777
      LOGICAL REJECT,FIRST,IMPLCT,BANDED,CALJAC,STARTN,CALHES
      LOGICAL INDEX1,INDEX2,INDEX3,LAST,PRED,CHANGE,UNEXP,UNEXN,VARIAB
      EXTERNAL FCN,JAC
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
C -------- CHECK THE INDEX OF THE PROBLEM -----
      INDEX1=NIND1.NE.0
      INDEX2=NIND2.NE.0
      INDEX3=NIND3.NE.0
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
      IF (IMPLCT) CALL MAS(NM1,FMAS,LDMAS,RPAR,IPAR)
      VARIAB=NSMIN.LT.NSMAX
C ---------- CONSTANTS ---------
      EXPO=1.D0/(NS+1.D0)
      SQ6=DSQRT(6.D0)
      C31=(4.D0-SQ6)/10.D0
      C32=(4.D0+SQ6)/10.D0
      C31M1=C31-1.D0
      C32M1=C32-1.D0
      C31MC2=C31-C32
      DD1=-(13.D0+7.D0*SQ6)/3.D0
      DD2=(-13.D0+7.D0*SQ6)/3.D0
      DD3=-1.D0/3.D0
      N2=2*N
      N3=3*N
      N4=4*N
      N5=5*N
      N6=6*N
      UNEXP=.FALSE.
      UNEXN=.FALSE.
      CHANGE=.FALSE.
      IKEEP=0
      ICHAN=0
      THETA=0.0D0
      THETAT=0.0D0
      NN=N
      NNS=N*NS
      NSCON=NS
      LRC=NNS+N
      CALL COERTV(NSMAX)
      CALL COERCV(NS,C,DD,U1,ALPH,BETA)
      IF (M1.GT.0) IJOB=IJOB+10
      POSNEG=SIGN(1.D0,XEND-X)
      HMAXN=MIN(ABS(HMAX),ABS(XEND-X))
      IF (ABS(H).LE.10.D0*UROUND) H=1.0D-6
      H=MIN(ABS(H),HMAXN)
      H=SIGN(H,POSNEG)
      HOLD=H
      REJECT=.FALSE.
      FIRST=.TRUE.
      LAST=.FALSE.
      IF ((X+H*1.0001D0-XEND)*POSNEG.GE.0.D0) THEN
         H=XEND-X
         LAST=.TRUE.
      END IF
      HOPT=H
      FACCON=1.D0
      NSING=0
      XOLD=X
      IF (IOUT.NE.0) THEN
          IRTRN=1
          NRSOL=1
          XOSOL=XOLD
          XSOL=X
          DO I=1,N
             CONT(I)=Y(I)
          END DO
          NSOLU=N
          HSOL=HOLD
          CALL SOLOUT(NRSOL,XOSOL,XSOL,Y,CONT,LRC,NSOLU,
     &                RPAR,IPAR,IRTRN)
          IF (IRTRN.LT.0) GOTO 179
      END IF
      MLE=MLJAC
      MUE=MUJAC
      MBJAC=MLJAC+MUJAC+1
      MBB=MLMAS+MUMAS+1
      MDIAG=MLE+MUE+1
      MDIFF=MLE+MUE-MUMAS
      MBDIAG=MUMAS+1
      EXPMNS=(NS+1.0D0)/(2.0D0*NS)
      QUOTT=ATOL(1)/RTOL(1)
      IF (ITOL.EQ.0) THEN
          RTOL1=0.1D0*RTOL(1)**EXPMNS
          ATOL1=RTOL1*QUOTT
          DO I=1,N
             SCAL(I)=ATOL1+RTOL1*ABS(Y(I))
          END DO
      ELSE
          DO I=1,N
             QUOTT=ATOL(I)/RTOL(I)
             RTOL1=0.1D0*RTOL(I)**EXPMNS
             ATOL1=RTOL1*QUOTT
             SCAL(I)=ATOL1+RTOL1*ABS(Y(I))
          END DO
      END IF
      HHFAC=H
      CALL FCN(N,X,Y,Y0,RPAR,IPAR)
      NFCN=NFCN+1
C --- BASIC INTEGRATION STEP
  10  CONTINUE
C *** *** *** *** *** *** ***
C  COMPUTATION OF THE JACOBIAN
C *** *** *** *** *** *** ***
      NJAC=NJAC+1
      IF (IJAC.EQ.0) THEN
C --- COMPUTE JACOBIAN MATRIX NUMERICALLY
         IF (BANDED) THEN
C --- JACOBIAN IS BANDED
            MUJACP=MUJAC+1
            MD=MIN(MBJAC,M2)
            DO MM=1,M1/M2+1
               DO K=1,MD
                  J=K+(MM-1)*M2
 12               FF(J)=Y(J)
                  FF(J+N)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J))))
                  Y(J)=Y(J)+FF(J+N)
                  J=J+MD
                  IF (J.LE.MM*M2) GOTO 12
                  CALL FCN(N,X,Y,CONT,RPAR,IPAR)
                  J=K+(MM-1)*M2
                  J1=K
                  LBEG=MAX(1,J1-MUJAC)+M1
 14               LEND=MIN(M2,J1+MLJAC)+M1
                  Y(J)=FF(J)
                  MUJACJ=MUJACP-J1-M1
                  DO L=LBEG,LEND
                     FJAC(L+MUJACJ,J)=(CONT(L)-Y0(L))/FF(J+N)
                  END DO
                  J=J+MD
                  J1=J1+MD
                  LBEG=LEND+1
                  IF (J.LE.MM*M2) GOTO 14
               END DO
            END DO
         ELSE
C --- JACOBIAN IS FULL
            DO I=1,N
               YSAFE=Y(I)
               DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE)))
               Y(I)=YSAFE+DELT
               CALL FCN(N,X,Y,CONT,RPAR,IPAR)
               DO J=M1+1,N
                 FJAC(J-M1,I)=(CONT(J)-Y0(J))/DELT
               END DO
               Y(I)=YSAFE
            END DO
         END IF
      ELSE
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
         CALL JAC(N,X,Y,FJAC,LDJAC,RPAR,IPAR)
      END IF
      CALJAC=.TRUE.
      CALHES=.TRUE.
  20  CONTINUE
C --- CHANGE THE ORDER HERE IF NECESSARY
      IF (VARIAB) THEN
         NSNEW=NS
         ICHAN=ICHAN+1
         HQUOT=H/HOLD
         THETAT=MIN(10.0D0,MAX(THETA,THETAT*0.5))
         IF (NEWT.GT.1.AND.THETAT.LE.VITU.AND.
     &      HQUOT.LT.HHOU.AND.HQUOT.GT.HHOD) NSNEW=MIN(NSMAX,NS+2)
         IF (THETAT.GE.VITD.OR.UNEXP) NSNEW=MAX(NSMIN,NS-2)
         IF (ICHAN.GE.1.AND.UNEXN) NSNEW=MAX(NSMIN,NS-2)
         IF (ICHAN.LE.10) NSNEW=MIN(NS,NSNEW)
         CHANGE=NS.NE.NSNEW
         UNEXN=.FALSE.
         UNEXP=.FALSE.
         IF (CHANGE) THEN
            NS=NSNEW
            ICHAN=1
            NNS=N*NS
            NSCON=NS
            LRC=NNS+N
            CALL COERCV(NS,C,DD,U1,ALPH,BETA)
            EXPO=1.D0/(NS+1.D0)
            EXPMNS=(NS+1.0D0)/(2.0D0*NS)
            RTOL1=0.1D0*RTOL(1)**EXPMNS
            ATOL1=RTOL1*QUOTT
            IF (ITOL.EQ.0) THEN
               DO I=1,N
                  SCAL(I)=ATOL1+RTOL1*ABS(Y(I))
               END DO
            ELSE
               DO I=1,N
                  QUOTT=ATOL(I)/RTOL(I)
                  RTOL1=0.1D0*RTOL(I)**EXPMNS
                  ATOL1=RTOL1*QUOTT
                  SCAL(I)=ATOL1+RTOL1*ABS(Y(I))
               END DO
            END IF
         END IF
      END IF
C --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS
      FAC1=U1/H
      CALL DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)
      IF (IER.NE.0) GOTO 78
      DO K=1,(NS-1)/2
         ALPHN(K)=ALPH(K)/H
         BETAN(K)=BETA(K)/H
         IAD=(K-1)*2*NM1+1
         CALL DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &         M1,M2,NM1,ALPHN(K),BETAN(K),EE2(1,IAD),EE2(1,IAD+NM1),
     &         LDE1,IP2((K-1)*NM1+1),IER,IJOB)
         IF (IER.NE.0) GOTO 78
      END DO
      NDEC=NDEC+1
  30  CONTINUE
      IF (VARIAB.AND.IKEEP.EQ.1) THEN
         ICHAN=ICHAN+1
         IKEEP=0
         IF (ICHAN.GE.10.AND.NS.LT.NSMAX) GOTO 20
      END IF
      NSTEP=NSTEP+1
      IF (NSTEP.GT.NMAX) GOTO 178
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND) GOTO 177
          IF (INDEX2) THEN
             DO I=NIND1+1,NIND1+NIND2
                SCAL(I)=SCAL(I)/HHFAC
             END DO
          END IF
          IF (INDEX3) THEN
             DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3
                SCAL(I)=SCAL(I)/(HHFAC*HHFAC)
             END DO
          END IF
      XPH=X+H
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      IF (NS.EQ.3) THEN
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN.OR.CHANGE) THEN
         DO I=1,NNS
            ZZ(I)=0.D0
            FF(I)=0.D0
         END DO
      ELSE
         HQUOT=H/HOLD
         C3Q=HQUOT
         C1Q=C31*C3Q
         C2Q=C32*C3Q
         DO I=1,N
            AK1=CONT(I+N)
            AK2=CONT(I+N2)
            AK3=CONT(I+N3)
            Z1I=C1Q*(AK1+(C1Q-C32M1)*(AK2+(C1Q-C31M1)*AK3))
            Z2I=C2Q*(AK1+(C2Q-C32M1)*(AK2+(C2Q-C31M1)*AK3))
            Z3I=C3Q*(AK1+(C3Q-C32M1)*(AK2+(C3Q-C31M1)*AK3))
            ZZ(I)=Z1I
            ZZ(I+N)=Z2I
            ZZ(I+N2)=Z3I
            FF(I)=TI311*Z1I+TI312*Z2I+TI313*Z3I
            FF(I+N)=TI321*Z1I+TI322*Z2I+TI323*Z3I
            FF(I+N2)=TI331*Z1I+TI332*Z2I+TI333*Z3I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
      NEWT=0
      NIT=NIT1
      EXPMI=1.0D0/EXPMNS
      FNEWT=MAX(10*UROUND/RTOL1,MIN(0.03D0,RTOL1**(EXPMI-1.0D0)))
      FACCON=MAX(FACCON,UROUND)**0.8D0
      THETA=ABS(THET)
  40  CONTINUE
      IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
      DO I=1,N
         CONT(I)=Y(I)+ZZ(I)
      END DO
      CALL FCN(N,X+C31*H,CONT,ZZ,RPAR,IPAR)
      DO I=1,N
         CONT(I)=Y(I)+ZZ(I+N)
      END DO
      CALL FCN(N,X+C32*H,CONT,ZZ(1+N),RPAR,IPAR)
      DO I=1,N
         CONT(I)=Y(I)+ZZ(I+N2)
      END DO
      CALL FCN(N,XPH,CONT,ZZ(1+N2),RPAR,IPAR)
      NFCN=NFCN+3
C ---     SOLVE THE LINEAR SYSTEMS
      DO I=1,N
         A1=ZZ(I)
         A2=ZZ(I+N)
         A3=ZZ(I+N2)
         ZZ(I)=TI311*A1+TI312*A2+TI313*A3
         ZZ(I+N)=TI321*A1+TI322*A2+TI323*A3
         ZZ(I+N2)=TI331*A1+TI332*A2+TI333*A3
      END DO
      CALL SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &      M1,M2,NM1,FAC1,ALPHN(1),BETAN(1),E1,EE2,EE2(1,1+NM1),LDE1,
     &          ZZ,ZZ(1+N),ZZ(1+N2),FF,FF(1+N),FF(1+N2),
     &          CONT,IP1,IP2,IPHES,IER,IJOB)
      NSOL=NSOL+1
      NEWT=NEWT+1
      DYNO=0.D0
      DO I=1,N
         DENOM=SCAL(I)
         DYNO=DYNO+(ZZ(I)/DENOM)**2+(ZZ(I+N)/DENOM)**2
     &          +(ZZ(I+N2)/DENOM)**2
      END DO
      DYNO=DSQRT(DYNO/NNS)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
      IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
         THQ=DYNO/DYNOLD
         IF (NEWT.EQ.2) THEN
            THETA=THQ
         ELSE
            THETA=SQRT(THQ*THQOLD)
         END IF
         THQOLD=THQ
         IF (THETA.LT.0.99D0) THEN
            FACCON=THETA/(1.0D0-THETA)
            DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
            IF (DYTH.GE.1.0D0) THEN
               QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
               HHFAC=0.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
               H=HHFAC*H
               REJECT=.TRUE.
               LAST=.FALSE.
               IF (HHFAC.LE.0.5D0) UNEXN=.TRUE.
               IF (CALJAC) GOTO 20
               GOTO 10
            END IF
         ELSE
            GOTO 78
         END IF
      END IF
      DYNOLD=MAX(DYNO,UROUND)
      DO I=1,N
         IN=I+N
         IN2=I+N2
         F1I=FF(I)+ZZ(I)
         F2I=FF(IN)+ZZ(IN)
         F3I=FF(IN2)+ZZ(IN2)
         FF(I)=F1I
         FF(IN)=F2I
         FF(IN2)=F3I
         ZZ(I)=T311*F1I+T312*F2I+T313*F3I
         ZZ(IN)=T321*F1I+T322*F2I+T323*F3I
         ZZ(IN2)=T331*F1I+    F2I
      END DO
      IF (FACCON*DYNO.GT.FNEWT) GOTO 40
C --- ERROR ESTIMATION
      CALL ESTRAD (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD1,DD2,DD3,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,
     &          E1,LDE1,ZZ,ZZ(1+N),ZZ(1+N2),CONT,FF,FF(1+N),IP1,
     &          IPHES,SCAL,ERR,FIRST,REJECT,FAC1,RPAR,IPAR)
C       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT
      IF (ERR.LT.1.D0) THEN
         DO I=1,N
            Y(I)=Y(I)+ZZ(I+N2)
            Z2I=ZZ(I+N)
            Z1I=ZZ(I)
            CONT(I+N)=(Z2I-ZZ(I+N2))/C32M1
            AK=(Z1I-Z2I)/C31MC2
            ACONT3=Z1I/C31
            ACONT3=(AK-ACONT3)/C32
            CONT(I+N2)=(AK-CONT(I+N))/C31M1
            CONT(I+N3)=CONT(I+N2)-ACONT3
         END DO
      END IF
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      ELSE
      IF (NS.EQ.5) THEN
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN.OR.CHANGE) THEN
         DO I=1,NNS
            ZZ(I)=0.D0
            FF(I)=0.D0
         END DO
      ELSE
         HQUOT=H/HOLD
         DO K=1,NS
            CCQ=C(K)*HQUOT
            DO I=1,N
               VAL=CONT(I+NS*N)
               DO L=NS-1,1,-1
                  VAL=CONT(I+L*N)+(CCQ-C(NS-L)+1.D0)*VAL
               END DO
               ZZ(I+(K-1)*N)=CCQ*VAL
            END DO
         END DO
         DO I=1,N
            Z1I=ZZ(I)
            Z2I=ZZ(I+N)
            Z3I=ZZ(I+N2)
            Z4I=ZZ(I+N3)
            Z5I=ZZ(I+N4)
            FF(I)=TI511*Z1I+TI512*Z2I+TI513*Z3I+TI514*Z4I+TI515*Z5I
            FF(I+N)=TI521*Z1I+TI522*Z2I+TI523*Z3I+TI524*Z4I+TI525*Z5I
            FF(I+N2)=TI531*Z1I+TI532*Z2I+TI533*Z3I+TI534*Z4I+TI535*Z5I
            FF(I+N3)=TI541*Z1I+TI542*Z2I+TI543*Z3I+TI544*Z4I+TI545*Z5I
            FF(I+N4)=TI551*Z1I+TI552*Z2I+TI553*Z3I+TI554*Z4I+TI555*Z5I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
      NEWT=0
      NIT=NIT1+5
      EXPMI=1.0D0/EXPMNS
      FNEWT=MAX(10*UROUND/RTOL1,MIN(0.03D0,RTOL1**(EXPMI-1.0D0)))
      FACCON=MAX(FACCON,UROUND)**0.8D0
      THETA=ABS(THET)
 140  CONTINUE
      IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
      DO K=0,NS-1
         IADD=K*N
         DO I=1,N
            CONT(I)=Y(I)+ZZ(IADD+I)
         END DO
         CALL FCN(N,X+C(K+1)*H,CONT,ZZ(IADD+1),RPAR,IPAR)
      END DO
      NFCN=NFCN+NS
C ---     SOLVE THE LINEAR SYSTEMS
      DO I=1,N
         Z1I=ZZ(I)
         Z2I=ZZ(I+N)
         Z3I=ZZ(I+N2)
         Z4I=ZZ(I+N3)
         Z5I=ZZ(I+N4)
         ZZ(I)=TI511*Z1I+TI512*Z2I+TI513*Z3I+TI514*Z4I+TI515*Z5I
         ZZ(I+N)=TI521*Z1I+TI522*Z2I+TI523*Z3I+TI524*Z4I+TI525*Z5I
         ZZ(I+N2)=TI531*Z1I+TI532*Z2I+TI533*Z3I+TI534*Z4I+TI535*Z5I
         ZZ(I+N3)=TI541*Z1I+TI542*Z2I+TI543*Z3I+TI544*Z4I+TI545*Z5I
         ZZ(I+N4)=TI551*Z1I+TI552*Z2I+TI553*Z3I+TI554*Z4I+TI555*Z5I
      END DO
      CALL SLVRAR(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E1,LDE1,ZZ,FF,IP1,IPHES,IER,IJOB)
      DO K=1,2
         IAD=(K-1)*2*NM1+1
         CALL SLVRAI(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &     M1,M2,NM1,ALPHN(K),BETAN(K),EE2(1,IAD),EE2(1,IAD+NM1),LDE1,
     &     ZZ(1+(2*K-1)*N),ZZ(1+2*K*N),FF(1+(2*K-1)*N),
     &     FF(1+2*K*N),CONT,IP2((K-1)*NM1+1),IPHES,IER,IJOB)
      END DO
      NSOL=NSOL+1
      NEWT=NEWT+1
      DYNO=0.D0
      DO I=1,N
         DENOM=SCAL(I)
         DO K=0,NS-1
            DYNO=DYNO+(ZZ(I+K*N)/DENOM)**2
         END DO
      END DO
      DYNO=DSQRT(DYNO/NNS)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
      IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
         THQ=DYNO/DYNOLD
         IF (NEWT.EQ.2) THEN
            THETA=THQ
         ELSE
            THETA=SQRT(THQ*THQOLD)
         END IF
         THQOLD=THQ
         IF (THETA.LT.0.99D0) THEN
            FACCON=THETA/(1.0D0-THETA)
            DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
            IF (DYTH.GE.1.0D0) THEN
               QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
               HHFAC=0.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
               H=HHFAC*H
               REJECT=.TRUE.
               LAST=.FALSE.
               IF (HHFAC.LE.0.5D0) UNEXN=.TRUE.
               IF (CALJAC) GOTO 20
               GOTO 10
            END IF
         ELSE
            GOTO 78
         END IF
      END IF
      DYNOLD=MAX(DYNO,UROUND)
      DO I=1,N
         Z1I=FF(I)+ZZ(I)
         Z2I=FF(I+N)+ZZ(I+N)
         Z3I=FF(I+N2)+ZZ(I+N2)
         Z4I=FF(I+N3)+ZZ(I+N3)
         Z5I=FF(I+N4)+ZZ(I+N4)
         FF(I)=Z1I
         FF(I+N)=Z2I
         FF(I+N2)=Z3I
         FF(I+N3)=Z4I
         FF(I+N4)=Z5I
         ZZ(I)=T511*Z1I+T512*Z2I+T513*Z3I+T514*Z4I+T515*Z5I
         ZZ(I+N)=T521*Z1I+T522*Z2I+T523*Z3I+T524*Z4I+T525*Z5I
         ZZ(I+N2)=T531*Z1I+T532*Z2I+T533*Z3I+T534*Z4I+T535*Z5I
         ZZ(I+N3)=T541*Z1I+T542*Z2I+T543*Z3I+T544*Z4I+T545*Z5I
         ZZ(I+N4)=T551*Z1I+     Z2I         +     Z4I
      END DO
      IF (FACCON*DYNO.GT.FNEWT) GOTO 140
C --- ERROR ESTIMATION
      CALL ESTRAV (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,NS,NNS,E1,LDE1,
     &          ZZ,CONT,FF,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
C       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT
      IF (ERR.LT.1.D0) THEN
         DO I=1,N
            Y(I)=Y(I)+ZZ(I+N4)
            CONT(I+N5)=ZZ(I)/C(1)
         END DO
         DO K=1,NS-1
            FACT=1.D0/(C(NS-K)-C(NS-K+1))
            DO I=1,N
               CONT(I+K*N)=(ZZ(I+(NS-K-1)*N)-ZZ(I+(NS-K)*N))*FACT
            END DO
         END DO
         DO J=2,NS
            DO K=NS,J,-1
               FACT=1.D0/(C(NS-K)-C(NS-K+J))
               DO I=1,N
                  CONT(I+K*N)=(CONT(I+K*N)-CONT(I+(K-1)*N))*FACT
               END DO
            END DO
         END DO
      END IF
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      ELSE
      IF (NS.EQ.7) THEN
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN.OR.CHANGE) THEN
         DO I=1,NNS
            ZZ(I)=0.D0
            FF(I)=0.D0
         END DO
      ELSE
         HQUOT=H/HOLD
         DO K=1,NS
            CCQ=C(K)*HQUOT
            DO I=1,N
               VAL=CONT(I+NS*N)
               DO L=NS-1,1,-1
                  VAL=CONT(I+L*N)+(CCQ-C(NS-L)+1.D0)*VAL
               END DO
               ZZ(I+(K-1)*N)=CCQ*VAL
            END DO
         END DO
         DO I=1,N
            Z1I=ZZ(I)
            Z2I=ZZ(I+N)
            Z3I=ZZ(I+N2)
            Z4I=ZZ(I+N3)
            Z5I=ZZ(I+N4)
            Z6I=ZZ(I+N5)
            Z7I=ZZ(I+N6)
            FF(I)=TI711*Z1I+TI712*Z2I+TI713*Z3I+TI714*Z4I+TI715*Z5I
     *               +TI716*Z6I+TI717*Z7I
            FF(I+N)=TI721*Z1I+TI722*Z2I+TI723*Z3I+TI724*Z4I+TI725*Z5I
     *               +TI726*Z6I+TI727*Z7I
            FF(I+N2)=TI731*Z1I+TI732*Z2I+TI733*Z3I+TI734*Z4I+TI735*Z5I
     *               +TI736*Z6I+TI737*Z7I
            FF(I+N3)=TI741*Z1I+TI742*Z2I+TI743*Z3I+TI744*Z4I+TI745*Z5I
     *               +TI746*Z6I+TI747*Z7I
            FF(I+N4)=TI751*Z1I+TI752*Z2I+TI753*Z3I+TI754*Z4I+TI755*Z5I
     *               +TI756*Z6I+TI757*Z7I
            FF(I+N5)=TI761*Z1I+TI762*Z2I+TI763*Z3I+TI764*Z4I+TI765*Z5I
     *               +TI766*Z6I+TI767*Z7I
            FF(I+N6)=TI771*Z1I+TI772*Z2I+TI773*Z3I+TI774*Z4I+TI775*Z5I
     *               +TI776*Z6I+TI777*Z7I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
      NEWT=0
      NIT=NIT1+10
      EXPMI=1.0D0/EXPMNS
      FNEWT=MAX(10*UROUND/RTOL1,MIN(0.03D0,RTOL1**(EXPMI-1.0D0)))
      FACCON=MAX(FACCON,UROUND)**0.8D0
      THETA=ABS(THET)
 240  CONTINUE
      IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
      DO K=0,NS-1
         IADD=K*N
         DO I=1,N
            CONT(I)=Y(I)+ZZ(IADD+I)
         END DO
         CALL FCN(N,X+C(K+1)*H,CONT,ZZ(IADD+1),RPAR,IPAR)
      END DO
      NFCN=NFCN+NS
C ---     SOLVE THE LINEAR SYSTEMS
      DO I=1,N
         Z1I=ZZ(I)
         Z2I=ZZ(I+N)
         Z3I=ZZ(I+N2)
         Z4I=ZZ(I+N3)
         Z5I=ZZ(I+N4)
         Z6I=ZZ(I+N5)
         Z7I=ZZ(I+N6)
         ZZ(I)=TI711*Z1I+TI712*Z2I+TI713*Z3I+TI714*Z4I+TI715*Z5I
     *            +TI716*Z6I+TI717*Z7I
         ZZ(I+N)=TI721*Z1I+TI722*Z2I+TI723*Z3I+TI724*Z4I+TI725*Z5I
     *            +TI726*Z6I+TI727*Z7I
         ZZ(I+N2)=TI731*Z1I+TI732*Z2I+TI733*Z3I+TI734*Z4I+TI735*Z5I
     *            +TI736*Z6I+TI737*Z7I
         ZZ(I+N3)=TI741*Z1I+TI742*Z2I+TI743*Z3I+TI744*Z4I+TI745*Z5I
     *            +TI746*Z6I+TI747*Z7I
         ZZ(I+N4)=TI751*Z1I+TI752*Z2I+TI753*Z3I+TI754*Z4I+TI755*Z5I
     *            +TI756*Z6I+TI757*Z7I
         ZZ(I+N5)=TI761*Z1I+TI762*Z2I+TI763*Z3I+TI764*Z4I+TI765*Z5I
     *            +TI766*Z6I+TI767*Z7I
         ZZ(I+N6)=TI771*Z1I+TI772*Z2I+TI773*Z3I+TI774*Z4I+TI775*Z5I
     *            +TI776*Z6I+TI777*Z7I
      END DO
      CALL SLVRAR(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &        M1,M2,NM1,FAC1,E1,LDE1,ZZ,FF,IP1,IPHES,IER,IJOB)
      DO K=1,3
         IAD=(K-1)*2*NM1+1
         CALL SLVRAI(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &      M1,M2,NM1,ALPHN(K),BETAN(K),EE2(1,IAD),EE2(1,IAD+NM1),LDE1,
     &      ZZ(1+(2*K-1)*N),ZZ(1+2*K*N),FF(1+(2*K-1)*N),
     &      FF(1+2*K*N),CONT,IP2((K-1)*NM1+1),IPHES,IER,IJOB)
      END DO
      NSOL=NSOL+1
      NEWT=NEWT+1
      DYNO=0.D0
      DO I=1,N
         DENOM=SCAL(I)
         DO K=0,NS-1
            DYNO=DYNO+(ZZ(I+K*N)/DENOM)**2
         END DO
      END DO
      DYNO=DSQRT(DYNO/NNS)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
      IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
         THQ=DYNO/DYNOLD
         IF (NEWT.EQ.2) THEN
            THETA=THQ
         ELSE
            THETA=SQRT(THQ*THQOLD)
         END IF
         THQOLD=THQ
         IF (THETA.LT.0.99D0) THEN
            FACCON=THETA/(1.0D0-THETA)
            DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
            IF (DYTH.GE.1.0D0) THEN
               QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
               HHFAC=0.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
               H=HHFAC*H
               REJECT=.TRUE.
               LAST=.FALSE.
               IF (HHFAC.LE.0.5D0) UNEXN=.TRUE.
               IF (CALJAC) GOTO 20
               GOTO 10
            END IF
         ELSE
            GOTO 78
         END IF
      END IF
      DYNOLD=MAX(DYNO,UROUND)
      DO I=1,N
         Z1I=FF(I)+ZZ(I)
         Z2I=FF(I+N)+ZZ(I+N)
         Z3I=FF(I+N2)+ZZ(I+N2)
         Z4I=FF(I+N3)+ZZ(I+N3)
         Z5I=FF(I+N4)+ZZ(I+N4)
         Z6I=FF(I+N5)+ZZ(I+N5)
         Z7I=FF(I+N6)+ZZ(I+N6)
         FF(I)=Z1I
         FF(I+N)=Z2I
         FF(I+N2)=Z3I
         FF(I+N3)=Z4I
         FF(I+N4)=Z5I
         FF(I+N5)=Z6I
         FF(I+N6)=Z7I
         ZZ(I)=T711*Z1I+T712*Z2I+T713*Z3I+T714*Z4I+T715*Z5I
     *            +T716*Z6I+T717*Z7I
         ZZ(I+N)=T721*Z1I+T722*Z2I+T723*Z3I+T724*Z4I+T725*Z5I
     *            +T726*Z6I+T727*Z7I
         ZZ(I+N2)=T731*Z1I+T732*Z2I+T733*Z3I+T734*Z4I+T735*Z5I
     *            +T736*Z6I+T737*Z7I
         ZZ(I+N3)=T741*Z1I+T742*Z2I+T743*Z3I+T744*Z4I+T745*Z5I
     *            +T746*Z6I+T747*Z7I
         ZZ(I+N4)=T751*Z1I+T752*Z2I+T753*Z3I+T754*Z4I+T755*Z5I
     *            +T756*Z6I+T757*Z7I
         ZZ(I+N5)=T761*Z1I+T762*Z2I+T763*Z3I+T764*Z4I+T765*Z5I
     *            +T766*Z6I+T767*Z7I
         ZZ(I+N6)=T771*Z1I+Z2I+Z4I+Z6I
      END DO
      IF (FACCON*DYNO.GT.FNEWT) GOTO 240
C --- ERROR ESTIMATION
      CALL ESTRAV (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,NS,NNS,E1,LDE1,
     &          ZZ,CONT,FF,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
C       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT
      IF (ERR.LT.1.D0) THEN
         DO I=1,N
            Y(I)=Y(I)+ZZ(I+(NS-1)*N)
            CONT(I+NS*N)=ZZ(I)/C(1)
         END DO
         DO K=1,NS-1
            FACT=1.D0/(C(NS-K)-C(NS-K+1))
            DO I=1,N
               CONT(I+K*N)=(ZZ(I+(NS-K-1)*N)-ZZ(I+(NS-K)*N))*FACT
            END DO
         END DO
         DO J=2,NS
            DO K=NS,J,-1
               FACT=1.D0/(C(NS-K)-C(NS-K+J))
               DO I=1,N
                  CONT(I+K*N)=(CONT(I+K*N)-CONT(I+(K-1)*N))*FACT
               END DO
            END DO
         END DO
      END IF
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      ELSE
CASE       (NS.EQ.1)
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      IF (FIRST.OR.STARTN.OR.CHANGE) THEN
         DO I=1,NS
            ZZ(I)=0.D0
            FF(I)=0.D0
         END DO
      ELSE
         HQUOT=H/HOLD
         DO I=1,N
            Z1I=HQUOT*CONT(I+N)
            ZZ(I)=Z1I
            FF(I)=Z1I
         END DO
      END IF
C *** *** *** *** *** *** ***
C  LOOP FOR THE SIMPLIFIED NEWTON ITERATION
C *** *** *** *** *** *** ***
      NEWT=0
      NIT=NIT1-3
      EXPMI=1.0D0/EXPMNS
      FNEWT=MAX(10*UROUND/RTOL1,0.03D0)
      FACCON=MAX(FACCON,UROUND)**0.8D0
      THETA=ABS(THET)
 440  CONTINUE
      IF (NEWT.GE.NIT) GOTO 78
C ---     COMPUTE THE RIGHT-HAND SIDE
      DO I=1,N
         CONT(I)=Y(I)+ZZ(I)
      END DO
      CALL FCN(N,XPH,CONT,ZZ,RPAR,IPAR)
      NFCN=NFCN+1
C ---     SOLVE THE LINEAR SYSTEMS
      CALL SLVRAR(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E1,LDE1,ZZ,FF,IP1,IPHES,IER,IJOB)
      NSOL=NSOL+1
      NEWT=NEWT+1
      DYNO=0.D0
      DO I=1,N
         DENOM=SCAL(I)
         DYNO=DYNO+(ZZ(I)/DENOM)**2
      END DO
      DYNO=DSQRT(DYNO/NNS)
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE
      IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN
         THQ=DYNO/DYNOLD
         IF (NEWT.EQ.2) THEN
            THETA=THQ
         ELSE
            THETA=SQRT(THQ*THQOLD)
         END IF
         THQOLD=THQ
         IF (THETA.LT.0.99D0) THEN
            FACCON=THETA/(1.0D0-THETA)
            DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT
            IF (DYTH.GE.1.0D0) THEN
               QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH))
               HHFAC=0.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT))
               H=HHFAC*H
               REJECT=.TRUE.
               LAST=.FALSE.
               IF (HHFAC.LE.0.5D0) UNEXN=.TRUE.
               IF (CALJAC) GOTO 20
               GOTO 10
            END IF
         ELSE
            GOTO 78
         END IF
      END IF
      DYNOLD=MAX(DYNO,UROUND)
      DO I=1,N
         F1I=FF(I)+ZZ(I)
         FF(I)=F1I
         ZZ(I)=F1I
      END DO
      IF (FACCON*DYNO.GT.FNEWT) GOTO 440
C --- ERROR ESTIMATION
      CALL ESTRAV (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,NS,NNS,E1,LDE1,
     &          ZZ,CONT,FF,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
C       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT
      IF (ERR.LT.1.D0) THEN
         DO I=1,N
            Y(I)=Y(I)+ZZ(I)
            CONT(I+N)=ZZ(I)
         END DO
      END IF
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
      END IF
      END IF
      END IF
C *** *** *** *** *** *** ***
C *** *** *** *** *** *** ***
C --- COMPUTATION OF HNEW
C --- WE REQUIRE .2<=HNEW/H<=8.
      FAC=MIN(SAFE,(1+2*NIT)*SAFE/(NEWT+2*NIT))
      QUOT=MAX(FACR,MIN(FACL,ERR**EXPO/FAC))
      HNEW=H/QUOT
C *** *** *** *** *** *** ***
C  IS THE ERROR SMALL ENOUGH ?
C *** *** *** *** *** *** ***
      IF (ERR.LT.1.D0) THEN
C --- STEP IS ACCEPTED
         FIRST=.FALSE.
         NACCPT=NACCPT+1
         IF (PRED.AND..NOT.CHANGE) THEN
C       --- PREDICTIVE CONTROLLER OF GUSTAFSSON
            IF (NACCPT.GT.1) THEN
               FACGUS=(HACC/H)*(ERR**2/ERRACC)**EXPO/SAFE
               FACGUS=MAX(FACR,MIN(FACL,FACGUS))
               QUOT=MAX(QUOT,FACGUS)
               HNEW=H/QUOT
            END IF
            HACC=H
            ERRACC=MAX(1.0D-2,ERR)
         END IF
         XOLD=X
         HOLD=H
         X=XPH
C       --- UPDATE SCALING
         IF (ITOL.EQ.0) THEN
            DO I=1,N
               SCAL(I)=ATOL1+RTOL1*ABS(Y(I))
            END DO
         ELSE
            DO I=1,N
               QUOTT=ATOL(I)/RTOL(I)
               RTOL1=0.1D0*RTOL(I)**EXPMNS
               ATOL1=RTOL1*QUOTT
               SCAL(I)=ATOL1+RTOL1*ABS(Y(I))
            END DO
         END IF
         IF (IOUT.NE.0) THEN
             NRSOL=NACCPT+1
             XSOL=X
             XOSOL=XOLD
             DO I=1,N
                CONT(I)=Y(I)
             END DO
             NSOLU=N
             HSOL=HOLD
             CALL SOLOUT(NRSOL,XOSOL,XSOL,Y,CONT,LRC,NSOLU,
     &                   RPAR,IPAR,IRTRN)
             IF (IRTRN.LT.0) GOTO 179
         END IF
         CALJAC=.FALSE.
         IF (LAST) THEN
            H=HOPT
            IDID=1
            RETURN
         END IF
         CALL FCN(N,X,Y,Y0,RPAR,IPAR)
         NFCN=NFCN+1
         HNEW=POSNEG*MIN(ABS(HNEW),HMAXN)
         HOPT=HNEW
         HOPT=MIN(H,HNEW)
         IF (REJECT) HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
         REJECT=.FALSE.
         IF ((X+HNEW/QUOT1-XEND)*POSNEG.GE.0.D0) THEN
            H=XEND-X
            LAST=.TRUE.
         ELSE
            QT=HNEW/H
            HHFAC=H
            IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) THEN
               IKEEP=1
               GOTO 30
            END IF
            H=HNEW
         END IF
         HHFAC=H
         IF (THETA.LE.THET) GOTO 20
         GOTO 10
      ELSE
C --- STEP IS REJECTED
         REJECT=.TRUE.
         LAST=.FALSE.
         IF (FIRST) THEN
             H=H*0.1D0
             HHFAC=0.1D0
         ELSE
             HHFAC=HNEW/H
             H=HNEW
         END IF
         IF (NACCPT.GE.1) NREJCT=NREJCT+1
         IF (CALJAC) GOTO 20
         GOTO 10
      END IF
C --- UNEXPECTED STEP-REJECTION
  78  CONTINUE
      UNEXP=.TRUE.
      IF (IER.NE.0) THEN
          NSING=NSING+1
          IF (NSING.GE.5) GOTO 176
      END IF
      H=H*0.5D0
      HHFAC=0.5D0
      REJECT=.TRUE.
      LAST=.FALSE.
      IF (CALJAC) GOTO 20
      GOTO 10
C --- FAIL EXIT
 176  CONTINUE
      WRITE(6,979)X
      WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER
      IDID=-4
      RETURN
 177  CONTINUE
      WRITE(6,979)X
      WRITE(6,*) ' STEP SIZE T0O SMALL, H=',H
      IDID=-3
      RETURN
 178  CONTINUE
      WRITE(6,979)X
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'
      IDID=-2
      RETURN
C --- EXIT CAUSED BY SOLOUT
 179  CONTINUE
      WRITE(6,979)X
 979  FORMAT(' EXIT OF RADAU AT X=',E18.4)
      IDID=2
      RETURN
      END
C
C     END OF SUBROUTINE RADCOV

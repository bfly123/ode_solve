      SUBROUTINE COERCV(NS,C,DD,U1,ALPH,BETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(0:NS),DD(NS),ALPH(NS),BETA(NS)
      C(0)=0.D0
      C(NS)=1.D0
      GOTO (1,11,3,11,5,11,7) NS
 11   CONTINUE
      RETURN
  1   CONTINUE
      C(1)=1.0D0
      U1=1.0D0
      DD(1)=-1.0D0
      RETURN
  3   CONTINUE
      SQ6=DSQRT(6.D0)
      C(1)=(4.D0-SQ6)/10.D0
      C(2)=(4.D0+SQ6)/10.D0
      ST9=9.D0**(1.D0/3.D0)
      U1=(6.D0+ST9*(ST9-1))/30.D0
      ALP=(12.D0-ST9*(ST9-1))/60.D0
      BET=ST9*(ST9+1)*DSQRT(3.D0)/60.D0
      CNO=ALP**2+BET**2
      U1=1.0D0/U1
      ALPH(1)=ALP/CNO
      BETA(1)=BET/CNO
      RETURN
  5   CONTINUE
      C(1)=  0.5710419611451768219312D-01
      C(2)=  0.2768430136381238276800D+00
      C(3)=  0.5835904323689168200567D+00
      C(4)=  0.8602401356562194478479D+00
      DD(1)= -0.2778093394406463730479D+02
      DD(2)=  0.3641478498049213152712D+01
      DD(3)= -0.1252547721169118720491D+01
      DD(4)=  0.5920031671845428725662D+00
      DD(5)= -0.2000000000000000000000D+00
      U1=  0.6286704751729276645173D+01
      ALPH(1)=  0.3655694325463572258243D+01
      BETA(1)=  0.6543736899360077294021D+01
      ALPH(2)=  0.5700953298671789419170D+01
      BETA(2)=  0.3210265600308549888425D+01
      RETURN
  7   CONTINUE
      C(1)=  0.2931642715978489197205D-01
      C(2)=  0.1480785996684842918500D+00
      C(3)=  0.3369846902811542990971D+00
      C(4)=  0.5586715187715501320814D+00
      C(5)=  0.7692338620300545009169D+00
      C(6)=  0.9269456713197411148519D+00
      DD(1)= -0.5437443689412861451458D+02
      DD(2)=  0.7000024004259186512041D+01
      DD(3)= -0.2355661091987557192256D+01
      DD(4)=  0.1132289066106134386384D+01
      DD(5)= -0.6468913267673587118673D+00
      DD(6)=  0.3875333853753523774248D+00
      DD(7)= -0.1428571428571428571429D+00
      U1=  0.8936832788405216337302D+01
      ALPH(1)=  0.4378693561506806002523D+01
      BETA(1)=  0.1016969328379501162732D+02
      ALPH(2)=  0.7141055219187640105775D+01
      BETA(2)=  0.6623045922639275970621D+01
      ALPH(3)=  0.8511834825102945723051D+01
      BETA(3)=  0.3281013624325058830036D+01
      RETURN
      END
C
C     END OF SUBROUTINE COERCV
C
C ***********************************************************


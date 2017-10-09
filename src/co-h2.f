c   pleas use it as such following example: 
c       IMPLICIT REAL*8 (A-H,O-Z)
c       OPEN(6,FILE='co-h2-pes-vmlrq-g.chk')
c       Do 300 XPHI=0.0D0, 90.0D0, 30.D0
c       Do 300 XTH1=0.0D0, 180.0D0, 30.0D0
c       DO 300 XTH2=0.0D0, 180.0D0, 30.0D0
c       DO 300   RX=3.0D0, 10.0D0, 0.2
c         CALL  COH2PES(XPHI,XTH1,XTH2,RX,V,0) 
c         WRITE(6,663)XPHI,XTH1,XTH2,RX,V
c  300  CONTINUE
c  663  FORMAT(1x,3f7.1,f8.2,4f20.4)
c       END

c******************************************************************************
      SUBROUTINE COH2PES(XPHI,XTH1,XTH2,XR,V,iv3)
c******************************************************************************
c** Subroutine to generate values of the vibrationally averaged 4D-MLR analyic
c  potential energy surfaces for complexes formed between H2 and CO in vibrational
c  level v= 0 or 1, as determined by Hui Li,Xiao-Long Zhang, Robert J. Le Roy and
c  Pierre-Nicholas Roy [JCP, submitted 2013]. On first call, input values of 
c  select appropriate 4D-MLR expansion parameters from DATA statements in subroutine 
c  PARAREAD;  subsequent calls generate additional potential function values for 
c  that same case.
c** Input variables:
c  XR  - distance between CO and H2 centre of mass in [Angst], pointing from
c  the center of mass of CO to the center of mass of H2.
c  XTH1  - Jacobi angular coordinate 'theta1' in degrees, which is the angle between 
c the vector XR pointing from the center of mass of CO to the center of mass of H2 and
c  the vector pointing from O atom to C. 
c  XTH2  - Jacobi angular coordinate 'theta2' in degrees, which is the angle between 
c the vector XR pointing from the center of mass of CO to the center of mass of H2 and
c  the vector pointing from H_2 atom to H_1. 
c  XPHI  - Jacobi dihedral angular coordinate 'phi' in degrees,dihedral angle between 
c  the two half planes extending from the vector R to C and the H1 atom.
c** Output:   V [cm-1]  is the calculated interaction energy '\Delta{V}'.
c-----------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER MXPARM,MMAX,id1,id2 
      PARAMETER (MXPARM=600,MMAX=16)
      PARAMETER (id1=2, id2=1)
      INTEGER iv3,NPARM,NP,NQ,NDE,NRE,NCN,MCM,NS,NL,NPOW
      INTEGER NC6L1,NC6L2,NC6LMAX,NC7L1,NC7L2,NC7LMAX
      INTEGER NC8L1,NC8L2,NC8LMAX,NQM1,NQM2
      INTEGER NBETA(0:50),NBETAL1(0:50),NBETAL2(0:50),NBETALMAX(0:50)
      REAL*8  BETA(0:50),Pn1(0:50,0:50), Pn2(0:50,0:50)
      REAL*8  PV(MXPARM),PD(MXPARM)
      REAL*8  C6(MXPARM),C7(MXPARM),C8(MXPARM)
      REAL*8  el(20)
      REAL*8  YC,Re,De,Vasy,RREF,AREF,AREFp,AREFq,Rep,CN,VLRe,
     2 dVLRedRe,phiINF,RTPp,RTPq,yp,yq,ype,dype,yPOW,XP,SUM,DSUM,VLR,
     3 XPW,DER,XDE,XRE,YTH1,YTH2,YPHI,CTH1,CTH2,STH1,STH2,CPHI,
     4 C5Sum,C6Sum,C8Sum,Qa,Qb,bDAMP,T0,DM(MMAX),DMP(MMAX),DMPP(MMAX)
      common/mulprod/q(MMAX,MMAX)
      COMMON /DATABLK/PV,PI,bDAMP,C6,C7,C8,RREF,CN,NBETAL1,NBETAL2,
     1 NBETALMAX,NDEL1,NDEL2,NDELMAX,NREL1,NREL2,NRELMAX,NCN,MCM,
     2 NP,NQ,NS,NL,NC6L1,NC6L2,NC6LMAX,NC7L1,NC7L2,NC7LMAX,
     3 NC8L1,NC8L2,NC8LMAX,NQM1,NQM2
c-----------------------------------------------------------------------
      DATA IPAR/0/
      SAVE IPAR
c-----------------------------------------------------------------------------
      IF(IPAR.EQ.0) THEN
          CALL PARAREAD(iv3)
           IPAR= 1
           ENDIF
        PI=DACOS(-1.0D0)
        RY=XR
        YTH1=(180.-XTH2)*PI/180.D0
        YTH2=(180.-XTH1)*PI/180.D0
        YPHI=XPHI*PI/180.D0
        CTH1=DCOS(YTH1)
        CTH2=DCOS(YTH2)
        STH1=DSIN(YTH1)
        STH2=DSIN(YTH2)
        CPHI=DCOS(YPHI)
        CALL plmrb(Pn1,CTH1,2)
        CALL plmrb(Pn2,CTH2,2)
c caculate the derivative of the parameters of De not including
c the coefficient before, only the three angle function A(TH1,TH2,PHI) 
c with associated Legendre expansion.
c     De expantion 
      IP=0
      De=0.0d0
      DO 201 L1=0,NDEL1,id1
      DO 201 L2=0,NDEL2,id2
       LMIN=IABS(L1-L2)
c        MM=MIN0(L1,L2)
       LLMAX=MIN0(NDELMAX,L1+L2)
        DO 201 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IP=IP+1
          PD(IP)=al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)
          De=De+PD(IP)*PV(IP)
         ENDIF
 201   CONTINUE
c caculate the derivative of the parameters of Re not including
c the coefficient before, only the three angle function A(TH1,TH2,PHI) 
c with associated Legendre expansion.
c     Re expantion 
      NDE=IP
      Re=0.0d0
      DO 202 L1=0,NREL1,id1
      DO 202 L2=0,NREL2,id2
       LMIN=IABS(L1-L2)
c        MM=MIN0(L1,L2)
       LLMAX=MIN0(NRELMAX,L1+L2)
        DO 202 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IP=IP+1
          PD(IP)=al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)
          Re=Re+PD(IP)*PV(IP)
         ENDIF
  202  CONTINUE
       NRE=IP-NDE
c C6,C7 and C8 coefficient which included the induce and dispersion parts reported 
c by Jankowski and Szalewicz J. Chem. Phys. 123 (2005) 104301
c and c J. Chem. Phys. 108 (1998) 3554  
       AREF= RREF*Re
c       IF(RREF.LE.0.d0) AREF= Re
       AREFp= AREF**NP
       AREFq= AREF**NQ
       Rep= Re**NP

      IK=0
      C6Sum=0.0d0
      DO 331 L1=0,NC6L1,2
      DO 331 L2=0,NC6L2,2
       LMIN=IABS(L1-L2)
c       MM=MIN0(L1,L2)
       LLMAX=MIN0(NC6LMAX,L1+L2)
        DO 331 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IK=IK+1
          C6Sum=C6Sum+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*C6(IK)
         ENDIF
  331  CONTINUE

      IK=0
      C7Sum=0.0d0
      DO 332 L1=0,NC7L1,2
      DO 332 L2=1,NC7L2,2
       LMIN=IABS(L1-L2)
c       MM=MIN0(L1,L2)
       LLMAX=MIN0(NC7LMAX,L1+L2)
        DO 332 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IK=IK+1
          C7Sum=C7Sum+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*C7(IK)
         ENDIF
  332  CONTINUE

      IK=0
      C8Sum=0.0d0
      DO 333 L1=0,NC8L1,2
      DO 333 L2=0,NC8L2,2
       LMIN=IABS(L1-L2)
c       MM=MIN0(L1,L2)
       LLMAX=MIN0(NC8LMAX,L1+L2)
        DO 333 L=LMIN,LLMAX
          LTOT=L1+L2+L
         IF(MOD(LTOT,2).EQ.0) THEN
          IK=IK+1
          C8Sum=C8Sum+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*C8(IK)
         ENDIF
  333  CONTINUE


c  V_elect accounts for the electrostatic interaction due to the 
c  molecular permanent dipole, quadrupoles etc. For CO-H2 system
c  main electrostatic interactions which reported by Jankowski and
c  Szalewicz J. Chem. Phys. 108, 3554 (1998) are the following 
c  products of multipole moments q(CO,H2)=q(L1,L2) (L1,L2=1,2
c  are relative the dipole, quadrupole moments of CO or H2) are:    

      IF(NQM1.GT.0.d0) CALL elst(YPHI,YTH1,YTH2,el)

      C4Sum=-el(4)
      C5Sum=-el(5)
      C6Sum=C6Sum-el(6)
      C7Sum=C7Sum-el(7)
      C8Sum=C8Sum-el(8)
c** For normal inverse-power sum MLR case, with or without damping
c  calculate damping coeffient for H2-CO as below
c Ip(H2)=15.45 ev W. Kolos and J. Rychlewski JCP 98,3960(1993) 
c Ip(CO)=14.01 ev=113031.3/8065.5447 by P. Erman et. al. Chem. Phys. Lett. 215, 173(1993) 
c Ip(H)=13.60 ev 
c Pd(H2)=[Ip(H2)/Ip(H)]^(2/3)=1.210833319 
c Pd(CO)=[Ip(CO)/Ip(H)]^(2/3)=1.019998386  
c bd(H2,CO)=2.78*[2*Pd(H2)*Pd(CO2)/(Pd(H2)+Pd(CO2))]=3.078164541  
       IF(bDAMP.GT.0.d0) CALL DAMPIG(Re,bDAMP,MMAX,DM,DMP,DMPP)
         T0= C6Sum/Re**NCN
         IF(bDAMP.GT.0.d0) T0= T0*DM(NCN)
         VLRe= T0
         dVLRedRe= -NCN*T0/Re
         IF(bDAMP.GT.0.d0)
     &    dVLRedRe= dVLRedRe + T0*DMP(NCN)/DM(NCN)
          
       IF(MCM.GT.NCN) THEN
        MMN= MCM - NCN
       IF(NP.LE.MMN) MMN= 0
       IF(MMN.GT.0) THEN
C  with and without damping function 
       IF(bDAMP.GT.0.d0) THEN
         VLRe= C6Sum*DM(6)/Re**6+C7Sum*DM(7)/Re**7+C8Sum*DM(8)/Re**8+
     &   C4Sum*DM(4)/Re**4+C5Sum*DM(5)/Re**5                  ! electrostatic C4,C5     
         dVLRedRe=dVLRedRe-7.0D0*C7Sum*DM(7)/Re**8+C7Sum*DMP(7)/Re**7
     1   -8.0D0*C8Sum*DM(8)/Re**9+C8Sum*DMP(8)/Re**8
     2  -4.0D0*C4Sum*DM(4)/Re**5+C4Sum*DMP(4)/Re**4
     3  -5.0D0*C5Sum*DM(5)/Re**6+C5Sum*DMP(5)/Re**5
         ELSE
         VLRe= C6Sum/Re**6+C7Sum/Re**7+C8Sum/Re**8+
     &   C4Sum/Re**4+C5Sum/Re**5                              ! electrostatic C4,C5     
         dVLRedRe=dVLRedRe-7.0D0*C7Sum/Re**8
     1   -8.0D0*C8Sum/Re**9-4.0D0*C4Sum/Re**5
     2   -5.0D0*C5Sum/Re**6
       ENDIF
       ENDIF
         phiINF= DLOG(2.d0*De/VLRe)
       ENDIF
       RTPp= RY**NP
       RTPq= RY**NQ
       yp= (RTPp - AREFp)/(RTPp + AREFp)
       yq= (RTPq - AREFq)/(RTPq + AREFq)
       ype= (RTPp - Rep)/(RTPp + Rep)
c caculate the derivative of the parameters of BETA(N) not including
c the coefficient before, only the three angle function A(TH1,TH2,PHI)
c with associated Legendre expansion.
       NPOW= NS
       IF(RY.GE.Re) NPOW= NL
        yPOW= 1.d0 - yp
        SUM=0.0
        DSUM=0.0
        IP=NDE+NRE
        DO 204 J=0,NPOW
           BETA(J)=0.0
           NPS=0
         DO 203 L1=0,NBETAL1(J),id1
         DO 203 L2=0,NBETAL2(J),id2
            LMIN=IABS(L1-L2)
c           MM=MIN0(L1,L2)
            LLMAX=MIN0(NBETALMAX(J),L1+L2)
         DO 203 L=LMIN,LLMAX
            LTOT=L1+L2+L
           IF(MOD(LTOT,2).EQ.0) THEN
            IP=IP+1
            NPS=NPS+1
            PD(IP)=al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*yq**(J)
            BETA(J)=BETA(J)+al1l2l0(L1,L2,L,YTH1,YTH2,YPHI)*PV(IP)
           ENDIF
  203  CONTINUE
         NBETA(J)=NPS
        SUM=SUM+BETA(J)*yq**(J)
c        IF(RREF.LE.0.D0) DSUM=DSUM+yPOW*BETA(J)*(J)*yp**(J-1)
        DSUM=DSUM+yPOW*BETA(J)*(J)*yq**(J-1)
  204  CONTINUE
c  caculate the derivative of the parameters of Vasy 
        PD(IP+1)=1.0D0
        Vasy=PD(IP+1)*PV(IP+1)
        XP= SUM*yPOW+ phiINF*yp

       IF(bDAMP.GT.0.d0) CALL DAMPIG(RY,bDAMP,MMAX,DM,DMP,DMPP)
c with and without damping function 
        IF(bDAMP.GT.0.d0)THEN
          VLR= C6Sum*DM(6)/RY**6+
     &    C4Sum*DM(4)/RY**4+C5Sum*DM(5)/RY**5        ! electrostatic C4,C5     
         ELSE
          VLR= C6Sum/RY**6+
     &    C4Sum/RY**4+C5Sum/RY**5                    ! electrostatic C4,C5     
         ENDIF
        IF(MMN.GT.0) THEN
          IF(bDAMP.GT.0.d0)THEN
            VLR= C6Sum*DM(6)/RY**6+
     1      C7Sum*DM(7)/RY**7+C8Sum*DM(8)/RY**8+
     2      C4Sum*DM(4)/RY**4+C5Sum*DM(5)/RY**5      ! electrostatic C4,C5     
          ELSE
            VLR= C6Sum/RY**6+
     1      C7Sum/RY**7+C8Sum/RY**8+
     2      C4Sum/RY**4+C5Sum/RY**5                  ! electrostatic C4,C5     
          ENDIF
       ENDIF
         XPW= DEXP(-XP*ype) * VLR/VLRe
         YC= De*(1.d0 - XPW)**2-De+Vasy
         V=YC
      RETURN
      END 
c-----------------------------------------------------------------------------
      SUBROUTINE PARAREAD(iv3)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MXPARM=600,MMAX=16,MXDE=96, MXRE=78, MXPHI=107,ISP=2)
      PARAMETER (MXC6=6,MXC7=7,MXC8=18,MQQ=32,isC=12,isO=16)
      INTEGER I,J,K,L1,L2
      INTEGER isv,iv3,IP,NDEL1,NDEL2,NDELMAX,NREL1,NREL2,NRELMAX
      INTEGER IFXP(MXPARM),NPARM,NP,NQ,NDE,NRE,NCN,MCM,NS,NL,NPOW
      INTEGER NC6L1,NC6L2,NC6LMAX,NC7L1,NC7L2,NC7LMAX
      INTEGER NC8L1,NC8L2,NC8LMAX,NQM1,NQM2
      INTEGER NBETA(0:50),NBETAL1(0:50),NBETAL2(0:50),NBETALMAX(0:50)
      INTEGER NB(0:4) 
      REAL*8  BETA(0:50),Pn1(0:50,0:50), Pn2(0:50,0:50),CCN(2)
      REAL*8  CC6(ISP,MXPARM),CC7(ISP,MXPARM),CC8(ISP,MXPARM)
      REAL*8  PV(MXPARM),C6(MXPARM),C7(MXPARM),C8(MXPARM)
      REAL*8  qq(ISP,MXPARM),el(20)
      REAL*8  DE(ISP,MXDE),RE(ISP,MXRE),PHIPHI(ISP,MXPHI)
      REAL*8  YC,RREF,bDAMP
      common/mulprod/q(MMAX,MMAX)
      COMMON /DATABLK/PV,PI,bDAMP,C6,C7,C8,RREF,CN,NBETAL1,NBETAL2,
     1 NBETALMAX,NDEL1,NDEL2,NDELMAX,NREL1,NREL2,NRELMAX,NCN,MCM,
     2 NP,NQ,NS,NL,NC6L1,NC6L2,NC6LMAX,NC7L1,NC7L2,NC7LMAX,
     3 NC8L1,NC8L2,NC8LMAX,NQM1,NQM2
c-----------------------------------------------------------------------
      CALL fct(40)
      CALL fill3j(13,13,13)
      DATA NDEL1/4/,NDEL2/12/,NDELMAX/12/
      DATA NREL1/4/,NREL2/10/,NRELMAX/10/
      DATA NP/5/,NQ/3/,NS/4/,NL/4/,RREF/1.2d0/,bDAMP/-3.08d0/
      DATA NPOW/4/
      DATA NB(0)/0/, NBETAL1(0)/6/, NBETAL2(0)/6/,NBETALMAX(0)/6/
      DATA NB(1)/1/, NBETAL1(1)/4/, NBETAL2(1)/4/,NBETALMAX(1)/4/
      DATA NB(2)/2/, NBETAL1(2)/2/, NBETAL2(2)/3/,NBETALMAX(2)/3/
      DATA NB(3)/3/, NBETAL1(3)/2/, NBETAL2(3)/2/,NBETALMAX(3)/2/
      DATA NB(4)/4/, NBETAL1(4)/2/, NBETAL2(4)/2/,NBETALMAX(4)/2/
      DATA NCN/6/, MCM/8/, CCN(1)/1.0d0/,CCN(2)/1.0d0/
      DATA NC6L1/2/,NC6L2/2/,NC6LMAX/4/
      DATA NC7L1/2/,NC7L2/3/,NC7LMAX/5/
      DATA NC8L1/4/,NC8L2/4/,NC8LMAX/6/
      DATA NQM1/8/,NQM2/8/
      DATA rbohr/0.5291772108d0/
      const=2.0d0*109737.31568525d0
      PI=DACOS(-1.0D0)
c    The products of multipole moments for electrostatic potential 
      DATA (qq(1,I),I=1,MQQ)/-0.2507499744D-01, -0.7055255778D+00,
     &-0.1673033366D+01, -0.3832278506D+01, -0.5488563844D+01,
     &-0.9625119280D+01, -0.9721102605D+01, -0.1321158135D+02,
     &-0.1727366451D-01, -0.4860224678D+00, -0.1152519244D+01,
     &-0.2639980059D+01, -0.4177607776D+01, -0.7327648561D+01,
     &-0.7329639285D+01, -0.9922143024D+01, -0.1090503880D-01,
     &-0.3068308909D+00, -0.7275970348D+00, -0.1666646065D+01,
     &-0.5775551813D+01, -0.1011914913D+02, -0.1065525948D+02,
     &-0.1472172089D+02, 0.3523532824D-02, 0.7257092676D+00,
     & 0.1864523280D+01, 0.4705015475D+01, 0.6431841208D+01,
     & 0.1126839311D+02, 0.1189450980D+02, 0.1644935614D+02/

c     read C6 expantion coefficients      
      DATA (CC6(1,I),I=1, MXC6)/0.3114216572D+02,0.8885711004D+01,
     & 0.8345203277D+01, 0.2427542512D+00, 0.6487880548D+00,
     & 0.7006910915D+01/ 
c     read C7 expantion coefficients      
      DATA (CC7(1,I),I=1, MXC7)/0.1609819683D+03,-0.7652324520D+01,
     &-0.6608207936D+01, 0.3296753267D+02, -0.1118541861D+00,
     &-0.2609930193D+00, -0.4372022481D+01/
c     read C8 expantion coefficients      
      DATA (CC8(1,I),I=1, MXC8)/0.8274331003D+03,0.1858830915D+04,
     &-0.5593655430D+02, 0.5060561936D+03, 0.1980847038D+02,
     &-0.5843802331D+02, 0.4425611116D+03, -0.7175184304D+00,
     &-0.1615538984D+01, -0.3925865625D+02, 0.9703792283D+01,
     & 0.1113875669D+00, -0.9339367353D-01, 0.7501673195D+01,
     & 0.0000000000D+00, 0.0000000000D+00, 0.0000000000D+00,
     & 0.0000000000D+00/

      DATA (DE(1,I),I=1,MXDE)/
     & 4.64050D+01,-4.43600D+00,-1.80500D+01,-3.16000D+00, 4.14000D+00,
     & 5.71000D+00, 0.00000D+00,-2.68000D+00,-5.20000D-01, 7.90000D-01,
     & 3.60000D-01,-1.40000D-01,-1.10000D-01, 2.79000D+00,-1.34000D+00,
     &-7.99000D+00, 1.24000D+00, 1.05600D+01, 1.81740D+02, 8.60000D-01,
     & 6.34000D+00, 7.37700D+01, 1.20000D-01, 7.09000D+00,-4.83000D+00,
     &-2.10000D-01, 7.50000D-01,-1.87400D+01, 2.60000D-01,-1.67000D+00,
     &-3.32000D+00, 2.10000D-01,-7.00000D-01, 5.40000D+00,-6.00000D-02,
     & 2.00000D-01, 2.80000D+00,-2.00000D-01, 3.00000D-01,-1.00000D+00,
     &-1.00000D-01, 0.00000D+00,-1.00000D+00, 1.00000D-01, 0.00000D+00,
     & 3.00000D-01, 0.00000D+00, 1.09000D+00, 5.00000D-02, 1.04000D+00,
     & 1.10000D-01, 1.97000D+00, 5.28000D+00, 0.00000D+00, 1.20000D-01,
     & 1.11000D+00, 3.51000D+00, 0.00000D+00, 1.30000D-01, 7.00000D-01,
     & 2.20000D+00, 1.66100D+01, 0.00000D+00, 1.00000D-01, 5.00000D-01,
     & 1.70000D+00, 1.40000D+01, 4.00000D-02, 0.00000D+00, 1.00000D-01,
     & 1.00000D+00, 4.20000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     & 3.00000D-01,-2.30000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     &-3.00000D-01,-2.20000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     &-2.00000D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     & 0.00000D+00, 0.00000D+00, 0.00000D+00,-1.00000D-01, 0.00000D+00,
     & 0.00000D+00/


      DATA (RE(1,I),I=1,MXRE)/
     & 3.76883D+00, 4.62130D-01, 9.15300D-01, 8.44000D-02,-7.36000D-02,
     &-7.71000D-02,-1.06000D-02, 3.01000D-02, 1.17000D-02,-6.00000D-03,
     &-3.40000D-03, 1.00800D-01, 2.49000D-02, 2.38000D-02,-8.30000D-03,
     &-1.33000D-02,-9.46300D-01, 4.10000D-03,-3.24000D-02,-5.63600D-01,
     & 8.20000D-03,-3.33000D-02,-1.15500D-01, 3.00000D-03,-6.00000D-03,
     & 1.00100D-01,-5.00000D-03, 1.00000D-02, 6.99000D-02,-4.00000D-03,
     & 5.00000D-03,-1.10000D-02, 0.00000D+00, 0.00000D+00,-2.20000D-02,
     & 1.00000D-03, 0.00000D+00, 1.00000D-03,-2.00000D-03, 7.00000D-04,
     &-2.10000D-03, 5.50000D-03, 3.10000D-03, 0.00000D+00, 3.26000D-02,
     & 0.00000D+00, 0.00000D+00,-4.00000D-03,-5.00000D-03, 5.00000D-04,
     & 0.00000D+00, 0.00000D+00, 1.30000D-02, 1.18400D-01, 1.00000D-03,
     & 0.00000D+00, 0.00000D+00, 1.30000D-02, 9.90000D-02, 0.00000D+00,
     & 0.00000D+00, 0.00000D+00, 7.00000D-03, 2.60000D-02,-1.20000D-03,
     & 0.00000D+00,-1.00000D-03, 2.00000D-03,-2.30000D-03, 0.00000D+00,
     & 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 2.00000D-03,
     & 0.00000D+00, 0.00000D+00, 0.00000D+00/

      DATA (PHIPHI(1,I),I=1,MXPHI)/
     &-5.7700D-02, 9.8700D-02, 3.4720D-01, 1.7400D-01, 5.9000D-02,
     &-9.4000D-02,-3.1000D-02, 0.0000D+00, 1.5000D-02, 1.4800D-01,
     &-6.0000D-03,-2.3000D-02, 2.6600D-01,-6.0000D-03,-9.0000D-03,
     &-2.7200D-01,-1.7000D-02,-2.3000D-02,-2.6400D-01, 1.5000D-02,
     &-2.1000D-02, 2.0000D-02,-1.6000D-02, 0.0000D+00, 0.0000D+00,
     &-1.9000D-02, 8.0000D-03, 0.0000D+00,-7.7000D-02, 3.0000D-03,
     &-1.3000D-02, 0.0000D+00,-2.0000D-03, 1.4000D-02,-1.0000D-02,
     & 0.0000D+00,-4.0000D-03, 4.0000D-03, 0.0000D+00, 2.0000D-03,
     &-4.0000D-03, 0.0000D+00,-3.0000D-03,-2.0000D-03, 7.0000D-03,
     & 0.0000D+00, 0.0000D+00, 3.0000D-03,-2.0000D-03, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 1.6850D-01,-1.3400D-01, 1.7900D-01,-3.4000D-02,-7.0000D-02,
     & 1.0700D-01, 5.2000D-02, 4.4200D-01, 0.0000D+00, 6.7000D-02,
     & 1.3000D-01, 0.0000D+00,-1.0000D-02, 0.0000D+00, 0.0000D+00,
     &-1.6000D-02, 0.0000D+00, 3.3000D-02, 0.0000D+00, 0.0000D+00,
     &-2.0000D-02, 0.0000D+00, 2.0000D-02,-3.000D-02,
     &-6.0000D-02,-1.5000D-01, 1.3000D-01, 3.0000D-02, 8.0000D-02,
     & 0.0000D+00, 5.1000D-01, 0.0000D+00, 6.0000D-02, 0.0000D+00,
     & 0.0000D+00,
     & 1.5000D-01, 1.0000D-01, 2.9000D-01, 0.0000D+00,-1.5000D-01,
     &-2.0000D-02, 0.0000D+00,
     &-6.7000D-01, 2.3000D-01, 4.0000D-01, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00/



c    isC=12, isO=16, v3=1

c    The products of multipole moments for electrostatic potential 
      DATA (qq(2,I),I=1,MQQ)/-0.2013547700D-01, -0.6912770104D+00,
     & -0.1676480968D+01, -0.3829557636D+01, -0.5488563844D+01,
     & -0.9625119280D+01, -0.9721102605D+01, -0.1321158135D+02,
     & -0.1387092762D-01, -0.4762069146D+00, -0.1154894228D+01,
     & -0.2638105706D+01, -0.4177607776D+01, -0.7327648561D+01,
     & -0.7329639285D+01, -0.9922143024D+01, -0.8756856647D-02,
     & -0.3006342330D+00, -0.7290963861D+00, -0.1665462767D+01,
     & -0.5775551813D+01, -0.1011914913D+02, -0.1065525948D+02,
     & -0.1472172089D+02, 0.3523532824D-02, 0.7257092676D+00,
     &  0.1864523280D+01, 0.4705015475D+01, 0.6431841208D+01,
     &  0.1126839311D+02, 0.1189450980D+02, 0.1644935614D+02/
c     read C6 expantion coefficients      
      DATA (CC6(2,I),I=1, MXC6)/0.3138057401D+02, 0.8953735407D+01,
     & 0.8409089833D+01, 0.2446126521D+00, 0.6537548403D+00,
     & 0.7060552197D+01/
c     read C7 expantion coefficients      
      DATA (CC7(2,I),I=1, MXC7)/0.1622143629D+03, -0.7710906755D+01,
     & -0.6658796954D+01, 0.3321991504D+02, -0.1127104839D+00,
     & -0.2629910467D+00, -0.4405492421D+01/
c     read C8 expantion coefficients      
      DATA (CC8(2,I),I=1, MXC8)/0.8337674997D+03, 0.1873061162D+04,
     &-0.5636477559D+02, 0.5099302978D+03, 0.1996011377D+02,
     &-0.5888539456D+02, 0.4459491303D+03, -0.7230113800D+00,
     &-0.1627906714D+01, -0.3955920019D+02, 0.9778079492D+01,
     & 0.1122402925D+00, -0.9410864714D-01, 0.7559102121D+01,
     & 0.0000000000D+00, 0.0000000000D+00, 0.0000000000D+00,
     & 0.0000000000D+00/   

      DATA (DE(2,I),I=1,MXDE)/
     & 4.66490D+01,-4.61500D+00,-1.80000D+01,-3.52000D+00, 3.88000D+00,
     & 5.86000D+00, 1.60000D-01,-2.71000D+00,-6.00000D-01, 7.60000D-01,
     & 3.90000D-01,-1.40000D-01,-1.30000D-01, 3.18000D+00,-1.57000D+00,
     &-1.21400D+01, 1.20000D+00, 1.04800D+01, 1.79620D+02, 8.30000D-01,
     & 6.12000D+00, 7.52100D+01, 9.00000D-02, 7.06000D+00, -4.26000D+00,
     &-2.20000D-01, 8.50000D-01, -1.89800D+01, 2.50000D-01,-1.60000D+00,
     &-3.68000D+00, 2.10000D-01, -7.00000D-01, 5.27000D+00,-6.00000D-02,
     & 1.00000D-01, 2.90000D+00, -1.40000D-01, 3.00000D-01,-9.00000D-01,
     &-1.00000D-01, 1.00000D-01, -1.00000D+00, 1.00000D-01,-1.00000D-01,
     & 3.00000D-01, -1.00000D-01, 1.06000D+00, 5.00000D-02, 9.50000D-01,
     & 1.00000D-01, 1.95000D+00, 5.33000D+00, 0.00000D+00, 1.50000D-01,
     & 1.03000D+00, 2.80000D+00, 0.00000D+00, 8.00000D-02, 7.60000D-01,
     & 2.10000D+00, 1.56500D+01, 0.00000D+00, 1.40000D-01, 5.00000D-01,
     & 1.70000D+00, 1.39400D+01, 0.00000D+00, 0.00000D+00, 1.00000D-01,
     & 1.00000D+00, 4.50000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     & 3.00000D-01, -2.10000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     &-2.00000D-01, -2.30000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     &-2.00000D-01, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     &-1.00000D-01, 0.00000D+00, 0.00000D+00, -1.00000D-01, 0.00000D+00,
     & 0.00000D+00/
     

      DATA (RE(2,I),I=1,MXRE)/
     & 3.76885D+00, 4.60710D-01, 9.17300D-01, 8.91000D-02, -7.06000D-02,
     &-7.75000D-02, -1.18000D-02, 2.99000D-02, 1.21000D-02,-5.70000D-03,
     &-3.20000D-03, 9.94000D-02, 2.60000D-02, 4.17000D-02, -8.40000D-03,
     &-1.38000D-02, -9.27400D-01, 3.80000D-03,-3.17000D-02,-5.68000D-01,
     & 6.90000D-03, -3.33000D-02,-1.19200D-01, 2.50000D-03,-6.00000D-03,
     & 9.86000D-02, -4.60000D-03, 1.00000D-02, 7.07000D-02,-3.80000D-03,
     & 5.00000D-03, -1.00000D-02, 0.00000D+00, 0.00000D+00,-2.20000D-02,
     & 1.00000D-03, 0.00000D+00, 1.00000D-03, -1.00000D-03, 0.00000D+00,
     &-2.20000D-03, 4.10000D-03, 2.70000D-03, -3.00000D-03, 3.16000D-02,
     & 0.00000D+00, 0.00000D+00,-5.00000D-03,-8.50000D-03, 6.00000D-04,
     & 0.00000D+00, 0.00000D+00, 1.30000D-02, 1.13700D-01, 1.10000D-03,
     & 0.00000D+00, 0.00000D+00, 1.30000D-02, 9.90000D-02, 0.00000D+00,
     & 0.00000D+00, 0.00000D+00, 7.00000D-03, 2.70000D-02, -1.20000D-03,
     & 0.00000D+00, 0.00000D+00, 2.00000D-03, -2.00000D-03, 0.00000D+00,
     & 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     & 0.00000D+00, 0.00000D+00, 0.00000D+00/


      DATA (PHIPHI(2,I),I=1,MXPHI)/
     &-5.7500D-02, 1.0060D-01, 3.5000D-01, 1.7900D-01, 6.2000D-02,
     &-9.5000D-02, -3.1000D-02, 0.0000D+00, 1.5000D-02, 1.3700D-01,
     &-4.0000D-03, -2.7000D-02, 2.4800D-01, -6.0000D-03, -5.0000D-03,
     &-2.6800D-01, -2.7000D-02, -2.3000D-02, -2.6000D-01, 1.5000D-02,
     &-2.0000D-02, 1.8000D-02, -1.5000D-02, 0.0000D+00, 0.0000D+00,
     &-1.6000D-02, 6.0000D-03, 0.0000D+00, -7.1000D-02, 6.0000D-03,
     &-4.0000D-03, 0.0000D+00, -2.0000D-03, 4.0000D-03, 0.0000D+00,
     & 0.0000D+00, -5.0000D-03, 3.0000D-03, 0.0000D+00, 0.0000D+00,
     &-4.0000D-03, 0.0000D+00, -3.0000D-03, -2.0000D-03, 6.0000D-03,
     & 2.0000D-03, 0.0000D+00, 0.0000D+00, -2.0000D-03, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00,
     & 1.7040D-01, -1.2500D-01, 1.7300D-01, -2.0000D-02, -6.3000D-02,
     & 1.0700D-01, 5.1000D-02, 4.3800D-01, 0.0000D+00, 6.0000D-02,
     & 1.0000D-01, 0.0000D+00, 0.0000D+00, -3.0000D-02, 0.0000D+00,
     &-1.7000D-02, 0.0000D+00, 2.7000D-02, 0.0000D+00, 1.0000D-02,
     & 0.0000D+00, 0.0000D+00, 0.0000D+00, -5.0000D-03,
     &-7.0000D-02, -1.3000D-01, 1.7000D-01, 5.0000D-02, 8.0000D-02,
     & 0.0000D+00, 5.0000D-01, 0.0000D+00, 1.8000D-01, 0.0000D+00,
     & 0.0000D+00,
     & 1.4000D-01, 9.0000D-02, 2.9000D-01, 0.0000D+00, -1.5000D-01,
     & 0.0000D+00, 2.0000D-01,
     &-6.6000D-01, 1.9000D-01, 3.0000D-01, 0.0000D+00, 0.0000D+00,
     & 0.0000D+00, 0.0000D+00/


        isv= 0
       IF((isC.EQ.12).and.(isO.EQ.16).and.(iv3.EQ.0)) isv= 1
       IF((isC.EQ.12).and.(isO.EQ.16).and.(iv3.EQ.1)) isv= 2

c    read the products of multipole moments (in atomic units) for electrostatic energy
      K=0
      DO  L1=2, NQM1, 2
      DO  L2=1, NQM2, 1 
         K=K+1
         q(L1,L2)=qq(isv,K)
      ENDDO
      ENDDO

      CN=CCN(isv) 

c    read C6 expantion coefficients
      DO I=1,MXC6
          C6(I)=CC6(isv,I)*rbohr**6*const
          ENDDO
c    read C7 expantion coefficients
      DO I=1,MXC7
          C7(I)=CC7(isv,I)*rbohr**7*const
          ENDDO
c    read C8 expantion coefficients
      DO I=1,MXC8
          C8(I)=CC8(isv,I)*rbohr**8*const
          ENDDO

      DO I=1,MXDE
c... first prepare well depth expansion parameters
          PV(I)=DE(isv,I)
          ENDDO
      IP=MXDE
      DO I=1,MXRE
c... next prepare potential minimum position expansion parameters
          IP=IP+1
          PV(IP)=RE(isv,I)
          ENDDO
      IP=MXDE+MXRE
c... then, prepare exponent coefficient \beta_i expansion parameters
      DO I=1,MXPHI
          IP=IP+1
          PV(IP)=PHIPHI(isv,I)
          ENDDO
c** Finally ... (if desired) read shift for definition of energy zero
      IP=MXDE+MXRE+MXPHI
      NPARM=IP+1
      PV(NPARM)=0.0D0
      RETURN
  600 FORMAT(/' Generate 4D-MLRR potential for {',i2,'}C{',i2,'}O(v3=',
     1  i1,')-H2')
  602 FORMAT(/' *** 4D Potential H2-CO Function selection parameters  i
     1sC=',i2,'   isO=',i3,'   v3=',i2/'   do not match defined cases.'
     2  '  Generate potential for {12}C{16}O(v=0)-H2 instead.')
      END
c***********************************************************************
c    
c This subroutine fills up the matrix w3j for C6 with values of 3-j coefficient
c     
      subroutine fill3j(l1max,l2max,lmax)
      implicit real*8 (a-h,o-z)
      common /w3jcg/w3j(0:50,0:30,0:30,0:60)
      dimension x(60)

      do l1=0,l1max
       do l2=0,l2max
        lmin=iabs(l1-l2)
        mm=min0(l1,l2)
        llmax=min0(lmax,l1+l2)

        do l=lmin,llmax
         do m=0,mm
          m1=m
          m2=-m
          mmm=0
          call cgc(l1,m1,l2,m2,l,mmm,c,1)
          w3j(m,l1,l2,l)=c
         end do
        end do
       end do
      end do

      return
      end
C ----------------------------------------------------------------------------
c
c Calculate the function Al1l2L for a given set of angles...
c It is assumed that the th1 and th2 angles are between the monomer bond
c and the "inner" part of the intermolecular axis.
c
       function al1l2l0(l1,l2,l,th1,th2,phi)
       implicit real*8 (a-h,o-z)
       dimension p1(0:50,0:50), p2(0:50,0:50)
       common /w3jcg/w3j(0:50,0:30,0:30,0:60)
       data izer/0/, ione/1/, pifact/12.566370614359d0/
c
       c1 = dcos(th1)
       c2 = dcos(th2)

       call plmrb(p1,c1,l1)
       call plmrb(p2,c2,l2)
       mmax = min(l1,l2)
       sum = 0.d0
       do m=1,mmax
        value=w3j(m,l1,l2,l)
        sum = sum + (-1)**m*value*p1(l1,m)*p2(l2,m)*dcos(m*phi)
       end do
       value=w3j(0,l1,l2,l)
       sum = 2*sum + value*p1(l1,0)*p2(l2,0)
c
       al1l2l0 = sum*pifact*(-1.d0)**(l1+l2+l)/dsqrt((2.d0*l1+1.d0)*
     1           (2.d0*l2+1.d0))
c
       return
       end

C --------------------------------------------------------------------------
c
c Compute the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials....
c
        subroutine plmrb(p,x,lmax)
        implicit real*8 (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
c inverse of dsqrt(2Pi)
        data twopinv /0.3989422804014d0/
c
c starting value
c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
c
c compute the diagonal elements
c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
c
c compute P_lm along the columns with fixed m

c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do
        end do
c
c Renormalize values...
c
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do
c
        return
        end

C -------------------------------------------------------------------------
c
c compute the matrix of N!
c
        subroutine fct(nmax)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        f(0) = 1.d0
        do i=1,nmax
         f(i) = f(i-1)*i
        end do
        return
        end
C -------------------------------------------------------------------------
c
Calculate the Clebsh-Gordan coefficient (or the 3-j symbol)
c The parameter ind3j.eq.1 indicates that the 3-J symbol is returned
c
        subroutine cgc(j1,m1,j2,m2,j,m,value,ind3j)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        d3jfact = 1.d0
        if(ind3j.eq.1) then
         d3jfact = ((-1.d0)**(j1-j2-m))/dsqrt(dfloat(2*j+1))
         m = -m
        endif

c
c Check the triangle conditions
c
        if(j.gt.(j1+j2)) write(6,*)'triangle violated'
        if(j.lt.abs(j1-j2)) write(6,*)'triangle violated'
        if((m1+m2).ne.m) then
          value = 0.d0
          return
        endif


c Calculation proper... the pre-sum factor....
c
        facn = (2*j+1)*f(j1+j2-j)*f(j1-m1)*f(j2-m2)*f(j+m)*f(j-m)
        facd = f(j1+j2+j+1)*f(j+j1-j2)*f(j+j2-j1)*f(j1+m1)*f(j2+m2)
        fac = dsqrt(facn/facd)

c
c determine the limit of k summation...
c
        kmax = min(j2+j-m1,j-m,j1-m1)
        if(kmax.lt.0) kmax = 0
        kmin = max(-j1-m1,-j2+j-m1,0)

c
c perform the summation (at least one cycle must be completed...
c
        sum = 0.d0
        do k = kmin,kmax
         facn = f(j1+m1+k)*f(j2+j-m1-k)
         facd = f(k)*f(j-m-k)*f(j1-m1-k)*f(j2-j+m1+k)
         sum = sum + (facn/facd)*(-1)**k
        end do
        value = d3jfact*fac*sum*(-1)**(j1-m1)
       return
       end
C -------------------------------------------------------------------------
       SUBROUTINE flush(nunit)
       endfile nunit
       backspace nunit
       end
C -------------------------------------------------------------------------

c***********************************************************************
c      SUBROUTINE DAMPIG(r,b,MMAX,DM,DMP,DMPP)
c** Subroutine to generate values DM(m) and the first and second radial
c  derivatives DMP(m) and DMPP(m) of normalized incomplete gamma 
c  functions of orders m=0 to MMAX, at radial distance 'r', for damping
c  parameter 'b'.  NOTE that  DM(m)= {Tang-Toennies function}(m+1).
c***********************************************************************
c     INTEGER MMAX,I,m
c     REAL*8 b,r,br,XP,SSm,SSm1,SSm2,TK,DM(MMAX),DMP(MMAX),DMPP(MMAX)
c     br= b*r
c     XP= DEXP(-br)
c     SSm= 0.d0
c     SSm1= 0.d0
c     SSm2= 0.d0
c     TK= 1.d0
c     DO  m=1, MMAX
c         SSm2= SSm1
c         SSm1= SSm
c         SSm= SSm+ TK
c         DM(m)= 1.d0 - XP*SSM
c         DMP(m)= b*XP*(SSm - SSm1)
c         DMPP(m)= b**2 *XP*(2.d0*SSm1 - SSm - SSm2)
c         TK= TK*br/DFLOAT(m)
c         ENDDO
c     RETURN
c     END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DAMPIG(r,b,MMAX,DM,DMP,DMPP)
c** Subroutine to generate values DM and the first and second radial
c  derivatives DMP and DMPP of the Tang-Toennies/incomplete gamma 
c  function damping function of all orders up to NMAX at radial distance
c   'r' for damping parameter 'b'
c***********************************************************************
      INTEGER MMAX,m,n
      REAL*8 b,r,br,XP,SSm,SSm1,SSm2,mfact,TK,DM(MMAX),DMP(MMAX),
     1       DMPP(MMAX),SM(-1:MMAX)

      br= b*r
      XP= DEXP(-br)
      SSm= 0.d0
      SSm1= 0.d0
      SSm2= 0.d0
      TK= 1.d0
      DO m=1,MMAX
         SSm2= SSm1
         SSm1= SSm
         SSm= SSm + TK
         DM(m)= 1.d0 - XP*SSm
         DMP(m)= b*XP*(SSm - SSm1)
         DMPP(m)= b**2 *XP* (2.d0*SSm1 - SSm - SSm2)
         TK= TK*br/DFLOAT(m)
      ENDDO
      IF(DM(MMAX).LT.1.0D-13) THEN
         mfact= 1.d0
         DO n=1,MMAX
            mfact= mfact*DFLOAT(n)
         ENDDO
         SSm= 0.d0
         TK= (br)**MMAX/mfact
         DO n=MMAX,MMAX+3
            SSm= SSm+TK
            TK= TK*br/DFLOAT(n+1)
         ENDDO
         SM(MMAX)= SSm
         TK= (br)**MMAX/mfact
         DO n=1,MMAX-1
            TK= TK*DFLOAT(MMAX+1-n)/br
            SM(MMAX-n)= TK+SM(MMAX-n+1)
         ENDDO
         SM(0)=1.d0 + SM(1)
         SM(-1)=SM(0)
         DO m=1,MMAX
            DM(m)= XP*SM(m)
            DMP(m)= -b*XP*(SM(m)-SM(m-1))
            DMPP(m)= b**2 *XP*(SM(m)-2.d0*SM(m-1)+SM(m-2))
         ENDDO
      ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c Subroutine for the calculation of the multipole
c part of the electrostatic energy for a given
c dimer conformation. th1 and th2 expressed
c in radians. th2 is in fact equal to Pi-theta2,
c where theta2 is in the same coordinate system as theta1.
c
        subroutine elst(XPHI,TH1,TH2,el)
        implicit real*8 (a-h,o-z)
        parameter (MMAX=16)
        dimension el(20)
        common/factorial/ f(0:40)
        common/mulprod/ q(MMAX,MMAX)
        data rbohr/0.5291772108d0/
        const=2.0d0*109737.31568525d0
c
        els = 0.d0
        do i=1,20
         el(i) = 0.d0
        end do
c
        do L1=2,8,2
        do L2=1,8,1
           LL = L1 + L2
           mola = (-1)**L1
           glam = mola*al1l2l0(L1,L2,LL,TH1,TH2,XPHI)*q(L1,L2)
           term = dsqrt(f(2*LL+1)/(f(2*L1)*f(2*L2)))*glam*rbohr**(LL+1)
           el(LL+1) = el(LL+1) + term
         end do
        end do
c
        do i=1,20
         el(i) = const*el(i)       !output unit is cm-1*angstrom^i  
        end do
        return
        end
C --------------------------------------------------------------------------
c    




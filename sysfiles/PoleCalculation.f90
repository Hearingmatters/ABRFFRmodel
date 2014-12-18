!this subroutine will calculate the polepositions of the next iteration
!depending on the current V and/or Y or envelope (depends on the choice)
SUBROUTINE PoleCalculation
  USE Declare
  IMPLICIT NONE
  
  REAL(8), PARAMETER :: THdB=30d0 
  REAL(8), PARAMETER :: Ax=1d0 !slopes of function for Q=0.03 dBvel
  REAL(8), PARAMETER :: Bx=-167.6014d0  !offset of function for Q=0.03 dBvel 
  REAL(8) PoleE(1:n)
  REAL(8) Yknee1(1:n) 
  REAL(8) Vknee1(1:n) 
  REAL(8) Vknee2(1:n)
  REAL(8) Yknee2(1:n)
  REAL(8) Yknee1CST(1:n)
  REAL(8) Yknee2CST(1:n)
  REAL(8) Vvect(1:n)
  REAL(8) Yvect(1:n)
  REAL(8) Yknee1F(1:n)
  REAL(8) Yknee2F(1:n)
  REAL(8) Rth(1:n)
  REAL(8) NDY1(1:n)
  REAL(8) NDY2(1:n)
  REAL(8) NDV1(1:n)
  REAL(8) NDV2(1:n)
  REAL(8) RthY1(1:n)
  REAL(8) RthY2(1:n)
  REAL(8) RthV1(1:n)
  REAL(8) RthV2(1:n)
  REAL(8) dY1(1:n)
  REAL(8) dY2(1:n)
  REAL(8) dV1(1:n)
  REAL(8) dV2(1:n)
  REAL(8),PARAMETER ::  factor=100d0
  REAL(8) BoffsetV(1:n)
  REAL(8) BoffsetY(1:n)
  REAL(8) PoleS(1:n)
  REAL(8) Theta0(1:n)
  REAL(8) Theta(1:n)
  REAL(8) Sfoc(1:n)
  REAL(8) Sa(1:n)
  REAL(8) Sb(1:n)
  REAL(8) Se(1:n)
  REAL(8) Sxp(1:n)
  REAL(8) Syp(1:n)
  REAL(8) Sy(1:n)
  REAL(8) Rand(1:n)
  REAL(8) Yint(1:n)
  REAL(8) Vint(1:n)
  INTEGER ctr
  INTEGER fhtel
  INTEGER onektel
  INTEGER, dimension(2) :: seed

PoleE(1:n)=0.3d0 !passive model pole
Yknee1(1:n)=6.9183d-10 !(Ybm at 30dB) for scaling symmetric model
Vknee1(1:n)=4.3652d-6 !(Vbm at 30dB) for scaling symmetric model

!Use for PoleE always the 0.3 value (max passive)
!From the 30 dB nonlinearity threshold, find the intersection (vbm) between the curves of the Q=0.3 pole (linear and passive)
!and the line starting from 30 dB and the starting pole with the desired compression slope

IF (useExtStartingPoles) THEN
	Vknee1(1:n)=Vthresholds;
	!Yknee1(1:n)=Ythresholds;
ENDIF

BoffsetV(1:n)=-compressionslope*THdB+20d0*LOG10(Vknee1(1:n)) !find the offset of the compression curves
BoffsetY(1:n)=-compressionslope*THdB+20d0*LOG10(Yknee1(1:n))
!find the intersection of this curve with the Q=0.03 curve
Vint(1:n)=(Bx-BoffsetV(1:n))/(compressionslope-Ax) !is the intersection in dB on xaxis
Vknee2(1:n)=Ax*Vint(1:n)+Bx !what it corresponds to in dB on y axis
Vknee2(1:n)=10d0 ** (Vknee2(1:n)/20d0) !is the intersection in m/s

Yint(1:n)=(Bx-BoffsetY(1:n))/(compressionslope-Ax) !is the intersection in dB xaxis
Yknee2(1:n)=Ax*Yint(1:n)+Bx !what it corresponds to in dB on y axis
Yknee2(1:n)=10d0 ** (Yknee2(1:n)/20d0) !is the intersection in m

 
!Introduce variability in the Vknee thresholds, kneepoint varies slightly from section
!to section
IF (useZweigIrregularity) THEN
   seed(1)=Subjectnr * 29
   seed(2)=Subjectnr + 2010 * 08                
   CALL RANDOM_SEED(PUT=seed) 
   CALL RANDOM_NUMBER(Rth(1:n))
   
   !section to section variability in the kneepoint threshold in dB
   IF (useKneeVar) THEN
      dY1(1:n)=20d0*LOG10(Yknee1(1:n))
      dY2(1:n)=20d0*LOG10(Yknee2(1:n))
      dV1(1:n)=20d0*LOG10(Vknee1(1:n))
      dV2(1:n)=20d0*LOG10(Vknee2(1:n))
      NDY1(1:n)= dY1(1:n) + ((Rth(1:n) - 0.5d0)* KneeVar/0.5d0) !normal dist in [dB]
      NDY2(1:n)= dY2(1:n) + ((Rth(1:n) - 0.5d0)* KneeVar/0.5d0) !normal dist in [dB]
      NDV1(1:n)= dV1(1:n) + ((Rth(1:n) - 0.5d0)* KneeVar/0.5d0) !normal dist in [dB]
      NDV2(1:n)= dV2(1:n) + ((Rth(1:n) - 0.5d0)* KneeVar/0.5d0) !normal dist in [dB]
      RthY1(1:n)= 10d0 ** (NDY1(1:n)/20d0) !log dist around kneepoint in disp [m]
      RthY2(1:n)= 10d0 ** (NDY2(1:n)/20d0) !log dist around kneepoint in disp [m]
      RthV1(1:n)= 10d0 ** (NDV1(1:n)/20d0) !log dist around kneepoint in vel [m]
      RthV2(1:n)= 10d0 ** (NDV2(1:n)/20d0) !log dist around kneepoint in vel [m]
   ELSE
      RthY1(1:n)= Yknee1(1:n)
      RthY2(1:n)= Yknee2(1:n)
      RthV1(1:n)= Vknee1(1:n)
      RthV2(1:n)= Vknee2(1:n)
   ENDIF
   !the main variability sits in the starting pole variation  
    DO i=1,n
         Rand(i)= IrrPct * (Rth(i) - 0.5d0)
	IF (useExtStartingPoles) THEN
	 PoleS(i)=(1d0+Rand(i)) * StartingPoles(i)
	ELSE
	 PoleS(i)=(1d0+Rand(i)) * SheraPo
	ENDIF
    ENDDO 
    !PoleS(1:n)=(1d0+Rand) * SheraPo

ELSE !no irregularities
   IF (useExtStartingPoles) THEN
	PoleS(1:n)=StartingPoles(1:n)
   ELSE
   	PoleS(1:n)= SheraPo
   ENDIF
   RthY1(1:n)= Yknee1(1:n)
   RthY2(1:n)= Yknee2(1:n)
   RthV1(1:n)= Vknee1(1:n)
   RthV2(1:n)= Vknee2(1:n)
ENDIF

!!Don't put irregularities in the low frequencies
!!If subject dependent irregularity patterns are needed, they need to go here
    ctr=0
    onektel=0
    DO i=1,n       
      IF (omega(i)/(2d0 *pi) .GE. 100) THEN
         ctr=ctr+1;
      ENDIF
      IF (omega(i)/(2d0 *pi) .GE. 1000) THEN
         onektel=onektel+1;   
         !find the section freq
      ENDIF
   ENDDO
    
   IF(.NOT.useLFirregularity) THEN
      IF (useExtStartingPoles) THEN
	PoleS(ctr:n)=StartingPoles(ctr:n)
      ELSE 		
        PoleS(ctr:n)= SheraPo
      ENDIF
      RthY1(ctr:n)= Yknee1(ctr:n)
      RthY2(ctr:n)= Yknee2(ctr:n)
      RthV1(ctr:n)= Vknee1(ctr:n)
      RthV2(ctr:n)= Vknee2(ctr:n)
   END IF

!! Here the poles are calculated for the nonlinear case
IF (Nonlinear) THEN
  IF (SheraNonlinearityType == DISP) THEN    
     !nonlinearity threshold for the displacement is proportional to 1/omega,
       Yknee1CST(1:n)= RthY1(1:n) * omega(onektel) !normalized to 1 kHz location 
       Yknee2CST(1:n)= RthY2(1:n) * omega(onektel)
       Yknee1F(1:n)= Yknee1CST(1:n) / omega(1:n) !kneepoint function ifo omega
       Yknee2F(1:n)= Yknee2CST(1:n) / omega(1:n) !so that at 1kHz, YkneeF = Yknee
       
       Yvect(1:n)=ABS(Y(1:n))/Yknee1F(1:n)
       Theta0(1:n)=ATAN(((PoleE(1:n)-PoleS(1:n))*factor)/((Yknee2F(1:n)/Yknee1F(1:n))-1d0)) 
       Theta(1:n)=Theta0/2d0
       Sfoc(1:n)=(PoleS(1:n)*factor)/(Yknee2F(1:n)/Yknee1F(1:n))
       Se(1:n)=1d0/COS((pi-Theta0(1:n))/2d0)
       Sb(1:n)=Sfoc(1:n)/Se(1:n)
       Sa(1:n)=Sfoc(1:n)*DSQRT(1d0-(1d0/(Se(1:n)**2d0)))
       Sxp(1:n)=(Yvect(1:n)-1d0)*COS(Theta(1:n))/COS(2d0*Theta(1:n))
       Syp(1:n)=Sb(1:n)*DSQRT(1d0+(Sxp(1:n)/Sa(1:n))**2d0)
       Sy(1:n)=Sxp(1:n)*SIN(Theta(1:n))+Syp*COS(Theta(1:n))
       SheraP(1:n)=PoleS(1:n)+Sy(1:n)/factor
  ENDIF !end for the displacement nonlinearity
  
  IF (SheraNonlinearityType == VEL) THEN 

     Vvect(1:n)=ABS(V(1:n))/RthV1(1:n)
     Theta0(1:n)=ATAN(((PoleE(1:n)-PoleS(1:n))*factor)/((RthV2(1:n)/RthV1(1:n))-1d0)) 
     Theta(1:n)=Theta0/2d0
     Sfoc(1:n)=(PoleS(1:n)*factor)/(RthV2(1:n)/RthV1(1:n))
     Se(1:n)=1d0/COS((pi-Theta0(1:n))/2d0)
     Sb(1:n)=Sfoc(1:n)/Se(1:n)
     Sa(1:n)=Sfoc(1:n)*DSQRT(1d0-(1d0/(Se(1:n)**2d0)))
     Sxp(1:n)=(Vvect(1:n)-1d0)*COS(Theta(1:n))/COS(2d0*Theta(1:n))
     Syp(1:n)=Sb(1:n)*DSQRT(1d0+(Sxp(1:n)/Sa(1:n))**2d0)
     Sy(1:n)=Sxp(1:n)*SIN(Theta(1:n))+Syp*COS(Theta(1:n))
     SheraP(1:n)=PoleS(1:n)+Sy(1:n)/factor

  ENDIF
     
  IF (SheraNonlinearityType == ENV) THEN
     SheraP(1:n)=0.061d0 !is not implemented yet, so is set to active,linear behaviour for now
  ENDIF

!to stop the poles from increasing after a certain value lead to linear behavior at high levels
  DO i=1,n
     IF(SheraP(i).GT.PoleE(i)) THEN
        SheraP(i)=PoleE(i)
     ENDIF
  ENDDO   

ELSE !The poles for the linear implementation
   SheraP(1:n)=PoleS(1:n)
ENDIF

debug(1:n)=PoleS(1:n) !SheraP(1:n)
END SUBROUTINE PoleCalculation

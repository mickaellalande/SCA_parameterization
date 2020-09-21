MODULE cloudth_mod

  IMPLICIT NONE

CONTAINS

       SUBROUTINE cloudth(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)


      IMPLICIT NONE


!===========================================================================
! Auteur : Arnaud Octavio Jam (LMD/CNRS)
! Date : 25 Mai 2010
! Objet : calcule les valeurs de qc et rneb dans les thermiques
!===========================================================================


#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"
#include "nuage.h"

      INTEGER itap,ind1,ind2
      INTEGER ngrid,klev,klon,l,ig 
      
      REAL ztv(ngrid,klev)
      REAL po(ngrid)
      REAL zqenv(ngrid)   
      REAL zqta(ngrid,klev)
          
      REAL fraca(ngrid,klev+1)
      REAL zpspsk(ngrid,klev)
      REAL paprs(ngrid,klev+1)
      REAL ztla(ngrid,klev)
      REAL zthl(ngrid,klev)

      REAL zqsatth(ngrid,klev)
      REAL zqsatenv(ngrid,klev)
      
      
      REAL sigma1(ngrid,klev)
      REAL sigma2(ngrid,klev)
      REAL qlth(ngrid,klev)
      REAL qlenv(ngrid,klev)
      REAL qltot(ngrid,klev) 
      REAL cth(ngrid,klev)  
      REAL cenv(ngrid,klev)   
      REAL ctot(ngrid,klev)
      REAL rneb(ngrid,klev)
      REAL t(ngrid,klev)
      REAL qsatmmussig1,qsatmmussig2,sqrt2pi,pi
      REAL rdd,cppd,Lv
      REAL alth,alenv,ath,aenv
      REAL sth,senv,sigma1s,sigma2s,xth,xenv
      REAL Tbef,zdelta,qsatbef,zcor
      REAL qlbef  
      REAL ratqs(ngrid,klev) ! determine la largeur de distribution de vapeur
      
      REAL zpdf_sig(ngrid),zpdf_k(ngrid),zpdf_delta(ngrid)
      REAL zpdf_a(ngrid),zpdf_b(ngrid),zpdf_e1(ngrid),zpdf_e2(ngrid) 
      REAL zqs(ngrid), qcloud(ngrid)
      REAL erf




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gestion de deux versions de cloudth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (iflag_cloudth_vert.GE.1) THEN
      CALL cloudth_vert(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)
      RETURN
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-------------------------------------------------------------------------------
! Initialisation des variables r?elles
!-------------------------------------------------------------------------------
      sigma1(:,:)=0.
      sigma2(:,:)=0.
      qlth(:,:)=0.
      qlenv(:,:)=0.  
      qltot(:,:)=0.
      rneb(:,:)=0.
      qcloud(:)=0.
      cth(:,:)=0.
      cenv(:,:)=0.
      ctot(:,:)=0.
      qsatmmussig1=0.
      qsatmmussig2=0.
      rdd=287.04
      cppd=1005.7
      pi=3.14159 
      Lv=2.5e6
      sqrt2pi=sqrt(2.*pi)



!-------------------------------------------------------------------------------
! Calcul de la fraction du thermique et des ?cart-types des distributions
!-------------------------------------------------------------------------------                 
      do ind1=1,ngrid

      if ((ztv(ind1,1).gt.ztv(ind1,2)).and.(fraca(ind1,ind2).gt.1.e-10)) then 

      zqenv(ind1)=(po(ind1)-fraca(ind1,ind2)*zqta(ind1,ind2))/(1.-fraca(ind1,ind2))


!      zqenv(ind1)=po(ind1)
      Tbef=zthl(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef




      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2)) 




      Tbef=ztla(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatth(ind1,ind2)=qsatbef
            
      alth=(0.622*Lv*zqsatth(ind1,ind2))/(rdd*ztla(ind1,ind2)**2)   
      ath=1./(1.+(alth*Lv/cppd))
      sth=ath*(zqta(ind1,ind2)-zqsatth(ind1,ind2)) 
      
      

!------------------------------------------------------------------------------
! Calcul des ?cart-types pour s
!------------------------------------------------------------------------------

!      sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+ratqs(ind1,ind2)*po(ind1)
!      sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.4+0.002*zqta(ind1,ind2) 
!       if (paprs(ind1,ind2).gt.90000) then
!       ratqs(ind1,ind2)=0.002
!       else 
!       ratqs(ind1,ind2)=0.002+0.0*(90000-paprs(ind1,ind2))/20000
!       endif
       sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+0.002*po(ind1)
       sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.01)**0.4+0.002*zqta(ind1,ind2) 
!       sigma1s=ratqs(ind1,ind2)*po(ind1)
!      sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.4+0.00003  
 
!------------------------------------------------------------------------------
! Calcul de l'eau condens?e et de la couverture nuageuse
!------------------------------------------------------------------------------
      sqrt2pi=sqrt(2.*pi)
      xth=sth/(sqrt(2.)*sigma2s)
      xenv=senv/(sqrt(2.)*sigma1s)
      cth(ind1,ind2)=0.5*(1.+1.*erf(xth))
      cenv(ind1,ind2)=0.5*(1.+1.*erf(xenv)) 
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)   

      qlth(ind1,ind2)=sigma2s*((exp(-1.*xth**2)/sqrt2pi)+xth*sqrt(2.)*cth(ind1,ind2))
      qlenv(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt(2.)*cenv(ind1,ind2))   
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ctot(ind1,ind2).lt.1.e-10) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else   
                
      ctot(ind1,ind2)=ctot(ind1,ind2) 
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqs(ind1)

      endif                           
      
          
!     print*,sth,sigma2s,qlth(ind1,ind2),ctot(ind1,ind2),qltot(ind1,ind2),'verif'


      else  ! gaussienne environnement seule
      
      zqenv(ind1)=po(ind1)
      Tbef=t(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef
      

!      qlbef=Max(po(ind1)-zqsatenv(ind1,ind2),0.)
      zthl(ind1,ind2)=t(ind1,ind2)*(101325/paprs(ind1,ind2))**(rdd/cppd)
      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2)) 
      

      sigma1s=ratqs(ind1,ind2)*zqenv(ind1)

      sqrt2pi=sqrt(2.*pi)
      xenv=senv/(sqrt(2.)*sigma1s)
      ctot(ind1,ind2)=0.5*(1.+1.*erf(xenv))
      qltot(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt(2.)*cenv(ind1,ind2))
      
      if (ctot(ind1,ind2).lt.1.e-3) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else   
                
      ctot(ind1,ind2)=ctot(ind1,ind2) 
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqsatenv(ind1,ind2)

      endif    
 
 
 
 
 
 
      endif   
      enddo
     
      return
!     end
END SUBROUTINE cloudth



!===========================================================================
     SUBROUTINE cloudth_vert(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)

!===========================================================================
! Auteur : Arnaud Octavio Jam (LMD/CNRS)
! Date : 25 Mai 2010
! Objet : calcule les valeurs de qc et rneb dans les thermiques
!===========================================================================


      USE ioipsl_getin_p_mod, ONLY : getin_p

      IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"
#include "nuage.h"
      
      INTEGER itap,ind1,ind2
      INTEGER ngrid,klev,klon,l,ig 
      
      REAL ztv(ngrid,klev)
      REAL po(ngrid)
      REAL zqenv(ngrid)   
      REAL zqta(ngrid,klev)
          
      REAL fraca(ngrid,klev+1)
      REAL zpspsk(ngrid,klev)
      REAL paprs(ngrid,klev+1)
      REAL ztla(ngrid,klev)
      REAL zthl(ngrid,klev)

      REAL zqsatth(ngrid,klev)
      REAL zqsatenv(ngrid,klev)
      
      
      REAL sigma1(ngrid,klev)                                                         
      REAL sigma2(ngrid,klev)
      REAL qlth(ngrid,klev)
      REAL qlenv(ngrid,klev)
      REAL qltot(ngrid,klev) 
      REAL cth(ngrid,klev)  
      REAL cenv(ngrid,klev)   
      REAL ctot(ngrid,klev)
      REAL rneb(ngrid,klev)
      REAL t(ngrid,klev)                                                                  
      REAL qsatmmussig1,qsatmmussig2,sqrt2pi,pi
      REAL rdd,cppd,Lv,sqrt2,sqrtpi
      REAL alth,alenv,ath,aenv
      REAL sth,senv,sigma1s,sigma2s,xth,xenv
      REAL xth1,xth2,xenv1,xenv2,deltasth, deltasenv
      REAL IntJ,IntI1,IntI2,IntI3,coeffqlenv,coeffqlth
      REAL Tbef,zdelta,qsatbef,zcor
      REAL qlbef  
      REAL ratqs(ngrid,klev) ! determine la largeur de distribution de vapeur
      ! Change the width of the PDF used for vertical subgrid scale heterogeneity
      ! (J Jouhaud, JL Dufresne, JB Madeleine)
      REAL,SAVE :: vert_alpha
      !$OMP THREADPRIVATE(vert_alpha)
      LOGICAL, SAVE :: firstcall = .TRUE.
      !$OMP THREADPRIVATE(firstcall)
      
      REAL zpdf_sig(ngrid),zpdf_k(ngrid),zpdf_delta(ngrid)
      REAL zpdf_a(ngrid),zpdf_b(ngrid),zpdf_e1(ngrid),zpdf_e2(ngrid) 
      REAL zqs(ngrid), qcloud(ngrid)
      REAL erf

!------------------------------------------------------------------------------
! Initialisation des variables r?elles
!------------------------------------------------------------------------------
      sigma1(:,:)=0.
      sigma2(:,:)=0.
      qlth(:,:)=0.
      qlenv(:,:)=0.  
      qltot(:,:)=0.
      rneb(:,:)=0.
      qcloud(:)=0.
      cth(:,:)=0.
      cenv(:,:)=0.
      ctot(:,:)=0.
      qsatmmussig1=0.
      qsatmmussig2=0.
      rdd=287.04
      cppd=1005.7
      pi=3.14159 
      Lv=2.5e6
      sqrt2pi=sqrt(2.*pi)
      sqrt2=sqrt(2.)
      sqrtpi=sqrt(pi)

      IF (firstcall) THEN
        vert_alpha=0.5
        CALL getin_p('cloudth_vert_alpha',vert_alpha)
        WRITE(*,*) 'cloudth_vert_alpha = ', vert_alpha
        firstcall=.FALSE.
      ENDIF

!-------------------------------------------------------------------------------
! Calcul de la fraction du thermique et des ?cart-types des distributions
!-------------------------------------------------------------------------------                 
      do ind1=1,ngrid

      if ((ztv(ind1,1).gt.ztv(ind1,2)).and.(fraca(ind1,ind2).gt.1.e-10)) then 

      zqenv(ind1)=(po(ind1)-fraca(ind1,ind2)*zqta(ind1,ind2))/(1.-fraca(ind1,ind2))


!      zqenv(ind1)=po(ind1)
      Tbef=zthl(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef




      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2)) 




      Tbef=ztla(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatth(ind1,ind2)=qsatbef
            
      alth=(0.622*Lv*zqsatth(ind1,ind2))/(rdd*ztla(ind1,ind2)**2)   
      ath=1./(1.+(alth*Lv/cppd))
      sth=ath*(zqta(ind1,ind2)-zqsatth(ind1,ind2)) 
      
      

!------------------------------------------------------------------------------
! Calcul des ?cart-types pour s
!------------------------------------------------------------------------------

      sigma1s=(0.92**0.5)*(fraca(ind1,ind2)**0.5)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+ratqs(ind1,ind2)*po(ind1)
      sigma2s=0.09*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.5+0.002*zqta(ind1,ind2) 
!       if (paprs(ind1,ind2).gt.90000) then
!       ratqs(ind1,ind2)=0.002
!       else 
!       ratqs(ind1,ind2)=0.002+0.0*(90000-paprs(ind1,ind2))/20000
!       endif
!       sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+0.002*po(ind1)
!       sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.01)**0.4+0.002*zqta(ind1,ind2) 
!       sigma1s=ratqs(ind1,ind2)*po(ind1)
!      sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.4+0.00003  
 
!------------------------------------------------------------------------------
! Calcul de l'eau condens?e et de la couverture nuageuse
!------------------------------------------------------------------------------
      sqrt2pi=sqrt(2.*pi)
      xth=sth/(sqrt(2.)*sigma2s)
      xenv=senv/(sqrt(2.)*sigma1s)
      cth(ind1,ind2)=0.5*(1.+1.*erf(xth))
      cenv(ind1,ind2)=0.5*(1.+1.*erf(xenv)) 
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)   

      qlth(ind1,ind2)=sigma2s*((exp(-1.*xth**2)/sqrt2pi)+xth*sqrt(2.)*cth(ind1,ind2))
      qlenv(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt(2.)*cenv(ind1,ind2))   
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)
     
       IF (iflag_cloudth_vert == 1) THEN
!-------------------------------------------------------------------------------
!  Version 2: Modification selon J.-Louis. On condense ?? partir de qsat-ratqs
!-------------------------------------------------------------------------------
!      deltasenv=aenv*ratqs(ind1,ind2)*po(ind1)
!      deltasth=ath*ratqs(ind1,ind2)*zqta(ind1,ind2)
      deltasenv=aenv*ratqs(ind1,ind2)*zqsatenv(ind1,ind2)
      deltasth=ath*ratqs(ind1,ind2)*zqsatth(ind1,ind2)
!      deltasenv=aenv*0.01*po(ind1)
!     deltasth=ath*0.01*zqta(ind1,ind2)    
      xenv1=(senv-deltasenv)/(sqrt(2.)*sigma1s)
      xenv2=(senv+deltasenv)/(sqrt(2.)*sigma1s)
      xth1=(sth-deltasth)/(sqrt(2.)*sigma2s)
      xth2=(sth+deltasth)/(sqrt(2.)*sigma2s)
      coeffqlenv=(sigma1s)**2/(2*sqrtpi*deltasenv)
      coeffqlth=(sigma2s)**2/(2*sqrtpi*deltasth)
      
      cth(ind1,ind2)=0.5*(1.+1.*erf(xth2))
      cenv(ind1,ind2)=0.5*(1.+1.*erf(xenv2)) 
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)   

      IntJ=sigma1s*(exp(-1.*xenv1**2)/sqrt2pi)+0.5*senv*(1+erf(xenv1))
      IntI1=coeffqlenv*0.5*(0.5*sqrtpi*(erf(xenv2)-erf(xenv1))+xenv1*exp(-1.*xenv1**2)-xenv2*exp(-1.*xenv2**2))
      IntI2=coeffqlenv*xenv2*(exp(-1.*xenv2**2)-exp(-1.*xenv1**2))
      IntI3=coeffqlenv*0.5*sqrtpi*xenv2**2*(erf(xenv2)-erf(xenv1))

      qlenv(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
!      qlenv(ind1,ind2)=IntJ 
!      print*, qlenv(ind1,ind2),'VERIF EAU'


      IntJ=sigma2s*(exp(-1.*xth1**2)/sqrt2pi)+0.5*sth*(1+erf(xth1))
!      IntI1=coeffqlth*((0.5*xth1-xth2)*exp(-1.*xth1**2)+0.5*xth2*exp(-1.*xth2**2))
!      IntI2=coeffqlth*0.5*sqrtpi*(0.5+xth2**2)*(erf(xth2)-erf(xth1))
      IntI1=coeffqlth*0.5*(0.5*sqrtpi*(erf(xth2)-erf(xth1))+xth1*exp(-1.*xth1**2)-xth2*exp(-1.*xth2**2))
      IntI2=coeffqlth*xth2*(exp(-1.*xth2**2)-exp(-1.*xth1**2))
      IntI3=coeffqlth*0.5*sqrtpi*xth2**2*(erf(xth2)-erf(xth1))
      qlth(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
!      qlth(ind1,ind2)=IntJ
!      print*, IntJ,IntI1,IntI2,IntI3,qlth(ind1,ind2),'VERIF EAU2'
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)

      ELSE IF (iflag_cloudth_vert == 2) THEN

!-------------------------------------------------------------------------------
!  Version 3: Modification Jean Jouhaud. On condense a partir de -delta s
!-------------------------------------------------------------------------------
!      deltasenv=aenv*ratqs(ind1,ind2)*po(ind1)
!      deltasth=ath*ratqs(ind1,ind2)*zqta(ind1,ind2)
!      deltasenv=aenv*ratqs(ind1,ind2)*zqsatenv(ind1,ind2)
!      deltasth=ath*ratqs(ind1,ind2)*zqsatth(ind1,ind2)
      deltasenv=aenv*vert_alpha*sigma1s
      deltasth=ath*vert_alpha*sigma2s
      
      xenv1=-(senv+deltasenv)/(sqrt(2.)*sigma1s)
      xenv2=-(senv-deltasenv)/(sqrt(2.)*sigma1s)
      xth1=-(sth+deltasth)/(sqrt(2.)*sigma2s)
      xth2=-(sth-deltasth)/(sqrt(2.)*sigma2s)
!     coeffqlenv=(sigma1s)**2/(2*sqrtpi*deltasenv)
!     coeffqlth=(sigma2s)**2/(2*sqrtpi*deltasth)
      
      cth(ind1,ind2)=0.5*(1.-1.*erf(xth1))
      cenv(ind1,ind2)=0.5*(1.-1.*erf(xenv1))
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)

      IntJ=0.5*senv*(1-erf(xenv2))+(sigma1s/sqrt2pi)*exp(-1.*xenv2**2)
      IntI1=(((senv+deltasenv)**2+(sigma1s)**2)/(8*deltasenv))*(erf(xenv2)-erf(xenv1))
      IntI2=(sigma1s**2/(4*deltasenv*sqrtpi))*(xenv1*exp(-1.*xenv1**2)-xenv2*exp(-1.*xenv2**2))
      IntI3=((sqrt2*sigma1s*(senv+deltasenv))/(4*sqrtpi*deltasenv))*(exp(-1.*xenv1**2)-exp(-1.*xenv2**2))

!      IntI1=0.5*(0.5*sqrtpi*(erf(xenv2)-erf(xenv1))+xenv1*exp(-1.*xenv1**2)-xenv2*exp(-1.*xenv2**2))
!      IntI2=xenv2*(exp(-1.*xenv2**2)-exp(-1.*xenv1**2))
!      IntI3=0.5*sqrtpi*xenv2**2*(erf(xenv2)-erf(xenv1))

      qlenv(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
!      qlenv(ind1,ind2)=IntJ 
!      print*, qlenv(ind1,ind2),'VERIF EAU'

      IntJ=0.5*sth*(1-erf(xth2))+(sigma2s/sqrt2pi)*exp(-1.*xth2**2)
      IntI1=(((sth+deltasth)**2+(sigma2s)**2)/(8*deltasth))*(erf(xth2)-erf(xth1))
      IntI2=(sigma2s**2/(4*deltasth*sqrtpi))*(xth1*exp(-1.*xth1**2)-xth2*exp(-1.*xth2**2))
      IntI3=((sqrt2*sigma2s*(sth+deltasth))/(4*sqrtpi*deltasth))*(exp(-1.*xth1**2)-exp(-1.*xth2**2))
      
      qlth(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
!      qlth(ind1,ind2)=IntJ
!      print*, IntJ,IntI1,IntI2,IntI3,qlth(ind1,ind2),'VERIF EAU2'
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)
      



      ENDIF ! of if (iflag_cloudth_vert==1 or 2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (cenv(ind1,ind2).lt.1.e-10.or.cth(ind1,ind2).lt.1.e-10) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else 
                
      ctot(ind1,ind2)=ctot(ind1,ind2) 
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqs(ind1)
!      qcloud(ind1)=fraca(ind1,ind2)*qlth(ind1,ind2)/cth(ind1,ind2) &
!    &             +(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)/cenv(ind1,ind2)+zqs(ind1) 

      endif  
                       
      
          
!     print*,sth,sigma2s,qlth(ind1,ind2),ctot(ind1,ind2),qltot(ind1,ind2),'verif'


      else  ! gaussienne environnement seule
      
      zqenv(ind1)=po(ind1)
      Tbef=t(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef
      

!      qlbef=Max(po(ind1)-zqsatenv(ind1,ind2),0.)
      zthl(ind1,ind2)=t(ind1,ind2)*(101325/paprs(ind1,ind2))**(rdd/cppd)
      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2)) 
      

      sigma1s=ratqs(ind1,ind2)*zqenv(ind1)

      sqrt2pi=sqrt(2.*pi)
      xenv=senv/(sqrt(2.)*sigma1s)
      ctot(ind1,ind2)=0.5*(1.+1.*erf(xenv))
      qltot(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt(2.)*cenv(ind1,ind2))
      
      if (ctot(ind1,ind2).lt.1.e-3) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else   
                
      ctot(ind1,ind2)=ctot(ind1,ind2) 
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqsatenv(ind1,ind2)

      endif    
 
 
 
 
 
 
      endif   
      enddo
     
      return
!     end
END SUBROUTINE cloudth_vert

       SUBROUTINE cloudth_v3(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,ctot_vol,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)


      IMPLICIT NONE


!===========================================================================
! Author : Arnaud Octavio Jam (LMD/CNRS)
! Date : 25 Mai 2010
! Objet : calcule les valeurs de qc et rneb dans les thermiques
!===========================================================================


#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"
#include "nuage.h"

      INTEGER itap,ind1,ind2
      INTEGER ngrid,klev,klon,l,ig 
      
      REAL ztv(ngrid,klev)
      REAL po(ngrid)
      REAL zqenv(ngrid)   
      REAL zqta(ngrid,klev)
          
      REAL fraca(ngrid,klev+1)
      REAL zpspsk(ngrid,klev)
      REAL paprs(ngrid,klev+1)
      REAL ztla(ngrid,klev)
      REAL zthl(ngrid,klev)

      REAL zqsatth(ngrid,klev)
      REAL zqsatenv(ngrid,klev)
      
      REAL sigma1(ngrid,klev)                                                         
      REAL sigma2(ngrid,klev)
      REAL qlth(ngrid,klev)
      REAL qlenv(ngrid,klev)
      REAL qltot(ngrid,klev) 
      REAL cth(ngrid,klev)
      REAL cenv(ngrid,klev)   
      REAL ctot(ngrid,klev)
      REAL cth_vol(ngrid,klev)
      REAL cenv_vol(ngrid,klev)
      REAL ctot_vol(ngrid,klev)
      REAL rneb(ngrid,klev)      
      REAL t(ngrid,klev)
      REAL qsatmmussig1,qsatmmussig2,sqrt2pi,sqrt2,sqrtpi,pi
      REAL rdd,cppd,Lv
      REAL alth,alenv,ath,aenv
      REAL sth,senv,sigma1s,sigma2s,xth,xenv, exp_xenv1, exp_xenv2,exp_xth1,exp_xth2
      REAL Tbef,zdelta,qsatbef,zcor
      REAL qlbef  
      REAL ratqs(ngrid,klev) ! Determine the width of the vapour distribution
      REAL zpdf_sig(ngrid),zpdf_k(ngrid),zpdf_delta(ngrid)
      REAL zpdf_a(ngrid),zpdf_b(ngrid),zpdf_e1(ngrid),zpdf_e2(ngrid) 
      REAL zqs(ngrid), qcloud(ngrid)
      REAL erf



      IF (iflag_cloudth_vert.GE.1) THEN
      CALL cloudth_vert_v3(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,ctot_vol,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)
      RETURN
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-------------------------------------------------------------------------------
! Initialisation des variables r?elles
!-------------------------------------------------------------------------------
      sigma1(:,:)=0.
      sigma2(:,:)=0.
      qlth(:,:)=0.
      qlenv(:,:)=0.  
      qltot(:,:)=0.
      rneb(:,:)=0.
      qcloud(:)=0.
      cth(:,:)=0.
      cenv(:,:)=0.
      ctot(:,:)=0.
      cth_vol(:,:)=0.
      cenv_vol(:,:)=0.
      ctot_vol(:,:)=0.
      qsatmmussig1=0.
      qsatmmussig2=0.
      rdd=287.04
      cppd=1005.7
      pi=3.14159 
      Lv=2.5e6
      sqrt2pi=sqrt(2.*pi)
      sqrt2=sqrt(2.)
      sqrtpi=sqrt(pi)


!-------------------------------------------------------------------------------
! Cloud fraction in the thermals and standard deviation of the PDFs
!-------------------------------------------------------------------------------                 
      do ind1=1,ngrid

      if ((ztv(ind1,1).gt.ztv(ind1,2)).and.(fraca(ind1,ind2).gt.1.e-10)) then 

      zqenv(ind1)=(po(ind1)-fraca(ind1,ind2)*zqta(ind1,ind2))/(1.-fraca(ind1,ind2))

      Tbef=zthl(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES*FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef


      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)     !qsl, p84
      aenv=1./(1.+(alenv*Lv/cppd))                                      !al, p84
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2))                          !s, p84

!po = qt de l'environnement ET des thermique
!zqenv = qt environnement
!zqsatenv = qsat environnement
!zthl = Tl environnement


      Tbef=ztla(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatth(ind1,ind2)=qsatbef
            
      alth=(0.622*Lv*zqsatth(ind1,ind2))/(rdd*ztla(ind1,ind2)**2)       !qsl, p84
      ath=1./(1.+(alth*Lv/cppd))                                        !al, p84
      sth=ath*(zqta(ind1,ind2)-zqsatth(ind1,ind2))                      !s, p84
      
!zqta = qt thermals
!zqsatth = qsat thermals
!ztla = Tl thermals

!------------------------------------------------------------------------------
! s standard deviations
!------------------------------------------------------------------------------

!     tests
!     sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+0.002*po(ind1)
!     sigma1s=(0.92*(fraca(ind1,ind2)**0.5)/(1-fraca(ind1,ind2))*(((sth-senv)**2)**0.5))+ratqs(ind1,ind2)*po(ind1)
!     sigma2s=(0.09*(((sth-senv)**2)**0.5)/((fraca(ind1,ind2)+0.02)**0.5))+0.002*zqta(ind1,ind2)
!     final option
      sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+ratqs(ind1,ind2)*po(ind1)
      sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.01)**0.4+0.002*zqta(ind1,ind2)
 
!------------------------------------------------------------------------------
! Condensed water and cloud cover
!------------------------------------------------------------------------------
      xth=sth/(sqrt2*sigma2s)
      xenv=senv/(sqrt2*sigma1s)
      cth(ind1,ind2)=0.5*(1.+1.*erf(xth))       !4.18 p 111, l.7 p115 & 4.20 p 119 thesis Arnaud Jam
      cenv(ind1,ind2)=0.5*(1.+1.*erf(xenv))     !4.18 p 111, l.7 p115 & 4.20 p 119 thesis Arnaud Jam
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)
      ctot_vol(ind1,ind2)=ctot(ind1,ind2)

      qlth(ind1,ind2)=sigma2s*((exp(-1.*xth**2)/sqrt2pi)+xth*sqrt2*cth(ind1,ind2))
      qlenv(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt2*cenv(ind1,ind2))
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)

      if (ctot(ind1,ind2).lt.1.e-10) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 
      else
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqs(ind1)
      endif

      else  ! Environnement only, follow the if l.110
      
      zqenv(ind1)=po(ind1)
      Tbef=t(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef

!     qlbef=Max(po(ind1)-zqsatenv(ind1,ind2),0.)
      zthl(ind1,ind2)=t(ind1,ind2)*(101325/paprs(ind1,ind2))**(rdd/cppd)
      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2))
      
      sigma1s=ratqs(ind1,ind2)*zqenv(ind1)

      xenv=senv/(sqrt2*sigma1s)
      ctot(ind1,ind2)=0.5*(1.+1.*erf(xenv))
      ctot_vol(ind1,ind2)=ctot(ind1,ind2)
      qltot(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt2*cenv(ind1,ind2))

      if (ctot(ind1,ind2).lt.1.e-3) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2)
      else   
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqsatenv(ind1,ind2)
      endif


      endif       ! From the separation (thermal/envrionnement) et (environnement) only, l.110 et l.183
      enddo       ! from the loop on ngrid l.108
      return
!     end
END SUBROUTINE cloudth_v3



!===========================================================================
     SUBROUTINE cloudth_vert_v3(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,ctot_vol,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)

!===========================================================================
! Auteur : Arnaud Octavio Jam (LMD/CNRS)
! Date : 25 Mai 2010
! Objet : calcule les valeurs de qc et rneb dans les thermiques
!===========================================================================


      USE ioipsl_getin_p_mod, ONLY : getin_p
      USE phys_output_var_mod, ONLY : cloudth_sth,cloudth_senv, &
     &                                cloudth_sigmath,cloudth_sigmaenv

      IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "thermcell.h"
#include "nuage.h"
      
      INTEGER itap,ind1,ind2
      INTEGER ngrid,klev,klon,l,ig 
      
      REAL ztv(ngrid,klev)
      REAL po(ngrid)
      REAL zqenv(ngrid)   
      REAL zqta(ngrid,klev)
          
      REAL fraca(ngrid,klev+1)
      REAL zpspsk(ngrid,klev)
      REAL paprs(ngrid,klev+1)
      REAL ztla(ngrid,klev)
      REAL zthl(ngrid,klev)

      REAL zqsatth(ngrid,klev)
      REAL zqsatenv(ngrid,klev)
      
      REAL sigma1(ngrid,klev)                                                         
      REAL sigma2(ngrid,klev)
      REAL qlth(ngrid,klev)
      REAL qlenv(ngrid,klev)
      REAL qltot(ngrid,klev) 
      REAL cth(ngrid,klev)
      REAL cenv(ngrid,klev)   
      REAL ctot(ngrid,klev)
      REAL cth_vol(ngrid,klev)
      REAL cenv_vol(ngrid,klev)
      REAL ctot_vol(ngrid,klev)
      REAL rneb(ngrid,klev)
      REAL t(ngrid,klev)                                                                  
      REAL qsatmmussig1,qsatmmussig2,sqrtpi,sqrt2,sqrt2pi,pi
      REAL rdd,cppd,Lv
      REAL alth,alenv,ath,aenv
      REAL sth,senv,sigma1s,sigma2s,sigma1s_fraca,sigma1s_ratqs
      REAL xth,xenv,exp_xenv1,exp_xenv2,exp_xth1,exp_xth2
      REAL xth1,xth2,xenv1,xenv2,deltasth, deltasenv
      REAL IntJ,IntI1,IntI2,IntI3,IntJ_CF,IntI1_CF,IntI3_CF,coeffqlenv,coeffqlth
      REAL Tbef,zdelta,qsatbef,zcor
      REAL qlbef  
      REAL ratqs(ngrid,klev) ! determine la largeur de distribution de vapeur
      ! Change the width of the PDF used for vertical subgrid scale heterogeneity
      ! (J Jouhaud, JL Dufresne, JB Madeleine)
      REAL,SAVE :: vert_alpha, vert_alpha_th
      !$OMP THREADPRIVATE(vert_alpha, vert_alpha_th)
      REAL,SAVE :: sigma1s_factor=1.1
      REAL,SAVE :: sigma1s_power=0.6
      REAL,SAVE :: cloudth_ratqsmin=-1.
      !$OMP THREADPRIVATE(sigma1s_factor,sigma1s_power,cloudth_ratqsmin)
      INTEGER, SAVE :: iflag_cloudth_vert_noratqs=0
      !$OMP THREADPRIVATE(iflag_cloudth_vert_noratqs)

      LOGICAL, SAVE :: firstcall = .TRUE.
      !$OMP THREADPRIVATE(firstcall)

      REAL zpdf_sig(ngrid),zpdf_k(ngrid),zpdf_delta(ngrid)
      REAL zpdf_a(ngrid),zpdf_b(ngrid),zpdf_e1(ngrid),zpdf_e2(ngrid) 
      REAL zqs(ngrid), qcloud(ngrid)
      REAL erf

!------------------------------------------------------------------------------
! Initialize
!------------------------------------------------------------------------------
      sigma1(:,:)=0.
      sigma2(:,:)=0.
      qlth(:,:)=0.
      qlenv(:,:)=0.  
      qltot(:,:)=0.
      rneb(:,:)=0.
      qcloud(:)=0.
      cth(:,:)=0.
      cenv(:,:)=0.
      ctot(:,:)=0.
      cth_vol(:,:)=0.
      cenv_vol(:,:)=0.
      ctot_vol(:,:)=0.
      qsatmmussig1=0.
      qsatmmussig2=0.
      rdd=287.04
      cppd=1005.7
      pi=3.14159 
      Lv=2.5e6
      sqrt2pi=sqrt(2.*pi)
      sqrt2=sqrt(2.)
      sqrtpi=sqrt(pi)

      IF (firstcall) THEN
        vert_alpha=0.5
        CALL getin_p('cloudth_vert_alpha',vert_alpha)
        WRITE(*,*) 'cloudth_vert_alpha = ', vert_alpha
        ! The factor used for the thermal is equal to that of the environment
        !   if nothing is explicitly specified in the def file
        vert_alpha_th=vert_alpha
        CALL getin_p('cloudth_vert_alpha_th',vert_alpha_th)
        WRITE(*,*) 'cloudth_vert_alpha_th = ', vert_alpha_th
        ! Factor used in the calculation of sigma1s
        CALL getin_p('cloudth_sigma1s_factor',sigma1s_factor)
        WRITE(*,*) 'cloudth_sigma1s_factor = ', sigma1s_factor
        ! Power used in the calculation of sigma1s
        CALL getin_p('cloudth_sigma1s_power',sigma1s_power)
        WRITE(*,*) 'cloudth_sigma1s_power = ', sigma1s_power
        ! Minimum value for the environmental air subgrid water distrib
        CALL getin_p('cloudth_ratqsmin',cloudth_ratqsmin)
        WRITE(*,*) 'cloudth_ratqsmin = ', cloudth_ratqsmin
        ! Remove the dependency to ratqs from the variance of the vertical PDF
        CALL getin_p('iflag_cloudth_vert_noratqs',iflag_cloudth_vert_noratqs)
        WRITE(*,*) 'iflag_cloudth_vert_noratqs = ', iflag_cloudth_vert_noratqs

        firstcall=.FALSE.
      ENDIF

!-------------------------------------------------------------------------------
! Calcul de la fraction du thermique et des ecart-types des distributions
!-------------------------------------------------------------------------------                 
      do ind1=1,ngrid

      if ((ztv(ind1,1).gt.ztv(ind1,2)).and.(fraca(ind1,ind2).gt.1.e-10)) then !Thermal and environnement

      zqenv(ind1)=(po(ind1)-fraca(ind1,ind2)*zqta(ind1,ind2))/(1.-fraca(ind1,ind2)) !qt = a*qtth + (1-a)*qtenv


      Tbef=zthl(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES*FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef


      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)     !qsl, p84
      aenv=1./(1.+(alenv*Lv/cppd))                                      !al, p84
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2))                          !s, p84

!zqenv = qt environnement
!zqsatenv = qsat environnement
!zthl = Tl environnement


      Tbef=ztla(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatth(ind1,ind2)=qsatbef
            
      alth=(0.622*Lv*zqsatth(ind1,ind2))/(rdd*ztla(ind1,ind2)**2)       !qsl, p84
      ath=1./(1.+(alth*Lv/cppd))                                        !al, p84
      sth=ath*(zqta(ind1,ind2)-zqsatth(ind1,ind2))                      !s, p84
      
      
!zqta = qt thermals
!zqsatth = qsat thermals
!ztla = Tl thermals

!------------------------------------------------------------------------------
! s standard deviation
!------------------------------------------------------------------------------

      sigma1s_fraca = (sigma1s_factor**0.5)*(fraca(ind1,ind2)**sigma1s_power) / &
     &                (1-fraca(ind1,ind2))*((sth-senv)**2)**0.5
!     sigma1s_fraca = (1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5
      IF (cloudth_ratqsmin>0.) THEN
         sigma1s_ratqs = cloudth_ratqsmin*po(ind1)
      ELSE
         sigma1s_ratqs = ratqs(ind1,ind2)*po(ind1)
      ENDIF
      sigma1s = sigma1s_fraca + sigma1s_ratqs
      sigma2s=(0.09*(((sth-senv)**2)**0.5)/((fraca(ind1,ind2)+0.02)**0.5))+0.002*zqta(ind1,ind2)
!      tests
!      sigma1s=(0.92**0.5)*(fraca(ind1,ind2)**0.5)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+ratqs(ind1,ind2)*po(ind1)
!      sigma1s=(0.92*(fraca(ind1,ind2)**0.5)/(1-fraca(ind1,ind2))*(((sth-senv)**2)**0.5))+0.002*zqenv(ind1)
!      sigma2s=0.09*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.5+0.002*zqta(ind1,ind2) 
!      sigma2s=(0.09*(((sth-senv)**2)**0.5)/((fraca(ind1,ind2)+0.02)**0.5))+ratqs(ind1,ind2)*zqta(ind1,ind2)
!       if (paprs(ind1,ind2).gt.90000) then
!       ratqs(ind1,ind2)=0.002
!       else 
!       ratqs(ind1,ind2)=0.002+0.0*(90000-paprs(ind1,ind2))/20000
!       endif
!       sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+0.002*po(ind1)
!       sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.01)**0.4+0.002*zqta(ind1,ind2) 
!       sigma1s=ratqs(ind1,ind2)*po(ind1)
!      sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.4+0.00003  
 
       IF (iflag_cloudth_vert == 1) THEN
!-------------------------------------------------------------------------------
!  Version 2: Modification from Arnaud Jam according to JL Dufrense. Condensate from qsat-ratqs
!-------------------------------------------------------------------------------

      deltasenv=aenv*ratqs(ind1,ind2)*zqsatenv(ind1,ind2)
      deltasth=ath*ratqs(ind1,ind2)*zqsatth(ind1,ind2)

      xenv1=(senv-deltasenv)/(sqrt(2.)*sigma1s)
      xenv2=(senv+deltasenv)/(sqrt(2.)*sigma1s)
      xth1=(sth-deltasth)/(sqrt(2.)*sigma2s)
      xth2=(sth+deltasth)/(sqrt(2.)*sigma2s)
      coeffqlenv=(sigma1s)**2/(2*sqrtpi*deltasenv)
      coeffqlth=(sigma2s)**2/(2*sqrtpi*deltasth)
      
      cth(ind1,ind2)=0.5*(1.+1.*erf(xth2))
      cenv(ind1,ind2)=0.5*(1.+1.*erf(xenv2)) 
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)   

      ! Environment
      IntJ=sigma1s*(exp(-1.*xenv1**2)/sqrt2pi)+0.5*senv*(1+erf(xenv1))
      IntI1=coeffqlenv*0.5*(0.5*sqrtpi*(erf(xenv2)-erf(xenv1))+xenv1*exp(-1.*xenv1**2)-xenv2*exp(-1.*xenv2**2))
      IntI2=coeffqlenv*xenv2*(exp(-1.*xenv2**2)-exp(-1.*xenv1**2))
      IntI3=coeffqlenv*0.5*sqrtpi*xenv2**2*(erf(xenv2)-erf(xenv1))

      qlenv(ind1,ind2)=IntJ+IntI1+IntI2+IntI3

      ! Thermal
      IntJ=sigma2s*(exp(-1.*xth1**2)/sqrt2pi)+0.5*sth*(1+erf(xth1))
      IntI1=coeffqlth*0.5*(0.5*sqrtpi*(erf(xth2)-erf(xth1))+xth1*exp(-1.*xth1**2)-xth2*exp(-1.*xth2**2))
      IntI2=coeffqlth*xth2*(exp(-1.*xth2**2)-exp(-1.*xth1**2))
      IntI3=coeffqlth*0.5*sqrtpi*xth2**2*(erf(xth2)-erf(xth1))
      qlth(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)

      ELSE IF (iflag_cloudth_vert >= 3) THEN

!-------------------------------------------------------------------------------
!  Version 3: Changes by J. Jouhaud; condensation for q > -delta s
!-------------------------------------------------------------------------------
!      deltasenv=aenv*ratqs(ind1,ind2)*po(ind1)
!      deltasth=ath*ratqs(ind1,ind2)*zqta(ind1,ind2)
!      deltasenv=aenv*ratqs(ind1,ind2)*zqsatenv(ind1,ind2)
!      deltasth=ath*ratqs(ind1,ind2)*zqsatth(ind1,ind2)
      IF (iflag_cloudth_vert == 3) THEN
        deltasenv=aenv*vert_alpha*sigma1s
        deltasth=ath*vert_alpha_th*sigma2s
      ELSE IF (iflag_cloudth_vert == 4) THEN
        IF (iflag_cloudth_vert_noratqs == 1) THEN
          deltasenv=vert_alpha*max(sigma1s_fraca,1e-10)
          deltasth=vert_alpha_th*sigma2s
        ELSE
          deltasenv=vert_alpha*sigma1s
          deltasth=vert_alpha_th*sigma2s
        ENDIF
      ENDIF
      
      xenv1=-(senv+deltasenv)/(sqrt(2.)*sigma1s)
      xenv2=-(senv-deltasenv)/(sqrt(2.)*sigma1s)
      exp_xenv1 = exp(-1.*xenv1**2)
      exp_xenv2 = exp(-1.*xenv2**2)
      xth1=-(sth+deltasth)/(sqrt(2.)*sigma2s)
      xth2=-(sth-deltasth)/(sqrt(2.)*sigma2s)
      exp_xth1 = exp(-1.*xth1**2)
      exp_xth2 = exp(-1.*xth2**2)
      
      !CF_surfacique
      cth(ind1,ind2)=0.5*(1.-1.*erf(xth1))
      cenv(ind1,ind2)=0.5*(1.-1.*erf(xenv1))
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2) 


      !CF_volumique & eau condense
      !environnement
      IntJ=0.5*senv*(1-erf(xenv2))+(sigma1s/sqrt2pi)*exp_xenv2
      IntJ_CF=0.5*(1.-1.*erf(xenv2))
      if (deltasenv .lt. 1.e-10) then 
      qlenv(ind1,ind2)=IntJ
      cenv_vol(ind1,ind2)=IntJ_CF
      else
      IntI1=(((senv+deltasenv)**2+(sigma1s)**2)/(8*deltasenv))*(erf(xenv2)-erf(xenv1))
      IntI2=(sigma1s**2/(4*deltasenv*sqrtpi))*(xenv1*exp_xenv1-xenv2*exp_xenv2)
      IntI3=((sqrt2*sigma1s*(senv+deltasenv))/(4*sqrtpi*deltasenv))*(exp_xenv1-exp_xenv2)
      IntI1_CF=((senv+deltasenv)*(erf(xenv2)-erf(xenv1)))/(4*deltasenv)
      IntI3_CF=(sqrt2*sigma1s*(exp_xenv1-exp_xenv2))/(4*sqrtpi*deltasenv)
      qlenv(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
      cenv_vol(ind1,ind2)=IntJ_CF+IntI1_CF+IntI3_CF
      endif

      !thermique
      IntJ=0.5*sth*(1-erf(xth2))+(sigma2s/sqrt2pi)*exp_xth2
      IntJ_CF=0.5*(1.-1.*erf(xth2))
      if (deltasth .lt. 1.e-10) then
      qlth(ind1,ind2)=IntJ
      cth_vol(ind1,ind2)=IntJ_CF
      else
      IntI1=(((sth+deltasth)**2+(sigma2s)**2)/(8*deltasth))*(erf(xth2)-erf(xth1))
      IntI2=(sigma2s**2/(4*deltasth*sqrtpi))*(xth1*exp_xth1-xth2*exp_xth2)
      IntI3=((sqrt2*sigma2s*(sth+deltasth))/(4*sqrtpi*deltasth))*(exp_xth1-exp_xth2)
      IntI1_CF=((sth+deltasth)*(erf(xth2)-erf(xth1)))/(4*deltasth)
      IntI3_CF=(sqrt2*sigma2s*(exp_xth1-exp_xth2))/(4*sqrtpi*deltasth)
      qlth(ind1,ind2)=IntJ+IntI1+IntI2+IntI3
      cth_vol(ind1,ind2)=IntJ_CF+IntI1_CF+IntI3_CF
      endif


      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)
      ctot_vol(ind1,ind2)=fraca(ind1,ind2)*cth_vol(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv_vol(ind1,ind2)

      ENDIF ! of if (iflag_cloudth_vert==1 or 3 or 4)


      if (cenv(ind1,ind2).lt.1.e-10.or.cth(ind1,ind2).lt.1.e-10) then
      ctot(ind1,ind2)=0.
      ctot_vol(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else 
                
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqs(ind1)
!      qcloud(ind1)=fraca(ind1,ind2)*qlth(ind1,ind2)/cth(ind1,ind2) &
!    &             +(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)/cenv(ind1,ind2)+zqs(ind1) 

      endif  

      else  ! Environment only
      
      zqenv(ind1)=po(ind1)
      Tbef=t(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef
      

!     qlbef=Max(po(ind1)-zqsatenv(ind1,ind2),0.)
      zthl(ind1,ind2)=t(ind1,ind2)*(101325/paprs(ind1,ind2))**(rdd/cppd)
      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2))
      sth=0.

      sigma1s=ratqs(ind1,ind2)*zqenv(ind1)
      sigma2s=0.

      xenv=senv/(sqrt2*sigma1s)
      ctot(ind1,ind2)=0.5*(1.+1.*erf(xenv))
      ctot_vol(ind1,ind2)=ctot(ind1,ind2)
      qltot(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt2*cenv(ind1,ind2))
      
      if (ctot(ind1,ind2).lt.1.e-3) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else   
                
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqsatenv(ind1,ind2)

      endif    
 
      endif       ! From the separation (thermal/envrionnement) et (environnement) only, l.335 et l.492
      ! Outputs used to check the PDFs
      cloudth_senv(ind1,ind2) = senv
      cloudth_sth(ind1,ind2) = sth
      cloudth_sigmaenv(ind1,ind2) = sigma1s
      cloudth_sigmath(ind1,ind2) = sigma2s

      enddo       ! from the loop on ngrid l.333
     
      return
!     end
END SUBROUTINE cloudth_vert_v3
!
END MODULE cloudth_mod

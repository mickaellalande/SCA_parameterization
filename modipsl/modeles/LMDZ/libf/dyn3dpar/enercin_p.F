      SUBROUTINE enercin_p ( vcov, ucov, vcont, ucont, ecin )
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur: P. Le Van
c   -------
c
c   Objet:
c   ------
c
c *********************************************************************
c .. calcul de l'energie cinetique aux niveaux s  ......
c *********************************************************************
c  vcov, vcont, ucov et ucont sont des arguments d'entree pour le s-pg .
c  ecin         est  un  argument de sortie pour le s-pg
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      REAL vcov( ip1jm,llm ),vcont( ip1jm,llm ),
     * ucov( ip1jmp1,llm ),ucont( ip1jmp1,llm ),ecin( ip1jmp1,llm )

      REAL ecinni( iip1 ),ecinsi( iip1 )

      REAL ecinpn, ecinps
      INTEGER     l,ij,i,ijb,ije

      EXTERNAL    SSUM
      REAL        SSUM



c                 . V
c                i,j-1

c      alpha4 .       . alpha1


c        U .      . P     . U
c       i-1,j    i,j      i,j

c      alpha3 .       . alpha2


c                 . V
c                i,j

c    
c  L'energie cinetique au point scalaire P(i,j) ,autre que les poles, est :
c       Ecin = 0.5 * U(i-1,j)**2 *( alpha3 + alpha4 )  +
c              0.5 * U(i  ,j)**2 *( alpha1 + alpha2 )  +
c              0.5 * V(i,j-1)**2 *( alpha1 + alpha4 )  +
c              0.5 * V(i,  j)**2 *( alpha2 + alpha3 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,llm
      
      ijb=ij_begin
      ije=ij_end+iip1
      
      IF (pole_nord) ijb=ij_begin+iip1
      IF (pole_sud)  ije=ij_end-iip1
      
      DO 1  ij = ijb, ije -1
      ecin( ij+1, l )  =    0.5  *
     * (   ucov( ij   ,l ) * ucont( ij   ,l ) * alpha3p4( ij +1 )   +
     *     ucov( ij+1 ,l ) * ucont( ij+1 ,l ) * alpha1p2( ij +1 )   +
     *     vcov(ij-iim,l ) * vcont(ij-iim,l ) * alpha1p4( ij +1 )   +
     *     vcov( ij+ 1,l ) * vcont( ij+ 1,l ) * alpha2p3( ij +1 )   )
   1  CONTINUE

c    ... correction pour  ecin(1,j,l)  ....
c    ...   ecin(1,j,l)= ecin(iip1,j,l) ...

CDIR$ IVDEP
      DO 2 ij = ijb, ije, iip1
      ecin( ij,l ) = ecin( ij + iim, l )
   2  CONTINUE

c     calcul aux poles  .......

      IF (pole_nord) THEN
    
        DO  i = 1, iim
         ecinni(i) = vcov(    i  ,  l) * 
     *               vcont(    i    ,l) * aire(   i   )
        ENDDO

        ecinpn = 0.5 * SSUM( iim,ecinni,1 ) / apoln

        DO ij = 1,iip1
          ecin(   ij     , l ) = ecinpn
        ENDDO
   
      ENDIF

      IF (pole_sud) THEN
    
        DO  i = 1, iim
         ecinsi(i) = vcov(i+ip1jmi1,l)* 
     *               vcont(i+ip1jmi1,l) * aire(i+ip1jm)
        ENDDO

        ecinps = 0.5 * SSUM( iim,ecinsi,1 ) / apols

        DO ij = 1,iip1
          ecin( ij+ ip1jm, l ) = ecinps
        ENDDO
   
      ENDIF

      
   5  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END

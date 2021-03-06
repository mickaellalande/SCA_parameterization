      SUBROUTINE divgrad2_loc ( klevel, h, deltapres, lh, divgra_out )
c
c     P. Le Van
c
c   ***************************************************************
c
c     .....   calcul de  (div( grad ))   de (  pext * h ) .....
c   ****************************************************************
c   h ,klevel,lh et pext  sont des arguments  d'entree pour le s-prg
c         divgra     est  un argument  de sortie pour le s-prg
c
      USE parallel_lmdz
      USE times
      USE mod_hallo
      USE divgrad2_mod
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "comdissipn.h"

c    .......    variables en arguments   .......
c
      INTEGER klevel
      REAL h( ijb_u:ije_u,klevel ), deltapres( ijb_u:ije_u,klevel )
      REAL divgra_out( ijb_u:ije_u,klevel)
c    .......    variables  locales    ..........
c
      REAL     signe, nudivgrs, sqrtps( ijb_u:ije_u,llm )
      INTEGER  l,ij,iter,lh
c    ...................................................................
      Type(Request),SAVE :: request_dissip
!$OMP THREADPRIVATE(request_dissip)
      INTEGER ijb,ije

c
c
      signe    = (-1.)**lh
      nudivgrs = signe * cdivh

c      CALL SCOPY ( ip1jmp1 * klevel, h, 1, divgra, 1 )
      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, klevel
        divgra(ijb:ije,l)=h(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
c
c$OMP BARRIER
       call Register_Hallo_u(divgra,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER

      CALL laplacien_loc( klevel, divgra, divgra )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)     
      DO l = 1, klevel
       DO ij = ijb, ije
        sqrtps( ij,l ) = SQRT( deltapres(ij,l) )
       ENDDO
      ENDDO
c$OMP END DO NOWAIT

c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)     
      DO l = 1, klevel
        DO ij = ijb, ije
         divgra(ij,l) = divgra(ij,l) * sqrtps(ij,l)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT
   
c    ........    Iteration de l'operateur  laplacien_gam    ........
c
      DO  iter = 1, lh - 2
c$OMP BARRIER
       call Register_Hallo_u(divgra,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)

c$OMP BARRIER


       CALL laplacien_gam_loc(klevel,cuvscvgam2,cvuscugam2,unsair_gam2,
     *                     unsapolnga2, unsapolsga2,  divgra, divgra )
      ENDDO
c
c    ...............................................................

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)     
      DO l = 1, klevel
        DO ij = ijb, ije
          divgra(ij,l) = divgra(ij,l) * sqrtps(ij,l)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT
c
c$OMP BARRIER
       call Register_Hallo_u(divgra,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER

      CALL laplacien_loc ( klevel, divgra, divgra )
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l  = 1,klevel
      DO ij = ijb,ije
      divgra_out(ij,l) =  nudivgrs * divgra(ij,l) / deltapres(ij,l)
      ENDDO
      ENDDO
c$OMP END DO NOWAIT

      RETURN
      END

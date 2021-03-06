       SUBROUTINE nxgraro2_loc(klevel,xcov,ycov,lr,grx_out,gry_out)
c
c      P.Le Van .
c   ***********************************************************
c                                 lr
c      calcul de  ( nxgrad (rot) )   du vect. v  ....
c
c       xcov et ycov  etant les compos. covariantes de  v
c   ***********************************************************
c     xcov , ycov et lr  sont des arguments  d'entree pour le s-prog
c      grx   et  gry     sont des arguments de sortie pour le s-prog
c
c
      USE write_Field_p
      USE parallel_lmdz
      USE times
      USE mod_hallo
      USE mod_filtreg_p
      USE nxgraro2_mod
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comdissipn.h"
c
c    ......  variables en arguments  .......
c
      INTEGER klevel
      REAL xcov( ijb_u:ije_u,klevel ), ycov( ijb_v:ije_v,klevel )
      REAL  grx_out( ijb_u:ije_u,klevel ),gry_out(ijb_v:ije_v,klevel)
c
c    ......   variables locales     ........
c
      REAL  signe, nugradrs
      INTEGER l,ij,iter,lr
      Type(Request),SAVE :: Request_dissip
!$OMP THREADPRIVATE(Request_dissip)
c    ........................................................
c
      INTEGER :: ijb,ije,jjb,jje
      
c
c
      signe    = (-1.)**lr
      nugradrs = signe * crot
c
c      CALL SCOPY ( ip1jmp1* klevel, xcov, 1, grx, 1 )
c      CALL SCOPY (  ip1jm * klevel, ycov, 1, gry, 1 )
 
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO    l = 1, klevel
        grx(ijb:ije,l)=xcov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT

c$OMP BARRIER
       call Register_Hallo_u(grx,llm,0,1,1,0,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER

      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO    l = 1, klevel
        gry(ijb:ije,l)=ycov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
 
c
      CALL     rotatf_loc ( klevel, grx, gry, rot )
c      call write_field3d_p('rot1',reshape(rot,(/iip1,jjm,llm/)))

c$OMP BARRIER
       call Register_Hallo_v(rot,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER
      
      CALL laplacien_rot_loc ( klevel, rot, rot,grx,gry      )
c       call write_field3d_p('rot2',reshape(rot,(/iip1,jjm,llm/)))
c
c    .....   Iteration de l'operateur laplacien_rotgam  .....
c
      DO  iter = 1, lr -2
c$OMP BARRIER
       call Register_Hallo_v(rot,llm,1,1,1,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER

	CALL laplacien_rotgam_loc( klevel, rot, rot )
      ENDDO
      
c       call write_field3d_p('rot3',reshape(rot,(/iip1,jjm,llm/)))
      
c
c
      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
       
      CALL filtreg_p( rot, jjb_v,jje_v, jjb,jje,jjm,
     &                klevel, 2,1, .FALSE.,1)
c$OMP BARRIER
       call Register_Hallo_v(rot,llm,1,0,0,1,Request_dissip)
       call SendRequest(Request_dissip)
c$OMP BARRIER
       call WaitRequest(Request_dissip)
c$OMP BARRIER

      CALL nxgrad_loc ( klevel, rot, grx, gry )

c
      ijb=ij_begin
      ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO    l = 1, klevel
        
         if(pole_sud) ije=ij_end-iip1
         DO  ij = ijb, ije
          gry_out( ij,l ) = gry( ij,l ) * nugradrs
         ENDDO
        
         if(pole_sud) ije=ij_end
         DO  ij = ijb, ije
          grx_out( ij,l ) = grx( ij,l ) * nugradrs
         ENDDO
     
      ENDDO
c$OMP END DO NOWAIT
c
      RETURN
      END

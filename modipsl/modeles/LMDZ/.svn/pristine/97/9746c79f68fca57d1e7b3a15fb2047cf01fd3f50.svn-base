      do l = 1, llm
!
! au dessus de 700hPa, on relaxe vers profil init
!      on fait l'hypothese que dans ce cas, il n'y a plus d'eau liq. au dessus 700hpa
!      donc la relaxation en thetal et qt devient relaxation en tempe et qv
       relax_u(l)=(u(l)-u_mod(l))/tau_sandu  ! pour u et v on relaxe sur tte la colonne
       relax_v(l)=(v(l)-v_mod(l))/tau_sandu
       relax_q(l,1)=0.
       relax_q(l,2)=0.
       relax_thl(l)=0.
!      print *,'nudge: l tau_sandu u u_mod',l,tau_sandu,u(l),u_mod(l)
!
       if (l.ge.llm700) then 
         relax_q(l,1)=(q(l,1)-q_mod(l))/tau_sandu
         relax_q(l,2)=0.
         relax_thl(l)=(temp(l)-t_mod(l))/tau_sandu
       endif
!     print *,'l dq1 relax dqadv',l,dq(l,1),relax_q(l,1),d_q_adv(l,1)
!     print *,'l dq2 relax dqadv',l,dq(l,2),relax_q(l,2),d_q_adv(l,2)
!     print *,'l dt  relax dtadv',l,dt_phys(l),relax_thl(l),d_t_adv(l)
      enddo
        u(1:mxcalc)=u(1:mxcalc) + timestep*(                                &          
     &              du_phys(1:mxcalc) - relax_u(1:mxcalc))
        v(1:mxcalc)=v(1:mxcalc) + timestep*(                                &
     &               dv_phys(1:mxcalc) - relax_v(1:mxcalc))
!       q(1:mxcalc,:)=q(1:mxcalc,:)+timestep*(
!    .              dq(1:mxcalc,:) - relax_q(1:mxcalc,:)+
!    .              d_q_adv(1:mxcalc,:))
        q(1:mxcalc,1)=q(1:mxcalc,1)+timestep*(                              &
     &      dq(1:mxcalc,1) - relax_q(1:mxcalc,1)+d_q_adv(1:mxcalc,1))
        q(1:mxcalc,2)=q(1:mxcalc,2)+timestep*(                              &
     &      dq(1:mxcalc,2) - relax_q(1:mxcalc,2)+d_q_adv(1:mxcalc,2))
        temp(1:mxcalc)=temp(1:mxcalc)+timestep*(                            &
     &      dt_phys(1:mxcalc)-relax_thl(1:mxcalc)+d_t_adv(1:mxcalc))


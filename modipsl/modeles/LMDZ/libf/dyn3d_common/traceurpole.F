!
! $Id: traceurpole.F 2622 2016-09-04 06:12:02Z emillour $
!
          subroutine traceurpole(q,masse)

          implicit none
      
      include "dimensions.h"
      include "paramet.h"
      include "comdissip.h"
      include "comgeom2.h"
      include "description.h"


c   Arguments
       integer iq
       real masse(iip1,jjp1,llm)
       real q(iip1,jjp1,llm)
       

c   Locals
      integer i,j,l
      real sommemassen(llm)
      real sommemqn(llm)
      real sommemasses(llm)
      real sommemqs(llm)
      real qpolen(llm),qpoles(llm)

    
c On impose une seule valeur au p�le Sud j=jjm+1=jjp1       
      sommemasses=0
      sommemqs=0
          do l=1,llm
             do i=1,iip1          
                 sommemasses(l)=sommemasses(l)+masse(i,jjp1,l)
                 sommemqs(l)=sommemqs(l)+masse(i,jjp1,l)*q(i,jjp1,l)
             enddo         
          qpoles(l)=sommemqs(l)/sommemasses(l)
          enddo

c On impose une seule valeur du traceur au p�le Nord j=1
      sommemassen=0
      sommemqn=0  
         do l=1,llm
           do i=1,iip1              
               sommemassen(l)=sommemassen(l)+masse(i,1,l)
               sommemqn(l)=sommemqn(l)+masse(i,1,l)*q(i,1,l)
           enddo
           qpolen(l)=sommemqn(l)/sommemassen(l) 
         enddo
    
c On force le traceur � prendre cette valeur aux p�les
        do l=1,llm
            do i=1,iip1
               q(i,1,l)=qpolen(l)
               q(i,jjp1,l)=qpoles(l)
             enddo
        enddo

      
      return
      end           

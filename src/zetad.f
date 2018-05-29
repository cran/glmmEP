cccccccccc FORTRAN subroutine zetad.f cccccccccc

c For computation of the first derivative of
c the "zeta" function as defined in the R package `sn'.

c Last changed: 14 AUG 2017

      subroutine zetad(x,ans)     
      double precision x,pi,rt2,ans,toler,tiny,ansprv,Cprv,Dprv,
     +                 aj,Ccur,Dcur,Delta   
      integer j

      if (x.gt.(-3.0)) then
         pi = 4.0*atan(1.0)
         rt2 = sqrt(2.0)      
         ans = 2.0*exp(-0.5*x*x)/(sqrt(2.0*pi)*erfc(-x/rt2))
      else
         toler = 1.0e-10
         tiny = 1.0e-30
         ansprv = tiny
         Cprv = tiny
         Dprv = 0.0
         Delta = 2.0 + toler
         j = 0
c
c        Top of iteration loop.
c      
10          j = j + 1
            aj = j - 1
            if (j.eq.1) then
               aj = 1.0
            endif
            Dcur = aj*Dprv - x
            if (abs(Dcur).lt.tiny) then
                Dcur = tiny
            endif
            Dcur = 1/Dcur
            Ccur = (aj/Cprv) - x
            if (abs(Ccur).lt.tiny) then
               Ccur = tiny
            endif
            Delta = Ccur*Dcur
            anscur = ansprv*Delta
            ansprv = anscur
            Cprv = Ccur
            Dprv = Dcur
            if (abs(Delta-1.0).ge.toler) then
               goto 10
            endif
         ans = 1/anscur
      endif

      return
      end

cccccccccc End of zetad.f cccccccccc

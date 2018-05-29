cccccccccc FORTRAN subroutine logphi.f cccccccccc

c For computation of the logarithm of the N(0,1) cumulative
c distribution function.

c Last changed: 07 FEB 2018

      subroutine logphi(x,ans)     
      double precision x,pi,zdv,ans
     
      if (x.gt.0.0) then
         ans = log(erfc(-x/sqrt(2.0))/2.0)
      else
         pi = 4.0*atan(1.0)
         call zetad(x,zdv)
         ans = -0.5*x*x - log(zdv) - 0.5*log(2*pi)
      endif

      return
      end

cccccccccc End of logphi.f cccccccccc

cccccccccc FORTRAN subroutine logdet.f cccccccccc

c For computation of the logarithm of the 
c determinant of a symmetric positive definite
c matrix.

c Last changed: 11 AUG 2017

      subroutine logdet(A,idmn,ipvt,work,det,ans)     
      double precision A(idmn,idmn),det(2),ans,work(idmn,idmn)
      integer idmn,ipvt(idmn),info     

      call dgefa(A,idmn,idmn,ipvt,info)
      call dgedi(A,idmn,idmn,ipvt,det,work,10)

      ans = log(abs(det(1))) + det(2)*log(10.0)
      
      return
      end

cccccccccc End of logdet.f cccccccccc

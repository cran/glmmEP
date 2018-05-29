cccccccccc FORTRAN subroutine cpbt.f cccccccccc

c Last changed: 24 MAY 2018

      subroutine cpbt(a1,a2,b1,b2,c0,c1,idmn,idmnsq,lena2,Dd,DdPlus,
     +                wka,wkb,A2ina1,B2inb1,A2inc1,ipvt,det,work,A2mat,
     +                A2neg,B2mat,B2neg,ans)     
      double precision a1(idmn),a2(lena2),b1(idmn),b2(lena2),c0,
     +                 c1(idmn),A2ina1(idmn),B2inb1(idmn),A2inc1(idmn),
     +                 wka(idmnsq),wkb(idmnsq),Dd(idmnsq,lena2),
     +                 DdPlus(lena2,idmnsq),A2mat(idmn,idmn),
     +                 A2neg(idmn,idmn),B2mat(idmn,idmn),
     +                 B2neg(idmn,idmn),a1A2a1,b1B2b1,c1A2a1,c1A2c1,
     +                 r2,xlgphr,det(2),work(idmn,idmn),xldmA2,xldmB2,
     +                 ans
      integer i,j,k,idmn,idmnsq,lena2,ipos,info,ipvt(idmn)

c Obtain D_d^+*a_2 and D_d^+*b_2 :

      do 10 j = 1,idmnsq
         wka(j) = 0.0
         wkb(j) = 0.0
         do 20 i = 1,lena2
            wka(j) = wka(j) + DdPlus(i,j)*a2(i)
            wkb(j) = wkb(j) + DdPlus(i,j)*b2(i)
20       continue
10    continue
      
c Obtain A_2 and B_2:

      ipos = 1
      do 30 j = 1,idmn
         do 40 i = 1,idmn
             A2mat(i,j) = wka(ipos)
             A2neg(i,j) = -wka(ipos)
             B2mat(i,j) = wkb(ipos)
             B2neg(i,j) = -wkb(ipos)
             ipos = ipos + 1
40       continue
30    continue

c Store a1 and b1 in A2ina1 and B2inb1:

      do 50 i = 1,idmn
         A2ina1(i) = a1(i)
         B2inb1(i) = b1(i)
         A2inc1(i) = c1(i)
50    continue

c Obtain A_2^(-1)a_1 and B_2^(-1)b_1 (noting storage
c of solution vectors in `a1' and `b1':

      call dgefa(A2mat,idmn,idmn,ipvt,info)
      call dgesl(A2mat,idmn,idmn,ipvt,A2ina1,0.0)
      call dgesl(A2mat,idmn,idmn,ipvt,A2inc1,0.0)

      call dgefa(B2mat,idmn,idmn,ipvt,info)
      call dgesl(B2mat,idmn,idmn,ipvt,B2inb1,0.0)

c Obtain a_1^TA_2^(-1)a_1 and b1_T^B_2^(-1)b1:

      a1A2a1 = 0.0
      b1B2b1 = 0.0
      c1A2a1 = 0.0
      c1A2c1 = 0.0
      do 60 i = 1,idmn
         a1A2a1 = a1A2a1 + a1(i)*A2ina1(i)
         b1B2b1 = b1B2b1 + b1(i)*B2inb1(i)
         c1A2a1 = c1A2a1 + a1(i)*A2inc1(i)
         c1A2c1 = c1A2c1 + c1(i)*A2inc1(i)
60    continue

c Obtain r2 and log{Phi(r2)}:

      r2 = (2*c0 - c1A2a1)/sqrt(2.0*(2.0-c1A2c1))
      call logPhi(r2,xlgphr) 

c Obtain log|-A2| and log|-B2|:

      call logdet(A2neg,idmn,ipvt,work,det,xldmA2)     
      call logdet(B2neg,idmn,ipvt,work,det,xldmB2)   

c Obtain answer:

      ans = xlgphr + 0.25*(b1B2b1 - a1A2a1) + 0.5*(xldmB2 - xldmA2)

      return
      end

cccccccccc End of cpbt.f cccccccccc

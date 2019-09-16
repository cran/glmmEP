cccccccccc FORTRAN subroutine kpbt.f cccccccccc

c Last changed: 28 MAY 2018

      subroutine kpbt(a1,a2,c0,c1,idmn,idmnsq,lena2,Dd,DdPlus,wk1,
     +                A2ina1,A2inc1,ipvt,A2mat,A2str,R2comp,wk2,
     +                R5,R5TA2,vR5TA2,ans1,ans2)     
      double precision a1(idmn),a2(lena2),c0,c1(idmn),A2ina1(idmn),
     +                 A2inc1(idmn),wk1(idmnsq),Dd(idmnsq,lena2),
     +                 DdPlus(lena2,idmnsq),A2mat(idmn,idmn),c1A2a1,
     +                 c1A2c1,r1,r2,r3,r4,zetdv,zetddv,
     +                 A2str(idmn,idmn),R2comp(idmn,idmn),wk2(idmn),
     +                 R5(idmn,idmn),R5TA2(idmn,idmn),vR5TA2(idmnsq),
     +                 ans1(idmn),ans2(lena2)
      integer i,j,k,idmn,idmnsq,lena2,ipos,info,ipvt(idmn)

c Obtain D_d^+*a_2:

      do 10 j = 1,idmnsq
         wk1(j) = 0.0
         do 20 i = 1,lena2
            wk1(j) = wk1(j) + DdPlus(i,j)*a2(i)
20       continue
10    continue
      
c Obtain A_2:

      ipos = 1
      do 30 j = 1,idmn
         do 40 i = 1,idmn
             A2mat(i,j) = wk1(ipos)
             A2str(i,j) = wk1(ipos)
             ipos = ipos + 1
40       continue
30    continue

c Store a_1 and c_1 in A2ina1 and A2inc1:

      do 50 i = 1,idmn
         A2ina1(i) = a1(i)
         A2inc1(i) = c1(i)
50    continue

c Obtain A_2^(-1)c_1 and A_2^(-1)a1 (noting storage
c of solution vectors in `c1' and `a1':

      call dgefa(A2mat,idmn,idmn,ipvt,info)
      call dgesl(A2mat,idmn,idmn,ipvt,A2inc1,0)
      call dgesl(A2mat,idmn,idmn,ipvt,A2ina1,0)

c Obtain c_1^TA_2^(-1)c_1 and c1_T^A_2^(-1)a1:

      c1A2c1 = 0.0
      c1A2a1 = 0.0
      do 60 i = 1,idmn
         c1A2c1 = c1A2c1 + c1(i)*A2inc1(i)
         c1A2a1 = c1A2a1 + c1(i)*A2ina1(i)
60    continue

c Obtain r_1 and r_2:
  
      r1 = sqrt(2.0*(2.0 - c1A2c1))
      r2 = (2.0*c0 - c1A2a1)/r1

c Obtain zetdv and zetddv:

      call zetad(r2,zetdv)
      zetddv = -zetdv*(r2 + zetdv) 

c Obtain r_3 and r_4:

      r3 = 2.0*zetdv/r1
      r4 = -2.0*zetddv/(r1*r1)

c Obtain A_2 + r_4*c1*c1^T:

      do 70 i = 1,idmn
         do 80 j = 1,idmn
            R2comp(i,j) = A2str(i,j) + r4*c1(i)*c1(j)
80       continue
70    continue

c Obtain R_5:

      call dgefa(R2comp,idmn,idmn,ipvt,info)
      do 90 j = 1,idmn
         do 100 i = 1,idmn
            wk2(i) = A2str(i,j)
100      continue
         call dgesl(R2comp,idmn,idmn,ipvt,wk2,0.0)
         do 110 i = 1,idmn
            R5(i,j) = wk2(i)
110      continue        
90    continue

c Obtain ans1 = R_5^T*(a_1 + r_3*c1):

      do 120 i = 1,idmn
         ans1(i) = 0.0
         do 130 j = 1,idmn
            ans1(i) = ans1(i) + R5(j,i)*(a1(j) + r3*c1(j))
130      continue
120   continue

c Obtain R_5^T*A_2:

      do 140 i = 1,idmn
         do 150 j = 1,idmn
            R5TA2(i,j) = 0.0
            do 160 k = 1,idmn
               R5TA2(i,j) = R5TA2(i,j) + R5(k,i)*A2str(k,j)
160         continue       
150      continue
140   continue

c Obtain vec(R_5^T*A_2):

      ipos = 1
      do 170 j = 1,idmn
         do 180 i = 1,idmn
            vR5TA2(ipos) = R5TA2(i,j)
            ipos = ipos + 1
180      continue
170   continue

c Obtain D_d^T vec(R_5^T*A_2):

      do 190 i = 1,lena2
         ans2(i) = 0.0
         do 200 j = 1,idmnsq
            ans2(i) = ans2(i) + Dd(j,i)*vR5TA2(j)   
200      continue
190   continue

      return
      end

cccccccccc End of kpbt.f cccccccccc

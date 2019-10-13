cccccccccc FORTRAN subroutine asn.f cccccccccc

c For computation of the log-partition function
c of the Multivariate Normal distribution.

c Last changed: 17 SEP 2019

      subroutine asn(a1,a2,A2ina1,idmn,idmnsq,lena2,A2mat,xm2A2,DdPlus,
     +               wkv,ipvt,det,work,ans)     
      double precision a1(idmn),a2(lena2),A2ina1(idmn),wkv(idmnsq),
     +                 DdPlus(lena2,idmnsq),A2mat(idmn,idmn),
     +                 xm2A2(idmn,idmn),det(2),work(idmn,idmn),xldm2A,
     +                 a1A2a1,ans
      integer i,j,idmn,idmnsq,lena2,ipos,info,ipvt(idmn)

c Obtain D_d^+*a_2:

      do 10 j = 1,idmnsq
         wkv(j) = 0.0
         do 20 i = 1,lena2
            wkv(j) = wkv(j) + DdPlus(i,j)*a2(i)
20       continue
10    continue

c Obtain A_2:

      ipos = 1
      do 30 j = 1,idmn
         do 40 i = 1,idmn
             A2mat(i,j) = wkv(ipos)
             xm2A2(i,j) = -2.0*wkv(ipos)
             ipos = ipos + 1
40       continue
30    continue

c Store a_1 in A2ina1:

      do 50 i = 1,idmn
         A2ina1(i) = a1(i)
50    continue

c Obtain A_2^(-1)a1 (noting storage
c of solution vector `a1'):

      call dgefa(A2mat,idmn,idmn,ipvt,info)
      call dgesl(A2mat,idmn,idmn,ipvt,A2ina1,0)

c Obtain a1_T^A_2^(-1)a1:

      a1A2a1 = 0.0
      do 60 i = 1,idmn
         a1A2a1 = a1A2a1 + a1(i)*A2ina1(i)
60    continue

c Obtain log|-2*A2|:

      call logdet(xm2A2,idmn,ipvt,work,det,xldm2A) 

c Return answer:

      ans = -0.25*a1A2a1 - 0.5*xldm2A

      return
      end

cccccccccc End of asn.f cccccccccc

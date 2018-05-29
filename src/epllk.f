cccccccccc FORTRAN subroutine epllk.f cccccccccc

c For computing the expectation propagation to the 
c log-likelihood for probit mixed models.

c Last changed: 24 MAY 2018

      subroutine epllk(beta,etaSg0,etaSg,m,nVec,nMax,numObs,indStt,
     +                 idF,idR,idRsq,lena,lena2,nlena,yDagg,Xf,Xr,
     +                 yDcur,XfCur,XrCur,uHat,Dd,DdPlus,
     +                 Xbeta,XuHat,etaFtS,etaStF,SUMlt,c1Cur,etaIN1,
     +                 etaIN2,eINa1,eINa2,eINb1,eINb2,etaPvF,etaPvS,
     +                 etaCrF,etaCrS,wk1,ipvt,A2ina1,A2inc1,A2mat,A2str,
     +                 R2comp,wk2,R5,R5TA2,vR5TA2,xkans1,xkans2,
     +                 B2inb1,work,A2neg,B2mat,B2str,B2neg,
     +                 xm2A2,det,wka,wkb,wkv,xmiscl,etaOut)
      double precision beta(idF),etaSg0,etaSg(lena),yDagg(numObs),
     +                 Xf(numObs,idF),Xr(numObs,idR),yDcur(nMax),
     +                 XfCur(nMax,idF),XrCur(nMax,idR),uHat(m,idR),
     +                 Dd(idRsq,lena2),DdPlus(lena2,idRsq),
     +                 Xbeta(nMax),XuHat(nMax),pi,aHat,zetdf,zetddf,
     +                 etaFtS(nMax,lena),etaStF(nMax,lena),SUMlt(lena),
     +                 c0Cur,c1Cur(idR),etaIN1(idR),etaIN2(lena2),
     +                 relErr,eINa1(idR),eINa2(lena2),eINb1(idR),
     +                 SUMlt0,eINb2(lena2),etaPvF(nlena),etaPvS(nlena),
     +                 etaCrF(nlena),etaCrS(nlena),wk1(idRsq),errCur,
     +                 A2ina1(idR),A2inc1(idR),A2mat(idR,idR),EPrltl,
     +                 A2str(idR,idR),R2comp(idR,idR),wk2(idR),
     +                 R5(idR,idR),R5TA2(idR,idR),vR5TA2(idRsq),
     +                 xkans1(idR),xkans2(lena2),wka(idRsq),wkb(idRsq),
     +                 b1str(idR),det(2),work(idR,idR),A2neg(idR,idR),
     +                 B2mat(idR,idR),B2str(idR,idR),B2neg(idR,idR),
     +                 xm2A2(idR,idR),wkv(idRsq),AsNans,eFtS0,xmiscl(3),
     +                 etaOut(m,lena)
      integer m,nVec(m),nMax,numObs,indStt(m),idR,idRsq,lena,lena2,
     +        nlena,idF,idFwk,i,j,k,kd,ipvt(idR),iPos,icnvgd,itNum,
     +        iEPmax


c     Set the expectation propagation control parameters:

      iEPmax = xmiscl(1)
      EPrltl = xmiscl(2)

c     Set 'idF' working variable:

      idFwk = idF
      if (idF.eq.0) then
         idFwk = 1
      endif
      
c     Set value of pi:

      pi = 4.0*atan(1.0)

c etaFtS: eta$"p(y|u;beta0,beta,x)->u"
c etaStF: eta$"u->p(y|u;beta0,beta,x)"

c     Add up each of the m terms in the
c     EP-approximate log-likelihood:

      xmiscl(3) = 0.0
      do 10 i=1,m    

c        Obtain XfCur, XrCur, Xbeta(j) and XuHat(j):

         do 20 j = 1,nVec(i)        
            yDCur(j) = yDagg(indStt(i)+j-1) 
            do 30 k = 1,idFwk 
               XfCur(j,k) = Xf(indStt(i)+j-1,k)
30          continue
            do 40 k = 1,idR
               XrCur(j,k) = Xr(indStt(i)+j-1,k)
40          continue
            Xbeta(j) = 0.0
            do 50 k = 1,idFwk
               Xbeta(j) = Xbeta(j) + XfCur(j,k)*beta(k)
50          continue
            XuHat(j) = 0.0
            do 60 k = 1,idR
               XuHat(j) = XuHat(j) + XrCur(j,k)*uHat(i,k)
60          continue
20       continue

c        Obtain starting values:

         do 70 j = 1,nVec(i)
            
            aHat = yDcur(j)*(Xbeta(j) + XuHat(j))
            
            call zetad(aHat,zetdf)
            zetddf = -zetdf*(aHat + zetdf)

            do 80 k = 1,idR
               etaFtS(j,k) = yDCur(j)*zetdf*XrCur(j,k)
     +                       - zetddf*XrCur(j,k)*XuHat(j) 
80          continue

            iPos = idR
            do 90 k = 1,idR
               do 100 kd = k,idR
                  iPos = iPos + 1
                  etaFtS(j,iPos) = 0.5*zetddf*XrCur(j,k)*XrCur(j,kd)
100            continue
90          continue
70       continue

c        Now obtain ith term of the EP-approximate log-likelihood:                

         icnvgd = 0
         itNum = 0
c
c           Top of iteration loop.
c
110         itNum = itNum + 1  

            do 120 k = 1,lena   
               SUMlt(k) = 0.0
               do 130 j = 1,nVec(i)
                  SUMlt(k) = SUMlt(k) + etaFtS(j,k)
130            continue
120         continue

            do 140 j = 1,nVec(i)

               do 150 k = 1,lena
                  etaStF(j,k) = etaSg(k) + SUMlt(k) - etaFtS(j,k)
150            continue

               c0Cur = yDcur(j)*Xbeta(j)
               do 160 k = 1,idR
                  c1Cur(k) = yDcur(j)*XrCur(j,k)
160            continue
                
               do 170 k = 1,idR
                  etaIN1(k) = etaStF(j,k)
170            continue

               do 180 k = 1,lena2
                  etaIN2(k) = etaStF(j,idR+k)
180            continue

               call kpbt(etaIN1,etaIN2,c0Cur,c1Cur,idR,idRsq,lena2,
     +                    Dd,DdPlus,wk1,A2ina1,A2inc1,ipvt,A2mat,A2str,
     +                    R2comp,wk2,R5,R5TA2,vR5TA2,xkans1,xkans2)    

               do 190 k = 1,idR
                  etaFtS(j,k) = xkans1(k) - etaStF(j,k)
190            continue

               do 200 k = 1,lena2
                  etaFtS(j,idR+k) = xkans2(k) - etaStF(j,idR+k)
200            continue       
140         continue     

            if (itNum.gt.1) then
          
c              Do set up for convergence check:

               iPos = 0
               do 210 j = 1,nVec(i)
                  do 220 k = 1,lena
                     iPos = iPos + 1
                     etaCrF(iPos) = etaFtS(j,k)
                     etaCrS(iPos) = etaStF(j,k)
220               continue
210            continue             

               relErr = -1.0
               do 230 k = 1,nlena
                  errCur = abs((etaCrF(k) - etaPvF(k))/etaPvF(k))
                  if (errCur.gt.relErr) then
                     relErr = errCur
                  endif
                  errCur = abs((etaCrS(k) - etaPvS(k))/etaPvS(k))
                  if (errCur.gt.relErr) then
                     relErr = errCur
                  endif
230            continue                  

c              Check for convergence:

               if (relErr.lt.EPrltl.or.itNum.ge.iEPmax) then
                  icnvgd = 1
               endif     
            endif

c           Set `etaPvF' and `etaPvS' (previous iteration) values:

            do 240 k = 1,nlena
               etaPvF(k) = etaCrF(k)
               etaPvS(k) = etaCrS(k)
240         continue            

            if (icnvgd.eq.0) then 
               goto 110
            endif
c
c Bottom of iteration loop.
c 

         if (itNum.ge.iEPmax) then
            call intpr("",0,itNum,0)
            call intpr("Warning: the number of expectation",34,itNum,0) 
            call intpr("propagation iterations reached",30,itNum,0)
            call intpr("its maximum value (EPmaxit).",28,itNum,0) 
            call intpr("",0,itNum,0)
         endif

c        Update the zeroth entries of the eta vectors:

         SUMlt0 = 0.0
         do 250 j = 1,nVec(i)
            do 260 k = 1,idR
               eINa1(k) = etaStF(j,k)
               eINb1(k) = etaFtS(j,k) + etaStF(j,k)
260         continue

            do 270 k = 1,lena2
               eINa2(k) = etaStF(j,idR+k)
               eINb2(k) = etaFtS(j,idR+k) + etaStF(j,idR+k)
270         continue

            c0Cur = yDcur(j)*Xbeta(j)
            do 280 k = 1,idR
                  c1Cur(k) = yDcur(j)*XrCur(j,k)
280         continue
           
            call cpbt(eINa1,eINa2,eINb1,eINb2,c0Cur,c1Cur,idR,idRsq,
     +                lena2,Dd,DdPlus,wka,wkb,A2ina1,B2inb1,A2inc1,ipvt,
     +                det,work,A2mat,A2neg,B2mat,B2neg,eFtS0) 
                       
            SUMlt0 = SUMlt0 + eFtS0

250      continue

c        Obtain the ith term of the approximate log-likelihood
c        and update:

         do 290 k = 1,idR
            etaIN1(k) = etaSg(k) + SUMlt(k)
            etaOut(i,k) = etaIN1(k)
290      continue

         do 300 k = 1,lena2
            etaIN2(k) = etaSg(idR+k) + SUMlt(idR+k)
            etaOut(i,idR+k) = etaIN2(k)
300      continue

         call asn(etaIN1,etaIN2,A2ina1,idR,idRsq,lena2,A2mat,xm2A2,
     +            DdPlus,wkv,ipvt,det,work,AsNans)

         xmiscl(3) = xmiscl(3) + 0.5*idR*log(2*pi) 
     +              + etaSg0 + SUMlt0 + AsNans    

10    continue

      return
      end

cccccccccc End of epllk.f cccccccccc

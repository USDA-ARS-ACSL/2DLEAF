*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      subroutine Break()
      include 'public.ins'
      Dimension SumFlux(3),S(3),iCav(NumELD)
      Common /Brake/ NumVL, Dx,Dz1,Dz2,Dair(10),jPrint(10),RA,BMult
         If(lInput.eq.0) goto 11
         Open(10,file='Param_br.dat')
         Read(10,*)
         Read(10,*)
         Read(10,*) NumVL, Dx, Dz1, Dz2, IAX
         Close(10)
         Open(10,file='Param_g.dat')
         Do i=1,8
           Read(10,*)
         Enddo
         Read(10,*) (Dair(jjj),jjj=1,NumG)
         Read(10,*)
         Read(10,*) BMult
         Close(10)
         Do jjj=1,NumG
           Dair(jjj)=Dair(jjj)*BMult
         Enddo
         Open(10,file='Output.dat')
         Do i=1,6
           Read(10,*)
         Enddo
         Read(10,*) (jPrint(i),i=1,NumG)
         Close(10)
         Open(18,file='g____.grp',status='unknown')
         Open(19,file='fluxnod.out',status='unknown')
         Open(20,file='fluxtot.out',status='unknown')
            RA=1.E+05
         Write(20,131) 
         Write(52,*)"Calculating..."
c         CALL WCOU@('Iter   ',RED+INVERSE_GREEN,WINDOW(3))
         Write(53,*)"Iter   "
         Iter=0       
         WRITE(LINE,100) Iter
c         CALL WCOU@(LINE,RED+INVERSE_GREEN,WINDOW(4))
         Write(54,*) Line
         lInput=0
         Return
11        A=0.
          Do jjj=1,NumG
          If(jPrint(jjj).ne.0) then
            s(jjj)=0.
            Do i=1,NumVL
               W=Dx
               if(i.eq.1.or.i.eq.NumVL) W=W/2
               s(jjj)=s(jjj)+Dair(jjj)*W*(g(i,jjj)-g(i+NumVL,jjj))/Dz1
            Enddo
CYAP 07/12
            If(IAX.eq.2) then
              Do i=1,NumVL
               W=Dx
               if(i.eq.1.or.i.eq.NumVL) W=W/2
               s(jjj)=
     & s(jjj)+Dair(jjj)*W*
     & (g(NumNP-(i-1),jjj)-g(NumNP-(i-1)-NumVL,jjj))/Dz2
              Enddo
            Endif
            SumFlux(jjj)=s(jjj)
            s(jjj)=0.
            Do i=1,NumBP
              n=KXB(i)
              k=CodeG(n)
              if(k.eq.-4) then
                If(jjj.ne.3) then 
                  W=Width(i)
                Else
                  W=1.
                Endif
                s(jjj)=s(jjj)-VarBG(i,jjj,3)*W
              Endif
            Enddo
            A=AMAX1(A,abs(SumFlux(jjj)/s(jjj)-1.))
          Endif
        Enddo
        A=abs(RA/s(1)-1.)
        RA=s(1)
        Write(20,132) Iter, s(1)
        If(A.lt.1.E-03.or.Iter.gt.20) then
          Do jjj=1,NumG
            If(jPrint(jjj).ne.0) then
              write(18,130) jjj
              Do i=1,NumNP
                write(18,140) x(i),y(i),g(i,jjj)
              Enddo
              write(19,130) jjj
              Write(19,150)
              Do i=1,NumBP
                n=KXB(i)
                k=CodeG(n)
                if(k.eq.-4) then
                  Write(19,160) n, -VarBG(i,jjj,3),Width(i),g(n,jjj)
                endif
              Enddo
              If(jjj.eq.1) then
                Write(20,170) jjj, s(jjj)
              else
                Write(20,170)  jjj, SumFlux(jjj)
              Endif
            Endif
          Enddo
c Here goes average concentration calculation
         Open(73,file='aver.dat')
         Read(73,'(30I2)') (iCav(n),n=1,NumEl)
         Read(73,*) iFirst,iLast
         Close(73)
         Open(21,file='avercon.out')
         Do jjj=1,3
         If(jPrint(jjj).ne.0) then
         SumC=0.
         SumA=0.
         Do n=1,NumEl
           If(iCav(n).eq.1) then
             NUS=4
             If(KX(n,3).eq.KX(n,4)) NUS=3
c                        Loop on subelements
             Do k=1,NUS-2
               i=KX(N,1)
               j=KX(N,k+1)
               l=KX(N,k+2)
               Ci=x(l)-x(j)
               Cj=x(i)-x(l)
               Ck=x(j)-x(i)
               Bi=y(j)-y(l)
               Bj=y(l)-y(i)
               Bk=y(i)-y(j)
               AE=(Ck*Bj-Cj*Bk)/2.
               gE=(g(i,jjj)+g(j,jjj)+g(l,jjj))/3.
               SumC=SumC+gE*AE
               SumA=SumA+AE
             Enddo
           Endif
         Enddo
         Write(21,180) jjj,SumC/SumA
         Endif
         Enddo
c        The whole 'between-box space
         Do jjj=1,3
         If(jPrint(jjj).ne.0) then
         SumC=0.
         SumA=0.
         Do n=iFirst,iLast
             NUS=4
             If(KX(n,3).eq.KX(n,4)) NUS=3
c                        Loop on subelements
             Do k=1,NUS-2
               i=KX(N,1)
               j=KX(N,k+1)
               l=KX(N,k+2)
               Ci=x(l)-x(j)
               Cj=x(i)-x(l)
               Ck=x(j)-x(i)
               Bi=y(j)-y(l)
               Bj=y(l)-y(i)
               Bk=y(i)-y(j)
               AE=(Ck*Bj-Cj*Bk)/2.
               gE=(g(i,jjj)+g(j,jjj)+g(j,jjj))/3.
               SumC=SumC+gE*AE
               SumA=SumA+AE
             Enddo
         Enddo
         Write(21,181) jjj,SumC/SumA
         Endif
         Enddo
c Surface
         SumC=0.
         SumA=0.
         Do i=1,NumBP
           n=KXB(i)
           k=CodeG(n)
           if(k.eq.-4) then
             SumC=SumC+g(n,1)
             SumA=SumA+1.
           endif
         Enddo
         Write(21,182) SumC/SumA
         Close(21)
c         call TEXT_MODE@
         stop
        else
          Iter=Iter+1
          WRITE(LINE,100) Iter
          write(54,*) Line
c          CALL WCOU@(LINE,RED+INVERSE_GREEN,WINDOW(4))
        endif
100    FORMAT(1x,I12)
130   format(//
     !         '      x              z              Gas # ',I3,/)
131   format('Iter      SumFlux')
132   format(I3,E14.4)
140   format(1x,3G17.6)
150   format(//,'Node      Flux            Width          Conc')
160   format(i5,3E17.6)
170   format('   Gas # ',I3,'   Total flux = ',E17.6)
180   format('   Gas # ',I3,' Average concentration in the',
     &      ' stomatal cavity = ', E13.4)
181   format('   Gas # ',I3,' Average concentration in the',
     &      ' between-box space = ', E13.4)
182   format('   Gas # 1       Average concentration on the',
     &      ' cell surface = ', E13.4)
      return
      end         


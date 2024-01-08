c||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      subroutine GasMover ()
      include 'public.ins'

      double precision A,B,C,P,Sum
      dimension A(2*MBandD-1,NumNPD),B(NumNPD),Qc(NumNPD),
     !                S(3,3),Ci(3),Bi(3),List(3),
     !                Diff(NumNPD),
     &                ALP(2*MBandD-1,NumNPD),BET(NumNPD)
      Integer gKod         
      Common /GasM/ NLevel,epsi,Dair(NumGD),BMult
        

      if (lInput.eq.0) goto 11
         
* --- Reading of the input files and initial calculations --------------
      im=130
      il=0       
      open(40,file='Param_G.dat', status='old',ERR=10)   
      write(52,*)"GasMover: computational parameters"
c      CALL WCOU@('GasMover: computational parameters',
c     !            RED+INVERSE_GREEN,WINDOW(2))
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10) NumG
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10) epsi
      Nlevel=1
      if(epsi.lt.0.999) NLevel=2      
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10) (Dair(jjj),jjj=1,NumG)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10) BMult
      close(40)
*
      im=155
      il=0
      open(41,file='Nodal_G.dat', status='old',ERR=10)
      Write(52,*) "Gas mover: nodal information"
c      CALL WCOU@('Gas mover: nodal information',
c    !            RED+INVERSE_GREEN,WINDOW(2))
      im=im+1
      il=il+1
      read(41,*,ERR=10)
      im=im+1
      il=il+1
      read(41,*,ERR=10)
      im=im+1
         do i=1,NumNP
      il=il+1
      read(41,*,ERR=10) n,(g(n,jjj),jjj=1,NumG)
         enddo
      close(41)
      return
         
* --- Beginning of calculations --------------------------------------

11    continue 

*     Calculate gas flow

      do 711 jjj=1,NumG 
      IF(Iter.gt.0. and. JJJ.eq.3) goto 711
      IF(JJJ.NE.2) THEN
      Do i=1,NumNP
        Diff(i)=Dair(jjj)
        If(LocInd(i).eq.0) Diff(i)=Diff(i)*BMult
      Enddo
      xMul=1.
      alf=1.-epsi
c        
      do 13 i=1,NumNP
        B(i) =0.
        Qc(i)=0.
        BET(i)=0.
          do 12 j=1,2*MBand-1
            A(j,i)=0.
            ALP(j,i)=0.
12        continue
13    continue
c==== Assembling of the matrixes
c
*       Loop on elements
        do 19 n=1,NumEl
          NUS=4
          if(KX(n,3).eq.KX(n,4)) NUS=3
*         Loop on subelements
          do 18 k=1,NUS-2
            i=KX(n,1)
            j=KX(n,k+1)
            l=KX(n,k+2)
            List(1)=i
            List(2)=j
            List(3)=l
            Ci(1)=x(l)-x(j)
            Ci(2)=x(i)-x(l)
            Ci(3)=x(j)-x(i)
            Bi(1)=y(j)-y(l)
            Bi(2)=y(l)-y(i)
            Bi(3)=y(i)-y(j)
            AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
            Ec=(Diff(i)+Diff(j)+Diff(l))/3.
            FMul=xMul*AE/4.
            SMul1=-1./AE/4.*xMul
            do 17 j1=1,3
              i1=List(j1)
              do 16 j2=1,3
                i2=List(j2)
                S(j1,j2)=SMul1*Ec*(Bi(j1)*Bi(j2)+Ci(j1)*Ci(j2))
                  iB=MBand+i2-i1
                  A(iB,i1)=A(iB,i1)+S(j1,j2)
                  Qc(i1)=Qc(i1)-S(j1,j2)*g(i2,jjj)
16            continue
17          continue
18        continue
19      continue
        do 21 i=1,NumNP
         If(JJJ.ne.3) then
          B(i)  =Qc(i)
         Endif
          do 20 j=1,MBandD
            ALP(j,i)=A(j,i)
20        continue
21    continue
c     Boundary conditions
      do 114 i=1,NumNP
      if(CodeG(i).ne.0) then
        do 111 j=1,NumBP
          if(KXB(j).eq.i) then
CYAP
             if(CodeG(i).gt.0.OR.JJJ.eq.3) then
                gKod=1
                if(CodeG(i).eq.1.or.jjj.eq.3) then 
                   gBnd=g(i,jjj)
                   VarBG(j,jjj,3)=Qc(i)
                else
                   gBnd=VarBG(j,jjj,1)
                endif
             else
                  gKod=3
                  Qg=VarBG(j,jjj,3)*Width(j)
                  Pg=VarBg(j,jjj,2)*Width(j)
             endif
             goto 112
          endif
111     continue

*     Dirichlet boundary condition
112       if(gKod.eq.1) then
            do 113 j=1,2*MBand-1
              A(j,i)=0.
CYAP??
              ALP(j,i)=0.
113         continue
            A(MBand,i)=1.
            B(i)=gBnd
            ALP(MBand,i)=1.
            BET(i)=0.
          end if
*     Cauchy boundary condition
          if(gKod.eq.3) then
               BET(i)=B(i)-Qg
               ALP(MBand,i)=ALP(MBand,i)-Pg
          end if
        end if
114    continue
      If(JJJ.ne.3) then
        Do i=1,NumNPD
          B(i)=BET(i)
        Do j=1,MBandD
          A(j,i)=ALP(j,i)
        Enddo
        Enddo
      Endif
c     Solve the global matrix equation for transport
c     Implicit formulation
      N1=NumNP-1
      do 212 k=1,N1
        P=1./A(MBand,k)
        kk=k+1
        kc=MBand
        do 211 i=kk,NumNP
          kc=kc-1
          if(kc.le.0) goto 212
          C=-P*A(kc,i)
          A(kc,i)=C
          ii=kc+1
          L=kc+MBand-1
          do 211 j=ii,L
            jj=j+MBand-kc
            A(j,i)=A(j,i)+C*A(jj,k)
211      continue
212    continue
      do 214 i=2,NumNP
        jj=MBand+1-i
        ii=1
        if(jj.le.0) then
          jj=1
          ii=i-MBand+1
        end if
        Sum=0.
        do 213 j=jj,MBand-1
          Sum=Sum+A(j,i)*B(ii)
          ii=ii+1
213      continue
        B(i)=B(i)+Sum
214    continue              
      B(NumNP)=B(NumNP)/A(MBand,NumNP)
      do 216 k=1,N1
        i=NumNP-k
        jj=i
        m=min0(2*MBand-1,MBand+k)
        Sum=0.
        do 215 j=MBand+1,m
          jj=jj+1
          Sum=Sum+A(j,i)*B(jj)
215      continue
        B(i)=(B(i)-Sum)/A(MBand,i)
216    continue  
*Final assignment
      do i=1,NumNP
       If(JJJ.ne.3)then
        g(i,jjj)=sngl(B(i))+g(i,jjj)
       Else
        g(i,jjj)=sngl(B(I))
*        if(g(i,jjj).lt.1.0E-30) g(i,jjj)=0.
       Endif
      enddo      
*       
      ENDIF
711   continue
*
      return
10    call errmes(im,il)
      return
      end


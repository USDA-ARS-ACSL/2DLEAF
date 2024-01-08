       Subroutine SetSurface()
       Include 'Public.ins'
       Real Kc,Ko,IPPFD
       COMMON /SurfPar/ Vmax,Tau,Rd,Kc,Ko,IPPFD,Pml,Alpha
       If (lInput.eq.0) goto 11
       im=160
       il=0
       Open (5,file='SetSurf.dat',status='old',ERR=10)
        im=im+1
        il=il+1
        READ (5,*,ERR=10)
        im=im+1
        il=il+1
        READ (5,*,ERR=10)
        im=im+1
        il=il+1
        READ (5,*,ERR=10) Vmax,Tau,Rd,Kc,Ko
        im=im+1
        il=il+1
        READ (5,*,ERR=10)
        im=im+1
        il=il+1
        READ (5,*,ERR=10) Pml, Alpha, IPPFD
c..................... Routine calculations
   11  Continue
C==== March along the boundary
      Do i=1,NumBP
      n=KXB(i)
      k=CodeG(n)
      if(k.eq.-4) then
C     Local concentrations
      C_CO2=g(n,1)
      C_O2 =g(n,2)
C
C Calculating CO2 - saturated rate of photosynthesis and CO2 compensation
C point
C
      z=Alpha*IPPFD
      Pm=z/sqrt(1+(z/Pml)**2)
      Gamma= C_O2/Tau/2.
C
C Calculating two possible rates of CO2 consumption
C
C
C Electron transport limiting
C
      Wj=Pm/(1.+2.*Gamma/C_CO2)
C Rubisco activity limiting
C
      Wc=Vmax/(1.+(Kc/C_CO2)*(1.+C_O2/Ko))
C
C True rate is the minimal one
C
      if(Wj.lt.Wc) then
         VarBG(i,1,2) = 3.* Gamma*Pm/(C_CO2+Gamma)**2.
         VarBG(i,1,3) =-(1.-Gamma/C_CO2)*Wj+Rd
      else
         VarBG(i,1,2) = Vmax*(Kc*(1+C_O2/Ko)+Gamma)
     &                  /(C_CO2+Kc*(1+C_O2/Ko))**2.
         VarBG(i,1,3) =-(1.-Gamma/C_CO2)*Wc+Rd
      endif
C Oxygen is the same for now
      VarBG(i,2,3)=0.
      Endif
      Enddo
      Return
10    call errmes(im,il)
      END


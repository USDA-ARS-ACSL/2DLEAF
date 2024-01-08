C This program calculates input parameters for the Harley-Tenhunen
C model using dependencies on temperature and cell area index
C Everything here is in mol, m, s, J
C Fig 2-4 and Table 2.1 are used
C 13.81=ln(10^6). 6.91=ln(10^3)
      Real Kc, Ko, IPPFD
      Open(10,File='p.dat')
      Read(10,*)
      Read(10,*)
      Read(10,*)
      Read(10,*) T,C30,IPPFD
      Close(10)
      Open(11,File='CAI.dat')
      Read(11,*)
      Read(11,*) Dummy,CAI
      Close(11)
      Pml =1.3*T1(29.20-13.81, 62.6E+03, 0.637E+03, 200E+03, T)/CAI
      Vmax=3.3*T1(34.63-13.81, 74.6E+03, 0.635E+03, 200E+03, T)/CAI
      Tau =T2(-3.95,      -29.0E+03,                     T)
      Kc  =T2(31.95-13.81, 65.0E+03,                     T)
      Ko  =T2(19.61 -6.91, 36.0E+03,                     T)
      Alpha=0.009
C What instead?
      Rd=0.025*Vmax
C
C HT have Kc and Ko in mkL L-1. We need them in mol L-1.
C At this point we have them in m^3 m-3 or in bars.
C Since M=pv/(RT) and p=10^Pa,
      B=1E+05/(8.31*T)
      Kc=Kc*B
      Ko=Ko*B
C
      Open(11,File='Setsurf.dat')
      Write(11,100)
      Write(11,200)
      Write(11,201) Vmax,Tau,Rd,Kc,Ko
      Write(11,300)
      Write(11,301) Pml,Alpha,IPPFD
      Close(11)
100   FORMAT('Data for the HT Setsurf module. Units are: mol, m, s')
200   FORMAT('Vmax       Tau         Rd         Kc         Ko')
201   FORMAT(5E12.4)
300   FORMAT('Pml        Alpha       IPPFD')
301   FORMAT(3E12.4)
      Stop
      End
      Function T1(c,DeltaH_a,DeltaS,DeltaH_d,T)
      U=8.31*T
      T1=exp(c-DeltaH_a/U)/(1.+exp((DeltaS*T-DeltaH_d)/U))
      Return
      End
      Function T2(c,DeltaH_a,T)
      U=8.31*T
      T2=exp(c-DeltaH_a/U)
      Return
      End


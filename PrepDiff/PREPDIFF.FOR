C This program calculates diffusion coefficients of CO2, O2 and
C water vapor in the air and puts the values into file PARAM_G.dat
      Character*80 LIN
      Open(10,file='param_g.dat',status='OLD',ERR=300)
      Open(12,file='paratemp.dat')
      Do i=1,8
        Read(10,'(A80)') LIN
        Write(12,'(A80)') LIN
      Enddo
      Open(11,file='p.dat')
      Read(11,*) 
      Read(11,*)
      Read(11,*)
      Read(11,*) T
      DCO2 =0.139E-04
      D_O2  =0.634E-04
      DH2O =0.239E-04
      TCORR=(T/273.15)**2
      Write(12,'(3E15.4)') DCO2*TCORR,D_O2*TCORR,DH2O*TCORR
      Read(11,*)
      Read(11,*) BMult
      Write(12,*)
      Write(12,*) BMult
      Stop
300   Write(*,*) 'No param_g.dat file'
      End


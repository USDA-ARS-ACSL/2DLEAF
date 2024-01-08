      Real IPPFD
      Dimension x(6500),y(6500),LocInd(6500)
      Character*80 A
      Open(10,file='grid_bnd.dat')
      Open(12,file='grid_new.dat')
      Open(11,file='p.dat')
      Open(17,file='prepgrid.dat')
      Read(17,*)
      Read(17,*)
      Read(17,*) NumGl, NumVl, dx, dy, IAX
      If(IAX.eq.2) then
        Read(11,*)
        Read(11,*) d_top,ylid_top,d_bottom,ylid_bottom
      Else
        Read(11,*)
        Read(11,*) d_top,ylid_top
      Endif
	 Do i=1,2
        Read(10,'(A80)') A
        Write(12,'(A80)') A
      Enddo
      Read(10,*)  NumNP, NumEl, NumBP, NumVL, IAX
      Write(12,'(i6,4i10)')  NumNP, NumEl, NumBP, NumVL,IAX
      Read(10,'(A80)') A
      Write(12,'(A80)') A
      Do i=1,NumNP
        Read(10,*) n,x(i),y(i),LocInd(i)
      Enddo
      d_top_old=y(1)-ylid_top
      If(IAX.eq.2)  d_bottom_old=ylid_bottom-y(NumNP)
      Do i=1,NumNP
        y_new=y(i)
        If(y(i).GE.ylid_top) 
     &     y_new=ylid_top+(d_top/d_top_old)*(y(i)-ylid_top)
        If(IAX.eq.2) then
          If(y(i).LE.ylid_bottom) 
     &     y_new=ylid_bottom+(d_bottom/d_bottom_old)*(y(i)-ylid_bottom)
        Endif
        Write(12,'(i4,2E13.3,i3)') i,x(i),y_new, LocInd(i)
      Enddo
1     Read(10,'(A80)',ERR=300,END=300) A
      Write(12,'(A80)') A
      Goto 1
300   Close(10)
      Close(12)
      Dy1=dy*(d_top/d_top_old)
      If(IAX.eq.2) then
        Dy2=dy*(d_bottom/d_bottom_old)
      Else
        Dy2=0.
      Endif
      Open(14,file='param_br')
      Write(14,17)
17    Format('Break data',/,'NumVL   Dx    Dy1     Dy2    IAX')
      Write(14,*) NumVl,Dx,Dy1,Dy2,IAX
      Close(14)
c
      Read(11,*)
      Read(11,*) T,C30,IPPFD
      Open(15,File='Nodal_g.dat')
      Open(16,File='Nodal_g')
      Read(15,*)
      Write(16,*)
      Read(15,*)
      Write(16,*)
      Do i=1,5000
        Read(15,120,Err=70,End=70) k,u1,u2,u3
        If(k.le.NumVl) u3=C30
        if(k.gt.NumNP-NumVl.and.IAX.EQ.2) u3=C30
        Write(16,120) k,u1,u2,u3
      Enddo
70    Close(11)
      Close(15)
      Close(16)
120   Format(I4,3E11.4)
      Stop
      End

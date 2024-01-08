      Parameter (NumGlD=301,NumVLD=80,NumNPD=NumGLD*NumVLD,
     &           NumClD=40,NumBPD=3000, NumGD=3)
c This program generates rectangular grid data for stomata2.for
      dimension x(NumVLD),y(NumGLD),xN(NumNPD),yN(NumNPD),
     &          xLU(NumCLD),yLU(NumCLD),xRB(NumCLD),yRB(NumCLD),
     &          xLUC(NumCLD),yLUC(NumCLD),xRBC(NumCLD),yRBC(NumCLD),
     &          iCav(NumNPD),
     &          MatNum(NumCLD),
     &          KX(NumNPD,4),U(4),
     &          KodeG(NumBPD),Width(NumBPD),KXB(NumBPD),
     &          gInit(NumGD),NodNum(NumGLD,NumVLD)
      integer e
      open(8,file='prepgrid.dat')
      read(8,*)
      read(8,*)
      read(8,*) NumGl,NumVL,dx,dy,IAX
      edx=dx/3.
      edy=dy/3.
c NumGL - number of transverse lines in the rectangular part
c NumVL - number of vertical lines in the rectangular part
      Do j=1,NumVL
        x(j)=float(j-1)*dx
      Enddo
c x coordinates of the vertical grid lines from the left to the right
      Do i=1,NumGL
        y(i)=float(NumGL-i)*dy
      Enddo
c y coordinates of transverse grid lines from the top to the bottom
      read(8,*)
      read(8,*) NumCl
c NumCl number of cells in the system
      read(8,*)
      do i=1,NumCl
      read(8,*) xLU(i),yLU(i),xRB(i),yRB(i),MatNum(i)
      enddo
c xLU transverse coordinate of the upper left corner of the cell
c yLU vertical coordinate of the upper left corner of the cell
c xRB transverse coordinate of the lower right corner of the cell
c yRB vertical coordinate of the lower right corner of the cell
c MatNum code of cell boundary: 0 - impermeable, -4 - sink/source
c 
C  Calculate widths of the boundary strips U(j)
C
c     j=1 Horizontal border
      U(1)=DX
c     j=2 Vertical border
      U(2)=DY
C     j=3 Corner hanging in the domain
      U(3)=0.5*(DX+dY)
C     j=4: Corner attached to the vertical wall
      U(4)=0.5*DX

c
c                Nodal coordinates and boundary points
c
      NumNP=0
      NumBP=0
      do i=1,NumGL
      do j=1,NumVL
        xP=x(j)
        yP=y(i)
C
        iNod=0
c       iNod will be 0 if the point is not within a cell
        iBound=0
c       iBound will be k if the point is on the border of a cell # k
        Do k=1,NumCL
C         M1=0 if the point is not within a cell
          M1=In_cell(xP,yP,k,NumCl,xLU,yLU,xRB,yRB,edx,edy)
C If the point is not within a cell then may be it is on the border?
          If(M1.eq.0) then
C           M2=0 if the point is not on boundary of the cell
C           M2 is the index of the position of the boundary point
C           see indexes for U;  M2=1,2 or 3.
            M2=In_boundary(xP,yP,k,NumCl,xLU,yLU,xRB,yRB,edx,edy)
C If the point is on the boundary then?
            If(M2.gt.0) then
c              M3=0 if the point is not on the walls
               M3=In_Segment(xP,yP,0.,y(1),0.,y(NumGl),edx,edy)+
     &    In_Segment(xP,yP,x(NumVL),y(1),x(NumVL),y(NumGL),edx,edy)
c              M3=1 if the point is on the external boundary
               If (M3.eq.1) then
                 If(M2.eq.3) then
                   N2=4
                   iBound=k
                 Else
                   M1=1
                 Endif
               Else
                 iBound=k
                 N2=M2
               Endif
            Endif
          Endif
        iNod=iNod+M1
        Enddo
        If(iNod.eq.0) then
          NumNP=NumNP+1
          xN(NumNP)=xP
          yN(NumNP)=yP
c         NodNum is the nodal number of grid point
          NodNum(i,j)=NumNP
          if(iBound.gt.0) then
            NumBP=NumBP+1
            KXB(NumBP)=NumNP
            Width(NumBP)=U(N2)
            KodeG(NumBP)=MatNum(iBound)
          Endif
        Endif
      Enddo
      Enddo
C
C     Before adding boundary nodes find external coordinates of 'lids'
C
      ymax=-1000.
      ymin=1000.
      Do j=1,NumBP
        If(yN(KXB(j)).gt.ymax)  ymax=yN(KXB(j))
        If(yN(KXB(j)).lt.ymin)  ymin=yN(KXB(j))
      Enddo
c
c     Add boundary nodes at the top and at the bottom of the domain
c
      Do n=1,NumNP
        if(abs(yN(n)-y(1)).lt.edy) then
          NumBP=NumBP+1
          KXB(NumBP)=n
          Width(NumBP)=Dx
          KodeG(NumBP)=1
          If(abs(xN(n)).lt.edx. or.abs(xN(n)-x(NumVL)).lt.edx) 
     &    Width(NumBP)=Dx/2.
        Endif
CYAP Two-sided leaf
        If(IAX.eq.2.and.abs(yN(n)-y(NumGL)).lt.edy) then
          NumBP=NumBP+1
          KXB(NumBP)=n
          Width(NumBP)=Dx
          KodeG(NumBP)=1
          If(abs(xN(n)).lt.edx. or.abs(xN(n)-x(NumVL)).lt.edx) 
     &    Width(NumBP)=Dx/2.
        Endif
      Enddo
c
c                       Elements
c
      NumEl=0
      do i=1,NumGL-1
      do j=1,NumVL-1
C Center of gravity
        xC=(x(j)+x(j+1))/2.
        yC=(y(i)+y(i+1))/2.
        iEl=0
        Do k=1,NumCL
C         M1=0 if the point is out of cell
          M1=In_cell(xC,yC,k,NumCl,xLU,yLU,xRB,yRB,edx,edy)
          iEl=iEl+M1
        Enddo
        If(iEl.eq.0) then
          NumEl=NumEl+1
          KX(NumEl,1)=NodNum(i,j)
          KX(NumEl,2)=NodNum(i+1,j)
          KX(NumEl,3)=NodNum(i+1,j+1)
          KX(NumEl,4)=NodNum(i,j+1)
        Endif
c  Check on errors
        Do k=1,4
          if(KX(NumEl,k).eq.0) write(*,*) 'Error in El # ', NumEl
        Enddo
      Enddo
      Enddo
C
C Total width of boundary
C
      Total=0.
      Do m=1,NumBP
        if(KodeG(m).eq.-4) Total=Total+Width(m)
      Enddo
C
C Width of boundaries per unit length of the leaf surface
C
      CAI=Total/Float(NumVL-1)/Dx
    
c
c                         Output to datafile
c
      open(9,file='grid_bnd.dat')
      write(9,101)
      write(9,102)
      write(9,103) NumNP,NumEl,NumBP,NumVL,IAX
      write(9,104)
101   format('*** GRID GENERATOR INFORMATION')
102   format('NumNP   NumEl   NumBP    NumVL    IAX')
103   format(1x,4I8)
104   format('   n           x           y  ')
      do n=1,NumNP
        LocInd=1
        If(yN(n).ge.ymax) LocInd=0
        If(IAX.eq.2.and.yN(n).le.ymin) LocInd=0	 
        write(9,105) n,xN(n),yN(n),LocInd
      enddo
105   format(i5,2E12.3,i8)
      write(9,106)
      write(9,107)
106   format('*** ELEMENT INFORMATION ')
107   format('          e          i          j          k',
     &       '          l')
      do e=1,NumEl
      write(9,108) e,(KX(e,k),k=1,4)
108   format(5I5)
      enddo
      write(9,109)
      write(9,110)
109   format('Boundary geometry information')
110   format('  n CodeG  Width')
      Do m=1,NumBP
        write(9,111) KXB(m),KodeG(m),Width(m)
      Enddo
111   format(2I5,2E10.3)
      close(9)
c
c     Initial distribution of concentrations
c
      Read(8,*)
      Read(8,*) NumG
      Read(8,*)
      Read(8,*) (gInit(j),j=1,NumG)
      Open(9,file='Nodal_G.dat')
      Write(9,*)
      Write(9,*)
      Do i=1,NumNP
        Write(9,120) i,(gInit(j),j=1,NumG)
      Enddo
120   Format(I4,3E11.4)
      Close(9)
      Open(7,file='CAI.dat')
      Write(7,*) 'Total length   Cell area index'
      Write(7,130) Total,CAI
      Close(7)
      Write(*,131) CAI
130   Format(2E15.4)
131   Format('Cell area index = ',F8.3,'; you have it in CAI.DAT')
      Read(8,*)
      Read(8,*) NumCavRect
c NumCavRect number of rectangles to cover the assumed stomatal cavity
      read(8,*)
      do i=1,NumCavRect
      read(8,*) xLUC(i),yLUC(i),xRBC(i),yRBC(i)
      enddo
c xLU transverse coordinate of the upper left corner of the rectangle
c yLU vertical coordinate of the upper left corner of the rectangle
c xRB transverse coordinate of the lower right corner of the rectangle
c yRB vertical coordinate of the lower right corner of the rectangle
      If(IAX.EQ.2) then
        Read(8,*)
        Read(8,*) yCav_top, yCav_bottom
      Else
        Read(8,*)
        Read(8,*) yCav_top
      Endif
c yCav_top    vertical coordivate of the bottom of the guard cell row on top
c yCav_bottom vertical coordivate of the top of the guard cell row on bottom
c
C
C Added to output average concentrations in stomatal cavity
C
      Do e=1,NumEl
        iflag=0
        Do l=1,NumCavRect
        Do k=1,4
          j=KX(e,k)
          if(xN(j).gt.XLUC(l)-edx.and.xN(j).lt.XRBC(l)+edx
     &     .and.yN(j).gt.YRBC(l)-edy.and.yN(j).lt.YLUC(l)+edy)
     &    iflag=1
        Enddo
        Enddo
      iCav(e)=iflag
      Enddo
      open(11,file='aver.dat')
      write(11,'(30I2)') (iCav(e),e=1,NumEl)
c
c First and last elements within a leaf
c    
      e=1
70    if(e.gt.NumEl) goto 71
        iflag=0
        Do k=1,4
          j=KX(e,k)
          if(yN(j).lt.ycav_top-edy) iflag=iflag+1
        Enddo
      if(iflag.ne.4) then
         e=e+1
         goto 70
      Else
         iFirst=e
      Endif
71    Write(11,'(I6)') iFirst
      iLast=NumEl
      If(IAX.eq.2) then
        e=iFirst
80      if(e.gt.NumEl) goto 81
        iflag=0
        Do k=1,4
            j=KX(e,k)
            if(yN(j).gt.ycav_bottom+edy) iflag=iflag+1
        Enddo
        if(iflag.eq.4) then
           e=e+1
           goto 80
        Else
           iLast=e
        Endif
81      Continue
      Endif
      Write(11,'(I6)') iLast
      Close(11)
      Stop
      end
*
C     In_Cell=0 if the point is out of cell
      INTEGER FUNCTION In_Cell(xP,yP,k,NumCl,xLU,yLU,xRB,yRB,edx,edy)
        Dimension xLU(NumCL),yLU(NumCL),xRB(NumCL),yRB(NumCL)
        iAux=1
        If(xP.lt.xLU(k)+edx. or. xP.gt.xRB(k)-edx. or. 
     &     yP.gt.yLU(k)-edy. or. yP.lt.yRB(k)+edy) iAux=0
        In_Cell=iAux
        Return
      End
*
C     In_Boundary=0 if the point is not on boudary of the cell
      INTEGER FUNCTION In_boundary(xP,yP,k,NumCl,xLU,yLU,xRB,yRB,
     &    edx,edy)
        Dimension xLU(NumCL),yLU(NumCL),xRB(NumCL),yRB(NumCL)
        iAux=0
        iX=0
        iY=0
        If(In_segment(xP,yP,xLU(k),yLU(k),xRB(k),yLU(k),edx,edy).eq.1)
     &  iX=iX+1
        If(In_segment(xP,yP,xLU(k),yRB(k),xRB(k),yRB(k),edx,edy).eq.1)
     &  iX=iX+1
        If(In_segment(xP,yP,xLU(k),yRB(k),xLU(k),yLU(k),edx,edy).eq.1)
     &  iY=iY+1
        If(In_segment(xP,yP,xRB(k),yRB(k),xRB(k),yLU(k),edx,edy).eq.1)
     &  iY=iY+1
        iAux=iX+2*iY
        If(IAux.gt.0) then
          Continue
        Endif
        In_Boundary=iAux
        Return
      End
*
c     In_Segment=0 if the point is not on the segment
      INTEGER FUNCTION In_segment(xP,yP,x1,y1,x2,y2,edx,edy)
        iAux=0
        If(abs(x1-x2).lt.edx) then
           If(abs(xP-x1).lt.edx) then
             If((yP-y1+edy)*(y2+edy-yP).GE.0.) then
               iAux=1
             Endif
           Endif
        Endif
        If(abs(y1-y2).lt.edy) then
          If(abs(yP-y1).lt.edy) then
             If((xP-x1+edx)*(x2+edx-xP).GE.0.) then
               iAux=1
             Endif
          endif
        Endif
        If(iAux.eq.1) then
          Continue
        Endif
        In_Segment=iAux
      Return
      End

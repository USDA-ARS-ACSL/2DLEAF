*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      subroutine Get_Grid_and_Boundary()
      include 'public.ins'
*
      integer e
      logical lConst
      im=20
      il=0
*
      open(40,file='Grid_Bnd.dat',status='old',ERR=10)
      write(52,*) "Grid and boundary data"
c      CALL WCOU@('Grid and boundary data',RED+INVERSE_GREEN,WINDOW(2))
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10) NumNP,NumEl,NumBP,IJ,IAX
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
* Reading nodal information
*
      do i=1,NumNP
       il=il+1
       read(40,*,ERR=10) n,x(n),y(n),LocInd(n)
      enddo
*  Reading element information
*
      Num=0
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      do e=1,NumEL
       il=il+1
       read(40,*,ERR=10) Num,(KX(Num,i),i=1,4)
      enddo
*  Reading boundary information
*
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      il=il+1
      read(40,*,ERR=10)
      im=im+1
      do i=1,NumBP
          il=il+1
          read(40,*,ERR=10) n,CodeG(n),Width(i)
          KXB(i)=n
      enddo
*       
      do n=1,NumNP
          ListN(n)=1
          ListNE(n)=0
      enddo
*  Area of elements and bandwidth MBand
*
      lConst=.true.
      MBand=1
      xMul=1.
      do e=1,NumEl
        Area(e)=0.
        NUS=4
        IF(KX(e,3).eq.KX(e,4)) NUS=3
        do k=1,NUS-2
          MB=1
          i=KX(e,1)
          j=KX(e,k+1)
          l=KX(e,k+2)
          CJ=x(i)-x(l)
          CK=x(j)-x(i)
          BJ=y(l)-y(i)
          BK=y(i)-y(j)
          AE=(CK*BJ-CJ*BK)/2.
          Area(e)=Area(e)+xMul*AE
          ListNE(i)=ListNE(i)+1
          ListNE(j)=ListNE(j)+1
          ListNE(l)=ListNE(l)+1
          if(abs(i-j).gt.MB) MB=abs(i-j)
          if(abs(i-l).gt.MB) MB=abs(i-l)
          if(abs(j-l).gt.MB) MB=abs(j-l)
          if(MB.gt.MBand) MBand=MB
          if(e.eq.1.and.k.eq.1) then
            MB1=MB
          else
            if(MB1.ne.MB) lConst=.false.
          end if
        enddo
      enddo
      MBand=MBand+1
      if(.not.lConst) IJ=NumNP
*
      return
10    call errmes(im,il)
      end


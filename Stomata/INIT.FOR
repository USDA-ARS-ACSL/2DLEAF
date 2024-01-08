      subroutine Initialize()  
      include 'public.ins'
      TYPE (QWINFO) winfo
      Integer status
    
c      CALL VGA@
c wcreate creates a text window
c      CALL WCREATE@(1,1,76,1,YELLOW+INVERSE_GREEN,WINDOW(1))
c      CALL WCREATE@(1,32,36,1,WHITE+INVERSE_GREEN,WINDOW(2))
c      CALL WCREATE@(54,32,8,1,Green+INVERSE_GREEN,WINDOW(3))
c      CALL WCREATE@(62,32,14,1,Green+INVERSE_GREEN,WINDOW(4))
c      CALL POPW@(WINDOW(1))
c      CALL POPW@(WINDOW(2))
c      CALL POPW@(WINDOW(3))
c      CALL POPW@(WINDOW(4))
c      CALL WCOU@(CENTRE@('STOMATA',74),RED+INVERSE_GREEN,WINDOW(1))
*
       open(unit=51,file="user",TITLE="window 51")
       winfo.X    = 1
       winfo.Y    = 1
       winfo.H    = 76
       winfo.W    = 1
       winfo.TYPE = QWIN$SET
       status     = SETWSIZEQQ(4, winfo)

       open(unit=52,file="user" ,TITLE="window 52")
        winfo.X    = 1
       winfo.Y    = 32
       winfo.H    = 36
       winfo.W    = 1
       winfo.TYPE = QWIN$SET
       status     = SETWSIZEQQ(4, winfo)
       open(unit=53,file="user" ,TITLE="window 53")
       winfo.X    = 54
       winfo.Y    = 32
       winfo.H    = 8
       winfo.W    = 1
       winfo.TYPE = QWIN$SET
       status     = SETWSIZEQQ(4, winfo)
       Open(unit=54,file="user", TITLE="window 54")
       winfo.X    = 62
       winfo.Y    = 32
       winfo.H    = 14
       winfo.W    = 1
       winfo.TYPE = QWIN$SET
       status     = SETWSIZEQQ(4, winfo)
       open(4,file='STOMATA.LOG')
       do i=1,NumNPD
       do j=1,NumGD
         g(i,j)=0.
       enddo
       CodeG(i)=0
       enddo
*
*
       NumBP=0
       do i=1,NumBPD
         KXB(i)=0
         Width(i)=0.
         FlTrans(i)=0.
         do j=1,NumGD
         do k=1,3
            VarBG(i,j,k)=0.
         enddo
         enddo
       enddo
*
       lInput=1
       do i=0,NumObjD
         tNext(i)=1.E+31
       enddo
       do i=1,4
       dtMx(i)=1.E+31
       tTDB(i)=1.E+31
       enddo
       tatm=1.E+31
      return
      end

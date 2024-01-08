c==== Error messages
      Subroutine ErrMes(N,M)
      include 'Public.ins'
      integer*2 i
      if(N.ne.0) then
      write(4,20) N,M
      write(LINE,10) N
      write(52,*) Line
c      CALL WCOU@(LINE,
c     !            RED+INVERSE_GREEN,WINDOW(2))
c      call get_key@(i)
c      call text_mode@
      stop
      endif
10    format('Error # ',i4)
20    format('Error # ',i4,',   Line # ',i4)
      return
      end
       

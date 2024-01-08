*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     Stomata  - numerical model of plant stomata                      *
*                                                                      *
*                      Release 01                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      include 'public.ins'
      call Initialize()
      call Get_Grid_and_Boundary()
1     continue
      call GasMover()
      call SetSurface()
c      call ScreenOutput()
      call Break()
      goto 1
      end         


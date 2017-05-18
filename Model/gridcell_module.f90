MODULE GRIDCELL_MOD

!  The Gridcell Module contains the subroutine gridcell which is used for 
!   determining if a point lies within an element.  Subroutine gridcell 
!   serves two purposes.  When passed an element number through the optional 
!   argument checkele, it determines if a particle is in that particular
!   element.  When the optional argument is omitted, the subroutine determines
!   in which element the particle is currently located.
!
! Originally created by:  Elizabeth North
! Modified by:            Zachary Schlag
! Created on:             2004
! Last Modified on:       19 Aug 2008

IMPLICIT NONE
PUBLIC

CONTAINS

  ! Subroutine gridcell determines if a particle is inside or outside
  ! an element (e.g., a grid cell). If the particle is on a boundary 
  ! segment or a grid point (node) of the element, then the particle
  ! is considered to be in the grid cell.  If the optional argument
  ! checkele is present, then it only checks one particular element, 
  ! otherwise it checks all the elements on the given grid.
  SUBROUTINE gridcell(elements,ele_y,ele_x,Xpos,Ypos,P_ele,triangle,checkele)
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: elements
    INTEGER, INTENT(OUT) :: P_ele,triangle
    INTEGER, INTENT(IN), OPTIONAL :: checkele
    DOUBLE PRECISION, INTENT(IN)  :: ele_x(4,elements),ele_y(4,elements), &
                                     Xpos,Ypos

    INTEGER :: i,p,counter(4),total,elestart,eleend
    DOUBLE PRECISION :: slope,xintersect,by1,by2,bx1,bx2,blowy,bhighy

    if( present(checkele) )then   !if the variable checkele is present:
      elestart = checkele         !  initialize elestart and eleend to checkele
      eleend = checkele           !  so checkele is the only element checked
    else                          !else:
      elestart = 1                !  initialize elestart to 1 and
      eleend = elements           !  eleend to elements
    endif                         !  so that all the elements are checked

    do i=elestart,eleend
      counter = 0

      ! 1 check if particle is to west or east of element 
      if ( (Ypos.LT.ele_y(1,i) .AND. Ypos.LT.ele_y(2,i) .AND.   &
            Ypos.LT.ele_y(3,i) .AND. Ypos.LT.ele_y(4,i) ) .OR.  &
           (Ypos.GT.ele_y(1,i) .AND. Ypos.GT.ele_y(2,i) .AND.   &
            Ypos.GT.ele_y(3,i) .AND. Ypos.GT.ele_y(4,i) ) ) cycle

      ! 2 check if particle is to north or south of element
      if ( (Xpos.LT.ele_x(1,i) .AND. Xpos.LT.ele_x(2,i) .AND.   &
            Xpos.LT.ele_x(3,i) .AND. Xpos.LT.ele_x(4,i) ) .OR.  &
           (Xpos.GT.ele_x(1,i) .AND. Xpos.GT.ele_x(2,i) .AND.   &
            Xpos.GT.ele_x(3,i) .AND. Xpos.GT.ele_x(4,i) ) ) cycle

      ! 3  check if particle on boundary point exactly, if so then in element
      if ((Xpos.EQ.ele_x(1,i) .AND. Ypos.EQ.ele_y(1,i)) .OR.    &
          (Xpos.EQ.ele_x(2,i) .AND. Ypos.EQ.ele_y(2,i)) .OR.    &
          (Xpos.EQ.ele_x(3,i) .AND. Ypos.EQ.ele_y(3,i)) .OR.    &
          (Xpos.EQ.ele_x(4,i) .AND. Ypos.EQ.ele_y(4,i)) ) then
        P_ele = i   !particle on boundary point
        triangle = 1
        exit
      endif

      ! 4 check if there is a horizontal boundary segment and if 
      !     particle is on it
      if (ele_y(1,i).EQ.ele_y(2,i) .OR. ele_y(1,i).EQ.ele_y(3,i) .OR.     &
          ele_y(1,i).EQ.ele_y(4,i) .OR. ele_y(2,i).EQ.ele_y(3,i) .OR.     &
          ele_y(2,i).EQ.ele_y(4,i) .OR. ele_y(3,i).EQ.ele_y(4,i))then

        if (ele_y(1,i).EQ.ele_y(2,i) .AND. Ypos.EQ.ele_y(1,i) ) then
          if( ( ele_x(1,i).GT.ele_x(2,i) .AND. Xpos.GT.ele_x(2,i) .AND.   &
              Xpos.LT.ele_x(1,i) ) .OR. ( ele_x(2,i).GT.ele_x(1,i) .AND.  &
              Xpos.GT.ele_x(1,i) .AND. Xpos.LT.ele_x(2,i) ) ) then
            P_ele = i   !particle on boundary segement
            triangle = 1
            exit
          else
            cycle  !particle not on bnd segement so not in element
          endif
        endif

        if (ele_y(1,i).EQ.ele_y(3,i) .AND. Ypos.EQ.ele_y(1,i) ) then
          if( ( ele_x(1,i).GT.ele_x(3,i) .AND. Xpos.GT.ele_x(3,i) .AND.   &
              Xpos.LT.ele_x(1,i) ) .OR. ( ele_x(3,i).GT.ele_x(1,i) .AND.  &
              Xpos.GT.ele_x(1,i) .AND. Xpos.LT.ele_x(3,i) ) ) then
            P_ele = i   !particle on boundary segement
            triangle = 1
            exit
          else
            cycle  !particle not on bnd segement so not in element
          endif
        endif

        if (ele_y(1,i).EQ.ele_y(4,i) .AND. Ypos.EQ.ele_y(1,i) ) then
          if( ( ele_x(1,i).GT.ele_x(4,i) .AND. Xpos.GT.ele_x(4,i) .AND.   &
              Xpos.LT.ele_x(1,i) ) .OR. ( ele_x(4,i).GT.ele_x(1,i) .AND.  &
              Xpos.GT.ele_x(1,i) .AND. Xpos.LT.ele_x(4,i) ) ) then
            P_ele = i   !particle on boundary segement
            triangle = 1
            exit
          else
            cycle  !particle not on bnd segement so not in element
          endif
        endif

        if (ele_y(2,i).EQ.ele_y(3,i) .AND. Ypos.EQ.ele_y(2,i) ) then
          if( ( ele_x(2,i).GT.ele_x(3,i) .AND. Xpos.GT.ele_x(3,i) .AND.   &
              Xpos.LT.ele_x(2,i) ) .OR. ( ele_x(3,i).GT.ele_x(2,i) .AND.  &
              Xpos.GT.ele_x(2,i) .AND. Xpos.LT.ele_x(3,i) ) ) then
            P_ele = i   !particle on boundary segement
            triangle = 1
            exit
          else
            cycle  !particle not on bnd segement so not in element
          endif
        endif

        if (ele_y(2,i).EQ.ele_y(4,i) .AND. Ypos.EQ.ele_y(2,i) ) then
          if( ( ele_x(2,i).GT.ele_x(4,i) .AND. Xpos.GT.ele_x(4,i) .AND.   &
              Xpos.LT.ele_x(2,i) ) .OR. ( ele_x(4,i).GT.ele_x(2,i) .AND.  &
              Xpos.GT.ele_x(2,i) .AND. Xpos.LT.ele_x(4,i) ) ) then
            P_ele = i   !particle on boundary segement
            triangle = 1
            exit
          else
            cycle  !particle not on bnd segement so not in element
          endif
        endif

        if (ele_y(3,i).EQ.ele_y(4,i) .AND. Ypos.EQ.ele_y(3,i) ) then
          if( ( ele_x(3,i).GT.ele_x(4,i) .AND. Xpos.GT.ele_x(4,i) .AND.   &
              Xpos.LT.ele_x(3,i) ) .OR. ( ele_x(4,i).GT.ele_x(3,i) .AND.  &
              Xpos.GT.ele_x(3,i) .AND. Xpos.LT.ele_x(4,i) ) ) then
            P_ele = i   !particle on boundary segement
            triangle = 1
            exit
          else
            cycle  !particle not on bnd segement so not in element
          endif
        endif

      endif

      ! 5 check if the particle Y-location is equal to boundary y-coordinate
      if( Ypos.EQ.ele_y(1,i) .OR. Ypos.EQ.ele_y(2,i) .OR.  &
          Ypos.EQ.ele_y(3,i) .OR. Ypos.EQ.ele_y(4,i) ) then
        !Find which y-bnd coordinate is northernmost
        bhighy = ele_y(1,i)
        if (ele_y(2,i).GT.bhighy) then
          bhighy = ele_y(2,i)
        endif
        if (ele_y(3,i).GT.bhighy) then
          bhighy = ele_y(3,i)
        endif
        if (ele_y(4,i).GT.bhighy) then
          bhighy = ele_y(4,i)
        endif

        !If Ypos equals northermost corner coordinate, then not in element
        if (Ypos.EQ.bhighy) cycle  !particle not in element

        !Find which y-bnd coordinate is southernmost
        blowy= ele_y(1,i)
        if (ele_y(2,i).LT.blowy) then
          blowy = ele_y(2,i)
        endif
        if (ele_y(3,i).LT.blowy) then
          blowy = ele_y(3,i)
        endif
        if (ele_y(4,i).LT.blowy) then
          blowy = ele_y(4,i)
        endif

        !if Ypos equals northermost corner coordinate, then not in element
        if (Ypos.EQ.blowy) cycle  !particle not in element
      endif

    ! 6 All other cases 
      do p = 1, 4
        if (p.EQ.1) then
          bx1 =  ele_x(1,i)
          by1 =  ele_y(1,i)
          bx2 =  ele_x(2,i)
          by2 =  ele_y(2,i)
        endif
        if (p.EQ.2) then
          bx1 =  ele_x(2,i)
          by1 =  ele_y(2,i)
          bx2 =  ele_x(3,i)
          by2 =  ele_y(3,i)
        endif
        if (p.EQ.3) then        
          bx1 =  ele_x(3,i)
          by1 =  ele_y(3,i)
          bx2 =  ele_x(4,i)
          by2 =  ele_y(4,i)
        endif
        if (p.EQ.4) then        
          bx1 =  ele_x(4,i)
          by1 =  ele_y(4,i)
          bx2 =  ele_x(1,i)
          by2 =  ele_y(1,i)
        endif


        !  if particle longitude (X position) is less than (to the west of)
        !  one boundary longitudes (X position)
        if (Xpos.LE.bx1 .OR. Xpos.LE.bx2) then
          !if particle Y-pos is inbetween bnd Y-positions
          if( ( by1.GT.by2 .AND. Ypos.GE.by2 .AND. Ypos.LE.by1 ) .OR.     &
              ( by2.GT.by1 .AND. Ypos.GE.by1 .AND. Ypos.LE.by2  ) ) then  
        
            if (bx1.EQ.bx2) then  !boundary line is vertical

              if (Xpos.EQ.bx1) then !particle on boundary segement
                P_ele = i
                triangle=1  !1 = in grid
                exit          
              else
                counter(p) = 1
                if (Ypos.EQ.by2) counter(p) = 0 
              endif

            else !boundary line is not vertical

              slope = ( (by1-by2)/(bx1-bx2) )
              xintersect = ( Ypos-by1+(slope*bx1) ) / slope
              if (xintersect.GT.Xpos) then
                counter(p) = 1 
                if (Ypos.EQ.by2) counter(p) = 0 
              endif
              if (xintersect.EQ.Xpos) then !particle on boundary segement
                P_ele = i
                triangle=1  !1 = in grid
                exit          
              endif  
            endif
          endif
        endif
      enddo

      total = counter(1) + counter(2) + counter(3) + counter(4)   
      if (mod(total,2).NE.0) then 
        P_ele = i
        triangle=1  !1 = in grid
        exit 
      endif
    enddo

  END SUBROUTINE gridcell


END MODULE GRIDCELL_MOD
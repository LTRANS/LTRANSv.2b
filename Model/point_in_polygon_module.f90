MODULE PIP_MOD  !Point-In-Polygon Module

!  The Point-in-Polygon Module contains one function, inpoly, which determines
!   if a point is inside or outside an irregularly shaped polygon using the 
!   'crossings method', a 'point-in-polygon' technique.
!
!  Created by:            Zachary Schlag
!  Created on:            07 May 2007
!  Last Modified on:      11 Aug 2008

IMPLICIT NONE
PUBLIC

CONTAINS

  !Function inpoly returns .TRUE. if the point (x,y) lies within 
  !  the polygon defined by the points in the variable e(n,2)
  !  It uses the point-in-polygon approach, shooting a ray out to the
  !  right of the point along the point's x coordinate, if an odd number 
  !  of boundaries are crossed then the point is in the polygon.
  !The variable onin, is an optional variable used to designate whether
  !  a point on a boundary line is considered as in the polygon or out
  !  if the variable is not included in the function call, the function
  !  defaults to on a line being in the polygon
  LOGICAL FUNCTION inpoly(x,y,n,e,onin)
  !Input Variables
  INTEGER, INTENT(in) :: n                  ! number of edge points in polygon
  LOGICAL, INTENT(in), OPTIONAL :: onin     ! define polygon edge as in or out
  DOUBLE PRECISION, INTENT(in) :: x,y, &    ! particle x & y location
                                  e(n,2)    ! polygon edge coordinates  

  !Additional Variables
  INTEGER :: i,j,hilo(n),crossed
  LOGICAL :: on,first,onout
  DOUBLE PRECISION :: m,b,ix

  if( present(onin) )then       !if the variable onin is present
    onout=.not.onin             !set onout to the opposite of onin
  else
    onout=.false.               !if onin is not present, default onout to false
  endif

  inpoly=.true.
  crossed=0

  !cycle through all edge points, determine if any are on the ray
  hilo=0        !initialize hilo to 0
  on=.false.    !initialize on to false, to indicate no edge points on the ray
  do i=1,n
    if( e(i,2) > y ) hilo(i)=1                  !if above ray, set hilo to  1
    if( e(i,2) < y ) hilo(i)=-1                 !if below ray, set hilo to -1
    if( e(i,2) == y .and. e(i,1)>x) on=.true.   !if on ray, set on true
    if( e(i,1)==x .and. e(i,2)==y)then          !if edge point on point:
      if(onout)inpoly=.false.                   ! if on means out inpoly=false
      return                                    ! return from function
    endif
  enddo

  !If any of the edge points are on the ray:
  if(on)then

    !initialize first to true in case the 'first' boundary point is on the ray
    first=.true.

    !initialize i to 1, used to iterate through boundary points
    i=1

    do
      if(i>n)exit                           !if all boundaries checked, exit
      if(hilo(i)==0 .and. e(i,1)>x)then     !if on ray:

        if(first)then                       !if first point:
          i=i+1                             !  increment i
          cycle                             !  and move on to next point
        endif

        !if the previous boundary point shares the point's x coordinate, but
        !  was not on the ray, then the boundary crossed the point:
        if(hilo(i-1)==0)then
          if(onout)inpoly=.false.           !if on means out set inpoly false
          return                            !return from the function
        endif

        !initialize j to 1, used to iterate through boundary points until
        !  the next boundary point is found which does not lie on the ray.
        !  if it iterates past the last boundary point, it cycles back to
        !  the beginning
        j=1
        do
          if((i+j)==(n+1))j=2-i     !if the end is reached, cycle to start
          if(hilo(i+j)/=0)exit      !if the boundary point is not on ray, exit
          !if the next boundary point not on the ray and shares the point's x
          !  coordinate, then the boundary crossed the point:
          if(e(i+j,1)<x)then
            if(onout)inpoly=.false.         !if on means out set inpoly false
            return                          !return from the function
          endif
          j=j+1     !increment j to check next boundary point
        enddo

        !if the boundary points on the ray are inbetween the preceeding and
        !  following boundary points (i.e., one of the points is above the ray
        !  and the other is below the ray), then increment crossed to indicate
        !  that the ray has been crossed
        if( (hilo(i-1)+hilo(i+j))==0 )crossed=crossed+1

        !if j<0 then it has iterated past the last boundary point and cycled
        !  back to the beginning, therefore all the boundaries have been
        !  checked and the loop can be exited
        if( j<0 ) exit

        !if it hasn't passed the last boundary point, move to 
        !  the next point not on the ray
        i=i+j
      endif
      first=.false. !if the point wasn't on the ray, set first to false
      i=i+1         !increment i to iterate to next boundary point
    enddo

  endif

  !Check if the ray crosses each boundary segment
  do i=1,n-1

    !if both boundary points are left of the point, cycle
    if(e(i,1)<=x .and. e(i+1,1)<=x)cycle

    !if both boundary points are below the point, cycle
    if(e(i,2)<=y .and. e(i+1,2)<=y)cycle

    !if both boundary points are above the point, cycle
    if(e(i,2)>=y .and. e(i+1,2)>=y)cycle

    !if both boundary points are to the right of the point, increment crossed
    !  and cycle (if they are not both above or both below, but are both
    !  right of, then it has to cross)
    if(e(i,1)> x .and. e(i+1,1)> x)then
      crossed=crossed+1
      cycle
    endif

    !if it reaches the following code then one boundary point is above, one is 
    !  below, and one is to the right of the point, while the other is left of 
    !  the point therefore it may or may not cross the ray, so the slope and 
    !  intercept must be calculated, then the intersection point of the ray 
    !  with the boundary line must be calculated using the y coordinate of the 
    !  point:

    m=  (e(i+1,2)-e(i,2))/(e(i+1,1)-e(i,1)) !slope
    b=  e(i,2)-m*e(i,1)                     !intercept
    ix= (y-b)/m                             !ix = intersection point

    if(ix==x)then                           !if x is on the intersection point:
      if(onout)inpoly=.false.               ! if on means out set inpoly false
      return                                ! return from the function
    endif

    !if x intersect is to the right of the point, then ray crosses boundary 
    !  line, therefore increment the variable crossed
    if(ix>x)crossed=crossed+1
  enddo

  !If the ray crosses an even # of times, then the point is not in the polygon
  if( mod(crossed,2)==0 )inpoly=.false.

  RETURN
  END FUNCTION inpoly

END MODULE PIP_MOD
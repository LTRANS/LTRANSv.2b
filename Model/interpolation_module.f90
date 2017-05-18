MODULE INT_MOD

!  The Interpolation Module contains two procedures that interpolate data:
!  Subroutine linint uses linear interpolation.
!  Subroutine polintd uses polynomial interpolation.
!
!  Created by:            Elizabeth North
!  Modified by:           Zachary Schlag
!  Created on:            2003
!  Last Modified on:      11 Aug 2008

IMPLICIT NONE
PUBLIC

CONTAINS

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                     SUBROUTINE linint                         ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE linint(xa,ya,n,x,y,m)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: xa(n),ya(n),x
    DOUBLE PRECISION, INTENT(OUT) :: y,m

    INTEGER :: jlo,jhi,k
    DOUBLE PRECISION :: b

    jlo=1  
    jhi=n

    !determine which two xa location that x lies between
    do
      k=(jhi+jlo)/2

      if(xa(k).gt.x)then
        jhi=k
      else
        jlo=k
      endif

      if (jhi-jlo == 1) exit
    enddo

    !calculate the slope and intersect
    m = ( ya(jlo) - ya(jhi) )/ ( xa(jlo) - xa(jhi) )
    b = ya(jlo) - m*xa(jlo)

    !linearly interpolate y
    y = m*x + b

    return

  END SUBROUTINE linint


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~                                                               ~~
! ~~                     FUNCTION polintd                          ~~
! ~~                                                               ~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  DOUBLE PRECISION FUNCTION polintd(xa,ya,n,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: xa(n),ya(n),x

    INTEGER :: i,ns
    DOUBLE PRECISION :: dif,dift,a,b,c

    ns=1
    dif=abs(x-xa(1))

    !Here we find the index, ns, of the closest table entry
    do i=1,n
      dift=abs(x-xa(i))
      if (dift.lt.dif) then
        ns=i
        dif=dift
      endif
    enddo 

    !calculate the value of c
    c = (xa(2)-x) * ( (ya(3)-ya(2)) / (xa(2)-xa(3)) )
    c = c - (xa(2)-x) * ( (ya(2)-ya(1)) / (xa(1)-xa(2)) )
    c = c / (xa(1) - xa(3))

    !calculate the values of a and b
    if (ns .EQ. 3) then
      a = (ya(3)-ya(2))/(xa(2)-xa(3))
      b = xa(3)-x
    else
      a = (ya(2)-ya(1))/(xa(1)-xa(2))
      b = xa(1)-x
    endif

    !calculate the value of polintd via polynomial interpolation
    polintd = ya(ns) + (xa(ns)-x)*a + b*c

  END FUNCTION polintd

END MODULE INT_MOD
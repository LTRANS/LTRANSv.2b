MODULE BOUNDARY_MOD

! This module contains variables and subroutines associated with the 
! creation of the land/sea boundaries.  The main purpose of this module 
! is to create the land/sea boundaries from a given masked rho grid.  
! The main subroutine in the module, createBounds, determines the number
! of boundary points, allocates an array of that size, and fills it with 
! the boundary points in order.

! Module and most subroutines created by:  Zachary Schlag
! Created on:           03 Apr 2008
! Last Modified on:        Feb 2011
!
! Subroutine intersect_reflect created by: Elizabeth North
! Created on:                  2005
! Last Modified on:     08 Aug 2008

IMPLICIT NONE
PRIVATE
SAVE

!*****************************************************************
!*                          VARIABLES                            *
!*****************************************************************

  TYPE :: boundary  !boundary point variable
    LOGICAL :: onU      !TRUE if point is on U grid, FALSE if on V Grid
    INTEGER :: ii       !i grid position of point
    INTEGER :: jj       !j grid position of point
    INTEGER :: poly     !polygon of point (1 = main body of water)
                        !                 (>1 = island)
  END TYPE boundary

  !initial holder of all boundary points, 
  !  deallocated after reformatted into hid,hx,hy,bx,by
  TYPE (boundary), ALLOCATABLE,  DIMENSION (:) :: bnds

  !final boundary variables, after reformatting from bnds
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: hx,hy,bx,by
  INTEGER, ALLOCATABLE, DIMENSION (:) :: hid
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: bnd_x,bnd_y

  !TRUE if the boundary is land, FALSE if it is open ocean
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: land

  INTEGER :: bnum,        & !total number of boundary points
             maxbound,    & !number of boundary points around water
             maxisland,   & !number of boundary points around islands
             numislands,  & !number of islands
             nbounds        !total number of boundary segments

  LOGICAL :: BND_SET = .false.  !Boundaries Set?

  !The following procedures have been made public:
  PUBLIC :: isBndSet,createBounds,mbounds,ibounds,intersect_reflect

CONTAINS

!*****************************************************************
!*                    FUNCTIONS & SUBROUTINES                    *
!*****************************************************************



    !*********************************************************
    !*                      isBndSet                         *
    !*********************************************************

LOGICAL FUNCTION isBndSet()
!This function simply returns TRUE if the boundaries have been created
!  and FALSE if they have not yet been created
  isBndSet = BND_SET
  RETURN
END FUNCTION isBndSet



    !*********************************************************
    !*                     createBounds                      *
    !*********************************************************

SUBROUTINE createBounds()
  USE PARAM_MOD, ONLY:ui,uj,vi,vj,BoundaryBLNs
  USE HYDRO_MOD, ONLY:getMask_Rho,getUVxy
!This subroutine creates boundaries based on the masking of the rho grid
!  The boundary points are U & V grid points directly between the rho points
!  This subroutine assumes there is only one body of water in the grid 
!    but there may exist multiple islands within the water
!  It creates the boundary points to go clockwise around the body of water
!    It then will create boundary points clockwise around any islands
!
!The subroutine does the following things:
!  1) Determines the number of boundary points in the given rho grid
!  2) Determines the form of each rho element
!    a) Elements consist of the rho nodes:  (i+1,j) --- (i+1,j+1)
!                                              |            |
!                                            (i,j)  ---  (i,j+1)
!
!    b) Element form is based on how the four nodes in the element are masked
!
!    c) If any elements exist where their form contains water and land crossing
!         diagonally, which are referred to as 'crosses' then these forms must
!         be solved as to which direction the boundarys are going through them
!
!  3) Starting at the first element it encounters with a water masked node, 
!       it goes around the body of water clockwise adding each boundary point
!       until it reaches the boundary point that it started on
!
!  4) If the number of boundary points used is not equal to the total number of
!       boundary points originally found in the grid, then find an unused
!       boundary element to begin on and go clockwise around the island it is 
!       part of, this is continued until all the boundary points are used


  LOGICAL :: U = .true.     !sent to subroutine ADD to indicate on U grid
  LOGICAL :: V = .false.    !sent to subroutine ADD to indicate on V grid
  LOGICAL :: wf = .false.   !short for waterfall, indicates path is going
                            !  through an element on the edge of rho grid
  INTEGER :: dir = 5        !direction boundary left the previous element
                            !  based on numeric keypad (8-up,6-right,etc)

  INTEGER :: ipos,jpos      !i and j position of current element when 
                            !  making bounds

  !The following variables are the same as the similarly named variables
  INTEGER :: ipos1,jpos1,ipos2,jpos2,dir1,dir2  ! above, only these are
  LOGICAL :: wf1, wf2                           ! used for solving crosses

  INTEGER :: i,j,STATUS,m,crossnum,oldcrossnum,count
  LOGICAL :: found, deadend1, deadend2

  LOGICAL :: polydone       !TRUE if current polygon is closed, else FALSE
  LOGICAL :: done = .false. !TRUE if all boundaries are used, else FALSE


  TYPE :: element       !element variable
    INTEGER :: form     !0-19 depends on masking of 4 nodes & situation
    LOGICAL :: used     !TRUE if element is already used, FALSE if not
                        !Note that in a 'cross' element, used will not
                        !  change to TRUE until it has been used twice
    LOGICAL :: unused   !TRUE if a 'cross' element has never been used,
                        !  FALSE if it has been used at least once
  END TYPE element

  TYPE (element), ALLOCATABLE, DIMENSION(:,:) :: ele        !elements

  !Used to reformat from bnds to bx,by,hx,hy,hid
  INTEGER :: k,c,numpoly,pstart,pend
  INTEGER, ALLOCATABLE, DIMENSION(:  ) :: polysizes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_u,y_u,x_v,y_v

  !Rho Mask Used to create boundaries
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: mask_rho

  !Keep track of boundary points on open ocean
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: oo

  ALLOCATE(x_u(ui,uj))
  ALLOCATE(y_u(ui,uj))
  ALLOCATE(x_v(vi,vj))
  ALLOCATE(y_v(vi,vj))
  ALLOCATE(mask_rho(vi,uj))

  CALL getMask_Rho(mask_rho)

  !Edit Rho Mask - remove nodes that don't have at least 2 neighbors
  !  This removes a situation that may occur where an area that is 
  !  within the created boundaries, is not covered by either the U 
  !  or V grid.
  do
    m=0
    do j=1,uj
      do i=1,vi
        if(mask_rho(i,j) == 1)then
          count = 0
          if(i> 1)then
            if(mask_rho(i-1,j)==1) count = count + 1
          endif
          if(i<vi)then
            if(mask_rho(i+1,j)==1) count = count + 1
          endif
          if(j> 1)then
            if(mask_rho(i,j-1)==1) count = count + 1
          endif
          if(j<uj)then
            if(mask_rho(i,j+1)==1) count = count + 1
          endif
          if(count < 2) then
            mask_rho(i,j) = 0
            m = m + 1
          endif
        endif
      enddo
    enddo
  if(m==0)exit
  enddo

    !*****************************************
    !*              SET BOUNDS               *
    !*****************************************

  !This section determines the number of boundary points in the given rho grid

!  write(*,*) 'Calculating Number of Boundary Points'
!  write(*,*) ' '

  bnum = 0  !initialize bnum to 0

  do i=2,vi-1
    do j=2,uj-1

      !Excluding the edges of the rho grid, increment bnum every time adjacent
      !  nodes are found that are of differing mask values
      if(i<vi-1)then
        if((mask_rho(i,j).EQ.0 .AND. mask_rho(i+1,j).EQ.1) .OR.                &
           (mask_rho(i,j).EQ.1 .AND. mask_rho(i+1,j).EQ.0)) bnum = bnum+1
      endif

      if(j<uj-1)then
        if((mask_rho(i,j).EQ.0 .AND. mask_rho(i,j+1).EQ.1) .OR.                &
           (mask_rho(i,j).EQ.1 .AND. mask_rho(i,j+1).EQ.0)) bnum = bnum+1
      endif

    enddo
  enddo

  !bnum must also be incremented where the second node from the edge of the 
  !  rho grid is masked as water because the boundary will need to go between
  !  this point and the edge point
  do i=2,vi-1
    if(mask_rho(i,2).EQ.1) bnum = bnum+1
    if(mask_rho(i,uj-1).EQ.1) bnum = bnum+1
  enddo

  do j=2,uj-1
    if(mask_rho(2,j).EQ.1) bnum = bnum+1
    if(mask_rho(vi-1,j).EQ.1) bnum = bnum+1
  enddo

  !Allocate the variable bnds to the total number of boundaries (bnum)
  ALLOCATE (bnds(bnum),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating bnds'

  !Initialize all the values in newly allocated bnds
  do i=1,bnum
    bnds(i)%onU  = .false.
    bnds(i)%ii   = 0
    bnds(i)%jj   = 0
    bnds(i)%poly = 0
  enddo


    !*****************************************
    !*              GET FORMS                *
    !*****************************************

  !This section determines the form of each rho element
  !
  !Remember that:
  !  a) Elements consist of the rho nodes:  (i+1,j) --- (i+1,j+1)
  !                                            |            |
  !                                          (i,j)  ---  (i,j+1)
  !
  !  b) Element form is based on how the four nodes in the element are masked
  !
  !  c) If any elements exist where their form contains water and land crossing
  !       diagonally, which are referred to as 'crosses' then these forms must
  !       be solved as to which direction the boundarys are going through them

  !Allocate the variable ele to the total number of elements
  ALLOCATE (ele(vi-1,uj-1),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating ele'

!  write(*,*) 'Getting Element Forms'
!  write(*,*) ' '

  crossnum = 0  !Initialize 'cross' element counter to 0

  do i=1,vi-1
    do j=1,uj-1
      ele(i,j)%used = .false.   !Initialize used to FALSE
      ele(i,j)%unused = .true.  !Initialize unused to TRUE (for 'crosses')

      !Determine the form of the current element
      ele(i,j)%form = 0
      if(mask_rho(i  ,j+1) == 0) ele(i,j)%form = ele(i,j)%form + 1
      if(mask_rho(i  ,j  ) == 0) ele(i,j)%form = ele(i,j)%form + 2
      if(mask_rho(i+1,j+1) == 0) ele(i,j)%form = ele(i,j)%form + 4
      if(mask_rho(i+1,j  ) == 0) ele(i,j)%form = ele(i,j)%form + 8

      SELECT CASE(ele(i,j)%form)

        CASE(6,9) !ud
        !Elements with forms 6 or 9 are 'crosses', increment the cross counter
          crossnum = crossnum + 1

        CASE(0,15) !0
        !Elements with forms 0 or 15 are all water or all land
        !Initialized used to TRUE because there are no boundaries through them
          ele(i,j)%used = .true.

      END SELECT

    enddo
  enddo

        !*********************************
        !*         SOLVE CROSSES         *
        !*********************************

  oldcrossnum = 0   !Initialize oldcrossnum to 0
                    !This variable is to know if an endless loop has occurred

  if(crossnum > 0) then
!    write(*,*) '  Crosses Exist... Now Solving'
!    write(*,*) ' '
  endif

  do 
    if(crossnum == 0) exit          !If no crosses exist unsolved, exit
    if(oldcrossnum == crossnum)then !If crosses still exist, but none were
!      write(*,*)'  Cross loop - using defaults'    !  solved in prior loop
!      write(*,*)' '                                !  use defaults of 16 & 19
      do i=1,vi-1
        do j=1,uj-1
          if(ele(i,j)%form == 6) ele(i,j)%form = 16
          if(ele(i,j)%form == 9) ele(i,j)%form = 19
        enddo
      enddo
    endif
    oldcrossnum = crossnum

    do i=1,vi-1
      if(crossnum == 0)exit         !If no crosses exist unsolved, exit

      do j=1,uj-1
        if(crossnum == 0)exit       !If no crosses exist unsolved, exit


        if(ele(i,j)%form == 6) then
          !Form 6 has water nodes in the top left and bottom right.
          !  Boundary edges travelling clockwise around water therefore
          !  enter the element from top or bottom and exit left or right


          if(i==1 .OR. j==uj-1 .OR. i==vi-1 .OR. j==1)then
            !if the cross is in an element on the edge of the rho grid,
            !  then the boundary will either enter the top and exit left,
            !  or enter the bottom and exit right
            ele(i,j)%form = 17 
            crossnum = crossnum - 1
          else


            !if its not on the edge of the grid it will have to be solved
            !  the hard way.  This is done by exiting the two possible
            !  exits (left & right) and following the boundaries until
            !  either, one returns to this cross or, both hit another
            !  unsolved cross.  If one returns to this cross then this 
            !  cross is solved.  If neither direction returns to this 
            !  cross, skip this cross.
            !  As other crosses are solved this one becomes more likely
            !  to find a path back to itself.

            ipos1 = i           !initialize two paths to the cross location
            jpos1 = j
            ipos2 = i
            jpos2 = j
            wf1 = .false.       !initialize wf & deadend variables to false
            wf2 = .false.
            deadend1 = .false.
            deadend2 = .false.

            ! Initialized ipos1,jpos1 to the element Left of the cross
            if(jpos1 == 2)then         !If element to left is on edge
              wf1 = .true.             !  switch wf to TRUE and move up
              jpos1 = jpos1 - 1
              ipos1 = ipos1 + 1
              dir1 = 8
              if(ipos1 == vi-1)then    !If element to left and up from
                jpos1 = jpos1 + 1      !  cross is the corner, move right
                dir1 = 6               !  from corner
              endif
            else
              jpos1 = jpos1 - 1        !Else just move left from cross
              dir1 = 4
            endif
          
            ! Initialized ipos2,jpos2 to the element Right of the cross
            if(jpos2 == uj-2)then      !If element to right is on edge
              wf2 = .true.             !  switch wf to TRUE and move down
              jpos2 = jpos + 1
              ipos2 = ipos - 1
              dir2 = 2
              if(ipos2 == 1)then       !If element to right and down from
                jpos2 = jpos2 - 1      !  cross is corner, move left from
                dir2 = 4               !  corner
              endif
            else
              jpos2 = jpos2 + 1        !Else just move right from cross
              dir2 = 6
            endif

            do
              if(ipos1 == i .AND. jpos1 == j)then   !If left path returned to
                if(dir1 == 2)then                   !  cross, set new form
                  ele(i,j)%form = 16                !  based on the direction
                elseif(dir1 == 8)then               !  the path returns from
                  ele(i,j)%form = 17
                else
                  write(*,*)'Problem Form 6 Left Solution'
                endif
                crossnum = crossnum - 1             !  and decrement crossnum
                exit

              elseif(ipos2 == i .AND. jpos2 == j)then  !If right path
                if(dir2 == 2)then                   !  returned to cross, set
                  ele(i,j)%form = 17                !  new form based on the
                elseif(dir2 == 8)then               !  direction the path
                  ele(i,j)%form = 16                !  returns from
                else
                  write(*,*)'Problem Form 6 Right Solution'
                endif
                crossnum = crossnum - 1             !  and decrement crossnum
                exit
              endif

              !if the left path has hit a dead end, switch deadend1 to TRUE
              if(ele(ipos1,jpos1)%form==6 .OR.             &
                 ele(ipos1,jpos1)%form==9) deadend1 = .true.

              !if the right path has hit a dead end, switch deadend2 to TRUE
              if(ele(ipos2,jpos2)%form==6 .OR.             &
                 ele(ipos2,jpos2)%form==9) deadend2 = .true.

              !if both paths have hit a dead end, this cross cannot be
              !  solved yet, so move on and come back to it later
              if(deadend1 .AND. deadend2) exit
            
              !if left path has not hit dead end, move to next point on path
              if(.NOT. deadend1) then
                CALL getNext(ipos1,jpos1,wf1,dir1,ele(ipos1,jpos1)%form)
              endif

              !if right path has not hit dead end, move to next point on path
              if(.NOT. deadend2) then
                CALL getNext(ipos2,jpos2,wf2,dir2,ele(ipos2,jpos2)%form)
              endif

            enddo

          endif


        elseif(ele(i,j)%form == 9)then
          !Form 9 has water nodes in the bottom left and top right.
          !  Boundary edges travelling clockwise around water therefore
          !  enter the element from left or right and exit top or bottom


          if(i==1 .OR. j==1 .OR. i==vi-1 .OR. j==uj-1)then
            !if the cross is in an element on the edge of the rho grid,
            !  then the boundary will either enter right and exit the top,
            !  or enter left and exit the bottom
            ele(i,j)%form = 18
            crossnum = crossnum - 1


            !if its not on the edge of the grid it will have to be solved
            !  the hard way.  This is done by exiting the two possible
            !  exits (up & down) and following the boundaries until
            !  either, one returns to this cross or, both hit another
            !  unsolved cross.  If one returns to this cross then this 
            !  cross is solved.  If neither direction returns to this 
            !  cross, skip this cross.
            !  As other crosses are solved this one becomes more likely
            !  to find a path back to itself.
          else
            ipos1 = i           !initialize two paths to the cross location
            jpos1 = j
            ipos2 = i
            jpos2 = j
            wf1 = .false.       !initialize wf & deadend variables to false
            wf2 = .false.
            deadend1 = .false.
            deadend2 = .false.


            ! Initialized ipos1,jpos1 to the element Above the cross
            if(ipos1 == vi-2)then      !If element above is on edge
              wf1 = .true.             !  switch wf to TRUE & move right
              ipos1 = ipos1 + 1
              jpos1 = jpos1 + 1
              dir1 = 6
              if(jpos1 == uj-1)then    !If element above and to right
                ipos1 = ipos1 - 1      !  of cross is the corner, move
                dir1 = 2               !  down from corner
              endif
            else
              ipos1 = ipos1 + 1        !Else just move above the cross
              dir1 = 8
            endif


            ! Initialized ipos1,jpos1 to the element Below the cross
            if(ipos2 == 2)then         !If element below is on edge
              wf2 = .true.             !  switch wf to TRUE & move left
              ipos2 = ipos2 - 1
              jpos2 = jpos2 - 1
              if(jpos2 == 1)then       !If element below and left of
                ipos2 = ipos2 + 1      !  cross is the corner, move
                dir2 = 8               !  up from corner
              else
                dir2 = 4
              endif
            else
              ipos2 = ipos2 - 1        !Else just move below the cross
              dir2 = 2
            endif


            do
              if(ipos1 == i .AND. jpos1 == j)then   !If up path returned to
                if(dir1 == 4)then                   !  cross, set new form
                  ele(i,j)%form = 19                !  based on the direction
                elseif(dir1 == 6)then               !  the path returns from
                  ele(i,j)%form = 18
                else
                  write(*,*)'Problem Form 9 Up Solution'
                endif
                crossnum = crossnum - 1             !  and decrement crossnum
                exit

              elseif(ipos2 == i .AND. jpos2 == j)then  !If down path returned
                if(dir2 == 4)then                   !  to cross, set new form
                  ele(i,j)%form = 18                !  based on the direction
                elseif(dir2 == 6)then               !  the path returns from
                  ele(i,j)%form = 19
                else
                  write(*,*)'Problem Form 9 Down Solution'
                endif
                crossnum = crossnum - 1             !  and decrement crossnum
                exit
              endif

              !if the up path has hit a dead end, switch deadend1 to TRUE
              if(ele(ipos1,jpos1)%form==6 .OR.             &
                 ele(ipos1,jpos1)%form==9) deadend1 = .true.

              !if the down path has hit a dead end, switch deadend2 to TRUE
              if(ele(ipos2,jpos2)%form==6 .OR.             &
                 ele(ipos2,jpos2)%form==9) deadend2 = .true.
 
              !if both paths have hit a dead end, this cross cannot be
              !  solved yet, so move on and come back to it later
              if(deadend1 .AND. deadend2) exit

              !if up path has not hit dead end, move to next point on path
              if(.NOT. deadend1) then
                CALL getNext(ipos1,jpos1,wf1,dir1,ele(ipos1,jpos1)%form)
              endif

              !if down path has not hit dead end, move to next point on path
              if(.NOT. deadend2) then
                CALL getNext(ipos2,jpos2,wf2,dir2,ele(ipos2,jpos2)%form)
              endif

            enddo

          endif
        endif
      enddo
    enddo

  enddo

  if(oldcrossnum > 0) then
!    write(*,*) '  Crosses Solved'
!    write(*,*) ' '
  endif


        !*********************************
        !*  FIND ELEMENT TO START FROM   *
        !*********************************

  found = .false.

  if(mask_rho(2,2) == 1)then    !check to see if the boundary will have to
    ipos = 1                    !  pass through the bottom left element
    jpos = 1
    found = .true.
  endif

  do i=2,vi-2
    if(found)exit               !if found, exit loop
    do j=2,uj-2
      if(ele(i,j)%form > 0 .AND. ele(i,j)%form < 15)then
        ipos = i                !if an element that isn't all water or all
        jpos = j                !  land has been found, store its location
        found = .true.          !  in ipos & jpos and set found to TRUE
        exit                    !  to exit the loops
      endif
    enddo
  enddo


        !*********************************
        !*       CREATE BOUNDARIES       *
        !*********************************

  polydone = .false.                   !Initialize polydone and wf to FALSE
  wf = .false.

!  write(*,*) 'Creating Boundaries'
!  write(*,*) ' '

  if(ipos == 1 .AND. jpos == 1)then    !If it starts in bottom left corner:
    CALL add(V,ipos+1,jpos,polydone,done)  !  Add V(2,1)
    wf = .true.                        !  Initialize wf to TRUE
    dir = 8                            !  Initialize dir to up (8)
    ipos = ipos + 1                    !  And move to element(2,1)
  endif

  do
    if(polydone)exit         !If water boundaries have been completed, exit


    if(ele(ipos,jpos)%form > 15 .AND. ele(ipos,jpos)%unused)then
      !if the element is a cross and has not been used, set unused to FALSE
      !  this indicates that it has been used once
      ele(ipos,jpos)%unused = .false.
    else
      !if the element is NOT a cross, or if it is & has been used once prior
      !  set used to TRUE to indicate that the element will never need to be
      !  used again
      ele(ipos,jpos)%used = .true.
    endif


    !NOTE: for more comprehensive explaination of following code see the
    !  subroutine getNext.  The following code is the same as the code
    !  found there, except this code contains calls to the subroutine
    !  add, which adds the boundaries to bnds.  getNext does not add
    !  boundaries so it does not contain calls to add.

    if(wf)then                              !if element on an edge
      if(dir == 2) then                     !right edge
        SELECT CASE(ele(ipos,jpos)%form)
          CASE(0,1,4,5)
            CALL add(V,ipos,jpos,polydone,done)
            ipos = ipos - 1
            if(ipos == 1)then
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              dir = 4
              jpos = jpos - 1
            endif
          CASE DEFAULT
            wf = .false.
        END SELECT
      elseif(dir == 4) then                 !bottom edge
        SELECT CASE(ele(ipos,jpos)%form)
          CASE(0,1,2,3)
            CALL add(U,ipos,jpos,polydone,done)
            jpos = jpos - 1
            if(jpos == 1)then
              if(.NOT. polydone)CALL add(V,ipos+1,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              dir = 8
              ipos = ipos + 1
            endif
          CASE DEFAULT
          wf = .false.
        END SELECT
      elseif(dir == 6) then                 !top edge
        SELECT CASE(ele(ipos,jpos)%form)
          CASE(0,4,8,12)
            CALL add(U,ipos,jpos+1,polydone,done)
            jpos = jpos + 1
            if(jpos == uj-1)then
              if(.NOT. polydone) CALL add(V,ipos,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              dir = 2
              ipos = ipos - 1
            endif
          CASE DEFAULT
            wf = .false.
        END SELECT
      elseif(dir == 8) then                 !left edge
        SELECT CASE(ele(ipos,jpos)%form)
          CASE(0,2,8,10)
            CALL add(V,ipos+1,jpos,polydone,done)
            ipos = ipos + 1
            if(ipos == vi-1)then
              if(.NOT. polydone)CALL add(U,ipos,jpos+1,polydone,done)
              ele(ipos,jpos)%used = .true.
              dir = 6
              jpos = jpos + 1
            endif
          CASE DEFAULT
            wf = .false.
        END SELECT
      else
        write(*,*) 'Error: wf direction not one of 2,4,6,8'
      endif


    else                                    !if element NOT on edge

      SELECT CASE(ele(ipos,jpos)%form)


        CASE(1,5,13)                        !These forms exit element bottom
          CALL add(V,ipos,jpos,polydone,done)
          if(ipos == 2)then
            wf = .true.
            ipos = ipos - 1
            if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done)
            ele(ipos,jpos)%used = .true.
            jpos = jpos - 1
            if(jpos == 1)then
              if(.NOT. polydone)CALL add(V,ipos+1,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              ipos = ipos + 1
              dir = 8
            else
              dir = 4
            endif
          else
            ipos = ipos - 1
            dir = 2
          endif


        CASE(2,3,7)                         !These forms exit element left
          CALL add(U,ipos,jpos,polydone,done)
          if(jpos == 2)then
            wf = .true.
            jpos = jpos - 1
            if(.NOT. polydone) CALL add(V,ipos+1,jpos,polydone,done)
            ele(ipos,jpos)%used = .true.
            ipos = ipos + 1
            if(ipos == vi-1)then
              if(.NOT. polydone)CALL add(U,ipos,jpos+1,polydone,done)
              ele(ipos,jpos)%used = .true.
              jpos = jpos + 1
              dir = 6
            else
              dir = 8
            endif
          else
            jpos = jpos - 1
            dir = 4
          endif


        CASE(4,12,14)                       !These forms exit element right
          CALL add(U,ipos,jpos+1,polydone,done)
          if(jpos == uj-2)then
            wf = .true.
            jpos = jpos + 1
            if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done)
            ele(ipos,jpos)%used = .true.
            ipos = ipos - 1
            if(ipos == 1)then
              if(.NOT. polydone) CALL add(V,ipos,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              jpos = jpos - 1
              dir = 4
            else
              dir = 2
            endif
          else
            jpos = jpos + 1
            dir = 6
          endif


        CASE(8,10,11)                       !These forms exit element top
          CALL add(V,ipos+1,jpos,polydone,done)
          if(ipos == vi-2)then
            wf = .true.
            ipos = ipos + 1
            if(.NOT. polydone) CALL add(U,ipos,jpos+1,polydone,done)
            ele(ipos,jpos)%used = .true.
            jpos = jpos + 1
            if(jpos == uj-1)then
              if(.NOT. polydone) CALL add(V,ipos,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              ipos = ipos - 1
              dir = 2
            else
              dir = 6
            endif
          else
            ipos = ipos + 1
            dir = 8
          endif


        CASE(16,17)                         !These forms exit right or left
                                            !16: u->r & d->l  17: d->r & u->l

          if((ele(ipos,jpos)%form == 16 .AND. dir == 2).OR.   & !if exit right
             (ele(ipos,jpos)%form == 17 .AND. dir == 8)) then

            CALL add(U,ipos,jpos+1,polydone,done)
            if(jpos == uj-2)then
              wf = .true.
              jpos = jpos + 1
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              ipos = ipos - 1
              if(ipos == 1)then
                if(.NOT. polydone)CALL add(V,ipos,jpos,polydone,done)
                ele(ipos,jpos)%used = .true.
                jpos = jpos - 1
                dir = 4
              else
                dir = 2
              endif
            else
              jpos = jpos + 1
              dir = 6
            endif
          
          else                                                  !if exit left

            CALL add(U,ipos,jpos,polydone,done)
            if(jpos == 2)then
              wf = .true.
              jpos = jpos - 1
              if(.NOT.polydone)CALL add(V,ipos+1,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              ipos = ipos + 1
              if(ipos == vi-1)then
                if(.NOT. polydone) CALL add(U,ipos,jpos+1,polydone,done)
                ele(ipos,jpos)%used = .true.
                jpos = jpos + 1
                dir = 6
              else
                dir = 8
              endif
            else
              jpos = jpos - 1
              dir = 4
            endif
          
          endif


        CASE(18,19)                         !These forms exit up or down
                                            !18: r->u & l->d  19: l->u & r->d

          if((ele(ipos,jpos)%form == 18 .AND. dir == 4).OR.   & !if exit top
             (ele(ipos,jpos)%form == 19 .AND. dir == 6)) then

            CALL add(V,ipos+1,jpos,polydone,done)
            if(ipos == vi-2)then
              wf = .true.
              ipos = ipos + 1
              if(.NOT.polydone)CALL add(U,ipos,jpos+1,polydone,done)
              ele(ipos,jpos)%used = .true.
              jpos = jpos + 1
              if(jpos == uj-1)then
                if(.NOT. polydone)CALL add(V,ipos,jpos,polydone,done)
                ele(ipos,jpos)%used = .true.
                ipos = ipos - 1
                dir = 2
              else
                dir = 6
              endif
            else
              ipos = ipos + 1
              dir = 8
            endif

          else                                                  !if exit down

            CALL add(V,ipos,jpos,polydone,done)
            if(ipos == 2)then
              wf = .true.
              ipos = ipos - 1
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done)
              ele(ipos,jpos)%used = .true.
              jpos = jpos - 1
              if(jpos == 1)then
                if(.NOT. polydone) CALL add(V,ipos+1,jpos,polydone,done)
                ele(ipos,jpos)%used = .true.
                ipos = ipos + 1
                dir = 8
              else
                dir = 4
              endif
            else
              ipos = ipos - 1
              dir = 2
            endif

          endif

      END SELECT
    endif

  enddo


        !*********************************
        !*         ISLAND SECTION        *
        !*********************************

  do
    if(done)exit                  !if not all boundary points have been used
                                  !  Islands must exist...


        !*********************************
        !*  FIND ELEMENT TO START ISLAND *
        !*********************************

    found = .false.               !initialize found to false

    do i=2,vi-2
      do j=2,uj-2
        if(ele(i,j)%form < 15 .AND. (ele(i,j)%used .EQV. .false.))then
          ipos = i                !find an element that hasn't been used yet
          jpos = j                !store its location in ipos,jpos
          found = .true.
          exit
        endif
      enddo
      if(found)exit               !if one has been found, exit
    enddo

        !*********************************
        !*    CREATE ISLAND BOUNDARIES   *
        !*********************************


    polydone = .false.            !initialize polydone to FALSE

    do

      if(polydone)exit            !if the current island is finished, exit


      if(ele(ipos,jpos)%form > 15 .AND. ele(ipos,jpos)%unused)then
        !if the element is a cross and has not been used, set unused to FALSE
        !  this indicates that it has been used once
        ele(ipos,jpos)%unused = .false.
      else
        !if the element is NOT a cross, or if it is & has been used once prior
        !  set used to TRUE to indicate that the element will never need to be
        !  used again
        ele(ipos,jpos)%used = .true.
      endif

      !NOTE: there is no wf section this time because islands cannot be on
      !  the edges of the grid... land on the edges of the grid would be
      !  already used in the previous section

      SELECT CASE(ele(ipos,jpos)%form)

        CASE(1,3,11)              !These forms exit element right
          CALL add(U,ipos,jpos+1,polydone,done)
          jpos = jpos + 1
          dir = 6

        CASE(2,10,14)             !These forms exit element bottom
          CALL add(V,ipos,jpos,polydone,done)
          ipos = ipos - 1
          dir = 2

        CASE(4,5,7)               !These forms exit element top
          CALL add(V,ipos+1,jpos,polydone,done)
          ipos = ipos + 1
          dir = 8

        CASE(8,12,13)             !These forms exit element left
          CALL add(U,ipos,jpos,polydone,done)
          jpos = jpos - 1
          dir = 4

        CASE(16,17)               !These forms exit up or down
                                  !16: r->u & l->d    17: l->u & r->d

          if((ele(ipos,jpos)%form == 16 .AND. dir == 4) .OR.  & !if exit up
             (ele(ipos,jpos)%form == 17 .AND. dir == 6)) then
            CALL add(V,ipos+1,jpos,polydone,done)
            ipos = ipos + 1
            dir = 8
          else                                                  !if exit down
            CALL add(V,ipos,jpos,polydone,done)
            ipos = ipos - 1
            dir = 2
          endif

        CASE(18,19)               !These forms exit right or left
                                  !18: u->r & d->l    19: d->r & u->l

          if((ele(ipos,jpos)%form == 18 .AND. dir == 2) .OR.  & !if exit right
             (ele(ipos,jpos)%form == 19 .AND. dir == 8)) then
            CALL add(U,ipos,jpos+1,polydone,done)
            jpos = jpos + 1
            dir = 6           
          else                                                  !if exit left
            CALL add(U,ipos,jpos,polydone,done)
            jpos = jpos - 1
            dir = 4           
          endif

      END SELECT

    enddo

  enddo

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ~              Create Boundary Blanking Files for Surfer                ~
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(BoundaryBLNs) then
  CALL output_llBounds()
  CALL output_xyBounds()
endif

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ~               Prepare Boundary and Island Coordinates                 ~
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Set up ROMS boundaries. Use location of u and v grid points to demarcate the 
! land boundary between rho points where mask_rho = 0 and mask_rho = 1.

  !get the x/y locations of the U and V grid nodes
  CALL getUVxy(x_u,y_u,x_v,y_v)

  numpoly = bnds(bnum)%poly       ! # of polygons
  numislands = numpoly - 1        ! # of islands

  !allocate polysizes to contain the # of boundary points in each polygon
  ALLOCATE(polysizes(numpoly),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate polysizes'

  polysizes = 0    !initialize to 0

  !iterate through boundary points, counting the # in each polygon
  do k=1,bnum
    polysizes(bnds(k)%poly) = polysizes(bnds(k)%poly) + 1
  enddo

  maxbound = polysizes(1) + 1    !set # of boundary points around main boundary
  ALLOCATE(bx(maxbound))         !  and allocate bx and by to this #
  ALLOCATE(by(maxbound))         !the +1 is to close the polygon
  ALLOCATE(oo(maxbound))         !Track open ocean boundary points

  !initialize oo
  oo = .FALSE.

  !initialize # of island boundaries to the # of islands
  !  this takes care of +1 for each island (see comment above about +1)
  maxisland = numislands          
  if(numislands>0)then
    !Sum up the total # of island boundary points
    do k=1,numislands
      maxisland = maxisland + polysizes(1+k)
    enddo

    !  and allocate hid,hx,hy to the # of island boundary points
    ALLOCATE(hid(maxisland))
    ALLOCATE(hx(maxisland))
    ALLOCATE(hy(maxisland))
  endif

  pstart = 1

  !iterate through the polygons adding the boundary points to bx,by,hid,hx,hy
  do m=1,numpoly
    pend = pstart + polysizes(m) - 1
    if(m<3) c = 0
    do k=pstart,pend
      c = c + 1
      i = bnds(k)%ii
      j = bnds(k)%jj
      if(m == 1) then             !if it is the main body of water boundary
        if( bnds(k)%onU ) then    !  if its on the U grid
          bx(c) = x_u(i,j)        !  add the U grid node location to bx and by
          by(c) = y_u(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (i==1   ) .AND. mask_rho( 1,j)==1) &
            .OR. ( (i==vi-1) .AND. mask_rho(vi,j)==1) ) oo(c) = .TRUE.
        else                      !  if its on the V grid
          bx(c) = x_v(i,j)        !  add the V grid node location to bx and by
          by(c) = y_v(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (j==1   ) .AND. mask_rho(i, 1)==1) &
            .OR. ( (j==uj-1) .AND. mask_rho(i,uj)==1) ) oo(c) = .TRUE.
        endif
      else                        !if it is an island boundary
        hid(c) = m+999           !  add island # to hid
                                  !  note: island IDs start at 1001
        if( bnds(k)%onU ) then    !  if its on the U grid
          hx(c) = x_u(i,j)        !  add the U grid node location to hx and hy
          hy(c) = y_u(i,j)
        else                      !  if its on the V grid
          hx(c) = x_v(i,j)        !  add the V grid node location to hx and hy
          hy(c) = y_v(i,j)
        endif
      endif
    enddo

    !At the end of each polygon, put the first boundary point on the end to 
    !  close the polygon
    c = c + 1
    k = pstart
    i = bnds(k)%ii
    j = bnds(k)%jj
    if(m == 1)then                !if its mainbay bound
      if( bnds(k)%onU ) then
          bx(c) = x_u(i,j)
          by(c) = y_u(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (i==1   ) .AND. mask_rho( 1,j)==1) &
            .OR. ( (i==vi-1) .AND. mask_rho(vi,j)==1) ) oo(c) = .TRUE.
      else
          bx(c) = x_v(i,j)
          by(c) = y_v(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (j==1   ) .AND. mask_rho(i, 1)==1) &
            .OR. ( (j==uj-1) .AND. mask_rho(i,uj)==1) ) oo(c) = .TRUE.
      endif
    else                          !if its island bound
      if( bnds(k)%onU ) then
          hid(c) = m+999
          hx(c) = x_u(i,j)
          hy(c) = y_u(i,j)
      else
          hid(c) = m+999
          hx(c) = x_v(i,j)
          hy(c) = y_v(i,j)
      endif
    endif

    pstart = pstart + polysizes(m)

  enddo


  ! Construct bnd_x and bnd_y matrices that contain boundary line segments.
  ! bnd_x(1,i),bnd_y(1,i) = coordinates of first point in line segment.
  ! bnd_x(2,i),bnd_y(2,i) = coordinates of second point in line segment.
  ! The bnd_x/bnd_y matrices are used to 
  !   1) determine if a particle intersects a boundary segment and 
  !   2) reflect the particle off the boundary segment.
  nbounds = maxbound + maxisland - numislands - 1
  ALLOCATE(bnd_x(2,nbounds))
  ALLOCATE(bnd_y(2,nbounds))


  !Allocate land to track land/ocean boundaries
  ALLOCATE (land(nbounds),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating land'

  !Initialize newly allocated land
  !TRUE if the boundary is land, FALSE if it is open ocean
  land = .TRUE.

  if(BoundaryBLNs) then
    OPEN(10,FILE='OpenOceanBoundaryMidpoints.csv',STATUS='REPLACE')
  endif

  ! First create the line segments for the land and open ocean boundaries by
  !   simply writing each x/y location and the next x/y location. The fact 
  !   that the first bx/by coordinate equals the last ensures that all line 
  !   segments surrounding the model are included.
  do i=1,maxbound-1
    bnd_x(1,i) = bx(i)   
    bnd_y(1,i) = by(i)
    bnd_x(2,i) = bx(i+1)
    bnd_y(2,i) = by(i+1)
    if(oo(i).AND.oo(i+1))then
      land(i) = .FALSE. !Not a land boundary
      if(BoundaryBLNs) WRITE(10,*) (bx(i)+bx(i+1))/2,',',(by(i)+by(i+1))/2
    endif
  enddo

  if(BoundaryBLNs) CLOSE(10)

  if(numislands>0)then
    ! Now create the line segments for the island boundaries. Same principle as
    !   above, but need to account for the fact that there are islands in the
    !   hx/hy matrices. Use hid matrix of island id numbers to guide process. 
    j=1                         
    k=hid(1)
    do i=1,maxisland-1
      if(hid(i) .NE. k)then       !if moving on to a new island
        k = hid(i)                !update k to new island id
        j = j + 1                 !and increment j
      endif
      bnd_x(1,maxbound+i-j) = hx(i)
      bnd_y(1,maxbound+i-j) = hy(i)
      bnd_x(2,maxbound+i-j) = hx(i+1)
      bnd_y(2,maxbound+i-j) = hy(i+1)
    enddo
  endif


  DEALLOCATE(bnds,polysizes)
  DEALLOCATE(x_u)
  DEALLOCATE(y_u)
  DEALLOCATE(x_v)
  DEALLOCATE(y_v)
  DEALLOCATE(mask_rho)
  DEALLOCATE(oo)

END SUBROUTINE createBounds



    !*********************************************************
    !*                         add                           *
    !*********************************************************

SUBROUTINE add(isU,iii,jjj,polydone,done)
!This subroutine is for adding boundary points to bnds
! INPUT:
!   isU - TRUE if point is on U grid, FALSE if it is on V grid
!   iii - i position on the grid
!   jjj - j position on the grid

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iii,jjj
  LOGICAL, INTENT(IN) :: isU
  LOGICAL, INTENT(INOUT) :: polydone,done

  INTEGER :: polystart = 1  !bnd location at start of current polygon
  INTEGER :: polynum   = 1  !# of the current polygon
  INTEGER :: b = 0          !# of boundary points used so far
  !Note that any local variable that is initialized in a type declaration
  !  statement, as done above, is automatically given the SAVE attribute,
  !  this means that the value in these variables will be retained in 
  !  subsequent calls to the subroutine

  !Check if the current polygon has returned to its starting location
! CRS commented out next two lines
!  if( (polystart /= (b+1)) .AND. (bnds(polystart)%onU == isU) .AND.       &
!      (bnds(polystart)%ii == iii).AND.(bnds(polystart)%jj == jjj) )then
! CRS replaced them with:
  if( (polystart .NE. (b+1)) .AND. (bnds(polystart)%onU .EQV. isU) .AND. &
      (bnds(polystart)%ii == iii).AND.(bnds(polystart)%jj == jjj) )then

    !switch polydone to TRUE to indicate the current polygon is complete
    polydone = .true.

    !since the current polygon is complete, increment polynum 
    polynum = polynum + 1
    !and set polystart to the first position of the subsequent polygon
    polystart = b + 1

    if(b == bnum)then
      done = .true.
      BND_SET = .true.
    endif

  else

    !If the current polygon isn't complete
    b = b+1             !increment b to the next boundary point and
    bnds(b)%onU  = isU  !store values for onU, ii, jj, and polynum
    bnds(b)%ii   = iii
    bnds(b)%jj   = jjj
    bnds(b)%poly = polynum

  endif

END SUBROUTINE add



    !*********************************************************
    !*                       getNext                         *
    !*********************************************************

SUBROUTINE getNext(i,j,wf,dir,form)
  USE PARAM_MOD, ONLY: vi,uj
  !This subroutine is for finding the next element, following the bounds
  !  clockwise around water
  !This is all based on the element location, whether or not its on an 
  !  edge, the direction the path is coming from, and the element form
  IMPLICIT NONE

  INTEGER :: i,j,dir,form
  LOGICAL :: wf

  if(wf)then                      !If the element is on the edge
    if(dir == 2) then             !If path is going down,
      SELECT CASE(form)           !  then its on right edge

        CASE(0,1,4,5)             !If path stays on edge
          i = i - 1
          if(i == 2)then          !If path reaches bottom,
            dir = 4               !  then go left on bottom edge
            j = j - 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    elseif(dir == 4) then         !If path is going left,
      SELECT CASE(form)           !  then its on bottom edge

        CASE(0,1,2,3)             !If path stays on edge
          j = j - 1
          if(j == 2)then          !If path reaches left edge,
            dir = 8               !  then go up on left edge
            i = i + 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    elseif(dir == 6) then         !If path is going right,
      SELECT CASE(form)           !  then its on top edge

        CASE(0,4,8,12)            !If path stays on edge
          j = j + 1
          if(j == uj-2)then       !If path reaches right edge,
            dir = 2               !  then go down on right edge
            i = i - 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    elseif(dir == 8) then         !If path is going up,
      SELECT CASE(form)           !  then its on left edge

        CASE(0,2,8,10)            !If path stays on edge
          i = i + 1
          if(i == vi-2)then       !If path reaches top edge,
            dir = 6               !  then go right on top edge
            j = j + 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    else
      write(*,*) 'Error: wf direction not one of 2,4,6,8'
    endif





  else                            !If element is not on edge
    SELECT CASE(form)

    CASE(1,5,13)                  !These forms exit element bottom
      if(i == 2)then              !  if moving into bottom edge
        wf = .true.               !    wf is true
        i = i - 1                 !    move down & left
        j = j - 1
        if(j == 1)then            !  if moving to down left corner
          i = i + 1               !    move up
          dir = 8
        else
          dir = 4
        endif
      else                        !  if not moving to bottom edge
        i = i - 1                 !    just move down
        dir = 2
      endif

    CASE(2,3,7)                   !These forms exit element left
      if(j == 2)then              !  if moving into left edge
        wf = .true.               !    wf is true
        j = j - 1                 !    move left and up
        i = i + 1
        if(i == vi-1)then         !  if moving to top left corner
          j = j + 1               !    move right
          dir = 6
        else
          dir = 8
        endif
      else                        !  if not moving to left edge
        j = j - 1                 !    just move left
        dir = 4
      endif

    CASE(4,12,14)                 !These forms exit element right
      if(j == uj-2)then           !  if moving to right edge
        wf = .true.               !    wf is true
        j = j + 1                 !    move right and down
        i = i - 1
        if(i == 1)then            !  moving to down right corner
          j = j - 1               !    move left
          dir = 4
        else
          dir = 2
        endif
      else                        !  if not moving to right edge
        j = j + 1                 !    just move right
        dir = 6
      endif

    CASE(8,10,11)                 !These forms exit element top
      if(i == vi-2)then           !  if moving to top edge
        wf = .true.               !    wf is true
        i = i + 1                 !    move up and right
        j = j + 1
        if(j == uj-1)then         !  if moving to top right corner
          i = i - 1               !    move down
          dir = 2
        else
          dir = 6
        endif
      else                        !  if not moving to top edge
        i = i + 1                 !    just move up
        dir = 8
      endif

    CASE(16,17)                   !These forms exit right or left
                                  !16: u->r & d->l    17: d->r & u->l

                                  !if exit right
      if((form == 16 .AND. dir == 2).OR. (form == 17 .AND. dir == 8))then

        if(j == uj-2)then         !  if moving to right edge
          wf = .true.             !    wf is true
          j = j + 1               !    move right and down
          i = i - 1
          if(i == 1)then          !  moving to down right corner
            j = j - 1             !    move left
            dir = 4
          else
            dir = 2
          endif
        else                      !  if not moving to right edge
          j = j + 1               !    just move right
          dir = 6
        endif
      
      else                        !if exit left

        if(j == 2)then            !  if moving to left edge
          wf = .true.             !    wf is true
          j = j - 1               !    move left and up
          i = i + 1
          if(i == vi-1)then       !  if moving to top left corner
            j = j + 1             !    move right
            dir = 6
          else
            dir = 8
          endif
        else                      !  if not moving to left edge
          j = j - 1               !    just move left
          dir = 4
        endif
      
      endif

    CASE(18,19)                   !These forms exit up or down
                                  !18: r->u & l->d    19: l->u & r->d

                                  !if exit up
      if((form == 18 .AND. dir == 4).OR. (form == 19 .AND. dir == 6))then

        if(i == vi-2)then         !  if moving to top edge
          wf = .true.             !    wf is true
          i = i + 1               !    move up and right
          j = j + 1
          if(j == uj-1)then       !  if moving to top right corner
            i = i - 1             !    move down
            dir = 2
          else
            dir = 6
          endif
        else                      !  if not moving to top edge
          i = i + 1               !    just move up
          dir = 8
        endif

      else                        !if exit down

        if(i == 2)then            !  if moving to bottom edge
          wf = .true.             !    wf is true
          i = i - 1               !    move down and left
          j = j - 1
          if(j == 1)then          !  if moving to down left corner
            i = i + 1             !    move up
            dir = 8
          else
            dir = 4
          endif
        else                      !  if not moving to bottom edge
          i = i - 1               !    just move down
          dir = 2
        endif

      endif

    END SELECT
  endif

END SUBROUTINE getNext

  ! This subroutine uses the point-in-polygon approach to determine if a
  !   centroid is inside the model domain 
  ! It returns the value 1 in the variable inbounds if the centroid is in
  !   bounds and 0 if it is outside the boundaries
  SUBROUTINE mbounds(Ypos,Xpos,inbounds)

    USE PIP_MOD, ONLY: inpoly
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: inbounds
    DOUBLE PRECISION, INTENT(IN) :: Ypos,Xpos

    INTEGER :: i
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: blatlon

    ALLOCATE(blatlon(maxbound,2))
    
    !read main boundaries into one variable for inpoly
    do i=1,maxbound
      blatlon(i,1) = bx(i)
      blatlon(i,2) = by(i)
    enddo
      
    inbounds = 0    !Initialize inbounds to 0
    !if the point is in the boundaries, change inbounds to 1
    if(inpoly(Xpos,Ypos,maxbound,blatlon)) inbounds = 1

    DEALLOCATE(blatlon)

  END SUBROUTINE mbounds



  ! This subroutine uses the point-in-polygon approach to determine if a
  ! centriod is inside one of the islands within the model domain 
  ! It returns the value 1 in the variable in_island if the centroid is in 
  !   island boundaries and 0 if it is outside the island boundaries
  SUBROUTINE ibounds(in_island,claty,clongx,island)

    USE PIP_MOD, ONLY: inpoly
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: in_island
    DOUBLE PRECISION, INTENT(IN) :: claty,clongx
    DOUBLE PRECISION, INTENT(OUT) :: island

    LOGICAL :: endIsle
    INTEGER :: i,j,start,count
    DOUBLE PRECISION :: isle
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: isbnds

    in_island = 0            !initialize in_island to 0
    island = 0.0             !initialize island to 0.0

    if(numislands>0)then

      i = 1                  !initialize i to 1
      isle = hid(i)          !initialize isle to first island id
      start = 0              !initialize start to 0
      do

        i = i + 1            !iterate through island boundary points

        !Determine if the end of an island has been reached by checking if it
        !  is either the very last point, or if not, if the next point is not
        !  part of the same island
        endIsle = .FALSE.
        if(i == maxisland)then
          endIsle = .TRUE.
        else
          if(hid(i+1) /= isle) endIsle = .TRUE.
        endif

        !if the end of an island has been reached:
        !  set count to the number of points around the island, allocate
        !  isbnds, then iterate through the island edge points storing
        !  the boundary locations into isbnds
        !  Next, call inpoly to see if the point is in the island
        !    If so, set island to the id of the island it is in, deallocate
        !      isbnds, and exit the subroutine
        !    If not, deallocate isbnds and repeat for all remaining islands
        if(endIsle)then
          count = i - start
          ALLOCATE(isbnds(count,2))
          do j=1,count
            isbnds(j,1) = hx(start+j)
            isbnds(j,2) = hy(start+j)
          enddo
          if(inpoly(clongx,claty,count,isbnds))then
            in_island = 1
            island = isle
            DEALLOCATE(isbnds)
            exit
          endif
          DEALLOCATE(isbnds)
          if(i == maxisland)exit
          start = i
          isle = hid(i+1)
        endif
      enddo

    endif

  END SUBROUTINE ibounds


  ! This subroutine calculates the intersection between the particle
  ! trajectory and the boundary line in a grid cell, and then calculates
  ! the reflection, returning the new particle location
  subroutine intersect_reflect(Xpos,Ypos,nXpos,nYpos,fintersectX,fintersectY,  &
    freflectX,freflectY,intersectf,skipbound,isWater)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: intersectf
    INTEGER, INTENT(INOUT) :: skipbound
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,nXpos,nYpos
    DOUBLE PRECISION, INTENT(OUT) :: fintersectX,fintersectY,freflectX,freflectY
    LOGICAL, OPTIONAL, INTENT(OUT) :: isWater
    INTEGER :: i,intersect,skipboundi
    DOUBLE PRECISION :: crossk,dPBC,mBCperp,rx1,rx2,ry1,ry2,Bp,distBC,dist1,   &
      dist2,intersctx,interscty,rPxyzX,rPxyzY,Mbc,Bbc,Mp,bcx1,bcy1,bcx2,bcy2,  &
      bBCperp,xhigh,xlow,yhigh,ylow,d_Pinter,dtest,bxhigh,bxlow,byhigh,bylow

    distBC=0.0
    Mbc = 0.0
    Bbc = 0.0
    Mp = 0.0
    Bp = 0.0
    intersect=0
    intersectf=0
    skipboundi = skipbound
    fintersectX = -999999.
    fintersectY = -999999.
    freflectX = -999999.
    freflectY = -999999.
    dtest = 999999.
    isWater = .FALSE.

    if (Xpos.GE.nXpos) then
      xhigh = Xpos
      xlow = nXpos
    else
      xhigh = nXpos
      xlow = Xpos
    endif  

    if (Ypos.GE.nYpos) then
      yhigh = Ypos
      ylow = nYpos
    else
      yhigh = nYpos
      ylow = Ypos
    endif

    do i=1,nbounds

      if (i == skipbound) cycle

        intersect = 0
        bcx1=bnd_x(1,i)
        bcy1=bnd_y(1,i)
        bcx2=bnd_x(2,i)
        bcy2=bnd_y(2,i)

        !If the boundary segment end points are both east, west, north, or 
        !  south of the particle's previous or new location, cycle to next 
        !  boundary
        if( ((bcx1 > xhigh) .AND. (bcx2 > xhigh)) .OR. &
            ((bcx1 < xlow ) .AND. (bcx2 < xlow )) .OR. &
            ((bcy1 > yhigh) .AND. (bcy2 > yhigh)) .OR. &
            ((bcy1 < ylow ) .AND. (bcy2 < ylow ))      ) cycle
        
        if (bcx1.GE.bcx2) then
          bxhigh = bcx1
          bxlow = bcx2
        else
          bxhigh = bcx2
          bxlow = bcx1
        endif  

        if (bcy1.GE.bcy2) then
          byhigh = bcy1
          bylow = bcy2
        else
          byhigh = bcy2
          bylow = bcy1
        endif

        !First determine if an undefined denominator is possible
        if (bcx1.EQ.bcx2 .OR. nXpos.EQ.Xpos ) then
          !test if they both vertical, if so cycle because they cannot intersect
          if (bcx1.EQ.bcx2 .AND. nXpos.EQ.Xpos ) cycle
          !test if perpendicular and parrallel to coordinate axes
          if (bcx1.EQ.bcx2 .AND. nYpos.EQ.Ypos ) then
            !undefined denominator, perp. & || to axes
            intersctx = bcx1
            interscty = nYpos
            if (intersctx.LE.xhigh  .AND. intersctx.GE.xlow .AND.              &
              interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
              intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
              interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
              rx1=nXpos+(DBLE(2.0)*dPBC)
              ry1=nYpos
              rx2=nXpos-(DBLE(2.0)*dPBC)
              ry2=nYpos
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              elseif(dist1.GT.dist2) then
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          elseif (nXpos.EQ.Xpos .AND. bcy1.EQ.bcy2 ) then
            !undefined denominator, perp. & || to axes
            intersctx = nXpos
            interscty = bcy1
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
              interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
              intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
              interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
              rx1=nXpos
              ry1=nYpos+(DBLE(2.0)*dPBC)
              rx2=nXpos
              ry2=nYpos-(DBLE(2.0)*dPBC)
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              elseif(dist1.GT.dist2) then
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          elseif (bcx1.EQ.bcx2 .AND. nYpos.NE.Ypos ) then
            !undefined denominator, not perpendicular
            Mp = (nYpos-Ypos)/(nXpos-Xpos)
            Bp = Ypos - Mp*Xpos
            intersctx = bcx1
            interscty = Mp*intersctx + Bp
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              dPBC = nXpos-intersctx
              rx1=nXpos+(DBLE(2.0)*dPBC)
              ry1=nYpos
              rx2=nXpos-(DBLE(2.0)*dPBC)
              ry2=nYpos
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              elseif(dist1.GT.dist2) then
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          elseif (nXpos.EQ.Xpos .AND. bcy1.NE.bcy2  ) then
            !undefined denominator, not perpendicular
            Mbc = (bcy2-bcy1)/(bcx2-bcx1)
            Bbc = bcy2 - Mbc*bcx2
            intersctx = nXpos
            interscty = Mbc*intersctx + Bbc
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                interscty.LE.byhigh .AND. interscty.GE.bylow ) then
              !Now use cross product to determine the distance of the particle 
              !  from the boundary
              distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
              crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
              dPBC = sqrt(crossk**2)/distBC
              !find line perpendicular to boundary
              mBCperp = DBLE(-1.0)/Mbc
              bBCperp = nYpos - mBCperp*nXpos
              !find two potential reflection points
              rx1 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) ) +nXpos
              ry1 = mBCperp*rx1 + bBCperp
              rx2 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) )       &
                  * DBLE(-1.0) + nXpos
              ry2 = mBCperp*rx2 + bBCperp
              !point closest to intersection of boundary and particle trajectory
              !  is the right one
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              elseif(dist1.GT.dist2) then
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          endif
        else

          if(intersect == 0)then

            Mbc = (bcy2-bcy1)/(bcx2-bcx1)
            Bbc = bcy2 - Mbc*bcx2
            Mp = (nYpos-Ypos)/(nXpos-Xpos)
            Bp = Ypos - Mp*Xpos
            intersctx = (Bbc - Bp)/(Mp - Mbc)
            interscty = Mp*intersctx + Bp

            !when bc parallel with x-axis, byhigh=bylow=intersecty
            if (Mbc.EQ.0.0) interscty = byhigh
        
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                interscty.LE.byhigh .AND. interscty.GE.bylow  ) then

              if (Mbc.EQ.0.0) then  !inverse slope denominator not OK
                dPBC = nYpos-bcy1
                rx1=nXpos
                ry1=nYpos+(DBLE(2.0)*dPBC)
                rx2=nXpos
                ry2=nYpos-(DBLE(2.0)*dPBC)
                dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
                dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 ) 
                if(dist1.LT.dist2) then
                  rPxyzX= rx1
                  rPxyzY= ry1
                elseif(dist1.GT.dist2) then
                  rPxyzX= rx2
                  rPxyzY= ry2
                endif
                  intersect=1
                endif  

              if(intersect == 0)then

                !Now use cross product to determine the distance of the
                !  particle from the boundary
                distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
                crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
                dPBC = sqrt(crossk**2)/distBC
                !find line perpendicular to boundary
                mBCperp = DBLE(-1.0)/Mbc
                bBCperp = nYpos - mBCperp*nXpos
                !find two potential reflection points
                rx1 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2)) +nXpos
                ry1 = mBCperp*rx1 + bBCperp
                rx2 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2))       &
                    * DBLE(-1.0) + nXpos
                ry2 = mBCperp*rx2 + bBCperp
                !point closest to intersection of boundary and particle 
                !  trajectory is the right one
                dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
                dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
                if(dist1.LT.dist2) then
                  rPxyzX= rx1
                  rPxyzY= ry1
                elseif(dist1.GT.dist2) then
                  rPxyzX= rx2
                  rPxyzY= ry2
                endif
                intersect=1
                
              endif
            endif  
          endif
        endif


        d_Pinter = sqrt( (Xpos-intersctx)**2 + (Ypos-interscty)**2 )
        if( (intersect .EQ. 1) .AND. (d_Pinter .LT. dtest) ) then
          fintersectX = intersctx
          fintersectY = interscty
          freflectX = rPxyzX
          freflectY = rPxyzY
          intersectf = 1
          dtest = d_Pinter
          skipboundi = i
          isWater = .NOT. land(i)
        endif

    enddo

    skipbound = skipboundi
  END SUBROUTINE intersect_reflect


SUBROUTINE output_xyBounds()
!This subroutine creates a blanking file for surfer of the model 
!  boundaries in meters
  USE PARAM_MOD, ONLY: ui,uj,vi,vj
  USE HYDRO_MOD, ONLY: getUVxy
  IMPLICIT NONE

  INTEGER :: i,j,k,m,numpoly,STATUS,pstart,pend
  INTEGER, ALLOCATABLE, DIMENSION(:) :: polysizes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_u,y_u,x_v,y_v

  ALLOCATE(x_u(ui,uj))
  ALLOCATE(y_u(ui,uj))
  ALLOCATE(x_v(vi,vj))
  ALLOCATE(y_v(vi,vj))

  write(*,*) 'output metric blanking file'

  !get the x/y locations of the U and V grid nodes
  CALL getUVxy(x_u,y_u,x_v,y_v)

  open(1,FILE='xybounds.bln',STATUS='REPLACE')

    numpoly = bnds(bnum)%poly       !numpoly = # of polygons

    !Allocate polysizes to the number of polygons
    allocate(polysizes(numpoly),STAT=STATUS)
    IF(STATUS /= 0) write(*,*) 'Problem allocate polysizes'

    polysizes = 0       !initialize all polygon sizes to 0

    !iterate through boundaries incrementing polysize depending on
    !   polygon each boundary is part of
    do k=1,bnum
      polysizes(bnds(k)%poly) = polysizes(bnds(k)%poly) + 1
    enddo

    pstart = 1

    do m=1,numpoly                          !iterate through the polygons
      write(1,*) polysizes(m)+1,',',1       !write polysize,1 for .bln 1st row
      pend = pstart + polysizes(m) - 1      !calculate last bnds location
      do k=pstart,pend                      !iterate through polygon bnds
        i = bnds(k)%ii                      !i position
        j = bnds(k)%jj                      !j position
        if( bnds(k)%onU ) then              !if on U Grid
          write(1,*) x_u(i,j),',',y_u(i,j)  !  write x & y of U node
        else                                !else
          write(1,*) x_v(i,j),',',y_v(i,j)  !  write x & y of V node
        endif
      enddo
      k = pstart                            !go back to polygon's first bnd
      i = bnds(k)%ii                        !  to close the polygon
      j = bnds(k)%jj
      if( bnds(k)%onU ) then
        write(1,*) x_u(i,j),',',y_u(i,j)
      else
        write(1,*) x_v(i,j),',',y_v(i,j)
      endif

      pstart = pstart + polysizes(m)        !calculate start of next polygon

    enddo

  close(1)

  DEALLOCATE(x_u,y_u,x_v,y_v)

END SUBROUTINE output_xyBounds


SUBROUTINE output_llBounds()
!This subroutine creates a blanking file for Surfer/Scripter of the model 
!  boundaries in longitude and latitude
  USE PARAM_MOD,   ONLY: ui,uj,vi,vj
  USE HYDRO_MOD,   ONLY: getUVxy
  USE CONVERT_MOD, ONLY: x2lon,y2lat
  IMPLICIT NONE

  INTEGER :: i,j,k,m,numpoly,STATUS,pstart,pend
  INTEGER, ALLOCATABLE, DIMENSION(:) :: polysizes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_u,y_u,x_v,y_v

  ALLOCATE(x_u(ui,uj))
  ALLOCATE(y_u(ui,uj))
  ALLOCATE(x_v(vi,vj))
  ALLOCATE(y_v(vi,vj))

  write(*,*) 'output lat/long blanking file'

  !get the x/y locations of the U and V grid nodes
  CALL getUVxy(x_u,y_u,x_v,y_v)

  open(1,FILE='llbounds.bln',STATUS='REPLACE')

    numpoly = bnds(bnum)%poly       !numpoly = # of polygons

    !Allocate polysizes to the number of polygons
    allocate(polysizes(numpoly),STAT=STATUS)
    IF(STATUS /= 0) write(*,*) 'Problem allocate polysizes'

    polysizes = 0       !initialize all polygon sizes to 0

    !iterate through boundaries incrementing polysize depending on
    !   polygon each boundary is part of
    do k=1,bnum
      polysizes(bnds(k)%poly) = polysizes(bnds(k)%poly) + 1
    enddo

    pstart = 1

    do m=1,numpoly                      !iterate through the polygons
      write(1,*) polysizes(m)+1,',',1   !write polysize,1 for .bln first row
      pend = pstart + polysizes(m) - 1  !calculate polygon last bnds location
      do k=pstart,pend                  !iterate through polygon bnds
        i = bnds(k)%ii                  !i position
        j = bnds(k)%jj                  !j position
        if( bnds(k)%onU ) then          !if on U Grid write lon & lat of U node
          write(1,*) x2lon(x_u(i,j),y_u(i,j)),',',y2lat(y_u(i,j))
        else                            !else write lon & lat of V node
          write(1,*) x2lon(x_v(i,j),y_v(i,j)),',',y2lat(y_v(i,j))
        endif
      enddo
      k = pstart                        !go back to polygon's first bnd
      i = bnds(k)%ii                    !  to close the polygon
      j = bnds(k)%jj
      if( bnds(k)%onU ) then
        write(1,*) x2lon(x_u(i,j),y_u(i,j)),',',y2lat(y_u(i,j))
      else
        write(1,*) x2lon(x_v(i,j),y_v(i,j)),',',y2lat(y_v(i,j))
      endif

      pstart = pstart + polysizes(m)    !calculate start of next polygon

    enddo

  close(1)

  DEALLOCATE(x_u,y_u,x_v,y_v)

END SUBROUTINE output_llBounds

END MODULE

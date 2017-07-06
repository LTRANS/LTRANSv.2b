MODULE BEHAVIOR_MOD

! The behavior module is used to assign biological or physical characteristics
!   to particles. 
! Currently particle movement is in the vertical direction. 
!
! Particle characteristics can include a swimming/sinking speed component and 
! a behavioral cue component that can depend upon particle age. The swimming/
! sinking speed component controls the speed of particle motion and can be 
! constant or set with a function. The behavioral cue component regulates the 
! direction of particle movement. For biological behaviors, a random component 
! is added to the swimming speed and direction to simulate random variation in
! the movements of individuals (in behavior types 1 - 5, see list below). 
! Physical characteristics can also be assigned to particles, like constant 
! sinking velocity, without the additional random movements (behavior type 6). 
! The following behavior types are currently available in LTRANS and are 
! specified using the Behavior parameter in the LTRANS.inc file:
!
!
! Passive (no behavior): Behavior = 0. In this case, the behavior module is not
!   executed. Particle motion is based on advection, and, if turned on, 
!   horizontal and vertical turbulence.
!
! Near-surface orientation: Behavior = 1. Particles swim up if they are deeper 
!   than 1 m from the surface.  
!
! Near-bottom orientation: Behavior = 2. Particles swim down if they are 
!   shallower than 1 m from the bottom.  
!
! Diurnal vertical migration: Behavior = 3. Particles swim down if light levels 
!   at the particle location exceed a predefined threshold value.  
!
! Crassostrea virginica oyster larvae: Behavior = 4. Swimming speeds and 
!   direction of motion vary depending upon age (stage) according to field and 
!   laboratory observations (see North et al. 2008). 
!
! C. ariakensis oyster larvae: Behavior = 5. Swimming speeds and direction of 
!   motion vary depending upon age (stage) according to field and laboratory 
!   observations (see North et al. 2008).
!
! Sinking velocity: Behavior = 6. Particles move up or down with constant
!   sinking (or floating) speeds without individual random motion. Code that 
!   calculates salinity and temperature at the particle location is included 
!   (but commented out) as a basis for calculating density-dependent sinking 
!   velocities.  
!
! Tidal Stream Transport: Behavior = 7.
!
!
! Behavior algorithms and code created by: Elizabeth North
! Module structure created by:             Zachary Schlag
! Created on:                              2004
! Last Modified on:                        22 March 2011
!

  IMPLICIT NONE
  PRIVATE
  SAVE

  !Timer for C. ariakensis downward swimming behavior
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: timer

  !Behavior of each particle
  INTEGER, ALLOCATABLE, DIMENSION(:) :: P_behave


  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:  ) ::    &
    P_pediage, & !Age at which the particle will settle (become a pediveliger)
    P_deadage, & !Age at which the particle will stop moving (die)
      !The following are for calculating salt gradient:
    P_Sprev,   & !Salinity at particle's previous location 
    P_zprev      !Particle's previous depth 


  !Swimming speed (age-dependent, linear increase unless constant)   
  !(n,1)slope, (n,2)intercept, (n,3) speed at current age
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: P_swim

  !For behavior 7, tracks if particle is on the bottom
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: bottom

  !Tracks if the particle is dead (TRUE) or alive (FALSE)
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: dead

  !Tracks if particles are Out Of Bounds (ie cross open ocean bound)
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: oob

  !The following procedures have been made public:
  PUBLIC :: initBehave,updateStatus,behave,getStatus,finBehave, &
            setOut,isOut,die,isDead

CONTAINS

  SUBROUTINE initBehave()    !Initialize the behavior module
    USE PARAM_MOD, ONLY: numpar,Behavior,swimfast,swimslow,swimstart,     &
                         pediage,deadage,Sgradient,settlementon 
    USE SETTLEMENT_MOD, ONLY: initSettlement
    USE NORM_MOD,   ONLY: norm
    IMPLICIT NONE
    INTEGER :: n

    write(*,*) 'initialize behavior'    

    !Allocate Behavior Variables
    ALLOCATE(timer(numpar))
    ALLOCATE(P_behave(numpar))
    ALLOCATE(P_pediage(numpar))
    ALLOCATE(P_deadage(numpar))
    ALLOCATE(P_Sprev(numpar))
    ALLOCATE(P_zprev(numpar))
    ALLOCATE(P_swim(numpar,3))
    ALLOCATE(dead(numpar))
    ALLOCATE(oob(numpar))

    IF(Behavior == 7) THEN
      ALLOCATE(bottom(numpar))
      bottom = .TRUE.
    ENDIF

    do n=1,numpar
      !Set behavior to the one specified in LTRANS.inc
      P_behave(n) = Behavior  !Behavior
      P_pediage(n) = pediage  !age at which particle reaches maximum swimming
                              !speed and can settle (becomes a pediveliger) (s)
      P_deadage(n) = deadage  !age at which particle stops moving (dies) (s)
      !Note: the following code assigns different veliger and pediveliger
      !  stage durations
      !P_pediage(n) = (14. + norm()*0.5)*24.*3600.
      !P_deadage(n) = P_pediage(n) + (7. + norm()*0.5)*24.*3600.

      !Calculate slope and intercept for age-dependent linear swimming speed
      P_swim(n,1) = (swimfast - swimslow)/(P_pediage(n) - swimstart) !slope
      P_swim(n,2) = swimfast - P_swim(n,1)*P_pediage(n)              !intercept
      P_swim(n,3) = 0.0                                   !swimming speed (m/s)
      !Note: P_swim(n,3) is updated at each time step in Subroutine behave
    enddo

    !The following variables are used by the C. virginica and C. ariakensis 
    !  behavior routines
    timer = DBLE(0.0)         !to count how long C. arikensis particles swim down

    ! Initialize salt storage matrices 
    P_Sprev = 0.0       !Initialized to 0.0
    P_zprev = 0.0       !Initialized to 0.0

    ! Initialize dead to .FALSE. i.e. all particles are initially alive
    dead = .FALSE.

    ! Initialize out of bounds tracker to .FALSE.
    !   (i.e. all particles start in bounds)
    oob = .FALSE.

    !if Settlement is turned on then inform Settlement module of the age at 
    !  which particle can settle (i.e., become pediveligers)
    if(settlementon)then
      CALL initSettlement(P_pediage)
    endif

  END SUBROUTINE initBehave


  SUBROUTINE updateStatus(P_age,n)  !Update particle status
    USE SETTLEMENT_MOD, ONLY: isSettled
    USE PARAM_MOD, ONLY: settlementon,mortality
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: P_age

    !Determine if particle dies from old age, if so kill it
    if ((P_age .GE. P_deadage(n)) .AND. mortality) then
      if(settlementon)then
        if(.NOT. isSettled(n)) call die(n)
      else
        call die(n)
      endif
    endif

  END SUBROUTINE updateStatus

  SUBROUTINE behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,P_zb,P_zc,P_zf,        &
                    P_zetac,P_age,P_depth,P_U,P_V,P_angle,n,it,ex,ix,          &
                    daytime,p,bott,XBehav,YBehav,ZBehav)
    USE PARAM_MOD, ONLY: us,dt,idt,twistart,twiend,Em,pi,daylength,Kd,thresh,  &
                         Sgradient,swimfast,swimstart,sink,Hswimspeed,         &
                         Swimdepth
    USE HYDRO_MOD, ONLY: WCTS_ITPI
    USE RANDOM_MOD, ONLY: genrand_real1
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: daytime
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,Pwc_zb(:),Pwc_zc(:),        &
                                    Pwc_zf(:),P_zb,P_zc,P_zf,P_zetac,P_age,    &
                                    P_depth,P_U,P_V,P_angle,ex(3),ix(3)
    INTEGER, INTENT(IN) :: n,it,p
    LOGICAL, INTENT(OUT) :: bott
    DOUBLE PRECISION, INTENT(OUT) :: XBehav,YBehav,ZBehav
    
    INTEGER :: btest,i,deplvl
    DOUBLE PRECISION :: negpos,dev1,devB,switch,switchslope
    DOUBLE PRECISION :: P_S,parBehav,Sslope,deltaS,deltaz
    ! DOUBLE PRECISION :: P_T !not needed unless temperature code below is enabled
    DOUBLE PRECISION :: dtime,tst,E0,P_light
    DOUBLE PRECISION :: currentspeed,Hdistance,theta,X,Y


    !   ***************** Initialize Return Values
    XBehav = 0.0
    YBehav = 0.0
    ZBehav = 0.0

    !   ***************** Update vertical swimming speeds based on particle age
    
    if(P_age .GE. swimstart) P_swim(n,3) = P_swim(n,1)*P_age+P_swim(n,2)
    if(P_age .GE. P_pediage(n)) P_swim(n,3) = swimfast
    
    !   ***************** Prepare for TYPE 4 & 5 (Oyster Larvae) Behaviors
    
    !Update pediveliger behavior/status and timer
    IF(P_behave(n) .EQ. 4 .OR. P_behave(n) .EQ. 5) THEN
       
      !Set behavior code for pediveligers
      if (P_age .GE. P_pediage(n) .AND. P_age .LT. P_deadage(n)) then
        P_behave(n) = 2
      endif
       
      !decrement timer
      timer(n) = max(DBLE(0.0), timer(n)-DBLE(dt))
    ENDIF
    
    !obtain salinity at particle location (P_S) to cue oyster larvae or tidal
    !  stream transport behavior
    IF ((P_behave(n).EQ.4) .OR. (P_behave(n).EQ.5 .AND. timer(n).EQ.0.0) .OR.  &
         P_behave(n).EQ.7) THEN 
       
      do i=3,us-2
        if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR.              &
            (Zpar .LT. Pwc_zf(i))) exit
      enddo
      deplvl = i-2   !depth level
       
      !Salinity at particle location
      P_S = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,    &
                      P_zc,P_zf,ex,ix,p,4)
       
    ENDIF

    
    !           *********************************************************
    !           *                                                       *
    !           *                   Behaviors                           *
    !           *                                                       *
    !           *********************************************************
    
    parBehav = 0.0
    
    !TYPE 1. Surface oriented. Particle swims up if deeper than 1 m.
    IF (P_behave(n).EQ.1) THEN 
       btest = 0   !switch to control behavior
       
       !particle has 80% chance of swimming up if deeper than 1.0 m of bottom
       if (P_zc .LT. (P_zetac-1.0)) then
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.80 
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3) 
          btest = 1
       end if
       
       !if within 1 m of surface, swim randomly (50% chance of swimming up)
       if (btest.EQ.0) then    
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.5 
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)  
       end if
       
    END IF
    
    
    !TYPE 2. Near-bottom. Particle swim down if not within 1 m of bottom.
    IF (P_behave(n).EQ.2 .OR. (P_behave(n).EQ.5 .AND. timer(n).GT.0.0)) THEN
       btest = 0   !switch to control behavior
       
       !particle has 80% change of swimming down if greater than 1.0 m 
       !  from bottom
       if (P_zc .GT. (P_depth+1.0)) then
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.20 
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)
          btest = 1
       end if
       
       !if within 1 m of bottom, just swim randomly
       if (btest.EQ.0) then    
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.5 
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)   
       end if
       
    END IF
    
    !TYPE 3: Diurnal Vertical Migration
    IF (P_behave(n).EQ.3) THEN
       
       !A. Find daytime in hrs since midnight (dtime)
       dtime = (daytime - aint(daytime))*DBLE(24.0)  !time of day 
       !This assumes that model simulations start at midnight
       
       !B. Calcluate irradiance at the water's surface (E0)
       tst = 0.0  !seconds since twilight start
       E0 = 0.0   !irradiance at the water's surface
       if (dtime.GT.twiStart .AND. dtime.LT.twiEnd) then
          tst=(dtime-twiStart)*DBLE(3600.0)
          E0= Em*SIN(PI*tst/(daylength*DBLE(3600.0)))*     &
                 SIN(PI*tst/(daylength*DBLE(3600.0)))
       else 
          E0 = 0.0
       end if
       
       !C. Calcluate irradiance at depth of the particle
       P_light = E0 * exp(Kd*P_zc)
       
       !If light at particle location is less than threshold, random swimming
       if (P_light.LT.thresh ) then
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.5 
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)   
       end if
       !If light at particle > threshold, then have 80% chance of swimming down
       if (P_light.GT.thresh ) then
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.20 
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)  
       end if
       
    END IF
    
    
    !TYPE 4. Crassostrea virginica -- above the halocline
    IF (P_behave(n).EQ.4) THEN 
       if (it.EQ.1) then
          P_Sprev(n) = P_S  !for first iteration
          P_zprev(n) = P_zc
       endif
       btest = 0   !switch to control behavior
       Sslope = 0.0  !salinity gradient that larvae swam through
       
       !determine if larva swam through salinity gradient large enough to 
       !  cue behavior;  if so, then 80% chance of swimming up                                            
       deltaS = P_Sprev(n) - P_S
       deltaz = P_zprev(n) - P_zc
       if (it.GT.1) Sslope = deltaS/deltaz
       if (abs(Sslope).GT.Sgradient) then
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.80 
          if (dev1.GT.switch) negpos = -1.0
          parBehav=negpos*P_swim(n,3)
          btest = 1
       endif
       
       !if no directed swimming, swim randomly with probabilities that result
       !  in particles moving up initially, then slowly moving toward bottom 
       !  with increasing age
       if (btest.EQ.0) then    
          negpos = 1.0
          dev1=genrand_real1()
          if (P_age .LT. 1.5*24.*3600.) then     !if Age < 1.5 Days
            switch = 0.1
          elseif (P_age .LT. 5.*24.*3600.) then  !if 1.5 Days <= Age < 5.0 Days
            switch = 0.49
          elseif (P_age .LT. 8.*24.*3600.) then  !if 5.0 Days <= Age < 8.0 Days
            switch = 0.50
          else                                   !if Age >= 8.0 Days
             switchslope = (DBLE(0.50)-DBLE(0.517)) /                          &
                           (DBLE(8.0)*DBLE(24.0)*DBLE(3600.0) - P_pediage(n))
             switch = switchslope*P_age + DBLE(0.50) -                         &
                      switchslope*DBLE(8.0)*DBLE(24.0)*DBLE(3600.0) 
             if (P_zc .LT. P_depth+1.) switch = 0.5
          endif
          if (dev1.GT.(1-switch)) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)  
       endif
       
       !update previous salt and depth matrix for next iteration
       P_Sprev(n) = P_S
       P_zprev(n) = P_zc
    ENDIF
    
    
    !TYPE 5. Crassostrea ariakensis -- below the halocline
    IF (P_behave(n).EQ.5 .AND. timer(n).EQ.0.0) THEN 
       if (it.EQ.1) then
          P_Sprev(n) = P_S  !for first iteration
          P_zprev(n) = P_zc
       endif
       btest = 0   !switch to control behavior
       Sslope = 0.0  !salinity gradient that larvae swam through
       
       !determine if larva swam through salinity gradient large enough to 
       !  cue behavior.  If so, then 80% chance of swimming down. Set timer
       !  to keep particle near bottom for 2 hrs
       deltaS = P_Sprev(n) - P_S
       deltaz = P_zprev(n) - P_zc
       if (it.GT.1) Sslope = deltaS/deltaz
       if (abs(Sslope).GT.Sgradient) then
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.20 
          btest = 1
          timer(n) = DBLE(2.0)*DBLE(3600.0)  !2 hr times 3600 s
          if (dev1.GT.switch) negpos = -1.0
          parBehav=negpos*P_swim(n,3) 
          !keep bottom oriented behavior from starting until after particle
          !  is 3.5 days old  
          if (P_age .LT. 3.5*24.*3600.) then  
             btest = 0
             timer(n) = 0.
          endif
       endif
       
       !if no directed swimming, just swim randomly with probabilities that 
       !  result in particles moving up initially, then moving toward bottom 
       !  with increasing age
       if (btest.EQ.0) then    
          negpos = 1.0
          dev1=genrand_real1()
          switch = 0.495 
          if (P_age .LT. 1.5*24.*3600.) switch = 0.9
          if (P_age .GT. 2.0*24.*3600. .AND. P_age .LT. 3.5*24.*3600.) then
             switchslope = (DBLE(0.3)-DBLE(0.495)) /                 &
                           (DBLE(2.0)*DBLE(24.0)*DBLE(3600.0) -      &
                            DBLE(3.5)*DBLE(24.0)*DBLE(3600.0))
             switch = switchslope*P_age+DBLE(0.3) -                  &
                      switchslope*DBLE(2.0)*DBLE(24.0)*DBLE(3600.0) 
          endif
          if (dev1.GT.switch) negpos = -1.0
          devB=genrand_real1()
          parBehav=negpos*devB*P_swim(n,3)  
       endif
       
       !update previous salt and depth matrix for next iteration        
       P_Sprev(n) = P_S
       P_zprev(n) = P_zc
    ENDIF
    
    !TYPE 6. Constant -- no random motion to vertical movement
    IF ((P_behave(n).EQ.6)) THEN               
       if(P_age .GE. swimstart) then
          parBehav = sink 
       else
          parBehav = P_swim(n,3)
       endif
       
       !Note: the code below is included if someone wants to calculate density
       ! ! To calculate salinity (P_S) and temperature (P_T) at particle location
       !  do i=3,us-2
       !    if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR.     & 
       !        (Zpar .LT. Pwc_zf(i))) exit
       !  enddo
       !  deplvl = i-2   !depth level
       !
       !  !Salinity at particle location
       !  P_S = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,     &
       !                  P_zb,P_zc,P_zf,ex,ix,p,4)
       ! 
       !  !Temperature at particle location 
       !  P_T = WCTS_ITPI("temp",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,     &
       !                  P_zb,P_zc,P_zf,ex,ix,p,4)
    ENDIF
    
    !Calculate movement due to behavior for all behaviors other than 7
    ZBehav = parBehav * idt
    
    !TYPE 7. Tidal Stream Transport: if flooding, then swim in direction of
    !  currents, else sit on bottom
    IF ((P_behave(n).EQ.7)) THEN 
       ! Set initial values for the first iteration
       if (it.EQ.1) then         
          P_Sprev(n) = P_S
          currentspeed = 0.0    
       endif
       !Find current speed at the particle location ( c = sqrt(a**2 + b**2) )
       currentspeed = sqrt( (P_U*cos(P_angle) - P_V*sin(P_angle))**2 +         &
                            (P_U*sin(P_angle) + P_V*cos(P_angle))**2 )
       if (bottom(n) .EQV. .TRUE.) then          !CRS
          !if particle is on bottom, test if salinity is increasing
          if (P_Sprev(n).LT.P_S) then       !if salinity is increasing:
             bottom(n) = .FALSE.            !  come off bottom 
             ZBehav = P_depth + Swimdepth   !  and swim to the swimming depth
          else
             ZBehav = -9999  !if salinity is not increasing, stay on bottom
          end if
       else        
          !if particle is not on bottom, test if currents are not slack 
          !  (defined as 0.05 m/s)
          if (currentspeed.GT.0.05) then
             !if the current speed is greater than 0.05 m/s, then swim in the 
             !  direction of the current
             Hdistance = Hswimspeed*idt
             !find theta of currents
             theta = atan( (P_U*sin(P_angle) + P_V*cos(P_angle)) /             &
                           (P_U*cos(P_angle) - P_V*sin(P_angle)) )
             X = (P_U*cos(P_angle) - P_V*sin(P_angle))
             Y = (P_U*sin(P_angle) + P_V*cos(P_angle))
             if(X.GT.0.0) then
                XBehav =  Hdistance*cos(theta)     
                YBehav =  Hdistance*sin(theta)     
             end if
             if(X.LT.0.0) then
                XBehav =  DBLE(-1.0)*Hdistance*cos(theta)     
                YBehav =  DBLE(-1.0)*Hdistance*sin(theta)     
             end if
             if(X.EQ.0 .AND. Y.GE.0.0) then
                XBehav =  0.0     
                YBehav =  Hdistance
             end if
             if(X.EQ.0 .AND. Y.LE.0.0) then
                XBehav =  0.0     
                YBehav =  DBLE(-1.0)*Hdistance
             end if
             !keep vertical position of particle at swim depth
             ZBehav = P_depth + Swimdepth
          else
             ZBehav = -9999       !if the current speed is less than 0.05 m/s,
             bottom(n) = .TRUE.   !  then swim to bottom
          end if
       end if
       bott = bottom(n)
    ENDIF

! ******************* End Particle Behavior ******************************
  END SUBROUTINE behave


  INTEGER FUNCTION getStatus(n)
    !This function returns an identification number that describes a particle's  
    !behavior type or status for use in visualization routines. It was
    !initially developed to contain the color code for plotting in Surfer.)                
    USE PARAM_MOD, ONLY: SETTLEMENTON,OPENOCEANBOUNDARY
    USE SETTLEMENT_MOD, ONLY: isSettled
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    getStatus = P_behave(n)          ! Set Status to behavior ID
                                     ! Change if Dead, Settled, or OutOfBounds

    if(dead(n)) getStatus = -1         ! -1 = Dead
    if(settlementon)then
      if(isSettled(n)) getStatus = -2  ! -2 = Settled
    endif
    if(OpenOceanBoundary)then
      if(oob(n)) getStatus = -3        ! -3 = Out of Bounds
    endif

  END FUNCTION getStatus


  LOGICAL FUNCTION isDead(n)
  !This function returns .TRUE. if the particle is "dead", and FALSE if not
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    isDead = dead(n)

  END FUNCTION isDead


  SUBROUTINE die(n)
  !This subroutine sets the value of dead(n) to TRUE, indicating 
  !  the particle is "dead"
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    dead(n) = .TRUE.

  END SUBROUTINE die


  SUBROUTINE setOut(n)
    !This subroutine changes particle n's status to out of bounds
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    oob(n) = .TRUE.

  END SUBROUTINE setOut

  LOGICAL FUNCTION isOut(n)
    !This function returns the value of oob for particle n
    ! i.e. Returns True if the particle is out of bounds, False if in bounds
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    isOut = oob(n)

  END FUNCTION isOut


  SUBROUTINE finBehave()    !Finish the behavior module
  USE PARAM_MOD     , ONLY: settlementon
  USE SETTLEMENT_MOD, ONLY: finSettlement
  IMPLICIT NONE

    !Deallocate Behavior Variables
    DEALLOCATE(timer)
    DEALLOCATE(P_behave)
    DEALLOCATE(P_pediage)
    DEALLOCATE(P_deadage)
    DEALLOCATE(P_Sprev)
    DEALLOCATE(P_zprev)
    DEALLOCATE(P_swim)
    DEALLOCATE(dead)
    DEALLOCATE(oob)

    if(ALLOCATED(bottom))DEALLOCATE(bottom)

    !If Settlement is on, Deallocate Settlement Variables
    if(settlementon) call finSettlement()

  END SUBROUTINE finBehave


END MODULE BEHAVIOR_MOD

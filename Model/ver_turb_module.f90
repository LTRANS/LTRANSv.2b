MODULE VTURB_MOD

! A random displacement model is implemented to simulate sub-grid scale
!   turbulent particle motion in the vertical (z) direction. 
!
! Created by:             Elizabeth North
! Modified by:            Zachary Schlag
! Created on:             2003
! Last Modified on:       18 Aug 2008

IMPLICIT NONE
PUBLIC

CONTAINS

!      ********************** Vertical Turbulence *************************
!     ***********************************************************************
!    ** Random Displacement Model (Visser 1997 MEPS 158:275-281) for        **
!    ** simulating displacement due to turbulent diffusion                  **
!    ** (vertical direction)                                                **
!    **    z(t+1)= z + K'(z)dt + R{ 2/r K[z+0.5K'(z)dt]dt}**0.5             **
!    **    where z = particle vertical location at time t                   **
!    **      K' = dK/dz (Kprime) and K = vertical diffusivity (KH from ROMS)**
!    **      dt = time step of RDM (deltat)                                 **
!    **      R  = random process with mean = 0 and standard deviation = r.  **
!    **                                                                     **
!    ** Programmed by EW North February 2003 UMCES HPL enorth@hpl.umces.edu **
!    *************************************************************************

  SUBROUTINE VTurb(P_zc,P_depth,P_zetac,p,ex,ix,Pwc_wzb,Pwc_wzc,Pwc_wzf,TurbV)
    USE PARAM_MOD,  ONLY: ws,idt
    USE HYDRO_MOD,  ONLY: getInterp
    USE INT_MOD,    ONLY: linint,polintd
    USE NORM_MOD,   ONLY: norm
    USE TENSION_MOD, ONLY: TSPSI,HVAL,HPVAL
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p
    DOUBLE PRECISION, INTENT(IN) :: P_zc,P_depth,P_zetac,ex(3),ix(3),     &
                                    Pwc_wzb(:),Pwc_wzc(:),Pwc_wzf(:)
    DOUBLE PRECISION, INTENT(OUT) :: TurbV

    !Background vertical diffusivity from ROMS
    DOUBLE PRECISION, PARAMETER:: background = 1.0E-6

    DOUBLE PRECISION :: deltat,DEV,r
    INTEGER :: i,j,k,jlo,loop
    DOUBLE PRECISION :: slopem,ParZc,Kprimec,KprimeZc,newZc,KH3rdc,Z3rdc, &
                        thisyc,ey(3)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: slopekb,slopekc,slopekf,    &
      interceptb,interceptc,interceptf,movexb,moveyb,movexc,moveyc,            &
      movexf,moveyf,ifitxb,ifityb,ifitxc,ifityc,ifitxf,ifityf,ifitx,ifity,     &
      Pwc_KHb,Pwc_KHc,Pwc_KHf,newxb,newyb,newxc,newyc,newxf,newyf

    !TSPACK Variables
    INTEGER :: IER,SigErr
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: YPKc,SIGMAKc

    !Number of values to proliferate to
    INTEGER :: p2

    p2 = ws*4

    !ALLOCATE VARIABLES
    ALLOCATE(Pwc_KHb(ws))
    ALLOCATE(Pwc_KHc(ws))
    ALLOCATE(Pwc_KHf(ws))
    ALLOCATE(slopekb(ws-1))
    ALLOCATE(slopekc(ws-1))
    ALLOCATE(slopekf(ws-1))
    ALLOCATE(interceptb(ws-1))
    ALLOCATE(interceptc(ws-1))
    ALLOCATE(interceptf(ws-1))
    ALLOCATE(movexb(p2+1))
    ALLOCATE(moveyb(p2+1))
    ALLOCATE(movexc(p2+1))
    ALLOCATE(moveyc(p2+1))
    ALLOCATE(movexf(p2+1))
    ALLOCATE(moveyf(p2+1))
    ALLOCATE(ifitxb(p2))
    ALLOCATE(ifityb(p2))
    ALLOCATE(ifitxc(p2))
    ALLOCATE(ifityc(p2))
    ALLOCATE(ifitxf(p2))
    ALLOCATE(ifityf(p2))
    ALLOCATE(ifitx(p2))
    ALLOCATE(ifity(p2))
    ALLOCATE(newxb(p2+7))
    ALLOCATE(newyb(p2+7))
    ALLOCATE(newxc(p2+7))
    ALLOCATE(newyc(p2+7))
    ALLOCATE(newxf(p2+7))
    ALLOCATE(newyf(p2+7))
    ALLOCATE(YPKc(p2))
    ALLOCATE(SIGMAKc(p2))
 
    !       *********************************************************
    !       *           Find kh in water column profile             *
    !       *********************************************************

    !i. find KH in water column profile at particle location
    do i=1,ws

        Pwc_KHb(i) = getInterp("khb",i)
        Pwc_KHc(i) = getInterp("khc",i)
        Pwc_KHf(i) = getInterp("khf",i)

    enddo

    !       *********************************************************
    !       *           Proliferate depth and kh values             *
    !       *********************************************************

    !ii. proliferate data points

    !  a. proliferate x(z) and y(KH) values
    newxb=0.0
    newxc=0.0
    newxf=0.0
    do j=1,p2+7
      newxb(j)=Pwc_wzb(1)+(float(j-4))*(Pwc_wzb(ws)-Pwc_wzb(1))/DBLE(p2)
      newxc(j)=Pwc_wzc(1)+(float(j-4))*(Pwc_wzc(ws)-Pwc_wzc(1))/DBLE(p2)
      newxf(j)=Pwc_wzf(1)+(float(j-4))*(Pwc_wzf(ws)-Pwc_wzf(1))/DBLE(p2)
    enddo

    do i=1,ws-1
      slopekb(i)=(Pwc_KHb(i)-Pwc_KHb(i+1))/(Pwc_wzb(i)-Pwc_wzb(i+1))
      interceptb(i) = Pwc_KHb(i) - slopekb(i)*Pwc_wzb(i)
      slopekc(i)=(Pwc_KHc(i)-Pwc_KHc(i+1))/(Pwc_wzc(i)-Pwc_wzc(i+1))
      interceptc(i) = Pwc_KHc(i) - slopekc(i)*Pwc_wzc(i)
      slopekf(i)=(Pwc_KHf(i)-Pwc_KHf(i+1))/(Pwc_wzf(i)-Pwc_wzf(i+1))
      interceptf(i) = Pwc_KHf(i) - slopekf(i)*Pwc_wzf(i)
    enddo

    newyb=0.0
    jlo=1
    do j=5,p2+3
      do
        if(Pwc_wzb(jlo+1).gt.newxb(j))exit
        jlo=jlo+1
      enddo
      newyb(j)=slopekb(jlo)*newxb(j) + interceptb(jlo)
    enddo

    newyc=0.0
    jlo=1
    do j=5,p2+3
      do
        if(Pwc_wzc(jlo+1).gt.newxc(j))exit
        jlo=jlo+1
      enddo
      newyc(j)=slopekc(jlo)*newxc(j) + interceptc(jlo)
    enddo

    newyf=0.0
    jlo=1
    do j=5,p2+3      
      do
        if(Pwc_wzf(jlo+1).gt.newxf(j))exit
        jlo=jlo+1
      enddo
      newyf(j)=slopekf(jlo)*newxf(j) + interceptf(jlo)
    enddo

    !       *********************************************************
    !       *          Create Extra Points for Moving Average       *
    !       *********************************************************

    !  b. create extra end point for moving average
    do i=1,4
      newyb(i)=Pwc_KHb(1) 
      newyc(i)=Pwc_KHb(1) 
      newyf(i)=Pwc_KHb(1) 
      newyb(i+p2+3)=Pwc_KHb(ws)
      newyc(i+p2+3)=Pwc_KHc(ws)
      newyf(i+p2+3)=Pwc_KHf(ws)
    enddo

    !           *********************************************************
    !           *               Take Moving Average                     *
    !           *********************************************************

    !iii. take 8-point moving average of proliferated values
    do i=2,p2-1
      moveyb(i)=( newyb(i)+newyb(i+1)+newyb(i+2)+newyb(i+3)+  &
                  newyb(i+4)+newyb(i+5)+newyb(i+6)+newyb(i+7) ) / DBLE(8.0)
      movexb(i)=newxb(i)+ (newxb(i+7)-newxb(i))/DBLE(2.0)
      moveyc(i)=( newyc(i)+newyc(i+1)+newyc(i+2)+newyc(i+3)+  &
                  newyc(i+4)+newyc(i+5)+newyc(i+6)+newyc(i+7) ) / DBLE(8.0)
      movexc(i)=newxc(i)+ (newxc(i+7)-newxc(i))/DBLE(2.0)
      moveyf(i)=( newyf(i)+newyf(i+1)+newyf(i+2)+newyf(i+3)+  &
                  newyf(i+4)+newyf(i+5)+newyf(i+6)+newyf(i+7) ) / DBLE(8.0)
      movexf(i)=newxf(i)+ (newxf(i+7)-newxf(i))/DBLE(2.0)
    enddo

    !iv. fix end points of series to true data values
    movexb(1)= Pwc_wzb(1)
    moveyb(1)= Pwc_KHb(1) 
    movexb(p2)=Pwc_wzb(ws)
    moveyb(p2)=Pwc_KHb(ws)

    movexc(1)= Pwc_wzc(1)
    moveyc(1)= Pwc_KHc(1) 
    movexc(p2)=Pwc_wzc(ws)
    moveyc(p2)=Pwc_KHc(ws)

    movexf(1)= Pwc_wzf(1)
    moveyf(1)= Pwc_KHf(1) 
    movexf(p2)=Pwc_wzf(ws)
    moveyf(p2)=Pwc_KHf(ws)

    !       *********************************************************
    !       *      Interpolate from external b,c,f to internal      *
    !       *********************************************************

    !v. fit polynomical and find internal time step b,c,f values at each depth

    !  a. Z profiles
    !    1. Prepare external time step values
    do k=1,p2
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = movexb(k)
        ey(2) = movexb(k)
        ey(3) = movexc(k)
      else
        ey=0.0
        ey(1) = movexb(k)
        ey(2) = movexc(k)
        ey(3) = movexf(k)
      endif
    !    2. Call polynomial interpolator
      ifitxb(k) = 0.0
      ifitxc(k) = 0.0
      ifitxf(k) = 0.0
      ifitxb(k) = polintd(ex,ey,3,ix(1))
      ifitxc(k) = polintd(ex,ey,3,ix(2))
      ifitxf(k) = polintd(ex,ey,3,ix(3))
    enddo

    !  b. KH profiles
    !    1. Prepare external time step values
    do k=1,p2
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = moveyb(k)
        ey(2) = moveyb(k)
        ey(3) = moveyc(k)
      else
        ey=0.0
        ey(1) = moveyb(k)
        ey(2) = moveyc(k)
        ey(3) = moveyf(k)
      endif
    !    2. Call Polynomial interpolator
      ifityb(k) = 0.0
      ifityc(k) = 0.0
      ifityf(k) = 0.0
      ifityb(k) = polintd(ex,ey,3,ix(1))
      ifityc(k) = polintd(ex,ey,3,ix(2))
      ifityf(k) = polintd(ex,ey,3,ix(3))
    enddo

    do k=1,p2
      if (ifityb(k).LT.0.0) ifityb(k)=0.0
      if (ifityc(k).LT.0.0) ifityc(k)=0.0
      if (ifityf(k).LT.0.0) ifityf(k)=0.0
    enddo


    !vi. Calculate profile of vertical diffusivity using R-K 4th order in time.
    do k=1,p2
      ifity(k)=( ifityb(k) + DBLE(4.0)*ifityc(k) + ifityf(k) )/DBLE(6.0)
      ifitx(k)=( ifitxb(k) + DBLE(4.0)*ifitxc(k) + ifitxf(k) )/DBLE(6.0)
    enddo

    !vii. fit a tension spline to water column profile of KH using TSPACK
    SigErr=0
    CALL TSPSI (p2,ifitx,ifity,YPKc,SIGMAKc,IER,SigErr)

    !viii. Initialize. Set deltat, number of iterations, initial z-coordinate
    deltat=2.0  ! dt= 1 sec
    loop= idt/int(deltat) ! number of iterations of RDM loop  
    ParZc = P_zc  !set initial particle depth  

    !       *********************************************************
    !       *           Random Displacement Model Loop              *
    !       *********************************************************

    !ix. Begin iterations
    do i=1,loop   
    !  a. Determine the second term in the equation: K'(z)deltat
    !    1. Find Kprime & solve for second term in RDM equation
    !       (K'(z)deltat = KprimeZ)
      Kprimec=0.0
      KprimeZc=0.0
      if (ParZc.LT.P_depth .OR. ParZc.GT.P_zetac) then
        Kprimec=0.0                                                
      else
        if (SigErr.EQ.0) then
          Kprimec= HPVAL (ParZc,p2,ifitx,ifity,YPKc,SIGMAKc,IER)
        else
          CALL linint(ifitx,ifity,p2,ParZc,thisyc,Kprimec)
        endif
      endif
      KprimeZc=DBLE(-1.0)*Kprimec*deltat    !need minus sign to switch axes   

    !  b. Determine the 3rd term in the RDM equation: 
    !     R{ 2/r K[z+0.5K'(z)dt]dt}**0.5
    !    i. Find K at location of [z+0.5K'(z)dt] = Z3rd
    !      1. calculate Z3rd and make sure within boudaries
      Z3rdc = ParZc + DBLE(0.5)*KprimeZc

    !      2. Find KH at the location Z3rd
      KH3rdc=0.0
      if (Z3rdc.LT.P_depth .OR. Z3rdc.GT.P_zetac) then
        KH3rdc=background                                
      else
        if (SigErr.EQ.0) then
          KH3rdc= HVAL (Z3rdc,p2,ifitx,ifity,YPKc,SIGMAKc,IER)
        else
          CALL linint(ifitx,ifity,p2,Z3rdc,KH3rdc,slopem)
        endif
        if (KH3rdc.LT.background) KH3rdc=background    
      endif

    !  c. Solve the entire equation
      DEV=norm()        ! the random deviate
      r=1.                ! the standard deviation of the random deviate
      newZc = ParZc + KprimeZc + DEV* (DBLE(2.0)/r * KH3rdc*deltat)**0.5

    !x. update particle z-coordinate
      ParZc = newZc

    enddo   !End RDM iterations

!           *********************************************************
!           *       Calculate displacement due to vertical turb     *
!           *********************************************************

    !xi. find vertical displacement of particle due to turbulent diffusion
    TurbV = P_zc - ParZc

! **************** End Vertical Turbulence (RDM) *************************
! ************************************************************************

    !DEALLOCATE VARIABLES
    DEALLOCATE(interceptb)
    DEALLOCATE(interceptc)
    DEALLOCATE(interceptf)
    DEALLOCATE(Pwc_KHb)
    DEALLOCATE(Pwc_KHc)
    DEALLOCATE(Pwc_KHf)
    DEALLOCATE(slopekb)
    DEALLOCATE(slopekc)
    DEALLOCATE(slopekf)
    DEALLOCATE(SIGMAKc)
    DEALLOCATE(movexb)
    DEALLOCATE(moveyb)
    DEALLOCATE(movexc)
    DEALLOCATE(moveyc)
    DEALLOCATE(movexf)
    DEALLOCATE(moveyf)
    DEALLOCATE(ifitxb)
    DEALLOCATE(ifityb)
    DEALLOCATE(ifitxc)
    DEALLOCATE(ifityc)
    DEALLOCATE(ifitxf)
    DEALLOCATE(ifityf)
    DEALLOCATE(ifitx)
    DEALLOCATE(ifity)
    DEALLOCATE(newxb)
    DEALLOCATE(newyb)
    DEALLOCATE(newxc)
    DEALLOCATE(newyc)
    DEALLOCATE(newxf)
    DEALLOCATE(newyf)
    DEALLOCATE(YPKc)

  END SUBROUTINE VTurb

END MODULE VTURB_MOD
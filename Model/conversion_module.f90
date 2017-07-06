MODULE CONVERT_MOD

! The Conversion Module contains all the procedures necessary to convert 
!  locations between latitude and longitude coordinates and metric (x and y) 
!  coordinates.  The default setting assumes a spherical earth 
!  (SphericalProjection = .TRUE. in LTRANS.data). This code is based on the 
!  equations in earthdist.m from Seagrid which was written by Rich Signel. 
!  These equations and the inverse of earthdist.m were coded in Fortran by 
!  Elizabeth North. 
!  
!  If the user specifies (SphericalProjection = .FALSE. in LTRANS.data), then
!  a Mercator projection is used. The equations are from the sg_mercator.m and
!  seagrid2roms.m matlab scripts that are found in Seagrid.  
!
!  This Conversion Module module contains four interface blocks to cover the 
!  four necessary conversions:  
!
!  lon2x = convert longitude to x-coordinate 
!  lat2y = convert latitude  to y-coordinate 
!  x2lon = convert x-coordinate to longitude 
!  y2lat = convert y-coordinate to latitude.
!
!  NOTE: The Interface blocks were created so that the functions can handle 
!        both REAL and DOUBLE PRECISION variables as input 
!  NOTE: Output is DOUBLE PRECISION regardless of input
!  NOTE: RCF = radian conversion factor (RCF = 180/pi)
!  NOTE: The values of lonmin and latmin which are specified in LTRANS.data by
!        the user are altered in the PARAMETER MODULE. 1 is subtracted from
!        each. This is done to ensure that there are no problems that could 
!        result from rounding by the user. 
!
!  The functions used to make these conversions are as follows:

!  Spherical Projection:
!   x&y = acos(ax*bx + ay*by + az*bz) * Earth_Radius
!   Note: both x & y are calculated as the arclength of the central angle
!     between the two vectors [ax,ay,az] and [bx,by,bz] multiplied by the
!     Earth's radius
!   lon = lonmin + RCF*acos( (cos(x/Earth_Radius) - sin(lat)*sin(lat)) / 
!                                                  (cos(lat)*cos(lat)) )
!   lat = y * RCF / Earth_Radius + latmin

!  Mercator Projection:
!   x   = lon / RCF * Earth_Radius
!   y   = log( tan( pi/4.0 + lat/(RCF*2.0) ))*Earth_Radius
!   lon = x / Earth_Radius * RCF
!   lat = 2.0*RCF * ( atan(exp(y/Earth_Radius)) - pi/4.0 )
!
! Created by:           Zachary Schlag
! Created on:           23 Jul 2008

! Last modified on:     08 Jul 2015
! Last modified by:     K.A.S. Mislan
! Last modification:    Added code to correctly convert grids that span >180 deg


USE PARAM_MOD, ONLY: PI,Earth_Radius,SphericalProjection,lonmin,latmin
IMPLICIT NONE
PRIVATE
SAVE


  !Return x location given longitude
  INTERFACE lon2x
    MODULE PROCEDURE rlon2x  !real input
    MODULE PROCEDURE dlon2x  !double precision input
  END INTERFACE lon2x


  !Return y location given latitude
  INTERFACE lat2y
    MODULE PROCEDURE rlat2y  !real input
    MODULE PROCEDURE dlat2y  !double precision input
  END INTERFACE lat2y


  !Return longitude given x location
  INTERFACE x2lon
    MODULE PROCEDURE rx2lon  !real input
    MODULE PROCEDURE dx2lon  !double precision input
  END INTERFACE x2lon


  !Return latitude given y location
  INTERFACE y2lat
    MODULE PROCEDURE ry2lat  !real input
    MODULE PROCEDURE dy2lat  !double precision input
  END INTERFACE y2lat

  !The following procedures have been made public for the use of other modules:
  PUBLIC :: lon2x,lat2y,x2lon,y2lat

CONTAINS

  DOUBLE PRECISION FUNCTION rlon2x(lon,lat)
    IMPLICIT NONE
    REAL, INTENT(IN) :: lon
    REAL, INTENT(IN), OPTIONAL :: lat
    DOUBLE PRECISION :: RCF  ! Radian conversion factor
    DOUBLE PRECISION :: c,ax,ay,az,bx,by,bz,alon,alat,blon,blat,x180,xtralon

    RCF = DBLE(180.0) / PI

    if(SphericalProjection)then
      if( present(lat) )then

        if((lon-lonmin) .gt. 180.0)then
          alon = 180.0
          alat = lat
          blon = 0
          blat = lat

          alon = alon / RCF
          alat = alat / RCF
          blon = blon / RCF
          blat = blat / RCF

          c  = cos(alat)
          ax = c * cos(alon)
          ay = c * sin(alon)
          az = sin(alat)

          c  = cos(blat)
          bx = c * cos(blon)
          by = c * sin(blon)
          bz = sin(blat)

          x180 = acos(ax*bx + ay*by + az*bz) * Earth_Radius

          alon = lon
          blon = lonmin+180.0

          alon = alon / RCF
          blon = blon / RCF

          ax = c * cos(alon)
          ay = c * sin(alon)
          bx = c * cos(blon)
          by = c * sin(blon)

          xtralon = acos(ax*bx + ay*by + az*bz) * Earth_Radius
          rlon2x = x180 + xtralon
        else
          alon = lon
          alat = lat
          blon = lonmin
          blat = lat

          alon = alon / RCF
          alat = alat / RCF
          blon = blon / RCF
          blat = blat / RCF

          c  = cos(alat)
          ax = c * cos(alon)
          ay = c * sin(alon)
          az = sin(alat)

          c  = cos(blat)
          bx = c * cos(blon)
          by = c * sin(blon)
          bz = sin(blat)

          rlon2x = acos(ax*bx + ay*by + az*bz) * Earth_Radius
        endif
      else
        write(*,*) "Problem lon2x: spherical projection without lat value"
      endif
    else
      rlon2x = lon / RCF * Earth_Radius
    endif

  END FUNCTION rlon2x


  DOUBLE PRECISION FUNCTION dlon2x(lon,lat)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: lon
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: lat
    DOUBLE PRECISION :: RCF  ! Radian conversion factor
    DOUBLE PRECISION :: c,ax,ay,az,bx,by,bz,alon,alat,blon,blat,x180,xtralon

    RCF = DBLE(180.0) / PI

    if(SphericalProjection)then
      if( present(lat) )then

	dlon2x = (lon-lonmin)/180.0 * Earth_Radius * pi * cos(lat/RCF)

      else
        write(*,*) "Problem lon2x: spherical projection without lat value"
      endif
    else
      dlon2x = lon / RCF * Earth_Radius
    endif

  END FUNCTION dlon2x


  DOUBLE PRECISION FUNCTION rlat2y(lat,lon)
    IMPLICIT NONE
    REAL, INTENT(IN) :: lat
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: lon
    DOUBLE PRECISION :: RCF  ! Radian conversion factor
    DOUBLE PRECISION :: c,ax,ay,az,bx,by,bz,alon,alat,blon,blat

    RCF = DBLE(180.0) / PI

    if(SphericalProjection) then
      if( present(lon) )then
        alon = lon
        alat = lat
        blon = lon
        blat = latmin
      else
	alon = lonmin
	alat = lat
        blon = lonmin
        blat = latmin   
      endif

      alon = alon / RCF
      alat = alat / RCF
      blon = blon / RCF
      blat = blat / RCF

      c = cos(alat)
      ax = c * cos(alon)
      ay = c * sin(alon)
      az = sin(alat)

      c = cos(blat)
      bx = c * cos(blon)
      by = c * sin(blon)
      bz = sin(blat)

      rlat2y = acos(ax*bx + ay*by + az*bz) * Earth_Radius

    else

      rlat2y = log( tan( pi/4.0 + lat/(RCF*2.0) ))*Earth_Radius

    endif

  END FUNCTION rlat2y


  DOUBLE PRECISION FUNCTION dlat2y(lat)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: lat
    DOUBLE PRECISION :: RCF  ! Radian conversion factor
    DOUBLE PRECISION :: c,ax,ay,az,bx,by,bz,alon,alat,blon,blat

    RCF = DBLE(180.0) / PI

    if(SphericalProjection) then

      dlat2y = (lat-latmin)*Earth_Radius*pi/180.0

    else

      dlat2y = log( tan( pi/4.0 + lat/(RCF*2.0) ))*Earth_Radius

    endif

  END FUNCTION dlat2y


  DOUBLE PRECISION FUNCTION rx2lon(x,y)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, INTENT(IN), OPTIONAL :: y
    DOUBLE PRECISION :: RCF  ! Radian conversion factor
    DOUBLE PRECISION :: lat,c,ax,ay,az,bx,by,bz,alon,alat,blon,blat,x180

    RCF = DBLE(180.0) / PI

    if(SphericalProjection) then

      if( PRESENT(y) )then
        lat = (y * RCF / Earth_Radius + latmin) / RCF

        alon = 0
        alat = lat
        blon = 180.0
        blat = lat

        alon = alon / RCF
        blon = blon / RCF

        c  = cos(alat)
        ax = c * cos(alon)
        ay = c * sin(alon)
        az = sin(alat)

        c  = cos(blat)
        bx = c * cos(blon)
        by = c * sin(blon)
        bz = sin(blat)

        x180 = acos(ax*bx + ay*by + az*bz) * Earth_Radius

        if(x .gt. x180)then

          rx2lon = lonmin + 180.0 + RCF*acos( (cos((x180-x)/Earth_Radius) - sin(lat)*sin(lat)) &
                                                          / (cos(lat)*cos(lat)) )
        else
          rx2lon = lonmin + RCF*acos( (cos(x/Earth_Radius) - sin(lat)*sin(lat)) &
                                                        / (cos(lat)*cos(lat)) )
        endif
      else
        write(*,*) "Problem x2lon: Spherical projection without y value"
      endif

    else
      rx2lon = x / Earth_Radius * RCF
    endif

  END FUNCTION rx2lon


  DOUBLE PRECISION FUNCTION dx2lon(x,y)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: x
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: y
    DOUBLE PRECISION :: RCF  ! Radian conversion factor
    DOUBLE PRECISION :: lat,c,ax,ay,az,bx,by,bz,alon,alat,blon,blat,x180

    RCF = DBLE(180.0) / PI

    if(SphericalProjection) then

      if( PRESENT(y) )then
        lat = y*180.0/(Earth_Radius*pi) + latmin

	dx2lon = x*180.0 / (Earth_Radius * pi * cos(lat/RCF)) + lonmin

      else
        write(*,*) "Problem x2lon: Spherical projection without y value"
      endif

    else
      dx2lon = x / Earth_Radius * RCF
    endif

  END FUNCTION dx2lon


  DOUBLE PRECISION FUNCTION ry2lat(y)
    IMPLICIT NONE
    REAL, INTENT(IN) :: y
    DOUBLE PRECISION :: RCF  ! Radian conversion factor

    RCF = DBLE(180.0) / PI

    if(SphericalProjection) then
      ry2lat = y * RCF / Earth_Radius + latmin
    else
      ry2lat = 2.0*RCF * ( atan(exp(y/Earth_Radius)) - pi/4.0 )
    endif

  END FUNCTION ry2lat


  DOUBLE PRECISION FUNCTION dy2lat(y)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: y
    DOUBLE PRECISION :: RCF  ! Radian conversion factor

    RCF = DBLE(180.0) / PI

    if(SphericalProjection) then
	dy2lat = y * RCF / Earth_Radius + latmin
    else
      dy2lat = 2.0*RCF * ( atan(exp(y/Earth_Radius)) - pi/4.0 )
    endif

  END FUNCTION dy2lat

END MODULE CONVERT_MOD

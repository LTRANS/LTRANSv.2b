MODULE PARAM_MOD 
 
!  The Parameter Module reads in the include file, LTRANS.data, making the 
!  parameters declared within available to all the other modules. It also 
!  reads in information from the NetCDF grid file and calculates values of
!  grid specific parameters, making them available to all the other modules.
! 
!  Created by:            Zachary Schlag 
!  Created on:            28 Jul 2008 
!  Last Modified on:         Feb 2013
 
IMPLICIT NONE 
PUBLIC 
SAVE 

  include 'LTRANS.h'

CONTAINS


  SUBROUTINE getParams()
  !Subroutine to read all input parameters from LTRANS.data 

    character(len=120) :: header
    integer :: istat,err

    err = 0

    OPEN(1,file='LTRANS.data')                  !--- read control variables:
      IF(err == 0) THEN
        READ(1,nml=numparticles ,IOSTAT=istat)  !--- number of particles
        IF(istat/=0)err = 10
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=timeparam    ,IOSTAT=istat)  !--- time info
        IF(istat/=0)err = 20
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=hydroparam   ,IOSTAT=istat)  !--- hydrodynamics info
        IF(istat/=0)err = 30
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=turbparam    ,IOSTAT=istat)  !--- turbulence info
        IF(istat/=0)err = 40
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=behavparam   ,IOSTAT=istat)  !--- behavior info
        IF(istat/=0)err = 50
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=dvmparam     ,IOSTAT=istat)  !--- diurnal vertical migration 
        IF(istat/=0)err = 60
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=settleparam  ,IOSTAT=istat)  !--- settlement info
        IF(istat/=0)err = 70
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=convparam    ,IOSTAT=istat)  !--- unit conversion
        IF(istat/=0)err = 80
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=romsgrid     ,IOSTAT=istat)  !--- roms grid
        IF(istat/=0)err = 90
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=romsoutput   ,IOSTAT=istat)  !--- roms history output file
        IF(istat/=0)err = 100
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=parloc       ,IOSTAT=istat)  !--- particle locations
        IF(istat/=0)err = 110
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=HabPolyLoc   ,IOSTAT=istat)  !--- habitat polygon info
        IF(istat/=0)err = 120
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=output       ,IOSTAT=istat)  !--- output related info
        IF(istat/=0)err = 130
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=other        ,IOSTAT=istat)  !--- other misc 
        IF(istat/=0)err = 140
      ENDIF
    CLOSE(1)

    IF(err == 0) THEN
      call gridData(IOSTAT=istat)
      IF(istat/=0)err = 150
    ENDIF

    SELECT CASE(err)
      CASE(0)
        header='No Errors'
      CASE(10)
        header='Error when reading numparticles, pls check LTRANS.data'
      CASE(20)
        header='Error when reading timeparam, pls check LTRANS.data'
      CASE(30)
        header='Error when reading hydroparam, pls check LTRANS.data'
      CASE(40)
        header='Error when reading turbparam, pls check LTRANS.data'
      CASE(50)
        header='Error when reading behavparam, pls check LTRANS.data'
      CASE(60)
        header='Error when reading behavdvm, pls check LTRANS.data'
      CASE(70)
        header='Error when reading settleparam, pls check LTRANS.data'
      CASE(80)
        header='Error when reading convparam, pls check LTRANS.data'
      CASE(90)
        header='Error when reading romsgrid, pls check LTRANS.data'
      CASE(100)
        header='Error when reading romsoutput, pls check LTRANS.data'
      CASE(110)
        header='Error when reading parloc, pls check LTRANS.data'
      CASE(120)
        header='Error when reading HabPolLoc, pls check LTRANS.data'
      CASE(130)
        header='Error when reading output, pls check LTRANS.data'
      CASE(140)
        header='Error when reading other, pls check LTRANS.data'
      CASE(150)
        header='Error when reading gridinfo, pls check GRID.data'
      CASE DEFAULT
        header='Error: unexpected err number'
    END SELECT

    IF(err/=0) CALL errorHandler(Header,-1)  !print the error message and stop

    !Subract 1 from lonmin and latmin to eliminate roundoff error
    lonmin = lonmin - 1
    latmin = latmin - 1

  END SUBROUTINE getParams


  SUBROUTINE errorHandler(header, flag)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: header
    INTEGER, INTENT(IN) :: flag
   
    IF (flag .eq. -1) THEN
      WRITE(*,"(A120)")header               !print error message in report.txt
      STOP
    ELSE
      WRITE(*,"('***** WARNING *****')")    !print warning message to screen
      WRITE(*,"(A120)")header
    ENDIF
   
  END SUBROUTINE errorHandler


  SUBROUTINE gridData(IOSTAT)
    ! This subroutine was originally a separate program: Grid_Generator.f90
    ! Created by:           Zachary Schlag
    ! Program Created:      28 Aug 2008
    ! Subroutine Created:   12 Nov 2010
    ! Last Modified on:        Feb 2011
    USE netcdf
    IMPLICIT NONE

    INTEGER, INTENT(OUT), OPTIONAL :: IOSTAT

    INCLUDE 'netcdf.inc'

    !NetCDF Variables
    INTEGER :: STATUS,GF_ID,VID,dimid,dimcount
    INTEGER :: xi_rho,xi_u,xi_v,eta_rho,eta_u,eta_v
    REAL, ALLOCATABLE, DIMENSION(:,:) :: mask_rho,mask_u,mask_v

    !Grid File Output Variables
    INTEGER :: nR,nU,nV,maxR,maxU,maxV,wetR,wetU,wetV

    !Iteration Variables
    INTEGER :: i,j,err

    err = 0

    ! *********************** GET GRID INFO ***********************

    ! OPEN NETCDF FILE - GET GF_ID VALUE

    STATUS = NF90_OPEN(NCgridfile,NF90_NOWRITE,GF_ID)
    if (STATUS .NE. NF90_NOERR) then
      write(*,*) 'Problem NF90_OPEN'
      err = 10
    endif

    ! GET VALUES FOR xi_rho,xi_u,xi_v,eta_rho,eta_u,eta_v

      STATUS = NF90_INQ_DIMID(GF_ID,'xi_rho',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_rho'
        err = 20 
      endif
      xi_rho = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'eta_rho',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid eta_rho'
        err = 20 
      endif
      eta_rho = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'xi_u',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_u'
        err = 20 
      endif
      xi_u = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'eta_u',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid eta_u'
        err = 20 
      endif
      eta_u = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'xi_v',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_v'
        err = 20 
      endif
      xi_v = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'eta_v',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid eta_v'
        err = 20 
      endif
      eta_v = dimcount

    ! ALLOCATE VARIABLE ARRAY DIMENSIONS

      ALLOCATE (mask_rho(xi_rho,eta_rho),STAT=STATUS)
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_rho'
        err = 30 
      endif
      ALLOCATE (mask_u(xi_u,eta_u),STAT=STATUS)
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_u'
        err = 30 
      endif
      ALLOCATE (mask_v(xi_v,eta_v),STAT=STATUS)
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_v'
        err = 30 
      endif

    ! READ IN DATA FROM NETCDF FILE TO VARIABLES

      ! rho grid mask
      STATUS = NF90_INQ_VARID(GF_ID,'mask_rho',VID)
      STATUS = NF90_GET_VAR(GF_ID,VID,mask_rho)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read mask_rho'
        err = 40 
      endif

      ! u grid mask
      STATUS = NF90_INQ_VARID(GF_ID,'mask_u',VID)
      STATUS = NF90_GET_VAR(GF_ID,VID,mask_u)
      if(STATUS .NE. NF90_NOERR) then
        write(*,*)'Problem read mask_u'
        err = 40 
      endif

      ! v grid mask
      STATUS = NF90_INQ_VARID(GF_ID,'mask_v',VID)
      STATUS = NF90_GET_VAR(GF_ID,VID,mask_v)
      if(STATUS .NE. NF90_NOERR) then
        write(*,*)'Problem read mask_v'
        err = 40 
      endif

    STATUS = NF90_CLOSE(GF_ID)
    if(STATUS /= NF90_NOERR) then
      write(*,*)'Problem closing GF_ID'
      err = 50
    endif




  ! ********************** MAKE GRID FILE **********************

    !Calculate number of nodes in each grid
    nR = xi_rho * eta_rho
    nU = xi_u * eta_u
    nV = xi_v * eta_v

    !Calculate Maximum number of elements in each grid
    maxR = (xi_rho-1)*(eta_rho-1)
    maxU = (xi_u-1)*(eta_u-1)
    maxV = (xi_v-1)*(eta_v-1)

    !Initialize number of wet elements to 0 for each grid
    wetR = 0
    wetU = 0
    wetV = 0

    !iterate through the elements of all three grids and sum up
    !  the number of elements with at least 1 wet node
    do i=1,xi_rho-1
      do j=1,eta_rho-1
        if(mask_rho(i+1,j+1)==1 .OR. mask_rho(i,j+1)==1 .OR.         &
           mask_rho(i+1,j)==1 .OR. mask_rho(i,j)==1) wetR = wetR + 1
      enddo
    enddo

    do i=1,xi_u-1
      do j=1,eta_u-1
        if(mask_u(i+1,j+1)==1 .OR. mask_u(i,j+1)==1 .OR.             &
           mask_u(i+1,j)==1 .OR. mask_u(i,j)==1) wetU = wetU + 1
      enddo
    enddo

    do i=1,xi_v-1
      do j=1,eta_v-1
        if(mask_v(i+1,j+1)==1 .OR. mask_v(i,j+1)==1 .OR.             &
           mask_v(i+1,j)==1 .OR. mask_v(i,j)==1) wetV = wetV + 1
      enddo
    enddo

    !Set Parameter Values
    ui = xi_u               ! u-grid dimension in x direction
    uj = eta_u              ! u-grid dimension in y direction
    vi = xi_v               ! v-grid dimension in x direction
    vj = eta_v              ! v-grid dimension in y direction
    rho_nodes = nR          ! number of rho points (vi*uj)
    u_nodes = nU            ! number of u points   (ui*uj)
    v_nodes = nV            ! number of v points   (vi*vj)
    max_rho_elements = maxR ! Number of elements using four rho points as corners
    max_u_elements = maxU   ! Number of elements using four  u  points as corners
    max_v_elements = maxV   ! Number of elements using four  v  points as corners
    rho_elements = wetR     ! rho elements with at least one corner with mask_rho=1
    u_elements = wetU       !  u  elements with at least one corner with mask_u=1
    v_elements = wetV       !  v  elements with at least one corner with mask_v=1


    !If IOSTAT is present, set return value to error code
    IF(PRESENT(IOSTAT)) IOSTAT = err
    !  0=No Errors                 30=Error allocating arrays
    ! 10=Error Opening NCgridfile  40=Error getting variables
    ! 20=Error getting dimensions  50=Error Closing NCgridfile

  END SUBROUTINE gridData

END MODULE PARAM_MOD 

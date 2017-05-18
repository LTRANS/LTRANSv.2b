MODULE HYDRO_MOD

!  This module handles all the input from the hydrodynamic NetCDF input files.
!  It is the only module that interacts with NetCDF input files.  It contains
!  all the variables read in from the NetCDF files.  It also contains all the
!  information and variables related to the grid elements.
!
!  Created by:            Zachary Schlag        
!  Created on:            07 Aug 2008
!  Last Modified on:         Feb 2013

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER :: iint, & !Keeps track of the input file, 0=file 1, 1=file 2, etc.
             stepf   !Keeps track of the forward time step

  !Used for reading in NetCDF variables one time step at a time
  INTEGER :: STARTr(4),COUNTr(4),STARTz(3),COUNTz(3)

  !Keeps track of the Rho, U, and V element that each particle is in
  INTEGER, ALLOCATABLE, DIMENSION(:) :: P_r_element,P_u_element,P_v_element

  !The Rho, U, and V nodes that make up the Rho, U, and V element that 
  !  the particle is in
  INTEGER :: rnode1,rnode2,rnode3,rnode4,unode1,unode2,unode3,unode4,vnode1,   &
             vnode2,vnode3,vnode4

  !These variables keep track of the interpolation method and weights
  INTEGER :: tOK
  DOUBLE PRECISION :: t,u,Wgt1,Wgt2,Wgt3,Wgt4

  !S-Level location variables
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SC,CS,SCW,CSW

  !Depth at each rho node location
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: depth
  
  !read in zeta,salinity,temperature,vertical diffusivity, and U,V,W velocities 
  !  at hydrodynamic back, center, and forward time
  INTEGER :: t_b,t_c,t_f,t_ijruv(12)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)  :: t_zeta
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: t_salt,t_temp,t_Kh,t_den,t_Uvel,      &
                                         t_Vvel,t_Wvel

  !Rho, U, and V grid wet elements(four node numbers that make up the element)
  !  (wet means at least 1 node is masked as water)
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: RE,UE,VE

  !For each element, a list containing itself and all the elements that share a 
  !  node with that element;  used to speed up determining which element the 
  !  particle has moved to, if it has moved at all
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: r_Adjacent,u_Adjacent,v_Adjacent

  !X/Y location of all the Rho,U,V grid nodes, and the angle between
  !  x-coordinate and true east direction (radian)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: rho_angle,rx,ry,ux,uy,vx,vy
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r_ele_x,r_ele_y,u_ele_x,    &
                                                   u_ele_y,v_ele_x,v_ele_y

  !U, and V grid metric node locations, and Rho grid masking
  !  These variables are shared with boundary_module.f90 to create the model 
  !  boundaries
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: x_u,y_u,x_v,y_v
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: m_r,m_u,m_v
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask_rho

  INTEGER, ALLOCATABLE, DIMENSION( : ) :: rho_mask,u_mask,v_mask   !ewn.v.2

  !Keeps track if the grid has been read in yet or not
  !  If the grid hasn't been read in, the boundaries can't be made
  LOGICAL :: GRD_SET = .FALSE.

  !The concatenated hydrodynamic input file name
  CHARACTER(len=200) :: filenm

  !Counters for NetCDF files
  INTEGER :: NCcount,NCstart,prcount

  !The following procedures have been made public:
  PUBLIC :: initGrid,initHydro,updateHydro,setEle,setEle_all,setInterp,        &
    getInterp,interp,WCTS_ITPI,getSlevel,getWlevel,getMask_Rho,getUVxy,        &
    getR_ele,getP_r_element,finHydro,initNetCDF,createNetCDF,writeNetCDF

CONTAINS

  SUBROUTINE initGrid()
    !This subroutine reads in the grid information and with it creates all the 
    !  element variables
    USE PARAM_MOD, ONLY: numpar,ui,vi,uj,vj,us,ws,rho_nodes,u_nodes,v_nodes,   &
        rho_elements,u_elements,v_elements,max_rho_elements,max_u_elements,    &
        max_v_elements,NCgridfile,prefix,suffix,filenum,numdigits
    USE CONVERT_MOD, ONLY: lon2x,lat2y
    USE netcdf
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

    INTEGER :: i,j,m,count,inele
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: romdepth,  &
                                                     x_rho,y_rho,angle
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask_u, mask_v
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: lon_rho,lat_rho,lon_u,    &
                                                     lat_u,lon_v,lat_v
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: r_ele,u_ele,v_ele
    !INTEGER, ALLOCATABLE, DIMENSION( : ) :: rho_mask,u_mask,v_mask  !ewn.v.2

    !ALLOCATE MODULE VARIABLES
    ALLOCATE(SC(us))
    ALLOCATE(CS(us))
    ALLOCATE(SCW(ws))
    ALLOCATE(CSW(ws))
    ALLOCATE(depth(rho_nodes))
    ALLOCATE(RE(4,rho_elements))
    ALLOCATE(UE(4,u_elements))
    ALLOCATE(VE(4,v_elements))
    ALLOCATE(r_Adjacent(rho_elements,10))
    ALLOCATE(u_Adjacent(u_elements,10))
    ALLOCATE(v_Adjacent(v_elements,10))
    ALLOCATE(rho_angle(rho_nodes))
    ALLOCATE(rx(rho_nodes))
    ALLOCATE(ry(rho_nodes))
    ALLOCATE(ux(u_nodes))
    ALLOCATE(uy(u_nodes))
    ALLOCATE(vx(v_nodes))
    ALLOCATE(vy(v_nodes))
    ALLOCATE(r_ele_x(4,rho_elements))
    ALLOCATE(r_ele_y(4,rho_elements))
    ALLOCATE(u_ele_x(4,u_elements))
    ALLOCATE(u_ele_y(4,u_elements))
    ALLOCATE(v_ele_x(4,v_elements))
    ALLOCATE(v_ele_y(4,v_elements))
    ALLOCATE(P_r_element(numpar))
    ALLOCATE(P_u_element(numpar))
    ALLOCATE(P_v_element(numpar))
    ALLOCATE(mask_rho(vi,uj))
    ALLOCATE(m_r(vi,uj))
    ALLOCATE(m_u(ui,uj))
    ALLOCATE(m_v(vi,vj))
    ALLOCATE(x_u(ui,uj))
    ALLOCATE(y_u(ui,uj))
    ALLOCATE(x_v(vi,vj))
    ALLOCATE(y_v(vi,vj))
    ALLOCATE(rho_mask(rho_nodes))  !ewn.v.2
    ALLOCATE(u_mask(u_nodes))
    ALLOCATE(v_mask(v_nodes))

    !ALLOCATE SUBROUTINE VARIABLES
    ALLOCATE(romdepth(vi,uj))
    ALLOCATE(mask_u(ui,uj))
    ALLOCATE(mask_v(vi,vj))
    ALLOCATE(x_rho(vi,uj))
    ALLOCATE(y_rho(vi,uj))
    ALLOCATE(lon_rho(vi,uj))
    ALLOCATE(lat_rho(vi,uj))
    ALLOCATE(lon_u(ui,uj))
    ALLOCATE(lat_u(ui,uj))
    ALLOCATE(lon_v(vi,vj))
    ALLOCATE(lat_v(vi,vj))
    ALLOCATE(angle(vi,uj))
    ALLOCATE(r_ele(4,max_rho_elements))
    ALLOCATE(u_ele(4,max_u_elements))
    ALLOCATE(v_ele(4,max_v_elements))
    !  ALLOCATE(rho_mask(rho_nodes))    !ewn.v.2
    !  ALLOCATE(u_mask(u_nodes))
    !  ALLOCATE(v_mask(v_nodes))


    WRITE(*,*) 'read-in grid information'

    ! *************************** READ IN GRID INFO **************************

    STATUS = NF90_OPEN(TRIM(NCgridfile),NF90_NOWRITE, NCID)
    if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
    if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! Depth (m)
      STATUS = NF90_INQ_VARID(NCID,'h',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find depth'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,romdepth)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read depth'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! longitude at rho (°)
      STATUS = NF90_INQ_VARID(NCID,'lon_rho',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lon_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lon_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! latitude at rho (°)
      STATUS = NF90_INQ_VARID(NCID,'lat_rho',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lat_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lat_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! longitude at u (°)
      STATUS = NF90_INQ_VARID(NCID,'lon_u',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lon_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lon_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! latitude at u (°)
      STATUS = NF90_INQ_VARID(NCID,'lat_u',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lat_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lat_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! longitude at v (°)
      STATUS = NF90_INQ_VARID(NCID,'lon_v',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lon_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lon_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! latitude at v (°)
      STATUS = NF90_INQ_VARID(NCID,'lat_v',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lat_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lat_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! mask on rho grid
      STATUS = NF90_INQ_VARID(NCID,'mask_rho',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find mask_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,mask_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! mask on u grid
      STATUS = NF90_INQ_VARID(NCID,'mask_u',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find mask_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,mask_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! mask on v grid
      STATUS = NF90_INQ_VARID(NCID,'mask_v',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find mask_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,mask_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! angle between x-coordinate and true east direction (radian)
      STATUS = NF90_INQ_VARID(NCID,'angle',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find angle'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,angle)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read angle'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

    STATUS = NF90_CLOSE(NCID)


    ! ************************* READ IN S LEVEL INFO *************************

    SELECT CASE(numdigits)
      CASE(1)
        WRITE(filenm,'(A,I1.1,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(2)
        WRITE(filenm,'(A,I2.2,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(3)
        WRITE(filenm,'(A,I3.3,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(4)
        WRITE(filenm,'(A,I4.4,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(5)
        WRITE(filenm,'(A,I5.5,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(6)
        WRITE(filenm,'(A,I6.6,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(7)
        WRITE(filenm,'(A,I7.7,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE(8)
        WRITE(filenm,'(A,I8.8,A)') TRIM(prefix),filenum,TRIM(suffix)
      CASE DEFAULT
        WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
        WRITE(*,*) 'Please use numdigit value from 1 to 8'
        WRITE(*,*) '  OR modify code in Hydrodynamic module'
        STOP
    END SELECT
!    write(*,*) TRIM(filenm)

    STATUS = NF90_OPEN(TRIM(filenm), NF90_NOWRITE, NCID)
    if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
    if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! s-coordinate on rho grid (sc_r)
      STATUS = NF90_INQ_VARID(NCID,'s_rho',VID)
      IF(STATUS /= NF90_NOERR)THEN
        STATUS = NF90_INQ_VARID(NCID,'sc_r',VID)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding SC'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      ENDIF
      STATUS = NF90_GET_VAR(NCID,VID,SC)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read SC'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! Cs value on rho grid (Cs_r)
      STATUS = NF90_INQ_VARID(NCID,'Cs_r',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find CS'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,CS)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read CS'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! s-coordinate on w grid (sc_w)
      STATUS = NF90_INQ_VARID(NCID,'s_w',VID)
      IF(STATUS /= NF90_NOERR)THEN
        STATUS = NF90_INQ_VARID(NCID,'sc_w',VID)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding SCW'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      ENDIF
      STATUS = NF90_GET_VAR(NCID,VID,SCW)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read SCW'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! Cs value on w grid (Cs_w)
      STATUS = NF90_INQ_VARID(NCID,'Cs_w',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find CSW'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,CSW)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read CSW'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

    !close the dataset and reassign the NCID
    STATUS = NF90_CLOSE(NCID)


    ! *************************** CREATE ELEMENTS *****************************

    !Store Mask Values to multiply by...
    m_r = mask_rho
    m_u = mask_u
    m_v = mask_v

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ~  4B. Prepare Elements (i.e., assign ID numbers to rectangular grids)  ~
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    write(*,*) 'create elements'

    ! Convert rho nodes lon/lat to x/y coordinates
    do j=1,uj
      do i=1,vi
        x_rho(i,j) = lon2x(lon_rho(i,j),lat_rho(i,j))
        y_rho(i,j) = lat2y(lat_rho(i,j))
      enddo
    enddo

    ! Convert u nodes lon/lat to x/y coordinates
    do j=1,uj
      do i=1,ui
        x_u(i,j) = lon2x(lon_u(i,j),lat_u(i,j))
        y_u(i,j) = lat2y(lat_u(i,j))
      enddo
    enddo

    ! Convert v nodes lon/lat to x/y coordinates
    do j=1,vj
      do i=1,vi
        x_v(i,j) = lon2x(lon_v(i,j),lat_v(i,j))
        y_v(i,j) = lat2y(lat_v(i,j))
      enddo
    enddo

    ! Assign mask values to rho nodes 
    count = 0
    do j=1,uj
      do i=1,vi
        count = count + 1   !move to next node number
            !cycles through each variable replacing the vi,uj part with count
            !  essentially giving it node numbers
        rho_mask(count) = mask_rho(i,j)
        rho_angle(count) = angle(i,j)
      enddo
    enddo

    ! Assign mask values to u nodes
    count = 0
    do j=1,uj
      do i=1,ui
        count = count + 1   !move to next node number
            !cycles through each variable replacing the vi,uj part with count
            !  essentially giving it node numbers
        u_mask(count) = mask_u(i,j)
      enddo
    enddo

    ! Assign mask values to v nodes
    count = 0
    do j=1,vj
      do i=1,vi
        count = count + 1   !move to next node number
            !cycles through each variable replacing the vi,uj part with count
            !  essentially giving it node numbers
        v_mask(count) = mask_v(i,j)
      enddo
    enddo


    ! Create matrix that contains the node numbers for each rho element
    count = 0
    do j=1,uj-1                         !z2v3.2
      do i=1,vi-1
        count = count + 1
        r_ele(1,count) = i + (j-1)*vi
        r_ele(2,count) = i + 1 + (j-1)*vi
        r_ele(3,count) = i + 1 + j*vi
        r_ele(4,count) = i + j*vi
      enddo
    enddo

    ! Create matrix that contains the node numbers for each u element
    count = 0
    do j=1,uj-1                         !z2v3.2
      do i=1,ui-1
        count = count + 1
        u_ele(1,count) = i + (j-1)*ui
        u_ele(2,count) = i + 1 + (j-1)*ui
        u_ele(3,count) = i + 1 + j*ui
        u_ele(4,count) = i + j*ui
      enddo
    enddo

    ! Create matrix that contains the node numbers for each v element
    count = 0
    do j=1,vj-1                         !z2v3.2
      do i=1,vi-1
        count = count + 1
        v_ele(1,count) = i + (j-1)*vi
        v_ele(2,count) = i + 1 + (j-1)*vi
        v_ele(3,count) = i + 1 + j*vi
        v_ele(4,count) = i + j*vi
      enddo
    enddo


    ! Create matrix that contains only the rho elements that contain a node
    ! whose mask value = 1 (i.e., it has at least one water point). 
    count = 0
    do i=1,max_rho_elements
      inele = 0
      !using the mask determine if any of the nodes for the current
      !  element are inbounds, if so set inele to 1
      if( rho_mask(r_ele(1,i)) .EQ. 1) inele=1
      if( rho_mask(r_ele(2,i)) .EQ. 1) inele=1
      if( rho_mask(r_ele(3,i)) .EQ. 1) inele=1
      if( rho_mask(r_ele(4,i)) .EQ. 1) inele=1              !z2v3.2
      !if inele = 1 then at least one of the three nodes for this element
      !  are in bounds.
      if( inele .EQ. 1 ) then
        count = count + 1
        !create array of elements that contain at least one node in bounds
        RE(1,count) = r_ele(1,i)
        RE(2,count) = r_ele(2,i)
        RE(3,count) = r_ele(3,i)
        RE(4,count) = r_ele(4,i)                            !z2v3.2
      endif
    enddo

    ! Create matrix that contains only the u elements that contain a node
    ! whose mask value = 1 (i.e., it has at least one water point). 
    count = 0
    do i=1,max_u_elements
      inele = 0
      !using the mask determine if any of the nodes for the current
      !  element are inbounds, if so set inele to 1
      if( u_mask(u_ele(1,i)) .EQ. 1) inele=1
      if( u_mask(u_ele(2,i)) .EQ. 1) inele=1
      if( u_mask(u_ele(3,i)) .EQ. 1) inele=1
      if( u_mask(u_ele(4,i)) .EQ. 1) inele=1                    !z2v3.2
      !if inele = 1 then at least one of the three nodes for this element
      !  are in bounds.
      if( inele .EQ. 1 ) then
        count = count + 1
        !create array of elements that contain at least one node in bounds
        UE(1,count) = u_ele(1,i)
        UE(2,count) = u_ele(2,i)
        UE(3,count) = u_ele(3,i)
        UE(4,count) = u_ele(4,i)                            !z2v3.2
      endif 
    enddo

    ! Create matrix that contains only the v elements that contain a node
    ! whose mask value = 1 (i.e., it has at least one water point). 
    count = 0
    do i=1,max_v_elements
      inele = 0
      !using the mask determine if any of the nodes for the current
      !  element are inbounds, if so set inele to 1
      if( v_mask(v_ele(1,i)) .EQ. 1) inele=1
      if( v_mask(v_ele(2,i)) .EQ. 1) inele=1
      if( v_mask(v_ele(3,i)) .EQ. 1) inele=1
      if( v_mask(v_ele(4,i)) .EQ. 1) inele=1                    !z2v3.2
      !if inele = 1 then at least one of the three nodes for this element
      !  are in bounds.
      if( inele .EQ. 1 ) then
        count = count + 1
        !create array of elements that contain at least one node in bounds
        VE(1,count) = v_ele(1,i)
        VE(2,count) = v_ele(2,i)
        VE(3,count) = v_ele(3,i)
        VE(4,count) = v_ele(4,i)                            !z2v3.2
      endif 
    enddo

    ! Create matrices of  x/y for rho nodes and depth values in rho node number
    !   format 
    count = 0
    do j=1,uj
      do i=1,vi
        count = count + 1   !move to next node number
        !cycles through each variable replacing the vi,uj part with count
        !  essentially giving it node numbers
        rx(count) = x_rho(i,j)
        ry(count) = y_rho(i,j)
        depth(count) = romdepth(i,j)
      enddo
    enddo

    ! Create matrices of x/y values for u nodes in u node number format 
    count = 0
    do j=1,uj
      do i=1,ui
        count = count + 1   !move to next node number
        !cycles through each variable replacing the ui,uj part with count
        !  essentially giving it node numbers
        ux(count) = x_u(i,j)
        uy(count) = y_u(i,j)
      enddo
    enddo

    ! Create matrices of x/y values for v nodes in v node number format
    count = 0
    do j=1,vj
      do i=1,vi
        count = count + 1   !move to next node number
        !cycles through each variable replacing the vi,vj part with count
        !  essentially giving it node numbers
        vx(count) = x_v(i,j)
        vy(count) = y_v(i,j)
      enddo
    enddo


    ! Create matrices of x/y node values for each rho, u, and v element
    do j=1,rho_elements
      do i=1,4                                          !z2v3.3
        r_ele_x(i,j) = rx(RE(i,j))
        r_ele_y(i,j) = ry(RE(i,j))
      enddo
    enddo

    do j=1,u_elements
      do i=1,4                                          !z2v3.3
        u_ele_x(i,j) = ux(UE(i,j))
        u_ele_y(i,j) = uy(UE(i,j))
      enddo
    enddo

    do j=1,v_elements
      do i=1,4                                          !z2v3.3
        v_ele_x(i,j) = vx(VE(i,j))
        v_ele_y(i,j) = vy(VE(i,j))
      enddo
    enddo


    ! ************************ FIND ADJACENT ELEMENTS *************************

    ! Create search restriction algorithms
    write(*,*) 'find adjacent elements'

    ! I. For each element, list all elements that are adjacent to it 
    write(*,*) ' - rho'
    do i=1,rho_elements                         !z2v4.2 - start
      r_Adjacent(i,1) = i
      m=1
      do j=max(i-(vi+2),1),min(i+vi+2,rho_elements)
        if(j.EQ.i) cycle
        if(  (RE(1,i) .EQ. RE(1,j)) .OR. (RE(1,i) .EQ. RE(2,j))                &
        .OR. (RE(1,i) .EQ. RE(3,j)) .OR. (RE(1,i) .EQ. RE(4,j))                &
        .OR. (RE(2,i) .EQ. RE(1,j)) .OR. (RE(2,i) .EQ. RE(2,j))                &
        .OR. (RE(2,i) .EQ. RE(3,j)) .OR. (RE(2,i) .EQ. RE(4,j))                &
        .OR. (RE(3,i) .EQ. RE(1,j)) .OR. (RE(3,i) .EQ. RE(2,j))                &
        .OR. (RE(3,i) .EQ. RE(3,j)) .OR. (RE(3,i) .EQ. RE(4,j))                &
        .OR. (RE(4,i) .EQ. RE(1,j)) .OR. (RE(4,i) .EQ. RE(2,j))                &
        .OR. (RE(4,i) .EQ. RE(3,j)) .OR. (RE(4,i) .EQ. RE(4,j))  )then
          m=m+1
          r_Adjacent(i,m) = j
        endif
      enddo 
    enddo 

    write(*,*) ' - u'
    do i=1,u_elements
      u_Adjacent(i,1) = i
      m=1
      do j=max(i-(ui+2),1),min(i+ui+2,u_elements)
        if(j.EQ.i) cycle
        if(  (UE(1,i) .EQ. UE(1,j)) .OR. (UE(1,i) .EQ. UE(2,j))                &
        .OR. (UE(1,i) .EQ. UE(3,j)) .OR. (UE(1,i) .EQ. UE(4,j))                &
        .OR. (UE(2,i) .EQ. UE(1,j)) .OR. (UE(2,i) .EQ. UE(2,j))                &
        .OR. (UE(2,i) .EQ. UE(3,j)) .OR. (UE(2,i) .EQ. UE(4,j))                &
        .OR. (UE(3,i) .EQ. UE(1,j)) .OR. (UE(3,i) .EQ. UE(2,j))                &
        .OR. (UE(3,i) .EQ. UE(3,j)) .OR. (UE(3,i) .EQ. UE(4,j))                &
        .OR. (UE(4,i) .EQ. UE(1,j)) .OR. (UE(4,i) .EQ. UE(2,j))                &
        .OR. (UE(4,i) .EQ. UE(3,j)) .OR. (UE(4,i) .EQ. UE(4,j))  )then
          m=m+1
          u_Adjacent(i,m) = j
        endif
      enddo 
    enddo 

    write(*,*) ' - v'
    do i=1,v_elements
      v_Adjacent(i,1) = i
      m=1
      do j=max(i-(vi+2),1),min(i+vi+2,v_elements)
        if(j.EQ.i) cycle
        if(  (VE(1,i) .EQ. VE(1,j)) .OR. (VE(1,i) .EQ. VE(2,j))                &
        .OR. (VE(1,i) .EQ. VE(3,j)) .OR. (VE(1,i) .EQ. VE(4,j))                &
        .OR. (VE(2,i) .EQ. VE(1,j)) .OR. (VE(2,i) .EQ. VE(2,j))                &
        .OR. (VE(2,i) .EQ. VE(3,j)) .OR. (VE(2,i) .EQ. VE(4,j))                &
        .OR. (VE(3,i) .EQ. VE(1,j)) .OR. (VE(3,i) .EQ. VE(2,j))                &
        .OR. (VE(3,i) .EQ. VE(3,j)) .OR. (VE(3,i) .EQ. VE(4,j))                &
        .OR. (VE(4,i) .EQ. VE(1,j)) .OR. (VE(4,i) .EQ. VE(2,j))                &
        .OR. (VE(4,i) .EQ. VE(3,j)) .OR. (VE(4,i) .EQ. VE(4,j))  )then     
          m=m+1
          v_Adjacent(i,m) = j
        endif
      enddo 
    enddo

    GRD_SET = .TRUE.

    !DEALLOCATE SUBROUTINE VARIABLES
    DEALLOCATE(romdepth,mask_u,mask_v,x_rho,y_rho,angle)
    ! DEALLOCATE(r_ele,u_ele,v_ele,rho_mask,u_mask,v_mask)
    DEALLOCATE(r_ele,u_ele,v_ele)                              !ewn.v.2


  END SUBROUTINE initGrid



  SUBROUTINE initHydro()
    !This Subroutine reads in the hydrodynamic information for the first 
    !  iteration
    USE PARAM_MOD, ONLY: numpar,ui,vi,uj,vj,us,ws,rho_nodes,u_nodes,v_nodes,   &
        prefix,suffix,filenum,numdigits,readZeta,constZeta,readSalt,constSalt, &
        readTemp,constTemp,readDens,constDens,readU,constU,readV,constV,readW, &
        constW,readAks,constAks
    USE netcdf
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

    INTEGER :: i,j,k,count,counter

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) :: romZ
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: romW,romKH,romS,romT,romD,romU,   &
                                             romV

    !ALLOCATE MODULE VARIABLES
    ALLOCATE(t_zeta(3,rho_nodes))
    ALLOCATE(t_salt(3,rho_nodes,us))
    ALLOCATE(t_temp(3,rho_nodes,us))
    ALLOCATE(t_Wvel(3,rho_nodes,ws))
    ALLOCATE(t_den (3,rho_nodes,us))
    ALLOCATE(t_Kh  (3,rho_nodes,ws))
    ALLOCATE(t_Uvel(3,  u_nodes,us))
    ALLOCATE(t_Vvel(3,  v_nodes,us))

    t_zeta = 0
    t_salt = 0
    t_temp = 0
    t_den  = 0
    t_Kh   = 0
    t_Uvel = 0
    t_Vvel = 0
    t_Wvel = 0

    !ALLOCATE SUBROUTINE VARIABLES
    ALLOCATE(romZ(vi,uj,3))
    ALLOCATE(romW(vi,uj,ws,3))
    ALLOCATE(romS(vi,uj,us,3))
    ALLOCATE(romT(vi,uj,us,3))
    ALLOCATE(romD(vi,uj,us,3))
    ALLOCATE(romU(ui,uj,us,3))
    ALLOCATE(romV(vi,vj,us,3))
    ALLOCATE(romKH(vi,uj,ws,3))


      !Open netCDF file
      counter=iint+filenum  !176 + 1 = 177 --> June 26,1995

      SELECT CASE(numdigits)
        CASE(1)
          WRITE(filenm,'(A,I1.1,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(2)
          WRITE(filenm,'(A,I2.2,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(3)
          WRITE(filenm,'(A,I3.3,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(4)
          WRITE(filenm,'(A,I4.4,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(5)
          WRITE(filenm,'(A,I5.5,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(6)
          WRITE(filenm,'(A,I6.6,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(7)
          WRITE(filenm,'(A,I7.7,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(8)
          WRITE(filenm,'(A,I8.8,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE DEFAULT
          WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
          WRITE(*,*) 'Please use numdigit value from 1 to 8'
          WRITE(*,*) '  OR modify code in Hydrodynamic module'
          STOP
      END SELECT
      write(*,*) TRIM(filenm)

      !t_ijruv = (/175,195,155,175,175,195,155,175,175,195,155,175/)

      t_b = 1    !Back step is 1st time step in arrays
      t_c = 2    !Center step is 2nd time step in arrays
      t_f = 3    !Forward step is 3rd time step in arrays
      stepf=3    !Forward step is 3rd time step of file

      !Get i/j max/min for rho/u/v
      call setijruv()

      ! Read in data for first three external time steps
      STATUS = NF90_OPEN(TRIM(filenm), NF90_NOWRITE, NCID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
!       write(*,*) 'ncid = ', ncid

      if(readZeta)then  
        ! **** Zeta ****
        startz(1)=t_ijruv(1)
        startz(2)=t_ijruv(3)
        startz(3)=1

        countz(1)=t_ijruv(2)-t_ijruv(1)+1
        countz(2)=t_ijruv(4)-t_ijruv(3)+1
        countz(3)=3

        STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find zeta'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romZ(t_ijruv(1):t_ijruv(2),             &
                              t_ijruv(3):t_ijruv(4),:),STARTz,COUNTz)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read zeta array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romZ = constZeta
      endif

      if(readSalt)then
        ! **** Salt ****
        startr(1)=t_ijruv(1)
        startr(2)=t_ijruv(3)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(2)-t_ijruv(1)+1
        countr(2)=t_ijruv(4)-t_ijruv(3)+1
        countr(3)=us
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'salt',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find salt'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romS(t_ijruv(1):t_ijruv(2),             &
                              t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read salt array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romS = constSalt
      endif

      if(readTemp)then  
        ! **** Temp ****
        startr(1)=t_ijruv(1)
        startr(2)=t_ijruv(3)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(2)-t_ijruv(1)+1
        countr(2)=t_ijruv(4)-t_ijruv(3)+1
        countr(3)=us
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'temp',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find temp'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romT(t_ijruv(1):t_ijruv(2),             &
                              t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read temp array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romT = constTemp
      endif

      if(readDens)then  
        ! **** Density ****
        startr(1)=t_ijruv(1)
        startr(2)=t_ijruv(3)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(2)-t_ijruv(1)+1
        countr(2)=t_ijruv(4)-t_ijruv(3)+1
        countr(3)=us
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'rho',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find rho'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romD(t_ijruv(1):t_ijruv(2),             &
                              t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read rho array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romD = constDens
      endif

      if(readU)then  
        ! **** U velocity ****
        startr(1)=t_ijruv(5)
        startr(2)=t_ijruv(7)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(6)-t_ijruv(5)+1
        countr(2)=t_ijruv(8)-t_ijruv(7)+1
        countr(3)=us
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'u',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find u'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romU(t_ijruv(5):t_ijruv(6),             &
                              t_ijruv(7):t_ijruv(8),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read u array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romU = constU
      endif

      if(readV)then  
        ! **** V velocity ****
        startr(1)=t_ijruv(9)
        startr(2)=t_ijruv(11)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(10)-t_ijruv(9)+1
        countr(2)=t_ijruv(12)-t_ijruv(11)+1
        countr(3)=us
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'v',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find v'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romV(t_ijruv(9):t_ijruv(10),            &
                              t_ijruv(11):t_ijruv(12),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read v array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romV = constV
      endif

      if(readW)then  
        ! **** W velocity ****
        startr(1)=t_ijruv(1)
        startr(2)=t_ijruv(3)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(2)-t_ijruv(1)+1
        countr(2)=t_ijruv(4)-t_ijruv(3)+1
        countr(3)=ws
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'w',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find w'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romW(t_ijruv(1):t_ijruv(2),             &
                              t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read w array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romW = constW
      endif

      if(readAks)then  
        ! **** Vertical diffusivity for salt (Aks) ****
        startr(1)=t_ijruv(1)
        startr(2)=t_ijruv(3)
        startr(3)=1
        startr(4)=1

        countr(1)=t_ijruv(2)-t_ijruv(1)+1
        countr(2)=t_ijruv(4)-t_ijruv(3)+1
        countr(3)=ws
        countr(4)=3
        STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem find AKs'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif

        STATUS = NF90_GET_VAR(NCID,VID,romKH(t_ijruv(1):t_ijruv(2),            &
                              t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
        if (STATUS .NE. NF90_NOERR) then
          write(*,*) 'Problem read AKs array'
          write(*,*) NF90_STRERROR(STATUS)
          stop
        endif
      else
        romKH = constAks
      endif

      !close the dataset and reassign the NCID
      STATUS = NF90_CLOSE(NCID)


      !Reshape input to fit node numbers assigned to elements
      do j=t_ijruv(3),t_ijruv(4)
        do i=t_ijruv(1),t_ijruv(2)
          count = (j-1)*vi + i
          do k=1,us
            t_salt(1,count,k) = romS(i,j,k,1) * m_r(i,j)
            t_salt(2,count,k) = romS(i,j,k,2) * m_r(i,j)
            t_salt(3,count,k) = romS(i,j,k,3) * m_r(i,j)
            t_temp(1,count,k) = romT(i,j,k,1) * m_r(i,j)
            t_temp(2,count,k) = romT(i,j,k,2) * m_r(i,j)
            t_temp(3,count,k) = romT(i,j,k,3) * m_r(i,j)
            t_Wvel(1,count,k) = romW(i,j,k,1) * m_r(i,j)
            t_Wvel(2,count,k) = romW(i,j,k,2) * m_r(i,j)
            t_Wvel(3,count,k) = romW(i,j,k,3) * m_r(i,j)
            t_den(1,count,k) = (romD(i,j,k,1) + DBLE(1000.0)) * m_r(i,j)
            t_den(2,count,k) = (romD(i,j,k,2) + DBLE(1000.0)) * m_r(i,j)
            t_den(3,count,k) = (romD(i,j,k,3) + DBLE(1000.0)) * m_r(i,j)
            t_KH(1,count,k) = romKH(i,j,k,1) * m_r(i,j)
            t_KH(2,count,k) = romKH(i,j,k,2) * m_r(i,j)
            t_KH(3,count,k) = romKH(i,j,k,3) * m_r(i,j)
          enddo
          t_Wvel(1,count,ws) = romW(i,j,ws,1) * m_r(i,j)
          t_Wvel(2,count,ws) = romW(i,j,ws,2) * m_r(i,j)
          t_Wvel(3,count,ws) = romW(i,j,ws,3) * m_r(i,j)
          t_KH(1,count,ws) = romKH(i,j,ws,1) * m_r(i,j)
          t_KH(2,count,ws) = romKH(i,j,ws,2) * m_r(i,j)
          t_KH(3,count,ws) = romKH(i,j,ws,3) * m_r(i,j)
          t_zeta(1,count) = romZ(i,j,1) * m_r(i,j)
          t_zeta(2,count) = romZ(i,j,2) * m_r(i,j)
          t_zeta(3,count) = romZ(i,j,3) * m_r(i,j)
        enddo
      enddo

      do j=t_ijruv(7),t_ijruv(8)
        do i=t_ijruv(5),t_ijruv(6)
          count = (j-1)*ui + i
          do k=1,us
            t_Uvel(1,count,k) = romU(i,j,k,1) * m_u(i,j)
            t_Uvel(2,count,k) = romU(i,j,k,2) * m_u(i,j)
            t_Uvel(3,count,k) = romU(i,j,k,3) * m_u(i,j)
          enddo
        enddo
      enddo

      do j=t_ijruv(11),t_ijruv(12)
        do i=t_ijruv(9),t_ijruv(10)
          count = (j-1)*vi + i
          do k=1,us
            t_Vvel(1,count,k) = romV(i,j,k,1) * m_v(i,j)
            t_Vvel(2,count,k) = romV(i,j,k,2) * m_v(i,j)
            t_Vvel(3,count,k) = romV(i,j,k,3) * m_v(i,j)
          enddo
        enddo    
      enddo


    !DEALLOCATE SUBROUTINE VARIABLES
    DEALLOCATE(romZ,romW,romD,romKH,romS,romT,romU,romV)

  END SUBROUTINE initHydro


  SUBROUTINE updateHydro()
    USE PARAM_MOD, ONLY: ui,vi,uj,vj,us,ws,tdim,rho_nodes,u_nodes,v_nodes,     &
        prefix,suffix,filenum,numdigits,readZeta,constZeta,readSalt,constSalt, &
        readTemp,constTemp,readDens,constDens,readU,constU,readV,constV,readW, &
        constW,readAks,constAks,startfile
    USE netcdf
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

    INTEGER :: i,j,k,count,counter
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) :: romZf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: romSf,romTf,romDf,romUf,romVf,    &
                                             romWf,romKHf

    !ALLOCATE SUBROUTINE VARIABLES
    ALLOCATE(romZf(vi,uj,1))
    ALLOCATE(romSf(vi,uj,us,1))
    ALLOCATE(romTf(vi,uj,us,1))
    ALLOCATE(romDf(vi,uj,us,1))
    ALLOCATE(romUf(ui,uj,us,1))
    ALLOCATE(romVf(vi,vj,us,1))
    ALLOCATE(romWf(vi,uj,ws,1))
    ALLOCATE(romKHf(vi,uj,ws,1))

    !Rotate Indices
    t_b = mod(t_b,3)+1  ! 1 -> 2 -> 3 -> 1
    t_c = mod(t_c,3)+1  ! 2 -> 3 -> 1 -> 2
    t_f = mod(t_f,3)+1  ! 3 -> 1 -> 2 -> 3


    !if the current input file is not yet finished, just increment stepf to 
    !  the next time step
    IF ( (startfile .AND. (iint==0) .AND. (stepf==tdim)) .OR.   &
         (stepf .LT. tdim)                                      ) THEN

      stepf=stepf+1

    ELSE
    !if the current input file is finished, update filenm to next input file,
    !  and reset stepf to 1

      !Open netCDF file
      iint = iint+1
      counter=iint+filenum  !176 + 1 = 177 --> June 26,1995

      SELECT CASE(numdigits)
        CASE(1)
          WRITE(filenm,'(A,I1.1,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(2)
          WRITE(filenm,'(A,I2.2,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(3)
          WRITE(filenm,'(A,I3.3,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(4)
          WRITE(filenm,'(A,I4.4,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(5)
          WRITE(filenm,'(A,I5.5,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(6)
          WRITE(filenm,'(A,I6.6,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(7)
          WRITE(filenm,'(A,I7.7,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE(8)
          WRITE(filenm,'(A,I8.8,A)') TRIM(prefix),counter,TRIM(suffix)
        CASE DEFAULT
          WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
          WRITE(*,*) 'Please use numdigit value from 1 to 8'
          WRITE(*,*) '  OR modify code in Hydrodynamic module'
          STOP
      END SELECT
      write(*,*) TRIM(filenm)

      stepf = 1

    ENDIF


    !Get i/j max/min for rho/u/v
    call setijruv()


    ! Read in forward time step data 
    STATUS = NF90_OPEN(TRIM(filenm), NF90_NOWRITE, NCID)
    if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
    if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
!   write(*,*) 'ncid = ', ncid

    if(readZeta)then
      ! **** Zeta ****
      startz(1)=t_ijruv(1)
      startz(2)=t_ijruv(3)
      startz(3)=stepf

      countz(1)=t_ijruv(2)-t_ijruv(1)+1
      countz(2)=t_ijruv(4)-t_ijruv(3)+1
      countz(3)=1

      STATUS = NF90_INQ_VARID(NCID,'zeta',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find zeta f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romZf(t_ijruv(1):t_ijruv(2),              &
                            t_ijruv(3):t_ijruv(4),:),STARTz,COUNTz)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read zeta array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romZf = constZeta
    endif

    if(readU)then
      ! **** U velocity ****
      startr(1)=t_ijruv(5)
      startr(2)=t_ijruv(7)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(6)-t_ijruv(5)+1
      countr(2)=t_ijruv(8)-t_ijruv(7)+1
      countr(3)=us
      countr(4)=1

      STATUS = NF90_INQ_VARID(NCID,'u',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find u f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romUf(t_ijruv(5):t_ijruv(6),              &
                            t_ijruv(7):t_ijruv(8),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read u array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romUf = constU
    endif

    if(readV)then
      ! **** V velocity ****
      startr(1)=t_ijruv(9)
      startr(2)=t_ijruv(11)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(10)-t_ijruv(9)+1
      countr(2)=t_ijruv(12)-t_ijruv(11)+1
      countr(3)=us
      countr(4)=1
      STATUS = NF90_INQ_VARID(NCID,'v',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find v f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romVf(t_ijruv(9):t_ijruv(10),             &
                            t_ijruv(11):t_ijruv(12),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read v array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romVf = constV
    endif


    if(readSalt)then
      ! **** Salt ****
      startr(1)=t_ijruv(1)
      startr(2)=t_ijruv(3)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(2)-t_ijruv(1)+1
      countr(2)=t_ijruv(4)-t_ijruv(3)+1
      countr(3)=us
      countr(4)=1
      STATUS = NF90_INQ_VARID(NCID,'salt',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find salt f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romSf(t_ijruv(1):t_ijruv(2),              &
                            t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read salt array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romSf = constSalt
    endif

    if(readTemp)then
      ! **** Temp ****
      startr(1)=t_ijruv(1)
      startr(2)=t_ijruv(3)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(2)-t_ijruv(1)+1
      countr(2)=t_ijruv(4)-t_ijruv(3)+1
      countr(3)=us
      countr(4)=1
      STATUS = NF90_INQ_VARID(NCID,'temp',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find temp f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romTf(t_ijruv(1):t_ijruv(2),              &
                            t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read temp array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romTf = constTemp
    endif

    if(readDens)then
      ! **** Density ****
      startr(1)=t_ijruv(1)
      startr(2)=t_ijruv(3)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(2)-t_ijruv(1)+1
      countr(2)=t_ijruv(4)-t_ijruv(3)+1
      countr(3)=us
      countr(4)=1
      STATUS = NF90_INQ_VARID(NCID,'rho',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find rho f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romDf(t_ijruv(1):t_ijruv(2),              &
                            t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read rho array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romDf = constDens
    endif


    if(readW)then
      ! **** W velocity ****
      startr(1)=t_ijruv(1)
      startr(2)=t_ijruv(3)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(2)-t_ijruv(1)+1
      countr(2)=t_ijruv(4)-t_ijruv(3)+1
      countr(3)=ws
      countr(4)=1
      STATUS = NF90_INQ_VARID(NCID,'w',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find w f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romWf(t_ijruv(1):t_ijruv(2),              &
                            t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read w array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romWf = constW
    endif

    if(readAks)then
      ! **** Vertical diffusivity for salt (Aks) ****
      startr(1)=t_ijruv(1)
      startr(2)=t_ijruv(3)
      startr(3)=1
      startr(4)=stepf

      countr(1)=t_ijruv(2)-t_ijruv(1)+1
      countr(2)=t_ijruv(4)-t_ijruv(3)+1
      countr(3)=ws
      countr(4)=1
      STATUS = NF90_INQ_VARID(NCID,'AKs',VID)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem find AKs f'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
      STATUS = NF90_GET_VAR(NCID,VID,romKHf(t_ijruv(1):t_ijruv(2),             &
                            t_ijruv(3):t_ijruv(4),:,:),STARTr,COUNTr)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read AKs array'
        write(*,*) NF90_STRERROR(STATUS)
        stop
      endif
    else
      romKHf = constAks
    endif

    !close the dataset and reassign the NCID
    STATUS = NF90_CLOSE(NCID)


    !Reshape input to fit node numbers assigned to elements
    do j=t_ijruv(3),t_ijruv(4)
      do i=t_ijruv(1),t_ijruv(2)
        count = (j-1)*vi + i
        do k=1,us    
          t_salt(t_f,count,k) = romSf(i,j,k,1) * m_r(i,j)
          t_temp(t_f,count,k) = romTf(i,j,k,1) * m_r(i,j)
          t_Wvel(t_f,count,k) = romWf(i,j,k,1) * m_r(i,j)
          t_den(t_f,count,k)  = (romDf(i,j,k,1) + DBLE(1000.0)) * m_r(i,j)
          t_KH(t_f,count,k)   = romKHf(i,j,k,1) * m_r(i,j)
        enddo
        t_Wvel(t_f,count,ws)  = romWf(i,j,ws,1) * m_r(i,j)
        t_KH(t_f,count,ws)    = romKHf(i,j,ws,1) * m_r(i,j)
        t_zeta(t_f,count)     = romZf(i,j,1) * m_r(i,j)
      enddo
    enddo

    do j=t_ijruv(7),t_ijruv(8)
      do i=t_ijruv(5),t_ijruv(6)
        count = (j-1)*ui + i
        do k=1,us    
          t_Uvel(t_f,count,k) = romUf(i,j,k,1) * m_u(i,j)
        enddo
      enddo
    enddo

    do j=t_ijruv(11),t_ijruv(12)
      do i=t_ijruv(9),t_ijruv(10)
        count = (j-1)*vi + i
        do k=1,us    
          t_Vvel(t_f,count,k) = romVf(i,j,k,1) * m_v(i,j)
        enddo
      enddo    
    enddo

    !DEALLOCATE SUBROUTINE VARIABLES
    DEALLOCATE(romZf,romSf,romTf,romUf,romVf,romWf,romKHf)

    write(*,*) 'existing matrix,stepf=',stepf

  END SUBROUTINE updateHydro



  SUBROUTINE setEle(Xpar,Ypar,n,err,first)
    !This Subroutine determines which Rho, U, and V grid elements contain 
    !  the given particle
    USE PARAM_MOD, ONLY: rho_elements,u_elements,v_elements,ui,vi,us,ws
    USE GRIDCELL_MOD, ONLY: gridcell
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT), OPTIONAL :: err
    LOGICAL, INTENT(IN), OPTIONAL :: first

    LOGICAL :: fst
    INTEGER :: i,triangle,checkele,P_r_ele,P_u_ele,P_v_ele,oP_ele,P_ele,error

    error = 0

    if( PRESENT(first) ) then
      fst = first
    else
      fst = .FALSE.
    endif

    if(fst) then !if the first iteration

      !Find rho element in which particle is located
      P_r_ele=0
      triangle=0
      call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar,Ypar,P_r_ele,triangle)
      if (triangle.EQ.0) error = 1
      P_r_element(n)=P_r_ele

      !Find u element in which particle is located
      P_u_ele=0
      triangle=0
      call gridcell(u_elements,u_ele_y,u_ele_x,Xpar,Ypar,P_u_ele,triangle)
      if (triangle.EQ.0) error = 2
      P_u_element(n)=P_u_ele

      !Find v element in which particle is located
      P_v_ele=0
      triangle=0
      call gridcell(v_elements,v_ele_y,v_ele_x,Xpar,Ypar,P_v_ele,triangle)
      if (triangle.EQ.0) error = 3
      P_v_element(n)=P_v_ele


    else !if not the first iteration 

      !Find rho element in which particle is located
      oP_ele = P_r_element(n)
      do i=1,10
        if(r_Adjacent(oP_ele,i).EQ.0) error = 4

        triangle = 0
        checkele = r_Adjacent(oP_ele,i)
        call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar,Ypar,P_ele,triangle,   &
                      checkele)
        if(triangle .NE. 0) then
          P_r_element(n) = P_ele
          exit
        endif

      enddo !r_singlecellloop


      !Find u element in which particle is located
      oP_ele = P_u_element(n)
      do i=1,10
        if(u_Adjacent(oP_ele,i).EQ.0) error = 5

        triangle = 0
        checkele = u_Adjacent(oP_ele,i)
        call gridcell(u_elements,u_ele_y,u_ele_x,Xpar,Ypar,P_ele,triangle,     &
                      checkele)
        if(triangle .NE. 0) then
          P_u_element(n) = P_ele
          exit
        endif

      enddo


      !Find v element in which particle is located
      oP_ele = P_v_element(n)
      do i=1,10
        if(v_Adjacent(oP_ele,i).EQ.0) error = 6

        triangle = 0
        checkele = v_Adjacent(oP_ele,i)
        call gridcell(v_elements,v_ele_y,v_ele_x,Xpar,Ypar,P_ele,triangle,     &
                      checkele)
        if(triangle .NE. 0) then
          P_v_element(n) = P_ele
          exit
        endif

      enddo

    endif

    !Assign node numbers for rho,u,v calculations
    rnode1 = RE(1,P_r_element(n))
    rnode2 = RE(2,P_r_element(n))
    rnode3 = RE(3,P_r_element(n))
    rnode4 = RE(4,P_r_element(n))

    unode1 = UE(1,P_u_element(n))
    unode2 = UE(2,P_u_element(n))
    unode3 = UE(3,P_u_element(n))
    unode4 = UE(4,P_u_element(n))

    vnode1 = VE(1,P_v_element(n))
    vnode2 = VE(2,P_v_element(n))
    vnode3 = VE(3,P_v_element(n))
    vnode4 = VE(4,P_v_element(n))

    if(PRESENT(err)) err = error

  END SUBROUTINE setEle



  SUBROUTINE setEle_all(Xpar,Ypar,err,par)
    !This Subroutine determines which Rho, U, and V grid elements contain
    !  each particle
    USE PARAM_MOD, ONLY: numpar,rho_elements,u_elements,v_elements,ui,vi,us,ws
    USE GRIDCELL_MOD, ONLY: gridcell
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN), DIMENSION(numpar) :: Xpar,Ypar
    INTEGER, INTENT(OUT), OPTIONAL :: err,par ! err = error code
                                              ! par = problem particle

    INTEGER :: n,p,i,error
    INTEGER :: triangle,checkele,P_r_ele,P_u_ele,P_v_ele,oP_ele,P_ele

    error = 0
    p = 0

    !Find rho element in which first particle is located
    P_r_ele=0
    triangle=0
    call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar(1),Ypar(1),P_r_ele,triangle)
    if (triangle.EQ.0) then
      error = 1
      p = 1
    endif
    P_r_element(1)=P_r_ele

    !Find u element in which first particle is located
    P_u_ele=0
    triangle=0
    call gridcell(u_elements,u_ele_y,u_ele_x,Xpar(1),Ypar(1),P_u_ele,triangle)
    if (triangle.EQ.0) then
      error = 2
      p = 1
    endif
    P_u_element(1)=P_u_ele

    !Find v element in which first particle is located
    P_v_ele=0
    triangle=0
    call gridcell(v_elements,v_ele_y,v_ele_x,Xpar(1),Ypar(1),P_v_ele,triangle)
    if (triangle.EQ.0) then
      error = 3
      p = 1
    endif
    P_v_element(1)=P_v_ele

    !Find rho, u, and v elements in which subsequent particles are located
    parloop: do n=2,numpar

      !Find rho element in which particle is located
      oP_ele = P_r_element(n-1)
      do i=1,10
        if(r_Adjacent(oP_ele,i).EQ.0) then
          !If selective search based on previous particle location fails, 
          !  search all
          P_r_ele=0
          triangle=0
          call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar(n),Ypar(n),P_r_ele,  &
                        triangle)
          if (triangle.EQ.0) then
            error = 1
            p = n
            exit parloop
          endif
          P_r_element(n)=P_r_ele
        else
          triangle = 0
          checkele = r_Adjacent(oP_ele,i)
          call gridcell(rho_elements,r_ele_y,r_ele_x,Xpar(n),Ypar(n),P_ele,    &
                        triangle,checkele)
          if(triangle .NE. 0) then
            P_r_element(n) = P_ele
            exit
          endif
        endif
      enddo


      !Find u element in which particle is located
      oP_ele = P_u_element(n-1)
      do i=1,10
        if(u_Adjacent(oP_ele,i).EQ.0) then
          !If selective search based on previous particle location fails, 
          !  search all
          P_u_ele=0
          triangle=0
          call gridcell(u_elements,u_ele_y,u_ele_x,Xpar(n),Ypar(n),P_u_ele,    &
                        triangle)
          if (triangle.EQ.0) then
            error = 2
            p = n
            exit parloop
          endif
          P_u_element(n)=P_u_ele
        else
          triangle = 0
          checkele = u_Adjacent(oP_ele,i)
          call gridcell(u_elements,u_ele_y,u_ele_x,Xpar(n),Ypar(n),P_ele,      &
                        triangle,checkele)
          if(triangle .NE. 0) then
            P_u_element(n) = P_ele
            exit
          endif
        endif
      enddo


      !Find v element in which particle is located
      oP_ele = P_v_element(n-1)
      do i=1,10
        if(v_Adjacent(oP_ele,i).EQ.0) then
          !If selective search based on previous particle location fails, 
          !  search all
          P_v_ele=0
          triangle=0
          call gridcell(v_elements,v_ele_y,v_ele_x,Xpar(n),Ypar(n),P_v_ele,    &
                        triangle)
          if (triangle.EQ.0) then
            error = 3
            p = n
            exit parloop
          endif
          P_v_element(n)=P_v_ele
        else
          triangle = 0
          checkele = v_Adjacent(oP_ele,i)
          call gridcell(v_elements,v_ele_y,v_ele_x,Xpar(n),Ypar(n),P_ele,      &
                        triangle,checkele)
          if(triangle .NE. 0) then
            P_v_element(n) = P_ele
            exit
          endif
        endif
      enddo

    enddo parloop

    if(PRESENT(err)) err = error
    if(PRESENT(par)) par = p

  END SUBROUTINE setEle_all



  SUBROUTINE setInterp(xp,yp)
    !This subroutine calculates and stores the interpolation method and values 
    !  for the current particle
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: xp,yp

    double precision x1,x2,x3,x4,y1,y2,y3,y4
    double precision Dis1,Dis2,Dis3,Dis4,TDis

    x1 = rx(rnode1)
    x2 = rx(rnode2)
    x3 = rx(rnode3)
    x4 = rx(rnode4)
    y1 = ry(rnode1)
    y2 = ry(rnode2)
    y3 = ry(rnode3)
    y4 = ry(rnode4)

    Wgt1 = 0
    Wgt2 = 0
    Wgt3 = 0
    Wgt4 = 0

    tOK = 0 !to store information on the interpolation method

    ! bilinear interpolation of first triangle
    t = ((xp-x1)*(y3-y1)+(y1-yp)*(x3-x1)) / ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
    u = ((xp-x1)*(y2-y1)+(y1-yp)*(x2-x1)) / ((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))
    tOK = 1 !first triangle

    ! if outside triangle, then do bilinear interpolation of other triangle
    if( t.LT.0. .OR. u.LT.0. .OR. (t+u).GT.1.0 ) then
      t = ((xp-x3)*(y1-y3)+(y3-yp)*(x1-x3)) / ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))
      u = ((xp-x3)*(y4-y3)+(y3-yp)*(x4-x3)) / ((x1-x3)*(y4-y3)-(y1-y3)*(x4-x3))
      tOK = 2 !second triangle

      !if bilinear techniques are undefined, then use inverse weighted distance
      if( t.LT.0. .OR. u.LT.0. .OR. (t+u).GT.1.0 ) then
        !if particle on node, then set equal to node value
        if(  (xp.EQ.x1 .AND. yp.EQ.y1).OR.(xp.EQ.x2 .AND. yp.EQ.y2)            &
         .OR.(xp.EQ.x3 .AND. yp.EQ.y3).OR.(xp.EQ.x4 .AND. yp.EQ.y4))then
          if (xp.EQ.x1 .AND. yp.EQ.y1) Wgt1 = 1.0
          if (xp.EQ.x2 .AND. yp.EQ.y2) Wgt2 = 1.0
          if (xp.EQ.x3 .AND. yp.EQ.y3) Wgt3 = 1.0
          if (xp.EQ.x4 .AND. yp.EQ.y4) Wgt4 = 1.0
        else !use inverse weighted distance instead  
          Dis1=1./( SQRT( (x1-xp)**2 + (y1-yp)**2 ) ) 
          Dis2=1./( SQRT( (x2-xp)**2 + (y2-yp)**2 ) ) 
          Dis3=1./( SQRT( (x3-xp)**2 + (y3-yp)**2 ) ) 
          Dis4=1./( SQRT( (x4-xp)**2 + (y4-yp)**2 ) ) 
          TDis = Dis1+Dis2+Dis3+Dis4
          Wgt1= Dis1/TDis
          Wgt2= Dis2/TDis
          Wgt3= Dis3/TDis
          Wgt4= Dis4/TDis
          tOK = 3 !no triangle - used inverse weighted distance
        endif
      endif
    endif     

  END SUBROUTINE setInterp


  DOUBLE PRECISION FUNCTION getInterp(var,i)
    !This Function returns the interpolated value at the particle's location 
    !  using the interpolation variables stored from function setInterp, and
    !  the hydrodynamic variables that have been read in
    USE PARAM_MOD, ONLY: FreeSlip
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER, INTENT(IN), OPTIONAL :: i

    INTEGER :: rmasksum   !ewn.v.2

    DOUBLE PRECISION :: v1,v2,v3,v4

    !Determine which data to interpolate from
    SELECT CASE(var)
      CASE("depth")
        v1 = depth(rnode1)
        v2 = depth(rnode2)
        v3 = depth(rnode3)
        v4 = depth(rnode4)
      CASE("angle")
        v1 = rho_angle(rnode1)
        v2 = rho_angle(rnode2)
        v3 = rho_angle(rnode3)
        v4 = rho_angle(rnode4)
      CASE("zetab")
        v1 = t_zeta(t_b,rnode1)
        v2 = t_zeta(t_b,rnode2)
        v3 = t_zeta(t_b,rnode3)
        v4 = t_zeta(t_b,rnode4)
      CASE("zetac")
        v1 = t_zeta(t_c,rnode1)
        v2 = t_zeta(t_c,rnode2)
        v3 = t_zeta(t_c,rnode3)
        v4 = t_zeta(t_c,rnode4)
      CASE("zetaf")
        v1 = t_zeta(t_f,rnode1)
        v2 = t_zeta(t_f,rnode2)
        v3 = t_zeta(t_f,rnode3)
        v4 = t_zeta(t_f,rnode4)
      CASE("saltb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_salt(t_b,rnode1,i)
        v2 = t_salt(t_b,rnode2,i)
        v3 = t_salt(t_b,rnode3,i)
        v4 = t_salt(t_b,rnode4,i)
      CASE("saltc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_salt(t_c,rnode1,i)
        v2 = t_salt(t_c,rnode2,i)
        v3 = t_salt(t_c,rnode3,i)
        v4 = t_salt(t_c,rnode4,i)
      CASE("saltf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_salt(t_f,rnode1,i)
        v2 = t_salt(t_f,rnode2,i)
        v3 = t_salt(t_f,rnode3,i)
        v4 = t_salt(t_f,rnode4,i)
      CASE("tempb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_temp(t_b,rnode1,i)
        v2 = t_temp(t_b,rnode2,i)
        v3 = t_temp(t_b,rnode3,i)
        v4 = t_temp(t_b,rnode4,i)
      CASE("tempc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_temp(t_c,rnode1,i)
        v2 = t_temp(t_c,rnode2,i)
        v3 = t_temp(t_c,rnode3,i)
        v4 = t_temp(t_c,rnode4,i)
      CASE("tempf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_temp(t_f,rnode1,i)
        v2 = t_temp(t_f,rnode2,i)
        v3 = t_temp(t_f,rnode3,i)
        v4 = t_temp(t_f,rnode4,i)
      CASE("denb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_den(t_b,rnode1,i)
        v2 = t_den(t_b,rnode2,i)
        v3 = t_den(t_b,rnode3,i)
        v4 = t_den(t_b,rnode4,i)
      CASE("denc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_den(t_c,rnode1,i)
        v2 = t_den(t_c,rnode2,i)
        v3 = t_den(t_c,rnode3,i)
        v4 = t_den(t_c,rnode4,i)
      CASE("denf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_den(t_f,rnode1,i)
        v2 = t_den(t_f,rnode2,i)
        v3 = t_den(t_f,rnode3,i)
        v4 = t_den(t_f,rnode4,i)
      CASE("khb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_KH(t_b,rnode1,i)
        v2 = t_KH(t_b,rnode2,i)
        v3 = t_KH(t_b,rnode3,i)
        v4 = t_KH(t_b,rnode4,i)
      CASE("khc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_KH(t_c,rnode1,i)
        v2 = t_KH(t_c,rnode2,i)
        v3 = t_KH(t_c,rnode3,i)
        v4 = t_KH(t_c,rnode4,i)
      CASE("khf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_KH(t_f,rnode1,i)
        v2 = t_KH(t_f,rnode2,i)
        v3 = t_KH(t_f,rnode3,i)
        v4 = t_KH(t_f,rnode4,i)
      CASE DEFAULT
        write(*,*) 'Problem interpolating ',var
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop
    END SELECT

     !Ensure there is no friction near land (the free slip condition) !ewn.v.2
     !by setting values on land nodes equal to nearby water nodes
       IF (FreeSlip) THEN
       ! determine if a land element is in present
       rmasksum = 0
       rmasksum = rho_mask(rnode1)+ rho_mask(rnode2)+ rho_mask(rnode3)         &
                + rho_mask(rnode4)
       ! rmasksum, an integer, is the sum of the mask values of each node
       if (rmasksum .LT. 4) then  
         ! determine how many land nodes are present and reassign values 
         if (rmasksum .EQ. 3) then   !one land node
            if (rho_mask(rnode1) .EQ. 0) then
               v1 = 0.5*(v2+v4)
            else if (rho_mask(rnode2) .EQ. 0) then 
               v2 = 0.5*(v1+v3)
            else if (rho_mask(rnode3) .EQ. 0) then 
               v3 = 0.5*(v2+v4)
            else if (rho_mask(rnode4) .EQ. 0) then 
               v4 = 0.5*(v1+v3)
            end if
         else if (rmasksum .EQ. 2) then   !two land nodes
            if (rho_mask(rnode1).EQ.0 .AND. rho_mask(rnode2).EQ.0) then
                v1 = v4 
                v2 = v3
            else if (rho_mask(rnode2).EQ.0 .AND. rho_mask(rnode3).EQ.0) then
                v2 = v1 
                v3 = v4
            else if (rho_mask(rnode3).EQ.0 .AND. rho_mask(rnode4).EQ.0) then
                v3 = v2 
                v4 = v1
            else if (rho_mask(rnode4).EQ.0 .AND. rho_mask(rnode1).EQ.0) then
                v4 = v3 
                v1 = v2
            else if (rho_mask(rnode1).EQ.0 .AND. rho_mask(rnode3).EQ.0) then
            v1 = v4
            v3 = v2
            else if (rho_mask(rnode4).EQ.0 .AND. rho_mask(rnode2).EQ.0) then
            v4 = v1
            v2 = v3
            end if
         else if (rmasksum .EQ. 1) then   !three land nodes
            if (rho_mask(rnode1) .EQ. 1) then
                 v2 = v1
                 v3 = v1  
                 v4 = v1
            else if (rho_mask(rnode2) .EQ. 1) then
                 v1 = v2
                 v3 = v2  
                 v4 = v2
            else if (rho_mask(rnode3) .EQ. 1) then
                 v1 = v3
                 v2 = v3  
                 v4 = v3
            else if (rho_mask(rnode4) .EQ. 1) then
                 v1 = v4
                 v2 = v4  
                 v3 = v4
            end if
         end if
       end if
     ENDIF

    !interpolate using the variables from setInterp
    if(tOK == 1) then 
      getInterp = v1 + (v2-v1)*t + (v3-v1)*u
    elseif(tOK == 2) then
      getInterp = v3 + (v4-v3)*t + (v1-v3)*u
    else 
      getInterp = Wgt1*v1 + Wgt2*v2 + Wgt3*v3 + Wgt4*v4 
    endif

  END FUNCTION getInterp


  DOUBLE PRECISION FUNCTION interp(xp,yp,var,i)
    !This Function determines the method of interpolation and returns the 
    !  interpolated value at the particle's location using the hydrodynamic 
    !  variables that have been read in
    USE PARAM_MOD, ONLY: FreeSlip
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: xp,yp
    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER, INTENT(IN), OPTIONAL :: i

    INTEGER :: rmasksum,vmasksum,umasksum   !ewn.v.2

    INTEGER :: RUV
    DOUBLE PRECISION :: tt,uu,x1,x2,x3,x4,y1,y2,y3,y4,v1,v2,v3,v4,vp,     &
                        Dis1,Dis2,Dis3,Dis4,TDis,Wt1,Wt2,Wt3,Wt4

    RUV = 1

    !determine which data to interpolate from
    SELECT CASE(var)
      CASE("depth")
        v1 = depth(rnode1)
        v2 = depth(rnode2)
        v3 = depth(rnode3)
        v4 = depth(rnode4)
      CASE("angle")
        v1 = rho_angle(rnode1)
        v2 = rho_angle(rnode2)
        v3 = rho_angle(rnode3)
        v4 = rho_angle(rnode4)
      CASE("zetab")
        v1 = t_zeta(t_b,rnode1)
        v2 = t_zeta(t_b,rnode2)
        v3 = t_zeta(t_b,rnode3)
        v4 = t_zeta(t_b,rnode4)
      CASE("zetac")
        v1 = t_zeta(t_c,rnode1)
        v2 = t_zeta(t_c,rnode2)
        v3 = t_zeta(t_c,rnode3)
        v4 = t_zeta(t_c,rnode4)
      CASE("zetaf")
        v1 = t_zeta(t_f,rnode1)
        v2 = t_zeta(t_f,rnode2)
        v3 = t_zeta(t_f,rnode3)
        v4 = t_zeta(t_f,rnode4)
      CASE("saltb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_salt(t_b,rnode1,i)
        v2 = t_salt(t_b,rnode2,i)
        v3 = t_salt(t_b,rnode3,i)
        v4 = t_salt(t_b,rnode4,i)
      CASE("saltc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_salt(t_c,rnode1,i)
        v2 = t_salt(t_c,rnode2,i)
        v3 = t_salt(t_c,rnode3,i)
        v4 = t_salt(t_c,rnode4,i)
      CASE("saltf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_salt(t_f,rnode1,i)
        v2 = t_salt(t_f,rnode2,i)
        v3 = t_salt(t_f,rnode3,i)
        v4 = t_salt(t_f,rnode4,i)
      CASE("tempb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_temp(t_b,rnode1,i)
        v2 = t_temp(t_b,rnode2,i)
        v3 = t_temp(t_b,rnode3,i)
        v4 = t_temp(t_b,rnode4,i)
      CASE("tempc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_temp(t_c,rnode1,i)
        v2 = t_temp(t_c,rnode2,i)
        v3 = t_temp(t_c,rnode3,i)
        v4 = t_temp(t_c,rnode4,i)
      CASE("tempf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_temp(t_f,rnode1,i)
        v2 = t_temp(t_f,rnode2,i)
        v3 = t_temp(t_f,rnode3,i)
        v4 = t_temp(t_f,rnode4,i)
      CASE("denb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_den(t_b,rnode1,i)
        v2 = t_den(t_b,rnode2,i)
        v3 = t_den(t_b,rnode3,i)
        v4 = t_den(t_b,rnode4,i)
      CASE("denc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_den(t_c,rnode1,i)
        v2 = t_den(t_c,rnode2,i)
        v3 = t_den(t_c,rnode3,i)
        v4 = t_den(t_c,rnode4,i)
      CASE("denf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_den(t_f,rnode1,i)
        v2 = t_den(t_f,rnode2,i)
        v3 = t_den(t_f,rnode3,i)
        v4 = t_den(t_f,rnode4,i)
      CASE("uvelb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Uvel(t_b,unode1,i)
        v2 = t_Uvel(t_b,unode2,i)
        v3 = t_Uvel(t_b,unode3,i)
        v4 = t_Uvel(t_b,unode4,i)
        RUV = 2
      CASE("uvelc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Uvel(t_c,unode1,i)
        v2 = t_Uvel(t_c,unode2,i)
        v3 = t_Uvel(t_c,unode3,i)
        v4 = t_Uvel(t_c,unode4,i)
        RUV = 2
      CASE("uvelf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Uvel(t_f,unode1,i)
        v2 = t_Uvel(t_f,unode2,i)
        v3 = t_Uvel(t_f,unode3,i)
        v4 = t_Uvel(t_f,unode4,i)
        RUV = 2
      CASE("vvelb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Vvel(t_b,vnode1,i)
        v2 = t_Vvel(t_b,vnode2,i)
        v3 = t_Vvel(t_b,vnode3,i)
        v4 = t_Vvel(t_b,vnode4,i)
        RUV = 3
      CASE("vvelc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Vvel(t_c,vnode1,i)
        v2 = t_Vvel(t_c,vnode2,i)
        v3 = t_Vvel(t_c,vnode3,i)
        v4 = t_Vvel(t_c,vnode4,i)
        RUV = 3
      CASE("vvelf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Vvel(t_f,vnode1,i)
        v2 = t_Vvel(t_f,vnode2,i)
        v3 = t_Vvel(t_f,vnode3,i)
        v4 = t_Vvel(t_f,vnode4,i)
        RUV = 3
      CASE("wvelb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Wvel(t_b,rnode1,i)
        v2 = t_Wvel(t_b,rnode2,i)
        v3 = t_Wvel(t_b,rnode3,i)
        v4 = t_Wvel(t_b,rnode4,i)
      CASE("wvelc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Wvel(t_c,rnode1,i)
        v2 = t_Wvel(t_c,rnode2,i)
        v3 = t_Wvel(t_c,rnode3,i)
        v4 = t_Wvel(t_c,rnode4,i)
      CASE("wvelf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_Wvel(t_f,rnode1,i)
        v2 = t_Wvel(t_f,rnode2,i)
        v3 = t_Wvel(t_f,rnode3,i)
        v4 = t_Wvel(t_f,rnode4,i)
      CASE("khb")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_KH(t_b,rnode1,i)
        v2 = t_KH(t_b,rnode2,i)
        v3 = t_KH(t_b,rnode3,i)
        v4 = t_KH(t_b,rnode4,i)
      CASE("khc")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_KH(t_c,rnode1,i)
        v2 = t_KH(t_c,rnode2,i)
        v3 = t_KH(t_c,rnode3,i)
        v4 = t_KH(t_c,rnode4,i)
      CASE("khf")
        if(.not. present(i))then
          write(*,*) 'Problem interpolating ',var
          write(*,*) 'Optional Argument (i) Required for this Variable'
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop
        endif
        v1 = t_KH(t_f,rnode1,i)
        v2 = t_KH(t_f,rnode2,i)
        v3 = t_KH(t_f,rnode3,i)
        v4 = t_KH(t_f,rnode4,i)
      CASE DEFAULT
        write(*,*) 'Problem interpolating ',var
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop
    END SELECT

    !determine which node locations to interpolate from
    SELECT CASE(RUV)
      CASE(1)
        x1 = rx(rnode1)
        x2 = rx(rnode2)
        x3 = rx(rnode3)
        x4 = rx(rnode4)
        y1 = ry(rnode1)
        y2 = ry(rnode2)
        y3 = ry(rnode3)
        y4 = ry(rnode4)
         !Ensure there is no friction near land (the free slip condition) ewn.v.2
          IF (FreeSlip) THEN 
            rmasksum = 0
            ! determine if a land element is in present
            rmasksum = rho_mask(rnode1)+ rho_mask(rnode2)+ rho_mask(rnode3)    &
                     + rho_mask(rnode4)
            ! rmasksum – integer - the sum of the mask values of each node
            if (rmasksum .LT. 4) then  
              ! determine how many land nodes are present and reassign values 
              if (rmasksum .EQ. 3) then   !one land node
                 if (rho_mask(rnode1) .EQ. 0) then
                    v1 = 0.5*(v2+v4)
                 else if (rho_mask(rnode2) .EQ. 0) then 
                    v2 = 0.5*(v1+v3)
                 else if (rho_mask(rnode3) .EQ. 0) then 
                    v3 = 0.5*(v2+v4)
                 else if (rho_mask(rnode4) .EQ. 0) then 
                    v4 = 0.5*(v1+v3)
                 end if
              else if (rmasksum .EQ. 2) then   !two land nodes
                 if (rho_mask(rnode1).EQ.0 .AND. rho_mask(rnode2).EQ.0) then
                     v1 = v4 
                     v2 = v3
                 else if (rho_mask(rnode2).EQ.0 .AND. rho_mask(rnode3).EQ.0) then
                     v2 = v1 
                     v3 = v4
                 else if (rho_mask(rnode3).EQ.0 .AND. rho_mask(rnode4).EQ.0) then
                     v3 = v2 
                     v4 = v1
                 else if (rho_mask(rnode4).EQ.0 .AND. rho_mask(rnode1).EQ.0) then
                     v4 = v3 
                     v1 = v2
                 else if (rho_mask(rnode1).EQ.0 .AND. rho_mask(rnode3).EQ.0) then
                 v1 = v4
                 v3 = v2
                 else if (rho_mask(rnode4).EQ.0 .AND. rho_mask(rnode2).EQ.0) then
                 v4 = v1
                 v2 = v3
                 end if
              else if (rmasksum .EQ. 1) then   !three land nodes
                 if (rho_mask(rnode1) .EQ. 1) then
                      v2 = v1
                      v3 = v1  
                      v4 = v1
                 else if (rho_mask(rnode2) .EQ. 1) then
                      v1 = v2
                      v3 = v2  
                      v4 = v2
                 else if (rho_mask(rnode3) .EQ. 1) then
                      v1 = v3
                      v2 = v3  
                      v4 = v3
                 else if (rho_mask(rnode4) .EQ. 1) then
                      v1 = v4
                      v2 = v4  
                      v3 = v4
                 end if
              end if
            end if
          END IF 

      CASE(2)
        x1 = ux(unode1)
        x2 = ux(unode2)
        x3 = ux(unode3)
        x4 = ux(unode4)
        y1 = uy(unode1)
        y2 = uy(unode2)
        y3 = uy(unode3)
        y4 = uy(unode4)
         !Ensure there is no friction near land (the free slip condition) ewn.v.2
          IF (FreeSlip) THEN 
            ! determine if a land element is in present
            umasksum = 0
            umasksum = u_mask(unode1)+ u_mask(unode2)+ u_mask(unode3)          &
                     + u_mask(unode4)
            ! rmasksum – integer - the sum of the mask values of each node
            if (umasksum .LT. 4) then  
              ! determine how many land nodes are present and reassign values 
              if (umasksum .EQ. 1) then   !one land node
                 if (u_mask(unode1) .EQ. 0) then
                    v1 = 0.5*(v2+v4)
                 else if (u_mask(unode2) .EQ. 0) then 
                    v2 = 0.5*(v1+v3)
                 else if (u_mask(unode3) .EQ. 0) then 
                    v3 = 0.5*(v2+v4)
                 else if (u_mask(unode4) .EQ. 0) then 
                    v4 = 0.5*(v1+v3)
                 end if
              else if (umasksum .EQ. 2) then   !two land nodes
                 if (u_mask(unode1).EQ.0 .AND. u_mask(unode2).EQ.0) then
                     v1 = v4 
                     v2 = v3
                 else if (u_mask(unode2).EQ.0 .AND. u_mask(unode3).EQ.0) then
                     v2 = v1 
                     v3 = v4
                 else if (u_mask(unode3).EQ.0 .AND. u_mask(unode4).EQ.0) then
                     v3 = v2 
                     v4 = v1
                 else if (u_mask(unode4).EQ.0 .AND. u_mask(unode1).EQ.0) then
                     v4 = v3 
                     v1 = v2
                 else if (u_mask(unode1).EQ.0 .AND. u_mask(unode3).EQ.0) then
                 v1 = v4
                 v3 = v2
                 else if (u_mask(unode4).EQ.0 .AND. u_mask(unode2).EQ.0) then
                 v4 = v1
                 v2 = v3
                 end if
              else if (umasksum .EQ. 1) then   !three land nodes
                 if (u_mask(unode1) .EQ. 1) then
                      v2 = v1
                      v3 = v1  
                      v4 = v1
                 else if (u_mask(unode2) .EQ. 1) then
                      v1 = v2
                      v3 = v2  
                      v4 = v2
                 else if (u_mask(unode3) .EQ. 1) then
                      v1 = v3
                      v2 = v3  
                      v4 = v3
                 else if (u_mask(unode4) .EQ. 1) then
                      v1 = v4
                      v2 = v4  
                      v3 = v4
                 end if
              end if
            end if
          END IF 
      CASE(3)
        x1 = vx(vnode1)
        x2 = vx(vnode2)
        x3 = vx(vnode3)
        x4 = vx(vnode4)
        y1 = vy(vnode1)
        y2 = vy(vnode2)
        y3 = vy(vnode3)
        y4 = vy(vnode4)
         !Ensure there is no friction near land (the free slip condition) ewn.v.2
          IF (FreeSlip) THEN 
            ! determine if a land element is in present
            vmasksum = 0
            vmasksum = v_mask(vnode1)+ v_mask(vnode2)+ v_mask(vnode3)           &                   
                     + v_mask(vnode4)
            ! rmasksum – integer - the sum of the mask values of each node
            if (vmasksum .LT. 4) then  
              ! determine how many land nodes are present and reassign values 
              if (vmasksum .EQ. 3) then   !one land node
                 if (v_mask(vnode1) .EQ. 0) then
                    v1 = 0.5*(v2+v4)
                 else if (v_mask(vnode2) .EQ. 0) then 
                    v2 = 0.5*(v1+v3)
                 else if (v_mask(vnode3) .EQ. 0) then 
                    v3 = 0.5*(v2+v4)
                 else if (v_mask(vnode4) .EQ. 0) then 
                    v4 = 0.5*(v1+v3)
                 end if
              else if (vmasksum .EQ. 2) then   !two land nodes
                 if (v_mask(vnode1).EQ.0 .AND. v_mask(vnode2).EQ.0) then
                     v1 = v4 
                     v2 = v3
                 else if (v_mask(vnode2).EQ.0 .AND. v_mask(vnode3).EQ.0) then
                     v2 = v1 
                     v3 = v4
                 else if (v_mask(vnode3).EQ.0 .AND. v_mask(vnode4).EQ.0) then
                     v3 = v2 
                     v4 = v1
                 else if (v_mask(vnode4).EQ.0 .AND. v_mask(vnode1).EQ.0) then
                     v4 = v3 
                     v1 = v2
                 else if (v_mask(unode1).EQ.0 .AND. v_mask(unode3).EQ.0) then
                 v1 = v4
                 v3 = v2
                 else if (v_mask(unode4).EQ.0 .AND. v_mask(unode2).EQ.0) then
                 v4 = v1
                 v2 = v3
                 end if
              else if (vmasksum .EQ. 1) then   !three land nodes
                 if (v_mask(vnode1) .EQ. 1) then
                      v2 = v1
                      v3 = v1  
                      v4 = v1
                 else if (v_mask(vnode2) .EQ. 1) then
                      v1 = v2
                      v3 = v2  
                      v4 = v2
                 else if (v_mask(vnode3) .EQ. 1) then
                      v1 = v3
                      v2 = v3  
                      v4 = v3
                 else if (v_mask(vnode4) .EQ. 1) then
                      v1 = v4
                      v2 = v4  
                      v3 = v4
                 end if
              end if
            end if
          END IF 
    END SELECT

    !determine the method of interpolation:

    ! bilinear interpolation of first triangle
    tt = ((xp-x1)*(y3-y1)+(y1-yp)*(x3-x1)) / ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
    uu = ((xp-x1)*(y2-y1)+(y1-yp)*(x2-x1)) / ((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))
    vp = v1 + (v2-v1)*tt + (v3-v1)*uu

    ! if outside triangle, then do bilinear interpolation of other triangle
    if( tt.LT.0. .OR. uu.LT.0. .OR. (tt+uu).GT.1.0 ) then
      tt = ((xp-x3)*(y1-y3)+(y3-yp)*(x1-x3))/((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))
      uu = ((xp-x3)*(y4-y3)+(y3-yp)*(x4-x3))/((x1-x3)*(y4-y3)-(y1-y3)*(x4-x3))
      vp = v3 + (v4-v3)*tt + (v1-v3)*uu

      !if bilinear techniques are undefined, then use inverse weighted distance
      if( tt.LT.0. .OR. uu.LT.0. .OR. (tt+uu).GT.1.0 ) then
        !if particle on node, then set equal to node value
        if(  (xp.EQ.x1 .AND. yp.EQ.y1).OR.(xp.EQ.x2 .AND. yp.EQ.y2)            &
         .OR.(xp.EQ.x3 .AND. yp.EQ.y3).OR.(xp.EQ.x4 .AND. yp.EQ.y4))then
          if (xp.EQ.x1 .AND. yp.EQ.y1) vp=v1
          if (xp.EQ.x2 .AND. yp.EQ.y2) vp=v2
          if (xp.EQ.x3 .AND. yp.EQ.y3) vp=v3
          if (xp.EQ.x4 .AND. yp.EQ.y4) vp=v4
        else !use inverse weighted distance instead  
          Dis1=1./( SQRT( (x1-xp)**2 + (y1-yp)**2 ) ) 
          Dis2=1./( SQRT( (x2-xp)**2 + (y2-yp)**2 ) ) 
          Dis3=1./( SQRT( (x3-xp)**2 + (y3-yp)**2 ) ) 
          Dis4=1./( SQRT( (x4-xp)**2 + (y4-yp)**2 ) ) 
          TDis = Dis1+Dis2+Dis3+Dis4
          Wt1= Dis1/TDis
          Wt2= Dis2/TDis
          Wt3= Dis3/TDis
          Wt4= Dis4/TDis
          vp = Wt1*v1 + Wt2*v2 + Wt3*v3 + Wt4*v4 
        endif
      endif
    endif

    interp = vp

  END FUNCTION interp

  !This function creates a Water Column Tension Spline at back, center, and 
  !  forward hydrodynamic time then uses Polynomial Interpolation to determine
  !  Internal Time values to finally get the value of the particle in space and
  !  time.  The name is derived from: Water Column Tension Spline, Internal 
  !  Time Polynomial Interpolation.  The final variable v is for version
  !  (ie what is to be returned): 1-back, 2-center, 3-forward, 4-(b+4c+f)/6
  DOUBLE PRECISION FUNCTION WCTS_ITPI(var,Xpos,Ypos,deplvl,Pwc_zb,Pwc_zc,      &
                                      Pwc_zf,slvls,P_zb,P_zc,P_zf,ex,ix,p,v)
    USE TENSION_MOD, ONLY: TSPSI,HVAL
    USE INT_MOD, ONLY: linint,polintd
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: deplvl,slvls,p,v
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,Pwc_zb(slvls),Pwc_zc(slvls),     &
                                    Pwc_zf(slvls),P_zb,P_zc,P_zf,ex(3),ix(3)
    CHARACTER(LEN=*), INTENT(IN) :: var

    INTEGER, PARAMETER :: nN = 4

    CHARACTER(LEN=LEN(var)+1) :: varb,varc,varf
    INTEGER :: i
    DOUBLE PRECISION :: abb_zb(nN),abb_zc(nN),abb_zf(nN),abb_vb(nN),           &
      abb_vc(nN),abb_vf(nN),P_vb,P_vc,P_vf,P_V,ey(3),vb,vc,vf,slope

    !TSPACK Variables
    INTEGER :: IER,SigErr
    DOUBLE PRECISION :: YP(nN),SIGM(nN)
      
    varb = var//"b"
    varc = var//"c"
    varf = var//"f"

    do i=1,nN
      abb_zb(i) = Pwc_zb(i+deplvl-1)
      abb_zc(i) = Pwc_zc(i+deplvl-1)
      abb_zf(i) = Pwc_zf(i+deplvl-1)
      abb_vb(i) = interp(Xpos,Ypos,varb,i+deplvl-1)
      abb_vc(i) = interp(Xpos,Ypos,varc,i+deplvl-1)
      abb_vf(i) = interp(Xpos,Ypos,varf,i+deplvl-1)
    enddo

!       *********************************************************
!       *       5Aiic6b.  Fit Tension Spline to WC Profile      *
!       *********************************************************

    !ii. call TSPACK to fit a tension spline to water column profile 
    !  of U,V,W velocities at particle x-y location and find value at particle

    P_vb=0.0
    SigErr=0
    CALL TSPSI (nN,abb_zb,abb_vb,YP,SIGM,IER,SigErr)
    IF (SigErr.EQ.0) THEN
      P_vb = HVAL (P_zb,nN,abb_zb,abb_vb,YP,SIGM,IER)
    ELSE
      CALL linint(abb_zb,abb_vb,nN,P_zb,P_vb,slope)
    ENDIF

    P_vc=0.0
    SigErr=0
    CALL TSPSI (nN,abb_zc,abb_vc,YP,SIGM,IER,SigErr)
    IF (SigErr.EQ.0) THEN
      P_vc = HVAL (P_zc,nN,abb_zc,abb_vc,YP,SIGM,IER)
    ELSE
      CALL linint(abb_zc,abb_vc,nN,P_zc,P_vc,slope)
    ENDIF

    P_vf=0.0
    SigErr=0
    CALL TSPSI (nN,abb_zf,abb_vf,YP,SIGM,IER,SigErr)
    IF (SigErr.EQ.0) THEN
      P_vf = HVAL (P_zf,nN,abb_zf,abb_vf,YP,SIGM,IER)
    ELSE
      CALL linint(abb_zf,abb_vf,nN,P_zf,P_vf,slope)
    ENDIF


!       *********************************************************
!       *               Find Internal b,c,f Values              *
!       *********************************************************

    !iii. fit polynomial to hydrodynamic model output and find 
    !  internal b,c,f values

    !    1. Prepare external time step values
    if (p .EQ. 1) then
      ey=0.0
      ey(1) = P_vb
      ey(2) = P_vb
      ey(3) = P_vc
    else
      ey=0.0
      ey(1) = P_vb
      ey(2) = P_vc
      ey(3) = P_vf
    endif

    !    2. Get value
    vb = polintd(ex,ey,3,ix(1))
    vc = polintd(ex,ey,3,ix(2))
    vf = polintd(ex,ey,3,ix(3))
    P_V = (vb + vc*4 + vf) / DBLE(6.0)

    SELECT CASE (v)
      CASE (1)
        WCTS_ITPI = vb
      CASE (2)
        WCTS_ITPI = vc
      CASE (3)
        WCTS_ITPI = vf
      CASE (4)
        WCTS_ITPI = P_V
      CASE DEFAULT
        write(*,*) 'ERROR: Illegal WCTS_ITPI version number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop
    END SELECT

  END FUNCTION WCTS_ITPI

  DOUBLE PRECISION FUNCTION getSlevel(zeta,depth,i)
    !This function returns the depth of the current s-level
    USE PARAM_MOD, ONLY: hc,Vtransform
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: zeta,depth

    DOUBLE PRECISION :: S,h

    ! convert negative depth to positive depth
    h = DBLE(-1.0) * depth


    SELECT CASE(Vtransform)

      CASE(1)  !Rutgers-ROMS formulation, eqn (1) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = hc*SC(i)+(h-hc)*CS(i)
        getSlevel = S+zeta*(DBLE(1.0)+S/h)

      CASE(2)  !UCLA-formulation, eqn(2) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = (hc*SC(i)+h*CS(i)) / (hc+h)
        getSlevel = zeta+(zeta+h)*S

      CASE(3)  !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit
        !ocean circulation model using a generalized topography-following 
        !coordinate system, J. Comp. Phys., 115 (1), 228-244.

        getSlevel = zeta*(DBLE(1.0)+SC(i))+hc*SC(i)+(h-hc)*CS(i)

      CASE DEFAULT
        write(*,*) 'ERROR: Illegal Vtransform number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop

    END SELECT

  END FUNCTION getSlevel



  DOUBLE PRECISION FUNCTION getWlevel(zeta,depth,i)
    !This function returns the depth of the current w s-level
    USE PARAM_MOD, ONLY: hc,Vtransform
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: zeta,depth

    DOUBLE PRECISION :: S,h

    ! convert negative depth to positive depth
    h = DBLE(-1.0) * depth


    SELECT CASE(Vtransform)

      CASE(1)  !Rutgers-ROMS formulation, eqn (1) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = hc*SCW(i)+(h-hc)*CSW(i)
        getWlevel = S+zeta*(DBLE(1.0)+S/h)

      CASE(2)  !UCLA-formulation, eqn(2) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = (hc*SCW(i)+h*CSW(i))/(hc+h)
        getWlevel = zeta+(zeta+h)*S

      CASE(3)  !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit
        !ocean circulation model using a generalized topography-following 
        !coordinate system, J. Comp. Phys., 115 (1), 228-244.

        getWlevel = zeta*(DBLE(1.0)+SCW(i))+hc*SCW(i)+(h-hc)*CSW(i)

      CASE DEFAULT
        write(*,*) 'ERROR: Illegal Vtransform number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop

    END SELECT

  END FUNCTION getWlevel


  SUBROUTINE getMask_Rho(mask)
    !This subroutine returns the values in the variable mask_rho
    !This is used by createBounds() in the boundary module to make the 
    !  boundaries based on mask_rho
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: mask(:,:)

    if(GRD_SET)then
      mask = mask_rho
    else
      write(*,*) 'ERROR: Cannot create boundaries, mask_rho not yet read in'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      stop
    endif

  END SUBROUTINE getMask_Rho

  SUBROUTINE getUVxy(ux,uy,vx,vy)
    !This subroutine returns the values in the variables x_u,y_u,x_v,y_v
    !This is used by createBounds() in the boundary module to make the 
    !  boundaries on the U & V node locations
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: ux(:,:),uy(:,:),vx(:,:),vy(:,:)

    if(GRD_SET)then
      ux = x_u
      uy = y_u
      vx = x_v
      vy = y_v
    else
      write(*,*) 'ERROR: Cannot create boundaries, x_u, y_u, x_v, or y_v is ', &
                 'not yet read in'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      stop
    endif

  END SUBROUTINE getUVxy

  SUBROUTINE getR_ele(ele_x,ele_y)
    !This subroutine returns the values in the variables r_ele_x, and r_ele_y
    !This is used by createPolySpecs() in the settlement module to determine 
    !  which habitat polygons are in each element
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: ele_x(:,:),ele_y(:,:)

    if(GRD_SET)then
      ele_x = r_ele_x
      ele_y = r_ele_y
    else
      write(*,*) 'ERROR: Cannot create Poly Specs, r_ele_x or r_ele_y is not', &
                 ' yet created'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      stop
    endif

  END SUBROUTINE getR_ele

  INTEGER FUNCTION getP_r_element(n)
    !This subroutine returns the id of the rho element the particle is 
    !  currently in
    !This is used by settlement() to determine which habitat polygons to check
    !  for settlement
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    getP_r_element = P_r_element(n)
  END FUNCTION getP_r_element


  SUBROUTINE setijruv()
    USE PARAM_MOD, ONLY: numpar,vi,ui,vj,uj,ijbuff
    IMPLICIT NONE

    INTEGER :: i,j,n

    !t_ijruv  1 - rho_imin
    !         2 - rho_imax
    !         3 - rho_jmin
    !         4 - rho_jmax
    !         5 - u_imin
    !         6 - u_imax
    !         7 - u_jmin
    !         8 - u_jmax
    !         9 - v_imin
    !        10 - v_imax
    !        11 - v_jmin
    !        12 - v_jmax

    !rho
    i = mod(RE(1,P_r_element(1))-1,vi)+1
    j = (RE(1,P_r_element(1))-1)/vi + 1
    t_ijruv(1) = max(i-ijbuff,1)
    t_ijruv(2) = min(i+ijbuff+1,vi)
    t_ijruv(3) = max(j-ijbuff,1)
    t_ijruv(4) = min(j+ijbuff+1,uj)

    do n=2,numpar
      i = mod(RE(1,P_r_element(n))-1,vi)+1
      j = (RE(1,P_r_element(n))-1)/vi + 1
      if((i-ijbuff  ) < t_ijruv(1)) t_ijruv(1) = max(i-ijbuff,1)
      if((i+ijbuff+1) > t_ijruv(2)) t_ijruv(2) = min(i+ijbuff+1,vi)
      if((j-ijbuff  ) < t_ijruv(3)) t_ijruv(3) = max(j-ijbuff,1)
      if((j+ijbuff+1) > t_ijruv(4)) t_ijruv(4) = min(j+ijbuff+1,uj)
    enddo

    !u
    i = mod(UE(1,P_u_element(1))-1,ui)+1
    j = (UE(1,P_u_element(1))-1)/ui + 1
    t_ijruv(5) = max(i-ijbuff,1)
    t_ijruv(6) = min(i+ijbuff+1,ui)
    t_ijruv(7) = max(j-ijbuff,1)
    t_ijruv(8) = min(j+ijbuff+1,uj)

    do n=2,numpar
      i = mod(UE(1,P_u_element(n))-1,ui)+1
      j = (UE(1,P_u_element(n))-1)/ui + 1
      if((i-ijbuff  ) < t_ijruv(5)) t_ijruv(5) = max(i-ijbuff,1)
      if((i+ijbuff+1) > t_ijruv(6)) t_ijruv(6) = min(i+ijbuff+1,ui)
      if((j-ijbuff  ) < t_ijruv(7)) t_ijruv(7) = max(j-ijbuff,1)
      if((j+ijbuff+1) > t_ijruv(8)) t_ijruv(8) = min(j+ijbuff+1,uj)
    enddo

    !v
    i = mod(VE(1,P_v_element(1))-1,vi)+1
    j = (VE(1,P_v_element(1))-1)/vi + 1
    t_ijruv( 9) = max(i-ijbuff,1)
    t_ijruv(10) = min(i+ijbuff+1,vi)
    t_ijruv(11) = max(j-ijbuff,1)
    t_ijruv(12) = min(j+ijbuff+1,vj)

    do n=2,numpar
      i = mod(VE(1,P_v_element(n))-1,vi)+1
      j = (VE(1,P_v_element(n))-1)/vi + 1
      if((i-ijbuff  ) < t_ijruv( 9)) t_ijruv( 9) = max(i-ijbuff,1)
      if((i+ijbuff+1) > t_ijruv(10)) t_ijruv(10) = min(i+ijbuff+1,vi)
      if((j-ijbuff  ) < t_ijruv(11)) t_ijruv(11) = max(j-ijbuff,1)
      if((j+ijbuff+1) > t_ijruv(12)) t_ijruv(12) = min(j+ijbuff+1,vj)
    enddo


  END SUBROUTINE setijruv


  SUBROUTINE finHydro()
    !This subroutine closes all the module's allocatable variables
    IMPLICIT NONE

    !ALLOCATE MODULE VARIABLES
    DEALLOCATE(r_Adjacent,u_Adjacent,v_Adjacent)
    DEALLOCATE(rho_angle,mask_rho,depth)
    DEALLOCATE(SC,CS,SCW,CSW,RE,UE,VE)
    DEALLOCATE(rx,ry,ux,uy,vx,vy)
    DEALLOCATE(r_ele_x,r_ele_y)
    DEALLOCATE(u_ele_x,u_ele_y)
    DEALLOCATE(v_ele_x,v_ele_y)
    DEALLOCATE(x_u,y_u,x_v,y_v)

    DEALLOCATE(P_r_element,P_u_element,P_v_element)
    DEALLOCATE(t_zeta,t_salt,t_temp,t_Wvel,t_Uvel,t_Vvel,t_Kh)

  END SUBROUTINE finHydro


  SUBROUTINE initNetCDF()
  
    !Initialize NetCDF Counters
    NCcount = 0
    NCstart = 0

  END SUBROUTINE initNetCDF

  SUBROUTINE createNetCDF(dob)
    USE PARAM_MOD, ONLY: numpar,NCOutFile,outpath,outpathGiven,NCtime,         &
        SVN_Version,RunName,ExeDir,OutDir,RunBy,Institution,StartedOn,         &
        TrackCollisions,SaltTempOn
    USE netcdf
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(numpar), OPTIONAL, INTENT(IN) :: dob

    INCLUDE 'netcdf.inc'

    CHARACTER(LEN=200) :: ncFile
    INTEGER :: STATUS,NCID,numparID,timeID,pageID,modtimeID,lonID,latID,       &
               depthID,colorID,hitBID,hitLID,dobID,saltID,tempID

    !   NF90_CREATE           ! create netCDF dataset: enter define mode
    !        ...
    !      NF90_DEF_DIM       ! define dimensions: from name and length
    !        ...
    !      NF90_DEF_VAR       ! define variables: from name, type, dims
    !        ...
    !      NF90_PUT_ATT       ! assign attribute values
    !        ...
    !   NF90_ENDDEF           ! end definitions: leave define mode
    !        ...
    !      NF90_PUT_VAR       ! provide values for variable
    !        ...
    !   NF90_CLOSE            ! close: save new netCDF dataset


    !NF90_CREATE

    !Reset Print Counter to 0
    prcount = 0

    IF(outpathGiven)THEN
      IF(NCtime == 0 ) THEN
        ncFile = TRIM(outpath) // TRIM(NCOutFile) // '.nc'
      ELSE
        NCcount = NCcount + 1
        write(ncFile,"(A,A,A,I3.3,A)")TRIM(outpath),TRIM(NCOutFile),'_',       &
                                      NCcount,'.nc'
      ENDIF
    ELSE
      IF(NCtime == 0 ) THEN
        ncFile = TRIM(NCOutFile) // '.nc'
      ELSE
        NCcount = NCcount + 1
        write(ncFile,"(A,A,I3.3,A)")TRIM(NCOutFile),'_',NCcount,'.nc'
      ENDIF
    ENDIF

    write(*,*)'Creating NetCDF Output File: ',TRIM(ncFile)

    STATUS = NF90_CREATE(TRIM(ncFile), NF90_CLOBBER, NCID)
    IF(STATUS /= NF90_NOERR) THEN
      WRITE(*,*) 'Problem creating NetCDF output file'
      WRITE(*,*) NF_STRERROR(STATUS)
      STOP
    ENDIF

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_DIM

        STATUS = NF90_DEF_DIM(NCID,'numpar',numpar,numparID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: numpar dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_DIM(NCID,'time',NF90_UNLIMITED,timeID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: time dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_VAR

        STATUS = NF90_DEF_VAR(NCID,'model_time',NF_DOUBLE,(/timeID/),modtimeID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: time var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        IF( PRESENT(dob) )THEN
          STATUS = NF90_DEF_VAR(NCID,'dob',NF_DOUBLE,(/numparID/),dobID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: dob var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        STATUS = NF90_DEF_VAR(NCID,'age',NF_DOUBLE,(/numparID,timeID/),pageID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: age var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'lon',NF_DOUBLE,(/numparID,timeID/),lonID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: lon var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'lat',NF_DOUBLE,(/numparID,timeID/),latID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: lat var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'depth',NF_DOUBLE,(/numparID,timeID/),      &
                              depthID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: depth var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'color',NF_DOUBLE,(/numparID,timeID/),      &
                              colorID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: color var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        IF(TrackCollisions)THEN
          STATUS =NF90_DEF_VAR(NCID,'hitBottom',NF_DOUBLE,(/numparID,timeID/), &
                               hitBID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Bottom ', &
                                              'Collision var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_DEF_VAR(NCID,'hitLand',NF_DOUBLE,(/numparID,timeID/),  &
                                hitLID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Land ',   &
                                              'Collision var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        IF(SaltTempOn)THEN
          STATUS = NF90_DEF_VAR(NCID,'salinity',NF_DOUBLE,(/numparID,timeID/), &
                                saltID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'Salinity var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_DEF_VAR(NCID,'temperature',NF_DOUBLE,                  &
                                (/numparID,timeID/),tempID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'Temperature var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_PUT_ATT

        !Particle Time
        STATUS = NF90_PUT_ATT(NCID, modtimeID, "long_name",                    &
                              "time that has passed thus far in the model")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, modtimeID, "units", "seconds")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, modtimeID, "field",                        &
                              "model_time, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        IF( PRESENT(dob) )THEN
          !Particle Date of Birth
          STATUS = NF90_PUT_ATT(NCID, dobID, "long_name",                      &
                   "Date of Birth of particles in seconds from model start")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, dobID, "units", "seconds")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, dobID, "field", "age, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        !Particle Age
        STATUS = NF90_PUT_ATT(NCID, pageID, "long_name", "age of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, pageID, "units", "seconds")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, pageID, "field", "age, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        !Longitude
        STATUS = NF90_PUT_ATT(NCID, lonID, "long_name",                        &
                              "longitude of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, lonID, "units", "decimal degrees E")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, lonID, "field", "lon, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        !Latitude
        STATUS =NF90_PUT_ATT(NCID, latID, "long_name", "latitude of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, latID, "units", "decimal degrees N")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, latID, "field", "lat, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        !Depth
        STATUS = NF90_PUT_ATT(NCID, depthID, "long_name", "depth of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, depthID, "units", "meters below surface")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, depthID, "field", "depth, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        !Color
        STATUS = NF90_PUT_ATT(NCID, colorID, "long_name",                      &
                 "identification number for particle behavior or status")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, colorID, "units",                          &
                              "nondimensional, see LTRANS User Guide")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, colorID, "field", "color, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        IF(TrackCollisions)THEN
          !hitBottom
          STATUS = NF90_PUT_ATT(NCID, hitBID, "long_name",                     &
                                "# of times Particle Collided with Bottom")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitBID, "units", "Number of Collisions")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitBID, "field",                         &
                                "hitBottom, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


          !hitLand
          STATUS = NF90_PUT_ATT(NCID, hitLID, "long_name",                     &
                                "# of times Particle Collided with Land")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitLID, "units", "Number of Collisions")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitLID, "field",                         &
                                "hitLand, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        IF(SaltTempOn)THEN
          !salt
          STATUS = NF90_PUT_ATT(NCID, saltID, "long_name",                     &
                                "Salinity at the particle's location")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS=NF90_PUT_ATT(NCID,saltID, "field", "salinity, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


          !temp
          STATUS = NF90_PUT_ATT(NCID, tempID, "long_name",                     &
                                "Temperature at the particle's location")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, tempID, "units", "° Celcius")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, tempID, "field",                         &
                                "temperature, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF


        !Global
        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "type",                       &
                              "Position and characteristics of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "title", "LTRANS output")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "author", "Zachary Schlag")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "svn", SVN_Version)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "run_name", RunName)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "executable", ExeDir)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "output_loc", OutDir)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "administrator", RunBy)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "institution", Institution)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "date", StartedOn)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF_ENDDEF

      STATUS = NF90_ENDDEF(NCID)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: EndDef'
      IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF_PUT_VAR

      !Particle Date of Birth
      IF( PRESENT(dob) )THEN
        STATUS = NF90_INQ_VARID(NCID, "dob", dobID)
        STATUS = NF90_PUT_VAR(NCID, dobID, dob)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put dob'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)  
      ENDIF

    !NF_CLOSE

    STATUS = NF_CLOSE(NCID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Close'
    IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

  END SUBROUTINE createNetCDF


  SUBROUTINE writeNetCDF(time,age,lon,lat,depth,colors,hitB,hitL,Salt,Temp)
    USE PARAM_MOD, ONLY: numpar,SaltTempOn,NCOutFile,outpath,outpathGiven,     &
        NCtime,TrackCollisions
    USE netcdf
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: time
    DOUBLE PRECISION, INTENT(IN) :: age(numpar),lon(numpar),lat(numpar),       &
                                    depth(numpar),colors(numpar)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: Salt(numpar),Temp(numpar)
    INTEGER, INTENT(IN), OPTIONAL :: hitB(numpar),hitL(numpar)

    INCLUDE 'netcdf.inc'

    CHARACTER(LEN=200) :: ncFile
    INTEGER :: STATUS,NCID,modtimeID,pageID,lonID,latID,depthID,hitBID,hitLID, &
               colorID,saltID,tempID
    INTEGER :: NCelapsed

    !If only one NetCDF output file is being written to:
    IF(NCtime == 0) THEN

      IF(outpathGiven)THEN
        ncFile = TRIM(outpath) // TRIM(NCOutFile) // '.nc'
      ELSE
        ncFile = TRIM(NCOutFile) // '.nc'
      ENDIF

    !If sequentially numbered NetCDF output files are being written to:
    ELSE

      NCelapsed = time - NCstart

      !If specified time interval has been reached, create new NetCDF file
      IF(NCelapsed >= NCtime) THEN
        NCstart = time
        call createNetCDF()
      ENDIF

      IF(outpathGiven)THEN
        write(ncFile,"(A,A,A,I3.3,A)")TRIM(outpath),TRIM(NCOutFile),'_',       &
                                      NCcount,'.nc'
      ELSE
        write(ncFile,"(A,A,I3.3,A)")TRIM(NCOutFile),'_',NCcount,'.nc'
      ENDIF

    ENDIF

    prcount = prcount + 1

    STATUS = NF90_OPEN(TRIM(ncFile), NF90_WRITE, NCID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Particle Time
      STATUS = NF90_INQ_VARID(NCID, "model_time", modtimeID)
      STATUS = NF90_PUT_VAR(NCID, modtimeID, DBLE(time), start = (/ prcount /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put model_time, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)  

      !Particle Age
      STATUS = NF90_INQ_VARID(NCID, "age", pageID)
      STATUS = NF90_PUT_VAR(NCID, pageID, age,             &
                            start = (/ 1, prcount /),     &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put age, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Longitude
      STATUS = NF90_INQ_VARID(NCID, "lon", lonID)
      STATUS = NF90_PUT_VAR(NCID, lonID, lon,             &
                            start = (/ 1, prcount /),     &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put lon, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Latitude
      STATUS = NF90_INQ_VARID(NCID, "lat", latID)
      STATUS = NF90_PUT_VAR(NCID, latID, lat,             &
                            start = (/ 1, prcount /),     &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put lat, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Depth
      STATUS = NF90_INQ_VARID(NCID, "depth", depthID)
      STATUS = NF90_PUT_VAR(NCID, depthID, depth,         &
                            start = (/ 1, prcount /),     &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put depth, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Color
      STATUS = NF90_INQ_VARID(NCID, "color", colorID)
      STATUS = NF90_PUT_VAR(NCID, colorID, colors,        &
                            start = (/ 1, prcount /),     &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put color, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !hitBottom
      IF( PRESENT(hitB) )THEN
        STATUS = NF90_INQ_VARID(NCID, "hitBottom", hitBID)
        STATUS = NF90_PUT_VAR(NCID, hitBID, hitB,       &
                              start = (/ 1, prcount /),   &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*)                                    &
          'Problem put # Bottom Collisions, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      !hitLand
      IF( PRESENT(hitL) )THEN
        STATUS = NF90_INQ_VARID(NCID, "hitLand", hitLID)
        STATUS = NF90_PUT_VAR(NCID, hitLID, hitL,      &
                              start = (/ 1, prcount /),   &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*)                                    &
          'Problem put # Land Collisions, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      !hitBottom
      IF( PRESENT(salt) )THEN
        STATUS = NF90_INQ_VARID(NCID, "salinity", saltID)
        STATUS = NF90_PUT_VAR(NCID, saltID, salt,       &
                              start = (/ 1, prcount /),   &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put salinity, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      !hitLand
      IF( PRESENT(temp) )THEN
        STATUS = NF90_INQ_VARID(NCID, "temperature", tempID)
        STATUS = NF90_PUT_VAR(NCID, tempID, temp,      &
                              start = (/ 1, prcount /),   &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put temperature, time: ', &
                                            time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

    STATUS = NF_CLOSE(NCID)

  END SUBROUTINE writeNetCDF

END MODULE HYDRO_MOD

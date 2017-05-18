!
!This file defines many variables that are read from the GRID.data and LTRANS.data files
!And also groups the input parameters/variables into namelists
!
!NOTE: variables in a namelist can NOT be dynamic variables!!!!! 
!      dynamic namelists are NOT yet supported in FORTRAN90/95 standard
!
!Basically, this file replaces GRID.inc and LTRANS.inc in previous code
!

!--these are the grid dimension definitions, they should read from GRID.data

  INTEGER :: ui                     ! number u points in x direction
  INTEGER :: uj                     !        u           y 
  INTEGER :: vi                     !        v           x
  INTEGER :: vj                     !        v           y

  INTEGER :: rho_nodes              ! number rho nodes
  INTEGER :: u_nodes                ! number of u nodes
  INTEGER :: v_nodes                ! number of v nodes

  INTEGER :: max_rho_elements       ! total number of rho elements 
  INTEGER :: max_u_elements         ! total number of u elements
  INTEGER :: max_v_elements         ! total           v

  INTEGER :: rho_elements           ! number of rho elements with at least 1 vertex is water
  INTEGER :: u_elements             !            u 
  INTEGER :: v_elements             !            v

!group the grid info section in a namelist:

  namelist/gridinfo/ ui, uj, vi, vj,rho_nodes, u_nodes, v_nodes,     &
     &               max_rho_elements,max_u_elements,max_v_elements, &
     &               rho_elements, u_elements, v_elements         





!The following used to be in LTRANS.inc:


!*** NUMBER OF PARTICLES ***

  INTEGER :: numpar

  namelist/numparticles/numpar      ! Number of particles



!*** TIME PARAMETERS ***

  REAL             :: days          ! Number of days to run the model
  INTEGER          :: iprint        ! Print interval for LTRANS output (s); 3600 = every hour
  INTEGER          :: dt            ! External time step (duration between hydro model predictions) (s) 
  INTEGER          :: idt           ! Internal (particle tracking) time step (s)

  namelist/timeparam/ days,iprint,dt,idt



!!*** ROMS HYDRODYNAMIC MODULE PARAMETERS ***
  INTEGER          :: us            ! Number of Rho grid s-levels in ROMS hydro model
  INTEGER          :: ws            ! Number of W grid s-levels in ROMS hydro model
  INTEGER          :: tdim          ! Number of time steps per ROMS hydro predictions file
  REAL             :: hc            ! Min Depth - used in ROMS S-level transformations
  DOUBLE PRECISION :: z0            ! ROMS roughness parameter
  INTEGER          :: Vtransform    ! 1-WikiRoms Eq. 1, 2-WikiRoms Eq. 2, 3-Song/Haidvogel 1994 Eq.
  LOGICAL          :: readZeta      ! If .TRUE. read in sea-surface height   (zeta) from NetCDF file, else use constZeta
  DOUBLE PRECISION :: constZeta     ! Constant value for Zeta if readZeta is .FALSE.
  LOGICAL          :: readSalt      ! If .TRUE. read in salinity             (salt) from NetCDF file, else use constSalt
  DOUBLE PRECISION :: constSalt     ! Constant value for Salt if readSalt is .FALSE.
  LOGICAL          :: readTemp      ! If .TRUE. read in temperature          (temp) from NetCDF file, else use constTemp
  DOUBLE PRECISION :: constTemp     ! Constant value for Temp if readTemp is .FALSE.
  LOGICAL          :: readU         ! If .TRUE. read in u-momentum component (U   ) from NetCDF file, else use constU
  DOUBLE PRECISION :: constU        ! Constant value for U if readU is .FALSE.
  LOGICAL          :: readV         ! If .TRUE. read in v-momentum component (V   ) from NetCDF file, else use constV
  DOUBLE PRECISION :: constV        ! Constant value for V if readV is .FALSE.
  LOGICAL          :: readW         ! If .TRUE. read in w-momentum component (W   ) from NetCDF file, else use constW
  DOUBLE PRECISION :: constW        ! Constant value for W if readW is .FALSE.
  LOGICAL          :: readAks       ! If .TRUE. read in salinity vertical diffusion coefficient (Aks ) from NetCDF file, else use constAks
  DOUBLE PRECISION :: constAks      ! Constant value for Aks if readAks is .FALSE.
  LOGICAL          :: readDens
  DOUBLE PRECISION :: constDens

!
  namelist/hydroparam/us,ws,tdim,hc,z0,Vtransform,readZeta,constZeta,readSalt,   &
                    & constSalt,readTemp,constTemp,readU,readU,constU,readV,     &
                    & constV,readW,constW,readAks,constAks,readDens,constDens



!*** TURBULENCE MODULE PARAMETERS ***
  LOGICAL          :: HTurbOn       ! Horizontal Turbulence on (.TRUE.) or off (.FALSE.)
  LOGICAL          :: VTurbOn       ! Vertical   Turbulence on (.TRUE.) or off (.FALSE.)
  DOUBLE PRECISION :: ConstantHTurb ! Constant value of horizontal turbulence (m2/s)

  namelist/turbparam/HTurbOn,VTurbOn,ConstantHTurb



!*** BEHAVIOR MODULE PARAMETERS ***
  INTEGER :: Behavior               ! Behavior type (specify a number)
                                    !   Note: The behavior types numbers are: 
                                    !     0 Passive, 1 near-surface, 2 near-bottom, 3 DVM, 
                                    !     4 C.virginica oyster larvae, 5 C.ariakensis oyster larvae, 
                                    !     6 constant, 7 Tidal Stream Transport)
  LOGICAL :: OpenOceanBoundary      ! Note: If you want to allow particles to "escape" via open ocean 
                                    !   boundaries, set this to TRUE; Escape means that the particle 
                                    !   will stick to the boundary and stop moving
  LOGICAL :: mortality              ! TRUE if particles can die; else FALSE
  DOUBLE PRECISION :: deadage       ! Age at which a particle stops moving (i.e., dies) (s)
                                    !   Note: deadage stops particle motion for all behavior types (0-6)
  DOUBLE PRECISION :: pediage       ! Age when particle reaches max swim speed and can settle (s)
                                    !   Note: for oyster larvae behavior (types 4 & 5):
                                    !     pediage = age at which a particle becomes a pediveliger
                                    !   Note: pediage does not cause particles to settle if the Settlement module is not on
  DOUBLE PRECISION :: swimstart     ! Age that swimming or sinking begins (s) 1 day = 1.*24.*3600.
  DOUBLE PRECISION :: swimslow      ! Swimming speed when particle begins to swim (m/s)
  DOUBLE PRECISION :: swimfast      ! Maximum swimming speed (m/s)  0.05 m/s for 5 mm/s
                                    !   Note: for constant swimming speed for behavior types 1,2 & 3: 
                                    !     set swimslow = swimfast = constant speed
  DOUBLE PRECISION :: Sgradient     ! Salinity gradient threshold that cues larval behavior (psu/m)
                                    !   Note: This parameter is only used if Behavior = 4 or 5. 
  DOUBLE PRECISION :: sink          ! Sinking velocity for behavior type 6
                                    !   Note: This parameter is only used if Behavior = 6.
! Tidal Stream Transport behavior type:
  DOUBLE PRECISION :: Hswimspeed    ! Horizontal swimming speed (m/s)
  DOUBLE PRECISION :: Swimdepth     ! Depth at which fish swims during flood time 
                                    ! in meters above bottom (this should be a positive value)


  namelist/behavparam/Behavior,OpenOceanBoundary,mortality,deadage,pediage,swimstart,swimslow,swimfast,Sgradient,sink,Hswimspeed,Swimdepth



!*** DVM. The following are parameters for the Diurnal Vertical Migration (DVM) behavior type:
  DOUBLE PRECISION :: twistart      ! Time of twilight start (hr) **
  DOUBLE PRECISION :: twiend        ! Time of twilight end (hr) **
  DOUBLE PRECISION :: daylength     ! Length of day (hr) **
  DOUBLE PRECISION :: Em            ! Irradiance at solar noon (microE m^-2 s^-1) **
  DOUBLE PRECISION :: Kd            ! Vertical attenuation coefficient
  DOUBLE PRECISION :: thresh        ! Light threshold that cues behavior (microE m^-2 s^-1)
  !  Note: These values were calculated for September 1 at the latitude of 37.0 (Chesapeake Bay mouth)
  !  Note: Variables marked with ** were calculated with light_v2BlueCrab.f (not included in LTRANS yet)
  !  Note: These parameters are only used if Behavior = 3 

  namelist/dvmparam/twistart,twiend,daylength,Em,Kd,thresh



!*** SETTLEMENT MODULE PARAMETERS ***
  LOGICAL :: settlementon           ! settlement module on (.TRUE.) or off (.FALSE.)
  !  Note: If settlement is off: set minholeid, maxholeid, minpolyid, maxpolyid, pedges, & hedges to 1
  !        to avoid both wasted variable space and errors due to arrays of size 0.
  !        If settlement is on and there are no holes: set minholeid, maxholeid, & hedges to 1
  LOGICAL :: holesExist             ! Are there holes in habitat? yes(TRUE) no(FALSE)
  INTEGER :: minpolyid              ! Lowest habitat polygon id number
  INTEGER :: maxpolyid              ! Highest habitat polygon id number
  INTEGER :: minholeid              ! Lowest hole id number
  INTEGER :: maxholeid              ! Highest hole id number
  INTEGER :: pedges                 ! Number of habitat polygon edge points (# of rows in habitat polygon file)
  INTEGER :: hedges                 ! Number of hole edge points (number of rows in holes file)

  namelist/settleparam/settlementon,holesExist,minpolyid,maxpolyid,minholeid,maxholeid,pedges,hedges



!*** CONVERSION MODULE PARAMETERS ***
  DOUBLE PRECISION :: PI            ! Pi
  DOUBLE PRECISION :: Earth_Radius  ! Equatorial radius
  LOGICAL :: SphericalProjection    ! Spherical Projection from ROMS (T/F)
  DOUBLE PRECISION :: lonmin        ! minimum longitude value, only used if SphericalProjection is .TRUE.
  DOUBLE PRECISION :: latmin        ! minimum  latitude value, only used if SphericalProjection is .TRUE.

  namelist/convparam/PI,Earth_Radius,SphericalProjection,lonmin,latmin



!*** INPUT FILE NAME AND LOCATION PARAMETERS ***; 

!  ** ROMS NetCDF Model Grid file **
  CHARACTER(LEN=200) :: NCgridfile
    !Note: the path to the file is necessary if the file is not in the same folder as the code
    !Note: if .nc file in separate folder in Linux, then include path. For example:
    !      CHARACTER(LEN=29), PARAMETER :: NCgridfile = '/share/enorth/CPB_GRID_wUV.nc' 
    !Note: if .nc file in separate folder in Windows, then include path. For example:
    !      CHARACTER(LEN=23), PARAMETER :: NCgridfile = 'D:\ROMS\CPB_GRID_wUV.nc'

  namelist/romsgrid/NCgridfile



!  ** ROMS Predictions NetCDF Input File **
!  Filename = prefix + filenum + suffix
  CHARACTER(LEN=200) :: prefix      ! NetCDF Input Filename prefix
  CHARACTER(LEN=200) :: suffix      ! NetCDF Input Filename suffix
  INTEGER :: filenum                ! Number in First NetCDF Input Filename
  INTEGER :: numdigits              ! Number of digits in number portion of file name (with leading zeros)
  LOGICAL :: startfile              ! .TRUE. means the first file has an additional time step
  !Note: the path to the file is necessary if the file is not in the same folder as the code
  !Note: if .nc file in separate folder in Windows, then include path in prefix. For example:
  !      CHARACTER(LEN=15), PARAMETER :: prefix='D:\ROMS\y95hdr_'   
  !      if .nc file in separate folder in Linux, then include path in prefix. For example:
  !      CHARACTER(LEN=26), PARAMETER :: prefix='/share/lzhong/1995/y95hdr_'   

  namelist/romsoutput/prefix,suffix,filenum,numdigits,startfile




!  ** Particle Location Input File **
  CHARACTER(LEN=200) :: parfile     ! Particle locations file
  !Note: the path to the file is necessary if the file is not in the same folder as the code

  namelist/parloc/parfile



!  ** Habitat Polygon Location Input Files **
  CHARACTER(LEN=200) :: habitatfile ! Habitat polygon file
  CHARACTER(LEN=200) :: holefile    ! Holes in habitat file
  !Note: the path to the file is necessary if the file is not in the same folder as the code

  namelist/HabPolyLoc/habitatfile,holefile



!  ** Output Related Variables **
  CHARACTER(LEN=200) :: outpath     ! Location to write output files
  CHARACTER(LEN=100) :: NCOutFile   ! Name of the NetCDF output file if outputting to NetCDF
  LOGICAL :: outpathGiven           ! If TRUE files are written to the path given in outpath
  LOGICAL :: writeCSV               ! If TRUE write CSV output files
  LOGICAL :: writeNC                ! If TRUE write .NC output files
  INTEGER :: NCtime                 ! Time interval between creation of new NetCDF output files

  !NetCDF Model Metadata:
  CHARACTER(LEN=200) :: SVN_Version ! SVN Repository and Version #
  CHARACTER(LEN=100) :: RunName     ! Unique Identifier for this particular model run
  CHARACTER(LEN=200) :: ExeDir      ! Location of the model run executable
  CHARACTER(LEN=200) :: OutDir      ! Location of the model run output files
  CHARACTER(LEN=100) :: RunBy       ! Name of person who setup/run the model
  CHARACTER(LEN=100) :: Institution ! Place the model is run
  CHARACTER(LEN=200) :: StartedOn   ! Date the model run began

  namelist/output/outpath,NCOutFile,outpathGiven,writeCSV,writeNC,NCtime, &
                  SVN_Version,RunName,ExeDir,OutDir,RunBy,Institution,StartedOn



!*** OTHER PARAMETERS *** 
  INTEGER :: seed                   ! Seed value for random number generator (Mersenne Twister)
  INTEGER :: ErrorFlag              ! What to do if an error is encountered: 0=stop, 1=return particle to previous location
                                    ! 2=kill particle & stop tracking, 3=set particle out of bounds & stop tracking
									! Note: Options 1-3 will output information to ErrorLog.txt
  LOGICAL :: BoundaryBLNs           ! Create Surfer Blanking Files of boundaries? .TRUE.=yes, .FALSE.=no
  LOGICAL :: SaltTempOn             ! Calculate salinity and temperature at particle 
                                    ! location: yes (.TRUE.) or no (.FALSE.)
  LOGICAL :: TrackCollisions        ! Write Bottom and Land Hit Files? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteHeaders           ! Write .txt files with column headers? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteModelTiming       ! Write .csv file with model timing data? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteProblemFile       ! Write a file with problem particles for debugging? .TRUE.=yes, .FALSE.=no

  INTEGER :: ijbuff                 ! number of extra elements to read in on every side of the particles

  LOGICAL :: FreeSlip               ! use free slip condition?

  namelist/other/seed,BoundaryBLNs,SaltTempOn,TrackCollisions,WriteHeaders, &
                 WriteModelTiming,WriteProblemFile,ijbuff,ErrorFlag,FreeSlip

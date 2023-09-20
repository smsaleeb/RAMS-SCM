!##############################################################################
! Implementation of 1D stand-alone RAMS 6.2.13 microphysics,aerosol,Harrington
! radiation Single Column Model (SCM).
!
! Requirements for the SCM and RAMS physics to match:
! 1.Number of land-surface patches must be 2 (NPATCH=2) which is used for
!   dust lofting and aerosol dry deposition.
!##############################################################################
program rams_ccpp

use io_params, only:frqstate
use node_mod, only:mzp,mmzp
use mem_grid, only:iyear1,imonth1,idate1,itime1,ngrid,ngrids,time,dtlt

implicit none

integer :: ng

 !##############################################################################
 ! The SCM needs the following variables as input in order to proceed by default.
 ! At time = 0, the following variables have non-zero values passed in from RAMS.
 ! The value in parentheses next to the variable name is the variable dimension.
 ! mzp=number-of-vert-levels and np=number-of-surface-patch-areas
 !
 ! 01.albedt(1)       - surface albedo (decimal fraction)
 ! 02.dn0(mzp)        - reference state density (kg/m3)
 ! 03.glat(1)         - grid cell latitude (degrees)
 ! 04.glon(1)         - grid cell longitude (degrees)
 ! 05.leaf_class(np)  - vegetation classes (number designation)
 ! 06.patch_area(np)  - fractional area of each vegetation class
 ! 07.pi0(mzp)        - reference state PI (J/(kg*K)) (Exner function * cp)
 ! 08.pp(mzp)         - perturbation PI (J/(kg*K))
 ! 09.pressure(mzp)   - air pressure (Pa)
 ! 10.rlongup(1)      - upwelling longwave radiation from surface (W/m2)
 ! 11.rtp(mzp)        - total water/ice mixing ratio (vapor + condensate) (kg/kg)
 ! 12.rv(mzp)         - water vapor mixing ratio (kg/kg)
 ! 13.soil_text(np)   - soil textural class (number designation)
 ! 14.soil_water(np)  - volumetric soil moisture (m3/m3)
 ! 15.theta(mzp)      - potential temperature (K)
 ! 16.thp(mzp)        - ice-liquide potential temperature (K)
 ! 17.topt(1)         - topographic height (m)
 ! 18.up(mzp)         - u-wind (m/s)
 ! 19.vp(mzp)         - v-wind (m/s)
 ! 20.wp(mzp)         - w-wind (m/s) (=0.0 m/s at time=0.0 seconds)
 ! 21.zm(mzp)         - grid level height on M (W level) (Lorenz vertical grid)
 ! 22.zt(mzp)         - grid level height on T (scalar level) (Lorenz vertical grid)
 !
 ! The following are needed for dust lofting and aerosol deposition. These are
 ! hard coded further down in this file, and then overwritten if data available.
 ! If no "soil_text"  is passed in, "sandy clay loam" is used.
 ! If no "leaf_class" is passed in, "desert-like" is used.
 ! If no "patch_area" is passed in, water=10%, land=90%
 ! if no "soil_water" is passed in, "0.10 m3/m3" is used.
 ! 
 !##############################################################################
 ! Some default values for idealized SCM testing without driver simulation info.
 ! These could be used for custom physics testing.
 ! These are over-written if reading the "parent_sim_info.txt" file.
 iyear1     = 1991     !Start year of the simulation (YYYY)
 imonth1    = 4        !Start month of the simulation
 idate1     = 26       !Start date of the simulation
 itime1     = 2100     !Start time of the simulation (UTC)
 dtlt       = 10.      !Current grid timestep
 time       = 0.       !Current model run time (seconds)
 ngrid      = 1        !Cureent model grid (1 = parent grid)
 mzp        = 40       
 ngrids     = 1        !Total number of grids including any nested grids
 do ng=1,ngrids
  mmzp(ng) = 40        !Number vertical levels on T-grid for each grid
  frqstate(ng) = 600.  !Frequency(/sec) writing output for each grid
 enddo
 !##############################################################################
 ! Read in simulation and grid specific variables in this call. The variables
 ! read-in here should be passed to the RAMS physics if not being read from a
 ! file. See the file "write_scm.f90" subroutine "siminfo" for details.
 CALL siminfo ("parent_sim_info.txt","scm.in/")

 !Initialization of RAMS for model start at time=0 or history restart
 CALL rams_initialize ()

 !Call RAMS physics each timestep (aerosols, radiation, microphysics)
 CALL rams_physics ()

End !finished main program

!##############################################################################
!##############################################################################
!##############################################################################
subroutine rams_initialize ()

use mem_other
use mem_leaf
use mem_basic
use mem_grid
use mem_radiate
use mem_micro
use micphys
use node_mod, only:mmxp,mmyp,mmzp,ia,iz,ja,jz,mxp,myp

implicit none

!##############################################################################
! Local variables for main
!##############################################################################
integer :: ng
integer, save :: firstcall = 0

!Do the following remainder of subroutine only once at model startup
!for initialization, allocation, etc.
if(firstcall==0)then
 firstcall=1

 !##############################################################################
 ! Namelist for microphysics (retain order as in namelist 'micphys.in')
 !##############################################################################
 namelist /mcphys_setting/ numsteps       &
                         ,iprntstatic     &
                         ,initcustom      &
                         ,inegadj         &
                         ,ithermo         &
                         ,ivapdiff        &
                         ,iautoaccret     &
                         ,icldrime        &
                         ,iselfcollragh   &
                         ,iselfcollps     &
                         ,icollps         &
                         ,icolliceice     &
                         ,icollicerain    &
                         ,imelting        &
                         ,icldnuc         &
                         ,iicenuc         &
                         ,isedimentation  &
                         ,ipress          &
                         ,itheta          &
                         ,irv             &
                         ,iprntstmt       &
                         ,iswrtyp         &
                         ,ilwrtyp         &
                         ,radfrq          &
                         ,icheckmic       &
                         ,imbudget        &
                         ,irime           &
                         ,iplaws          &
                         ,isedim          &
                         ,icloud          & 
                         ,idriz           &
                         ,irain           &
                         ,ipris           &
                         ,isnow           &
                         ,iaggr           &
                         ,igraup          &
                         ,ihail           &
                         ,cparm           &
                         ,dparm           &
                         ,rparm           &
                         ,pparm           &
                         ,sparm           &
                         ,aparm           &
                         ,gparm           &
                         ,hparm           &
                         ,gnu             &
                         ,iaerosol        &
                         ,idust           &
                         ,isalt           &
                         ,iabcarb         &
                         ,idustloft       &
                         ,dustfile        &
                         ,iccnlev         &
                         ,iifn            &
                         ,iifn_formula    &
                         ,iaerorad        &
                         ,iaerodep        &
                         ,iaeroprnt       &
                         ,cin_max         &
                         ,ccn1_max        &
                         ,ccn2_max        &
                         ,dust1_max       &
                         ,dust2_max       &
                         ,saltf_max       &
                         ,saltj_max       &
                         ,salts_max       &
                         ,abc1_max        &
                         ,abc2_max        &
                         ,iaero_chem      &
                         ,aero_epsilon    &
                         ,aero_medrad     &
                         ,itrkepsilon     &
                         ,itrkdust        &
                         ,itrkdustifn

 !##############################################################################
 ! Assign some local values used in physics
 !##############################################################################
 !Surface info settings for setting soil, water, land class dimensions needed
 !for surface aerosol deposition, dust lofting, sea-salt aerosol generation. 
 nzg = 1    ! soil levels (default to 1 here)
 nzs = 1    ! surface water/snow levels (default to 1 here)
 npatch = 2 ! Forcing 2 surface patches. 1=water, 2=land
 !Microphysics and radiation flags common to RAMS. Keep the same here for
 !consistency. Could carefully changes these later, but not necessary.
 print_msg  = .true.   ! print flag associated with parallel Node-1 in RAMS
 level      = 3        ! Microphysics level to always equal 3
 writescm   = 1        ! 0 = no SCM output files, 1= yes SCM output files

 !##############################################################################
 ! Read microphysics namelist
 !##############################################################################
 if(iprntstmt>=1 .and. print_msg) then
 write(*,*)
 write(*,*)'***************************************************************'
 write(*,*) 'Reading namelist...'
 write(*,*)
 endif
 open(10, file = 'RAMSIN', status = 'old' )
 read(10, nml = mcphys_setting )
 close(10)

 !##############################################################################
 ! MICROPHYSICS AND AEROSOL FLAG CHECKING SECTION
 !##############################################################################
 IF( ((icloud .GE. 2 .AND. icloud .LE. 4 .AND. cparm .LE. 0.)  &
  .OR.(idriz  .GE. 2 .AND. idriz  .LE. 4 .AND. dparm .LE. 0.)  &
  .OR.(irain  .GE. 2 .AND. irain  .LE. 4 .AND. rparm .LE. 0.)  &
  .OR.(ipris  .GE. 2 .AND. ipris  .LE. 4 .AND. pparm .LE. 0.)  &
  .OR.(isnow  .GE. 2 .AND. isnow  .LE. 4 .AND. Sparm .LE. 0.)  &
  .OR.(igraup .GE. 2 .AND. igraup .LE. 4 .AND. gparm .LE. 0.)  &
  .OR.(iaggr  .GE. 2 .AND. iaggr  .LE. 4 .AND. Aparm .LE. 0.)  &
  .OR.(ihail  .GE. 2 .AND. ihail  .LE. 4 .AND. hparm .LE. 0.)) &
  .AND. level <=3) THEN
   print*,' '
   print*,'FATAL - Microphysics - xPARM must be positive'  &
             ,'if micro flags are set to 2, 3, or 4:'
   print*,'cparm:',cparm
   print*,'dparm:',dparm
   print*,'rparm:',rparm
   print*,'pparm:',pparm
   print*,'sparm:',sparm
   print*,'aparm:',aparm
   print*,'gparm:',gparm
   print*,'hparm:',hparm
   print*,' '
   stop
 endif

 if (iaerosol .lt. 0 .or. iaerosol .gt. 1) THEN
    print*,'FATAL - IAEROSOL OUT OF RANGE: MUST BE 0-1'
    stop
 endif
 if (idust .lt. 0 .or. idust .gt. 2) THEN
    print*,'FATAL - IDUST OUT OF RANGE: MUST BE 0-2'
    stop
 endif
 if (isalt .lt. 0 .or. isalt .gt. 2) THEN
    print*,'FATAL - ISALT OUT OF RANGE: MUST BE 0-2'
    stop
 endif
 if (iabcarb .lt. 0 .or. iabcarb .gt. 1) THEN
    print*,'FATAL - IABCARB OUT OF RANGE: MUST BE 0-1'
    stop
 endif
 if (iaerorad .lt. 0 .or. iaerorad .gt. 1) THEN
    print*,'FATAL - IAERORAD OUT OF RANGE: MUST BE 0-1'
    stop
 endif
 if ((iswrtyp.ne.0 .and. iswrtyp.ne.3) .or. &
     (ilwrtyp.ne.0 .and. ilwrtyp.ne.3)) then
    print*,'Fatal - ISWRTYP and ILWRTYP must = 0 or 3'
    stop
 endif
 if (iaerorad .eq. 1 .and. (iswrtyp.ne.3 .or. ilwrtyp.ne.3)) THEN
    print*,'FATAL - Aerosol radiation turned on but Harrington'
    print*,'        Radiation scheme is not on.'
    stop
 endif
 !Guidance on hydrometeor shape parameter GNU
 if (icloud .gt. 0 .and. gnu(1) < 4) then
  print*,'FATAL - GNU for Cloud droplets is strongly recommended >= 4'
  stop
 endif
 if (idriz .gt. 0 .and. gnu(8) < 4) then
  print*,'FATAL - GNU for Drizzle droplets is strongly recommended >= 4'
  stop
 endif
 !ALLOWABLE RANGE OF MICRO / RADIATION FLAGS
 if (icloud .lt. 0 .or. icloud .gt. 5) THEN
  print*,'FATAL - ICLOUD OUT OF RANGE'
  stop
 endif
 if (irain .lt. 0 .or. irain .gt. 5) THEN
  print*,'FATAL - IRAIN OUT OF RANGE'
  stop
 endif
 if (ipris .lt. 0 .or. ipris .gt. 5) THEN
  print*,'FATAL - IPRIS OUT OF RANGE'
  stop
 endif
 if (ipris .ge. 1 .and. ipris .le. 4) THEN
  print*,'FATAL - IPRIS MUST BE 0 or 5'
  stop
 endif
 if (isnow .lt. 0 .or. isnow .gt. 5) THEN
  print*,'FATAL - ISNOW OUT OF RANGE'
  stop
 endif
 if (iaggr .lt. 0 .or. iaggr .gt. 5) THEN
  print*,'FATAL - IAGGR OUT OF RANGE'
  stop
 endif
 if (igraup .lt. 0 .or. igraup .gt. 5) THEN
  print*,'FATAL - IGRAUP OUT OF RANGE'
  stop
 endif
 if (ihail .lt. 0 .or. ihail .gt. 5) THEN
  print*,'FATAL - IHAIL OUT OF RANGE'
  stop
 endif
 if (idriz .lt. 0 .or. idriz .gt. 5) THEN
  print*,'FATAL - IDRIZ OUT OF RANGE'
  stop
 endif
 if (iccnlev .lt. 0 .or. iccnlev .gt. 2) THEN
  print*,'FATAL - ICCNLEV OUT OF RANGE: MUST BE 0-2'
  stop
 endif
 if (imbudget .lt. 0 .or. imbudget .gt. 3) THEN
  print*,'FATAL - IMBUDGET OUT OF RANGE'
  stop
 endif
 if (iifn .lt. 0 .or. iifn .gt. 3) THEN
  print*,'FATAL - IIFN OUT OF RANGE: MUST BE 0-3'
  stop
 endif
 if (isedim .lt. 0 .or. isedim .gt. 1) THEN
  print*,'FATAL - ISEDIM OUT OF RANGE: MUST BE 0-1'
  stop
 endif
 if (itrkepsilon .lt. 0 .or. itrkepsilon .gt. 1 .or. &
     itrkdust    .lt. 0 .or. itrkdust    .gt. 1 .or. &
     itrkdustifn .lt. 0 .or. itrkdustifn .gt. 1) THEN
  print*,'FATAL - AEROSOL TRACKING FLAGS MUST BE 0-1'
  stop
 endif

 !One-moment cloud water checking
 if (icloud .lt. 5) then
  if (idriz .ge. 5) then
    print*,'FATAL - Microphysics - IDRIZ must be < 5 if ICLOUD < 5'
    stop
  endif
  if (iccnlev .gt. 0) THEN
    print*,'WARNING - iccnlev>0 and icloud<5.'
    print*,' No aerosol nucleation for these settings.'
    print*,' Aerosol removal/tracking now only active via deposition.'
  endif
  if (iifn .eq. 3) THEN
    print*,'FATAL - IIFN=3 requires ICLOUD>=5 for there to be any'
    print*,' primary heterogeneous ice nucleation.'
    stop
  endif
 endif

 !Two-moment cloud water checking
 if (icloud .eq. 5) then
  if (iaerosol.eq.0 .and. idust.eq.0 .and. isalt.eq.0 .and. iabcarb.eq.0) then
    print*,'Warning - If icloud==5, then either IAEROSOL, IDUST, IABCARB, or ISALT'
    print*,' must be turned on so that some aerosols are present in order'
    print*,' for cloud nucleation to occur.'
  endif
 endif

 if (imbudget .ge. 3 .and. iccnlev .eq. 0) THEN
   print*,'FATAL - If IMBUDGET=3 then ICCNLEV must be > 0.'
   print*,' Cannot track specific aerosol type nucleation'
   print*,' unless ICCNLEV>0 for nucleation scavenging.'
   stop
 endif
 if (imbudget .ge. 3 .and. idust .eq. 0) THEN
   print*,'FATAL - If IMBUDGET=3 then IDUST must be > 0.'
   print*,' Cannot track dust if there is no dust turned on.'
   stop
 endif
 if ((itrkepsilon.gt.0 .or. itrkdust.gt.0 .or. itrkdustifn.gt.0) &
      .and. iccnlev.lt.2) THEN
   print*,'FATAL - Aerosol tracking = 0 since ICCNLEV < 2.'
   print*,'Either set ICCNLEV>=2 or set ITRKEPSILON=ITRKDUST=ITRKDUSTIFN=0'
   stop
 endif
 if ((itrkdust.gt.0 .or. itrkdustifn.gt.0) .and. idust.eq.0) THEN
   print*,'FATAL - Dust tracking = 0 since IDUST=0.'
   print*,'Either set IDUST>0 or set ITRKDUST=ITRKDUSTIFN=0'
   stop
 endif

 !##############################################################################
 !Assign hydrometeor modes. This needs to be called before allocating main
 !variables since we need namelist flags and hydrometeor mode numbers to
 !determine which micro variables to allocate.
 if(iprntstmt>=1 .and. print_msg) then
 write(*,*)
 write(*,*)'***************************************************************'
 write(*,*) 'Assign hydrometeor modes so we know how to allocate memory...'
 write(*,*)
 endif
 CALL jnmbinit ()

 !##############################################################################
 !Allocate main variables
 if(iprntstmt>=1 .and. print_msg) then
 write(*,*)
 write(*,*)'***************************************************************'
 write(*,*) 'Dynamic memory allocation based on flags and grid size...'
 write(*,*)
 endif
 !Set X and Y dimensions to 1 grid cell for each model domain grid 
 !(parent + nested grids). This is set to work as a single column, so all
 !3D RAMS allocated physics arrays set as (z,x,y) will be allocated as (z,1,1) 
 !where z can vary in the column between domain grids. Also need to make the x,y 
 !dimensions equal to 1 in the event that the parent model does not use a z,x,y 
 !set of dimensions, but rather, something unique to a global gridded domain 
 !such as an icosahedral grid.
 do ng=1,ngrids
  mmyp(ng) = 1  !Set to 1 Y-direction grid cell
  mmxp(ng) = 1  !Set to 1 X-direction grid cell
 enddo
 mxp = 1 !mxp = current grid number X,I points
 myp = 1 !myp = current grid number Y,J points
 !Max value of vertical grid cells from among parent and nested grids, 
 !used in "McLatchy" radiation sounding for adding radiation levels if domain
 !top is lower than about 24km.
 maxnzp = 0
 do ng=1,ngrids
  maxnzp = max(maxnzp,mmzp(ng))
 enddo
 !These are used in RAMS for parallel processing on sub-domains over X & Y.
 !For true single column implementation (SCM) these should = 1.
 ia = 1     ! X-points start in sub-domain to loop over (= 1 for true SCM)
 iz = mxp   ! X-points stop  in sub-domain to loop over (= 1 for true SCM)
 ja = 1     ! Y-points start in sub-domain to loop over (= 1 for true SCM)
 jz = myp   ! Y-points stop  in sub-domain to loop over (= 1 for true SCM)
 !Allocating grid datatypes
 if(iprntstmt>=1 .and. print_msg)print*,'start other alloc'
 allocate(other_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_other (other_g(ng))
    CALL alloc_other (other_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),npatch)
 enddo
 !Allocating land/ocean surface info datatypes
 if(iprntstmt>=1 .and. print_msg)print*,'start leaf/land alloc'
 allocate(leaf_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_leaf (leaf_g(ng))
    CALL alloc_leaf (leaf_g(ng),mmxp(ng),mmyp(ng),nzg,nzs,npatch)
 enddo
 if(iprntstmt>=1 .and. print_msg)print*,'start grid alloc'
 allocate(grid_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_grid (grid_g(ng))
    CALL alloc_grid (grid_g(ng),mmxp(ng),mmyp(ng))
 enddo
 !Allocating dynamics datatypes
 if(iprntstmt>=1 .and. print_msg)print*,'start basic alloc'
 allocate(basic_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_basic (basic_g(ng))
    CALL alloc_basic (basic_g(ng),mmzp(ng),mmxp(ng),mmyp(ng))
 enddo
 !Allocating micro datatypes
 if(iprntstmt>=1 .and. print_msg)print*,'start micro alloc'
 allocate(micro_g(ngrids))
 allocate(pcp_tab(ngrids))
 do ng=1,ngrids
    CALL dealloc_micro (micro_g(ng))
    CALL dealloc_sedim (pcp_tab(ng))
    CALL alloc_micro (micro_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),1)
    CALL alloc_sedim (pcp_tab(ng),mmzp(ng))
 enddo
 !Allocating radiation datatypes
 if(iprntstmt>=1 .and. print_msg)print*,'start radiate alloc'
 allocate(radiate_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_radiate (radiate_g(ng))
    CALL alloc_radiate (radiate_g(ng),mmzp(ng),mmxp(ng),mmyp(ng))
 enddo

 !##############################################################################
 !Initialize surface vegetation class, soil class, and water/land class area. 
 !These are over-written below if read in or passed in from parent model.
 !These are needed for dust lofting, sea-salt production, and aerosol dry
 !deposition onto the surface. We assume 2 land patches. Patch 1 = water and
 !Patch 2 = land.
 !##############################################################################
 do ng=1,ngrids
  ! Water patch fraction (real)
  leaf_g(ng)%patch_area(:,:,1) = 0.10
  ! Land patch fraction (real)
  leaf_g(ng)%patch_area(:,:,2) = 0.90
  ! Land patch type = 3 = desert-like [real; uses nint(value)]
  leaf_g(ng)%leaf_class(:,:,1) = 0. !water patch
  leaf_g(ng)%leaf_class(:,:,2) = 3.
  ! Volumetric soil moisture (m3/m3) (real)
  leaf_g(ng)%soil_water(nzg,:,:,1) = 1.00 !water patch
  leaf_g(ng)%soil_water(nzg,:,:,2) = 0.10
  ! Soil textural class = 6 = sandy clay loam [real; uses nint(value)]
  leaf_g(ng)%soil_text(nzg,:,:,1) = 0. !water patch
  leaf_g(ng)%soil_text(nzg,:,:,2) = 6.
 enddo

 !##############################################################################
 ! Write out the namelist variables (retain order as in namelist)
 !##############################################################################
 if(iprntstmt>=1 .and. print_msg) then
 write(*,*)
 write(*,*)'***************************************************************'
 write(*,*)'                    NAMELIST VARIABLES                         '
 write(*,*)'See documentation file on namelist parameters & RAMS variables '
 write(*,*)'***************************************************************'
 write(*,'(1X,A24,I6)')    ,'NUMSTEPS                ',numsteps
 write(*,'(1X,A24,I6)')    ,'IPRNTSTATIC             ',iprntstatic
 write(*,'(1X,A24,I6)')    ,'INITCUSTOM              ',initcustom
 write(*,'(1X,A24,I6)')    ,'INEGADJ                 ',inegadj
 write(*,'(1X,A24,I6)')    ,'ITHERMO                 ',ithermo
 write(*,'(1X,A24,I6)')    ,'IVAPDIFF                ',ivapdiff
 write(*,'(1X,A24,I6)')    ,'IAUTOACCRET             ',iautoaccret
 write(*,'(1X,A24,I6)')    ,'ICLDRIME                ',icldrime
 write(*,'(1X,A24,I6)')    ,'ISELFCOLLRAGH           ',iselfcollragh
 write(*,'(1X,A24,I6)')    ,'ISELFCOLLPS             ',iselfcollps
 write(*,'(1X,A24,I6)')    ,'ICOLLPS                 ',icollps
 write(*,'(1X,A24,I6)')    ,'ICOLLICEICE             ',icolliceice
 write(*,'(1X,A24,I6)')    ,'ICOLLICERAIN            ',icollicerain
 write(*,'(1X,A24,I6)')    ,'IMELTING                ',imelting
 write(*,'(1X,A24,I6)')    ,'ICLDNUC                 ',icldnuc
 write(*,'(1X,A24,I6)')    ,'IICENUC                 ',iicenuc
 write(*,'(1X,A24,I6)')    ,'ISEDIMENTATION          ',isedimentation
 write(*,'(1X,A24,I6)')    ,'ISWRTYP                 ',iswrtyp
 write(*,'(1X,A24,I6)')    ,'ILWRTYP                 ',ilwrtyp
 write(*,'(1X,A24,I6)')    ,'IPRESS                  ',ipress
 write(*,'(1X,A24,I6)')    ,'ITHETA                  ',itheta
 write(*,'(1X,A24,I6)')    ,'IRV                     ',irv
 write(*,'(1X,A24,I6)')    ,'IPRNTSTMT               ',iprntstmt
 write(*,'(1X,A24,F11.4)') ,'RADFRQ                  ',radfrq
 write(*,'(1X,A24,I6)')    ,'ICHECKMIC               ',icheckmic
 write(*,'(1X,A24,I6)')    ,'IMBUDGET                ',imbudget
 write(*,'(1X,A24,I6)')    ,'IRIME                   ',irime
 write(*,'(1X,A24,I6)')    ,'IPLAWS                  ',iplaws
 write(*,'(1X,A24,I6)')    ,'ISEDIM                  ',isedim
 write(*,'(1X,A24,I6)')    ,'ICLOUD                  ',icloud
 write(*,'(1X,A24,I6)')    ,'IDRIZ                   ',idriz
 write(*,'(1X,A24,I6)')    ,'IRAIN                   ',irain
 write(*,'(1X,A24,I6)')    ,'IPRIS                   ',ipris
 write(*,'(1X,A24,I6)')    ,'ISNOW                   ',isnow
 write(*,'(1X,A24,I6)')    ,'IAGGR                   ',iaggr
 write(*,'(1X,A24,I6)')    ,'IGRAUP                  ',igraup
 write(*,'(1X,A24,I6)')    ,'IHAIL                   ',ihail
 write(*,'(1X,A24,E15.4)') ,'CPARM                   ',cparm
 write(*,'(1X,A24,E15.4)') ,'DPARM                   ',dparm
 write(*,'(1X,A24,E15.4)') ,'PPARM                   ',pparm
 write(*,'(1X,A24,E15.4)') ,'SPARM                   ',sparm
 write(*,'(1X,A24,E15.4)') ,'APARM                   ',aparm
 write(*,'(1X,A24,E15.4)') ,'GPARM                   ',gparm
 write(*,'(1X,A24,E15.4)') ,'HPARM                   ',hparm
 write(*,'(1X,A24,F8.1)')  ,'GNU-cloud               ',gnu(1)
 write(*,'(1X,A24,F8.1)')  ,'GNU-driz                ',gnu(8)
 write(*,'(1X,A24,F8.1)')  ,'GNU-rain                ',gnu(2)
 write(*,'(1X,A24,F8.1)')  ,'GNU-pris                ',gnu(3)
 write(*,'(1X,A24,F8.1)')  ,'GNU-snow                ',gnu(4)
 write(*,'(1X,A24,F8.1)')  ,'GNU-aggr                ',gnu(5)
 write(*,'(1X,A24,F8.1)')  ,'GNU-graup               ',gnu(6)
 write(*,'(1X,A24,F8.1)')  ,'GNU-hail                ',gnu(7)
 write(*,'(1X,A24,I6)')    ,'IAEROSOL                ',iaerosol
 write(*,'(1X,A24,I6)')    ,'IDUST                   ',idust
 write(*,'(1X,A24,I6)')    ,'ISALT                   ',isalt
 write(*,'(1X,A24,I6)')    ,'IABCARB                 ',iabcarb
 write(*,'(1X,A24,I6)')    ,'IDUSTLOFT               ',idustloft
 write(*,'(1X,A24,5X,A)')  ,'DUSTFILE                ',trim(dustfile)
 write(*,'(1X,A24,I6)')    ,'ICCNLEV                 ',iccnlev
 write(*,'(1X,A24,I6)')    ,'IIFN                    ',iifn
 write(*,'(1X,A24,I6)')    ,'IIFN_FORMULA            ',iifn_formula
 write(*,'(1X,A24,I6)')    ,'IAERORAD                ',iaerorad
 write(*,'(1X,A24,I6)')    ,'IAERODEP                ',iaerodep
 write(*,'(1X,A24,I6)')    ,'IAEROPRNT               ',iaeroprnt
 write(*,'(1X,A24,F15.8)') ,'CIN_MAX                 ',cin_max
 write(*,'(1X,A24,F15.8)') ,'CCN1_MAX                ',ccn1_max
 write(*,'(1X,A24,F15.8)') ,'CCN2_MAX                ',ccn2_max
 write(*,'(1X,A24,F15.8)') ,'DUST1_MAX               ',dust1_max
 write(*,'(1X,A24,F15.8)') ,'DUST2_MAX               ',dust2_max
 write(*,'(1X,A24,F15.8)') ,'SALTF_MAX               ',saltf_max
 write(*,'(1X,A24,F15.8)') ,'SALTJ_MAX               ',saltj_max
 write(*,'(1X,A24,F15.8)') ,'SALTS_MAX               ',salts_max
 write(*,'(1X,A24,F15.8)') ,'ABC1_MAX                ',abc1_max
 write(*,'(1X,A24,F15.8)') ,'ABC2_MAX                ',abc2_max
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(ccn-1)       ',iaero_chem(1)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(ccn-2)       ',iaero_chem(2)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(dust-1)      ',iaero_chem(3)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(dust-2)      ',iaero_chem(4)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(salt-film)   ',iaero_chem(5)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(salt-jet)    ',iaero_chem(6)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(salt-spume)  ',iaero_chem(7)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(dust-1)      ',iaero_chem(8)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(dust-2)      ',iaero_chem(9)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(abcarb-1)    ',iaero_chem(10)
 write(*,'(1X,A24,I6)')    ,'IAERO_CHEM(abcarb-2)    ',iaero_chem(11)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(ccn-1)     ',aero_epsilon(1)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(ccn-2)     ',aero_epsilon(2)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(dust-1)    ',aero_epsilon(3)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(dust-2)    ',aero_epsilon(4)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(salt-film) ',aero_epsilon(5)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(salt-jet)  ',aero_epsilon(6)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(salt-spume)',aero_epsilon(7)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(abcarb-1)  ',aero_epsilon(8)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(abcarb-2)  ',aero_epsilon(9)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(regen-1)   ',aero_epsilon(10)
 write(*,'(1X,A24,F9.2)')  ,'AERO_EPSILON(regen-2)   ',aero_epsilon(11)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(ccn-1)      ',aero_medrad(1)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(ccn-2)      ',aero_medrad(2)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(dust-1)     ',aero_medrad(3)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(dust-2)     ',aero_medrad(4)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(salt-film)  ',aero_medrad(5)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(salt-jet)   ',aero_medrad(6)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(salt-spume) ',aero_medrad(7)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(abcarb-1)   ',aero_medrad(8)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(abcarb-2)   ',aero_medrad(9)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(regen-1)    ',aero_medrad(10)
 write(*,'(1X,A24,E15.4)') ,'AERO_MEDRAD(regen-2)    ',aero_medrad(11)
 write(*,'(1X,A24,I6)')    ,'ITRKEPSILON             ',itrkepsilon
 write(*,'(1X,A24,I6)')    ,'ITRKDUST                ',itrkdust
 write(*,'(1X,A24,I6)')    ,'ITRKDUSTIFN             ',itrkdustifn
 write(*,*)
 write(*,*)'***************************************************************'
 write(*,*)'            HARDCODED SETTINGS NEEDED FOR SCM                  '
 write(*,*)'***************************************************************'
 write(*,'(1X,A24,I6)')    ,'LEVEL                   ',level
 write(*,'(1X,A24,I6)')    ,'NPATCH                  ',npatch
 write(*,'(1X,A24,I6)')    ,'NZG                     ',nzg
 write(*,'(1X,A24,I6)')    ,'NZS                     ',nzs
 endif

 !##############################################################################
 !Initialize aerosol density and vanthoff factors if they are used.
 if(iaerosol>0 .or. idust>0 .or. isalt>0 .or. iabcarb>0) then
  if(iprntstmt>=1 .and. print_msg) then
  write(*,*)
  write(*,*)'***************************************************************'
  write(*,*) 'Initializing aerosol characteristic factors'
  write(*,*)'***************************************************************'
  write(*,*)
  endif
  CALL aerosol_init ()
 endif

 !##############################################################################
 !Need this for full microphysics so that we have cfmas and pwmas in
 !analysis files for estimating hydrometeor diameters from mass
 !and number. Level=3 needs the full init for bulk microphysics.
 if(iprntstmt>=1 .and. print_msg) then
 write(*,*)
 write(*,*)'***************************************************************'
 write(*,*) 'Initializing some microphysics basic parameters and tables'
 write(*,*)'***************************************************************'
 write(*,*)
 endif
 CALL micro_master ()

endif !end if firstcall only at initialization of simulation

return
End Subroutine rams_initialize

!##############################################################################
!##############################################################################
!##############################################################################
subroutine rams_physics ()

use mem_grid, only:time,dtlt,iprntstmt,print_msg
use micphys, only:numsteps,writescm

implicit none

!##############################################################################
! Local variables for main
!##############################################################################
integer :: ng,i,j,k,istp
istp=0

!##############################################################################
!Read in single column variables from ASCII files. In couple CCPP model, this
!would be where variables would be passed from parent model to SCM.
if(writescm==1) then
 CALL system ('rm -f scm.out/*')
 CALL system ('rm -f scm.gem/out.*')
 CALL readwrite_scm (1,"scm.in/",istp)
endif

!##############################################################################
!Compute some grid information that RAMS would normally do at initialization
CALL init_rtgt_dz ()

!##############################################################################
!Initialize microphysics and aerosols and custom settings.
!Setting theta=thp and rtp=rv (or vice versa) depending on parent model input.
!Setting aerosol profiles depending on input sources.
!Setting any customizations for physics testing (for INITCUSTOM==1).
if(time == 0.0) then
 if(iprntstmt>=1 .and. print_msg) then
  write(*,*)
  write(*,*)'***************************************************************'
  write(*,*) 'Initializing thermo, aerosols, dust-frac, customizations'
  write(*,*)'***************************************************************'
 endif
 CALL init_aero_custom ()
endif !if time zero, do initialization

!##############################################################################
!Call radiation, aerosols, thermo, negative adjustment, and microphysics
!##############################################################################
do while (istp < numsteps)

 istp = istp + 1
 CALL predthp    () !Update theta-il from radiative temperature tendency.
 CALL radiate    () !Run condensate- and aerosol-sensitive Harrington radiation.
 CALL negadj1    () !Make sure moisture, condensate, aerosols are in bounds.
 CALL thermo     () !Update theta & rv following moisture adjustment in negadj1.
 CALL aerosols   () !Update aerosols via sources for dust and sea-salt.
 CALL micro      () !Run Level-3 RAMS microphysics.
 CALL trsets_ns  () !Set non-scalar top and bottom boundary conditions.
 CALL negadj1    () !Re-run moisture, condensate, aerosol bounds adjustment.
 CALL thermo     () !Recompute theta and rv in case adjustments were done.
 CALL checkmicro () !Check to make sure micro variables are ok.

 time = time + dtlt !increment time for offline version

 if(iprntstmt>=1 .and. print_msg) then
 write(*,*)
 write(*,*) 'Timestep-',istp,'Time-',time
 write(*,*)
 endif

 !Write out single column variables to ASCII files. In coupled CCPP model, this
 !would be where variables would be passed from SCM to parent model.
 if(writescm==1) CALL readwrite_scm (2,"scm.out/",istp)

 if(iprntstmt>=1 .and. print_msg) then
 if(istp >= numsteps)then
  write(*,*)
  write(*,*) 'Finished aerosols, microphysics, and radiation.'
  write(*,*)
 endif
 endif

enddo ! do while (istp < numsteps)

return
End Subroutine rams_physics

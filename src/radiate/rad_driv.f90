!##############################################################################
Subroutine radiate ()

use mem_grid
use mem_radiate
use mem_basic
use mem_micro
use rconstants
use rrad3
use micphys
use node_mod, only:mzp,mxp,myp,ia,iz,ja,jz

implicit none

integer, save :: ncall=0

real, dimension(nzpmax) ::  vctr1 ,vctr2 ,vctr3 ,vctr4 ,vctr5 ,vctr6   &
                           ,vctr7 ,vctr8 ,vctr9 ,vctr10,vctr11,vctr12

if (ilwrtyp + iswrtyp .eq. 0) return

CALL tend_copy (mzp,mxp,myp,ia,iz,ja,jz,radiate_g(ngrid)%fthrdp(1,1,1)  &
   ,radiate_g(ngrid)%fthrd(1,1,1))

if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. 0.001) then

   if(iprntstmt>=1 .and. print_msg) &
      print 90,time,time/3600.+(itime1/100+mod(itime1,100)/60.)
90      format(' Radiation Tendencies Updated Time =',F10.1,  &
        '  UTC TIME (HRS) =',F9.4)

   ! Compute solar zenith angle, multiplier for solar constant, sfc albedo,
   ! and surface upward longwave radiation. radprep in RAMS coupled to the 
   ! LEAF-3 model provides zenith angle, surface downward shortwave and 
   ! longwave for land-surface model in call to "sfcrad". Then "sfcrad" 
   ! returns albedo, upward longwave, and upward shortwave from surface.
   ! "sfcrad" is not called in this SCM version since we are not using LEAF-3
   ! land surface, but will end up using whatever CCPP provides surface model.

   CALL radprep (mxp,myp,ia,iz,ja,jz,jday       &
      ,solfac                                   &
      ,grid_g(ngrid)%glat            (1,1)      &
      ,grid_g(ngrid)%glon            (1,1)      &
      ,radiate_g(ngrid)%rlongup      (1,1)      &
      ,radiate_g(ngrid)%albedt       (1,1)      &
      ,radiate_g(ngrid)%cosz         (1,1))

   ! Using Harrington radiation

   if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then

      ! If first call for this node, initialize several quantities & Mclatchy
      ! sounding data.

      if (ncall .eq. 0) then
         CALL radinit (ng,nb,nsolb,npsb,nuum,prf,alpha,trf,beta  &
            ,xp,wght,wlenlo,wlenhi,solar0,ralcs,a0,a1,a2,a3  &
            ,exptabc,ulim,npartob,npartg,ncog,ncb  &
            ,ocoef,bcoef,gcoef,gnu)

         CALL mclatchy (1,mzp  &
            ,grid_g(ngrid)%glat       (1,1)  &
            ,grid_g(ngrid)%rtgt       (1,1)  &
            ,grid_g(ngrid)%topt       (1,1)  &
            ,radiate_g(ngrid)%rlongup (1,1)  &
            ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
            ,vctr8,vctr9,vctr10,vctr11,vctr12 &
            )

         ncall = ncall + 1
      endif

     ! For any call, interpolate the mclatchy sounding data by latitude and
     ! season.

     CALL mclatchy (2,mzp  &
         ,grid_g(ngrid)%glat       (1,1)  &
         ,grid_g(ngrid)%rtgt       (1,1)  &
         ,grid_g(ngrid)%topt       (1,1)  &
         ,radiate_g(ngrid)%rlongup (1,1)  &
         ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
         ,vctr8,vctr9,vctr10,vctr11,vctr12 &
         )

   endif

endif

return
END SUBROUTINE radiate

!##############################################################################
Subroutine tend_copy (m1,m2,m3,ia,iz,ja,jz,at,at2)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: at,at2

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         at(k,i,j) = at2(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE tend_copy

!##############################################################################
Subroutine radprep (m2,m3,ia,iz,ja,jz,jday   &
   ,solfac,glat,glon,rlongup,albedt,cosz)
   
use rconstants

implicit none

integer :: m2,m3,ia,iz,ja,jz,jday
real :: solfac
real, dimension(m2,m3) :: glat,glon,rlongup,albedt,cosz

integer :: ip,i,j

! Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].

CALL zen (m2,m3,ia,iz,ja,jz,jday,glat,glon,cosz,solfac)

return
END SUBROUTINE radprep

!##############################################################################
Subroutine zen (m2,m3,ia,iz,ja,jz,jday,glat,glon,cosz,solfac)

use mem_grid
use mem_radiate
use rconstants

implicit none

integer :: m2,m3,ia,iz,ja,jz,jday,i,j,oyr,omn,ody,otm
integer, external :: julday

real :: solfac,sdec,cdec,declin,d0,d02,dayhr,radlat,cslcsd,snlsnd  &
   ,gglon,dayhrr,hrangl,pisolar,eqt_julian_nudge,otmf
real, dimension(m2,m3) :: glat,glon,cosz

! Find the hour angle and julian day, then get cosine of zenith angle.
CALL date_add_to (iyear1,imonth1,idate1,itime1*100,time,'s',oyr,omn,ody,otm)

!print*,'Init-Start:(yr,mo,day,utc): ',iyear1,imonth1,idate1,itime1
!print*,'Elapsed time(sec): ',time

! Fix the julian day of the year for RCE simulations for const radiation
!SaleebySCM:no RCE flags being set
!if(irce == 1) then
! CALL date_add_to (iyear1,imonth1,idate1,itime1*100,0.,'s',oyr,omn,ody,otm)
!endif

otmf=float(otm)
dayhr = (otmf-mod(otmf,10000.))/10000. + (mod(otmf,10000.)/6000.)
jday = julday(omn,ody,oyr)

!print*,'Current (yr,mo,day,day-hour):',oyr,omn,ody
!print*,'Current Julian-Day, Day-hour:',jday,dayhr

!      sdec - sine of declination, cdec - cosine of declination
declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
sdec = sin(declin)
cdec = cos(declin)

! Find the factor, solfac, to multiply the solar constant to correct
! for Earth's varying distance to the sun.

d0 = 6.2831853 * float(jday-1) / 365.
d02 = d0 * 2.
solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
   + 0.000719 * cos(d02) + 0.000077 * sin(d02)

! Find day of year equation of time adjustment (Saleeby2008: improve accuracy)
pisolar=3.1415927
eqt_julian_nudge=0.0
if(jday>=1 .and. jday<=106) then
 eqt_julian_nudge = -14.2 * sin(pisolar * (jday +   7) / 111.)
elseif(jday>=107 .and. jday<=166) then
 eqt_julian_nudge =   4.0 * sin(pisolar * (jday - 106) /  59.)
elseif(jday>=167 .and. jday<=246) then
 eqt_julian_nudge =  -6.5 * sin(pisolar * (jday - 166) /  80.)
elseif(jday>=247 .and. jday<=366) then
 eqt_julian_nudge =  16.4 * sin(pisolar * (jday - 247) / 113.)
endif
eqt_julian_nudge = eqt_julian_nudge / 60.0

do j = ja,jz
   do i = ia,iz
      radlat = glat(i,j) * pi180
!SaleebySCM:always using lonrad=1
      !if (lonrad .eq. 0) radlat = centlat(1) * pi180
      if (radlat .eq. declin) radlat = radlat + 1.e-5
      cslcsd = cos(radlat) * cdec
      snlsnd = sin(radlat) * sdec
      gglon = glon(i,j)
!SaleebySCM:always using lonrad=1
      !if (lonrad .eq. 0) gglon = centlon(1)

      !Saleeby(2009): Add eqt_julian_nudge to improve hour angle accuracy
      dayhrr = mod(dayhr+gglon/15.+24.,24.) + eqt_julian_nudge
      hrangl = 15. * (dayhrr - 12.) * pi180

      !Cosine of zenith angle (angle in radians)
      cosz(i,j) = snlsnd + cslcsd * cos(hrangl)

      !Hold zenith angle constant for RCE simulations
!SaleebySCM:not RCE flags being set
      !if(irce == 1) cosz(i,j) = cos(rce_szen * pi180)

   enddo
enddo

return
END SUBROUTINE zen

!##############################################################################
Subroutine radcalc3 (m1,i,j,ngrid,maxnzp,mcat,iswrtyp,ilwrtyp,zm,zt &
   ,glat,rtgt,topt,rv,albedt,cosz,rlongup,rshort,rlong,aodt &
   ,fthrd,bext,swup,swdn,lwup,lwdn &
   ,pressure &
   ,dn0 &
   )

!-----------------------------------------------------------------------------
! radcalc3: column driver for two-stream radiation code
! variables used within routine radcalc3:
! ==================================================================
! Variables in rrad3 parameter statement
!  mb               : maximum allowed number of bands [=8]
!  mg               : maximum allowed number of gases [=3]
!  mk               : maximum number of pseudobands allowed for any gas [=7]
!  ncog             : number of fit coefficients (omega and asym) [=5]
!  ncb              : number of fit coefficients (extinction) [=2]
!  npartob          : number of hydrometeor categories (including different habits)
!  npartg           : number of hydrometeor categories used for gc coefficients [=7]
!  nrad             : total number of radiation levels used (m1 - 1 + narad)
!  narad            : number of radiation levels added above model
!  nsolb            : active number of solar bands
!  nb               : active number of bands
!  ng               : active number of gases
!  jday             : julian day
!  solfac           : solar constant multiplier for variable E-S distance
!  ralcs (mb)       : rayleigh scattering integration constants
!  solar1 (mb)      : solar fluxes at top of atmosphere - corrected for ES distance
!  solar0 (mb)      : solar fluxes at top of atmosphere - uncorrected for ES distance
!  nuum (mb)        :    continuum flags
!  a0,a1,a2,a3 (mb) : Planck func fit coefficients
!  npsb (mg,mb)     : number of pseudo bands
!  trf (mg,mb)      : reference temperature for xp and wght coefficients
!  prf (mg,mb)      : reference pressure for xp and wght coefficients
!  ulim (mg,mb)     : upper bound on pathlength for gases
!  xp (mg,mk,mb)    : coefficient used in computing gaseous optical depth
!  alpha (mg,mk,mb) : pressure correction factor exponent
!  beta (mg,mk,mb)  : temperature correction factor exponent
!  wght (mg,mk,mb)  : pseudo band weight
!  exptabc (150)    : table of exponential func values
!  ocoef(ncog,mb,npartob)  : fit coefficients for hyd. single scatter.
!  bcoef(ncb,mb ,npartob)  : fit coefficients for hyd. extinction coefficient.
!  gcoef(ncog,mb,npartg)   : fit coefficients for hyd. asymmetry parameter.
!
! Input variables from model
!
!  m1               : number of vertical levels in model grid
!  ncat             : max number of hydrometeor categories [=7]
!  mcat             : actual number of hydrometeor categories [= 0, 1, or 7]
!  nhcat            : number of hydrometeor categories including ice habits [=15]
!  iswrtyp          : shortwave radiation parameterization selection flag
!  ilwrtyp          : longwave radiation parameterization selection flag
!  glat             : latitude
!  rtgt             : terrain-following coordinate metric factor
!  topt             : topography height
!  albedt           : surface albedo
!  cosz             : solar zenith angle
!  aodt             : total aerosol optical depth (band=3)
!  rlongup          : upward longwave radiation at surface (W/m^2)
!  rshort           : downward shortwave radiation at surface (W/m^2)
!  rlong            : downward longwave radiation at surface (W/m^2)
!  jnmb (ncat)      : microphysics category flag
!  dnfac (nhcat)    : factor for computing dn from emb
!  pwmasi (nhcat)   : inverse of mass power law exponent for hydrometeors
!  zm (m1)          : model physical heights of W points (m)
!  zt (m1)          : model physical heights of T points (m)
!  press (nzpmax)   : model pressure (Pa)
!  tair (nzpmax)    : model temperature (K)
!  rv (m1)          : model vapor mixing ratio (kg/kg)
!  dn0 (m1)         : model air density (kg/m^3)
!  fthrd (m1)       : theta_il tendency from radiation
!  jhcat (nzpmax,ncat)  : microphysics category array
!  cx (nzpmax,ncat) : hydrometeor number concentration (#/kg)
!  emb (nzpmax,ncat): hydrometeor mean mass (kg)
!
! Variables input from model scratch space (redefined internally on each call)
!
!  zml (nrad)       : physical heights of W points of all radiation levels (m)
!  ztl (nrad)       : physical heights of T points of all radiation levels (m)
!  dzl (nrad)       : delta-z (m) of all radiation levels
!  pl (nrad)        : pressure (Pa)
!  tl (nrad)        : temperature (K)
!  dl (nrad)        : air density of all radiation levels (kg/m^3)
!  rl (nrad)        : vapor density of all radiation levels (kg/m^3)
!  vp (nrad)        : vapor pressure (Pa)
!  o3l (nrad)       : stores the calculated ozone profile (g/m^3)
!  flxu (nrad)      : Total upwelling flux (W/m^2)
!  flxd (nrad)      : Total downwelling flux (W/m^2)
!  t (nrad)         : layer transmission func
!  r (nrad)         : layer reflection func
!  tc (nrad)        : cumulative optical depth
!  sigu (nrad)      : upwelling layer source func
!  sigd (nrad)      : downwelling layer source func
!  re (nrad)        : cumulative reflection func
!  vd (nrad)        : multi-scat. diffuse downwelling contributions
!                         from source func
!  td (nrad)        : inverse of cumulative transmission fnct
!  vu (nrad)        : multi-scat. diffuse upwelling contributions
!                         from source func
!  tg (nrad)        : gaseous optical depth
!  tcr (nrad)       : continuum/Rayleigh optical depth
!  src (nrad)       : Planck func source for each band
!  fu (nrad*6)      : upwelling fluxes for pseudo-bands (W/m^2)
!  fd (nrad*6)      : downwelling fluxes for pseudo-bands (W/m^2)
!  u (nrad*mg)      : path-length for gases (H_2O, CO_2, O_3)  (Pa)
!  tp (nrad*mb)     : optical depth of hydrometeors (m^-1)
!  omgp (nrad*mb)   : Single scatter albedo of hydrometeors
!  gp (nrad*mb)     : Asymmetry factor of hydrometeors
!
! Locally-defined variables
!
!  ngass (mg) : flags indicating if H20,CO2,O3 are active for solar wavelengths
!  ngast (mg) : flags indicating if H20,CO2,O3 are active for long wavelengths
!
! Additional variables used only within routine mclatchy:
! ==================================================================
! namax : maximum allowed number of added rad levels above model top[=10]
!                       used for oc and bc coefficients [=13]
! mcdat (33,9,6)    : Mclatchy sounding data (33 levels, 9 soundings, 6 vars)
! mclat (33,9,6)    : mcdat interpolated by season to latitude bands
! mcol (33,6)       : mclat interpolated to lat-lon of grid column
!
! Additional variables used only within routine cloud_opt:
! ==================================================================
!  ib .......... band number
!  dn .......... characteristic diameter (m)
!  oc .......... scattering albedo fit coefficients
!  bc .......... extinction fit coefficients
!  gc .......... asymmetery fit coefficients
!  kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                   numbers as a func of 15 microphysics category numbers
!
! Particle Numbers describe the following particles:
!
!     Harrington radiation code             RAMS microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   rain                        2.  rain
!  3:   pristine ice columns        3.  pristine ice columns
!  4:   pristine ice rosettes       4.  snow columns
!  5:   pristine ice plates         5.  aggregates
!  6:   snow columns                6.  graupel
!  7:   snow rosettes               7.  hail
!  8:   snow plates                 8.  pristine ice hexagonal plates
!  9:   aggregates columns          9.  pristine ice dendrites
!  10:  aggregates rosettes        10.  pristine ice needles
!  11:  aggregates plates          11.  pristine ice rosettes
!  12:  graupel                    12.  snow hexagonal plates
!  13:  hail                       13.  snow dendrites
!                                  14.  snow needles
!                                  15.  snow rosettes
!
! for the asymmetery parameter, since we only have spherical
! particles, there are only 7 particle types...
!  1:   cloud drops
!  2:   rain
!  3:   pristine ice
!  4:   snow
!  5:   aggregates
!  6:   graupel
!  7:   hail
!---------------------------------------------------------------------------

use rconstants
use rrad3
use micphys
!use node_mod

implicit none

integer m1,maxnzp,mcat,ngrid
integer :: iswrtyp,ilwrtyp
integer i,j,k
integer, save :: ncall = 0,nradmax
integer, save :: ngass(mg)=(/1, 1, 1/),ngast(mg)=(/1, 1, 1/)
!     one can choose the gases of importance here,
!       ngas = 1    gas active
!            = 0    gas not active
!
!       ngas(1) = H2O
!       ngas(2) = CO2
!       ngas(3) =  O3

real, save :: eps=1.e-15
real :: glat,rtgt,topt,cosz,aodt,albedt,rlongup,rshort,rlong
real :: zm(m1),zt(m1),dn0(m1),rv(m1),fthrd(m1)
real :: pressure(m1)
real :: bext(m1),swup(m1),swdn(m1),lwup(m1),lwdn(m1)

real, allocatable, save, dimension(:)     :: zml,ztl,dzl,pl,tl,dl,rl,o3l  &
                                      ,vp,flxus,flxds,tg,tcr,src,t,r,tc  &
                                      ,flxul,flxdl  &
                                      ,sigu,sigd,re,vd,td,vu  &
                                      ,u,fu,fd,tp,omgp,gp

real :: exner(m1) !non-dimensional pressure

!Saleeby(2011):Variables for radiatively active aerosols
real :: relh(m1)
real, external :: rslf

if (ncall == 0) then
   ncall = 1
   nradmax = maxnzp + namax
   allocate(zml  (nradmax) ,ztl  (nradmax) ,dzl  (nradmax) ,pl (nradmax)  &
           ,tl   (nradmax) ,dl   (nradmax) ,rl   (nradmax) ,o3l(nradmax)  &
           ,vp   (nradmax) ,flxus(nradmax) ,flxds(nradmax) ,tg (nradmax)  &
           ,flxul(nradmax) ,flxdl(nradmax)                                &
           ,tcr  (nradmax) ,src  (nradmax) ,t    (nradmax) ,r  (nradmax)  &
           ,tc   (nradmax) ,sigu (nradmax) ,sigd (nradmax) ,re (nradmax)  &
           ,vd   (nradmax) ,td   (nradmax) ,vu   (nradmax))
   allocate(u(nradmax*mg),fu(nradmax*6),fd(nradmax*6))
   allocate(tp(nradmax*mb),omgp(nradmax*mb),gp(nradmax*mb))
   tg=0.
endif

nrad = m1 - 1 + narad

! rlongup used to set tl(1): stephan*tl^4=rlongup
 CALL mclatchy (3,m1  &
   ,glat,rtgt,topt,rlongup  &
   ,zm,zt,press,tair,dn0,rv,zml,ztl,pl,tl,dl,rl,o3l,dzl &
   )

! calculate non-dimensional pressure
do k=1,m1
  exner(k) = (press(k)*p00i)**rocp
enddo

! zero out scratch arrays
 CALL azero (nrad*mg,u)
 CALL azero (nrad*6,fu)
 CALL azero (nrad*6,fd)
 CALL azero (nrad*mb,tp)
 CALL azero (nrad*mb,omgp)
 CALL azero (nrad*mb,gp)
 CALL azero (nrad   ,tg)

! Saleeby(2011): Aerosol radiative impacts section
! This must be run before cloud_opt
! Only run this for level=3 microphysics
if(iaerorad==1 .and. level .ne. 4) then
 !Only do levels 1 through m1-1. If there are no additional radiation levels,
 !then pl and tl only get inititalized (mclatchy call with iaction == 3) up 
 !to level m1-1. Also, aerorad always computes from level 2 to m1-1.
 do k=1,m1-1
   relh(k) = rv(k)/rslf(pl(k),tl(k))
 enddo
 CALL aerorad (i,j,mb,nb,nrad,m1,dzl,relh,tp,omgp,gp,dn0,aodt)
endif

! Compute hydrometeor optical properties
 CALL cloud_opt (mb,nb,nrad,m1,mcat,dzl  &
   ,dn0,tp,omgp,gp &
   ,ocoef,bcoef,gcoef,ncog,ncb,npartob,npartg)

! Sum up attenutation by aerosols and hydrometeors
 CALL sum_opt (mb,nrad,nb,m1,tp,omgp,gp,bext,dzl)

! Get the path lengths for the various gases...
 CALL path_lengths (nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

do k = 1,nrad
   if (rl(k) <   0. .or.  &
       dl(k) <   0. .or.  &
       pl(k) <   0. .or.  &
      o3l(k) <   0.) then
      print*, 'Negative value of density, vapor, pressure, or ozone'
      print*, 'when calling Harrington radiation'
!      print*, 'at k,i,j = ',k,i+mi0(ngrid),j+mj0(ngrid)
      print*, 'at k,i,j = ',k,i,j
      print*, 'ngrid=',ngrid
      print*, 'stopping model'
      print*, 'rad: rl(k), dl(k), pl(k), o3l(k)'
      print*, rv(k), dl(k), pl(k), o3l(k)
      stop
   endif
enddo
do k = 1,nrad
   if (tl(k) < 160.) then
      print*, 'Temperature too low when calling Harrington radiation' 
!      print*, 'at k,i,j = ',k,i+mi0(ngrid),j+mj0(ngrid)
      print*, 'at k,i,j = ',k,i,j
      print*, 'ngrid,tempk=',ngrid,tl(k)
      !stop
   endif
enddo   

! call shortwave and longwave schemes...

if (iswrtyp == 3 .and. cosz > 0.03) then
   CALL azero2 (nrad,flxus,flxds)

   CALL swrad (nrad,nb,nsolb,npsb,       & !counters
      u,pl,tl,dzl,                       & !model variables
      xp,alpha,beta,wght,prf,trf,ralcs,  & !band specifics
      solar1,ngass,                      & !coefficients
      albedt,cosz,                       & !boundaries
      tp,omgp,gp,fu,fd,                  & !optical properties  
      flxus,flxds,ulim)                    !sw fluxes

   rshort = flxds(1)

   do k = 2,m1-1
      !divide by exner to get potential temp heating rate
      fthrd(k) = fthrd(k)  &
         + (flxds(k) - flxds(k-1) + flxus(k-1) - flxus(k)) &
            / (dl(k) * dzl(k) * cp * exner(k))
      swup(k) = flxus(k)
      swdn(k) = flxds(k)
    enddo
    !lower and upper boundary conditions on swup and swdn
    swup(1) = flxus(1)
    swup(m1) = flxus(nrad) ! use the top radiation value rather than m1 value
    swdn(1) = flxds(1)
    swdn(m1) = flxds(nrad) ! use the top radiation value rather than m1 value

else

   swup   = 0.
   swdn   = 0.
   rshort = 0.

endif

if (ilwrtyp == 3) then
   CALL azero2 (nrad,flxul,flxdl)

   CALL lwrad (i,j,nrad,nb,nsolb,npsb,nuum,   & !counters
      u,pl,tl,vp,                             & !model variables
      xp,alpha,beta,wght,prf,trf,ralcs,       & !band specifics
      a0,a1,a2,a3,                            & !coefficients
      exptabc,ngast,                          & !boundaries
      tp,omgp,gp,fu,fd,flxul,flxdl,ulim)        !fluxes, output

   !Set rlong to surface level downward longwave flux.
   rlong = flxdl(1)

   !Make lowest level upward longwave flux equal to rlongup
   !produced from land surface models (LEAF,SiB).
   flxul(1) = rlongup

   do k = 2,m1-1
      !divide by exner to get potential temp heating rate
      fthrd(k) = fthrd(k)  &
         + (flxdl(k) - flxdl(k-1) + flxul(k-1) - flxul(k)) &
            / (dl(k) * dzl(k) * cp * exner(k))
      lwup(k) = flxul(k)
      lwdn(k) = flxdl(k)
   enddo
   !lower and upper boundary conditions on lwup and lwdn
   lwup(1) = flxul(1)
   lwup(m1) = flxul(nrad) ! use the top radiation value rather than m1 value
   lwdn(1) = flxdl(1)
   lwdn(m1) = flxdl(nrad) ! use the top radiation value rather than m1 value

endif

return
END SUBROUTINE radcalc3

!##############################################################################
Subroutine cloud_opt (mb,nb,nrad,m1,mcat,dzl  &
   ,dn0,tp,omgp,gp,oc,bc,gc,ncog,ncb,npartob,npartg)

! computing properties of spherical liquid water and irregular ice
! using fits to adt theory
!
! ib .......... band number
! mb .......... maximum number of bands
! nb .......... total number of bands
! m1 .......... number of vertical levels
! dzl .......... delta z in each level (m)
! dn .......... characteristic diameter (m)
! emb ......... mean hydrometeor mass (kg)
! cx .......... hydrometeor concentration (#/kg)
! tp .......... optical depth
! omgp ........ scattering albedo
! gp .......... asymmetry parameter
! oc .......... scattering albedo fit coefficients
! bc .......... extinction fit coefficients
! gc .......... asymmetry fit coefficients
! ncog ........ number of fit coefficients (omega and asym)
! ncb ......... number of fit coefficients (extinction)
! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                 numbers as a func of 15 microphysics category numbers
! dn0 ......... model air density (kg/m^3)
! dnfac ....... factor for computing dn from emb
! pwmasi ...... inverse of power used in mass power law
! npartob ..... number of hydrometeor categories (including different habits)
!                 used for oc and bc coefficients
! npartg ...... number of hydrometeor categories used for gc coefficients
!
! Saleeby(2008): would need to modify dnmin,dnmax,kradcat for drizzle

use micphys

implicit none

integer mb,nb,ib,nrad,m1,ncog,ncb,krc,npartob,npartg
integer icat,mcat,k,ihcat

integer kradcat(15)
real dzl(nrad),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),dn0(m1)
real oc(ncog,mb,npartob),bc(ncb,mb,npartob),gc(ncog,mb,npartg)
real ext,om,gg,dn,gnu_tab

real dnmin(7),dnmax(7)
data dnmin /   1.,   10.,   1.,  125.,   10.,   10.,   10./
data dnmax /1000.,10000., 125.,10000.,10000.,10000.,10000./

data kradcat/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/

if (level <=3) then
   do icat = 1,mcat 
      if (jnmb(icat) .gt. 0) then
 
         do k = 2,m1-1

            if (cx(k,icat) .gt. 1.e-9) then
               ihcat = jhcat(k,icat)
               krc = kradcat(ihcat)
               dn = dnfac(ihcat) * emb(k,icat) ** pwmasi(ihcat) * 1.e6

               !Modification by Adele Igel (Oct 2020).
               !Use an effective Dn here constrainted by effective diameter.
               !Gnu from microphysics can be different from Gnu allowed in 
               !radiation due to lookup tables available only for Gnu=1,2.
               !The two lines below help mitigate this difference.
               !See Harrington's thesis for more detail.
               gnu_tab = real(max(1,min(2,nint(gnu(icat)))))
               dn = dn * (gnu(icat)+2.) / (gnu_tab + 2.)

               dn = max(dnmin(icat),min(dnmax(icat),dn))

               do ib = 1,nb

                  ext = cx(k,icat) * dn0(k) * dzl(k)  &
                     * bc(1,ib,krc) * dn ** bc(2,ib,krc)
                  om = oc(1,ib,krc)  &
                     + oc(2,ib,krc) * exp(oc(3,ib,krc) * dn)  &
                     + oc(4,ib,krc) * exp(oc(5,ib,krc) * dn)
                  gg = gc(1,ib,icat)  &
                     + gc(2,ib,icat) * exp(gc(3,ib,icat) * dn)  &
                     + gc(4,ib,icat) * exp(gc(5,ib,icat) * dn)

                  tp(k,ib)   = tp(k,ib)   + ext
                  omgp(k,ib) = omgp(k,ib) + om * ext
                  gp(k,ib)   = gp(k,ib)   + gg * om * ext
  
               enddo

            endif 
         enddo
      endif
   enddo
endif

return
END SUBROUTINE cloud_opt

!##############################################################################
Subroutine sum_opt (mb,nrad,nb,m1,tp,omgp,gp,bext,dzl)

implicit none

integer nb,ib,m1,k,mb,nrad
real tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),bext(m1),dzl(m1)

! Combine the optical properties....

do ib = 1,nb
   do k = 2,m1-1

      if (tp(k,ib) .gt. 0.0) then
         gp(k,ib) = gp(k,ib) / omgp(k,ib)
         !Saleeby(2010): Need this min/max func to prevent 'gp'
         ! from being unphysical. Needed for RCE simulations.
         gp(k,ib) = MAX(MIN(gp(k,ib),1.),0.)
         omgp(k,ib) = omgp(k,ib) / tp(k,ib)
         !Saleeby(2010): Need this min/max func to prevent 'omgp'
         ! from being unphysical. Needed for RCE simulations.
         omgp(k,ib) = MAX(MIN(omgp(k,ib),1.),0.)
      else
         omgp(k,ib) = 0.0
         gp(k,ib) = 0.0
      endif

      !Check for validity of opt values before calling radiation
      if (tp(k,ib) .lt. 0) then
         print*, 'tp(k,ib) less than zero for k,ib = ',k,ib
         print*, 'tp(k,ib) = ',tp(k,ib)
         stop 'opt1'
      endif
      if (omgp(k,ib) .lt. 0. .or. omgp(k,ib) .gt. 1.) then
         print*, 'omgp(k,ib) out of range [0,1] for k,ib = ',k,ib
         print*, 'omgp(k,ib) = ',omgp(k,ib)
         stop 'opt2'
      endif
      if (gp(k,ib) .lt. 0. .or. gp(k,ib) .gt. 1.) then
         print*, 'gp(k,ib) out of range [0,1] for k,ib = ',k,ib
         print*, 'gp(k,ib) = ',gp(k,ib)
         stop 'opt3'
      endif

   enddo
enddo

! Calculating visual range (km) using the Koschmeider equation
! Units: tp(k,ib) : optical thickness in level k, band ib (dimensionless)
!        dzl(k) : [m]
!        final bext(k) : [km]
! Consider only over band #3 (245-700 nm)
 do k = 2,m1-1
  bext(k) = 0.
  bext(k) = tp(k,3) / dzl(k) !Compute extinction coefficient
  !Prevent infinite visibility: constrain max vis to 1000km
  if(bext(k) .lt. 3.912e-6) bext(k) = 3.912e-6
  bext(k) = 3.912/bext(k)/1000.
 enddo

return
END SUBROUTINE sum_opt

!##############################################################################
Subroutine path_lengths (nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

! Get the path lengths for the various gases...

implicit none

integer :: nrad
real, dimension(nrad) :: rl,dzl,dl,o3l,vp,pl
real, dimension(nrad,3) :: u
real :: rvk0,rvk1,dzl9,rmix,eps
integer :: k

u(1,1) = .5 * (rl(2) + rl(1)) * 9.81 * dzl(1)
u(1,2) = .5 * (dl(2) + dl(1)) * .45575e-3 * 9.81 * dzl(1)
u(1,3) = o3l(1) * 9.81 * dzl(1)

rvk0 = rl(1)
do k = 2,nrad
   rvk1 = (rl(k) + 1.e-6)
   dzl9 = 9.81 * dzl(k)
   rmix = rvk1 / dl(k)
   vp(k) = pl(k) * rmix / (.622 + rmix)
   u(k,1) = (rvk1 - rvk0) / (log(rvk1 / rvk0) + eps) * dzl9
   u(k,2) = (dl(k) - dl(k-1)) / (log(dl(k) / dl(k-1)) + eps)  &
       * dzl9 * 0.45575e-3
   u(k,3) = 0.5 * dzl9 * (o3l(k) + o3l(k-1))
   rvk0 = rvk1
enddo

return
END SUBROUTINE path_lengths

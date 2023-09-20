!##############################################################################
Subroutine init_aero_custom ()

use mem_grid, only:ngrid
use mem_basic, only:basic_g
use mem_micro, only:micro_g
use micphys, only:ipris,iifn,level,iaerosol,idust,isalt,iabcarb,initcustom
use mem_other, only:other_g
use node_mod, only:mzp,mxp,myp,ia,iz,jz,jz

implicit none

   CALL init_thermo (mzp,mxp,myp     &
      ,basic_g(ngrid)%thp   (1,1,1)  &
      ,basic_g(ngrid)%theta (1,1,1)  &
      ,basic_g(ngrid)%rtp   (1,1,1)  &
      ,basic_g(ngrid)%rv    (1,1,1))

   if(iaerosol > 0) CALL init_ccn1 (mzp,mxp,myp    &
      ,micro_g(ngrid)%cn1np (1,1,1)  &
      ,micro_g(ngrid)%cn1mp (1,1,1)  &
      ,basic_g(ngrid)%dn0   (1,1,1),ngrid)

   if(iaerosol > 0) CALL init_ccn2 (mzp,mxp,myp   &
      ,micro_g(ngrid)%cn2np (1,1,1)  &
      ,micro_g(ngrid)%cn2mp (1,1,1)  &
      ,basic_g(ngrid)%dn0   (1,1,1),ngrid)

   if(idust    > 0) CALL init_dust (mzp,mxp,myp   &
      ,micro_g(ngrid)%md1np (1,1,1)  &
      ,micro_g(ngrid)%md2np (1,1,1)  &
      ,micro_g(ngrid)%md1mp (1,1,1)  &
      ,micro_g(ngrid)%md2mp (1,1,1),ngrid)

   if(iabcarb  > 0) CALL init_absorbing_carbon (mzp,mxp,myp   &
      ,micro_g(ngrid)%abc1np (1,1,1)  &
      ,micro_g(ngrid)%abc2np (1,1,1)  &
      ,micro_g(ngrid)%abc1mp (1,1,1)  &
      ,micro_g(ngrid)%abc2mp (1,1,1),ngrid)

   if(isalt    > 0) CALL init_salt (mzp,mxp,myp   &
      ,micro_g(ngrid)%salt_film_np (1,1,1)  &
      ,micro_g(ngrid)%salt_jet_np  (1,1,1)  &
      ,micro_g(ngrid)%salt_spum_np (1,1,1)  &
      ,micro_g(ngrid)%salt_film_mp (1,1,1)  &
      ,micro_g(ngrid)%salt_jet_mp  (1,1,1)  &
      ,micro_g(ngrid)%salt_spum_mp (1,1,1),ngrid)

   if(ipris >= 5 .and. (iifn==1.or.iifn==2)) then
       CALL init_ifn (mzp,mxp,myp    &
        ,micro_g(ngrid)%cifnp (1,1,1)  &
        ,basic_g(ngrid)%dn0   (1,1,1),ngrid)
   endif

   if(idust == 2) then
    CALL init_dustfrac_read ()
    CALL init_dustfrac (mxp,myp)
   endif

   if(initcustom == 1) CALL init_custom (mzp,mxp,myp)

return
END SUBROUTINE init_aero_custom

!##############################################################################
Subroutine init_rtgt_dz ()

use mem_grid, only:grid_g,ztop,zm,zt,dzt,dzm,ngrid
use node_mod, only:mzp,ia,iz,ja,jz

implicit none

integer :: i,j,k

do j = ja,jz
 do i = ia,iz

  !Compute some grid information that RAMS would normally do at initialization
  !topt = topography(m) at scalar horizontal points (ie. T)
  !topm = topography(m) at momentum horizontal points (ie. U,V)
  !ztop = top prognostic vertical level altitude AGL (m)
  !rtgt = vertical level topographic adjustment factor on horizontal T grid
  !rtgm = vertical level topographic adjustment factor on horizontal M grid
  !dzt  = inverse of vertical grid spacing on T grid
  !dzm  = inverse of vertical grid spacing on M grid
  ztop = zm(mzp-1)
  grid_g(ngrid)%rtgt(i,j) = 1. - grid_g(ngrid)%topt(i,j) / ztop
  do k = 2,mzp
   dzt(k) = 1. / (zm(k) - zm(k-1))
  enddo
  dzt(1) = dzt(2) * dzt(2) / dzt(3)
  do k = 1,mzp-1
   dzm(k) = 1. / (zt(k+1) - zt(k))
  enddo
  dzm(mzp) = dzm(mzp-1) * dzm(mzp-1) / dzm(mzp-1)

 enddo
enddo

return
END SUBROUTINE init_rtgt_dz

!##############################################################################
Subroutine init_dustfrac (n2,n3)

use micphys
use mem_micro
use mem_grid

implicit none

integer :: n2,n3,isource,jsource,i,j

do j = 1,n3
   do i = 1,n2

  !Idealized lofting with every grid cell able to loft with 100% fraction
  if(idustloft==0) then
   if(iprntstmt>=1 .and. print_msg) &
    print*,' '
    print*,'Setting idealized dust erodible fraction of 1.0'
   micro_g(ngrid)%dustfrac(i,j) = 1.0

  !Ginoux 1.0 degree dust erodible fraction lofting database
  elseif(idustloft==1) then
   if(iprntstmt>=1 .and. print_msg) &
    print*,' '
    print*,'Setting Ginoux 1.0 degree dust erodible fraction'
   isource = nint(grid_g(ngrid)%glon(i,j)-(-179.50))+1
   jsource = nint(grid_g(ngrid)%glat(i,j)-(-089.50))+1
   micro_g(ngrid)%dustfrac(i,j) = ginoux_1p00_deg_dustfrac(isource,jsource)

  endif

 enddo
enddo

return
END SUBROUTINE init_dustfrac

!##############################################################################
Subroutine init_dustfrac_read ()

use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

character(len=strl1) :: cname
real, allocatable, dimension(:) :: lat_source,lon_source
integer :: ny_source,nx_source,i,j

!Standard DUSTFILE name from RAMSIN. This is used when not needing to
!create grid-specific dust file or reading in pre-created grid-specific
!dust file.
cname=dustfile(1:len_trim(dustfile))

!Ginoux 1.0 degree dust erodible fraction
if(idustloft==1) then
 nx_source=360
 ny_source=180
 allocate(lat_source(ny_source))
 allocate(lon_source(nx_source))
 open(91,file=cname,form='formatted',status='old')
 if(iprntstmt>=1 .and. print_msg) &
  print*,' '
  print*,'GINOUX 1-DEGREE DUST SOURCE ERODIBLE FRACTION DATA READ IN'
 do i = 1,nx_source
  do j = 1,ny_source
    read(91,*)lat_source(j),lon_source(i),ginoux_1p00_deg_dustfrac(i,j)
  enddo
 enddo
 close(91)
 deallocate(lat_source,lon_source) 

endif

return
END SUBROUTINE init_dustfrac_read

!##############################################################################
Subroutine init_custom (n1,n2,n3)

use mem_grid, only:ngrid
use mem_micro, only:micro_g
use micphys

implicit none

integer :: n1,n2,n3,i,j,k

! Initialize thermo and moisture for RAMS micro from parent model
! We assume there is no condensate or supersaturation to start at time zero
! such that theta-il = theta and rtp = rv.

do j = 1,n3
   do i = 1,n2
      do k = 34,34

       !Setup for 1mm rain drop diameters
       micro_g(ngrid)%rrp(k,i,j) = 5.24e-3   ! kg/kg rain water
       micro_g(ngrid)%crp(k,i,j) = 10000.    ! 10000 /kg

      enddo
   enddo
enddo

return
END SUBROUTINE init_custom

!##############################################################################
Subroutine init_thermo (n1,n2,n3,thp,theta,rtp,rv)

use micphys, only:itheta,irv

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: thp,theta,rtp,rv

! Initialize thermo and moisture for RAMS micro from parent model
! We assume there is no condensate or supersaturation to start at time zero
! such that theta-il = theta and rtp = rv.

do j = 1,n3
   do i = 1,n2
      do k = 1,n1

        !Initialize THP assuming THETA is the input thermodynamic variable.
        if(itheta == 1) thp(k,i,j) = theta(k,i,j)
        !Initialize RTP assuming RV (water vapor) is the input moisture variable.
        if(irv    == 1) rtp(k,i,j) = rv(k,i,j)

      enddo
   enddo
enddo

return
END SUBROUTINE init_thermo

!##############################################################################
Subroutine aerosols ()

use mem_basic
use mem_micro
use mem_grid
use mem_leaf
use micphys
use node_mod, only:mzp,mxp,myp,ia,iz,ja,jz

implicit none

!Run the SEASALT and DUST Source model before the call to Micro
if(idust==2) then
  CALL dust_sources (mzp,mxp,myp,ia,iz,ja,jz      &
                    ,grid_g(ngrid)%rtgt           &
                    ,grid_g(ngrid)%glat           &
                    ,grid_g(ngrid)%glon           &
                    ,basic_g(ngrid)%up            &
                    ,basic_g(ngrid)%vp            &
                    ,basic_g(ngrid)%dn0           &
                    ,micro_g(ngrid)%md1np         &
                    ,micro_g(ngrid)%md2np         &
                    ,micro_g(ngrid)%md1mp         &
                    ,micro_g(ngrid)%md2mp         &
                    ,micro_g(ngrid)%dustfrac      &
                    ,leaf_g(ngrid)%soil_water     &
                    ,leaf_g(ngrid)%patch_area     &
                    ,leaf_g(ngrid)%leaf_class     &
                    ,leaf_g(ngrid)%soil_text      &
                   !,leaf_g(ngrid)%veg_rough      &
                    )
endif

if(isalt==2) then
  CALL salt_sources (mzp,mxp,myp,ia,iz,ja,jz      &
                    ,grid_g(ngrid)%rtgt           &
                    ,basic_g(ngrid)%up            &
                    ,basic_g(ngrid)%vp            &
                    ,basic_g(ngrid)%dn0           &
                    ,micro_g(ngrid)%salt_film_np  &
                    ,micro_g(ngrid)%salt_jet_np   &
                    ,micro_g(ngrid)%salt_spum_np  &
                    ,micro_g(ngrid)%salt_film_mp  &
                    ,micro_g(ngrid)%salt_jet_mp   &
                    ,micro_g(ngrid)%salt_spum_mp  &
                    ,leaf_g(ngrid)%patch_area     &
                    ,leaf_g(ngrid)%leaf_class     &
                    )
endif


return
END SUBROUTINE aerosols

!##############################################################################
Subroutine aerosol_init ()

use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

real :: weightfac

! Set aerosol density depending on chemistry and soluble fraction
! Pure quantity densities (kg/m3) are:
! NH42S04 = 1769. (ammonium sulfate)
! Clay Dust (smaller) = 2500.
! Silt Dust (larger) = 2650.
! NaCl = 2165. (sodium chloride)

! Set Aerosol density (kg/m3) based on weighted mixture of soluble
! and insoluble material. Assume insoluble core to be like that of 
! silt dust with density = 2650 kg/m3, except for acat=3 which is
! already set to small sized clay dust
! Also set vanthoff factors for given chemistry

if(iprntstmt>=1 .and. print_msg) print*,''
if(iprntstmt>=1 .and. print_msg) print*,'Setting up default aerosol densities:'

do acat=1,aerocat
 
 aero_rhosol(acat)   = 0.0
 aero_vanthoff(acat) = 0.0

 if(acat==3) then !if small dust (clay)
   weightfac = 2500. * (1.0-aero_epsilon(acat)) !clay dust core
 else !all other
   weightfac = 2650. * (1.0-aero_epsilon(acat)) !silt dust core
 endif

 if(iaero_chem(acat)==1) then !NH42S04
   aero_rhosol(acat) = 1769. * aero_epsilon(acat) + weightfac
   aero_vanthoff(acat) = 3
 elseif(iaero_chem(acat)==2) then !NaCl
   aero_rhosol(acat) = 2165. * aero_epsilon(acat) + weightfac
   aero_vanthoff(acat) = 2
 endif

 if(iprntstmt>=1 .and. print_msg) print*,'acat,rg,rho,i:',acat &
     ,aero_medrad(acat),aero_rhosol(acat),aero_vanthoff(acat)

enddo

if(iprntstmt>=1 .and. print_msg) print*,''

return
END SUBROUTINE aerosol_init

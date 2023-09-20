!##############################################################################
Subroutine siminfo (fname,fln)

use mem_grid
use node_mod, only:mmzp,mzp
use io_params, only:frqstate

implicit none

integer :: i
character(len=*) :: fname,fln
character(len=64) :: fullname
logical :: file_exists

fullname=trim(fln)//trim(fname)

!READ SIMULATION INFO RELATED TO PARENT MODEL.
!IF COUPLED TO OTHER MODEL, THIS INFO WOULD BE PASSED TO THE DRIVER
!AND A FILE WOULD NOT BE READ IN.
inquire(file=trim(fullname),EXIST=file_exists)
if(file_exists)then
 print*,"Opening parent model info file"
 open(50,file=trim(fullname),status='old')
 read(50,*)iyear1  !Start year of the simulation (YYYY)
 read(50,*)imonth1 !Start month of the simulation
 read(50,*)idate1  !Start date of the simulation
 read(50,*)itime1  !Start time of the simulation (UTC)
 read(50,*)
 read(50,*)dtlt    !Current grid timestep
 read(50,*)time    !Current model run time (seconds)
 read(50,*)ngrid   !Cureent model grid
 read(50,*)mzp     !Current grid number of vertical levels
 read(50,*)
 read(50,*)ngrids  !Total number of grids including any nested grids
 do i=1,ngrids     
 read(50,*)mmzp(i) !Number vertical levels on T-grid for each grid
 enddo
 do i=1,ngrids
 read(50,*)frqstate(i) !Frequency(/sec) writing output for each grid
 enddo
 close(50)
else
 print*," "
 print*,"Input file does not exist: ",trim(fullname)
 print*,"Using default idealized input from main.f90"
 print*," "
endif

return
END SUBROUTINE siminfo


!##############################################################################
Subroutine iofil (fname,x,ix,wf,fln,istatic)

use mem_grid, only:time,iprntstatic

implicit none

character(len=*) :: fname,fln
character(len=64) :: fullname
integer :: ix                    !number vertical levels
integer :: wf                    !read or write flag
integer :: istatic               !flag to designate changing variable
real, dimension(ix) :: x         !variable to read/write
logical :: file_exists
integer k,l

fullname=trim(fln)//trim(fname)

!READ SCM DATA
if(wf == 1) then 
  inquire(file=trim(fullname),EXIST=file_exists)
  if(file_exists)then
   open(50,file=trim(fullname),status='old')
   read(50,*)
   do k=1,ix
     read(50,*) l,x(k)
   enddo
   close(50)
  else
   print*,"Input file does not exist: ",trim(fullname)
  endif
!WRITE SCM DATA
elseif (wf == 2) then
 if(iprntstatic==1 .or. (iprntstatic==0 .and. istatic==1))then
  open(50,file=trim(fullname),status='unknown',access='append')
  write(50,*)'Elapsed time',time
  do k=1,ix
    write(50,'(i5,e26.18)') k,x(k)
  enddo
  close(50)
 endif
endif

return
END SUBROUTINE iofil

!##############################################################################
Subroutine readwrite_scm (wf,fln,istp)

use mem_other
use mem_basic
use mem_micro
use mem_grid
use mem_radiate
use mem_leaf
use micphys
use rconstants
use node_mod, only:mzp

implicit none
character(len=*) :: fln
integer :: m1,k,i,j,wf,ng,istp

if(iprntstmt>=1 .and. print_msg) then
 if(wf==1)then
  write(*,*)
  write(*,*)'***************************************************************'
  write(*,*) 'Reading in variables at time:',time
  write(*,*)'***************************************************************'
  write(*,*)
 elseif(wf==2)then
  write(*,*) 'Writing out variables at time:',time
 endif
endif

ng = ngrid
m1 = mzp
i  = 1
j  = 1

!Air pressure here since RAMS uses Exner(PI)
call iofil('pressure.txt',other_g(ng)%pressure(:),m1,wf,fln,0)

!Land surface variables read in here if available. Defaults in main.
!These variables do not change.
call iofil('soil_water.txt',leaf_g(ng)%soil_water(nzg,i,j,:),npatch,wf,fln,0)
call iofil('patch_area.txt',leaf_g(ng)%patch_area(i,j,:),npatch,wf,fln,0)
call iofil('leaf_class.txt',leaf_g(ng)%leaf_class(i,j,:),npatch,wf,fln,0)
call iofil('soil_text.txt' ,leaf_g(ng)%soil_text(nzg,i,j,:),npatch,wf,fln,0)

!Model grid point information (not changing)
call iofil('zm.txt',zm,m1,wf,fln,0)
call iofil('zt.txt',zt,m1,wf,fln,0)
call iofil('glat.txt',grid_g(ng)%glat(i,j),1,wf,fln,0)
call iofil('glon.txt',grid_g(ng)%glon(i,j),1,wf,fln,0)
call iofil('topt.txt',grid_g(ng)%topt(i,j),1,wf,fln,0)

!Basic dynamic input fields (not updated)
call iofil('up.txt'   ,basic_g(ng)%up   (:,i,j),m1,wf,fln,0)
call iofil('vp.txt'   ,basic_g(ng)%vp   (:,i,j),m1,wf,fln,0)
call iofil('wp.txt'   ,basic_g(ng)%wp   (:,i,j),m1,wf,fln,0)
call iofil('pp.txt'   ,basic_g(ng)%pp   (:,i,j),m1,wf,fln,0)
call iofil('dn0.txt'  ,basic_g(ng)%dn0  (:,i,j),m1,wf,fln,0)
call iofil('pi0.txt'  ,basic_g(ng)%pi0  (:,i,j),m1,wf,fln,0)

!Thermodynamic fields (these are updated)
call iofil('thp.txt'  ,basic_g(ng)%thp  (:,i,j),m1,wf,fln,1)
call iofil('theta.txt',basic_g(ng)%theta(:,i,j),m1,wf,fln,1)
call iofil('rtp.txt'  ,basic_g(ng)%rtp  (:,i,j),m1,wf,fln,1)
call iofil('rv.txt'   ,basic_g(ng)%rv   (:,i,j),m1,wf,fln,1)

!Radiation variables I/O if radiation turned on
if(iswrtyp == 3 .or. ilwrtyp == 3) then
 !These are computed in radation (yes updated)
 call iofil('fthrd.txt'  ,radiate_g(ng)%fthrd  (:,i,j),m1,wf,fln,1)
 call iofil('fthrdp.txt' ,radiate_g(ng)%fthrdp (:,i,j),m1,wf,fln,1)
 call iofil('bext.txt'   ,radiate_g(ng)%bext   (:,i,j),m1,wf,fln,1)
 call iofil('swdn.txt'   ,radiate_g(ng)%swdn   (:,i,j),m1,wf,fln,1)
 call iofil('swup.txt'   ,radiate_g(ng)%swup   (:,i,j),m1,wf,fln,1)
 call iofil('lwdn.txt'   ,radiate_g(ng)%lwdn   (:,i,j),m1,wf,fln,1)
 call iofil('lwup.txt'   ,radiate_g(ng)%lwup   (:,i,j),m1,wf,fln,1)
 call iofil('rshort.txt' ,radiate_g(ng)%rshort (i,j),1,wf,fln,1)
 call iofil('rlong.txt'  ,radiate_g(ng)%rlong  (i,j),1,wf,fln,1)
 call iofil('cosz.txt'   ,radiate_g(ng)%cosz   (i,j),1,wf,fln,1)
 call iofil('aodt.txt'   ,radiate_g(ng)%aodt   (i,j),1,wf,fln,1)
 !These are computed in land surface model (not updated)
 call iofil('rlongup.txt',radiate_g(ng)%rlongup(i,j),1,wf,fln,0)
 call iofil('albedt.txt' ,radiate_g(ng)%albedt (i,j),1,wf,fln,0)
endif

!CCN and GCCN
if(iaerosol > 0) then
 call iofil('cn1np.txt',micro_g(ng)%cn1np(:,i,j),m1,wf,fln,1)
 call iofil('cn1mp.txt',micro_g(ng)%cn1mp(:,i,j),m1,wf,fln,1)
 call iofil('cn2np.txt',micro_g(ng)%cn2np(:,i,j),m1,wf,fln,1)
 call iofil('cn2mp.txt',micro_g(ng)%cn2mp(:,i,j),m1,wf,fln,1)
endif
!Dust aerosols
if(idust > 0) then
 call iofil('md1np.txt',micro_g(ng)%md1np(:,i,j),m1,wf,fln,1)
 call iofil('md2np.txt',micro_g(ng)%md2np(:,i,j),m1,wf,fln,1)
 call iofil('md1mp.txt',micro_g(ng)%md1mp(:,i,j),m1,wf,fln,1)
 call iofil('md2mp.txt',micro_g(ng)%md2mp(:,i,j),m1,wf,fln,1)
 if(idust == 2) call iofil('dustfrac.txt',micro_g(ng)%dustfrac(i,j),1,wf,fln,1)
endif
!Absorbing carbon aerosols
if(iabcarb > 0) then
 call iofil('abc1np.txt',micro_g(ng)%abc1np(:,i,j),m1,wf,fln,1)
 call iofil('abc2np.txt',micro_g(ng)%abc2np(:,i,j),m1,wf,fln,1)
 call iofil('abc1mp.txt',micro_g(ng)%abc1mp(:,i,j),m1,wf,fln,1)
 call iofil('abc2mp.txt',micro_g(ng)%abc2mp(:,i,j),m1,wf,fln,1)
endif
!Sea salt aerosols
if(isalt > 0) then
 call iofil('salt_film_np.txt',micro_g(ng)%salt_film_np(:,i,j),m1,wf,fln,1)
 call iofil('salt_jet_np.txt' ,micro_g(ng)%salt_jet_np (:,i,j),m1,wf,fln,1)
 call iofil('salt_spum_np.txt',micro_g(ng)%salt_spum_np(:,i,j),m1,wf,fln,1)
 call iofil('salt_film_mp.txt',micro_g(ng)%salt_film_mp(:,i,j),m1,wf,fln,1)
 call iofil('salt_jet_mp.txt' ,micro_g(ng)%salt_jet_mp (:,i,j),m1,wf,fln,1)
 call iofil('salt_spum_mp.txt',micro_g(ng)%salt_spum_mp(:,i,j),m1,wf,fln,1)
endif

!Cloud water
if(level == 2 .or. level == 3) then
 if(icloud >= 1) call iofil('rcp.txt',micro_g(ng)%rcp(:,i,j),m1,wf,fln,1)
endif

!Full level=3 microphysics variables 
if(level == 3) then
 !Surface precip type variables
 call iofil('pcpg.txt' ,micro_g(ng)%pcpg (i,j),1,wf,fln,1)
 call iofil('qpcpg.txt',micro_g(ng)%qpcpg(i,j),1,wf,fln,1)
 call iofil('dpcpg.txt',micro_g(ng)%dpcpg(i,j),1,wf,fln,1)
 !Ice nuclei
 if(ipris>=5 .and. (iifn==1.or.iifn==2)) then
  call iofil('cifnp.txt',micro_g(ng)%cifnp(:,i,j),m1,wf,fln,1)
 endif
 !Drizzle droplet vars
 if(idriz >= 1) then
  call iofil('rdp.txt'  ,micro_g(ng)%rdp  (:,i,j),m1,wf,fln,1)
  call iofil('pcpvd.txt',micro_g(ng)%pcpvd(:,i,j),m1,wf,fln,1)
  call iofil('pcprd.txt',micro_g(ng)%pcprd(i,j),1,wf,fln,1)
  call iofil('accpd.txt',micro_g(ng)%accpd(i,j),1,wf,fln,1)
 endif
 !Rain hydrometeor vars
 if(irain >= 1) then 
  call iofil('rrp.txt'  ,micro_g(ng)%rrp  (:,i,j),m1,wf,fln,1)
  call iofil('pcpvr.txt',micro_g(ng)%pcpvr(:,i,j),m1,wf,fln,1)
  call iofil('q2.txt'   ,micro_g(ng)%q2   (:,i,j),m1,wf,fln,1)
  call iofil('pcprr.txt',micro_g(ng)%pcprr(i,j),1,wf,fln,1)
  call iofil('accpr.txt',micro_g(ng)%accpr(i,j),1,wf,fln,1)
 endif
 !Pristine ice hydrometeor  vars
 if(ipris >= 1) then
  call iofil('rpp.txt'  ,micro_g(ng)%rpp  (:,i,j),m1,wf,fln,1)
  call iofil('pcpvp.txt',micro_g(ng)%pcpvp(:,i,j),m1,wf,fln,1)
  call iofil('pcprp.txt',micro_g(ng)%pcprp(i,j),1,wf,fln,1)
  call iofil('accpp.txt',micro_g(ng)%accpp(i,j),1,wf,fln,1)
 endif
 !Snow hydrometeor vars
 if(isnow >= 1) then
  call iofil('rsp.txt'  ,micro_g(ng)%rsp  (:,i,j),m1,wf,fln,1)
  call iofil('pcpvs.txt',micro_g(ng)%pcpvs(:,i,j),m1,wf,fln,1)
  call iofil('pcprs.txt',micro_g(ng)%pcprs(i,j),1,wf,fln,1)
  call iofil('accps.txt',micro_g(ng)%accps(i,j),1,wf,fln,1)
 endif
 !Aggregate hydrometeor vars
 if(iaggr >= 1) then
  call iofil('rap.txt'  ,micro_g(ng)%rap  (:,i,j),m1,wf,fln,1)
  call iofil('pcpva.txt',micro_g(ng)%pcpva(:,i,j),m1,wf,fln,1)
  call iofil('pcpra.txt',micro_g(ng)%pcpra(i,j),1,wf,fln,1)
  call iofil('accpa.txt',micro_g(ng)%accpa(i,j),1,wf,fln,1)
 endif
 !Graupel hydrometeor vars
 if(igraup >= 1) then
  call iofil('rgp.txt'  ,micro_g(ng)%rgp  (:,i,j),m1,wf,fln,1)
  call iofil('pcpvg.txt',micro_g(ng)%pcpvg(:,i,j),m1,wf,fln,1)
  call iofil('q6.txt'   ,micro_g(ng)%q6   (:,i,j),m1,wf,fln,1)
  call iofil('pcprg.txt',micro_g(ng)%pcprg(i,j),1,wf,fln,1)
  call iofil('accpg.txt',micro_g(ng)%accpg(i,j),1,wf,fln,1)
 endif
 !Hail hydrometeor vars
 if(ihail >=1) then
  call iofil('rhp.txt'  ,micro_g(ng)%rhp  (:,i,j),m1,wf,fln,1)
  call iofil('pcpvh.txt',micro_g(ng)%pcpvh(:,i,j),m1,wf,fln,1)
  call iofil('q7.txt'   ,micro_g(ng)%q7   (:,i,j),m1,wf,fln,1)
  call iofil('pcprh.txt',micro_g(ng)%pcprh(i,j),1,wf,fln,1)
  call iofil('accph.txt',micro_g(ng)%accph(i,j),1,wf,fln,1)
 endif
 !Two-moment number concentration vars
 if(jnmb(1) >= 5) call iofil('ccp.txt',micro_g(ng)%ccp(:,i,j),m1,wf,fln,1)
 if(jnmb(8) >= 5) call iofil('cdp.txt',micro_g(ng)%cdp(:,i,j),m1,wf,fln,1)
 if(jnmb(2) >= 5) call iofil('crp.txt',micro_g(ng)%crp(:,i,j),m1,wf,fln,1)
 if(jnmb(3) >= 5) call iofil('cpp.txt',micro_g(ng)%cpp(:,i,j),m1,wf,fln,1)
 if(jnmb(4) >= 5) call iofil('csp.txt',micro_g(ng)%csp(:,i,j),m1,wf,fln,1)
 if(jnmb(5) >= 5) call iofil('cap.txt',micro_g(ng)%cap(:,i,j),m1,wf,fln,1)
 if(jnmb(6) >= 5) call iofil('cgp.txt',micro_g(ng)%cgp(:,i,j),m1,wf,fln,1)
 if(jnmb(7) >= 5) call iofil('chp.txt',micro_g(ng)%chp(:,i,j),m1,wf,fln,1)

 !Variables related to aerosol tracking in hydrometeors
 if(iccnlev >= 2) then
  !Regenerated aerosol vars
  call iofil('regen_aero1_np.txt',micro_g(ng)%regen_aero1_np(:,i,j),m1,wf,fln,1)
  call iofil('regen_aero1_mp.txt',micro_g(ng)%regen_aero1_mp(:,i,j),m1,wf,fln,1)
  call iofil('regen_aero2_np.txt',micro_g(ng)%regen_aero2_np(:,i,j),m1,wf,fln,1)
  call iofil('regen_aero2_mp.txt',micro_g(ng)%regen_aero2_mp(:,i,j),m1,wf,fln,1)
  !Tracking ALL aerosol mass in hydrometeor vars
  if(jnmb(1)>=1) call iofil('cnmcp.txt',micro_g(ng)%cnmcp(:,i,j),m1,wf,fln,1)
  if(jnmb(8)>=1) call iofil('cnmdp.txt',micro_g(ng)%cnmdp(:,i,j),m1,wf,fln,1)
  if(jnmb(2)>=1) call iofil('cnmrp.txt',micro_g(ng)%cnmrp(:,i,j),m1,wf,fln,1)
  if(jnmb(3)>=1) call iofil('cnmpp.txt',micro_g(ng)%cnmpp(:,i,j),m1,wf,fln,1)
  if(jnmb(4)>=1) call iofil('cnmsp.txt',micro_g(ng)%cnmsp(:,i,j),m1,wf,fln,1)
  if(jnmb(5)>=1) call iofil('cnmap.txt',micro_g(ng)%cnmap(:,i,j),m1,wf,fln,1)
  if(jnmb(6)>=1) call iofil('cnmgp.txt',micro_g(ng)%cnmgp(:,i,j),m1,wf,fln,1)
  if(jnmb(7)>=1) call iofil('cnmhp.txt',micro_g(ng)%cnmhp(:,i,j),m1,wf,fln,1)
  !Accumulated ALL aerosol mass at surface via sedimentation in hydrometeors
  call iofil('accpaero.txt',micro_g(ng)%accpaero(i,j),1,wf,fln,1)
  call iofil('pcpraero.txt',micro_g(ng)%pcpraero(i,j),1,wf,fln,1)
  !Tracking DUST in hydrometeors vars
  if(itrkdust==1 .and. idust>0) then
   !Accumulated DUST  mass at surface via sedimentation in hydrometeors
   call iofil('accpdust.txt',micro_g(ng)%accpdust(i,j),1,wf,fln,1)
   call iofil('pcprdust.txt',micro_g(ng)%pcprdust(i,j),1,wf,fln,1)
   !Tracking dust as CCN + IN in hydrometeors vars
   if(jnmb(1)>=1) call iofil('dnmcp.txt',micro_g(ng)%dnmcp(:,i,j),m1,wf,fln,1)
   if(jnmb(8)>=1) call iofil('dnmdp.txt',micro_g(ng)%dnmdp(:,i,j),m1,wf,fln,1)
   if(jnmb(2)>=1) call iofil('dnmrp.txt',micro_g(ng)%dnmrp(:,i,j),m1,wf,fln,1)
   if(jnmb(3)>=1) call iofil('dnmpp.txt',micro_g(ng)%dnmpp(:,i,j),m1,wf,fln,1)
   if(jnmb(4)>=1) call iofil('dnmsp.txt',micro_g(ng)%dnmsp(:,i,j),m1,wf,fln,1)
   if(jnmb(5)>=1) call iofil('dnmap.txt',micro_g(ng)%dnmap(:,i,j),m1,wf,fln,1)
   if(jnmb(6)>=1) call iofil('dnmgp.txt',micro_g(ng)%dnmgp(:,i,j),m1,wf,fln,1)
   if(jnmb(7)>=1) call iofil('dnmhp.txt',micro_g(ng)%dnmhp(:,i,j),m1,wf,fln,1)
  endif
  !Tracking dust as Ice Nuclei in hydrometeors vars
  if(itrkdustifn==1 .and. idust>0) then
   if(jnmb(1)>=1) call iofil('dincp.txt',micro_g(ng)%dincp(:,i,j),m1,wf,fln,1)
   if(jnmb(8)>=1) call iofil('dindp.txt',micro_g(ng)%dindp(:,i,j),m1,wf,fln,1)
   if(jnmb(2)>=1) call iofil('dinrp.txt',micro_g(ng)%dinrp(:,i,j),m1,wf,fln,1)
   if(jnmb(3)>=1) call iofil('dinpp.txt',micro_g(ng)%dinpp(:,i,j),m1,wf,fln,1)
   if(jnmb(4)>=1) call iofil('dinsp.txt',micro_g(ng)%dinsp(:,i,j),m1,wf,fln,1)
   if(jnmb(5)>=1) call iofil('dinap.txt',micro_g(ng)%dinap(:,i,j),m1,wf,fln,1)
   if(jnmb(6)>=1) call iofil('dingp.txt',micro_g(ng)%dingp(:,i,j),m1,wf,fln,1)
   if(jnmb(7)>=1) call iofil('dinhp.txt',micro_g(ng)%dinhp(:,i,j),m1,wf,fln,1)
  endif
  !Tracking solubility fraction for regenerated aerosols vars
  if(itrkepsilon==1) then
   call iofil('resol_aero1_mp.txt',micro_g(ng)%resol_aero1_mp(:,i,j),m1,wf,fln,1)
   call iofil('resol_aero2_mp.txt',micro_g(ng)%resol_aero2_mp(:,i,j),m1,wf,fln,1)
   if(jnmb(1)>=1) call iofil('snmcp.txt',micro_g(ng)%snmcp(:,i,j),m1,wf,fln,1)
   if(jnmb(8)>=1) call iofil('snmdp.txt',micro_g(ng)%snmdp(:,i,j),m1,wf,fln,1)
   if(jnmb(2)>=1) call iofil('snmrp.txt',micro_g(ng)%snmrp(:,i,j),m1,wf,fln,1)
   if(jnmb(3)>=1) call iofil('snmpp.txt',micro_g(ng)%snmpp(:,i,j),m1,wf,fln,1)
   if(jnmb(4)>=1) call iofil('snmsp.txt',micro_g(ng)%snmsp(:,i,j),m1,wf,fln,1)
   if(jnmb(5)>=1) call iofil('snmap.txt',micro_g(ng)%snmap(:,i,j),m1,wf,fln,1)
   if(jnmb(6)>=1) call iofil('snmgp.txt',micro_g(ng)%snmgp(:,i,j),m1,wf,fln,1)
   if(jnmb(7)>=1) call iofil('snmhp.txt',micro_g(ng)%snmhp(:,i,j),m1,wf,fln,1)
  endif
 endif
 !Vars for tracking Ice nucleation of aerosols (D > 0.5 microns) (DeMott method)
 if(iifn==3 .and. iccnlev>=1) then
  if(jnmb(1)>=5)call iofil('ifnnucp.txt',micro_g(ng)%ifnnucp(:,i,j),m1,wf,fln,1)
  if(jnmb(1)>=5)call iofil('immercp.txt',micro_g(ng)%immercp(:,i,j),m1,wf,fln,1)
  if(jnmb(8)>=5)call iofil('immerdp.txt',micro_g(ng)%immerdp(:,i,j),m1,wf,fln,1)
  if(jnmb(2)>=5)call iofil('immerrp.txt',micro_g(ng)%immerrp(:,i,j),m1,wf,fln,1)
 endif
 !Microphysics budget variables
 if(imbudget>=1) then !18
  call iofil('latheatvap.txt' ,micro_g(ng)%latheatvap (:,i,j),m1,wf,fln,1)
  call iofil('latheatfrz.txt' ,micro_g(ng)%latheatfrz (:,i,j),m1,wf,fln,1)
  call iofil('nuccldrt.txt'   ,micro_g(ng)%nuccldrt   (:,i,j),m1,wf,fln,1)
  call iofil('cld2raint.txt'  ,micro_g(ng)%cld2raint  (:,i,j),m1,wf,fln,1)
  call iofil('ice2raint.txt'  ,micro_g(ng)%ice2raint  (:,i,j),m1,wf,fln,1)
  call iofil('nucicert.txt'   ,micro_g(ng)%nucicert   (:,i,j),m1,wf,fln,1)
  call iofil('vapliqt.txt'    ,micro_g(ng)%vapliqt    (:,i,j),m1,wf,fln,1)
  call iofil('vapicet.txt'    ,micro_g(ng)%vapicet    (:,i,j),m1,wf,fln,1)
  call iofil('evapliqt.txt'   ,micro_g(ng)%evapliqt   (:,i,j),m1,wf,fln,1)
  call iofil('evapicet.txt'   ,micro_g(ng)%evapicet   (:,i,j),m1,wf,fln,1)
  call iofil('freezingt.txt'  ,micro_g(ng)%freezingt  (:,i,j),m1,wf,fln,1)
  call iofil('meltingt.txt'   ,micro_g(ng)%meltingt   (:,i,j),m1,wf,fln,1)
  call iofil('melticet.txt'   ,micro_g(ng)%melticet   (:,i,j),m1,wf,fln,1)
  call iofil('rimecldt.txt'   ,micro_g(ng)%rimecldt   (:,i,j),m1,wf,fln,1)
  call iofil('rain2icet.txt'  ,micro_g(ng)%rain2icet  (:,i,j),m1,wf,fln,1)
  call iofil('aggregatet.txt' ,micro_g(ng)%aggregatet (:,i,j),m1,wf,fln,1)
  call iofil('latheatvapt.txt',micro_g(ng)%latheatvapt(:,i,j),m1,wf,fln,1)
  call iofil('latheatfrzt.txt',micro_g(ng)%latheatfrzt(:,i,j),m1,wf,fln,1)
 endif
 if(imbudget>=2) then !37
  call iofil('inuchomrt.txt'    ,micro_g(ng)%inuchomrt    (:,i,j),m1,wf,fln,1)
  call iofil('inuccontrt.txt'   ,micro_g(ng)%inuccontrt   (:,i,j),m1,wf,fln,1)
  call iofil('inucifnrt.txt'    ,micro_g(ng)%inucifnrt    (:,i,j),m1,wf,fln,1)
  call iofil('inuchazrt.txt'    ,micro_g(ng)%inuchazrt    (:,i,j),m1,wf,fln,1)
  call iofil('vapcldt.txt'      ,micro_g(ng)%vapcldt      (:,i,j),m1,wf,fln,1)
  call iofil('vapraint.txt'     ,micro_g(ng)%vapraint     (:,i,j),m1,wf,fln,1)
  call iofil('vapprist.txt'     ,micro_g(ng)%vapprist     (:,i,j),m1,wf,fln,1)
  call iofil('vapsnowt.txt'     ,micro_g(ng)%vapsnowt     (:,i,j),m1,wf,fln,1)
  call iofil('vapaggrt.txt'     ,micro_g(ng)%vapaggrt     (:,i,j),m1,wf,fln,1)
  call iofil('vapgraut.txt'     ,micro_g(ng)%vapgraut     (:,i,j),m1,wf,fln,1)
  call iofil('vaphailt.txt'     ,micro_g(ng)%vaphailt     (:,i,j),m1,wf,fln,1)
  call iofil('vapdrizt.txt'     ,micro_g(ng)%vapdrizt     (:,i,j),m1,wf,fln,1)
  call iofil('evapcldt.txt'     ,micro_g(ng)%evapcldt     (:,i,j),m1,wf,fln,1)
  call iofil('evapraint.txt'    ,micro_g(ng)%evapraint    (:,i,j),m1,wf,fln,1)
  call iofil('evapprist.txt'    ,micro_g(ng)%evapprist    (:,i,j),m1,wf,fln,1)
  call iofil('evapsnowt.txt'    ,micro_g(ng)%evapsnowt    (:,i,j),m1,wf,fln,1)
  call iofil('evapaggrt.txt'    ,micro_g(ng)%evapaggrt    (:,i,j),m1,wf,fln,1)
  call iofil('evapgraut.txt'    ,micro_g(ng)%evapgraut    (:,i,j),m1,wf,fln,1)
  call iofil('evaphailt.txt'    ,micro_g(ng)%evaphailt    (:,i,j),m1,wf,fln,1)
  call iofil('evapdrizt.txt'    ,micro_g(ng)%evapdrizt    (:,i,j),m1,wf,fln,1)
  call iofil('meltprist.txt'    ,micro_g(ng)%meltprist    (:,i,j),m1,wf,fln,1)
  call iofil('meltsnowt.txt'    ,micro_g(ng)%meltsnowt    (:,i,j),m1,wf,fln,1)
  call iofil('meltaggrt.txt'    ,micro_g(ng)%meltaggrt    (:,i,j),m1,wf,fln,1)
  call iofil('meltgraut.txt'    ,micro_g(ng)%meltgraut    (:,i,j),m1,wf,fln,1)
  call iofil('melthailt.txt'    ,micro_g(ng)%melthailt    (:,i,j),m1,wf,fln,1)
  call iofil('rimecldsnowt.txt' ,micro_g(ng)%rimecldsnowt (:,i,j),m1,wf,fln,1)
  call iofil('rimecldaggrt.txt' ,micro_g(ng)%rimecldaggrt (:,i,j),m1,wf,fln,1)
  call iofil('rimecldgraut.txt' ,micro_g(ng)%rimecldgraut (:,i,j),m1,wf,fln,1)
  call iofil('rimecldhailt.txt' ,micro_g(ng)%rimecldhailt (:,i,j),m1,wf,fln,1)
  call iofil('rain2prt.txt'     ,micro_g(ng)%rain2prt     (:,i,j),m1,wf,fln,1)
  call iofil('rain2snt.txt'     ,micro_g(ng)%rain2snt     (:,i,j),m1,wf,fln,1)
  call iofil('rain2agt.txt'     ,micro_g(ng)%rain2agt     (:,i,j),m1,wf,fln,1)
  call iofil('rain2grt.txt'     ,micro_g(ng)%rain2grt     (:,i,j),m1,wf,fln,1)
  call iofil('rain2hat.txt'     ,micro_g(ng)%rain2hat     (:,i,j),m1,wf,fln,1)
  call iofil('aggrselfprist.txt',micro_g(ng)%aggrselfprist(:,i,j),m1,wf,fln,1)
  call iofil('aggrselfsnowt.txt',micro_g(ng)%aggrselfsnowt(:,i,j),m1,wf,fln,1)
  call iofil('aggrprissnowt.txt',micro_g(ng)%aggrprissnowt(:,i,j),m1,wf,fln,1)
 endif
 if(imbudget==3 .and. idust>0) then !4
  call iofil('dust1cldrt.txt',micro_g(ng)%dust1cldrt(:,i,j),m1,wf,fln,1)
  call iofil('dust2cldrt.txt',micro_g(ng)%dust2cldrt(:,i,j),m1,wf,fln,1)
  call iofil('dust1drzrt.txt',micro_g(ng)%dust1drzrt(:,i,j),m1,wf,fln,1)
  call iofil('dust2drzrt.txt',micro_g(ng)%dust2drzrt(:,i,j),m1,wf,fln,1)
 endif
endif !if level == 3

return
END SUBROUTINE readwrite_scm

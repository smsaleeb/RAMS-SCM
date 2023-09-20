!##############################################################################
Subroutine trsets_ns ()

use mem_grid
use mem_basic
use mem_radiate
use mem_micro
use micphys
use node_mod, only:mzp,mxp,myp

implicit none

integer :: m1,m2,m3,i,j

!Here we set top and bottom boundaries for non-advected
!3D vars or "non-scalars" since our definition of scalar refers
!to 3D vars that are advected and diffused.

!Note that in SCM version we add all the advected scalars here.
!In standard RAMS, this is done for advected scalars in another
!part of the code that we do not need for SCM.

CALL trsetns (mzp,mxp,myp,basic_g(ngrid)%thp,2)
CALL trsetns (mzp,mxp,myp,basic_g(ngrid)%rtp,2)

if(ilwrtyp+iswrtyp > 0) then
   CALL trsetns (mzp,mxp,myp,radiate_g(ngrid)%fthrd,1)
endif

if(ilwrtyp == 3 .or. iswrtyp == 3) then
   CALL trsetns (mzp,mxp,myp,radiate_g(ngrid)%bext,2)
endif

if (ipris >=5 .and. (iifn == 1 .or. iifn == 2)) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cifnp,2)
endif
if (iaerosol > 0) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cn1np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cn1mp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cn2np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cn2mp,2)
endif
if (idust > 0) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%md1np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%md1mp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%md2np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%md2mp,2)
endif
if (isalt > 0) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%salt_film_np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%salt_film_mp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%salt_jet_np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%salt_jet_mp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%salt_spum_np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%salt_spum_mp,2)
endif
if (iabcarb > 0) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%abc1np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%abc1mp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%abc2np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%abc2mp,2)
endif

if (iccnlev>=2) then
   if (jnmb(1) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmcp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmcp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmcp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dincp,2)
   endif
   if (jnmb(2) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmrp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmrp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmrp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dinrp,2)
   endif
   if (jnmb(3) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmpp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmpp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmpp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dinpp,2)
   endif
   if (jnmb(4) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmsp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmsp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmsp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dinsp,2)
   endif
   if (jnmb(5) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmap,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmap,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmap,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dinap,2)
   endif
   if (jnmb(6) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmgp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmgp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmgp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dingp,2)
   endif
   if (jnmb(7) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmhp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmhp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmhp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dinhp,2)
   endif
   if (jnmb(8) >= 1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cnmdp,2)
    if(itrkepsilon==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%snmdp,2)
    if(itrkdust==1)    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dnmdp,2)
    if(itrkdustifn==1) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%dindp,2)
   endif
   !Regenerated aerosol variables
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%regen_aero1_np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%regen_aero1_mp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%regen_aero2_np,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%regen_aero2_mp,2)
   if(itrkepsilon==1) then
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%resol_aero1_mp,2)
    CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%resol_aero2_mp,2)
   endif
endif

if (iifn==3 .and. iccnlev>=1) then
   if (jnmb(1) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%ifnnucp,2)
   if (jnmb(1) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%immercp,2)
   if (jnmb(2) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%immerrp,2)
   if (jnmb(8) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%immerdp,2)
endif

if (jnmb(1) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rcp,2)
   if (jnmb(1) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%ccp,2)
endif

if (jnmb(2) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rrp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%q2,1)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvr,2)
   if (jnmb(2) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%crp,2)
endif

if (jnmb(3) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rpp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvp,2)
   if (jnmb(3) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cpp,2)
endif

if (jnmb(4) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rsp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvs,2)
   if (jnmb(4) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%csp,2)
endif

if (jnmb(5) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rap,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpva,2)
   if (jnmb(5) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cap,2)
endif

if (jnmb(6) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rgp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%q6,1)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvg,2)
   if (jnmb(6) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cgp,2)
endif

if (jnmb(7) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rhp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%q7,1)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvh,2)
   if (jnmb(7) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%chp,2) 
endif

if (jnmb(8) >= 1) then
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%rdp,2)
   CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvd,2)
   if (jnmb(8) >= 5) CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%cdp,2)  
endif

return
END SUBROUTINE trsets_ns

!##############################################################################
Subroutine trsetns (m1,m2,m3,ap,flg)

use mem_grid

implicit none

real, dimension(m1,m2,m3) :: ap
integer :: m1,m2,m3,i,j,flg
real :: dzmr

dzmr = dzm(m1-2) / dzm(m1-1)

!Top and bottom boundary setting
do j = 1,m3
 do i = 1,m2

  if(flg==1) &
   ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))

  if(flg==2) &
   ap(m1,i,j) = max(0., ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j)))

  ap(1,i,j)  = ap(2,i,j)

 enddo
enddo

return
END SUBROUTINE trsetns


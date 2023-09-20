!##############################################################################
Module mem_micro

implicit none

   Type micro_vars
   
   ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          rcp,rdp,rrp,rpp,rsp,rap,rgp,rhp           &
                         ,ccp,cdp,crp,cpp,csp,cap,cgp,chp           &
                         ,cifnp,q2,q6,q7                            &
                         ,pcpvd,pcpvr,pcpvp,pcpvs,pcpva,pcpvg,pcpvh &
                         !Bin precip vars
                         ,pcpvic,pcpvip,pcpvid                      &
                         !Aerosol categories mass and number
                         ,cn1np,cn2np,cn1mp,cn2mp                   &
                         ,md1np,md2np,md1mp,md2mp                   &
                         ,salt_film_np,salt_jet_np,salt_spum_np     &
                         ,salt_film_mp,salt_jet_mp,salt_spum_mp     &
                         ,abc1np,abc2np,abc1mp,abc2mp               &
                         ,regen_aero1_np,regen_aero1_mp             &
                         ,regen_aero2_np,regen_aero2_mp             &
                         !Immersion freezing nuclei tracking
                         ,immercp,immerdp,immerrp,ifnnucp           &
                         !Aerosol tracking variables
                         ,cnmcp,cnmdp,cnmrp,cnmpp,cnmsp             &
                         ,cnmap,cnmgp,cnmhp                         &
                         ,dnmcp,dnmdp,dnmrp,dnmpp,dnmsp             &
                         ,dnmap,dnmgp,dnmhp                         &
                         ,dincp,dindp,dinrp,dinpp,dinsp             &
                         ,dinap,dingp,dinhp                         &
                         ,snmcp,snmdp,snmrp,snmpp,snmsp             &
                         ,snmap,snmgp,snmhp                         &
                         ,resol_aero1_mp,resol_aero2_mp &
           ! MICRO BUDGET PROCESSES (imbudget >=1)
           ,latheatvap,latheatfrz,nuccldrt,cld2raint,ice2raint,nucicert    &
           ,vapliqt,vapicet,evapliqt,evapicet,freezingt,meltingt           &
           ,melticet,rimecldt,rain2icet,aggregatet                         &
           ,latheatvapt,latheatfrzt                                        &
           ! MICRO BUDGET PROCESSES (imbudget >=2)
           ,inuchomrt,inuccontrt,inucifnrt,inuchazrt,vapcldt,vapraint      &
           ,vapprist,vapsnowt,vapaggrt,vapgraut,vaphailt,vapdrizt          &
           ,evapcldt,evapraint,evapprist,evapsnowt,evapaggrt,evapgraut     &
           ,evaphailt,evapdrizt                                            &
           ,meltprist,meltsnowt,meltaggrt,meltgraut,melthailt              &
           ,rimecldsnowt,rimecldaggrt,rimecldgraut,rimecldhailt            &
           ,rain2prt,rain2snt,rain2agt,rain2grt,rain2hat                   &
           ,aggrselfprist,aggrselfsnowt,aggrprissnowt                      &
           ! MICRO BUDGET PROCESSES (imbudget >=3)
           ,dust1cldrt,dust2cldrt,dust1drzrt,dust2drzrt

   ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          accpr,accpp,accps,accpa,accpg,accph,accpd &
                         ,pcprr,pcprp,pcprs,pcpra,pcprg,pcprh,pcprd &
                         ,pcpg,qpcpg,dpcpg                          &
                         !Accumulated aerosols and accumulation rate
                         ,accpdust,pcprdust,accpaero,pcpraero       &
                         !Dust erodible fraction
                         ,dustfrac
                          
   End Type               

   type (micro_vars), allocatable :: micro_g(:)

   !Sedimentation table variables
   Type pcp_tab_type
   real, allocatable, dimension(:,:,:,:,:,:) :: pcpfillc,pcpfillr
   real, allocatable, dimension(:,:,:,:,:)   :: sfcpcp
   real, allocatable, dimension(:,:,:,:,:)   :: allpcp
   End Type
   type (pcp_tab_type), allocatable :: pcp_tab(:)
                                                 
Contains                  

!##############################################################################
Subroutine alloc_sedim (pcp_tab,n1)

use micphys

implicit none

   type (pcp_tab_type) :: pcp_tab
   integer, intent(in) :: n1

! Micro Level=3 sedimentation tables
   if (level == 3) then
    allocate (pcp_tab%pcpfillc(n1,maxkfall,nembfall,nhcat,ndensrtgt,nband))
    allocate (pcp_tab%pcpfillr(n1,maxkfall,nembfall,nhcat,ndensrtgt,nband))
    allocate (pcp_tab%sfcpcp(maxkfall,nembfall,nhcat,ndensrtgt,nband))
    allocate (pcp_tab%allpcp(n1,nembfall,nhcat,ndensrtgt,nband))
   endif

return
END SUBROUTINE alloc_sedim

!##############################################################################
Subroutine dealloc_sedim (pcp_tab)

implicit none

   type (pcp_tab_type) :: pcp_tab

   if (allocated(pcp_tab%pcpfillc)) deallocate(pcp_tab%pcpfillc)
   if (allocated(pcp_tab%pcpfillr)) deallocate(pcp_tab%pcpfillr)
   if (allocated(pcp_tab%sfcpcp))   deallocate(pcp_tab%sfcpcp)
   if (allocated(pcp_tab%allpcp))   deallocate(pcp_tab%allpcp)

return
END SUBROUTINE dealloc_sedim

!##############################################################################
Subroutine alloc_micro (micro,n1,n2,n3,n4)

use micphys

implicit none          

   type (micro_vars) :: micro
   integer, intent(in) :: n1,n2,n3,n4
   integer :: numpts3d,numpts2d

   numpts3d=n1*n2*n3
   numpts2d=n2*n3

! Allocate arrays based on options (if necessary)
      if (level >= 0 .and. level .ne. 4) then
         if(iaerosol > 0) then
            allocate (micro%cn1np(n1,n2,n3))
            allocate (micro%cn1mp(n1,n2,n3))
            allocate (micro%cn2np(n1,n2,n3))
            allocate (micro%cn2mp(n1,n2,n3))
            call azero (numpts3d,micro%cn1np)
            call azero (numpts3d,micro%cn1mp)
            call azero (numpts3d,micro%cn2np)
            call azero (numpts3d,micro%cn2mp)
         endif
         if(idust > 0) then
            allocate (micro%md1np(n1,n2,n3))
            allocate (micro%md2np(n1,n2,n3))
            allocate (micro%md1mp(n1,n2,n3))
            allocate (micro%md2mp(n1,n2,n3))
            call azero (numpts3d,micro%md1np)
            call azero (numpts3d,micro%md2np)
            call azero (numpts3d,micro%md1mp)
            call azero (numpts3d,micro%md2mp)
            if(idust == 2) then
              allocate (micro%dustfrac(n2,n3))
              call azero (numpts2d,micro%dustfrac)
            endif
         endif
         if(isalt > 0) then
            allocate (micro%salt_film_np(n1,n2,n3))
            allocate (micro%salt_jet_np(n1,n2,n3))
            allocate (micro%salt_spum_np(n1,n2,n3))
            allocate (micro%salt_film_mp(n1,n2,n3))
            allocate (micro%salt_jet_mp(n1,n2,n3))
            allocate (micro%salt_spum_mp(n1,n2,n3))
            call azero (numpts3d,micro%salt_film_np)
            call azero (numpts3d,micro%salt_jet_np)
            call azero (numpts3d,micro%salt_spum_np)
            call azero (numpts3d,micro%salt_film_mp)
            call azero (numpts3d,micro%salt_jet_mp)
            call azero (numpts3d,micro%salt_spum_mp)
         endif
         if(iabcarb > 0) then
            allocate (micro%abc1np(n1,n2,n3))
            allocate (micro%abc2np(n1,n2,n3))
            allocate (micro%abc1mp(n1,n2,n3))
            allocate (micro%abc2mp(n1,n2,n3))
            call azero (numpts3d,micro%abc1np)
            call azero (numpts3d,micro%abc2np)
            call azero (numpts3d,micro%abc1mp)
            call azero (numpts3d,micro%abc2mp)
         endif
      endif

      if (level >= 2 .and. level .ne. 4) then
         if(icloud >= 1) then
           allocate (micro%rcp(n1,n2,n3))
           call azero (numpts3d,micro%rcp)
         endif
      endif

      if (level == 3) then
         if(ipris>=5 .and. (iifn==1.or.iifn==2)) allocate (micro%cifnp(n1,n2,n3))
         if(idriz >= 1)  then
            allocate (micro%rdp(n1,n2,n3))
            allocate (micro%accpd(n2,n3))
            allocate (micro%pcprd(n2,n3))
            allocate (micro%pcpvd(n1,n2,n3))
            call azero (numpts3d,micro%rdp)
            call azero (numpts2d,micro%accpd)
            call azero (numpts2d,micro%pcprd)
            call azero (numpts3d,micro%pcpvd)
         endif
         if(irain >= 1)  then
            allocate (micro%rrp(n1,n2,n3))
            allocate (micro%accpr(n2,n3))
            allocate (micro%pcprr(n2,n3))
            allocate (micro%pcpvr(n1,n2,n3))
            allocate (micro%q2(n1,n2,n3))
            call azero (numpts3d,micro%rrp)
            call azero (numpts2d,micro%accpr)
            call azero (numpts2d,micro%pcprr)
            call azero (numpts3d,micro%pcpvr)
            call azero (numpts3d,micro%q2)
         endif
         if(ipris >= 1)  then
            allocate (micro%rpp(n1,n2,n3))
            allocate (micro%accpp(n2,n3))
            allocate (micro%pcprp(n2,n3))
            allocate (micro%pcpvp(n1,n2,n3))
            call azero (numpts3d,micro%rpp)
            call azero (numpts2d,micro%accpp)
            call azero (numpts2d,micro%pcprp)
            call azero (numpts3d,micro%pcpvp)
         endif
         if(isnow >= 1)  then
            allocate (micro%rsp(n1,n2,n3))
            allocate (micro%accps(n2,n3))
            allocate (micro%pcprs(n2,n3))
            allocate (micro%pcpvs(n1,n2,n3))
            call azero (numpts3d,micro%rsp)
            call azero (numpts2d,micro%accps)
            call azero (numpts2d,micro%pcprs)
            call azero (numpts3d,micro%pcpvs)
         endif
         if(iaggr >= 1)  then
            allocate (micro%rap(n1,n2,n3))
            allocate (micro%accpa(n2,n3))
            allocate (micro%pcpra(n2,n3))
            allocate (micro%pcpva(n1,n2,n3))
            call azero (numpts3d,micro%rap)
            call azero (numpts2d,micro%accpa)
            call azero (numpts2d,micro%pcpra)
            call azero (numpts3d,micro%pcpva)
         endif
         if(igraup >= 1) then
            allocate (micro%rgp(n1,n2,n3))
            allocate (micro%accpg(n2,n3))
            allocate (micro%pcprg(n2,n3))
            allocate (micro%pcpvg(n1,n2,n3))
            allocate (micro%q6(n1,n2,n3))
            call azero (numpts3d,micro%rgp)
            call azero (numpts2d,micro%accpg)
            call azero (numpts2d,micro%pcprg)
            call azero (numpts3d,micro%pcpvg)
            call azero (numpts3d,micro%q6)
         endif
         if(ihail >= 1)  then
            allocate (micro%rhp(n1,n2,n3))
            allocate (micro%accph(n2,n3))
            allocate (micro%pcprh(n2,n3))
            allocate (micro%pcpvh(n1,n2,n3))
            allocate (micro%q7(n1,n2,n3))
            call azero (numpts3d,micro%rhp)
            call azero (numpts2d,micro%accph)
            call azero (numpts2d,micro%pcprh)
            call azero (numpts3d,micro%pcpvh)
            call azero (numpts3d,micro%q7)
         endif

         if(jnmb(1) >= 5) allocate (micro%ccp(n1,n2,n3))
         if(jnmb(8) >= 5) allocate (micro%cdp(n1,n2,n3))
         if(jnmb(2) >= 5) allocate (micro%crp(n1,n2,n3))
         if(jnmb(3) >= 5) allocate (micro%cpp(n1,n2,n3))
         if(jnmb(4) >= 5) allocate (micro%csp(n1,n2,n3))
         if(jnmb(5) >= 5) allocate (micro%cap(n1,n2,n3))
         if(jnmb(6) >= 5) allocate (micro%cgp(n1,n2,n3))
         if(jnmb(7) >= 5) allocate (micro%chp(n1,n2,n3))

         if(jnmb(1) >= 5) call azero (numpts3d,micro%ccp)
         if(jnmb(8) >= 5) call azero (numpts3d,micro%cdp)
         if(jnmb(2) >= 5) call azero (numpts3d,micro%crp)
         if(jnmb(3) >= 5) call azero (numpts3d,micro%cpp)
         if(jnmb(4) >= 5) call azero (numpts3d,micro%csp)
         if(jnmb(5) >= 5) call azero (numpts3d,micro%cap)
         if(jnmb(6) >= 5) call azero (numpts3d,micro%cgp)
         if(jnmb(7) >= 5) call azero (numpts3d,micro%chp)

         allocate (micro%pcpg(n2,n3))
         allocate (micro%qpcpg(n2,n3))
         allocate (micro%dpcpg(n2,n3))
         call azero (numpts2d,micro%pcpg)
         call azero (numpts2d,micro%qpcpg)
         call azero (numpts2d,micro%dpcpg)

         if(iccnlev>=2) then
           allocate (micro%regen_aero1_np(n1,n2,n3))
           allocate (micro%regen_aero1_mp(n1,n2,n3))
           allocate (micro%regen_aero2_np(n1,n2,n3))
           allocate (micro%regen_aero2_mp(n1,n2,n3))
           call azero (numpts3d,micro%regen_aero1_np)
           call azero (numpts3d,micro%regen_aero1_mp)
           call azero (numpts3d,micro%regen_aero2_np)
           call azero (numpts3d,micro%regen_aero2_mp)

           if(jnmb(1) >= 1) allocate (micro%cnmcp(n1,n2,n3))
           if(jnmb(8) >= 1) allocate (micro%cnmdp(n1,n2,n3))
           if(jnmb(2) >= 1) allocate (micro%cnmrp(n1,n2,n3))
           if(jnmb(3) >= 1) allocate (micro%cnmpp(n1,n2,n3))
           if(jnmb(4) >= 1) allocate (micro%cnmsp(n1,n2,n3))
           if(jnmb(5) >= 1) allocate (micro%cnmap(n1,n2,n3))
           if(jnmb(6) >= 1) allocate (micro%cnmgp(n1,n2,n3))
           if(jnmb(7) >= 1) allocate (micro%cnmhp(n1,n2,n3))
           if(jnmb(1) >= 1) call azero (numpts3d,micro%cnmcp)
           if(jnmb(8) >= 1) call azero (numpts3d,micro%cnmdp)
           if(jnmb(2) >= 1) call azero (numpts3d,micro%cnmrp)
           if(jnmb(3) >= 1) call azero (numpts3d,micro%cnmpp)
           if(jnmb(4) >= 1) call azero (numpts3d,micro%cnmsp)
           if(jnmb(5) >= 1) call azero (numpts3d,micro%cnmap)
           if(jnmb(6) >= 1) call azero (numpts3d,micro%cnmgp)
           if(jnmb(7) >= 1) call azero (numpts3d,micro%cnmhp)

           allocate (micro%accpaero(n2,n3))
           allocate (micro%pcpraero(n2,n3))
           call azero (numpts2d,micro%accpaero)
           call azero (numpts2d,micro%pcpraero)

           if(itrkdust==1 .and. idust>0) then
            allocate (micro%accpdust(n2,n3))
            allocate (micro%pcprdust(n2,n3))
            call azero (numpts2d,micro%accpdust)
            call azero (numpts2d,micro%pcprdust)
            if(jnmb(1) >= 1) allocate (micro%dnmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%dnmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%dnmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%dnmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%dnmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%dnmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%dnmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%dnmhp(n1,n2,n3))
            if(jnmb(1) >= 1) call azero (numpts3d,micro%dnmcp)
            if(jnmb(8) >= 1) call azero (numpts3d,micro%dnmdp)
            if(jnmb(2) >= 1) call azero (numpts3d,micro%dnmrp)
            if(jnmb(3) >= 1) call azero (numpts3d,micro%dnmpp)
            if(jnmb(4) >= 1) call azero (numpts3d,micro%dnmsp)
            if(jnmb(5) >= 1) call azero (numpts3d,micro%dnmap)
            if(jnmb(6) >= 1) call azero (numpts3d,micro%dnmgp)
            if(jnmb(7) >= 1) call azero (numpts3d,micro%dnmhp)
           endif
           if(itrkdustifn==1 .and. idust>0) then
            if(jnmb(1) >= 1) allocate (micro%dincp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%dindp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%dinrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%dinpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%dinsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%dinap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%dingp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%dinhp(n1,n2,n3))
            if(jnmb(1) >= 1) call azero (numpts3d,micro%dincp)
            if(jnmb(8) >= 1) call azero (numpts3d,micro%dindp)
            if(jnmb(2) >= 1) call azero (numpts3d,micro%dinrp)
            if(jnmb(3) >= 1) call azero (numpts3d,micro%dinpp)
            if(jnmb(4) >= 1) call azero (numpts3d,micro%dinsp)
            if(jnmb(5) >= 1) call azero (numpts3d,micro%dinap)
            if(jnmb(6) >= 1) call azero (numpts3d,micro%dingp)
            if(jnmb(7) >= 1) call azero (numpts3d,micro%dinhp)
           endif
           if(itrkepsilon==1) then
            allocate (micro%resol_aero1_mp(n1,n2,n3))
            allocate (micro%resol_aero2_mp(n1,n2,n3))
            call azero (numpts3d,micro%resol_aero1_mp)
            call azero (numpts3d,micro%resol_aero2_mp)
            if(jnmb(1) >= 1) allocate (micro%snmcp(n1,n2,n3))
            if(jnmb(8) >= 1) allocate (micro%snmdp(n1,n2,n3))
            if(jnmb(2) >= 1) allocate (micro%snmrp(n1,n2,n3))
            if(jnmb(3) >= 1) allocate (micro%snmpp(n1,n2,n3))
            if(jnmb(4) >= 1) allocate (micro%snmsp(n1,n2,n3))
            if(jnmb(5) >= 1) allocate (micro%snmap(n1,n2,n3))
            if(jnmb(6) >= 1) allocate (micro%snmgp(n1,n2,n3))
            if(jnmb(7) >= 1) allocate (micro%snmhp(n1,n2,n3))
            if(jnmb(1) >= 1) call azero (numpts3d,micro%snmcp)
            if(jnmb(8) >= 1) call azero (numpts3d,micro%snmdp)
            if(jnmb(2) >= 1) call azero (numpts3d,micro%snmrp)
            if(jnmb(3) >= 1) call azero (numpts3d,micro%snmpp)
            if(jnmb(4) >= 1) call azero (numpts3d,micro%snmsp)
            if(jnmb(5) >= 1) call azero (numpts3d,micro%snmap)
            if(jnmb(6) >= 1) call azero (numpts3d,micro%snmgp)
            if(jnmb(7) >= 1) call azero (numpts3d,micro%snmhp)
           endif
         endif

         if(iifn==3 .and. iccnlev>=1) then
           if(jnmb(1) >= 5) allocate (micro%ifnnucp(n1,n2,n3))
           if(jnmb(1) >= 5) allocate (micro%immercp(n1,n2,n3))
           if(jnmb(8) >= 5) allocate (micro%immerdp(n1,n2,n3))
           if(jnmb(2) >= 5) allocate (micro%immerrp(n1,n2,n3))
           if(jnmb(1) >= 5) call azero (numpts3d,micro%ifnnucp)
           if(jnmb(1) >= 5) call azero (numpts3d,micro%immercp)
           if(jnmb(8) >= 5) call azero (numpts3d,micro%immerdp)
           if(jnmb(2) >= 5) call azero (numpts3d,micro%immerrp)
         endif

         !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
         !Do not need to zero these since this is done in micro anyway.
         if(imbudget>=1) then
           allocate (micro%latheatvap(n1,n2,n3))
           allocate (micro%latheatfrz(n1,n2,n3))
           allocate (micro%nuccldrt(n1,n2,n3))
           allocate (micro%cld2raint(n1,n2,n3))
           allocate (micro%ice2raint(n1,n2,n3))
           allocate (micro%nucicert(n1,n2,n3))
           allocate (micro%vapliqt(n1,n2,n3))
           allocate (micro%vapicet(n1,n2,n3))
           allocate (micro%evapliqt(n1,n2,n3))
           allocate (micro%evapicet(n1,n2,n3))
           allocate (micro%freezingt(n1,n2,n3))
           allocate (micro%meltingt(n1,n2,n3))
           allocate (micro%melticet(n1,n2,n3))
           allocate (micro%rimecldt(n1,n2,n3))
           allocate (micro%rain2icet(n1,n2,n3))
           allocate (micro%aggregatet(n1,n2,n3))
           allocate (micro%latheatvapt(n1,n2,n3))
           allocate (micro%latheatfrzt(n1,n2,n3))
         endif
         if(imbudget>=2) then
           allocate (micro%inuchomrt(n1,n2,n3))
           allocate (micro%inuccontrt(n1,n2,n3))
           allocate (micro%inucifnrt(n1,n2,n3))
           allocate (micro%inuchazrt(n1,n2,n3))
           allocate (micro%vapcldt(n1,n2,n3))
           allocate (micro%vapraint(n1,n2,n3))
           allocate (micro%vapprist(n1,n2,n3))
           allocate (micro%vapsnowt(n1,n2,n3))
           allocate (micro%vapaggrt(n1,n2,n3))
           allocate (micro%vapgraut(n1,n2,n3))
           allocate (micro%vaphailt(n1,n2,n3))
           allocate (micro%vapdrizt(n1,n2,n3))
           allocate (micro%evapcldt(n1,n2,n3))
           allocate (micro%evapraint(n1,n2,n3))
           allocate (micro%evapprist(n1,n2,n3))
           allocate (micro%evapsnowt(n1,n2,n3))
           allocate (micro%evapaggrt(n1,n2,n3))
           allocate (micro%evapgraut(n1,n2,n3))
           allocate (micro%evaphailt(n1,n2,n3))
           allocate (micro%evapdrizt(n1,n2,n3))
           allocate (micro%meltprist(n1,n2,n3))
           allocate (micro%meltsnowt(n1,n2,n3))
           allocate (micro%meltaggrt(n1,n2,n3))
           allocate (micro%meltgraut(n1,n2,n3))
           allocate (micro%melthailt(n1,n2,n3))
           allocate (micro%rimecldsnowt(n1,n2,n3))
           allocate (micro%rimecldaggrt(n1,n2,n3))
           allocate (micro%rimecldgraut(n1,n2,n3))
           allocate (micro%rimecldhailt(n1,n2,n3))
           allocate (micro%rain2prt(n1,n2,n3))
           allocate (micro%rain2snt(n1,n2,n3))
           allocate (micro%rain2agt(n1,n2,n3))
           allocate (micro%rain2grt(n1,n2,n3))
           allocate (micro%rain2hat(n1,n2,n3))
           allocate (micro%aggrselfprist(n1,n2,n3))
           allocate (micro%aggrselfsnowt(n1,n2,n3))
           allocate (micro%aggrprissnowt(n1,n2,n3))
         endif
         if(imbudget==3 .and. idust>0) then
           allocate (micro%dust1cldrt(n1,n2,n3))
           allocate (micro%dust2cldrt(n1,n2,n3))
           allocate (micro%dust1drzrt(n1,n2,n3))
           allocate (micro%dust2drzrt(n1,n2,n3))
         endif
      endif

return
END SUBROUTINE alloc_micro

!##############################################################################
Subroutine dealloc_micro (micro)

implicit none

   type (micro_vars) :: micro
   
   if (allocated(micro%rcp))     deallocate (micro%rcp)
   if (allocated(micro%rdp))     deallocate (micro%rdp)
   if (allocated(micro%rrp))     deallocate (micro%rrp)
   if (allocated(micro%rpp))     deallocate (micro%rpp)
   if (allocated(micro%rsp))     deallocate (micro%rsp)
   if (allocated(micro%rap))     deallocate (micro%rap)
   if (allocated(micro%rgp))     deallocate (micro%rgp)
   if (allocated(micro%rhp))     deallocate (micro%rhp)
   if (allocated(micro%ccp))     deallocate (micro%ccp)
   if (allocated(micro%cdp))     deallocate (micro%cdp)
   if (allocated(micro%crp))     deallocate (micro%crp)
   if (allocated(micro%cpp))     deallocate (micro%cpp)
   if (allocated(micro%csp))     deallocate (micro%csp)
   if (allocated(micro%cap))     deallocate (micro%cap)
   if (allocated(micro%cgp))     deallocate (micro%cgp)
   if (allocated(micro%chp))     deallocate (micro%chp)
   if (allocated(micro%cifnp))   deallocate (micro%cifnp)
   if (allocated(micro%q2))      deallocate (micro%q2)
   if (allocated(micro%q6))      deallocate (micro%q6)
   if (allocated(micro%q7))      deallocate (micro%q7)

   if (allocated(micro%cn1np))   deallocate (micro%cn1np)
   if (allocated(micro%cn2np))   deallocate (micro%cn2np)
   if (allocated(micro%cn1mp))   deallocate (micro%cn1mp)
   if (allocated(micro%cn2mp))   deallocate (micro%cn2mp)
   if (allocated(micro%md1np))   deallocate (micro%md1np)
   if (allocated(micro%md2np))   deallocate (micro%md2np)
   if (allocated(micro%md1mp))   deallocate (micro%md1mp)
   if (allocated(micro%md2mp))   deallocate (micro%md2mp)
   if (allocated(micro%dustfrac))deallocate (micro%dustfrac)
   if (allocated(micro%salt_film_np))  deallocate (micro%salt_film_np)
   if (allocated(micro%salt_jet_np))   deallocate (micro%salt_jet_np)
   if (allocated(micro%salt_spum_np))  deallocate (micro%salt_spum_np)
   if (allocated(micro%salt_film_mp))  deallocate (micro%salt_film_mp)
   if (allocated(micro%salt_jet_mp))   deallocate (micro%salt_jet_mp)
   if (allocated(micro%salt_spum_mp))  deallocate (micro%salt_spum_mp)
   if (allocated(micro%abc1np))   deallocate (micro%abc1np)
   if (allocated(micro%abc2np))   deallocate (micro%abc2np)
   if (allocated(micro%abc1mp))   deallocate (micro%abc1mp)
   if (allocated(micro%abc2mp))   deallocate (micro%abc2mp)
   if (allocated(micro%regen_aero1_np)) deallocate (micro%regen_aero1_np)
   if (allocated(micro%regen_aero1_mp)) deallocate (micro%regen_aero1_mp)
   if (allocated(micro%regen_aero2_np)) deallocate (micro%regen_aero2_np)
   if (allocated(micro%regen_aero2_mp)) deallocate (micro%regen_aero2_mp)

   if (allocated(micro%immercp)) deallocate (micro%immercp)
   if (allocated(micro%immerdp)) deallocate (micro%immerdp)
   if (allocated(micro%immerrp)) deallocate (micro%immerrp)
   if (allocated(micro%ifnnucp)) deallocate (micro%ifnnucp)

   if (allocated(micro%cnmcp))   deallocate (micro%cnmcp)
   if (allocated(micro%cnmdp))   deallocate (micro%cnmdp)
   if (allocated(micro%cnmrp))   deallocate (micro%cnmrp)
   if (allocated(micro%cnmpp))   deallocate (micro%cnmpp)
   if (allocated(micro%cnmsp))   deallocate (micro%cnmsp)
   if (allocated(micro%cnmap))   deallocate (micro%cnmap)
   if (allocated(micro%cnmgp))   deallocate (micro%cnmgp)
   if (allocated(micro%cnmhp))   deallocate (micro%cnmhp)
   if (allocated(micro%accpaero))deallocate (micro%accpaero)
   if (allocated(micro%pcpraero))deallocate (micro%pcpraero)

   if (allocated(micro%dnmcp))   deallocate (micro%dnmcp)
   if (allocated(micro%dnmdp))   deallocate (micro%dnmdp)
   if (allocated(micro%dnmrp))   deallocate (micro%dnmrp)
   if (allocated(micro%dnmpp))   deallocate (micro%dnmpp)
   if (allocated(micro%dnmsp))   deallocate (micro%dnmsp)
   if (allocated(micro%dnmap))   deallocate (micro%dnmap)
   if (allocated(micro%dnmgp))   deallocate (micro%dnmgp)
   if (allocated(micro%dnmhp))   deallocate (micro%dnmhp)
   if (allocated(micro%accpdust))deallocate (micro%accpdust)
   if (allocated(micro%pcprdust))deallocate (micro%pcprdust)

   if (allocated(micro%dincp))   deallocate (micro%dincp)
   if (allocated(micro%dindp))   deallocate (micro%dindp)
   if (allocated(micro%dinrp))   deallocate (micro%dinrp)
   if (allocated(micro%dinpp))   deallocate (micro%dinpp)
   if (allocated(micro%dinsp))   deallocate (micro%dinsp)
   if (allocated(micro%dinap))   deallocate (micro%dinap)
   if (allocated(micro%dingp))   deallocate (micro%dingp)
   if (allocated(micro%dinhp))   deallocate (micro%dinhp)

   if (allocated(micro%snmcp))   deallocate (micro%snmcp)
   if (allocated(micro%snmdp))   deallocate (micro%snmdp)
   if (allocated(micro%snmrp))   deallocate (micro%snmrp)
   if (allocated(micro%snmpp))   deallocate (micro%snmpp)
   if (allocated(micro%snmsp))   deallocate (micro%snmsp)
   if (allocated(micro%snmap))   deallocate (micro%snmap)
   if (allocated(micro%snmgp))   deallocate (micro%snmgp)
   if (allocated(micro%snmhp))   deallocate (micro%snmhp)
   if (allocated(micro%resol_aero1_mp)) deallocate (micro%resol_aero1_mp)
   if (allocated(micro%resol_aero2_mp)) deallocate (micro%resol_aero2_mp)

   if (allocated(micro%pcpvr))   deallocate (micro%pcpvr)
   if (allocated(micro%pcpvp))   deallocate (micro%pcpvp)
   if (allocated(micro%pcpvs))   deallocate (micro%pcpvs)
   if (allocated(micro%pcpva))   deallocate (micro%pcpva)
   if (allocated(micro%pcpvg))   deallocate (micro%pcpvg)
   if (allocated(micro%pcpvh))   deallocate (micro%pcpvh)
   if (allocated(micro%pcpvd))   deallocate (micro%pcpvd)

   if (allocated(micro%accpr))   deallocate (micro%accpr)
   if (allocated(micro%accpp))   deallocate (micro%accpp)
   if (allocated(micro%accps))   deallocate (micro%accps)
   if (allocated(micro%accpa))   deallocate (micro%accpa)
   if (allocated(micro%accpg))   deallocate (micro%accpg)
   if (allocated(micro%accph))   deallocate (micro%accph)
   if (allocated(micro%accpd))   deallocate (micro%accpd)
   if (allocated(micro%pcprr))   deallocate (micro%pcprr)
   if (allocated(micro%pcprp))   deallocate (micro%pcprp)
   if (allocated(micro%pcprs))   deallocate (micro%pcprs)
   if (allocated(micro%pcpra))   deallocate (micro%pcpra)
   if (allocated(micro%pcprg))   deallocate (micro%pcprg)
   if (allocated(micro%pcprh))   deallocate (micro%pcprh)
   if (allocated(micro%pcprd))   deallocate (micro%pcprd)
   if (allocated(micro%pcpg))    deallocate (micro%pcpg)
   if (allocated(micro%qpcpg))   deallocate (micro%qpcpg)
   if (allocated(micro%dpcpg))   deallocate (micro%dpcpg)

    !COMPUTE AND OUTPUT MICRO BUDGET PROCESSES
    if (allocated(micro%latheatvap))    deallocate (micro%latheatvap)
    if (allocated(micro%latheatfrz))    deallocate (micro%latheatfrz)
    if (allocated(micro%nuccldrt))      deallocate (micro%nuccldrt)
    if (allocated(micro%cld2raint))     deallocate (micro%cld2raint)
    if (allocated(micro%ice2raint))     deallocate (micro%ice2raint)
    if (allocated(micro%nucicert))      deallocate (micro%nucicert)
    if (allocated(micro%vapliqt))       deallocate (micro%vapliqt)
    if (allocated(micro%vapicet))       deallocate (micro%vapicet)
    if (allocated(micro%evapliqt))      deallocate (micro%evapliqt)
    if (allocated(micro%evapicet))      deallocate (micro%evapicet)
    if (allocated(micro%freezingt))     deallocate (micro%freezingt)
    if (allocated(micro%meltingt))      deallocate (micro%meltingt)
    if (allocated(micro%melticet))      deallocate (micro%melticet)
    if (allocated(micro%rimecldt))      deallocate (micro%rimecldt)
    if (allocated(micro%rain2icet))     deallocate (micro%rain2icet)
    if (allocated(micro%aggregatet))    deallocate (micro%aggregatet)
    if (allocated(micro%latheatvapt))   deallocate (micro%latheatvapt)
    if (allocated(micro%latheatfrzt))   deallocate (micro%latheatfrzt)

    if (allocated(micro%inuchomrt))     deallocate (micro%inuchomrt)
    if (allocated(micro%inuccontrt))    deallocate (micro%inuccontrt)
    if (allocated(micro%inucifnrt))     deallocate (micro%inucifnrt)
    if (allocated(micro%inuchazrt))     deallocate (micro%inuchazrt)
    if (allocated(micro%vapcldt))       deallocate (micro%vapcldt)
    if (allocated(micro%vapraint))      deallocate (micro%vapraint)
    if (allocated(micro%vapprist))      deallocate (micro%vapprist)
    if (allocated(micro%vapsnowt))      deallocate (micro%vapsnowt)
    if (allocated(micro%vapaggrt))      deallocate (micro%vapaggrt)
    if (allocated(micro%vapgraut))      deallocate (micro%vapgraut)
    if (allocated(micro%vaphailt))      deallocate (micro%vaphailt)
    if (allocated(micro%vapdrizt))      deallocate (micro%vapdrizt)
    if (allocated(micro%evapcldt))      deallocate (micro%evapcldt)
    if (allocated(micro%evapraint))     deallocate (micro%evapraint)
    if (allocated(micro%evapprist))     deallocate (micro%evapprist)
    if (allocated(micro%evapsnowt))     deallocate (micro%evapsnowt)
    if (allocated(micro%evapaggrt))     deallocate (micro%evapaggrt)
    if (allocated(micro%evapgraut))     deallocate (micro%evapgraut)
    if (allocated(micro%evaphailt))     deallocate (micro%evaphailt)
    if (allocated(micro%evapdrizt))     deallocate (micro%evapdrizt)
    if (allocated(micro%meltprist))     deallocate (micro%meltprist)
    if (allocated(micro%meltsnowt))     deallocate (micro%meltsnowt)
    if (allocated(micro%meltaggrt))     deallocate (micro%meltaggrt)
    if (allocated(micro%meltgraut))     deallocate (micro%meltgraut)
    if (allocated(micro%melthailt))     deallocate (micro%melthailt)
    if (allocated(micro%rimecldsnowt))  deallocate (micro%rimecldsnowt)
    if (allocated(micro%rimecldaggrt))  deallocate (micro%rimecldaggrt)
    if (allocated(micro%rimecldgraut))  deallocate (micro%rimecldgraut)
    if (allocated(micro%rimecldhailt))  deallocate (micro%rimecldhailt)
    if (allocated(micro%rain2prt))      deallocate (micro%rain2prt)
    if (allocated(micro%rain2snt))      deallocate (micro%rain2snt)
    if (allocated(micro%rain2agt))      deallocate (micro%rain2agt)
    if (allocated(micro%rain2grt))      deallocate (micro%rain2grt)
    if (allocated(micro%rain2hat))      deallocate (micro%rain2hat)
    if (allocated(micro%aggrselfprist)) deallocate (micro%aggrselfprist)
    if (allocated(micro%aggrselfsnowt)) deallocate (micro%aggrselfsnowt)
    if (allocated(micro%aggrprissnowt)) deallocate (micro%aggrprissnowt)

    if (allocated(micro%dust1cldrt))        deallocate (micro%dust1cldrt)
    if (allocated(micro%dust2cldrt))        deallocate (micro%dust2cldrt)
    if (allocated(micro%dust1drzrt))        deallocate (micro%dust1drzrt)
    if (allocated(micro%dust2drzrt))        deallocate (micro%dust2drzrt)

return
END SUBROUTINE dealloc_micro

END MODULE mem_micro

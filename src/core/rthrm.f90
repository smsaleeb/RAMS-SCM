!##############################################################################
Subroutine thermo ()

use mem_grid
use mem_basic
use mem_micro
use micphys, only:ithermo,level
use mem_other
use node_mod, only:mzp,mxp,myp

implicit none

if(ithermo == 1) then

if (level .eq. 3) then

   CALL wetthrm3 (mzp,mxp,myp,1,mxp,1,myp  &
     ,basic_g(ngrid)%thp(1,1,1) ,basic_g(ngrid)%theta(1,1,1)  &
     ,basic_g(ngrid)%rtp(1,1,1) ,basic_g(ngrid)%rv   (1,1,1)  &
     ,micro_g(ngrid)%rcp(1,1,1) ,micro_g(ngrid)%rrp  (1,1,1)  &
     ,micro_g(ngrid)%rpp(1,1,1) ,micro_g(ngrid)%rsp  (1,1,1)  &
     ,micro_g(ngrid)%rap(1,1,1) ,micro_g(ngrid)%rgp  (1,1,1)  &
     ,micro_g(ngrid)%rhp(1,1,1) ,micro_g(ngrid)%q6   (1,1,1)  &
     ,micro_g(ngrid)%q7 (1,1,1) ,micro_g(ngrid)%rdp  (1,1,1)  &
     ,basic_g(ngrid)%pi0(1,1,1) ,basic_g(ngrid)%pp   (1,1,1)  &
     ,other_g(ngrid)%pressure(1)                              &
     )

endif

endif !(if ithermo)

return
END SUBROUTINE thermo

!##############################################################################
Subroutine wetthrm3 (m1,m2,m3,ia,iz,ja,jz                   &
   ,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rap,rgp,rhp,q6,q7,rdp  &
   ,pi0,pp                                                  &
   ,pressure                                                &
   )

! This routine calculates theta and rv for "level 3 microphysics"
! given prognosed theta_il, cloud, rain, pristine ice, snow, aggregates,
! graupel, hail, q6, and q7.

use rconstants
use micphys, only:tair,til,qhydm,rliq,rice,jnmb,ipress

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: tcoal,fracliq,tairstr
real, dimension(m1) :: picpi,pressure
real, dimension(m1,m2,m3) :: pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp,rap &
                            ,rgp,rhp,q6,q7,rdp

do j = ja,jz
   do i = ia,iz

      do k = 1,m1
!SaleebySCM:Replace pressure and use this instead of Exner(PI) since
!other models use pressure rather than Exner function (PI0 + PP)
         if(ipress==0)then
          picpi(k) = (pi0(k,i,j) + pp(k,i,j)) * cpi
         elseif(ipress==1)then
          picpi(k) = (pressure(k)*p00i)**rocp
         endif
         tair(k) = theta(k,i,j) * picpi(k)
         til(k) = thp(k,i,j) * picpi(k)
         rliq(k) = 0.
         rice(k) = 0.
      enddo

      if (jnmb(1) .ge. 1) then
         do k = 1,m1
            rliq(k) = rliq(k) + rcp(k,i,j)
         enddo
      endif

      if (jnmb(2) .ge. 1) then
         do k = 1,m1
            rliq(k) = rliq(k) + rrp(k,i,j)
         enddo
      endif

      if (jnmb(3) .ge. 1) then
         do k = 1,m1
            rice(k) = rice(k) + rpp(k,i,j)
         enddo
      endif

      if (jnmb(4) .ge. 1) then
         do k = 1,m1
            rice(k) = rice(k) + rsp(k,i,j)
         enddo
      endif

      if (jnmb(5) .ge. 1) then
         do k = 1,m1
            rice(k) = rice(k) + rap(k,i,j)
         enddo
      endif

      if (jnmb(6) .ge. 1) then
         do k = 1,m1
            CALL qtc (q6(k,i,j),tcoal,fracliq)
            rliq(k) = rliq(k) + rgp(k,i,j) * fracliq
            rice(k) = rice(k) + rgp(k,i,j) * (1. - fracliq)
         enddo
      endif

      if (jnmb(7) .ge. 1) then
         do k = 1,m1
            CALL qtc (q7(k,i,j),tcoal,fracliq)
            rliq(k) = rliq(k) + rhp(k,i,j) * fracliq
            rice(k) = rice(k) + rhp(k,i,j) * (1. - fracliq)
         enddo
      endif

      if (jnmb(8) .ge. 1) then
         do k = 1,m1
            rliq(k) = rliq(k) + rdp(k,i,j)
         enddo
      endif

      do k = 1,m1
         qhydm(k) = alvl * rliq(k) + alvi * rice(k)
         rv(k,i,j) = rtp(k,i,j) - rliq(k) - rice(k)
      enddo

      do k = 1,m1
         if (tair(k) .gt. 253.) then
            tairstr = 0.5 * (til(k)  &
               + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
         else
            tairstr = til(k) * (1. + qhydm(k) * cp253i)
         endif
         theta(k,i,j) = tairstr / picpi(k)
      enddo

   enddo
enddo

return
END SUBROUTINE wetthrm3

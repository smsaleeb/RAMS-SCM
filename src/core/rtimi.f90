!##############################################################################
Subroutine predthp ()

use mem_grid, only:ngrid,dtlt
use node_mod, only:mzp
use mem_basic
use mem_radiate

implicit none


! Step thermodynamic variable from  t  to  t+1 with radiative tendendcy

! SaleebySCM:Here use mzp instead of mxyzp in full RAMS

if (ilwrtyp + iswrtyp > 0) &
 CALL update (mzp,basic_g(ngrid)%thp,radiate_g(ngrid)%fthrdp,dtlt)

return
END SUBROUTINE predthp


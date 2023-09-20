!##############################################################################
Module mem_basic

implicit none

   Type basic_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          up,vp,wp,pp  &
                         ,rv,theta,thp,rtp &
                         ,pi0,dn0

   End Type
   
   type (basic_vars), allocatable :: basic_g(:)
  
Contains

!##############################################################################
Subroutine alloc_basic (basic,n1,n2,n3)

use micphys

implicit none

   type (basic_vars) :: basic
   integer, intent(in) :: n1,n2,n3
   integer :: numpts

! Allocate arrays based on options (if necessary)

      allocate (basic%up(n1,n2,n3))
      allocate (basic%vp(n1,n2,n3))
      allocate (basic%wp(n1,n2,n3))
      allocate (basic%pp(n1,n2,n3))
      allocate (basic%thp(n1,n2,n3))
      allocate (basic%rtp(n1,n2,n3))
      allocate (basic%rv(n1,n2,n3))
      allocate (basic%theta(n1,n2,n3))
      allocate (basic%pi0(n1,n2,n3))
      allocate (basic%dn0(n1,n2,n3))

      !SaleebySCM: Doing this just to zero out data for SCM
      numpts=n1*n2*n3
      call azero (numpts,basic%up)
      call azero (numpts,basic%vp)
      call azero (numpts,basic%wp)
      call azero (numpts,basic%pp)
      call azero (numpts,basic%thp)
      call azero (numpts,basic%rtp)
      call azero (numpts,basic%rv)
      call azero (numpts,basic%theta)
      call azero (numpts,basic%pi0)
      call azero (numpts,basic%dn0)

return
END SUBROUTINE alloc_basic

!##############################################################################
Subroutine dealloc_basic (basic)

implicit none

   type (basic_vars) :: basic
   
   if (allocated(basic%up   ))    deallocate (basic%up   )
   if (allocated(basic%vp   ))    deallocate (basic%vp   )
   if (allocated(basic%wp   ))    deallocate (basic%wp   )
   if (allocated(basic%pp   ))    deallocate (basic%pp   )
   if (allocated(basic%thp  ))    deallocate (basic%thp  )
   if (allocated(basic%rtp  ))    deallocate (basic%rtp  )
   if (allocated(basic%rv   ))    deallocate (basic%rv   )
   if (allocated(basic%theta))    deallocate (basic%theta)
   if (allocated(basic%pi0  ))    deallocate (basic%pi0  )
   if (allocated(basic%dn0  ))    deallocate (basic%dn0  )

return
END SUBROUTINE dealloc_basic

END MODULE mem_basic

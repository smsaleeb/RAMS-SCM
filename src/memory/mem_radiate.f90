!##############################################################################
Module mem_radiate

implicit none

   Type radiate_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          fthrd,fthrdp,bext,swup,swdn,lwup,lwdn
                          
      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          rshort,rlong,rlongup,albedt,cosz,aodt

   End Type
   
   type (radiate_vars), allocatable :: radiate_g(:)
   
   integer :: ilwrtyp,iswrtyp
   real    :: radfrq

Contains

!##############################################################################
Subroutine alloc_radiate (radiate,n1,n2,n3)

implicit none

   type (radiate_vars) :: radiate
   integer, intent(in) :: n1,n2,n3
   integer :: numpts3d,numpts2d

   numpts3d=n1*n2*n3
   numpts2d=n2*n3

! Allocate arrays based on options (if necessary)
      
      if(ilwrtyp+iswrtyp > 0)  then
                         allocate (radiate%fthrd(n1,n2,n3))
                         allocate (radiate%fthrdp(n1,n2,n3))
                         allocate (radiate%rshort(n2,n3))
                         allocate (radiate%rlong(n2,n3))
                         allocate (radiate%rlongup(n2,n3))
                         allocate (radiate%albedt(n2,n3))
                         allocate (radiate%cosz(n2,n3))
                         allocate (radiate%aodt(n2,n3))
        call azero (numpts3d,radiate%fthrd)
        call azero (numpts3d,radiate%fthrdp)
        call azero (numpts2d,radiate%rshort)
        call azero (numpts2d,radiate%rlong)
        call azero (numpts2d,radiate%rlongup)
        call azero (numpts2d,radiate%albedt)
        call azero (numpts2d,radiate%cosz)
        call azero (numpts2d,radiate%aodt)
      endif
      if(ilwrtyp == 3 .or. iswrtyp == 3) then
         allocate (radiate%bext(n1,n2,n3))
         call azero (numpts3d,radiate%bext)
      endif
      if(ilwrtyp == 3)  then
         allocate (radiate%lwup(n1,n2,n3))
         allocate (radiate%lwdn(n1,n2,n3))
         call azero (numpts3d,radiate%lwup)
         call azero (numpts3d,radiate%lwdn)
      endif
      if(iswrtyp == 3)  then
         allocate (radiate%swup(n1,n2,n3))
         allocate (radiate%swdn(n1,n2,n3))
         call azero (numpts3d,radiate%swup)
         call azero (numpts3d,radiate%swdn)
      endif
   
return
END SUBROUTINE alloc_radiate

!##############################################################################
Subroutine dealloc_radiate (radiate)

implicit none

   type (radiate_vars) :: radiate

   if (allocated(radiate%fthrd))    deallocate (radiate%fthrd)
   if (allocated(radiate%fthrdp))   deallocate (radiate%fthrdp)
   if (allocated(radiate%rshort))   deallocate (radiate%rshort)
   if (allocated(radiate%rlong))    deallocate (radiate%rlong)
   if (allocated(radiate%rlongup))  deallocate (radiate%rlongup)
   if (allocated(radiate%albedt))   deallocate (radiate%albedt)
   if (allocated(radiate%cosz))     deallocate (radiate%cosz)
   if (allocated(radiate%aodt))     deallocate (radiate%aodt)
   if (allocated(radiate%bext))     deallocate (radiate%bext)
   if (allocated(radiate%swup))     deallocate (radiate%swup)
   if (allocated(radiate%swdn))     deallocate (radiate%swdn)
   if (allocated(radiate%lwup))     deallocate (radiate%lwup)
   if (allocated(radiate%lwdn))     deallocate (radiate%lwdn)

return
END SUBROUTINE dealloc_radiate

END MODULE mem_radiate

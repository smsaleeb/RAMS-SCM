!##############################################################################
Module mem_leaf

use grid_dims

implicit none

   Type leaf_vars
   
      ! Variables to be dimensioned by (nzg,nxp,nyp,npatch)
   real, allocatable, dimension(:,:,:,:) :: &
                  soil_water,soil_text

      ! Variables to be dimensioned by (nxp,nyp,npatch)
   real, allocatable, dimension(:,:,:) :: &
                  patch_area,leaf_class

   End Type
   
   type (leaf_vars), allocatable :: leaf_g(:)

Contains

!##############################################################################
Subroutine alloc_leaf (leaf,nx,ny,nzg,nzs,np)

implicit none

   type (leaf_vars) :: leaf
   integer, intent(in) :: nx,ny,nzg,nzs,np


! Allocate arrays based on options (if necessary)

   allocate (leaf%soil_water     (nzg,nx,ny,np))
   allocate (leaf%soil_text      (nzg,nx,ny,np))
   allocate (leaf%patch_area   (nx,ny,np))
   allocate (leaf%leaf_class   (nx,ny,np))

   call azero (nzg*nx*ny*np,leaf%soil_water)
   call azero (nzg*nx*ny*np,leaf%soil_text)
   call azero (nx*ny*np,leaf%patch_area)
   call azero (nx*ny*np,leaf%leaf_class)
   
return
END SUBROUTINE alloc_leaf

!##############################################################################
Subroutine dealloc_leaf (leaf)

implicit none

  type (leaf_vars) :: leaf

  if(allocated(leaf%soil_water))      deallocate (leaf%soil_water)
  if(allocated(leaf%soil_text))       deallocate (leaf%soil_text)
  if(allocated(leaf%patch_area))      deallocate (leaf%patch_area)
  if(allocated(leaf%leaf_class))      deallocate (leaf%leaf_class)
 
return
END SUBROUTINE dealloc_leaf

END MODULE mem_leaf

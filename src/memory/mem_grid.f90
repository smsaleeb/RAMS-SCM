!##############################################################################
Module mem_grid

use grid_dims

implicit none

   Type grid_vars
                            
      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                 topt,rtgt,glat,glon
              
   End Type
   
   type (grid_vars), allocatable :: grid_g(:)

   real,    dimension(nzpmax) :: zt,zm,dzt,dzm

   integer :: iprntstmt,iprntstatic,nzg,nzs,npatch,ngrid,ngrids,maxnzp &
             ,iyear1,imonth1,idate1,ihour1,itime1
   real    :: time,ztop,dtlt

   logical :: print_msg

Contains

!##############################################################################
Subroutine alloc_grid (grid,n2,n3)

implicit none

   type (grid_vars) :: grid
   integer, intent(in) :: n2,n3
   integer :: numpts

! Allocate arrays based on options (if necessary)

     numpts=n2*n3

     allocate (grid%topt(n2,n3))
     allocate (grid%rtgt(n2,n3))
     allocate (grid%glat(n2,n3))
     allocate (grid%glon(n2,n3))
     call azero (numpts,grid%topt)
     call azero (numpts,grid%rtgt)
     call azero (numpts,grid%glat)
     call azero (numpts,grid%glon)

return
END SUBROUTINE alloc_grid

!##############################################################################
Subroutine dealloc_grid (grid)

implicit none

   type (grid_vars) :: grid

! Deallocate arrays
      
     if(allocated(grid%topt   ))    deallocate (grid%topt)
     if(allocated(grid%rtgt   ))    deallocate (grid%rtgt)
     if(allocated(grid%glat   ))    deallocate (grid%glat)
     if(allocated(grid%glon   ))    deallocate (grid%glon)
    
return
END SUBROUTINE dealloc_grid

END MODULE mem_grid

!##############################################################################
Module mem_other

use grid_dims

implicit none

   Type other_vars
                            
      ! Variables to be dimensioned by (nzp)
   real, allocatable, dimension(:) :: &
       !Add new pressure variables here since RAMS uses Exner(PI)
       pressure

              
   End Type
   
   type (other_vars), allocatable :: other_g(:)

Contains

!##############################################################################
Subroutine alloc_other (other,n1,n2,n3,npatch)

implicit none

   type (other_vars) :: other
   integer, intent(in) :: n1,n2,n3,npatch

! Allocate arrays based on options (if necessary)
     allocate (other%pressure(n1))
     call azero (n1,other%pressure)
    
return
END SUBROUTINE alloc_other

!##############################################################################
Subroutine dealloc_other (other)

implicit none

   type (other_vars) :: other

! Deallocate arrays
      
if(allocated(other%pressure))   deallocate (other%pressure)
    
return
END SUBROUTINE dealloc_other

END MODULE mem_other

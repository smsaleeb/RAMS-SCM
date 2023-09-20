!##############################################################################
Module grid_dims

implicit none

! This contains very basic specification of grid dimensions and other 
!    parameters that will be used to dimension arrays and allocate memory.

! Set maximum values of parameters:

integer, parameter ::  &
  strl1        = 128   & ! Long character string length
 ,maxgrds      = 8     & ! Max # of grids
 ,nzpmax       = 332     ! Max # of points in z-direction

END MODULE grid_dims


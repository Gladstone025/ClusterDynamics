module prec_mod

  use iso_c_binding, only: dp => c_double, c_long, c_int

implicit none
real(dp), parameter :: Pi = 3.14159265359
real(dp), parameter :: evJ = 1.602176487e-19 
real(dp), parameter :: kb   = 1.38064852e-23

end module prec_mod

module para_Ov

use prec_mod
use cluster_mod

implicit none


!! ----------------------- Parametres Ovcharenko ---------------------- !!
real(dp), parameter :: Temp = 823
real(dp), parameter :: Vat = 1.205e-29
real(dp), parameter :: ww = (48.*Pi**2/Vat**2)**(1./3.)
real(dp), parameter :: gam = 1.
real(dp), parameter :: Evf = 1.77*evJ
real(dp), parameter :: Evm = 1.1*evJ
real(dp), parameter :: Dv = 1.e-6*exp(-Evm/(kb*Temp))
real(dp), parameter :: Cq = 1.e-7

real(dp), parameter :: G1 = 0._dp



CONTAINS


function beta(x)
	implicit none
	real(dp) :: beta, x
	beta = ww*x**(1./3.)*Dv          
end function

function alpha(x)
	implicit none
	real(dp) :: alpha, x, Evb, r
	r = (3*x*Vat/(4*Pi))**(1./3.)
	Evb = Evf - 2.*gam*Vat/r
	alpha = ww*x**(1./3.)*Dv*exp(-Evb/(kb*Temp))
end function

function F_scal(x,C1)
	implicit none
	real(dp) :: F_scal,x,C1
	F_scal = beta(x)*C1 - alpha(x)     
end function

function D_scal(x,C1)
	implicit none
	real(dp) :: D_scal,x,C1
	D_scal = beta(x)*C1 + alpha(x)     
end function


!! ---- En mode vectoriel ---- !!

function betav(x)
	implicit none
	real(dp),dimension(:) :: x
	real(dp),dimension(1:size(x)) :: betav
	betav = ww*x**(1./3.)*Dv          
end function

function alphav(x)
	implicit none
	real(dp),dimension(:) :: x
	real(dp),dimension(1:size(x)) :: alphav, r, Evb
	r = (3*x*Vat/(4*Pi))**(1./3.)
	Evb = Evf - 2.*gam*Vat/r
	alphav = ww*x**(1./3.)*Dv*exp(-Evb/(kb*Temp))         
end function

function F_vec(x,Conc1)
	implicit none
	real(dp),dimension(:) :: x
	real(dp),dimension(1:size(x)):: F_vec
	real(dp) :: Conc1
	F_vec = betav(x)*Conc1 - alphav(x)    
end function

function D_vec(x,Conc1)
	implicit none
	real(dp),dimension(:) :: x
	real(dp),dimension(1:size(x)):: D_vec
	real(dp) :: Conc1
	D_vec = betav(x)*Conc1 + alphav(x)     
end function









end module para_Ov

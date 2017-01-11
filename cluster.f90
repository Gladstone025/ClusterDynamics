module cluster_mod

use prec_mod

implicit none

integer :: Nv, Ni, Ns
integer :: mv, mi, ms

type :: cluster
	real(dp) :: defect
	real(dp) :: solute
	logical :: mobile = .False.
	integer :: ind = 0
end type cluster

interface operator (+)	
		module procedure c_add
end interface
interface operator (-)	
		module procedure c_sub
end interface
	
CONTAINS


function C2I(c)
	implicit none
	type(cluster) :: c
	real(dp) :: C2I, p1, p2
	p1 = c%defect
	p2 = c%solute
	C2I = 1 + Nv + p2*(Nv+Ni+1) + p1
end function


function IsMobile(n,s)
	implicit none
	real(dp) :: n, s
	logical :: IsMobile
	if (n.eq.0 .and. s.eq.1) then
		IsMobile = .True.
	else if (n.eq.1 .and. s.eq.0) then
		IsMobile = .True.
	else if (n.eq.-1 .and. s.eq.0) then
		IsMobile = .True.
	else
		IsMobile = .False.
	end if		
end function


function IsInDomain(c)
	implicit none
	type(cluster) :: c
	integer :: ind, Nmax
	logical :: IsInDomain
	ind = C2I(c)
	Nmax = (1+Nv+Ni)*(1+Ns)
	if ((ind.ge.1) .and. (ind.le.Nmax)) then
		IsInDomain = .True.
	else
		IsInDomain = .False.
	end if
	if ((c%defect.eq.0) .and. (c%solute.eq.0)) then
		IsInDomain = .False.
	end if
end function


type(cluster) function c_add(c1,c2)
	implicit none
	type(cluster), intent(in) :: c1, c2
	c_add%defect = c1%defect + c2%defect
	c_add%solute = c1%solute + c2%solute
	c_add%mobile = IsMobile(c_add%defect,c_add%solute)
	c_add%ind = C2I(c_add)
end function c_add
		
type(cluster) function c_sub(c1,c2)
	implicit none
	type(cluster), intent(in) :: c1, c2
	c_sub%defect = c1%defect - c2%defect
	c_sub%solute = c1%solute - c2%solute
	c_sub%mobile = IsMobile(c_sub%defect,c_sub%solute)
	c_sub%ind = C2I(c_sub)
end function c_sub	


		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Utile ?
function I2C(ind)
	implicit none
	integer :: ind, p1, p2
	type(cluster) :: I2C
	p2 = int(ind/(1+Nv+Ni))
	p1 = ind - p2*(1+Nv+Ni) - (Nv+1)
	I2C%defect = p1
	I2C%solute = p2
end function
		
		
		

end module cluster_mod

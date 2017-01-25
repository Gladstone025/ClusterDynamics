module cluster_mod

use prec_mod

implicit none

integer :: Nv, Ni, Ns
integer :: mv, mi, ms

type :: cluster
	real(dp) :: defect
	real(dp) :: solute
	logical :: mobile = .False.
	real(dp) :: ind = 0
end type cluster

interface operator (+)	
		module procedure c_add
end interface

interface operator (-)	
		module procedure c_sub
end interface

interface operator (==)
		module procedure c_eq
end interface

interface operator (/=)
		module procedure c_ne
end interface

interface assignment (=)
		module procedure c_assign
end interface


type(cluster), parameter :: C_null = cluster(0._dp,0._dp,.False.,0)
	
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
end function


function C_constructor(n,s)
	implicit none
	type(cluster) :: C_constructor
	real(dp) :: n, s
	integer :: ind
	logical :: mob
	C_constructor%defect = n
	C_constructor%solute = s
	C_constructor%mobile = IsMobile(n,s)
	C_constructor%ind = C2I(C_constructor)
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

logical function c_eq(c1,c2)
	implicit none
	type(cluster), intent(in) :: c1, c2
	if ((c1%defect.eq.c2%defect) .and. (c1%solute.eq.c2%solute) .and. &
		&(c1%mobile.eqv.c2%mobile) .and. (c1%ind.eq.c2%ind)) then
		c_eq = .True.
	else
		c_eq = .False.
	end if
end function c_eq	

logical function c_ne(c1,c2)
	implicit none
	type(cluster), intent(in) :: c1, c2
	if (.not.(c_eq(c1,c2))) then
		c_ne = .True.
	else
		c_ne = .False.
	end if
end function c_ne



subroutine c_assign(c1,c2)
	implicit none
	type(cluster), intent(out) :: c1
	type(cluster), intent(in) :: c2
	c1%defect = c2%defect
	c1%solute = c2%solute
	c1%mobile = c2%mobile
	c1%ind = c2%ind
end subroutine c_assign


function I2C(ind)
	implicit none
	integer :: ind
	real(dp) :: p1, p2
	type(cluster) :: I2C
	p2 = int(ind/(1+Nv+Ni))
	p1 = ind - p2*(1+Nv+Ni) - (Nv+1)
	I2C%defect = p1
	I2C%solute = p2
	I2C%mobile = IsMobile(p1,p2)
	I2C%ind = ind
end function
		
		
		

end module cluster_mod

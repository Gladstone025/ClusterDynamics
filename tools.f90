module tools

use prec_mod
use parametres, only: Pi, iloop, jloop, kloop

implicit none

contains 




subroutine output(C,T)
implicit none
	real(dp) :: T
	real(dp), dimension(:) :: C
	character(50) sTime
	character(60) fic

	write( sTime, '(f8.1)' )  T
	write(*,*) sTime

	fic = '../Data_Fe/EKMC_Distribution_'//trim(adjustl(sTime))
	open(unit=1,file=fic,status='replace')
	do jloop = 1, size(C)
		write(1,'(2(E15.6E3))') float(jloop), C(jloop)   
	enddo
	close(1)
end subroutine output


function approx_lineaire(x,ymin,ymax,xmin,xmax)
implicit none
	real(dp):: approx_lineaire,x,ymin,ymax,xmin,xmax
	approx_lineaire = ymin + (ymax-ymin)*(x- xmin)/(xmax-xmin) 
end function

function Noyau(x)
implicit none
	real(dp) :: x
	real(dp) :: Noyau
	Noyau = 1.d0/sqrt(2.d0*Pi)*exp(-x**2/2.d0)
end function

function Noyau_v(x)
implicit none
	real(dp), dimension(:) :: x
	real(dp), dimension(size(x)) :: Noyau_v
	Noyau_v = 1.d0/sqrt(2.d0*Pi)*exp(-x**2/2.d0)
end function

subroutine gen_rand_normal(W, Npart)
implicit none
	integer :: Npart
	real(dp) :: pi
	real(dp), dimension(1:Npart) :: W, U1, U2
	call random_number(U1)
	call random_number(U2)
	pi = 4.*atan(1.)
	W = sqrt(-2.*log(U1))*cos(2.*pi*U2)
end subroutine gen_rand_normal



subroutine tirage_exp(ee,nu)
implicit none
	real(dp) :: ee, nu, u1
	call random_number(u1)
	ee = -(1._dp/nu)*log(u1)
end subroutine

subroutine init_random_seed()
use iso_fortran_env, only: int64
implicit none
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid
	integer(int64) :: t

	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(newunit=un, file="/dev/urandom", access="stream", &
		form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
		read(un) seed
		close(un)
	else
		! Fallback to XOR:ing the current time and pid. The PID is
		! useful in case one launches multiple instances of the same
		! program in parallel.
		call system_clock(t)
		if (t == 0) then
			call date_and_time(values=dt)
			t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
		       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
		       + dt(3) * 24_int64 * 60 * 60 * 1000 &
		       + dt(5) * 60 * 60 * 1000 &
		       + dt(6) * 60 * 1000 + dt(7) * 1000 &
		       + dt(8)
		end if
		pid = getpid()
		t = ieor(t, int(pid, kind(t)))
		do i = 1, n
			seed(i) = lcg(t)
		end do
	end if
	call random_seed(put=seed)
	contains
	! This simple PRNG might not be good enough for real work, but is
	! sufficient for seeding a better PRNG.
	function lcg(s)
		integer :: lcg
		integer(int64) :: s
		if (s == 0) then
			s = 104729
		else
			s = mod(s, 4294967296_int64)
		end if
		s = mod(s * 279470273_int64, 4294967291_int64)
		lcg = int(mod(s, int(huge(0), int64)), kind(0))
	end function lcg
end subroutine init_random_seed


end module tools





module fonctioncvode


use prec_mod
use parametres, only: alpha, beta, iloop, rloop, quasi, Neq
!use data_fun, only: Neq

implicit none
contains

!! ---- fonction F pour solver ---- !!

subroutine FCVFUN(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun_')
	implicit none
	real(dp) :: T, C1
	integer(c_int) :: IER
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	ConcP(1:Neq) = 0.d0
	IER = 0
	C1 = Conc(1)
	do iloop = 2, Neq-1
		rloop = real(iloop,8)
		ConcP(iloop) = beta(rloop-1)*C1*Conc(iloop-1) - (beta(rloop)*C1+alpha(rloop))*Conc(iloop) + alpha(rloop+1)*Conc(iloop+1)
		!ConcP(1) = ConcP(1) - beta(rloop)*Conc(iloop)*C1 + alpha(rloop)*Conc(iloop)
	enddo
	if (quasi) then
		ConcP(1) = 0
	else
		ConcP(1) = -2.*beta(1.d0)*C1*C1 + alpha(2.d0)*Conc(2)
		do iloop = 2, Neq-1
			rloop = real(iloop,8)
			ConcP(1) = ConcP(1) - beta(rloop)*Conc(iloop)*C1 + alpha(rloop)*Conc(iloop)
		end do
	end if
	ConcP(Neq) = beta(rloop-1)*C1*Conc(iloop-1) - (beta(rloop)*C1+alpha(rloop))*Conc(iloop)
end subroutine


end module fonctioncvode

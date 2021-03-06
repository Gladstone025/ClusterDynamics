module parametres


use prec_mod

implicit none
real(dp), parameter :: Pi = 3.14159265359


!! ------------------------- Boundary and model parameters ------------------------- !!
integer, parameter :: N_front = 200
integer, parameter :: N_buff = 100

integer, parameter :: Nf_Inter = 400
integer, parameter :: Nb_Inter = 400

integer, parameter :: Nf_Vac = 200
integer, parameter :: Nb_Vac = 50

integer, parameter :: Nf_EL_Inter = 50
integer, parameter :: Nf_EL_Vac = 50

real(dp), parameter :: R_min = -1.e10

real(dp) :: M_queue = 0 
real(dp), dimension(1) :: Cvac
real(dp) :: CmobVac, CmobInter
logical :: quasi, quasi1
integer :: etape, methode, cas_physique, Nmax
logical :: Coupling_Inter = .False. 
logical :: Coupling_Vac = .False.

!! --------------------------- Parametres CVODE ----------------------------- !!
real(dp) :: T, T0, TF, DeltaT, TF1, TF2, T2
integer(c_int) :: IER
real(dp), dimension(:), allocatable :: C, C_init, C_stotodis, C_init1, C_init2
real(dp), dimension(:), allocatable :: C_Inter, C_Vac, C_Mob, bCn_sto, aCmob_sto
real(dp), dimension(1) :: RPAR
integer(c_long), dimension(1) :: IPAR
real(dp) :: tt, tout
real(dp) :: abstol, reltol

integer, parameter :: isolver = 1 
integer, parameter :: meth = 2
integer, parameter :: itmeth = 2
integer, parameter :: iatol = 1
integer, parameter :: itask = 1
integer(c_long), dimension(21) :: iout
real(dp), dimension(6) :: rout
integer(c_long) :: maxerrfail, nitermax
integer(c_long) :: Neq
integer :: Nv, Ni, mv, mi, Nmaxv, Nmaxi
real(dp), dimension(:,:), allocatable :: Alpha_tab, Beta_tab

!! ----------------------- Parametres stochastiques ------------------------- !!
integer, parameter :: Taille = 250000
real(dp), dimension(Taille) :: Xpart = 0._dp
real(dp), dimension(Taille) :: XpartInter = 0._dp
real(dp), dimension(Taille) :: XpartVac = 0._dp
real(dp), dimension(:), allocatable :: C_sto

real(dp) :: bCnCi_sto, aCi_sto, Si_sto
real(dp) :: bCnCv_sto, aCv_sto, Sv_sto

real(dp) :: N_0 = real(N_front+N_buff/5._dp,8)
real(dp) :: alpha_m = real(N_buff/10._dp,8)
real(dp) :: dt_sto
real(dp) :: MStochastique, MDiscret
real(dp) :: MStoInter, MStoVac, MDisInter, MDisVac

!! ---------------------------- Parametres MPI ------------------------------ !!
integer :: rank, numproc, ierror
real(dp) :: MPI_Cvac
real(dp), dimension(:), allocatable :: MPI_C_Mob
real(dp), dimension(:), allocatable :: MPI_C_stotodis

!! --------------------------- Parametres utiles ---------------------------- !!
real(dp) :: rloop
real(dp), parameter :: h = 0.4_dp

!! ----------------------- Parametres Ovcharenko ---------------------- !!
real(dp), parameter :: evJ = 1.602176487e-19 
real(dp), parameter :: kb   = 1.38064852e-23
real(dp), parameter :: Temp = 823
real(dp), parameter :: Vat = 1.205e-29
real(dp), parameter :: ww = (48.*Pi**2/Vat**2)**(1./3.)
real(dp), parameter :: gam = 1.
real(dp), parameter :: Evf = 1.77*evJ
real(dp), parameter :: Evm = 1.1*evJ
real(dp), parameter :: Dv = 1.e-6*exp(-Evm/(kb*Temp))
real(dp), parameter :: Cq = 1.e-7

real(dp), parameter :: G1 = 0._dp

!! -------------- Paramètres physiques (Fe) irradiation --------------- !!
real(dp), parameter :: a      = 2.87e-10
real(dp), parameter :: V_at   = (a**3)/2.d0
real(dp), parameter :: b      = a
real(dp), parameter :: K_1v   = 1.41e-10 !(3.d0/(8.d0*pi))**(1./3.)*a
real(dp), parameter :: K_2v   = 1.74d0
real(dp), parameter :: s_1v   = 1./3.
real(dp), parameter :: s_2v   = 0.d0
real(dp), parameter :: K_1i   = 5.7248e-11 !sqrt(a/(8.d0*pi*b))*a
real(dp), parameter :: K_2i   = 5.98149 !35.d0*sqrt(b/(8.d0*pi*a))-1.d0
real(dp), parameter :: s_1i   = 1./2.
real(dp), parameter :: s_2i   = 0.15d0
real(dp), parameter :: r_iv   = 0.65e-9
real(dp), parameter :: D_0    = 8.2e-7
real(dp), parameter :: Em_v   = 0.83*evJ
real(dp), parameter :: Em_2v  = 0.62*evJ
real(dp), parameter :: Em_3v  = 0.35*evJ
real(dp), parameter :: Em_4v  = 0.48*evJ
real(dp), parameter :: Em_i   = 0.34*evJ
real(dp), parameter :: Em_2i  = 0.42*evJ
real(dp), parameter :: Em_3i  = 0.43*evJ
real(dp), parameter :: Ef_v   = 2.12*evJ
real(dp), parameter :: Ef_i   = 3.77*evJ
real(dp), parameter :: Eb_2v  = 0.30*evJ
real(dp), parameter :: Eb_3v  = 0.37*evJ
real(dp), parameter :: Eb_4v  = 0.62*evJ
real(dp), parameter :: Eb_2i  = 0.80*evJ
real(dp), parameter :: Eb_3i  = 0.92*evJ
real(dp), parameter :: Eb_4i  = 1.64*evJ
real(dp), parameter :: rho_c  = a/2.d0
real(dp), parameter :: Z_v    = 1.d0
real(dp), parameter :: Z_i    = 1.1d0
real(dp), parameter :: mu     = 82._dp
real(dp), parameter :: nu     = 0.29_dp
real(dp), parameter :: TempFe = 563
real(dp), parameter :: T_max  = 0.1/(3.87e-9)
real(dp), parameter :: rho_d  = 1.e12
real(dp), parameter :: l_gb   = 100.e-6
real(dp), parameter :: k_b    = 1.38064852e-23
real(dp), parameter :: G_v    = 2.322e-10/V_at
real(dp), parameter :: G_i    = 1.1606517e-9/V_at
real(dp), parameter :: G_8v   = 1.161e-10/V_at
real(dp), parameter :: G_4i   = 8.7075e-14/V_at

real(dp), parameter :: D_v    = D_0*exp(-Em_v/(k_b*TempFe))
real(dp), parameter :: D_i    = D_0*exp(-Em_i/(k_b*TempFe))
real(dp), parameter :: Cv_eq  = 1./V_at*exp(-Ef_v/(k_b*TempFe))
real(dp), parameter :: Ci_eq  = 1./V_at*exp(-Ef_i/(k_b*TempFe))


contains
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


!! ---- Lacunes/interstitiels ---- !!

function r_n(n)
	implicit none
	real(dp) :: n, r_n
	if (n < 0) then             ! -- cas lacune
		r_n = K_1v*((abs(n)**s_1v-1./2.)+K_2v*(abs(n)**s_2v-1./2.))
	else						! -- cas interstitiel
		r_n = K_1i*((n**s_1i-1./2.)+K_2i*(n**s_2i-1./2.))
	endif
end function

function Z_nm(n, m)
	implicit none
	real(dp) :: n, m, Z_nm
	if (n > 0 .and. m > 0) then
		Z_nm = Z_i
	else
		Z_nm = Z_v
	endif
end function

function D_m(m)
	implicit none
	real(dp) :: m, D_m
	D_m = 0._dp
	select case (int(m))
		case(-4)
			D_m = D_0*exp(-Em_4v/(k_b*TempFe))
		case(-3)
			D_m = D_0*exp(-Em_3v/(k_b*TempFe))
		case(-2)
			D_m = D_0*exp(-Em_2v/(k_b*TempFe))
		case(-1)
			D_m = D_0*exp(-Em_v/(k_b*TempFe))
		case(1)
			D_m = D_0*exp(-Em_i/(k_b*TempFe))
		case(2) 
			D_m = D_0*exp(-Em_2i/(k_b*TempFe))
		case(3)
			D_m = D_0*exp(-Em_3i/(k_b*TempFe))
	end select
end function

function beta_nm(n,m)
	implicit none
	real(dp) :: n, m, beta_nm
	if ((n.eq.0) .or. (m.eq.0)) then
		beta_nm = 0._dp
	elseif ((n*m < 0) .and. (Z_v*(r_n(n)+r_n(m)) < r_iv)) then
		beta_nm = 4.*pi*r_iv*D_m(m)
	else
		if (m < 0) then
			beta_nm = 4.*pi*Z_v*(r_n(n)+r_n(m))*D_m(m)
		else
			beta_nm = 4.*pi*Z_nm(n,m)*(r_n(n)+r_n(m))*D_m(m)
		endif
	endif
end function

! -- liaison lacune - amas lacunaire
function Eb_vv(n)                  
	implicit none
	real(dp) :: n, Eb_vv
	if (n .eq. -2) then
		Eb_vv = Eb_2v
	elseif (n .eq. -3) then
		Eb_vv = Eb_3v
	elseif (n .eq. -4) then
		Eb_vv = Eb_4v 
	else 
		Eb_vv = Ef_v + (Eb_2v - Ef_v)/(2.d0**(2./3.)-1.d0)*(abs(n)**(2./3.) - abs(n+1)**(2./3.))	
	endif
end function

! -- liaison interstitiel - amas lacunaire
function Eb_iv(n)                  
	implicit none
	real(dp) :: n, Eb_iv
	Eb_iv = Ef_v + Ef_i - Eb_vv(n-1)
end function

! -- liaison interstitiel - amas interstitiel
function Eb_ii(n)                  
	implicit none
	real(dp) :: n, Eb_ii, rn, rn_1
	if (n .eq. 2) then
		Eb_ii = Eb_2i
	elseif (n .eq. 3) then
		Eb_ii = Eb_3i
	elseif (n .eq. 4) then
		Eb_ii = Eb_4i
	else
		rn = sqrt(n*V_at/(Pi*b))
		rn_1 = sqrt((n-1)*V_at/(Pi*b))
		Eb_ii = Ef_i + 1e9*(((rn_1*mu*b*b)/(2*(1-nu)))*(log(4*rn_1/rho_c)-1) - ((rn*mu*b*b)/(2*(1-nu)))*(log(4*rn/rho_c)-1)) 
	endif
end function

! -- liaison lacune - amas interstitiel
function Eb_vi(n)                  
	implicit none
	real(dp) :: n, Eb_vi
	Eb_vi = Ef_i + Ef_v - Eb_ii(n+1)
end function


function alpha_nm(n, m)
	implicit none
	real(dp) :: n, m, alpha_nm
	! cas non physiques
	if (n.eq.m .or. n.eq.0 .or. m.eq.0) then
		alpha_nm = 0._dp
	! cas lacune-lacune
	else if (n < 0 .and. m.eq.(-1)) then 
		alpha_nm = beta_nm(n-m,m)/V_at*exp(-Eb_vv(n)/(k_b*TempFe))
	! cas lacune-interstitiel
	else if (n > 0 .and. m.eq.(-1)) then
		alpha_nm = beta_nm(n-m,m)/V_at*exp(-Eb_vi(n)/(k_b*TempFe))
	! cas interstitiel-interstitiel
	else if (n > 0 .and. m.eq.(1)) then
		alpha_nm = beta_nm(n-m,m)/V_at*exp(-Eb_ii(n)/(k_b*TempFe))
	! cas interstitiel-lacune
	else if (n < 0 .and. m.eq.(1)) then 
		alpha_nm = beta_nm(n-m,m)/V_at*exp(-Eb_iv(n)/(k_b*TempFe))
	else
		alpha_nm = 0._dp
	end if
end function




function tab_fe(pos,debut)
	implicit none
	integer :: tab_fe, pos, debut
	tab_fe = 1 + debut + pos
end function


function F_fe(x,ConcMob)
	implicit none
	real(dp) :: x
	real(dp) :: F_fe
	real(dp), dimension(:) :: ConcMob
	integer :: mloop
	F_fe = 0._dp
	do mloop = -mv, mi
		rloop = real(mloop,8)
		F_fe = F_fe + rloop*(beta_nm(x,rloop)*ConcMob(tab_fe(mloop,mv))-alpha_nm(x,rloop)) 
	end do
end function

function D_fe(x,ConcMob)
	implicit none
	real(dp) :: x
	real(dp) :: D_fe
	real(dp), dimension(:) :: ConcMob
	integer :: mloop
	D_fe = 0._dp
	do mloop = -mv, mi
		rloop = real(mloop,8)
		D_fe = D_fe + 0.5_dp*rloop*rloop*(beta_nm(x,rloop)*ConcMob(tab_fe(mloop,mv)) + alpha_nm(x,rloop))
	end do
end function





end module parametres




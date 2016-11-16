module fonctioncvode


use prec_mod
use parametres
use tools


implicit none
contains

!! ---- fonction F pour solver ---- !!

subroutine FCVFUN(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun_')
	implicit none
	real(dp) :: T, Conc1, Concvac, Concinter
	real(dp) :: bCnCi, bCnCv, aCi, aCv, Sv, Si
	integer(c_int) :: IER
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	ConcP(1:Neq) = 0.d0
	IER = 0
	if (cas_physique.eq.1) then 
		Conc1 = Conc(1)
		do iloop = 2, Neq-1
			rloop = real(iloop,8)
			ConcP(iloop) = beta(rloop-1)*Conc1*Conc(iloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(iloop) + alpha(rloop+1)*Conc(iloop+1)
		enddo
		if (quasi) then
			ConcP(1) = 0
		else
			ConcP(1) = -2.*beta(1.d0)*Conc1*Conc1 + alpha(2.d0)*Conc(2)
			do iloop = 2, Neq-1
				rloop = real(iloop,8)
				ConcP(1) = ConcP(1) - beta(rloop)*Conc(iloop)*Conc1 + alpha(rloop)*Conc(iloop)
			end do
		end if
		ConcP(Neq) = beta(rloop-1)*Conc1*Conc(iloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(iloop)
	else if (cas_physique.eq.2) then
		bCnCi = 0._dp
		bCnCv = 0._dp
		aCi = 0._dp
		aCv = 0._dp
		Si = 0._dp
		Sv = 0._dp
		ConcP(1:Neq) = 0.d0
		IER = 0
		Concvac = Conc(tab_fe(-1,Nv))
		Concinter = Conc(tab_fe(1,Nv))
		do iloop = -Nv+1,-2
			ConcP(tab_fe(iloop,Nv)) = Beta_tab(iloop-1,1)*Conc(tab_fe(iloop-1,Nv))*Concinter - Alpha_tab(iloop,1)*Conc(tab_fe(iloop,Nv)) - &
									&Beta_tab(iloop,1)*Conc(tab_fe(iloop,Nv))*Concinter + Alpha_tab(iloop+1,1)*Conc(tab_fe(iloop+1,Nv)) + &
									&Beta_tab(iloop+1,-1)*Conc(tab_fe(iloop+1,Nv))*Concvac - Alpha_tab(iloop,-1)*Conc(tab_fe(iloop,Nv)) - &
									&Beta_tab(iloop,-1)*Conc(tab_fe(iloop,Nv))*Concvac + Alpha_tab(iloop-1,-1)*Conc(tab_fe(iloop-1,Nv))
			if (iloop .eq. -8) then
				ConcP(tab_fe(iloop,Nv)) = ConcP(tab_fe(iloop,Nv)) + G_8v
			end if
		end do
		do iloop = 2, Ni-1
			ConcP(tab_fe(iloop,Nv)) = Beta_tab(iloop-1,1)*Conc(tab_fe(iloop-1,Nv))*Concinter - Alpha_tab(iloop,1)*Conc(tab_fe(iloop,Nv)) - &
									&Beta_tab(iloop,1)*Conc(tab_fe(iloop,Nv))*Concinter + Alpha_tab(iloop+1,1)*Conc(tab_fe(iloop+1,Nv)) + &
									&Beta_tab(iloop+1,-1)*Conc(tab_fe(iloop+1,Nv))*Concvac - Alpha_tab(iloop,-1)*Conc(tab_fe(iloop,Nv)) - &
									&Beta_tab(iloop,-1)*Conc(tab_fe(iloop,Nv))*Concvac + Alpha_tab(iloop-1,-1)*Conc(tab_fe(iloop-1,Nv))
			if (iloop .eq. 4) then
				ConcP(tab_fe(iloop,Nv)) = ConcP(tab_fe(iloop,Nv)) + G_4i
			endif
		end do
		ConcP(tab_fe(Ni,Nv)) = Beta_tab(Ni-1,1)*Conc(tab_fe(Ni-1,Nv))*Concinter - Alpha_tab(Ni,1)*Conc(tab_fe(Ni,Nv)) - &
									&Beta_tab(Ni,1)*Conc(tab_fe(Ni,Nv))*Concinter - &
									&Alpha_tab(Ni,-1)*Conc(tab_fe(Ni,Nv)) - &
									&Beta_tab(Ni,-1)*Conc(tab_fe(Ni,Nv))*Concvac + Alpha_tab(Ni-1,-1)*Conc(tab_fe(Ni-1,Nv))
		ConcP(tab_fe(-Nv,Nv)) = 0._dp - Alpha_tab(-Nv,1)*Conc(tab_fe(-Nv,Nv)) - &
									&Beta_tab(-Nv,1)*Conc(tab_fe(-Nv,Nv))*Concinter + Alpha_tab(-Nv+1,1)*Conc(tab_fe(-Nv+1,Nv)) + &
									&Beta_tab(-Nv+1,-1)*Conc(tab_fe(-Nv+1,Nv))*Concvac - Alpha_tab(-Nv,-1)*Conc(tab_fe(-Nv,Nv)) - &
									&Beta_tab(-Nv,-1)*Conc(tab_fe(-Nv,Nv))*Concvac
		if (quasi) then
			ConcP(tab_fe(1,Nv)) = 0
			ConcP(tab_fe(-1,Nv)) = 0
		else
			! somme sur les amas immobiles
			do iloop = 2, Ni - 1
				bCnCi = bCnCi - Beta_tab(iloop,1)*Conc(tab_fe(iloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(iloop,-1)*Conc(tab_fe(iloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(iloop+1,1)*Conc(tab_fe(iloop+1,Nv))
				aCv = aCv + Alpha_tab(iloop-1,-1)*Conc(tab_fe(iloop-1,Nv))
				Si = Si + Beta_tab(iloop,1)*Conc(tab_fe(iloop,Nv))
				Sv = Sv + Beta_tab(iloop,-1)*Conc(tab_fe(iloop,Nv))
			enddo
			do iloop = -Nv + 2, -2 
				bCnCi = bCnCi - Beta_tab(iloop,1)*Conc(tab_fe(iloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(iloop,-1)*Conc(tab_fe(iloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(iloop+1,1)*Conc(tab_fe(iloop+1,Nv))
				aCv = aCv + Alpha_tab(iloop-1,-1)*Conc(tab_fe(iloop-1,Nv))
				Si = Si + Beta_tab(iloop,1)*Conc(tab_fe(iloop,Nv))
				Sv = Sv + Beta_tab(iloop,-1)*Conc(tab_fe(iloop,Nv))
			enddo
			! contribution amas mobiles
			bCnCi = bCnCi - 2._dp*Beta_tab(1,1)*Conc(tab_fe(1,Nv))*Conc(tab_fe(1,Nv)) - Beta_tab(-1,1)*Conc(tab_fe(-1,Nv))*Concinter + &
					&Beta_tab(2,-1)*Conc(tab_fe(2,Nv))*Conc(tab_fe(-1,Nv)) - Beta_tab(1,-1)*Conc(tab_fe(1,Nv))*Conc(tab_fe(-1,Nv))
			aCi = aCi + 2._dp*Alpha_tab(2,1)*Conc(tab_fe(2,Nv)) - Alpha_tab(1,-1)*Conc(tab_fe(1,Nv)) 
			bCnCv = bCnCv - 2._dp*Beta_tab(-1,-1)*Conc(tab_fe(-1,Nv))*Conc(tab_fe(-1,Nv)) - Beta_tab(1,-1)*Conc(tab_fe(1,Nv))*Concvac + &
					&Beta_tab(-2,1)*Conc(tab_fe(-2,Nv))*Conc(tab_fe(1,Nv)) - Beta_tab(-1,1)*Conc(tab_fe(-1,Nv))*Conc(tab_fe(1,Nv))
			aCv = aCv - Alpha_tab(-1,1)*Conc(tab_fe(-1,Nv)) + 2._dp*Alpha_tab(-2,-1)*Conc(tab_fe(-2,Nv))
			Si = Si + Beta_tab(-1,1)*(Conc(tab_fe(-1,Nv))+Conc(tab_fe(1,Nv)))
			Sv = Sv + Beta_tab(1,-1)*(Conc(tab_fe(1,Nv))+Conc(tab_fe(-1,Nv)))
			! contribution totale
			ConcP(tab_fe(1,Nv)) = G_i + bCnCi + aCi - Z_i*rho_d*D_i*(Concinter-Ci_eq) - &
										&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Concinter-Ci_eq)
			ConcP(tab_fe(-1,Nv)) = G_v + bCnCv + aCv - Z_v*rho_d*D_v*(Concvac-Cv_eq) - &
										&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Concvac-Cv_eq)
		end if
	else
		print *, "Cas physique inexistant"
		IER = 1
	end if
end subroutine


end module fonctioncvode

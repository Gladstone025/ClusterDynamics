module fonctioncvode


use prec_mod
use parametres
use tools


implicit none
contains

!! ---- fonction F pour solver ---- !!

subroutine FCVFUN1(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun1_')
	implicit none
	real(dp) :: T, Conc1, Concvac, Concinter
	real(dp) :: bCnCi, bCnCv, aCi, aCv, Sv, Si
	integer(c_int) :: IER
	integer :: nloop
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	ConcP(1:Neq) = 0.d0
	IER = 0
	if (cas_physique.eq.1) then 
		Conc1 = Conc(1)
		do nloop = 2, Neq-1
			rloop = real(nloop,8)
			ConcP(nloop) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop) + alpha(rloop+1)*Conc(nloop+1)
		enddo
		if (quasi) then
			ConcP(1) = 0
		else
			ConcP(1) = -2.*beta(1.d0)*Conc1*Conc1 + alpha(2.d0)*Conc(2)
			do nloop = 2, Neq-1
				rloop = real(nloop,8)
				ConcP(1) = ConcP(1) - beta(rloop)*Conc(nloop)*Conc1 + alpha(rloop)*Conc(nloop)
			end do
		end if
		ConcP(Neq) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop)
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
		do nloop = -Nv+1,-2
			ConcP(tab_fe(nloop,Nv)) = Beta_tab(nloop-1,1)*Conc(tab_fe(nloop-1,Nv))*Concinter - Alpha_tab(nloop,1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv)) + &
									&Beta_tab(nloop+1,-1)*Conc(tab_fe(nloop+1,Nv))*Concvac - Alpha_tab(nloop,-1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
			if (nloop .eq. -8) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_8v
			end if
		end do
		do nloop = 2, Ni-1
			ConcP(tab_fe(nloop,Nv)) = Beta_tab(nloop-1,1)*Conc(tab_fe(nloop-1,Nv))*Concinter - Alpha_tab(nloop,1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv)) + &
									&Beta_tab(nloop+1,-1)*Conc(tab_fe(nloop+1,Nv))*Concvac - Alpha_tab(nloop,-1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
			if (nloop .eq. 4) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_4i
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
			! somme sur les amas immobiles deterministes
			do nloop = 2, Nf_Inter - 1
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
			enddo
			do nloop = -Nf_Vac + 2, -2 
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
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
			! contribution totale (amas stochastiques, forces de puits et creation de defauts)
			bCnCi = bCnCi + bCnCi_sto
			bCnCv = bCnCv + bCnCv_sto
			aCi = aCi + aCi_sto
			aCv = aCv + aCv_sto
			Si = Si + Si_sto
			Sv = Sv + Sv_sto			
			ConcP(tab_fe(1,Nv)) = G_i + bCnCi + aCi - Z_i*rho_d*D_i*(Concinter-Ci_eq) - &
										&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Concinter-Ci_eq)
			ConcP(tab_fe(-1,Nv)) = G_v + bCnCv + aCv - Z_v*rho_d*D_v*(Concvac-Cv_eq) - &
										&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Concvac-Cv_eq)
		else
			! somme sur les amas immobiles
			do nloop = 2, Ni - 1
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
			enddo
			do nloop = -Nv + 2, -2 
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
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
			!do nloop = -20, 20
			!	print *, nloop, ConcP(tab_fe(nloop,Nv))
			!end do
		end if
	else
		print *, "Cas physique inexistant"
		IER = 1
	end if
end subroutine

subroutine FCVFUN2(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun2_')
	implicit none
	real(dp) :: T, Conc1, Concvac, Concinter
	real(dp) :: bCnCi, bCnCv, aCi, aCv, Sv, Si
	integer(c_int) :: IER
	integer :: nloop
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	ConcP(1:Neq) = 0.d0
	IER = 0
	if (cas_physique.eq.1) then 
		Conc1 = Conc(1)
		do nloop = 2, Neq-1
			rloop = real(nloop,8)
			ConcP(nloop) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop) + alpha(rloop+1)*Conc(nloop+1)
		enddo
		if (quasi) then
			ConcP(1) = 0
		else
			ConcP(1) = -2.*beta(1.d0)*Conc1*Conc1 + alpha(2.d0)*Conc(2)
			do nloop = 2, Neq-1
				rloop = real(nloop,8)
				ConcP(1) = ConcP(1) - beta(rloop)*Conc(nloop)*Conc1 + alpha(rloop)*Conc(nloop)
			end do
		end if
		ConcP(Neq) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop)
	else if (cas_physique.eq.2) then
		bCnCi = 0._dp
		bCnCv = 0._dp
		aCi = 0._dp
		aCv = 0._dp
		Si = 0._dp
		Sv = 0._dp
		ConcP(1:Neq) = 0.d0
		IER = 0
		if (quasi1) then 
			Concvac = CmobVac !Conc(tab_fe(-1,Nv))
			Concinter = CmobInter !Conc(tab_fe(1,Nv))
		else
			Concvac = Conc(tab_fe(-1,Nv))
			Concinter = Conc(tab_fe(1,Nv))
		end if
		do nloop = -Nv+1,-2
			ConcP(tab_fe(nloop,Nv)) = Beta_tab(nloop-1,1)*Conc(tab_fe(nloop-1,Nv))*Concinter - Alpha_tab(nloop,1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv)) + &
									&Beta_tab(nloop+1,-1)*Conc(tab_fe(nloop+1,Nv))*Concvac - Alpha_tab(nloop,-1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
			if ((nloop .eq. -8).and.(.not.quasi1)) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_8v
			end if
		end do
		do nloop = 2, Ni-1
			ConcP(tab_fe(nloop,Nv)) = Beta_tab(nloop-1,1)*Conc(tab_fe(nloop-1,Nv))*Concinter - Alpha_tab(nloop,1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv)) + &
									&Beta_tab(nloop+1,-1)*Conc(tab_fe(nloop+1,Nv))*Concvac - Alpha_tab(nloop,-1)*Conc(tab_fe(nloop,Nv)) - &
									&Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
			if ((nloop .eq. 4).and.(.not.quasi1)) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_4i
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
			! somme sur les amas immobiles
			do nloop = 2, Nf_Inter - 1
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
			enddo
			do nloop = -Nf_Vac + 2, -2 
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
			enddo
			! Somme sur les amas stochastiques
			!if (Coupling_Inter) then
			!	print *, "FCVfun Coupling Inter"
			!	do nloop = 1, Taille
			!		bCnCi = bCnCi - MStoInter/Taille*beta_nm(XpartInter(nloop),1._dp)*Concinter
			!		bCnCv = bCnCv - MStoInter/Taille*beta_nm(XpartInter(nloop),-1._dp)*Concvac
			!		aCi = aCi + MStoInter/Taille*alpha_nm(XpartInter(nloop)+1._dp,1._dp)
			!		aCv = aCv + MStoInter/Taille*alpha_nm(XpartInter(nloop)-1._dp,-1._dp)
			!		Si = Si + MStoInter/Taille*beta_nm(XpartInter(nloop),1._dp)
			!		Sv = Sv + MStoInter/Taille*beta_nm(XpartInter(nloop),-1._dp)
			!	enddo
			!end if
			!if (Coupling_Vac) then
			!	do nloop = 1, Taille
			!		bCnCi = bCnCi - MStoVac/Taille*beta_nm(XpartVac(nloop),1._dp)*Concinter
			!		bCnCv = bCnCv - MStoVac/Taille*beta_nm(XpartVac(nloop),-1._dp)*Concvac
			!		aCi = aCi + MStoVac/Taille*alpha_nm(XpartVac(nloop)+1._dp,1._dp)
			!		aCv = aCv + MStoVac/Taille*alpha_nm(XpartVac(nloop)-1._dp,-1._dp)
			!		Si = Si + MStoVac/Taille*beta_nm(XpartVac(nloop),1._dp)
			!		Sv = Sv + MStoVac/Taille*beta_nm(XpartVac(nloop),-1._dp)
			!	enddo
			!end if
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
			bCnCi = bCnCi + bCnCi_sto
			bCnCv = bCnCv + bCnCv_sto
			aCi = aCi + aCi_sto
			aCv = aCv + aCv_sto
			Si = Si + Si_sto
			Sv = Sv + Sv_sto			
			ConcP(tab_fe(1,Nv)) = G_i + bCnCi + aCi - Z_i*rho_d*D_i*(Concinter-Ci_eq) - &
										&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Concinter-Ci_eq)
			ConcP(tab_fe(-1,Nv)) = G_v + bCnCv + aCv - Z_v*rho_d*D_v*(Concvac-Cv_eq) - &
										&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Concvac-Cv_eq)
			if (quasi1) then
				ConcP(tab_fe(1,Nv)) = 0
				ConcP(tab_fe(-1,Nv)) = 0
			end if
		else
			! somme sur les amas immobiles
			do nloop = 2, Ni - 1
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
			enddo
			do nloop = -Nv + 2, -2 
				bCnCi = bCnCi - Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))*Concinter
				bCnCv = bCnCv - Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))*Concvac
				aCi = aCi + Alpha_tab(nloop+1,1)*Conc(tab_fe(nloop+1,Nv))
				aCv = aCv + Alpha_tab(nloop-1,-1)*Conc(tab_fe(nloop-1,Nv))
				Si = Si + Beta_tab(nloop,1)*Conc(tab_fe(nloop,Nv))
				Sv = Sv + Beta_tab(nloop,-1)*Conc(tab_fe(nloop,Nv))
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


! ---- Amas mobiles ---- !

subroutine FCVFUN3(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun3_')
	implicit none
	real(dp) :: T, Conc1, Concvac, Concinter
	real(dp) :: bCnCi, bCnCv, aCi, aCv, Sv, Si
	integer(c_int) :: IER
	integer :: mloop, nloop
	real(dp), dimension(-mv:mi) :: Cmob, bCnCmob, aCmob
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	ConcP(1:Neq) = 0.d0
	IER = 0
	if (cas_physique.eq.1) then 
		Conc1 = Conc(1)
		do nloop = 2, Neq-1
			rloop = real(nloop,8)
			ConcP(nloop) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop) + alpha(rloop+1)*Conc(nloop+1)
		enddo
		if (quasi) then
			ConcP(1) = 0
		else
			ConcP(1) = -2.*beta(1.d0)*Conc1*Conc1 + alpha(2.d0)*Conc(2)
			do nloop = 2, Neq-1
				rloop = real(nloop,8)
				ConcP(1) = ConcP(1) - beta(rloop)*Conc(nloop)*Conc1 + alpha(rloop)*Conc(nloop)
			end do
		end if
		ConcP(Neq) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop)
	else if (cas_physique.eq.2) then
		Cmob = Conc(tab_fe(-mv,Nv):tab_fe(mi,Nv))
		bCnCmob = 0._dp
		aCmob = 0._dp
		Si = 0._dp
		Sv = 0._dp
		ConcP(1:Neq) = 0.d0
		IER = 0
		do nloop = -Nv+1,-mv-1
			ConcP(tab_fe(nloop,Nv)) = 0._dp
			if (nloop .eq. -8) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_8v
			end if
			do mloop = -mv, mi
				if (((nloop+mloop) > -Nv) .and. ((nloop-mloop) > -Nv)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
											& (Beta_tab(nloop,mloop)*Cmob(mloop) + Alpha_tab(nloop,mloop))*Conc(tab_fe(nloop,Nv)) + &
											& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv))
				end if
			end do			
		end do
		do nloop = mi+1, Ni-1
			ConcP(tab_fe(nloop,Nv)) = 0._dp
			if (nloop .eq. 4) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_4i
			endif
			do mloop = -mv, mi
				if (((nloop+mloop) < Ni) .and. ((nloop-mloop) < Ni)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
											& (Beta_tab(nloop,mloop)*Cmob(mloop) + Alpha_tab(nloop,mloop))*Conc(tab_fe(nloop,Nv)) + &
											& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv))
				end if
			end do				
		end do
		! Gestion des amas mobiles
		if (quasi) then
			! somme sur les amas immobiles deterministes
			do nloop = -mv,mi
				do mloop = mi+1, Nf_Inter-mi
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
				do mloop = -Nf_Vac+mv, -mv-1 
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution amas mobiles
			do nloop = -mv, mi
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) 
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution totale (amas stochastiques, forces de puits et creation de defauts)
			do nloop = -mv, mi
				bCnCmob(nloop) = bCnCmob(nloop) + bCn_sto(nloop)*Cmob(nloop)
				aCmob(nloop) = aCmob(nloop) + aCmob_sto(nloop) 
			end do
			Si = Si + Si_sto
			Sv = Sv + Sv_sto	
			! Finalement on calcule la valeur de dC/dt
			do nloop = -mv, mi
				ConcP(tab_fe(nloop,Nv)) = bCnCmob(nloop) + aCmob(nloop)
				if (nloop.eq.(1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_i - Z_i*rho_d*D_i*(Cmob(nloop)-Ci_eq) - &
											&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Cmob(nloop)-Ci_eq)
				else if (nloop.eq.(-1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_v - Z_v*rho_d*D_v*(Cmob(nloop)-Cv_eq) - &
											&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Cmob(nloop)-Cv_eq)
				end if
			end do		
		else
			! contribution des amas immobiles
			do nloop = -mv,mi
				do mloop = mi+1, Ni-mi
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
				do mloop = -Nv+mv, -mv-1 
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution amas mobiles
			do nloop = -mv,-1
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) 
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			do nloop = 1, mi
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) 
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! Finalement on calcule la valeur de dC/dt
			do nloop = -mv, mi
				ConcP(tab_fe(nloop,Nv)) = bCnCmob(nloop) + aCmob(nloop)
				if (nloop.eq.(1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_i - Z_i*rho_d*D_i*(Cmob(nloop)-Ci_eq) - &
											&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Cmob(nloop)-Ci_eq)
				else if (nloop.eq.(-1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_v - Z_v*rho_d*D_v*(Cmob(nloop)-Cv_eq) - &
											&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Cmob(nloop)-Cv_eq)
				end if
			end do	
		end if
		!do nloop = -20, 20
		!	print *, nloop, ConcP(tab_fe(nloop,Nv))
		!end do
		ConcP(tab_fe(0,Nv)) = 0._dp
	else
		print *, "Cas physique inexistant"
		IER = 1
	end if
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




subroutine FCVFUN(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun_')
	implicit none
	real(dp) :: T, Conc1, Concvac, Concinter
	real(dp) :: bCnCi, bCnCv, aCi, aCv, Sv, Si
	integer(c_int) :: IER
	integer :: nloop, mloop
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	real(dp), dimension(-mv:mi) :: Cmob, bCnCmob, aCmob
	ConcP(1:Neq) = 0.d0
	IER = 0
	if (cas_physique.eq.1) then 
		Conc1 = Conc(1)
		do nloop = 2, Neq-1
			rloop = real(nloop,8)
			ConcP(nloop) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop) + alpha(rloop+1)*Conc(nloop+1)
		enddo
		if (quasi) then
			ConcP(1) = 0
		else
			ConcP(1) = -2.*beta(1.d0)*Conc1*Conc1 + alpha(2.d0)*Conc(2)
			do nloop = 2, Neq-1
				rloop = real(nloop,8)
				ConcP(1) = ConcP(1) - beta(rloop)*Conc(nloop)*Conc1 + alpha(rloop)*Conc(nloop)
			end do
		end if
		ConcP(Neq) = beta(rloop-1)*Conc1*Conc(nloop-1) - (beta(rloop)*Conc1+alpha(rloop))*Conc(nloop)
	else if (cas_physique.eq.2) then
		bCnCi = 0._dp
		bCnCv = 0._dp
		aCi = 0._dp
		aCv = 0._dp
		Si = 0._dp
		Sv = 0._dp
		bCnCmob = 0._dp
		aCmob = 0._dp
		ConcP(1:Neq) = 0.d0
		IER = 0
		Concvac = Conc(tab_fe(-1,Nv))
		Concinter = Conc(tab_fe(1,Nv))
		Cmob = Conc(tab_fe(-mv,Nv):tab_fe(mi,Nv))
		do nloop = -Nv+1,-mv-1
			ConcP(tab_fe(nloop,Nv)) = 0._dp
			if (nloop .eq. -8) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_8v
			end if
			do mloop = -mv, mi
				if (((nloop+mloop) > -Nv) .and. ((nloop-mloop) > -Nv)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
											& (Beta_tab(nloop,mloop)*Cmob(mloop) + Alpha_tab(nloop,mloop))*Conc(tab_fe(nloop,Nv)) + &
											& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv))
				end if
			end do			
		end do
		do nloop = mi+1, Ni-1
			ConcP(tab_fe(nloop,Nv)) = 0._dp
			if (nloop .eq. 4) then
				ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_4i
			endif
			do mloop = -mv, mi
				if (((nloop+mloop) < Ni) .and. ((nloop-mloop) < Ni)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
											& (Beta_tab(nloop,mloop)*Cmob(mloop) + Alpha_tab(nloop,mloop))*Conc(tab_fe(nloop,Nv)) + &
											& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv))
				end if
			end do				
		end do
		if (quasi) then
			! contribution des amas immobiles
			do nloop = -mv,mi
				do mloop = mi+1, Nf_Inter-mi
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
				do mloop = -Nf_Vac+mv, -mv-1 
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution amas mobiles 
			do nloop = -mv,-1
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop) - &
									& Beta_tab(mloop,nloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) + &
									& Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(nloop+mloop,Nv))
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			do nloop = 1, mi
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop) - &
									& Beta_tab(mloop,nloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) + &
									& Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(nloop+mloop,Nv)) 
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution totale (amas stochastiques, forces de puits et creation de defauts)
			do nloop = -mv, mi
				bCnCmob(nloop) = bCnCmob(nloop) + bCn_sto(nloop)*Cmob(nloop)
				aCmob(nloop) = aCmob(nloop) + aCmob_sto(nloop) 
			end do
			Si = Si + Si_sto
			Sv = Sv + Sv_sto			
			! contribution totale
			do nloop = -mv, mi
				ConcP(tab_fe(nloop,Nv)) = bCnCmob(nloop) + aCmob(nloop)
				if (nloop.eq.(1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_i - Z_i*rho_d*D_i*(Cmob(nloop)-Ci_eq) - &
											&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Cmob(nloop)-Ci_eq)
				else if (nloop.eq.(-1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_v - Z_v*rho_d*D_v*(Cmob(nloop)-Cv_eq) - &
											&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Cmob(nloop)-Cv_eq)
				end if
			end do
		else
			! contribution des amas immobiles
			do nloop = -mv,mi
				do mloop = mi+1, Ni-mi
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
				do mloop = -Nv+mv, -mv-1 
					bCnCmob(nloop) = bCnCmob(nloop) - Beta_tab(mloop,nloop)*Conc(tab_fe(mloop,Nv))*Cmob(nloop)
					aCmob(nloop) = aCmob(nloop) + Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(mloop+nloop,Nv))
					if (nloop.eq.1) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution amas mobiles 
			do nloop = -mv,-1
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop) - &
									& Beta_tab(mloop,nloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) + &
									& Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(nloop+mloop,Nv))
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			do nloop = 1, mi
				do mloop = -mv, mi
					bCnCmob(nloop) = bCnCmob(nloop) + Beta_tab(nloop-mloop,mloop)*Conc(tab_fe(nloop-mloop,Nv))*Cmob(mloop) - &
									& Beta_tab(nloop,mloop)*Cmob(nloop)*Cmob(mloop) - &
									& Beta_tab(mloop,nloop)*Cmob(nloop)*Cmob(mloop)
					aCmob(nloop) = aCmob(nloop) - Alpha_tab(nloop,mloop)*Cmob(nloop) + &
									& Alpha_tab(nloop+mloop,mloop)*Conc(tab_fe(nloop+mloop,Nv)) + &
									& Alpha_tab(mloop+nloop,nloop)*Conc(tab_fe(nloop+mloop,Nv)) 
					if (nloop.eq.(1)) then
						Si = Si + Beta_tab(mloop,1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(-1,mloop)*Conc(tab_fe(mloop,Nv))
					else if (nloop.eq.(-1)) then
						Sv = Sv + Beta_tab(mloop,-1)*Conc(tab_fe(mloop,Nv)) + Beta_tab(1,mloop)*Conc(tab_fe(mloop,Nv))
					end if
				end do
			end do
			! contribution totale
			do nloop = -mv, mi
				ConcP(tab_fe(nloop,Nv)) = bCnCmob(nloop) + aCmob(nloop)
				if (nloop.eq.(1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_i - Z_i*rho_d*D_i*(Cmob(nloop)-Ci_eq) - &
											&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Cmob(nloop)-Ci_eq)
				else if (nloop.eq.(-1)) then
					ConcP(tab_fe(nloop,Nv)) = ConcP(tab_fe(nloop,Nv)) + G_v - Z_v*rho_d*D_v*(Cmob(nloop)-Cv_eq) - &
											&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Cmob(nloop)-Cv_eq)
				end if
			end do	
		end if
	else
		print *, "Cas physique inexistant"
		IER = 1
	end if
end subroutine























end module fonctioncvode



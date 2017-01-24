module fonctioncvode


use prec_mod
use parametres
use tools


implicit none
contains
!!!
!
! A REVOIR : la distinction quasi/pas quasi (pas sur que necessaire !!)
!
!! ---- fonction F pour solver ---- !!


subroutine FCVFUN(T, Conc, ConcP, IPAR, RPAR, IER) bind(c,name='fcvfun_')
	implicit none
	real(dp) :: T, Conc1, Concvac, Concinter
	real(dp) :: Sv, Si
	integer(c_int) :: IER
	integer :: nloop, mloop, sloop
	real(dp), dimension(1:Neq) :: Conc, ConcP
	real(dp), dimension(:) :: RPAR
	integer(c_long), dimension(:) :: IPAR
	real(dp), dimension(-mv:mi) :: Cmob, bCnCmob, aCmob
	real(dp), dimension(-mv:mi,0:ms) :: CmobSol, bCnCmobSol, aCmobSol
	integer :: index_mob
	type(cluster) :: C_nu, C_mu, C_mob, C_diff
	type(cluster) :: MonoVac, MonoInter
	MonoVac = cluster(-1._dp, 0._dp,.True.,0)
	MonoVac%ind = C2I(MonoVac)
	MonoInter = cluster(1._dp, 0._dp,.True.,0)
	MonoInter%ind = C2I(MonoInter)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
			Si = Si - 2._dp*Beta_tab(1,1)*Cmob(1)
			Sv = Sv - 2._dp*Beta_tab(-1,-1)*Cmob(-1)
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
		print *, T
	else if (cas_physique.eq.3) then
		!! definition aCmob, bCnCmob
		Si = 0._dp
		Sv = 0._dp
		!bCnCmob = 0._dp
		!aCmob = 0._dp
		ConcP(1:Neq) = 0._dp
		IER = 0
		! Dynamique amas immobiles
		do nloop = 1, Neq
			C_nu = Det(nloop)
			if (.not.(C_nu%mobile)) then
				do mloop = 1, Nmob
					C_mob = Mob(mloop)
					! Flux nu-mu to nu
					C_diff = C_nu-C_mob
					if (IsInDomain(C_diff)) then
						ConcP(nint(C_nu%ind)) = ConcP(nint(C_nu%ind)) + &
										& Beta_clust(C_diff,C_mob)*Conc(nint(C_diff%ind))*Conc(nint(C_mob%ind)) - &
										& Alpha_clust(C_nu,C_mob)*Conc(nint(C_nu%ind))
					end if
					! Flux nu to nu+mu
					C_diff = C_nu+C_mob
					if (IsInDomain(C_diff)) then
						ConcP(nint(C_nu%ind)) = ConcP(nint(C_nu%ind)) - &
										& Beta_clust(C_nu,C_mob)*Conc(nint(C_nu%ind))*Conc(nint(C_mob%ind)) + &
										& Alpha_clust(C_diff,C_mob)*Conc(nint(C_diff%ind))
					end if
				end do
			end if
		end do
		! Dynamique amas mobiles
		do nloop = 1, Nmob
			C_nu = Mob(nloop)
			ConcP(nint(C_nu%ind)) = 0._dp
			do mloop = 1, Nmob
				C_mob = Mob(mloop)
				! Flux nu-mu to nu sur mobiles M_nu
				C_diff = C_nu-C_mob
				if (IsInDomain(C_diff)) then
					ConcP(nint(C_nu%ind)) = ConcP(nint(C_nu%ind)) + &
									& Beta_clust(C_diff,C_mob)*Conc(nint(C_diff%ind))*Conc(nint(C_mob%ind)) - &
									& Alpha_clust(C_nu,C_mob)*Conc(nint(C_nu%ind)) 
				end if	
				! Flux nu to nu+mu sur mobiles M
				C_diff = C_nu+C_mob	
				if (IsInDomain(C_diff)) then		 
					ConcP(nint(C_nu%ind)) = ConcP(nint(C_nu%ind)) - &
									 & Beta_clust(C_nu,C_mob)*Conc(nint(C_nu%ind))*Conc(nint(C_mob%ind)) + &
									 & Alpha_clust(C_diff,C_mob)*Conc(nint(C_diff%ind))				
				end if
				! Forces de puits
				if (C_nu%defect.eq.1) then
					Si = Si + Beta_Clust(MonoVac,C_mob)*Conc(nint(C_mob%ind))
				else if (C_nu%defect.eq.-1) then
					Sv = Sv + Beta_Clust(MonoInter,C_mob)*Conc(nint(C_mob%ind))
				end if
			end do
			do mloop = 1, Neq
				! Flux mu to nu+mu sur total Omega
				C_mu = Det(mloop)
				C_diff = C_nu+C_mu	
				if (IsInDomain(C_diff)) then
					ConcP(nint(C_nu%ind)) = ConcP(nint(C_nu%ind)) - &
									& Beta_clust(C_mu,C_nu)*Conc(nint(C_nu%ind))*Conc(nint(C_mu%ind)) + &
									& Alpha_clust(C_diff,C_nu)*Conc(nint(C_diff%ind))
				end if
				! Forces de puits
				if (C_nu%defect.eq.1) then
					Si = Si + Beta_Clust(C_mu,MonoInter)*Conc(nint(C_mu%ind))
				else if (C_nu%defect.eq.-1) then
					Sv = Sv + Beta_Clust(C_mu,MonoVac)*Conc(nint(C_mu%ind))
				end if
			end do
		end do
		! Gestion des termes sources
		do nloop = 1, Neq
			ConcP(nloop) = ConcP(nloop) + G_source(nloop)
		end do
		! Ajout des forces de puits
		Si = Si - 2._dp*Beta_Clust(MonoInter,MonoInter)*Conc(nint(MonoInter%ind))
		Sv = Sv - 2._dp*Beta_Clust(MonoVac,MonoVac)*Conc(nint(MonoVac%ind))
		ConcP(nint(MonoInter%ind)) = ConcP(nint(MonoInter%ind)) - Z_i*rho_d*D_i*(Conc(nint(MonoInter%ind)) - Ci_eq) - &
											&6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Conc(nint(MonoInter%ind)) - Ci_eq)
		ConcP(nint(MonoVac%ind)) = ConcP(nint(MonoVac%ind)) - Z_v*rho_d*D_v*(Conc(nint(MonoVac%ind)) - Cv_eq) - &
											&6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Conc(nint(MonoVac%ind)) - Cv_eq)
		! Gestion de (0,0)
		ConcP(nint(C_null%ind)) = 0._dp
		print *, T
	else
		print *, "Cas physique inexistant"
		IER = 1
	end if
end subroutine























end module fonctioncvode



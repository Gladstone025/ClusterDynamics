module stochastique

use prec_mod
use parametres
use tools
!use data_fun
!use mpi

implicit none

contains 

function Calcul_Cvac(Conc,HLangevin)
	implicit none
	integer :: Nconc, Npart
	real(dp) :: Calcul_Cvac
	real(dp) :: bCn, aCn
	real(dp), dimension(:) :: Conc
	real(dp), dimension(:) :: HLangevin
	bCn = 2.*beta(1.d0)*Conc(1)
	aCn = alpha(2.d0)*Conc(2)
	Nconc = size(Conc)
	Npart = size(HLangevin)
	do iloop = 2, Nconc-1 
		rloop = float(iloop)
		aCn = aCn + alpha(rloop)*Conc(iloop)
		bCn = bCn + beta(rloop)*Conc(iloop)
	enddo
	do iloop = 1, Npart
		bCn = bCn + MStochastique/Npart*beta(HLangevin(iloop))
		aCn = aCn + MStochastique/Npart*alpha(HLangevin(iloop))
	enddo
	Calcul_Cvac = (G1+aCn)/bCn
end function 


function Calcul_Cvac_strang(Conc,HLangevin,dts,Ts)
	implicit none
	integer :: Nconc, Npart, Ks
	real(dp) :: Calcul_Cvac_strang, Cvac1, Cvac2
	real(dp) :: dts, dts2, Ts
	real(dp) :: bCn, aCn
	real(dp), dimension(:) :: Conc
	real(dp), dimension(:) :: HLangevin
	Nconc = size(Conc)
	Npart = size(HLangevin)
	Ks = int(Ts/dts)
	dts2 = dts/2.
	Cvac1 = Conc(1)
	Cvac2 = Conc(1)
	bCn = 0 !beta(1.d0)*C(1)
	aCn = alpha(2.d0)*Conc(2)
	do iloop = 2, Nconc-1 
		rloop = float(iloop)
		aCn = aCn + alpha(rloop)*Conc(iloop)
		bCn = bCn + beta(rloop)*Conc(iloop)
	enddo
	do iloop = 1, Npart
		bCn = bCn + MStochastique/Npart*beta(HLangevin(iloop))
		aCn = aCn + MStochastique/Npart*alpha(HLangevin(iloop))
	enddo
	do kloop = 1, Ks
		Cvac1 = (Cvac2 - aCn/bCn)*exp(-bCn*dts2) + aCn/bCn
		Cvac2 = Cvac1/(1.+2.*beta(1.d0)*dts*Cvac1)
		Cvac1 = (Cvac2 - aCn/bCn)*exp(-bCn*dts2) + aCn/bCn
	enddo
	Calcul_Cvac_strang = Cvac1
end function


function Calcul_ConcMob(Conc,HInter,HVac,dts,Ts)
	implicit none
	integer :: Nconc, NpartInter, NpartVac, Ks
	real(dp) :: Cmob1, Cmob2
	real(dp) :: dts, dts2, Ts
	real(dp) :: bCmob, aCmob, Si, Sv
	real(dp), dimension(1:mv+mi+1) :: Calcul_ConcMob
	real(dp), dimension(:) :: Conc
	real(dp), dimension(:) :: HInter, Hvac
	Nconc = size(Conc)
	NpartInter = size(HInter)
	NpartVac = size(Hvac)
	Ks = int(Ts/dts)
	dts2 = dts/2.
	! A enrichir avec forces de puits et les quelques termes qui manquent dus au amas mobiles
	do iloop = -mv, mi
		rloop = real(iloop,8)
		Cmob1 = Conc(tab_fe(iloop,Nv))
		Cmob2 = Conc(tab_fe(iloop,Nv))
		aCmob = 0._dp
		bCmob = 0._dp
		Si = 0._dp
		Sv = 0._dp
		if (iloop.eq.-mv) then 
			aCmob = G_v - Alpha_tab(-1,1)*Conc(tab_fe(-1,Nv)) + 2._dp*Alpha_tab(-2,-1)*Conc(tab_fe(-2,Nv))
		else if (iloop.eq.mi) then
			aCmob = G_i + 2._dp*Alpha_tab(2,1)*Conc(tab_fe(2,Nv)) - Alpha_tab(1,-1)*Conc(tab_fe(1,Nv)) 
		end if
		do jloop = max(-Nv-iloop,-Nv), -2
			aCmob = aCmob + Alpha_tab(jloop+iloop,iloop)*Conc(tab_fe(jloop+iloop,Nv))
			bCmob = bCmob + Beta_tab(jloop,iloop)*Conc(tab_fe(jloop,Nv))
			Si = Si + Beta_tab(jloop,1)*Conc(tab_fe(jloop,Nv))
			Sv = Sv + Beta_tab(jloop,-1)*Conc(tab_fe(jloop,Nv))
		enddo
		do jloop = 2, min(Ni-iloop,Ni) 
			aCmob = aCmob + Alpha_tab(jloop+iloop,iloop)*Conc(tab_fe(jloop+iloop,Nv))
			bCmob = bCmob + Beta_tab(jloop,iloop)*Conc(tab_fe(jloop,Nv))
			Si = Si + Beta_tab(jloop,1)*Conc(tab_fe(jloop,Nv))
			Sv = Sv + Beta_tab(jloop,-1)*Conc(tab_fe(jloop,Nv))
		enddo
		if (Coupling_Inter) then
			do jloop = 1, NpartInter
				aCmob = aCmob + MStoInter/NpartInter*alpha_nm(HInter(jloop),rloop)
				bCmob = bCmob + MStoInter/NpartInter*beta_nm(HInter(jloop),rloop)
				Si = Si + MStoInter/NpartInter*beta_nm(HInter(jloop),1._dp)
				Sv = Sv + MStoInter/NpartInter*beta_nm(HInter(jloop),-1._dp)
			enddo
		end if
		if (Coupling_Vac) then
			do jloop = 1, NpartVac
				aCmob = aCmob + MStoVac/NpartVac*alpha_nm(HVac(jloop),rloop)
				bCmob = bCmob + MStoVac/NpartVac*beta_nm(HVac(jloop),rloop)	
				Si = Si + MStoInter/NpartInter*beta_nm(HInter(jloop),1._dp)
				Sv = Sv + MStoInter/NpartInter*beta_nm(HInter(jloop),-1._dp)		
			enddo
		end if
		if (iloop.eq.-mv) then 
			aCmob = aCmob - Z_v*rho_d*D_v*(Cmob1-Cv_eq) - 6._dp/l_gb*sqrt(Z_v*rho_d+Sv/D_v)*D_v*(Cmob1-Cv_eq)
		else if (iloop.eq.mi) then
			aCmob = aCmob - Z_i*rho_d*D_i*(Cmob1-Ci_eq) - 6._dp/l_gb*sqrt(Z_i*rho_d+Si/D_i)*D_i*(Cmob1-Ci_eq)
		end if
		do kloop = 1, Ks
			Cmob1 = (Cmob2 - aCmob/bCmob)*exp(-bCmob*dts2) + aCmob/bCmob
			Cmob2 = Cmob1/(1.+2.*Beta_tab(iloop,iloop)*dts*Cmob1)
			Cmob1 = (Cmob2 - aCmob/bCmob)*exp(-bCmob*dts2) + aCmob/bCmob
		enddo
		if (iloop.eq.0) then
			Cmob1 = 0._dp
		end if
		Calcul_ConcMob(tab_fe(iloop,mv)) = Cmob1
	end do
end function

function Langevin1part(x,ConcMob,Tfinal)
	implicit none
	integer :: K
	real(dp) :: x, Langevin1part, Tfinal
	real(dp) :: C1, sqrtdt
	real(dp), dimension(1) :: G
	real(dp), dimension(:) :: ConcMob
	K = int(Tfinal/dt_sto)
	sqrtdt = sqrt(dt_sto)
	Langevin1part = x
	if (cas_physique.eq.1) then
		C1 = ConcMob(1)
		do iloop = 1, K
			call gen_rand_normal(G,1)
			Langevin1part = Langevin1part + dt_sto*(F_scal(Langevin1part,C1)) + sqrtdt*sqrt(D_scal(Langevin1part,C1))*G(1)
		end do
	else if (cas_physique.eq.2) then
		do iloop = 1, K
			call gen_rand_normal(G,1)
			Langevin1part = Langevin1part + dt_sto*(F_fe(Langevin1part,ConcMob)) + &
								&sqrtdt*sqrt(D_fe(Langevin1part,ConcMob))*G(1)
		end do
	end if
end function Langevin1part


function EKMC1part(x,ConcMob,Tfinal)
	implicit none
	real(dp), dimension(:) :: ConcMob
	real(dp) :: x, EKMC1part
	real(dp) :: Tfinal, C1
	real(dp) :: dte, nue, expe, pos, u, s
	real(dp), dimension(:), allocatable :: p
	if (cas_physique.eq.1) then
		allocate(p(1))
		C1 = ConcMob(1)
		EKMC1part = x
		dte = 0.d0
		do while (dte < Tfinal)
			pos = EKMC1part
			nue = beta(pos)*C1+alpha(pos)
			call tirage_exp(expe,nue)
			dte = dte + expe
			if (dte < Tfinal) then
				p(1) = beta(pos)*C1/(beta(pos)*C1+alpha(pos))
				call random_number(u)
				if (u.le.p(1)) then 
					EKMC1part = EKMC1part+1.
				else 
					EKMC1part = EKMC1part-1.
				endif
			endif
		enddo
		deallocate(p)
	else if (cas_physique.eq.2) then
		allocate(p(-mv:mi))
		EKMC1part = x
		dte = 0.d0
		do while (dte < Tfinal)
			pos = EKMC1part
			nue = 0._dp
			do jloop = -mv, mi
				rloop = real(jloop,8)
				nue = nue + beta_nm(pos,rloop)*ConcMob(tab_fe(jloop,mv))+alpha_nm(pos,-rloop)
			end do
			call tirage_exp(expe,nue)
			dte = dte + expe
			if (dte < Tfinal) then
				! Pour les probas, on fait gaffe : beta_(n,i/v)C_i/v + alpha_(n,v/i) 
				! Penser a computer Alpha_tab et Beta_tab sur -max(mv,mi),max(mv,mi) avec zero pour les inexistants
				do jloop = -mv, mi
					rloop = real(jloop,8)
					p(jloop) = (beta_nm(pos,rloop)*ConcMob(tab_fe(jloop,mv))+alpha_nm(pos,-rloop))/nue
				end do
				!print *, "p + p = ", p(-mv)!+p(0)+p(mi) 
				call random_number(u)
				s = 0
				do jloop = -mv, mi
					s = s + p(jloop)
					if (s > u) then 
						EKMC1part = EKMC1part + real(jloop,8)
						exit
					endif
				end do
			endif
		enddo
		deallocate(p)
	end if
end function EKMC1part

subroutine PropagationSto(HLangevin,ConcMob,Tfinal)
	implicit none
	real(dp) :: Tfinal
	real(dp), dimension(:) :: ConcMob
	real(dp), dimension(:) :: HLangevin
	do jloop = 1, size(HLangevin)
		if ((HLangevin(jloop) > Nf_Inter+Nf_EL_Inter) .or. (HLangevin(jloop) < -(Nf_Vac+Nf_EL_Vac))) then
			HLangevin(jloop) = Langevin1part(HLangevin(jloop),ConcMob,Tfinal)
		else
			HLangevin(jloop) = EKMC1part(real(nint(HLangevin(jloop)),8),ConcMob,Tfinal)
		end if
	end do
end subroutine PropagationSto


! revoir toutes les boucles
subroutine Langevin(HLangevin,ConcMob,Tfinal)
	implicit none
	integer :: K
	real(dp) :: Tfinal
	real(dp) :: C1, sqrtdt
	real(dp), dimension(1) :: G
	real(dp), dimension(:) :: ConcMob
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(size(HLangevin)) :: G_vec
	K = int(Tfinal/dt_sto)
	sqrtdt = sqrt(dt_sto)
	if (cas_physique.eq.1) then
		C1 = ConcMob(1)
		do iloop = 1, K
			call gen_rand_normal(G_vec,Taille)
			HLangevin(:) = HLangevin(:) + dt_sto*(F_vec(HLangevin(:),C1)) + sqrtdt*sqrt(D_vec(HLangevin(:),C1))*G_vec
		end do
	else if (cas_physique.eq.2) then
		do iloop = 1, K
			print *, iloop
			do jloop = 1, size(HLangevin)
				call gen_rand_normal(G,1)
				HLangevin(jloop) = HLangevin(jloop) + dt_sto*(F_fe(HLangevin(jloop),ConcMob)) + &
									&sqrtdt*sqrt(D_fe(HLangevin(jloop),ConcMob))*G(1)
			end do
		end do
	end if
end subroutine Langevin


function Metropolis(Conc,N_0,alpha_m,Npart)
	implicit none
	integer :: Npart, Ntot
	real(dp) :: N_0, n0, n1, uniform, alea, U, p, acc, q0, q1, alpha_m
	real(dp), dimension(:) :: Conc
	real(dp), dimension(1:Npart) :: Metropolis
	n0 = N_0
	Ntot = size(Conc)
	print *, size(Conc)
	do iloop = 1, Npart
		call random_number(uniform)
		alea = -alpha_m + 2.*alpha_m*uniform
		n1 = min(real(Ntot-2,8),n0 + alea)
		call random_number(U)
		q1 = approx_lineaire(n1,Conc(int(n1)),Conc(int(n1)+1),real(int(n1),8),real(int(n1)+1,8))
		q0 = approx_lineaire(n0,Conc(int(n0)),Conc(int(n0)+1),real(int(n0),8),real(int(n0)+1,8))
		p = min(1.d0,q1/q0)
		if (U <= p .and. n1 >= N_front) then
			n0 = n1
			acc = acc + 1.d0
		endif
		Metropolis(iloop) = n0
	end do
	!print *, "taux d'acceptation : ", acc/float(Npart)
end function

subroutine SSA(HLangevin,ConcMob,Tfinal)
	implicit none
	integer :: Npart
	real(dp), dimension(:) :: ConcMob
	real(dp), dimension(:) :: HLangevin
	real(dp) :: Tfinal, C1
	real(dp) :: dte, nue, expe, pos, u, s
	real(dp), dimension(:), allocatable :: p
	Npart = size(HLangevin)
	if (cas_physique.eq.1) then
		allocate(p(1))
		C1 = ConcMob(1)
		do iloop = 1, Npart
			dte = 0.d0
			do while (dte < Tfinal)
				pos = HLangevin(iloop)
				nue = beta(pos)*C1+alpha(pos)
				call tirage_exp(expe,nue)
				dte = dte + expe
				if (dte < Tfinal) then
					p(1) = beta(pos)*C1/(beta(pos)*C1+alpha(pos))
					call random_number(u)
					if (u.le.p(1)) then 
						HLangevin(iloop) = HLangevin(iloop)+1.
					else 
						HLangevin(iloop) = HLangevin(iloop)-1.
					endif
				endif
			enddo
		enddo
		deallocate(p)
	else if (cas_physique.eq.2) then
		allocate(p(-mv:mi))
		do iloop = 1, Npart
			dte = 0.d0
			do while (dte < Tfinal)
				pos = HLangevin(iloop)
				nue = 0._dp
				do jloop = -mv, mi
					rloop = real(jloop,8)
					nue = nue + beta_nm(pos,rloop)*ConcMob(tab_fe(jloop,mv))+alpha_nm(pos,-rloop)
				end do
				call tirage_exp(expe,nue)
				dte = dte + expe
				if (dte < Tfinal) then
					! Pour les probas, on fait gaffe : beta_(n,i/v)C_i/v + alpha_(n,v/i) 
					! Penser a computer Alpha_tab et Beta_tab sur -max(mv,mi),max(mv,mi) avec zero pour les inexistants
					do jloop = -mv, mi
						rloop = real(jloop,8)
						p(jloop) = (beta_nm(pos,rloop)*ConcMob(tab_fe(jloop,mv))+alpha_nm(pos,-rloop))/nue
					end do
					!print *, "p + p = ", p(-mv)!+p(0)+p(mi) 
					call random_number(u)
					s = 0
					do jloop = -mv, mi
						s = s + p(jloop)
						if (s > u) then 
							HLangevin(iloop) = HLangevin(iloop) + real(jloop,8)
							exit
						endif
					end do
				endif
			enddo
		enddo
		deallocate(p)
	end if
end subroutine

function Multinomial(Conc,Nmin,Nmax,Npart)
	implicit none
	integer :: Npart, Ntot, Nmin, Nmax
	real(dp) :: s, u
	real(dp), dimension(:) :: Conc
	real(dp), dimension(size(Conc)) :: P
	real(dp), dimension(Npart) :: Multinomial
	Ntot = size(Conc)
	P = Conc/sum(Conc)
	do iloop = 1, Npart
		call random_number(u)
		s = 0
		do kloop = Nmin, Nmax
			s = s + P(kloop)
			if (s > u) then
				Multinomial(iloop) = real(kloop,8)
				exit
			endif
		end do	
	end do
end function


function Inboundary(x,InterVac,a)
	implicit none
	real(dp) :: x, a ! a sert a prendre en compte la masse qui depasse un peu
					 ! du fait de la methode des noyaux
	integer :: InterVac
	logical :: Inboundary
	Inboundary = .False.
	
	select case (InterVac)
	! - Cas simple
	case(0)
		if (x < N_front+a) then
			Inboundary = .True.
		end if
	! - Cas lacunes
	case(-1)
		if (x > -Nf_Vac-a) then
			Inboundary = .True.
		end if
	! - Cas Interstitiel
	case(1)
		if (x < Nf_Inter+a) then
			Inboundary = .True.
		end if	
	end select
end function

! a reecrire
function StotoDiscret(HLangevin,InterVac)
	implicit none
	integer :: InterVac
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(:), allocatable :: Conc, StotoDiscret
	integer :: Npart, n, ninf
	ninf = 0
	Npart = size(HLangevin)
	! Penser a ecrire une fonction du type "Inboundary" qui renvoit un booleen
	select case (InterVac)
	case(0)
		allocate(Conc(1:N_front+N_buff))
		allocate(StotoDiscret(1:N_front+N_buff))
	case(1)
		allocate(Conc(1:Nf_Inter+Nb_Inter))
		allocate(StotoDiscret(1:Nf_Inter+Nb_Inter))
	case(-1)
		allocate(Conc(-Nf_Vac-Nb_Vac:-1))
		allocate(StotoDiscret(1:Nf_Vac+Nb_Vac))
	end select
	Conc = 0._dp
	do iloop = 1, Npart
		if (Inboundary(HLangevin(iloop),InterVac,2._dp)) then
			n = nint(HLangevin(iloop))
			if (Inboundary(HLangevin(iloop),InterVac,0._dp)) then
				ninf = ninf + 1
			end if
			if (methode.eq.1) then
				Conc(n) = Conc(n) + 1._dp	
			else if (methode.eq.2) then
				do kloop = -2, 2
					Conc(n+kloop) = Conc(n+kloop) + 1./h*Noyau(1./h*(n+kloop-HLangevin(iloop)))
				end do
			else if (methode.eq.3) then
				if ((n > Nf_Inter+Nf_EL_Inter) .or. (n < -(Nf_Vac+Nf_EL_Vac))) then
					do kloop = -2, 2
						Conc(n+kloop) = Conc(n+kloop) + 1./h*Noyau(1./h*(n+kloop-HLangevin(iloop)))
					end do
				else
					Conc(n) = Conc(n) + 1._dp
				end if
			end if
		end if
	end do
	if (InterVac.eq.0) then
		MDiscret = MStochastique*(real(ninf,8)/real(Npart,8))
		Conc(N_front:N_front+N_buff) = 0._dp
		if (sum(Conc) > 0._dp) then
			Conc = MDiscret*Conc/sum(Conc)
		end if
		MStochastique = MStochastique - MDiscret
		StotoDiscret = Conc
	else if (InterVac.eq.1) then
		MDisInter = MStoInter*(real(ninf,8)/real(Npart,8))
		Conc(Nf_Inter:Nf_Inter+Nb_Inter) = 0._dp
		if (sum(Conc) > 0._dp) then
			Conc = MDisInter*Conc/sum(Conc)
		end if
		MStoInter = MStoInter - MDisInter
		StotoDiscret = Conc
	else if (InterVac.eq.(-1)) then
		MDisVac = MStoVac*(real(ninf,8)/real(Npart,8))
		Conc(-Nf_Vac-Nb_Vac:-Nf_Vac) = 0._dp
		if (sum(Conc) > 0._dp) then
			Conc = MDisVac*Conc/sum(Conc)
		end if
		MStoVac = MStoVac - MDisVac
		StotoDiscret(1:Nf_Vac+Nb_Vac) = Conc(-1:-Nf_Vac-Nb_Vac:-1)
	end if
end function



function Sampling(Conc,HLangevin,InterVac)
	implicit none
	integer :: InterVac
	real(dp), dimension(:) :: Conc
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(:), allocatable :: Conc_buff
	real(dp), dimension(size(HLangevin)) :: Sampling
	real(dp), dimension(:), allocatable :: Htot, Hsto
	real(dp), dimension(:), allocatable :: Sampling_inter
	integer :: Nconc, Npart, Ninf, Ndiff
	integer :: compt, nn, nsto
	real(dp) :: Mconc, Mlangevin
	real(dp) :: a, u
	! -- On initialise
	Ninf = 0
	compt = 0
	Npart = size(HLangevin)
	a = 0.5
	if (InterVac.eq.0) then
		allocate(Conc_buff(1:Neq))
		Conc_buff = 0._dp
		Conc_buff(N_front:N_front+N_buff) = Conc(N_front:N_front+N_buff)
		Npart = size(HLangevin)
		Mconc = sum(Conc(N_front:N_front+N_buff))
		Nconc = int(Mconc*Npart/MStochastique)
		MStochastique = MStochastique + Mconc
	else if (InterVac.eq.1) then
		allocate(Conc_buff(1:Ni))
		Conc_buff = 0._dp
		Conc_buff(Nf_Inter:Nf_Inter+Nb_Inter) = Conc(tab_fe(Nf_Inter,Nv):tab_fe(Nf_Inter+Nb_Inter,Nv))
		Npart = size(HLangevin)
		print *, Npart
		Mconc = sum(Conc_buff(Nf_Inter:Nf_Inter+Nb_Inter))
		print *, Mconc, MStoInter
		Nconc = int(Mconc*Npart/MStoInter)
		print *, Nconc
		MStoInter = MStoInter + Mconc	
	else if (InterVac.eq.(-1)) then
		allocate(Conc_buff(1:Nv))
		Conc_buff = 0._dp
		Conc_buff(Nf_Vac:Nf_Vac+Nb_Vac) = Conc(tab_fe(-Nf_Vac,Nv):tab_fe(-Nb_Vac-Nf_Vac,Nv):-1)
		Npart = size(HLangevin)
		print *, Npart
		Mconc = sum(Conc_buff(Nf_Vac:Nf_Vac+Nb_Vac))
		print *, Mconc, MStoVac
		Nconc = int(Mconc*Npart/MStoVac)
		print *, Nconc
		MStoVac = MStoVac + Mconc		
	end if
	!if (Nconc > 2*Npart) then
	!	DeltaT  = DeltaT/4._dp
	!	dt_sto = dt_sto/4._dp
	!end if
	allocate(Hsto(Nconc))
	if (methode.eq.1) then
		select case (InterVac)
		case(0)
			Hsto = Multinomial(Conc_buff,N_front,N_front+N_buff,Nconc)
		case(1)
			Hsto = Multinomial(Conc_buff,Nf_Inter,Nf_Inter+Nb_Inter,Nconc)
		case(-1)
			Hsto = -1._dp*Multinomial(Conc_buff,Nf_Vac-Nb_Vac,Nf_Vac,Nconc)
		end select
	else 
		select case (InterVac)
		case(0)
			Hsto = Metropolis(Conc_buff,real(N_front+N_buff/5._dp,8),real(N_buff/10._dp,8),Nconc)
		case(1)
			Hsto = Metropolis(Conc_buff,real(Nf_Inter+Nb_Inter/5._dp,8),real(Nb_Inter/10._dp,8),Nconc)
		case(-1)
			Hsto = -1._dp*Metropolis(Conc_buff,real(Nf_Vac+Nb_Vac/5._dp,8),real(Nb_Vac/10._dp,8),Nconc)
		end select
	end if
	print *, Npart+Nconc
	allocate(Htot(Npart+Nconc))
	Htot(1:Npart) = HLangevin(1:Npart)
	Htot(Npart+1:Npart+Nconc) = Hsto(1:Nconc)
	! -- Htot contient *toutes* les particules
	! -- On compte celles qui contribuent a la partie deterministe
	do iloop = 1, Npart+Nconc
		if (Inboundary(Htot(iloop),InterVac,0._dp)) then
			Ninf = Ninf+1
		endif
	end do
	Ndiff = Nconc - Ninf
	print *, "Ninf : ", Ninf, " / Nconc : ", Nconc
	print *, Inboundary(Htot(1),InterVac,0._dp), Inboundary(Htot(Npart+Nconc-10),InterVac,0._dp), Htot(1), Htot(Npart+Nconc-10)
	print *, InterVac, " / Ndiff : ", Ndiff
	! -- Puis on gere le resampling en fonction de Ndiff
	if (Ndiff < 0) then ! -- il faudra dupliquer des particules
	print *, "DUPLICATION"
		do iloop = 1, Npart+Nconc
			if (.not.Inboundary(Htot(iloop),InterVac,0._dp)) then
				compt = compt+1
				Sampling(compt) = Htot(iloop)
				nsto = nint(Htot(iloop))
				if (methode.eq.1) then
					C_sto(nsto) = C_sto(nsto) + 1._dp	
				else if (methode.eq.2) then
					do kloop = -2, 2
						C_sto(nsto+kloop) = C_sto(nsto+kloop) + 1./h*Noyau(1./h*(nsto+kloop-Htot(iloop)))
					end do
				else if (methode.eq.3) then
					if ((nsto > Nf_Inter+Nf_EL_Inter) .or. (nsto < -(Nf_Vac+Nf_EL_Vac))) then
						do kloop = -2, 2
							Conc(nsto+kloop) = Conc(nsto+kloop) + 1./h*Noyau(1./h*(nsto+kloop-HLangevin(iloop)))
						end do
					else
						Conc(nsto) = Conc(nsto) + 1._dp
					end if
				end if
			endif
		end do
		do iloop = compt+1,Npart
			call random_number(u)
			nn = int(u*compt)+1
			Sampling(iloop) = Sampling(nn)
			nsto = nint(Sampling(iloop))
			if (methode.eq.1) then
				C_sto(nsto) = C_sto(nsto) + 1._dp	
			else if (methode.eq.2) then
				do kloop = -2, 2
					C_sto(nsto+kloop) = C_sto(nsto+kloop) + 1./h*Noyau(1./h*(nsto+kloop-Htot(iloop)))
				end do
			else if (methode.eq.3) then
				if ((nsto > Nf_Inter+Nf_EL_Inter) .or. (nsto < -(Nf_Vac+Nf_EL_Vac))) then
					do kloop = -2, 2
						Conc(nsto+kloop) = Conc(nsto+kloop) + 1./h*Noyau(1./h*(nsto+kloop-HLangevin(iloop)))
					end do
				else
					Conc(nsto) = Conc(nsto) + 1._dp
				end if
			end if
		enddo
	else ! -- il faut supprimer des particules
	print *, "SUPPRESSION"
		allocate(Sampling_inter(Npart+Ndiff))
		Sampling_inter = 0._dp
		compt = 0
		do iloop = 1, Npart+Nconc
			if (.not.Inboundary(Htot(iloop),InterVac,0._dp)) then
				compt = compt+1
				Sampling_inter(compt) = Htot(iloop)
			endif
		end do	
		compt = 0
		do while (compt < Ndiff)
			call random_number(u)
			nn = int(u*(Npart+Ndiff))+1
			if (Sampling_inter(nn) > R_min) then
				Sampling_inter(nn) = R_min
				compt = compt+1
			end if
		end do
		compt = 0
		do iloop = 1, Npart+Ndiff
			if (Sampling_inter(iloop) > R_min) then
				compt = compt + 1
				Sampling(compt) = Sampling_inter(iloop)
				nsto = nint(Sampling_inter(iloop))
				if (methode.eq.1) then
					C_sto(nsto) = C_sto(nsto) + 1._dp	
				else if (methode.eq.2) then
					do kloop = -2, 2
						C_sto(nsto+kloop) = C_sto(nsto+kloop) + 1./h*Noyau(1./h*(nsto+kloop-Htot(iloop)))
					end do
				else if (methode.eq.3) then
					if ((nsto > Nf_Inter+Nf_EL_Inter) .or. (nsto < -(Nf_Vac+Nf_EL_Vac))) then
						do kloop = -2, 2
							Conc(nsto+kloop) = Conc(nsto+kloop) + 1./h*Noyau(1./h*(nsto+kloop-HLangevin(iloop)))
						end do
					else
						Conc(nsto) = Conc(nsto) + 1._dp
					end if
				end if
			end if
		end do
		deallocate(Sampling_inter)
	endif
	if (InterVac.eq.0) then
		C_sto = MStochastique*C_sto/sum(C_sto)
	else if (InterVac.eq.1) then
		C_sto(1:1200) = C_sto(1:1200)/sum(C_sto(1:1200))
	else if (InterVac.eq.(-1)) then
		C_sto(-1200:-1) = C_sto(-1200:-1)/sum(C_sto(-1200:-1))			
	end if
	deallocate(Hsto)
	deallocate(Htot)
end function


function Mtail(Conc,ninf,nmax)
	implicit none
	integer :: ninf, nmax
	real(dp) :: Mtail
	real(dp), dimension(:) :: Conc
	Mtail = 0._dp
	do iloop = ninf, nmax
		Mtail = Mtail + iloop*Conc(iloop)
	end do
end function

end module stochastique



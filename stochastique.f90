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
	integer :: Nconc, Npart, Ks
	real(dp) :: Cmob1, Cmob2
	real(dp) :: dts, dts2, Ts
	real(dp) :: bCmob, aCmob
	real(dp), dimension(-mv:mi) :: Calcul_ConcMob
	real(dp), dimension(:) :: Conc
	real(dp), dimension(:) :: HInter, Hvac
	Nconc = size(Conc)
	Npart = size(HInter)
	Ks = int(Ts/dts)
	dts2 = dts/2.
	Cmob1 = Conc(1)
	Cmob2 = Conc(1)
	bCmob = 0
	aCmob = 0
	! A enrichir avec forces de puits et les quelques termes qui manquent dus au amas mobiles
	do iloop = -mv, mi
		do jloop = max(-Nv-iloop,-Nv), -2
			aCmob = aCmob + Alpha_tab(jloop+iloop,iloop)*Conc(tab_fe(jloop+iloop,Nv))
			bCmob = bCmob + Beta_tab(jloop,iloop)*Conc(tab_fe(jloop,Nv))
		enddo
		do jloop = 2, min(Ni-iloop,Ni) 
			aCmob = aCmob + Alpha_tab(jloop+iloop,iloop)*Conc(tab_fe(jloop+iloop,Nv))
			bCmob = bCmob + Beta_tab(jloop,iloop)*Conc(tab_fe(jloop,Nv))
		enddo
		do jloop = 1, Npart
			bCmob = bCmob + MStoInter/Npart*beta(HInter(jloop))
			bCmob = bCmob + MStoVac/Npart*beta(HVac(jloop))			
			aCmob = aCmob + MStoInter/Npart*alpha(HInter(jloop))
			aCmob = aCmob + MStoVac/Npart*alpha(HVac(jloop))
		enddo
		do kloop = 1, Ks
			Cmob1 = (Cmob2 - aCmob/bCmob)*exp(-bCmob*dts2) + aCmob/bCmob
			Cmob2 = Cmob1/(1.+2.*Beta_tab(iloop,iloop)*dts*Cmob1)
			Cmob1 = (Cmob2 - aCmob/bCmob)*exp(-bCmob*dts2) + aCmob/bCmob
		enddo
		Calcul_ConcMob(iloop) = Cmob1
	end do
end function



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
		allocate(p(mv+mi+1))
		do iloop = 1, Npart
			dte = 0.d0
			do while (dte < Tfinal)
				pos = HLangevin(iloop)
				nue = 0._dp
				do jloop = -mv, mi
					rloop = real(jloop,8)
					nue = nue + beta_nm(pos,rloop)*ConcMob(jloop)+alpha_nm(pos,rloop)
				end do
				call tirage_exp(expe,nue)
				dte = dte + expe
				if (dte < Tfinal) then
					! Pour les probas, on fait gaffe : beta_(n,i/v)C_i/v + alpha_(n,v/i) 
					! Penser a computer Alpha_tab et Beta_tab sur -max(mv,mi),max(mv,mi) avec zero pour les inexistants
					do jloop = -mv, mi
						p(jloop) = (beta_nm(pos,rloop)*ConcMob(jloop)+alpha_nm(pos,-rloop))/nue
					end do
					call random_number(u)
					s = 0
					do jloop = -mv, mi
						s = s + p(jloop)
						if (s > u) then 
							HLangevin(jloop) = HLangevin(jloop) + real(jloop,8)
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



function StotoDiscret(HLangevin,InterVac)
	implicit none
	integer :: InterVac
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(N_front+N_buff) :: Conc, StotoDiscret
	integer :: Npart, n, ninf
	Conc = 0._dp
	ninf = 0
	Npart = size(HLangevin)
	! Penser a ecrire une fonction du type "Inboundary" qui renvoit un booleen
	do iloop = 1, Npart
		if (HLangevin(iloop) < N_front+2.) then
			n = nint(HLangevin(iloop))
			if (HLangevin(iloop) < N_front) then
				ninf = ninf + 1
			end if
			do kloop = -5, 5
				Conc(n+kloop) = Conc(n+kloop) + 1./h*Noyau(1./h*(n+kloop-HLangevin(iloop)))
			end do
		end if
	end do
	if (InterVac.eq.0) then
		MDiscret = MStochastique*(real(ninf,8)/real(Npart,8))
		Conc(N_front:N_front+N_buff) = 0._dp
		Conc = MDiscret*Conc/sum(Conc)
		MStochastique = MStochastique - MDiscret
	else if (InterVac.eq.1) then
		MDisInter = MStoInter*(real(ninf,8)/real(Npart,8))
		Conc(N_front:N_front+N_buff) = 0._dp
		Conc = MDisInter*Conc/sum(Conc)
		MStoInter = MStoInter - MDisInter
	else if (InterVac.eq.(-1)) then
		MDisVac = MStoVac*(real(ninf,8)/real(Npart,8))
		Conc(N_front:N_front+N_buff) = 0._dp
		Conc = MDisVac*Conc/sum(Conc)
		MStoVac = MStoVac - MDisVac
	end if
	StotoDiscret = Conc
end function



function Sampling(Conc,HLangevin,InterVac)
	implicit none
	integer :: InterVac
	real(dp), dimension(:) :: Conc
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(size(Conc)) :: Conc_buff
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
	Conc_buff = 0._dp
	Conc_buff(N_front:N_front+N_buff) = Conc(N_front:N_front+N_buff)
	Npart = size(HLangevin)
	Mconc = sum(Conc(N_front:N_front+N_buff))
	if (InterVac.eq.0) then
		Nconc = int(Mconc*Npart/MStochastique)
		MStochastique = MStochastique + Mconc
	else if (InterVac.eq.1) then
		Nconc = int(Mconc*Npart/MStoInter)
		MStoInter = MStoInter + Mconc	
	else if (InterVac.eq.(-1)) then
		Nconc = int(Mconc*Npart/MStoVac)
		MStoVac = MStoVac + Mconc		
	end if
	allocate(Hsto(Nconc))
	if (methode.eq.1) then
		Hsto = Multinomial(Conc_buff,N_front,N_front+N_buff,Nconc)
	else 
		Hsto = Metropolis(Conc_buff,real(N_front+N_buff/5._dp,8),real(N_buff/10._dp,8),Nconc)
	end if
	allocate(Htot(Npart+Nconc))
	Htot(1:Npart) = HLangevin(1:Npart)
	Htot(Npart+1:Npart+Nconc) = Hsto(1:Nconc)
	! -- Htot contient *toutes* les particules
	! -- On compte celles qui contribuent a la partie deterministe
	do iloop = 1, Npart+Nconc
		if (Htot(iloop) < N_front) then
			Ninf = Ninf+1
		endif
	end do
	Ndiff = Nconc - Ninf
	! -- Puis on gere le resampling en fonction de Ndiff
	if (Ndiff < 0) then ! -- il faudra dupliquer des particules
	print *, "DUPLICATION"
		do iloop = 1, Npart+Nconc
			if (Htot(iloop).ge.N_front) then
				compt = compt+1
				Sampling(compt) = Htot(iloop)
				nsto = nint(Htot(iloop))
				do kloop = -5, 5
					C_sto(nsto+kloop) = C_sto(nsto+kloop) + 1./(real(Npart,8)*h)*Noyau(1./h*(nsto+kloop-Htot(iloop)))
				end do
			endif
		end do
		do iloop = compt+1,Npart
			call random_number(u)
			nn = int(u*compt)+1
			Sampling(iloop) = Sampling(nn)
			nsto = nint(Sampling(iloop))
			do kloop = -5, 5
				C_sto(nsto+kloop) = C_sto(nsto+kloop) + 1./(real(Npart,8)*h)*Noyau(1./h*(nsto+kloop-Sampling(iloop)))
			end do
		enddo
	else ! -- il faut supprimer des particules
	print *, "SUPPRESSION"
		allocate(Sampling_inter(Npart+Ndiff))
		compt = 0
		do iloop = 1, Npart+Nconc
			if (Htot(iloop).ge.N_front) then
				compt = compt+1
				Sampling_inter(compt) = Htot(iloop)
			endif
		end do	
		compt = 0
		do while (compt < Ndiff)
			call random_number(u)
			nn = int(u*(Npart+Ndiff))+1
			if (Sampling_inter(nn) > 0) then
				Sampling_inter(nn) = -1.
				compt = compt+1
			end if
		end do
		compt = 0
		do iloop = 1, Npart+Ndiff
			if (Sampling_inter(iloop) > 0) then
				compt = compt + 1
				Sampling(compt) = Sampling_inter(iloop)
				nsto = nint(Sampling_inter(iloop))
				do kloop = -5, 5
					C_sto(nsto+kloop) = C_sto(nsto+kloop) + 1./(real(Npart,8)*h)*Noyau(1./h*(nsto+kloop-Sampling_inter(iloop)))
				end do
			end if
		end do
		deallocate(Sampling_inter)
	endif
	if (InterVac.eq.0) then
		C_sto = MStochastique*C_sto/sum(C_sto)
	else if (InterVac.eq.1) then
		C_sto = MStoInter*C_sto/sum(C_sto)	
	else if (InterVac.eq.(-1)) then
		C_sto = MStoVac*C_sto/sum(C_sto)			
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



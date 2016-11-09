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



subroutine Langevin(HLangevin,C1,Tfinal)
	implicit none
	integer :: K
	real(dp) :: Tfinal
	real(dp) :: C1, sqrtdt
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(size(HLangevin)) :: W
	K = int(Tfinal/dt_sto)
	sqrtdt = sqrt(dt_sto)
	do iloop = 1, K
		!print *, "Sto ... ", iloop
		call gen_rand_normal(W,Taille)
		HLangevin(:) = HLangevin(:) + dt_sto*(F_vec(HLangevin(:),C1)) + sqrtdt*sqrt(D_vec(HLangevin(:),C1))*W
	end do
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

subroutine SSA(HLangevin,C1,Tfinal)
	implicit none
	integer :: Npart
	real(dp), dimension(:) :: HLangevin
	real(dp) :: Tfinal, C1
	real(dp) :: dte, nue, expe, pos, u, p
	Npart = size(HLangevin)
	do iloop = 1, Npart
		dte = 0.d0
		do while (dte < Tfinal)
			pos = HLangevin(iloop)
			nue = beta(pos)*C1+alpha(pos)
			call tirage_exp(expe,nue)
			dte = dte + expe
			if (dte < Tfinal) then
				p = beta(pos)*C1/(beta(pos)*C1+alpha(pos))
				call random_number(u)
				if (u.le.p) then 
					HLangevin(iloop) = HLangevin(iloop)+1.
				else 
					HLangevin(iloop) = HLangevin(iloop)-1.
				endif
			endif
		enddo
	enddo
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



function StotoDiscret(HLangevin)
	implicit none
	real(dp), dimension(:) :: HLangevin
	real(dp), dimension(N_front+N_buff) :: Conc, StotoDiscret
	integer :: Npart, n, ninf
	Conc = 0._dp
	ninf = 0
	Npart = size(HLangevin)
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
	MDiscret = MStochastique*(real(ninf,8)/real(Npart,8))
	print *, "Mdiscret = ", MDiscret
	Conc(N_front:N_front+N_buff) = 0._dp
	Conc = MDiscret*Conc/sum(Conc)
	print *, "Mdiscret = ", sum(Conc)
	MStochastique = MStochastique - MDiscret
	StotoDiscret = Conc
end function



function Sampling(Conc,HLangevin)
	implicit none
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
	Nconc = int(Mconc*Npart/MStochastique)
	print *, "Nconc = ", Nconc
	print *, "Mconc = ", Mconc
	print *, "Npart = ", Npart
	print *, "MStochastique = ", MStochastique
	MStochastique = MStochastique + Mconc
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
	C_sto = MStochastique*C_sto/sum(C_sto)
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



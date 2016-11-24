program Couplage

use prec_mod
use parametres
use init_mod
use tools
use stochastique
use fonctioncvode
use mpi

implicit none

integer :: io
real(dp) :: DeltaT = 10._dp
real(dp) :: Mtemp

namelist /parameter/etape,dt_sto,methode,cas_physique

open(8,file="parametre_entree")!, status='OLD', recl=80, delim='APOSTROPHE')

read (8, nml=parameter)


Nmax = 2000

select case (etape)

! Methode explicite !
case(1)

	Neq = 1000
	
	allocate(C(Neq))
	allocate(C_init(Neq))

	call fnvinits(isolver,Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fnvinits'
		stop
	end if

	T = 0._dp
	T0 = 0._dp
	TF = 100._dp
	abstol = 1e-18_dp
	reltol = 1e-5_dp
	C_init(1) = Cq
	C_init(2:Neq) = 0._dp

	call fcvmalloc(T,C_init,meth,itmeth,iatol,reltol,abstol,Iout,Rout,IPAR,RPAR,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvmalloc', IER
		stop
	end if

	maxerrfail = 15
	nitermax = 10000

	call fcvsetiin('MAX_NSTEPS',nitermax,IER)
	call fcvsetiin('MAX_ERRFAIL',maxerrfail,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvsetiin'
		stop
	end if

	call fcvdense(Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvdense', IER
		stop
	end if

	do mainloop = 1,10
		call fcvode(TF,T,C,itask,IER)
		call output(C,T)
		TF = TF + 100._dp
	enddo

	deallocate(C)
	


case(2)


call MPI_INIT(ierror)

call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


	! --- Initialisation du générateur aléatoire 
	call init_random_seed

	! --- Initialisation du solver
	Neq = 250
	allocate(C(Neq))
	allocate(C_init(Neq))
	allocate(C_stotodis(Neq))
	allocate(MPI_C_stotodis(Neq))
	allocate(C_sto(1200))
	
	call fnvinits(isolver,Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fnvinits'
		stop
	end if

	T = 0._dp
	T0 = 0._dp
	TF = 100._dp
	abstol = 1e-18_dp
	reltol = 1e-5_dp
	C_init(1) = Cq
	C_init(2:Neq) = 0._dp

	call fcvmalloc(T0,C_init,meth,itmeth,iatol,reltol,abstol,Iout,Rout,IPAR,RPAR,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvmalloc', IER
		stop
	end if

	maxerrfail = 15
	nitermax = 10000

	call fcvsetiin('MAX_NSTEPS',nitermax,IER)
	call fcvsetiin('MAX_ERRFAIL',maxerrfail,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvsetiin'
		stop
	end if

	call fcvdense(Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvdense', IER
		stop
	end if


	! --- Boucle principale solver+stochastique
	
	! --- On propage en full EDO jusqu'au critere d'arret
	quasi = .False.
	do while (M_queue < Cq/1000._dp) ! 0.1% de la masse totale se trouve dans la queue
		call fcvode(TF,T,C,itask,IER)
		call output(C,T)
		M_queue = Mtail(C,N_front,N_front+N_buff)
		TF = TF + DeltaT
	end do
	
	! --- On genere la partie stochastique et on restart avec la nouvelle CI sur C
	quasi = .True.
	C_init = 0._dp
	C_init(N_front:N_front+N_buff) = C(N_front:N_front+N_buff)
	MStochastique = sum(C_init)

	M_queue = Mtail(C_init,N_front,N_front+N_buff)
	if (methode.eq.1) then 
		Xpart = Multinomial(C_init,N_front,N_front+N_buff,Taille) 
	else 
		Xpart = Metropolis(C_init,N_0,alpha_m,Taille)
	end if
	C_init = 0._dp
	C_init(1:N_front-1) = C(1:N_front-1)
	M_queue = 0._dp
	Cvac(1) = C(1)
	call FCVREINIT(T, C_init, iatol, reltol, abstol, IER)

	! --- Boucle de couplage
	do io = 1, 100
		C_sto = 0._dp

		if (methode.eq.1) then 
			call SSA(Xpart,Cvac,DeltaT)
		else 
			call Langevin(Xpart,Cvac,DeltaT) 
		end if

		call fcvode(TF,T,C,itask,IER)

		DeltaT = 500._dp
		TF = TF + DeltaT
		
		Cvac(1) = Calcul_Cvac(C,Xpart)
		
		call MPI_REDUCE(Cvac,MPI_Cvac,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
		MPI_Cvac = MPI_Cvac/real(numproc,8)
		call MPI_BCAST(MPI_Cvac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
		Cvac = MPI_Cvac


		C_stotodis = StotoDiscret(Xpart,0)
		Xpart = Sampling(C,Xpart,0)
		
		call MPI_REDUCE(C_stotodis,MPI_C_stotodis,Neq,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
		MPI_C_stotodis = MPI_C_stotodis/real(numproc,8)
		call MPI_BCAST(MPI_C_stotodis,Neq,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
		C_stotodis = MPI_C_stotodis
		
		C_init = 0._dp
		C_init(2:N_front-1) = C(2:N_front-1)+C_stotodis(2:N_front-1)
		C_init(1) = Cvac(1)
		C_sto(2:N_front-1) = C_init(2:N_front-1)
		
		if (rank.eq.0) then
			call output(C_sto,T)
		end if
		call FCVREINIT(T, C_init, iatol, reltol, abstol, IER)
	end do
	
	
	deallocate(C)
	deallocate(C_init)
	deallocate(C_stotodis)
	deallocate(MPI_C_stotodis)
	deallocate(C_sto)


	
call MPI_FINALIZE(ierror)	




case(3)

	Ni = 5000
	Nv = 1000
	Neq = 1 + Ni + Nv
	Nmaxi = 5000
	Nmaxv = 1000
	mv = 1
	mi = 1
	! Allocation des principaux tableaux
	allocate(C(Neq))
	allocate(C_init(Neq))
	!allocate(C_stotodis(Neq))
	!allocate(MPI_C_stotodis(Neq))
	!allocate(C_sto(Nmax))
	
	! Allocation et calcul des coefficients alpha/beta
	allocate(Alpha_tab(-Nmaxv:Nmaxi,-mv:mi))
	allocate(Beta_tab(-Nmaxv:Nmaxi,-mv:mi))
	
	do iloop = -Nmaxv,Nmaxi
		do jloop = -mv,mi
			Alpha_tab(iloop,jloop) = alpha_nm(real(iloop,8),real(jloop,8))
			Beta_tab(iloop,jloop) = beta_nm(real(iloop,8),real(jloop,8))
		end do 
	end do
	do iloop = -10,10
		print *, iloop, Alpha_tab(iloop,-1), Alpha_tab(iloop,0), Alpha_tab(iloop,1), &
		&Beta_tab(iloop,-1), Beta_tab(iloop,0), Beta_tab(iloop,1)
	end do

	print *, "init ok"

	call fnvinits(isolver,Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fnvinits'
		stop
	end if
	print *, "init ok 2"
	T = 0._dp
	T0 = 0._dp
	TF = 100._dp
	abstol = 1e-20_dp
	reltol = 1e-5_dp
	C_init(1) = Cq
	C_init(2:Neq) = 0._dp

	call fcvmalloc(T,C_init,meth,itmeth,iatol,reltol,abstol,Iout,Rout,IPAR,RPAR,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvmalloc', IER
		stop
	end if
	print *, "alloc ok"
	maxerrfail = 15
	nitermax = 10000

	call fcvsetiin('MAX_NSTEPS',nitermax,IER)
	call fcvsetiin('MAX_ERRFAIL',maxerrfail,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvsetiin'
		stop
	end if
	print *, "setting ok"
	call fcvdense(Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvdense', IER
		stop
	end if
	print *, "dense ok"
	!TF = 0.001_dp
	TF = 2.58397_dp
	do mainloop = 1,5
		print *, mainloop
		call fcvode(TF,T,C,itask,IER)
		call output(C,T)
		TF = TF*10._dp
	enddo

	deallocate(C)
	
	



case(4)


call MPI_INIT(ierror)

call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)


	! --- Initialisation du générateur aléatoire 
	call init_random_seed
	print *, "Initialisation"
	! --- Initialisation du solver
	Ni = Nf_Inter+Nb_Inter
	Nv = Nf_Vac+Nb_Vac
	Neq = 1 + Ni + Nv
	Nmaxi = 5000
	Nmaxv = 1000
	mv = 1
	mi = 1
	allocate(C(Neq))
	allocate(C_init(Neq))
	allocate(C_stotodis(-Nv:Ni))
	allocate(MPI_C_stotodis(-Nv:Ni))
	allocate(C_sto(-1200:1200))
	allocate(C_sto_Inter(1200))
	allocate(C_sto_Vac(1200))
	allocate(C_Inter(1:Ni))
	allocate(C_Vac(1:Nv))
	allocate(C_Mob(-mv:mi))
	allocate(MPI_C_mob(-mv:mi))
	
	allocate(Alpha_tab(-Nmaxv:Nmaxi,-mv:mi))
	allocate(Beta_tab(-Nmaxv:Nmaxi,-mv:mi))
	
	do iloop = -Nmaxv,Nmaxi
		do jloop = -mv,mi
			Alpha_tab(iloop,jloop) = alpha_nm(real(iloop,8),real(jloop,8))
			Beta_tab(iloop,jloop) = beta_nm(real(iloop,8),real(jloop,8))
		end do 
	end do
	
	do iloop = -20,20
			print *, Alpha_tab(iloop,-1), Alpha_tab(iloop,1), Beta_tab(iloop,-1), Beta_tab(iloop,1)
	end do
	
	call fnvinits(isolver,Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fnvinits'
		stop
	end if

	T = 0._dp
	T0 = 0._dp
	TF = 100._dp
	DeltaT = 100._dp
	abstol = 1e-10_dp
	reltol = 1e-5_dp
	C_init(1) = 0._dp
	C_init(2:Neq) = 0._dp

	call fcvmalloc(T0,C_init,meth,itmeth,iatol,reltol,abstol,Iout,Rout,IPAR,RPAR,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvmalloc', IER
		stop
	end if

	maxerrfail = 15
	nitermax = 10000

	call fcvsetiin('MAX_NSTEPS',nitermax,IER)
	call fcvsetiin('MAX_ERRFAIL',maxerrfail,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvsetiin'
		stop
	end if

	call fcvdense(Neq,IER)
	if (IER /= 0) then
		write(*,*) 'ERROR in fcvdense', IER
		stop
	end if

	print *, "Initialisation : ok"

	! --- Boucle principale solver+stochastique
	
	print *, "Boucle deterministe"
	! --- On propage en full EDO jusqu'au critere d'arret
	quasi = .False.
	do while ((C(tab_fe(Nf_Inter,Nv)) < 1.e3) .and. (C(tab_fe(-Nf_Vac,Nv)) < 1.e7))
		call fcvode(TF,T,C,itask,IER)
		call output(C,T)
		TF = TF + DeltaT
	end do
	if (C(tab_fe(Nf_Inter,Nv)) > 1.e3) then
		print *, "INTERSTITIELS"
		Coupling_Inter = .True.
	end if
	if (C(tab_fe(-Nf_Vac,Nv)) > 1.e7) then
		Coupling_Vac = .True.
		print *, "LACUNES"
	end if
	print *, "Boucle deterministe : ok"
	
	! --- On genere la partie stochastique et on restart avec la nouvelle CI sur C
	
	print *, "Initialisation boucle stochastique"
	quasi = .True.
	
	if (Coupling_Inter) then
		C_Inter = 0._dp
		C_Inter(Nf_Inter:Nf_Inter+Nb_Inter) = C(tab_fe(Nf_Inter,Nv):tab_fe(Nf_Inter+Nb_Inter,Nv))
		MStoInter = sum(C_Inter)
		if (methode.eq.1) then 
			XpartInter = Multinomial(C_Inter,Nf_Inter,Nf_Inter+Nb_Inter,Taille) 
		else 
			N_0 = real(Nf_Inter+Nb_Inter/5._dp,8)
			alpha_m = real(Nb_Inter/10._dp,8)
			XpartInter = Metropolis(C_Inter,N_0,alpha_m,Taille)
		end if
	end if
	if (Coupling_Vac) then		
		C_Vac = 0._dp
		C_Vac(Nf_Vac:Nf_Vac+Nb_Vac) = C(tab_fe(-Nf_Vac,Nv):tab_fe(-Nf_Vac-Nb_Vac,Nv):-1)
		MStoVac = sum(C_Vac)
		if (methode.eq.1) then 
			XpartVac = -1._dp*Multinomial(C_Vac,Nf_Vac,Nf_Vac+Nb_Vac,Taille) 
		else 
			N_0 = real(Nf_Vac+Nb_Vac/5._dp,8)
			alpha_m = real(Nb_Vac/10._dp,8)
			XpartVac = -1._dp*Metropolis(C_vac,N_0,alpha_m,Taille)
		end if
	end if
	
	C_init = 0._dp
	C_Mob = 0._dp
	C_Mob(-mv:mi) = C(tab_fe(-mv,Nv):tab_fe(mi,Nv))
	C_init(tab_fe(-Nf_Vac+1,Nv):tab_fe(Nf_Inter-1,Nv)) = C(tab_fe(-Nf_Vac+1,Nv):tab_fe(Nf_Inter-1,Nv))
	
	call FCVREINIT(T, C_init, iatol, reltol, abstol, IER)

	DeltaT = 100._dp
	print *, "Initialisation boucle stochastique : ok"
	
	! --- Boucle de couplage
	do io = 1, 100
		C_sto = 0._dp
		
		print *, "Calcul Langevin"
		if (Coupling_Inter) then
			if (methode.eq.1) then 
				call SSA(XpartInter,C_Mob,DeltaT)
			else 
				call Langevin(XpartInter,C_Mob,DeltaT) 
			end if
		end if
		if (Coupling_Vac) then
			if (methode.eq.1) then 
				call SSA(XpartVac,C_Mob,DeltaT)
			else 
				call Langevin(XpartVac,C_Mob,DeltaT) 
			end if
		end if
		print *, "Calcul Langevin : ok"

		print *, "Calcul deterministe"
		print *, TF, T
		call fcvode(TF,T,C,itask,IER)
		print *, TF, T
		print *, "Calcul deterministe : ok"
		
		DeltaT = 100._dp
		TF = TF + DeltaT
		
		print *, "Calcul amas mobiles"
		
		C_Mob = Calcul_ConcMob(C,XpartInter,XpartVac,DeltaT/1000._dp,DeltaT)
		
		call MPI_REDUCE(C_Mob,MPI_C_Mob,mv+mi+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
		MPI_C_Mob = MPI_C_Mob/real(numproc,8)
		call MPI_BCAST(MPI_C_Mob,mv+mi+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
		C_Mob = MPI_C_Mob

		print *, "Calcul amas mobiles : ok"

		C_stotodis = 0._dp
	
		if (Coupling_Inter) then
			C_stotodis(1:Nf_Inter+Nb_Inter) = StotoDiscret(XpartInter,1)
			XpartInter = Sampling(C,XpartInter,1)
			C_sto(1:1200) = MStoInter*C_sto(1:1200)
		end if
		if (Coupling_Vac) then
			C_stotodis(-Nf_Vac-Nb_Vac:-1) = StotoDiscret(XpartVac,-1)
			XpartVac = Sampling(C,XpartVac,-1)
			C_sto(-1200:-1) = MStoVac*C_sto(-1200:-1)
		end if

			
		
		call MPI_REDUCE(C_stotodis,MPI_C_stotodis,1+Nf_Inter+Nb_Inter+Nf_Vac+Nb_Vac,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
		MPI_C_stotodis = MPI_C_stotodis/real(numproc,8)
		call MPI_BCAST(MPI_C_stotodis,1+Nf_Inter+Nb_Inter+Nf_Vac+Nb_Vac,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
		C_stotodis = MPI_C_stotodis
		
		C_init = 0._dp
		C_init(tab_fe(mi+1,Nv):tab_fe(Nf_Inter-1,Nv)) = C(tab_fe(mi+1,Nv):tab_fe(Nf_Inter-1,Nv)) + C_stotodis(mi+1:Nf_Inter-1)
		C_init(tab_fe(-Nf_Vac+1,Nv):tab_fe(-mv-1,Nv)) = C(tab_fe(-Nf_Vac+1,Nv):tab_fe(-mv-1,Nv)) + C_stotodis(-Nf_Vac+1:-mv-1)
		C_init(tab_fe(-mv,Nv):tab_fe(mi,Nv)) = C_Mob(-mv:mi)
		C_sto(-Nf_Vac+1:Nf_Inter-1) = C_init(tab_fe(-Nf_Vac+1,Nv):tab_fe(Nf_Inter-1,Nv))
		

		if (C(tab_fe(Nf_Inter,Nv)) > 1.e3) then
			Coupling_Inter = .True.
			print *, "INTERSTITIELS"
		end if
		if (C(tab_fe(-Nf_Vac,Nv)) > 1.e7) then
			Coupling_Vac = .True.
			print *, "LACUNES"
		end if
		
		
		if (rank.eq.0) then
			call output(C_sto,T)
		end if
		call FCVREINIT(T, C_init, iatol, reltol, abstol, IER)
	end do
	
	
	deallocate(C)
	deallocate(C_init)
	deallocate(C_stotodis)
	deallocate(MPI_C_stotodis)
	deallocate(C_sto)
	deallocate(C_Inter)
	deallocate(C_Vac)
	deallocate(C_Mob)
	deallocate(MPI_C_Mob)

	
call MPI_FINALIZE(ierror)	






end select


end program

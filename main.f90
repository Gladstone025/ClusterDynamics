program Couplage

use prec_mod
use parametres
use tools
use stochastique
use fonctioncvode
!use data_fun
!use deterministe
use mpi

implicit none

integer :: io
real(dp) :: DeltaT = 100._dp
real(dp) :: Mtemp

namelist /parameter/etape,dt_sto,methode,model

open(8,file="parametre_entree")!, status='OLD', recl=80, delim='APOSTROPHE')

read (8, nml=parameter)

Neq = 250
Nmax = 2000

select case (etape)

! Methode explicite !
case(1)

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
	Cvac = C(1)
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
		
		Cvac = Calcul_Cvac(C,Xpart)
		
		call MPI_REDUCE(Cvac,MPI_Cvac,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
		MPI_Cvac = MPI_Cvac/real(numproc,8)
		call MPI_BCAST(MPI_Cvac,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
		Cvac = MPI_Cvac


		C_stotodis = StotoDiscret(Xpart)
		Xpart = Sampling(C,Xpart)
		
		call MPI_REDUCE(C_stotodis,MPI_C_stotodis,Neq,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierror)
		MPI_C_stotodis = MPI_C_stotodis/real(numproc,8)
		call MPI_BCAST(MPI_C_stotodis,Neq,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
		C_stotodis = MPI_C_stotodis
		
		C_init = 0._dp
		C_init(2:N_front-1) = C(2:N_front-1)+C_stotodis(2:N_front-1)
		C_init(1) = Cvac
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


end select


end program

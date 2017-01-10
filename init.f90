module init_mod

use prec_mod
use cluster_mod
use parametres

implicit none

CONTAINS


!function tab_fe(pos,debut)
!	implicit none
!	integer :: tab_fe, pos, debut
!	tab_fe = 1 + debut + pos
!end function


subroutine init()
	implicit none
	integer :: nloop, mloop, mobloop, sloop
	type(cluster) :: MonoInter, MonoVac, MonoSol, clust
	Ni = Nf_Inter+Nb_Inter
	Nv = Nf_Vac+Nb_Vac
	Ns = Nf_Sol+Nb_Sol
	
	mi = 1
	mv = 1
	ms = 1
	
	Neq = (1+Ni+Nv)*(1+Ns)
	Nmob = mi+mv+ms
	Nimmob = Neq-Nmob

	allocate(Mob(Nmob))
	MonoInter = cluster(1,0,.True.)
	MonoVac = cluster(-1,0,.True.)
	MonoSol = cluster(0,1,.True.)	
	Mob(1) = MonoInter
	Mob(2) = MonoVac
	Mob(3) = MonoSol
	
	allocate(Det(Neq))
	do nloop = -Nv, Ni
		do sloop = 0, Ns
			clust = cluster(nloop,sloop,.False.)
			if (IsMobile(nloop,sloop)) then
				clust%mobile = .True.
			end if
			Det(C2I(clust)) = clust
		end do	
	end do


end subroutine init





subroutine finalize()
	implicit none
	
	deallocate(C)
	deallocate(C_init)
	!deallocate(C_stotodis)
	!deallocate(MPI_C_stotodis)
	!deallocate(C_sto)


end subroutine finalize





end module init_mod

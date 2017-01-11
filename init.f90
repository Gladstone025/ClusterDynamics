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
	real(dp) :: rn, rm, rmob, rs
	type(cluster) :: MonoInter, MonoVac, MonoSol, clust
	!Ni = Nf_Inter+Nb_Inter
	!Nv = Nf_Vac+Nb_Vac
	!Ns = Nf_Sol+Nb_Sol
	
	!mi = 1
	!mv = 1
	!ms = 1
	
	Neq = (1+Ni+Nv)*(1+Ns)
	Nmob = mi+mv+ms
	Nimmob = Neq-Nmob

	allocate(Mob(Nmob))
	MonoInter = cluster(1,0,.True.)
	MonoInter%ind = C2I(MonoInter)
	MonoVac = cluster(-1,0,.True.)
	MonoVac%ind = C2I(MonoVac)
	!MonoSol = cluster(0,1,.True.)	
	!MonoSol%ind = C2I(MonoSol)
	Mob(1) = MonoInter
	Mob(2) = MonoVac
	!Mob(3) = MonoSol
	
	allocate(Det(Neq))
	allocate(G_source(Neq))
	G_source = 0._dp
	do nloop = -Nv, Ni
		do sloop = 0, Ns
			rn = real(nloop,8)
			rs = real(sloop,8)
			clust = cluster(rn,rs,.False.)
			if (IsMobile(rn,rs)) then
				clust%mobile = .True.
			end if
			clust%ind = C2I(clust)
			Det(clust%ind) = clust
			if (nloop.eq.-8) then
				G_source(clust%ind) = G_8v
			else if (nloop.eq.-1) then
				G_source(clust%ind) = G_v
			else if (nloop.eq.1) then
				G_source(clust%ind) = G_i
			else if (nloop.eq.4) then
				G_source(clust%ind) = G_4i
			end if
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

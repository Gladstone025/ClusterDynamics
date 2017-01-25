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
	integer :: nloop, mloop, mobloop, sloop, iloop
	real(dp) :: rn, rm, rmob, rs
	type(cluster) :: MonoInter, MonoVac, MonoSol, clust
	Ni = Nf_Inter+Nb_Inter
	Nv = Nf_Vac+Nb_Vac
	NS = 0
	!Ns = Nf_Sol+Nb_Sol
	
	mi = 1
	mv = 1
	ms = 0
	!ms = 1
	
	Neq = (1+Ni+Nv)*(1+Ns)
	Nmob = mi+mv+ms
	Nimmob = Neq-Nmob

	allocate(Mob(Nmob))
	allocate(Mob_index(Nmob))
	MonoInter = cluster(1,0,.True.)
	MonoInter%ind = C2I(MonoInter)
	MonoVac = cluster(-1,0,.True.)
	MonoVac%ind = C2I(MonoVac)
	!MonoSol = cluster(0,1,.True.)	
	!MonoSol%ind = C2I(MonoSol)
	Mob(1) = MonoInter
	Mob(2) = MonoVac
	!Mob(3) = MonoSol
	Mob_index(1) = MonoInter%ind
	Mob_index(2) = MonoVac%ind
	!Mob_index(3) = MonoSol%ind
	
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
			Det(nint(clust%ind)) = clust
			if (nloop.eq.-8) then
				G_source(nint(clust%ind)) = G_8v
			else if (nloop.eq.-1) then
				G_source(nint(clust%ind)) = G_v
			else if (nloop.eq.1) then
				G_source(nint(clust%ind)) = G_i
			else if (nloop.eq.4) then
				G_source(nint(clust%ind)) = G_4i
			end if
		end do	
	end do
	

	NbuffI = Nb_Inter+1
	allocate(BuffI_Index(NbuffI))
	iloop = 0
	do nloop = 0, Ni
		rn = real(nloop,8)
		clust = c_constructor(rn,0._dp)
		if (rn.le.Ni .and. rn.ge.Nf_Inter) then
			iloop = iloop+1
			BuffI_Index(iloop) = nint(clust%ind)
		end if
	end do
	print *, "NbuffI", iloop, NbuffI

	NbuffV = Nb_Vac+1 
	allocate(BuffV_Index(NbuffV))
	iloop = 0
	do nloop = -Nv, -1
		rn = real(nloop,8)
		clust = c_constructor(rn,0._dp)
		if (rn.ge.-Nv .and. rn.le.-Nf_Vac) then
			iloop = iloop+1
			BuffV_Index(iloop) = nint(clust%ind)
		end if
	end do
	print *, "NbuffV", iloop, NbuffV
		
	
	Nbound = 1
	allocate(FrontI_Index(Nbound))
	allocate(FrontV_Index(Nbound))
	clust = c_constructor(real(Nf_Inter,8),0._dp)
	FrontI_Index(Nbound) = int(clust%ind)
	clust = c_constructor(real(Nf_Vac,8),0._dp)
	FrontV_Index(Nbound) = int(clust%ind)


end subroutine init





subroutine init_FeHe

implicit none
	integer :: nloop, mloop, mobloop, sloop, iloop
	real(dp) :: rn, rm, rmob, rs
	type(cluster) :: MonoInter, MonoVac, MonoSol, clust
	NbuffI = Nb_Inter+1 !(Nb_Inter+1)*Nf_Sol + (Nb_Sol+1)*(Ni+1)
		allocate(BuffI_Index(NbuffI))
		iloop = 0
		do nloop = 0, Ni
			do sloop = 0, Ns
				rn = real(nloop,8)
				rs = real(sloop,8)
				clust = cluster(rn,rs,.False.)
				if (IsMobile(rn,rs)) then
					clust%mobile = .True.
				end if
				clust%ind = C2I(clust)
				if (rn.le.Ni .and. rn.ge.Nf_Inter .and. rs < Nf_Sol) then
					iloop = iloop+1
					BuffI_Index(iloop) = nint(clust%ind)
				else if (rs.le.Ns .and. rs.ge.Nf_Sol .and. rn.le.Ni) then
					iloop = iloop+1
					BuffI_Index(iloop) = nint(clust%ind)
				end if
			end do	
		end do
		print *, "NbuffI", iloop, NbuffI

		NbuffV = Nb_Vac+1 !(Nb_Vac+1)*Nf_Sol + (Nb_Sol+1)*Nv
		allocate(BuffV_Index(NbuffV))
		iloop = 0
		do nloop = -Nv, -1
			do sloop = 0, Ns
				rn = real(nloop,8)
				rs = real(sloop,8)
				clust = cluster(rn,rs,.False.)
				if (IsMobile(rn,rs)) then
					clust%mobile = .True.
				end if
				clust%ind = C2I(clust)
				if (rn.ge.-Nv .and. rn.le.Nf_Vac .and. rs < Nf_Sol) then
					iloop = iloop+1
					BuffV_Index(iloop) = nint(clust%ind)
				else if (rs.le.Ns .and. rs.ge.Nf_Sol .and. rn.ge.Nv) then
					iloop = iloop+1
					BuffV_Index(iloop) = nint(clust%ind)
				end if
			end do	
		end do
		print *, "NbuffV", iloop, NbuffV


end subroutine init_FeHe



subroutine finalize()
	implicit none
	
	deallocate(C)
	deallocate(C_init)
	!deallocate(C_stotodis)
	!deallocate(MPI_C_stotodis)
	!deallocate(C_sto)


end subroutine finalize





end module init_mod

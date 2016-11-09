module init_mod

use prec_mod
use parametres

implicit none

CONTAINS


function tab_fe(pos,debut)
	implicit none
	integer :: tab_fe, pos, debut
	tab_fe = 1 + debut + pos
end function


subroutine init()
	implicit none
	
	Ni = 3000
	Nv = 100
	Neq = 1 + Ni + Nv
	Nmaxi = 3000
	Nmaxv = 100
	mv = 1
	mi = 1
	! Allocation des principaux tableaux
	allocate(C(Neq))
	allocate(C_init(Neq))
	allocate(C_stotodis(Neq))
	allocate(MPI_C_stotodis(Neq))
	allocate(C_sto(Nmax))
	
	! Allocation et calcul des coefficients alpha/beta
	allocate(Alpha_tab(-Nmaxv:Nmaxi,-mv:mi))
	allocate(Beta_tab(-Nmaxv:Nmaxi,-mv:mi))
	
	do iloop = -Nmaxv,Nmaxi
		do jloop = -mv,mi
			Alpha_tab(iloop,jloop) = alpha_nm(real(iloop,8),real(jloop,8))
			Beta_tab(iloop,jloop) = beta_nm(real(iloop,8),real(jloop,8))
		end do 
	end do


end subroutine init





subroutine finalize()
	implicit none
	
	deallocate(C)
	deallocate(C_init)
	deallocate(C_stotodis)
	deallocate(MPI_C_stotodis)
	deallocate(C_sto)


end subroutine finalize





end module init_mod

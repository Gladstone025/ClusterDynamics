module parametres


use prec_mod
use cluster_mod
use para_Ov
use para_Fe

implicit none


!! ------------------------- Boundary and model parameters ------------------------- !!
integer, parameter :: N_front = 200
integer, parameter :: N_buff = 100

integer, parameter :: Nf_Inter = 400
integer, parameter :: Nb_Inter = 400

integer, parameter :: Nf_Vac = 200
integer, parameter :: Nb_Vac = 50

integer, parameter :: Nf_Sol = 200
integer, parameter :: Nb_Sol = 50

type(cluster), dimension(:), allocatable :: Mob, Det, Immob
real(dp), dimension(:), allocatable :: G_source

integer, parameter :: Nf_EL_Inter = 50
integer, parameter :: Nf_EL_Vac = 50
integer, parameter :: Nf_EL_Sol = 50

real(dp), parameter :: R_min = -1.e10

real(dp), dimension(:), allocatable :: Mob_i, Det_i, Immob_i

real(dp) :: M_queue = 0 
real(dp), dimension(1) :: Cvac
real(dp) :: CmobVac, CmobInter
logical :: quasi, quasi1
integer :: etape, methode, cas_physique, Nmax
logical :: Coupling_Inter = .False. 
logical :: Coupling_Vac = .False.

!! --------------------------- Parametres CVODE ----------------------------- !!
real(dp) :: T, T0, TF, DeltaT, TF1, TF2, T2
integer(c_int) :: IER
real(dp), dimension(:), allocatable :: C, C_init, C_stotodis, C_init1, C_init2
real(dp), dimension(:), allocatable :: C_Inter, C_Vac, C_Mob, bCn_sto, aCmob_sto
real(dp), dimension(1) :: RPAR
integer(c_long), dimension(1) :: IPAR
real(dp) :: tt, tout
real(dp) :: abstol, reltol

integer, parameter :: isolver = 1 
integer, parameter :: meth = 2
integer, parameter :: itmeth = 2
integer, parameter :: iatol = 1
integer, parameter :: itask = 1
integer(c_long), dimension(21) :: iout
real(dp), dimension(6) :: rout
integer(c_long) :: maxerrfail, nitermax
integer(c_long) :: Neq, Nmob, Nimmob
!integer :: Nv, Ni, Ns
!integer :: mv, mi, ms
integer :: Nmaxv, Nmaxi, Nmaxs
real(dp), dimension(:,:), allocatable :: Alpha_tab, Beta_tab
real(dp), dimension(:,:,:), allocatable :: Alpha_tab_Sol, Beta_tab_Sol

!! ----------------------- Parametres stochastiques ------------------------- !!
integer, parameter :: Taille = 50000
real(dp), dimension(Taille) :: Xpart = 0._dp
real(dp), dimension(Taille) :: XpartInter = 0._dp
real(dp), dimension(Taille) :: XpartVac = 0._dp
real(dp), dimension(:), allocatable :: C_sto

real(dp) :: bCnCi_sto, aCi_sto, Si_sto
real(dp) :: bCnCv_sto, aCv_sto, Sv_sto

real(dp) :: N_0 = real(N_front+N_buff/5._dp,8)
real(dp) :: alpha_m = real(N_buff/10._dp,8)
real(dp) :: dt_sto
real(dp) :: MStochastique, MDiscret
real(dp) :: MStoInter, MStoVac, MDisInter, MDisVac

!! ---------------------------- Parametres MPI ------------------------------ !!
integer :: rank, numproc, ierror
real(dp) :: MPI_Cvac
real(dp), dimension(:), allocatable :: MPI_C_Mob
real(dp), dimension(:), allocatable :: MPI_C_stotodis

!! --------------------------- Parametres utiles ---------------------------- !!
real(dp) :: rloop
real(dp), parameter :: h = 0.4_dp


CONTAINS


function Alpha_Clust(c_nu,c_mu)
	type(cluster) :: c_nu, c_mu
	real(dp) :: Alpha_Clust
	!select case(cas_physique)
	!case(3)
		Alpha_Clust = Alpha_Clust_Fe(c_nu,c_mu)
	!end select
end function

function Beta_Clust(c_nu,c_mu)
	type(cluster) :: c_nu, c_mu
	real(dp) :: Beta_Clust
	!select case(cas_physique)
	!case(3)
		Beta_Clust = Beta_Clust_Fe(c_nu,c_mu)
	!end select
end function



end module parametres




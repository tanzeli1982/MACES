module data_buffer_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the buffer of large temporary data in the model
!
!---------------------------------------------------------------------------------
   use data_type_mod

   implicit none
   public

   !! allocatable arrays for RungeKutta4 
   real(kind=8), allocatable, dimension(:,:) :: rk4_K1
   real(kind=8), allocatable, dimension(:,:) :: rk4_K2
   real(kind=8), allocatable, dimension(:,:) :: rk4_K3
   real(kind=8), allocatable, dimension(:,:) :: rk4_K4
   real(kind=8), allocatable, dimension(:,:) :: rk4_K5
   real(kind=8), allocatable, dimension(:,:) :: rk4_K6
   real(kind=8), allocatable, dimension(:,:) :: rk4_nxt4th
   real(kind=8), allocatable, dimension(:,:) :: rk4_nxt5th
   real(kind=8), allocatable, dimension(:,:) :: rk4_interim
   real(kind=8), allocatable, dimension(:,:) :: rk4_rerr
   
   !! allocatable arrays for Hydrodynamic model
   ! hydrodynamic state variables [h, U*h, N, Css, Cj]
   real(kind=8), allocatable, dimension(:,:) :: m_uhydro
   ! platform coordinate and elevation (m) 
   real(kind=8), allocatable, dimension(:)   :: m_X
   real(kind=8), allocatable, dimension(:)   :: m_dX
   real(kind=8), allocatable, dimension(:)   :: m_Zh
   real(kind=8), allocatable, dimension(:)   :: m_dZh
   ! current speed (m/s)
   real(kind=8), allocatable, dimension(:)   :: m_U
   ! significant wave height (m) and wave energy (W)
   real(kind=8), allocatable, dimension(:)   :: m_Hwav
   real(kind=8), allocatable, dimension(:)   :: m_Ewav
   real(kind=8), allocatable, dimension(:)   :: m_Uwav
   ! bottom shear stress (Pa)
   real(kind=8), allocatable, dimension(:)   :: m_tau
   ! sediment source and sink (kg/m2/s)
   real(kind=8), allocatable, dimension(:)   :: m_Esed
   real(kind=8), allocatable, dimension(:)   :: m_Dsed
   real(kind=8), allocatable, dimension(:)   :: m_Cz
   real(kind=8), allocatable, dimension(:)   :: m_kwav
   real(kind=8), allocatable, dimension(:)   :: m_Qb
   real(kind=8), allocatable, dimension(:)   :: m_Swg
   real(kind=8), allocatable, dimension(:)   :: m_Sbf
   real(kind=8), allocatable, dimension(:)   :: m_Swc
   real(kind=8), allocatable, dimension(:)   :: m_Sbrk
   ! temporary variables
   real(kind=8), allocatable, dimension(:,:) :: tmp_uhydroL
   real(kind=8), allocatable, dimension(:,:) :: tmp_uhydroR
   real(kind=8), allocatable, dimension(:,:) :: tmp_phi
   real(kind=8), allocatable, dimension(:,:) :: tmp_FL
   real(kind=8), allocatable, dimension(:,:) :: tmp_FR
   real(kind=8), allocatable, dimension(:,:) :: tmp_P
   real(kind=8), allocatable, dimension(:,:) :: tmp_SRC
   real(kind=8), allocatable, dimension(:,:) :: tmp_eigval
   real(kind=8), allocatable, dimension(:)   :: tmp_aL
   real(kind=8), allocatable, dimension(:)   :: tmp_aR
   real(kind=8), allocatable, dimension(:)   :: tmp_U
   real(kind=8), allocatable, dimension(:)   :: tmp_Cg
   real(kind=8), allocatable, dimension(:)   :: tmp_Nmax
   real(kind=8), allocatable, dimension(:)   :: tmp_Qb

   !! user-defined data
   type(ModelParams) :: m_params
   type(ForcingData) :: m_forcings

end module data_buffer_mod

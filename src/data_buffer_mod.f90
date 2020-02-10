module data_buffer_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the buffer of large temporary data in the model
!
!---------------------------------------------------------------------------------
   implicit none
   public

   integer, parameter :: Wss = 1, Wsal = 2
   integer, parameter :: EQM_WAVE = 1, TSNT_WAVE = 2

   !! allocatable arrays for Hydrodynamic model
   ! hydrodynamic state variables [h, U*h, Css, Cj]
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
   ! wave source and sink
   real(kind=8), allocatable, dimension(:)   :: m_Cz
   real(kind=8), allocatable, dimension(:)   :: m_kwav
   real(kind=8), allocatable, dimension(:)   :: m_Qb
   real(kind=8), allocatable, dimension(:)   :: m_Swg
   real(kind=8), allocatable, dimension(:)   :: m_Sbf
   real(kind=8), allocatable, dimension(:)   :: m_Swc
   real(kind=8), allocatable, dimension(:)   :: m_Sbrk
   ! sediment and constituents
   real(kind=8), allocatable, dimension(:,:) :: m_Cs
   ! temporary variables
   real(kind=8), allocatable, dimension(:,:) :: tmp_uhydro
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
   real(kind=8), allocatable, dimension(:)   :: tmp_B
   ! forcing variables
   real(kind=8), allocatable, dimension(:,:) :: Cs_source
   real(kind=8), allocatable, dimension(:,:) :: Cs_sink
   real(kind=8), allocatable, dimension(:)   :: fctr_wave
   real(kind=8) :: frc_Twav
   real(kind=8) :: frc_U10
   ! parameter variables
   real(kind=8), allocatable, dimension(:)   :: par_alphaA
   real(kind=8), allocatable, dimension(:)   :: par_betaA
   real(kind=8), allocatable, dimension(:)   :: par_alphaD
   real(kind=8), allocatable, dimension(:)   :: par_betaD
   real(kind=8), allocatable, dimension(:)   :: par_cD0
   real(kind=8), allocatable, dimension(:)   :: par_ScD
   real(kind=8) :: par_d50
   real(kind=8) :: par_Cz0
   real(kind=8) :: par_Kdf
   real(kind=8) :: par_cbc
   real(kind=8) :: par_cwc
   real(kind=8) :: par_fr
   ! rungekutta temporary arrays
   real(kind=8), allocatable, dimension(:,:) :: rk4_K1(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_K2(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_K3(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_K4(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_K5(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_K6(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_nxt4th(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_nxt5th(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_interim(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rk4_rerr(:,:)

end module data_buffer_mod

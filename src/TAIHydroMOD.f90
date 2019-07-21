module TAIHydroMOD
!---------------------------------------------------------------------------------
! Purpose:
!
! This module implements the 1-D transect-based hydrodynamic model
!
!---------------------------------------------------------------------------------
   use RungeKutta,   only : InitRKDataBuffer, DestructRKDataBuffer 

   type ParamData
      real(kind=8) :: d50     ! sediment median diameter (m)
      real(kind=8) :: C0      ! flow conductance on tidal flat
      real(kind=8) :: K       ! sediment dispersion coefficient (m2/s)
      real(kind=8) :: cbc     ! wave dissipation coefficient
      real(kind=8) :: fr      ! maximum allowed wave height fraction
      real(kind=8) :: cb      !

   end type

   type ForcingData
      real(kind=8) :: T       ! wave period (s)
      real(kind=8) :: U10     ! wind speed (m/s)
      real(kind=8) :: U0      ! seaward current speed (m/s)
      real(kind=8) :: h0      ! seaward water depth (m)
      real(kind=8) :: Uw0     ! seaward wave speed (m/s)
   end type

   implicit none
   integer, parameter :: NVAR = 5
   ! hydrodynamic state variables
   real(kind=8), allocatable, dimension(:,:) :: m_uhydro
   ! 
   real(kind=8), allocatable, dimension(:)   :: m_X
   real(kind=8), allocatable, dimension(:)   :: m_dX
   real(kind=8), allocatable, dimension(:)   :: m_Zh
   real(kind=8), allocatable, dimension(:)   :: m_U
   real(kind=8), allocatable, dimension(:)   :: m_Hwav
   real(kind=8), allocatable, dimension(:)   :: m_Ewav
   real(kind=8), allocatable, dimension(:)   :: m_tau
   real(kind=8), allocatable, dimension(:)   :: m_Cf
   real(kind=8), allocatable, dimension(:)   :: m_kwav
   real(kind=8), allocatable, dimension(:)   :: m_Qb
   real(kind=8), allocatable, dimension(:)   :: m_Swg
   real(kind=8), allocatable, dimension(:)   :: m_Sbf
   real(kind=8), allocatable, dimension(:)   :: m_Swc
   real(kind=8), allocatable, dimension(:)   :: m_Sbrk
   ! forcing data
   type(ForcingData) :: m_forcings
   ! parameters
   type(ParamData) :: m_params
   ! other variables
   integer :: NX

contains
   subroutine InitializeTAIMOD(X, Zh)
      implicit none
      real(kind=8), intent(in) :: X(:)    ! grid cell coordinate (m) 
      real(kind=8), intent(in) :: Zh(:)   ! grid cell elevation (m)
      integer :: ii

      NX = size(X)
      call InitRKDataBuffer(NVAR, NX)
      allocate(m_uhydro(NVAR,NX))
      allocate(m_X(NX))
      allocate(m_dX(NX))
      allocate(m_Zh(NX))
      allocate(m_U(NX))
      allocate(m_Hwav(NX))
      allocate(m_Ewav(NX))
      allocate(m_tau(NX))
      allocate(m_Qb(NX))
      allocate(m_Swg(NX))
      allocate(m_Sbf(NX))
      allocate(m_Swc(NX))
      allocate(m_Sbrk(NX))
      
      m_uhydro = 0.0
      do ii = 1, NX, 1
         if (ii==1) then
            m_dX(ii) = 0.5 * (m_X(ii) + m_X(ii+1))
         else if (ii==NX) then
            m_dX(ii) = 0.5 * (m_X(ii-1) + m_X(ii))
         else
            m_dX(ii) = 0.5 * (m_X(ii+1) - m_X(ii-1))
         end if
         m_uhydro(3,ii) = max(-Zh(ii), 0.0)
      end do
   end subroutine

   subroutine DestructTAIMOD()
      implicit none

      call DestructRKDataBuffer()
      deallocate(m_uhydro)
      deallocate(m_X)
      deallocate(m_dX)
      deallocate(m_Zh)
      deallocate(m_U)
      deallocate(m_Hwav)
      deallocate(m_Ewav)
      deallocate(m_tau)
      deallocate(m_Qb)
      deallocate(m_Swg)
      deallocate(m_Sbf)
      deallocate(m_Swc)
      deallocate(m_Sbrk)
   end subroutine

   subroutine ModelSetup(Zh, Bag, T, U0, h0, Uw0, U10)
      implicit none
      real(kind=8), intent(in) :: Zh(:)   ! platform elevation (m)
      real(kind=8), intent(in) :: Bag(:)  ! aboveground biomass (kg/m2)
      real(kind=8), intent(in) :: T
      real(kind=8), intent(in) :: U0
      real(kind=8), intent(in) :: h0
      real(kind=8), intent(in) :: Uw0
      real(kind=8), intent(in) :: U10

      m_Zh = Zh
      m_forcings%T = T
      m_forcings%U0 = U0
      m_forcings%h0 = h0
      m_forcings%Uw0 = Uw0
      m_forcings%U10 = U10
      call UpdateGroundRoughness(Bag)
      call UpdateWaveNumber()
      call UpdateWaveBrkPSB()
      call UpdateWaveBtmFriction()
      call UpdateWaveWhiteBrking()
      call UpdateWaveDepthBrking()
   end subroutine

   subroutine ModelCallback()
      implicit none

      call UpdateShearStress()
   end subroutine

   subroutine UpdateGroundRoughness(Bag)
      implicit none
      real(kind=8), intent(in) :: Bag(:)  ! aboveground biomass (kg/m2)
      real(kind=8) :: cb, C0, cD0, ScD
      real(kind=8) :: alphaA, alphaD
      real(kind=8) :: betaA, betaD
      real(kind=8) :: asb, dsb, cD, h
      integer :: ii

      cb = m_params%cb
      C0 = m_params%C0
      cD0 = m_params%cD0
      ScD = m_params%ScD
      alphaA = m_params%alphaA
      alphaB = m_params%alphaB
      betaA = m_params%betaA
      betaB = m_params%betaB
      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h<0.1) then
            m_Cf(ii) = 1d-20
         else
            asb = alphaA * Bag(ii)**betaA
            dsb = alphaD * Bag(ii)**betaD
            cD = cD0 + ScD * Bag(ii)
            m_Cf(ii) = C0*sqrt(2.0/(cD*(cb**2)*asb*h+2.0*(1-asb*dsb)))
         end if
      end do
   end subroutine

   subroutine UpdateWaveNumber()
      implicit none

   end subroutine

end module TAIHydroMOD

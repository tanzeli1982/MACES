module TAIHydroMOD
!---------------------------------------------------------------------------------
! Purpose:
!
! This module implements the 1-D transect-based hydrodynamic model
!
!---------------------------------------------------------------------------------
   use RungeKutta4,  only : InitRKDataBuffer, DestructRKDataBuffer 
   use RungeKutta4,  only : NonLRBrents, BinarySearch 
   use RungeKutta4,  only : FVSKT_Superbee, FVSKT_celledge

   implicit none
   integer, parameter :: NVAR = 5
   integer, parameter :: NPFT = 4
   real(kind=8), parameter :: PI = 3.14159265d+0
   real(kind=8), parameter :: e = 2.71828183d+0
   real(kind=8), parameter :: Roul = 1028.0
   real(kind=8), parameter :: Roua = 1.225
   real(kind=8), parameter :: Karman = 0.41
   real(kind=8), parameter :: G = 9.8
   ! model parameters
   real(kind=8) :: par_d50       ! sediment median diameter (m)
   real(kind=8) :: par_Cz0       ! the Chézy friction coefficient (m^0.5/s)
   real(kind=8) :: par_cD0(NPFT) ! bulk plant drag coefficient baseline (unitless)
   real(kind=8) :: par_ScD(NPFT) ! the slope between cD and Bag (unknown)
   real(kind=8) :: par_alphaA(NPFT) ! empirical coefficient for projected plant area
   real(kind=8) :: par_betaA(NPFT)  ! empirical coefficient for projected plant area
   real(kind=8) :: par_alphaD(NPFT) ! empirical coefficient for stem diameter
   real(kind=8) :: par_betaD(NPFT)  ! empirical coefficient for stem diameter
   real(kind=8) :: par_Kdf       ! sediment dispersion coefficient (m^2/s)
   real(kind=8) :: par_cbc       ! wave dissipation coefficient (unknown)
   real(kind=8) :: par_fr        ! maximum allowed wave height fraction
   ! model inputs
   real(kind=8) :: force_T       ! wave period (s)
   real(kind=8) :: force_U10     ! wind speed (m/s)
   ! hydrodynamic state variables
   real(kind=8), allocatable, dimension(:,:) :: m_uhydro
   ! platform coordinate and elevation (m) 
   real(kind=8), allocatable, dimension(:)   :: m_X
   real(kind=8), allocatable, dimension(:)   :: m_dX
   real(kind=8), allocatable, dimension(:)   :: m_Zh
   real(kind=8), allocatable, dimension(:)   :: m_dZh
   integer, allocatable, dimension(:)        :: m_pft
   ! current speed (m/s)
   real(kind=8), allocatable, dimension(:)   :: m_U
   ! significant wave height (m) and wave energy (W)
   real(kind=8), allocatable, dimension(:)   :: m_Hwav
   real(kind=8), allocatable, dimension(:)   :: m_Ewav
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
   ! constants
   real(kind=8), allocatable, dimension(:)   :: const_Qb
   ! other variables
   integer :: NX

contains
   subroutine Initialize(X, Zh)
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
      allocate(m_dZh(NX))
      allocate(m_pft(NX))
      allocate(m_U(NX))
      allocate(m_Hwav(NX))
      allocate(m_Ewav(NX))
      allocate(m_Esed(NX))
      allocate(m_Dsed(NX))
      allocate(m_tau(NX))
      allocate(m_Qb(NX))
      allocate(m_Swg(NX))
      allocate(m_Sbf(NX))
      allocate(m_Swc(NX))
      allocate(m_Sbrk(NX))
      allocate(m_Cz(NX))
      allocate(tmp_uhydroL(NVAR,NX))
      allocate(tmp_uhydroR(NVAR,NX))
      allocate(tmp_phi(NVAR,NX))
      allocate(tmp_FL(NVAR,NX))
      allocate(tmp_FR(NVAR,NX))
      allocate(tmp_P(NVAR,NX))
      allocate(tmp_SRC(NVAR,NX))
      allocate(tmp_eigval(NVAR,NX))
      allocate(tmp_aL(NX))
      allocate(tmp_aR(NX))
      allocate(tmp_U(NX))
      allocate(tmp_Cg(NX))
      allocate(tmp_Nmax(NX))
      allocate(const_Qb(101))
      
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
      const_Qb = (/0.0,0.4637,0.5005,0.5260,0.5461,0.5631,0.5780,0.5914, &
         0.6035,0.6147,0.6252,0.6350,0.6442,0.6530,0.6614,0.6694,0.6770, &
         0.6844,0.6915,0.6984,0.7050,0.7115,0.7177,0.7238,0.7298,0.7355, &
         0.7412,0.7467,0.7521,0.7573,0.7625,0.7676,0.7725,0.7774,0.7822, &
         0.7869,0.7915,0.7960,0.8005,0.8049,0.8092,0.8135,0.8177,0.8218, &
         0.8259,0.8299,0.8339,0.8378,0.8417,0.8455,0.8493,0.8531,0.8568, &
         0.8604,0.8640,0.8676,0.8711,0.8746,0.8781,0.8815,0.8849,0.8883, &
         0.8916,0.8949,0.8981,0.9014,0.9046,0.9078,0.9109,0.9140,0.9171, &
         0.9202,0.9232,0.9262,0.9292,0.9322,0.9352,0.9381,0.9410,0.9439, &
         0.9467,0.9496,0.9524,0.9552,0.9580,0.9607,0.9635,0.9662,0.9689, &
         0.9716,0.9742,0.9769,0.9795,0.9821,0.9847,0.9873,0.9899,0.9924, &
         0.9950,0.9975,1.0/)
   end subroutine

   subroutine Destruct()
      implicit none

      call DestructRKDataBuffer()
      deallocate(m_uhydro)
      deallocate(m_X)
      deallocate(m_dX)
      deallocate(m_Zh)
      deallocate(m_dZh)
      deallocate(m_pft)
      deallocate(m_U)
      deallocate(m_Hwav)
      deallocate(m_Ewav)
      deallocate(m_Esed)
      deallocate(m_Dsed)
      deallocate(m_tau)
      deallocate(m_Qb)
      deallocate(m_Swg)
      deallocate(m_Sbf)
      deallocate(m_Swc)
      deallocate(m_Sbrk)
      deallocate(m_Cz)
      deallocate(tmp_uhydroL)
      deallocate(tmp_uhydroR)
      deallocate(tmp_phi)
      deallocate(tmp_FL)
      deallocate(tmp_FR)
      deallocate(tmp_P)
      deallocate(tmp_SRC)
      deallocate(tmp_eigval)
      deallocate(tmp_aL)
      deallocate(tmp_aR)
      deallocate(tmp_U)
      deallocate(tmp_Cg)
      deallocate(tmp_Nmax)
      deallocate(const_Qb)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Set model parameters.
   ! cbc should be related to 1/Cz0
   ! d50 = 0.00025   ! m
   ! Cz0 = 65.0      ! m^0.5/s
   ! Kdf = 100.0     ! m^2/s
   ! cbc = 0.015     ! unknown 
   ! fr = 0.78       ! fraction
   ! alphaA(:) = 8.0
   ! betaA(:) = 0.5
   ! alphaD(:) = 0.005
   ! betaD(:) = 0.3
   ! cD0(:) = 1.1
   ! ScD(:) = -0.3
   !
   !------------------------------------------------------------------------------
   subroutine SetModelParameters(d50, Cz0, Kdf, cbc, fr, alphaA, betaA, &
                                 alphaD, betaD, cD0, ScD)
      implicit none
      real(kind=8), intent(in) :: d50
      real(kind=8), intent(in) :: Cz0
      real(kind=8), intent(in) :: Kdf
      real(kind=8), intent(in) :: cbc
      real(kind=8), intent(in) :: fr
      real(kind=8), intent(in) :: alphaA(:)
      real(kind=8), intent(in) :: betaA(:)
      real(kind=8), intent(in) :: alphaD(:)
      real(kind=8), intent(in) :: betaD(:)
      real(kind=8), intent(in) :: cD0(:)
      real(kind=8), intent(in) :: ScD(:)

      par_d50 = d50
      par_Cz0 = Cz0
      par_Kdf = Kdf
      par_cbc = cbc
      par_fr = fr
      par_alphaA = alphaA
      par_betaA = betaA
      par_alphaD = alphaD
      par_betaD = betaD
      par_cD0 = cD0
      par_ScD = ScD
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Pre-run model intermediate variable updates and read boundary 
   !          conditions.
   !
   !------------------------------------------------------------------------------
   subroutine ModelSetup(Zh, pft, Bag, Ero, Dep, T, U10, h0, U0, Hwav0, &
                         Css0, Cj0)
      implicit none
      real(kind=8), intent(in) :: Zh(:)   ! platform elevation (m)
      integer, intent(in) :: pft(:)       ! platform vegetation type
      real(kind=8), intent(in) :: Bag(:)  ! aboveground biomass (kg/m2)
      real(kind=8), intent(in) :: Ero(:)  ! sediment suspension (kg/m2/s)
      real(kind=8), intent(in) :: Dep(:)  ! sediment deposition (kg/m2/s)
      real(kind=8), intent(in) :: T       ! wave period (s)
      real(kind=8), intent(in) :: U10     ! wind speed (m/s)
      real(kind=8), intent(in) :: h0      ! seaward water depth (m)
      real(kind=8), intent(in) :: U0      ! seaward current speed (m/s)
      real(kind=8), intent(in) :: Hwav0   ! seaward significant wave height (m)
      real(kind=8), intent(in) :: Css0    ! seaward sediment conc (kg/m3)
      real(kind=8), intent(in) :: Cj0     ! seaward salinity (PSU)
      integer :: ii

      m_Zh = Zh
      m_Esed = Ero
      m_Dsed = Dep
      ! a typical wave period is 2s but increase greatly with the increase
      ! of wave speed (https://en.wikipedia.org/wiki/Wind_wave)
      force_T = T
      force_U10 = U10
      m_uhydro(1,1) = h0
      m_uhydro(2,1) = U0 * h0
      m_uhydro(3,1) = 0.125*Roul*G*(Hwav0**2)/(2.0*PI/T)
      m_uhydro(4,1) = Css0
      m_uhydro(5,1) = Cj0 
      do ii = 1, NX, 1
         if (ii==1) then
            m_dZh(ii) = 0.5 * (m_Zh(ii+1) - m_Zh(ii))
         else if (ii==NX) then
            m_dZh(ii) = 0.5 * (m_Zh(ii) - m_Zh(ii-1))
         else
            m_dZh(ii) = 0.5 * (m_Zh(ii+1) - m_Zh(ii-1))
         end if
      end do
      call UpdateGroundRoughness(pft, Bag)
      !call UpdateWaveNumber()
      call UpdateWaveNumber2()   ! much more efficient but less accurate
      !call UpdateWaveBrkProb()
      call UpdateWaveBrkProb2()  ! much more efficient but less accurate
      call UpdateWaveGeneration()
      call UpdateWaveBtmFriction()
      call UpdateWaveWhiteCapping()
      call UpdateWaveDepthBrking()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Post-run model update.
   !          At depths smaller than 0.1 m, the flow velocity and sediment 
   !          transport are taken to zero (linearly interpolated to zero)??
   !
   !------------------------------------------------------------------------------
   subroutine ModelCallback()
      implicit none
      real(kind=8) :: sigma
      integer :: ii
      
      sigma = 2.0*PI/force_T
      do ii = 1, NX, 1
         if (m_uhydro(1,ii)<=0) then
            m_uhydro(2:NVAR,ii) = 0.0d0
         end if
      end do
      m_U = m_uhydro(2,:) / max(0.1,m_uhydro(1,:))
      m_Ewav = sigma * m_uhydro(3,:)
      m_Hwav = sqrt(8.0*m_Ewav/G/Roul)
      call UpdateShearStress()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate ground roughness and bottom shear stress
   ! Chézy's coefficient vs manning's coefficient: Cz = h^(1/6) / n
   ! Chézy's coefficient vs Darcy-Weisbach friction factor: Cz = sqrt(8*G/f)
   ! From an empirical equation of Froehlich (2012; J. Irrigation and Drainage
   ! Engineering), Cz0 = sqrt(G)*(5.62*log10(h/1.7/D50) + 6.25
   ! From van Rijn's formula, Cz0 = 18*log10(12*h/3/D90)
   !
   !------------------------------------------------------------------------------
   subroutine UpdateGroundRoughness(pft, Bag)
      implicit none
      integer, intent(in) :: pft(:)          ! vegetation type
      real(kind=8), intent(in) :: Bag(:)     ! aboveground biomass (kg/m2)
      real(kind=8), parameter :: Cb = 2.5d-3 ! bed drag coefficient (unitless)
      real(kind=8) :: cD0, ScD, alphaA, alphaD
      real(kind=8) :: betaA, betaD
      real(kind=8) :: asb, dsb, cD, h
      integer :: ii

      do ii = 1, NX, 1
         cD0 = par_cD0(pft(ii))
         ScD = par_ScD(pft(ii))
         alphaA = par_alphaA(pft(ii))
         betaA = par_betaA(pft(ii))
         alphaD = par_alphaD(pft(ii))
         betaD = par_betaD(pft(ii))
         h = m_uhydro(1,ii)
         asb = alphaA + Bag(ii)**betaA
         dsb = alphaD + Bag(ii)**betaD
         cD = cD0 + ScD * Bag(ii)
         m_Cz(ii) = par_Cz0*sqrt(2.0/(cD*asb*h+2.0*(1.0-asb*dsb)*Cb))
      end do
   end subroutine

   subroutine UpdateShearStress()
      implicit none
      real(kind=8) :: fcurr, fwave
      real(kind=8) :: h, U, Uwav, Hwav
      real(kind=8) :: tau_curr, tau_wave
      integer :: ii

      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         U = m_U(ii)
         Hwav = m_Hwav(ii)
         if (h>0) then
            ! bottom shear stress by currents
            fcurr = 0.24/(log(4.8*h/par_d50))**2
            tau_curr = 0.125*Roul*fcurr*U**2
            ! bottom shear stress by wave
            Uwav = PI*Hwav/force_T/sinh(Karman*h)
            fwave = 1.39*(6.0*Uwav*force_T/PI/par_d50)**(-0.52)
            tau_wave = 0.5*fwave*Roul*Uwav**2
            m_tau(ii) = tau_curr*(1.0+1.2*(tau_wave/(tau_curr+tau_wave))**3.2)
         else
            m_tau(ii) = 0.0d0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave number.
   !
   !------------------------------------------------------------------------------
   subroutine WaveNumberEQ(kwav, coefs, fval)
      implicit none
      real(kind=8), intent(in) :: kwav
      real(kind=8), intent(in) :: coefs(2)
      real(kind=8), intent(out) :: fval
      real(kind=8) :: sigma, T, h

      T = coefs(1) 
      h = coefs(2)
      sigma = 2.0*PI/T
      fval = sqrt(G*kwav*tanh(kwav*h)) - sigma
   end subroutine

   subroutine UpdateWaveNumber()
      implicit none
      real(kind=8) :: coefs(2), xbounds(2)
      real(kind=8) :: h, sigma
      integer :: ii
      
      sigma = 2.0*PI/force_T     ! wave frequency (dispersion)
      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h>0) then
            xbounds = (/sigma**2/G, sigma/sqrt(G*h)/)
            coefs = (/force_T, h/)
            call NonLRBrents(WaveNumberEQ, coefs, xbounds, 1d-4, m_kwav(ii))
         else
            m_kwav(ii) = 0.0d0
         end if
      end do
   end subroutine

   subroutine UpdateWaveNumber2()
      implicit none
      real(kind=8) :: coefs(2), xbounds(2)
      real(kind=8) :: h, h1, sigma
      integer :: ii

      sigma = 2.0*PI/force_T     ! wave frequency (dispersion)
      h1 = m_uhydro(1,1)
      if (h1>0) then
         xbounds = (/sigma**2/G, sigma/sqrt(G*h1)/)
         coefs = (/force_T, h1/)
         call NonLRBrents(WaveNumberEQ, coefs, xbounds, 1d-4, m_kwav(1))
         do ii = 2, NX, 1
            h = m_uhydro(1,ii)
            if (h>0) then
               m_kwav(ii) = m_kwav(1)*sqrt(h1/max(0.1,h))
            else
               m_kwav(ii) = 0.0d0
            end if
         end do
      else
         m_kwav = 0.0d0
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave depth breaking possibility.
   !
   !------------------------------------------------------------------------------
   subroutine BreakProbEQ(Qb, coefs, fval)
      implicit none
      real(kind=8), intent(in) :: Qb
      real(kind=8), intent(in) :: coefs(2)
      real(kind=8), intent(out) :: fval
      real(kind=8) :: Hrms, Hmax
      
      Hmax = coefs(1)
      Hrms = coefs(2)
      fval = (1-Qb)/log(Qb) + (Hrms/Hmax)**2
   end subroutine

   subroutine UpdateWaveBrkProb()
      implicit none
      real(kind=8) :: xbounds(2), coefs(2)
      real(kind=8) :: h, Hrms, Hmax
      integer :: ii

      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h>0) then
            Hmax = par_fr * h
            Hrms = m_Hwav(ii)
            xbounds = (/1d-10, 1.0-1d-10/)
            coefs = (/Hmax, Hrms/)
            call NonLRBrents(BreakProbEQ, coefs, xbounds, 1d-4, m_Qb(ii))
         else
            m_Qb(ii) = 1.0d0
         end if
      end do
   end subroutine

   subroutine UpdateWaveBrkProb2()
      implicit none
      real(kind=8) :: h, fHrms
      integer :: ii, nQb, idx 

      nQb = size(const_Qb)
      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h>0) then
            fHrms = m_Hwav(ii) / (par_fr*h) 
            if (fHrms<=const_Qb(1)) then
               m_Qb(ii) = 0.0d0
            else if (fHrms>=const_Qb(nQb)) then
               m_Qb(ii) = 1.0d0
            else
               call BinarySearch(const_Qb, fHrms, idx)
               m_Qb(ii) = 0.01*DBLE(idx) - 0.005
            end if
         else
            m_Qb(ii) = 1.0d0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave sources and sinks.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateWaveGeneration()
      implicit none
      real(kind=8), parameter :: Cd = 1.3d-3    ! drag coefficient for U10
      real(kind=8) :: sigma, alpha, beta
      real(kind=8) :: kwav, h
      integer :: ii

      sigma = 2.0*PI/force_T
      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h>0) then
            kwav = m_kwav(ii)
            alpha = 80.0*sigma*(Roua*Cd*force_U10/Roul/G/kwav)**2
            beta = 5.0*Roua/Roul/force_T*(force_U10*kwav/sigma-0.9)
            m_Swg(ii) = alpha + beta * m_Ewav(ii)
         else
            m_Swg(ii) = 0.0d0
         end if
      end do
   end subroutine

   subroutine UpdateWaveBtmFriction()
      implicit none
      real(kind=8) :: Hwav, kwav, Ewav
      real(kind=8) :: Cf, Qb, h
      integer :: ii

      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h>0) then
            Hwav = m_Hwav(ii)
            kwav = m_kwav(ii)
            Ewav = m_Ewav(ii)
            Qb = m_Qb(ii)
            Cf = 2.0*par_cbc*PI*Hwav/force_T/sinh(kwav*h)
            m_Sbf(ii) = (1-Qb)*2.0*Cf*kwav*Ewav/sinh(2.0*kwav*h)
         else
            m_Sbf(ii) = 0.0d0
         end if
      end do
   end subroutine
   
   subroutine UpdateWaveWhiteCapping()
      implicit none
      real(kind=8), parameter :: gammaPM = 4.57d-3
      real(kind=8), parameter :: m = 2.0d0
      real(kind=8), parameter :: cwc = 3.33d-5
      real(kind=8) :: sigma

      sigma = 2.0*PI/force_T
      m_Swc = cwc*sigma*((m_Ewav*(sigma**4)/G**2/gammaPM)**m)*m_Ewav
   end subroutine
   
   subroutine UpdateWaveDepthBrking()
      implicit none
      real(kind=8), parameter :: Cd = 1.3d-3    ! drag coefficient for U10
      real(kind=8) :: sigma, Qb, Hmax, alpha
      real(kind=8) :: Hwav, Ewav, kwav, h
      integer :: ii

      sigma = 2.0*PI/force_T
      do ii = 1, NX, 1
         h = m_uhydro(1,ii)
         if (h>0) then
            Qb = m_Qb(ii)
            kwav = m_kwav(ii)
            Hwav = m_Hwav(ii)
            Ewav = m_Ewav(ii)
            Hmax = par_fr * h
            alpha = 80.0*sigma*(Roua*Cd*force_U10/Roul/G/kwav)**2
            m_Sbrk(ii) = 2.0*alpha/force_T*Qb*((Hmax/Hwav)**2)*Ewav
         else
            m_Sbrk(ii) = 0.0d0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge convection flux. 
   !          Significant wave height is limited by fr times water depth.
   !
   !------------------------------------------------------------------------------
   subroutine CalcEdgeConvectionFlux(uhydro, fluxes)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(inout) :: fluxes(:,:)
      real(kind=8) :: sigma
      integer :: indx

      sigma = 2.0*PI/force_T
      tmp_U = uhydro(2,:) / max(0.1,uhydro(1,:))
      tmp_Nmax = 0.125*Roul*G*(par_fr*uhydro(1,:))**2/sigma
      indx = count(uhydro(1,:)>0)
      tmp_Cg(1:indx) = 0.5*sigma*(1.0+2.0*m_kwav(1:indx)*uhydro(1,1:indx)/ &
         sinh(2.0*m_kwav(1:indx)*uhydro(1,1:indx)))/m_kwav(1:indx)
      tmp_Cg(indx+1:NX) = 0.0d0
      fluxes(1,:) = uhydro(2,:)
      fluxes(2,:) = uhydro(1,:)*(tmp_U**2) + 0.5*G*uhydro(1,:)**2
      fluxes(3,:) = tmp_Cg*min(uhydro(3,:),tmp_Nmax)
      fluxes(4,:) = uhydro(2,:)*uhydro(4,:)
      fluxes(5,:) = uhydro(2,:)*uhydro(5,:)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge maximum gradient.
   !          Ignore the non-diagonal Jacobi terms for N, Css and Cj
   !          
   !          ZGEEV('N', 'N', NVAR, jacobi, NVAR, b, DUMMY, 1, DUMMY, 1, &
   !                WORK, 2*NVAR, WORK, err)
   !          complex(kind=8) :: jacobi(NVAR,NVAR), b(NVAR)
   !          complex(kind=8) :: DUMMY(1,1), WORK(2*NVAR)
   !          integer :: err
   !
   !------------------------------------------------------------------------------
   subroutine CalcEdgeMaxGradient(uhydro, gradient)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(inout) :: gradient(:)
      real(kind=8) :: sigma
      integer :: indx

      sigma = 2.0*PI/force_T
      indx = count(uhydro(1,:)>0)
      tmp_Cg(1:indx) = 0.5*sigma*(1.0+2.0*m_kwav(1:indx)*uhydro(1,1:indx)/ &
         sinh(2.0*m_kwav(1:indx)*uhydro(1,1:indx)))/m_kwav(1:indx)
      tmp_Cg(indx+1:NX) = 0.0d0
      tmp_U = uhydro(2,:) / max(0.1,uhydro(1,:))
      tmp_eigval(1,:) = 3.0*tmp_U + sqrt(5.0*(tmp_U**2)+G*uhydro(1,:))
      tmp_eigval(2,:) = 3.0*tmp_U - sqrt(5.0*(tmp_U**2)+G*uhydro(1,:))
      tmp_eigval(3,:) = tmp_Cg
      tmp_eigval(4,:) = uhydro(2,:)
      tmp_eigval(5,:) = uhydro(2,:)
      gradient = maxval(abs(tmp_eigval), dim=1)
   end subroutine

   subroutine CalcCellDiffusionFlux(uhydro, fluxes)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(inout) :: fluxes(:,:)

      fluxes = 0.0d0
      fluxes(4,2:NX-1) = 0.5*par_Kdf*(uhydro(1,2:NX-1)+uhydro(1,3:NX))* &
         (uhydro(4,3:NX)-uhydro(4,2:NX-1))/m_dX(2:NX-1)
      fluxes(5,2:NX-1) = 0.5*par_Kdf*(uhydro(1,2:NX-1)+uhydro(1,3:NX))* &
         (uhydro(5,3:NX)-uhydro(5,2:NX-1))/m_dX(2:NX-1)
   end subroutine

   subroutine CalcCellStateSources(uhydro, sources)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(inout) :: sources(:,:)
      real(kind=8) :: sigma

      sigma = 2.0*PI/force_T
      tmp_U = uhydro(2,:) / max(0.1,uhydro(1,:))
      sources(1,:) = 0.0d0
      sources(2,:) = -tmp_U*abs(tmp_U)*G/m_Cz**2 - G*uhydro(1,:)*m_dZh/m_dX
      sources(3,:) = (m_Swg - (m_Sbf+m_Swc+m_Sbrk)) / sigma 
      sources(4,:) = m_Esed - m_Dsed*uhydro(4,:)/(m_uhydro(4,:)+1d-3)
      sources(5,:) = 0.0d0
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Use the 1-D Finite Volume Semi-discrete KT central scheme to 
   !          discretize partial differential equation in the space domain.
   !
   !------------------------------------------------------------------------------
   subroutine TAIHydroEquations(uhydro, duhydro)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(inout) :: duhydro(:,:)
      real(kind=8) :: Fminus(NVAR), Fplus(NVAR)
      real(kind=8) :: ap, am, dx
      integer :: ii

      ! calculate slope limiter
      call FVSKT_Superbee(uhydro, tmp_phi)
      ! calculate cell edge variable values
      call FVSKT_celledge(uhydro, tmp_phi, tmp_uhydroL, tmp_uhydroR) 
      ! calculate cell edge convective fluxes 
      call CalcEdgeConvectionFlux(tmp_uhydroL, tmp_FL)
      call CalcEdgeConvectionFlux(tmp_uhydroR, tmp_FR)
      ! calculate cell edge maximum gradient
      call CalcEdgeMaxGradient(tmp_uhydroL, tmp_aL)
      call CalcEdgeMaxGradient(tmp_uhydroR, tmp_aR)
      ! calculate cell diffusion flux 
      call CalcCellDiffusionFlux(uhydro, tmp_P)
      ! calculate cell state sources
      call CalcCellStateSources(uhydro, tmp_SRC) 
      ! calculate temporal gradients
      do ii = 1, NX, 1
         dx = m_dX(ii)
         ap = max(0.0, tmp_aR(ii), tmp_aL(ii+1))
         am = max(0.0, tmp_aR(ii-1), tmp_aL(ii))
         if (ii==0) then
            ! seaward boundary condition
            duhydro(:,ii) = 0.0d0
         else if (ii==NX) then
            ! landward boundary condition
            Fminus = 0.5*(tmp_FR(:,ii-1)+tmp_FL(:,ii)) - &
               0.5*am*(tmp_uhydroL(:,ii)-tmp_uhydroR(:,ii-1))
            duhydro(:,ii) = Fminus/dx - tmp_P(:,ii-1)/dx + tmp_SRC(:,ii)
         else
            Fminus = 0.5*(tmp_FR(:,ii-1)+tmp_FL(:,ii)) - &
               0.5*am*(tmp_uhydroL(:,ii)-tmp_uhydroR(:,ii-1))
            Fplus = 0.5*(tmp_FR(:,ii)+tmp_FL(:,ii+1)) - &
               0.5*ap*(tmp_uhydroL(:,ii+1)-tmp_uhydroR(:,ii))
            duhydro(:,ii) = -(Fplus - Fminus) / dx + (tmp_P(:,ii) - &
               tmp_P(:,ii-1)) / dx + tmp_SRC(:,ii)
         end if
      end do
   end subroutine

end module TAIHydroMOD

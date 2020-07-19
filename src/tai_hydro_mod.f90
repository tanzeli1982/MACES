module tai_hydro_mod 
!---------------------------------------------------------------------------------
! Purpose:
!
! This module implements the 1-D transect-based hydrodynamic model
!
!---------------------------------------------------------------------------------
   use data_buffer_mod
   use hydro_utilities_mod 

   implicit none
   public
   !f2py real(kind=8), allocatable, dimension(:) :: sim_h
   !f2py real(kind=8), allocatable, dimension(:) :: sim_U
   !f2py real(kind=8), allocatable, dimension(:) :: sim_Hwav
   !f2py real(kind=8), allocatable, dimension(:) :: sim_Uwav
   !f2py real(kind=8), allocatable, dimension(:) :: sim_tau
   !f2py real(kind=8), allocatable, dimension(:) :: sim_Css
   real(kind=8), allocatable, dimension(:) :: sim_h
   real(kind=8), allocatable, dimension(:) :: sim_U
   real(kind=8), allocatable, dimension(:) :: sim_Hwav
   real(kind=8), allocatable, dimension(:) :: sim_Uwav
   real(kind=8), allocatable, dimension(:) :: sim_tau
   real(kind=8), allocatable, dimension(:) :: sim_Css

contains
   subroutine InitHydroMod(xin, zhin, fetchin, Cs0, nvar, npft, nx)
      implicit none
      !f2py real(kind=8), intent(in) :: xin, zhin, fetchin
      !f2py real(kind=8), intent(in) :: Cs0
      !f2py integer, intent(in) :: nvar, npft
      !f2py integer, intent(hide), depend(xin) :: nx = len(xin)
      real(kind=8), dimension(nx) :: xin     ! platform x coordinate (m) 
      real(kind=8), dimension(nx) :: zhin    ! platform surface elevation (msl)
      real(kind=8), dimension(nx) :: fetchin ! platform fetch length (m)
      real(kind=8) :: Cs0
      integer :: nvar               ! state variable number
      integer :: npft               ! pft number
      integer :: nx                 ! grid cell number
      ! local variables
      integer :: ii

      ! hydrodynamics state allocatable arrays
      allocate(m_uhydro(nx,nvar))      ; m_uhydro = 0.0d0
      allocate(m_X(nx))                ; m_X = xin
      allocate(m_dX(nx))               ; m_dX = 0.0d0
      allocate(m_Zh(nx))               ; m_Zh = zhin
      allocate(m_dZh(nx))              ; m_dZh = 0.0d0
      allocate(m_xfetch(nx))           ; m_xfetch = fetchin
      allocate(m_U(nx))                ; m_U = 0.0d0
      allocate(m_Hwav(nx))             ; m_Hwav = 0.0d0
      allocate(m_kwav(nx))             ; m_kwav = INFNT
      allocate(m_Ewav(nx))             ; m_Ewav = 0.0d0
      allocate(m_Uwav(nx))             ; m_Uwav = 0.0d0
      allocate(m_Twav(nx))             ; m_Twav = 2.0d0
      allocate(m_Qb(nx))               ; m_Qb = 0.0d0
      allocate(m_tau(nx))              ; m_tau = 0.0d0
      allocate(m_Swg(nx))              ; m_Swg = 0.0d0
      allocate(m_Sbf(nx))              ; m_Sbf = 0.0d0
      allocate(m_Swc(nx))              ; m_Swc = 0.0d0
      allocate(m_Sbrk(nx))             ; m_Sbrk = 0.0d0
      allocate(m_Cz(nx))               ; m_Cz = 0.0d0
      allocate(m_Cs(nx))               ; m_Cs = 0.0d0
      ! hydrodynamics temporary allocatable arrays
      allocate(tmp_uhydro(nx,nvar))    ; tmp_uhydro = 0.0d0
      allocate(tmp_uhydroL(nx,nvar))   ; tmp_uhydroL = 0.0d0
      allocate(tmp_uhydroR(nx,nvar))   ; tmp_uhydroR = 0.0d0
      allocate(tmp_phi(nx,nvar))       ; tmp_phi = 0.0d0
      allocate(tmp_FL(nx,nvar))        ; tmp_FL = 0.0d0
      allocate(tmp_FR(nx,nvar))        ; tmp_FR = 0.0d0
      allocate(tmp_P(nx,nvar))         ; tmp_P = 0.0d0
      allocate(tmp_SRC(nx,nvar))       ; tmp_SRC = 0.0d0
      allocate(tmp_eigval(nx,nvar))    ; tmp_eigval = 0.0d0
      allocate(tmp_aL(nx))             ; tmp_aL = 0.0d0
      allocate(tmp_aR(nx))             ; tmp_aR = 0.0d0
      allocate(tmp_U(nx))              ; tmp_U = 0.0d0
      allocate(tmp_B(nx))              ; tmp_B = 0.0d0
      ! output variables
      allocate(sim_h(nx))              ; sim_h = 0.0d0
      allocate(sim_U(nx))              ; sim_U = 0.0d0
      allocate(sim_Hwav(nx))           ; sim_Hwav = 0.0d0
      allocate(sim_Uwav(nx))           ; sim_Uwav = 0.0d0
      allocate(sim_tau(nx))            ; sim_tau = 0.0d0
      allocate(sim_Css(nx))            ; sim_Css = 0.0d0
      ! sources and sinks
      allocate(Cs_source(nx))          ; Cs_source = 0.0d0
      allocate(Cs_sink(nx))            ; Cs_sink = 0.0d0
      allocate(fctr_wave(nx))          ; fctr_wave = 1.0d0
      ! user-defined allocatable arrays
      allocate(rk4_K1(nx,nvar))        ; rk4_K1 = 0.0d0
      allocate(rk4_K2(nx,nvar))        ; rk4_K2 = 0.0d0
      allocate(rk4_K3(nx,nvar))        ; rk4_K3 = 0.0d0
      allocate(rk4_K4(nx,nvar))        ; rk4_K4 = 0.0d0
      allocate(rk4_K5(nx,nvar))        ; rk4_K5 = 0.0d0
      allocate(rk4_K6(nx,nvar))        ; rk4_K6 = 0.0d0
      allocate(rk4_nxt4th(nx,nvar))    ; rk4_nxt4th = 0.0d0
      allocate(rk4_nxt5th(nx,nvar))    ; rk4_nxt5th = 0.0d0
      allocate(rk4_interim(nx,nvar))   ; rk4_interim = 0.0d0
      allocate(rk4_rerr(nx,nvar))      ; rk4_rerr = 0.0d0
      allocate(par_cD0(npft))          ; par_cD0 = 0.0d0
      allocate(par_ScD(npft))          ; par_ScD = 0.0d0
      allocate(par_alphaA(npft))       ; par_alphaA = 0.0d0
      allocate(par_betaA(npft))        ; par_betaA = 0.0d0
      allocate(par_alphaD(npft))       ; par_alphaD = 0.0d0
      allocate(par_betaD(npft))        ; par_betaD = 0.0d0

      do ii = 1, nx, 1
         if (ii==1) then
            m_dX(ii) = 0.5 * (m_X(ii) + m_X(ii+1))
         else if (ii==nx) then
            m_dX(ii) = 0.5 * (m_X(ii-1) + m_X(ii))
         else
            m_dX(ii) = 0.5 * (m_X(ii+1) - m_X(ii-1))
         end if
         m_uhydro(ii,1) = max(-m_Zh(ii), 0.0)
      end do
      m_uhydro(:,3) = Cs0
   end subroutine

   subroutine FinalizeHydroMod()
      implicit none

      ! deallocate hydrodynamics state arrays
      deallocate(m_uhydro)
      deallocate(m_X)
      deallocate(m_dX)
      deallocate(m_Zh)
      deallocate(m_dZh)
      deallocate(m_U)
      deallocate(m_Hwav)
      deallocate(m_kwav)
      deallocate(m_Ewav)
      deallocate(m_Uwav)
      deallocate(m_Twav)
      deallocate(m_Qb)
      deallocate(m_tau)
      deallocate(m_Swg)
      deallocate(m_Sbf)
      deallocate(m_Swc)
      deallocate(m_Sbrk)
      deallocate(m_Cz)
      deallocate(m_Cs)
      ! deallocate hydrodynamics temporary arrays
      deallocate(tmp_uhydro)
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
      deallocate(tmp_B)
      deallocate(sim_h)
      deallocate(sim_U)
      deallocate(sim_Hwav)
      deallocate(sim_Uwav)
      deallocate(sim_tau)
      deallocate(sim_Css)
      deallocate(Cs_source)
      deallocate(Cs_sink)
      deallocate(fctr_wave)
      ! deallocate user-defined arrays
      deallocate(rk4_K1)
      deallocate(rk4_K2)
      deallocate(rk4_K3)
      deallocate(rk4_K4)
      deallocate(rk4_K5)
      deallocate(rk4_K6)
      deallocate(rk4_nxt4th)
      deallocate(rk4_nxt5th)
      deallocate(rk4_interim)
      deallocate(rk4_rerr)
      deallocate(par_cD0)
      deallocate(par_ScD)
      deallocate(par_alphaA)
      deallocate(par_betaA)
      deallocate(par_alphaD)
      deallocate(par_betaD)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Set model parameters.
   !
   !------------------------------------------------------------------------------
   subroutine SetModelParams(d50, Cz0, Kdf, cbc, cwc, fr, alphaA, betaA, &
                             alphaD, betaD, cD0, ScD, n)
      implicit none
      !f2py real(kind=8), intent(in) :: d50, Cz0, Kdf, cbc, cwc, fr
      !f2py real(kind=8), intent(in) :: alphaA, betaA, alphaD, betaD
      !f2py real(kind=8), intent(in) :: cD0, ScD
      !f2py integer, intent(hide), depend(alphaA) :: n = len(alphaA)
      real(kind=8) :: d50, Cz0, Kdf, cbc, cwc, fr
      real(kind=8), dimension(n) :: alphaA, betaA
      real(kind=8), dimension(n) :: alphaD, betaD
      real(kind=8), dimension(n) :: cD0, ScD
      integer :: n

      par_d50 = d50
      par_Cz0 = Cz0
      par_Kdf = Kdf
      par_cbc = cbc
      par_cwc = cwc
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
   subroutine ModelSetup(sources, sinks, zh, pft, Bag, xref, Twav, &
                         h0, U10, Cs0, n)
      implicit none
      !f2py real(kind=8), intent(in) :: sources, sinks
      !f2py real(kind=8), intent(in) :: zh, Bag
      !f2py integer, intent(in) :: pft
      !f2py real(kind=8), intent(in) :: xref
      !f2py real(kind=8), intent(in) :: Twav, h0, U10
      !f2py real(kind=8), intent(in) :: U10, Cs0
      !f2py integer, intent(hide), depend(zh) :: n = len(zh)
      real(kind=8), dimension(n) :: sources, sinks
      real(kind=8), dimension(n) :: zh, Bag
      integer, dimension(n) :: pft
      real(kind=8) :: xref
      real(kind=8) :: Twav, h0, U10
      real(kind=8) :: Cs0
      integer :: n
      ! local variables
      integer :: ii

      Cs_source = sources
      Cs_sink = sinks
      ! a typical wave period is 2s but increase greatly with the increase
      ! of wave speed (https://en.wikipedia.org/wiki/Wind_wave)
      frc_Twav = Twav
      frc_U10 = U10

      m_Zh = zh
      do ii = 1, n, 1
         if (ii==1) then
            m_dZh(ii) = 0.5 * (m_Zh(ii+1) - m_Zh(ii))
         else if (ii==n) then
            m_dZh(ii) = 0.5 * (m_Zh(ii) - m_Zh(ii-1))
         else
            m_dZh(ii) = 0.5 * (m_Zh(ii+1) - m_Zh(ii-1))
         end if
      end do
      call CalcWaveReductionByVeg(m_X, m_dX, Bag, xref, fctr_wave)
      call UpdateGroundRoughness(pft, Bag, m_uhydro(:,1), m_Cz)

      ! boundary conditions
      m_uhydro(1,1) = h0
      !m_uhydro(1,2) = h0*U0
      m_uhydro(1,3) = h0*Cs0
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Post-run model update.
   !          At depths smaller than 0.1 m, the flow velocity and sediment 
   !          transport are taken to zero (linearly interpolated to zero)??
   !
   !------------------------------------------------------------------------------
   subroutine ModelCallback(wave_mod)
      implicit none
      !f2py integer, intent(in) :: wave_mod
      integer :: wave_mod     ! wave mode
      ! local variables
      real(kind=8) :: h, kwav, Twav
      integer :: ii, n, m
      
      n = size(m_uhydro,1)
      m = size(m_uhydro,2)
      do ii = 1, n, 1
         if (m_uhydro(ii,1)<=TOL_REL) then
            m_uhydro(ii:n,1) = 0.0
            exit
         end if
      end do
      do ii = 3, m, 1
         where (m_uhydro(:,ii)<0) m_uhydro(:,ii) = 0.0
      end do

      do ii = 1, n, 1
         h = m_uhydro(ii,1)
         if (h<=TOL_REL) then
            m_uhydro(ii,2:m) = 0.0
            m_U(ii) = 0.0
            m_Cs(ii) = 0.0
         else
            m_U(ii) = m_uhydro(ii,2) / max(0.1,h)
            m_Cs(ii) = m_uhydro(ii,3) / max(0.1,h)
         end if 
      end do

      ! update wave dynamics
      if (wave_mod==EQM_WAVE) then
         do ii = 1, n, 1
            h = m_uhydro(ii,1)
            if (h<=TOL_REL) then
               m_Hwav(ii) = 0.0
               m_Uwav(ii) = 0.0
               m_kwav(ii) = INFNT
               m_Twav(ii) = frc_Twav
            else
               call UpdateSgnftWaveHeight(frc_U10, m_xfetch(ii), &
                                          h, m_Hwav(ii), Twav)
               call UpdateWaveNumber(Twav, h, kwav)
               m_Twav(ii) = Twav
               m_kwav(ii) = kwav
               m_Uwav(ii) = 2*PI*m_Hwav(ii)/Twav/sinh(kwav*max(0.1,h))
            end if
         end do
      else
         Twav = frc_Twav
         m_Twav = frc_Twav
         call UpdateWaveNumber(Twav, m_uhydro(:,1), m_kwav)
         call UpdateSgnftWaveHeight(Twav, frc_U10, m_uhydro(:,1), &
                                    m_kwav, m_Ewav)
         do ii = 1, n, 1
            h = m_uhydro(ii,1)
            kwav = m_kwav(ii)
            if (h<=TOL_REL) then
               m_Hwav(ii) = 0.0
               m_Uwav(ii) = 0.0
            else
               !m_Hwav(ii) = fctr_wave(ii) * sqrt(8.0*m_Ewav(ii)/G/Roul)
               m_Hwav(ii) = sqrt(8.0*m_Ewav(ii)/G/Roul)
               m_Uwav(ii) = 2*PI*m_Hwav(ii)/frc_Twav/sinh(kwav*max(0.1,h))
            end if
         end do
         !call UpdateWaveBrkProb(m_uhydro(:,1), m_Hwav, m_Qb)
         !call UpdateWaveBrkProb2(m_uhydro(:,1), m_Hwav, m_Qb)
         !call UpdateWaveGeneration(frc_Twav, frc_U10, m_uhydro(:,1), &
         !                          m_kwav, m_Ewav, m_Swg)
         !call UpdateWaveBtmFriction(frc_Twav, m_uhydro(:,1), m_Hwav, &
         !                           m_kwav, m_Ewav, m_Qb, m_Sbf)
         !call UpdateWaveWhiteCapping(frc_Twav, m_Ewav, m_Swc)
         !call UpdateWaveDepthBrking(frc_Twav, frc_U10, m_uhydro(:,1), &
         !                           m_Hwav, m_kwav, m_Ewav, m_Qb, m_Sbrk)
      end if
      call UpdateShearStress(m_Twav, m_uhydro(:,1), m_U, m_Uwav, m_tau)

      sim_h = m_uhydro(:,1)
      sim_U = m_U
      sim_Css = m_Cs
      sim_Hwav = m_Hwav
      sim_Uwav = m_Uwav
      sim_tau = m_tau
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Solve the TAI equations using the 4th-order Runge-Kutta-Fehlberg
   !          method.
   !
   !------------------------------------------------------------------------------
   subroutine ModelRun(mode, tol, dyncheck, curstep, ncurstep, &
                       nextstep, error, n)
      implicit none
      !f2py integer, intent(in) :: mode
      !f2py logical, intent(in) :: dyncheck
      !f2py real(kind=8), intent(in) :: tol, curstep
      !f2py real(kind=8), intent(out) :: ncurstep, nextstep
      !f2py integer, intent(out) :: error
      !f2py integer, intent(hide), depend(tol) :: n = len(tol)
      integer :: mode, error
      logical, dimension(n) :: dyncheck
      real(kind=8), dimension(n) :: tol
      real(kind=8) :: curstep, ncurstep, nextstep
      integer :: n

      ncurstep = curstep
      call RK4Fehlberg(TAIHydroEquations, m_uhydro, mode, tol, &
                       dyncheck, tmp_uhydro, ncurstep, nextstep, error)
      if (error==0) then
         m_uhydro = tmp_uhydro
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge convection flux. 
   !          Significant wave height is limited by fr times water depth.
   !
   !------------------------------------------------------------------------------
   subroutine CalcEdgeConvectionFlux(uhydro, fluxes, n, m)
      implicit none
      !f2py real(kind=8), intent(in) :: uhydro
      !f2py real(kind=8), intent(out) :: fluxes
      !f2py integer, intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py integer, intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, fluxes
      integer :: n, m
      ! local variables
      real(kind=8) :: U
      integer :: ii

      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            U = uhydro(ii,2) / max(0.1,uhydro(ii,1))
         else
            U = 0.0
         end if
         fluxes(ii,1) = uhydro(ii,2)
         fluxes(ii,2) = uhydro(ii,1)*(U**2)
         fluxes(ii,3:m) = U*uhydro(ii,3:m)
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge maximum gradient.
   !          Ignore the non-diagonal Jacobi terms for N, Css and Cj
   !          
   !          ZGEEV('N', 'N', nvar, jacobi, nvar, b, DUMMY, 1, DUMMY, 1, &
   !                WORK, 2*nvar, WORK, err)
   !          complex(kind=8) :: jacobi(nvar,nvar), b(nvar)
   !          complex(kind=8) :: DUMMY(1,1), WORK(2*nvar)
   !          integer :: err
   !
   !------------------------------------------------------------------------------
   subroutine CalcEdgeMaxGradient(uhydro, gradient, n, m)
      implicit none
      !f2py real(kind=8), intent(in) :: uhydro
      !f2py real(kind=8), intent(out) :: gradient
      !f2py integer, intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py integer, intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro
      real(kind=8), dimension(n) :: gradient
      integer :: n, m
      ! local variables
      real(kind=8) :: U
      integer :: ii

      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            U = uhydro(ii,2) / max(0.1,uhydro(ii,1))
         else
            U = 0.0
         end if
         tmp_eigval(ii,1) = 1.5*U + sqrt(1.25*(U**2))
         tmp_eigval(ii,2) = 1.5*U - sqrt(1.25*(U**2))      
         tmp_eigval(ii,3:m) = U
      end do
      gradient = maxval(abs(tmp_eigval), dim=2)
   end subroutine

   subroutine CalcCellDiffusionFlux(uhydro, fluxes, n, m)
      implicit none
      !f2py real(kind=8), intent(in) :: uhydro
      !f2py real(kind=8), intent(out) :: fluxes
      !f2py integer, intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py integer, intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, fluxes
      integer :: n, m
      ! local variables
      real(kind=8) :: Cs1, Cs2
      real(kind=8) :: h1, h2
      integer :: ii

      do ii = 1, n, 1
         fluxes(ii,1:2) = 0.0d0
         if (ii<n) then
            h1 = uhydro(ii,1)
            h2 = uhydro(ii+1,1)
            if (h1>TOL_REL) then
               Cs1 = uhydro(ii,3) / max(0.1,h1)
            else
               Cs1 = 0.0
            end if
            if (h2>TOL_REL) then
               Cs2 = uhydro(ii+1,3) / max(0.1,h2)
            else
               Cs2 = 0.0
            end if
            fluxes(ii,3) = 0.5*par_Kdf*(h1+h2)*(Cs2-Cs1)/m_dX(ii)
         else
            fluxes(ii,3) = 0.0d0
         end if
      end do
   end subroutine

   subroutine CalcCellStateSources(uhydro, sources, n, m)
      implicit none
      !f2py real(kind=8), intent(in) :: uhydro
      !f2py real(kind=8), intent(out) :: sources
      !f2py integer, intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py integer, intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, sources
      integer :: n, m
      ! local variables
      real(kind=8) :: scaler, scalerL, scalerR
      integer :: ii

      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            tmp_U(ii) = uhydro(ii,2) / max(0.1,uhydro(ii,1))
         else
            tmp_U(ii) = 0.0
         end if
         if (ii==1) then
            tmp_B(ii) = 0.5*(m_Zh(ii+1) + uhydro(ii+1,1) - &
               m_Zh(ii) - uhydro(ii,1))
         else if (ii==n) then
            tmp_B(ii) = 0.5*(m_Zh(ii) + uhydro(ii,1) - &
               m_Zh(ii-1) - uhydro(ii-1,1))
         else
            tmp_B(ii) = 0.5*(m_Zh(ii+1) + uhydro(ii+1,1) - &
               m_Zh(ii-1) - uhydro(ii-1,1))
         end if
      end do
      
      sources(:,1) = 0.0d0
      do ii = 1, n, 1
         scaler = max(0.0,uhydro(ii,3))/(m_uhydro(ii,3)+TOL_REL)
         if (ii==1) then
            scalerR = max(0.0,uhydro(ii+1,3))/(m_uhydro(ii+1,3)+TOL_REL)
            sources(ii,2) = -(0.75*tmp_U(ii)*abs(tmp_U(ii))*G*m_Cz(ii)+ &
               0.25*tmp_U(ii+1)*abs(tmp_U(ii+1))*G*m_Cz(ii+1)) - &
               G*(0.75*uhydro(ii,1)+0.25*uhydro(ii+1,1))*tmp_B(ii)/m_dX(ii)
            sources(ii,3) = (0.75*Cs_source(ii)+0.25*Cs_source(ii+1)) - &
               (0.75*Cs_sink(ii)*scaler+0.25*Cs_sink(ii+1)*scalerR)
         else if (ii==n) then
            scalerL = max(0.0,uhydro(ii-1,3))/(m_uhydro(ii-1,3)+TOL_REL)
            sources(ii,2) = -(0.75*tmp_U(ii)*abs(tmp_U(ii))*G*m_Cz(ii)+ &
               0.25*tmp_U(ii-1)*abs(tmp_U(ii-1))*G*m_Cz(ii-1)) - &
               G*(0.75*uhydro(ii,1)+0.25*uhydro(ii-1,1))*tmp_B(ii)/m_dX(ii)
            sources(ii,3) = (0.25*Cs_source(ii-1)+0.75*Cs_source(ii)) - &
               (0.25*Cs_sink(ii-1)*scalerL+0.75*Cs_sink(ii)*scaler)
         else
            scalerR = max(0.0,uhydro(ii+1,3))/(m_uhydro(ii+1,3)+TOL_REL)
            scalerL = max(0.0,uhydro(ii-1,3))/(m_uhydro(ii-1,3)+TOL_REL)
            sources(ii,2) = -(0.5*tmp_U(ii)*abs(tmp_U(ii))*G*m_Cz(ii)+ &
               0.25*tmp_U(ii-1)*abs(tmp_U(ii-1))*G*m_Cz(ii-1)+ &
               0.25*tmp_U(ii+1)*abs(tmp_U(ii+1))*G*m_Cz(ii+1)) - &
               G*(0.5*uhydro(ii,1)+0.25*uhydro(ii-1,1)+0.25*uhydro(ii+1,1))* &
               tmp_B(ii)/m_dX(ii)
            sources(ii,3) = (0.25*Cs_source(ii-1)+0.5*Cs_source(ii)+ &
               0.25*Cs_source(ii+1)) - (0.25*Cs_sink(ii-1)*scalerL+ &
               0.5*Cs_sink(ii)*scaler+0.25*Cs_sink(ii+1)*scalerR)
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Use the 1-D Finite Volume Semi-discrete KT central scheme to 
   !          discretize partial differential equation in the space domain.
   !
   !------------------------------------------------------------------------------
   subroutine TAIHydroEquations(uhydro, duhydro, n, m)
      implicit none
      !f2py real(kind=8), intent(in) :: uhydro
      !f2py real(kind=8), intent(out) :: duhydro
      !f2py integer, intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py integer, intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, duhydro
      integer :: n, m
      ! local variables
      real(kind=8) :: Fminus(m), Fplus(m)
      real(kind=8) :: ap, am, dx
      integer :: ii

      ! calculate slope limiter
      call FVSKT_Superbee(uhydro, tmp_phi, n, m)
      ! calculate cell edge variable values
      call FVSKT_celledge(uhydro, tmp_phi, tmp_uhydroL, tmp_uhydroR, n, m)
      ! calculate cell edge convective fluxes 
      call CalcEdgeConvectionFlux(tmp_uhydroL, tmp_FL, n, m)
      call CalcEdgeConvectionFlux(tmp_uhydroR, tmp_FR, n, m)
      ! calculate cell edge maximum gradient
      call CalcEdgeMaxGradient(tmp_uhydroL, tmp_aL, n, m)
      call CalcEdgeMaxGradient(tmp_uhydroR, tmp_aR, n, m)
      ! calculate cell diffusion flux 
      call CalcCellDiffusionFlux(uhydro, tmp_P, n, m)
      ! calculate cell state sources
      call CalcCellStateSources(uhydro, tmp_SRC, n, m) 
      ! calculate temporal gradients
      do ii = 1, n, 1
         dx = m_dX(ii)
         if (ii==1) then
            ! seaward boundary condition
            duhydro(ii,:) = 0.0
         else if (ii==n) then
            ! landward boundary condition
            am = max(0.0, tmp_aR(ii-1), tmp_aL(ii))
            Fminus = 0.5*(tmp_FR(ii-1,:)+tmp_FL(ii,:)) - &
               0.5*am*(tmp_uhydroL(ii,:)-tmp_uhydroR(ii-1,:))
            duhydro(ii,:) = Fminus/dx - tmp_P(ii-1,:)/dx + tmp_SRC(ii,:)
         else
            ! inner grid cells
            ap = max(0.0, tmp_aR(ii), tmp_aL(ii+1))
            am = max(0.0, tmp_aR(ii-1), tmp_aL(ii))
            Fminus = 0.5*(tmp_FR(ii-1,:)+tmp_FL(ii,:)) - &
               0.5*am*(tmp_uhydroL(ii,:)-tmp_uhydroR(ii-1,:))
            Fplus = 0.5*(tmp_FR(ii,:)+tmp_FL(ii+1,:)) - &
               0.5*ap*(tmp_uhydroL(ii+1,:)-tmp_uhydroR(ii,:))
            duhydro(ii,:) = -(Fplus - Fminus) / dx + (tmp_P(ii,:) - &
               tmp_P(ii-1,:)) / dx + tmp_SRC(ii,:)
         end if
      end do
      where (uhydro(:,3)<=0 .and. duhydro(:,3)<0) duhydro(:,3) = 0
   end subroutine

end module tai_hydro_mod 

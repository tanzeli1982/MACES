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

contains
   subroutine InitHydroMod(xin, zhin, nvar, npft, nx)
      implicit none
      !f2py real(kind=8), intent(in) :: xin, zhin
      !f2py integer, intent(in) :: nvar, npft
      !f2py integer, intent(hide), depend(xin) :: nx = len(xin)
      real(kind=8), dimension(nx) :: xin     ! platform x coordinate (m) 
      real(kind=8), dimension(nx) :: zhin    ! platform surface elevation (msl)
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
      allocate(m_U(nx))                ; m_U = 0.0d0
      allocate(m_Hwav(nx))             ; m_Hwav = 0.0d0
      allocate(m_kwav(nx))             ; m_kwav = 0.0d0
      allocate(m_Ewav(nx))             ; m_Ewav = 0.0d0
      allocate(m_Uwav(nx))             ; m_Uwav = 0.0d0
      allocate(m_tau(nx))              ; m_tau = 0.0d0
      allocate(m_Qb(nx))               ; m_Qb = 0.0d0
      allocate(m_Swg(nx))              ; m_Swg = 0.0d0
      allocate(m_Sbf(nx))              ; m_Sbf = 0.0d0
      allocate(m_Swc(nx))              ; m_Swc = 0.0d0
      allocate(m_Sbrk(nx))             ; m_Sbrk = 0.0d0
      allocate(m_Cz(nx))               ; m_Cz = 0.0d0
      allocate(m_Cg(nx))               ; m_Cg = 0.0d0
      allocate(m_Nmax(nx))             ; m_Nmax = 0.0d0
      allocate(m_Css(nx))              ; m_Css = 0.0d0
      allocate(m_Cj(nx))               ; m_Cj = 0.0d0
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
      allocate(tmp_Css(nx))            ; tmp_Css = 0.0d0
      allocate(tmp_Cj(nx))             ; tmp_Cj = 0.0d0
      allocate(tmp_Qb(101))            ; tmp_Qb = 0.0d0
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
      allocate(force_pft(nx))          ; force_pft = -1
      allocate(force_Bag(nx))          ; force_Bag = 0.0d0
      allocate(force_Esed(nx))         ; force_Esed = 0.0d0
      allocate(force_Dsed(nx))         ; force_Dsed = 0.0d0
      allocate(force_Fminus(nvar))     ; force_Fminus = 0.0d0

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
      tmp_Qb = (/0.0,0.4637,0.5005,0.5260,0.5461,0.5631,0.5780,0.5914, &
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
      deallocate(m_tau)
      deallocate(m_Qb)
      deallocate(m_Swg)
      deallocate(m_Sbf)
      deallocate(m_Swc)
      deallocate(m_Sbrk)
      deallocate(m_Cz)
      deallocate(m_Cg)
      deallocate(m_Nmax)
      deallocate(m_Css)
      deallocate(m_Cj)
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
      deallocate(tmp_Qb)
      deallocate(tmp_Css)
      deallocate(tmp_Cj)
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
      deallocate(force_Esed)
      deallocate(force_Dsed)
      deallocate(force_Bag)
      deallocate(force_pft)
      deallocate(force_Fminus)
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
   subroutine SetModelParams(d50, Cz0, Kdf, cbc, fr, alphaA, betaA, &
                             alphaD, betaD, cD0, ScD, n)
      implicit none
      !f2py real(kind=8), intent(in) :: d50, Cz0, Kdf, cbc, fr
      !f2py real(kind=8), intent(in) :: alphaA, betaA, alphaD, betaD
      !f2py real(kind=8), intent(in) :: cD0, ScD
      !f2py integer, intent(hide), depend(alphaA) :: n = len(alphaA)
      real(kind=8) :: d50, Cz0, Kdf, cbc, fr
      real(kind=8), dimension(n) :: alphaA, betaA
      real(kind=8), dimension(n) :: alphaD, betaD
      real(kind=8), dimension(n) :: cD0, ScD
      integer :: n

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
   subroutine ModelSetup(zh, pft, Bag, Esed, Dsed, Twav, U10, h0, & 
                         U0, Hwav0, Css0, Cj0, n)
      implicit none
      !f2py real(kind=8), intent(in) :: zh, Bag, Esed, Dsed
      !f2py integer, intent(in) :: pft
      !f2py real(kind=8), intent(in) :: Twav, U10, h0, U0, Hwav0
      !f2py real(kind=8), intent(in) :: Css0, Cj0
      !f2py integer, intent(hide), depend(zh) :: n = len(zh)
      real(kind=8), dimension(n) :: zh, Bag
      real(kind=8), dimension(n) :: Esed, Dsed
      integer, dimension(n) :: pft
      real(kind=8) :: Twav, U10, h0, U0, Hwav0
      real(kind=8) :: Css0, Cj0
      real(kind=8) :: sigma, Nwav
      integer :: n
      integer :: ii

      force_Esed = Esed
      force_Dsed = Dsed
      force_pft = pft
      force_Bag = Bag
      ! a typical wave period is 2s but increase greatly with the increase
      ! of wave speed (https://en.wikipedia.org/wiki/Wind_wave)
      force_Twav = Twav
      force_U10 = U10

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

      call UpdateGroundRoughness(pft, Bag, m_uhydro(:,1), m_Cz)
      !call UpdateWaveNumber(Twav, m_uhydro(:,1), m_kwav)
      call UpdateWaveNumber2(Twav, m_uhydro(:,1), m_kwav)
      !call UpdateWaveBrkProb(m_uhydro(:,1), m_Hwav, m_Qb)
      call UpdateWaveBrkProb2(m_uhydro(:,1), m_Hwav, tmp_Qb, m_Qb)
      call UpdateWaveGeneration(Twav, U10, m_uhydro(:,1), m_kwav, m_Ewav, m_Swg)
      call UpdateWaveBtmFriction(Twav, m_uhydro(:,1), m_Hwav, &
                                 m_kwav, m_Ewav, m_Qb, m_Sbf)
      call UpdateWaveWhiteCapping(Twav, m_Ewav, m_Swc)
      call UpdateWaveDepthBrking(Twav, U10, m_uhydro(:,1), m_Hwav, &
                                 m_kwav, m_Ewav, m_Qb, m_Sbrk)

      sigma = 2.0*PI/Twav
      do ii = 1, n, 1
         if (m_uhydro(ii,1)>TOL_REL .and. m_kwav(ii)>TOL_REL) then
            m_Cg(ii) = 0.5*sigma*(1.0+2.0*m_kwav(ii)*m_uhydro(ii,1)/ &
               sinh(2.0*m_kwav(ii)*m_uhydro(ii,1)))/m_kwav(ii)
         else
            m_Cg(ii) = 0.0
         end if
         if (m_uhydro(ii,1)>TOL_REL) then
            m_Nmax(ii) = 0.125*Roul*G*(par_fr*m_uhydro(ii,1))**2/sigma
         else
            m_Nmax(ii) = 0.0
         end if
      end do
      Nwav = 0.125*Roul*G*(Hwav0**2)/sigma
      ! boundary conditions
      force_Fminus(1) = h0 * U0
      force_Fminus(2) = h0*U0**2 + 0.5*G*h0**2
      force_Fminus(3) = m_Cg(1)*Nwav
      force_Fminus(4) = U0*h0*Css0
      force_Fminus(5) = U0*h0*Cj0
      m_uhydro(1,1) = h0
      m_uhydro(1,2) = h0*U0
      m_uhydro(1,3) = Nwav
      m_uhydro(1,4) = h0*Css0
      m_uhydro(1,5) = h0*Cj0
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
      ! local variables
      real(kind=8) :: sigma, h
      integer :: ii, n, m
      
      n = size(m_uhydro,1)
      m = size(m_uhydro,2)
      sigma = 2.0*PI/force_Twav
      do ii = 1, n, 1
         if (m_uhydro(ii,1)<=TOL_REL) then
            m_uhydro(ii:n,1) = 0.0
            exit
         end if
      end do
      !where (m_uhydro(:,1)<0) m_uhydro(:,1) = 0.0
      where (m_uhydro(:,3)<0) m_uhydro(:,3) = 0.0
      where (m_uhydro(:,4)<0) m_uhydro(:,4) = 0.0
      where (m_uhydro(:,5)<0) m_uhydro(:,5) = 0.0
      do ii = 1, n, 1
         h = m_uhydro(ii,1)
         if (h<=TOL_REL) then
            m_uhydro(ii,2:m) = 0.0
            m_U(ii) = 0.0
            m_Css(ii) = 0.0
            m_Cj(ii) = 0.0
            m_Ewav(ii) = 0.0
            m_Hwav(ii) = 0.0
            m_Uwav(ii) = 0.0
         else
            m_U(ii) = m_uhydro(ii,2) / max(0.1,h)
            m_Css(ii) = m_uhydro(ii,4) / max(0.1,h)
            m_Cj(ii) = m_uhydro(ii,5) / max(0.1,h)
            m_Ewav(ii) = sigma * m_uhydro(ii,3)
            m_Hwav(ii) = sqrt(8.0*m_Ewav(ii)/G/Roul)
            m_Uwav(ii) = min(PI*m_Hwav(ii)/force_Twav/sinh(Karman*h), 20.0)
         end if 
      end do
      call UpdateShearStress(force_Twav, m_uhydro(:,1), m_U, m_Hwav, &
                             m_Uwav, m_tau)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Solve the TAI equations using the 4th-order Runge-Kutta-Fehlberg
   !          method.
   !
   !------------------------------------------------------------------------------
   subroutine ModelRun(mode, tol, dyncheck, curstep, ncurstep, nextstep, error, n)
      implicit none
      !f2py integer, intent(in) :: mode, dyncheck
      !f2py real(kind=8), intent(in) :: tol, curstep
      !f2py real(kind=8), intent(out) :: ncurstep, nextstep
      !f2py integer, intent(out) :: error
      !f2py integer, intent(hide), depend(tol) :: n = len(tol)
      integer :: mode, error
      integer, dimension(n) :: dyncheck
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
   ! Purpose: Return model simulations
   !
   !------------------------------------------------------------------------------
   subroutine GetModelSims(n, h, U, Hwav, tau, Css, Cj)
      implicit none
      !f2py intent(in) :: n
      !f2py intent(out) :: h, U, Hwav, tau, Css, Cj
      real(kind=8), dimension(n) :: h, U, Hwav
      real(kind=8), dimension(n) :: tau, Css, Cj
      integer :: n

      h = m_uhydro(:,1)
      U = m_U
      Hwav = m_Hwav
      tau = m_tau
      Css = m_Css
      Cj = m_Cj
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
      real(kind=8) :: sigma, U, B
      integer :: ii

      sigma = 2.0*PI/force_Twav
      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            U = uhydro(ii,2) / max(0.1,uhydro(ii,1))
            B = 0.5*G*uhydro(ii,1)**2
         else
            U = 0.0
            B = 0.0
         end if
         fluxes(ii,1) = uhydro(ii,2)
         fluxes(ii,2) = uhydro(ii,1)*(U**2) + B
         fluxes(ii,3) = m_Cg(ii)*min(uhydro(ii,3),m_Nmax(ii))
         fluxes(ii,4) = U*uhydro(ii,4)
         fluxes(ii,5) = U*uhydro(ii,5)
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
      real(kind=8) :: sigma, U, B
      integer :: ii

      sigma = 2.0*PI/force_Twav
      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            U = uhydro(ii,2) / max(0.1,uhydro(ii,1))
            B = G*uhydro(ii,1)
         else
            U = 0.0
            B = 0.0
         end if
         tmp_eigval(ii,1) = 1.5*U + sqrt(1.25*(U**2)+B)
         tmp_eigval(ii,2) = 1.5*U - sqrt(1.25*(U**2)+B)      
         tmp_eigval(ii,3) = m_Cg(ii)
         tmp_eigval(ii,4) = U
         tmp_eigval(ii,5) = U
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
      integer :: ii

      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            tmp_Css(ii) = uhydro(ii,4) / max(0.1,uhydro(ii,1))
            tmp_Cj(ii) = uhydro(ii,5) / max(0.1,uhydro(ii,1))
         else
            tmp_Css(ii) = 0.0
            tmp_Cj(ii) = 0.0
         end if
      end do
      fluxes = 0.0d0
      fluxes(2:n-1,4) = 0.5*par_Kdf*(uhydro(2:n-1,1)+uhydro(3:n,1))* &
         (tmp_Css(3:n)-tmp_Css(2:n-1))/m_dX(2:n-1)
      fluxes(2:n-1,5) = 0.5*par_Kdf*(uhydro(2:n-1,1)+uhydro(3:n,1))* &
         (tmp_Cj(3:n)-tmp_Cj(2:n-1))/m_dX(2:n-1)
   end subroutine

   subroutine CalcCellStateSources(uhydro, sources, n, m)
      implicit none
      !f2py real(kind=8), intent(in) :: uhydro
      !f2py real(kind=8), intent(out) :: sources
      !f2py integer, intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py integer, intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, sources
      integer :: n, m
      real(kind=8) :: sigma, scaler
      integer :: ii

      sigma = 2.0*PI/force_Twav
      do ii = 1, n, 1
         if (uhydro(ii,1)>TOL_REL) then
            tmp_U(ii) = uhydro(ii,2) / max(0.1,uhydro(ii,1))
            !tmp_B(ii) = m_dZh(ii)
         else
            tmp_U(ii) = 0.0
            !tmp_B(ii) = 0.0
         end if
      end do
      
      sources(:,1) = 0.0d0
      sources(:,5) = 0.0d0
      do ii = 1, n, 1
         scaler = max(0.0,uhydro(ii,4))/(m_uhydro(ii,4)+TOL_REL)
         if (ii==1) then
            sources(ii,2) = -(0.75*tmp_U(ii)*abs(tmp_U(ii))*G/m_Cz(ii)**2+ &
               0.25*tmp_U(ii+1)*abs(tmp_U(ii+1))*G/m_Cz(ii+1)**2) - &
               G*(0.75*uhydro(ii,1)+0.25*uhydro(ii+1,1))*m_dZh(ii)/m_dX(ii)
            sources(ii,3) = (0.75*(m_Swg(ii)-m_Sbf(ii)-m_Swc(ii)-m_Sbrk(ii)) + &
               0.25*(m_Swg(ii+1)-m_Sbf(ii+1)-m_Swc(ii+1)-m_Sbrk(ii+1))) / sigma
            sources(ii,4) = (0.75*force_Esed(ii)+0.25*force_Esed(ii+1)) - &
               (0.75*force_Dsed(ii)+0.25*force_Dsed(ii+1))*scaler
         else if (ii==n) then
            sources(ii,2) = -(0.75*tmp_U(ii)*abs(tmp_U(ii))*G/m_Cz(ii)**2+ &
               0.25*tmp_U(ii-1)*abs(tmp_U(ii-1))*G/m_Cz(ii-1)**2) - &
               G*(0.75*uhydro(ii,1)+0.25*uhydro(ii-1,1))*m_dZh(ii)/m_dX(ii)
            sources(ii,3) = (0.75*(m_Swg(ii)-m_Sbf(ii)-m_Swc(ii)-m_Sbrk(ii)) + &
               0.25*(m_Swg(ii-1)-m_Sbf(ii-1)-m_Swc(ii-1)-m_Sbrk(ii-1))) / sigma
            sources(ii,4) = (0.25*force_Esed(ii-1)+0.75*force_Esed(ii)) - &
               (0.25*force_Dsed(ii-1)+0.75*force_Dsed(ii))*scaler
         else
            sources(ii,2) = -(0.5*tmp_U(ii)*abs(tmp_U(ii))*G/m_Cz(ii)**2+ &
               0.25*tmp_U(ii+1)*abs(tmp_U(ii+1))*G/m_Cz(ii+1)**2+ &
               0.25*tmp_U(ii-1)*abs(tmp_U(ii-1))*G/m_Cz(ii-1)**2) - &
               G*(0.5*uhydro(ii,1)+0.25*uhydro(ii-1,1)+0.25*uhydro(ii+1,1))* &
               m_dZh(ii)/m_dX(ii)
            sources(ii,3) = (0.5*(m_Swg(ii)-m_Sbf(ii)-m_Swc(ii)-m_Sbrk(ii)) + &
               0.25*(m_Swg(ii+1)-m_Sbf(ii+1)-m_Swc(ii+1)-m_Sbrk(ii+1)) + &
               0.25*(m_Swg(ii-1)-m_Sbf(ii-1)-m_Swc(ii-1)-m_Sbrk(ii-1))) / sigma
            sources(ii,4) = (0.25*force_Esed(ii-1)+0.5*force_Esed(ii)+ &
               0.25*force_Esed(ii+1)) - (0.25*force_Dsed(ii-1)+ &
               0.5*force_Dsed(ii)+0.25*force_Dsed(ii+1))*scaler
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
            !ap = max(0.0, tmp_aR(ii), tmp_aL(ii+1))
            !Fplus = 0.5*(tmp_FR(ii,:)+tmp_FL(ii+1,:)) - &
            !   0.5*ap*(tmp_uhydroL(ii+1,:)-tmp_uhydroR(ii,:))
            duhydro(ii,:) = 0.0 !-(Fplus - force_Fminus) / dx + tmp_SRC(ii,:)
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
   end subroutine

end module tai_hydro_mod 

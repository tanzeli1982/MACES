module TAIHydroMOD
!---------------------------------------------------------------------------------
! Purpose:
!
! This module implements the 1-D transect-based hydrodynamic model
!
!---------------------------------------------------------------------------------
   use data_type_mod
   use data_buffer_mod
   use hydro_utilities_mod 

   implicit none

contains
   subroutine InitHydroMod(xin, zhin, nvar, npft, nx)
      implicit none
      !f2py intent(in) :: xin, zhin, nvar, npft
      !f2py intent(hide), depend(xin) :: nx = len(xin)
      real(kind=8), dimension(nx) :: xin     ! platform x coordinate (m) 
      real(kind=8), dimension(nx) :: zhin    ! platform surface elevation (msl)
      integer :: nvar               ! state variable number
      integer :: npft               ! pft number
      integer :: nx                 ! grid cell number
      ! local variables
      integer :: ii

      ! RungeKutta4 allocatable arrays
      allocate(rk4_K1(nvar,nx))        ; rk4_K1 = 0.0d0
      allocate(rk4_K2(nvar,nx))        ; rk4_K2 = 0.0d0
      allocate(rk4_K3(nvar,nx))        ; rk4_K3 = 0.0d0
      allocate(rk4_K4(nvar,nx))        ; rk4_K4 = 0.0d0
      allocate(rk4_K5(nvar,nx))        ; rk4_K5 = 0.0d0
      allocate(rk4_K6(nvar,nx))        ; rk4_K6 = 0.0d0
      allocate(rk4_nxt4th(nvar,nx))    ; rk4_nxt4th = 0.0d0
      allocate(rk4_nxt5th(nvar,nx))    ; rk4_nxt5th = 0.0d0
      allocate(rk4_interim(nvar,nx))   ; rk4_interim = 0.0d0
      allocate(rk4_rerr(nvar,nx))      ; rk4_rerr = 0.0d0
      ! hydrodynamics state allocatable arrays
      allocate(m_uhydro(nvar,nx))      ; m_uhydro = 0.0d0
      allocate(m_X(nx))                ; m_X = xin
      allocate(m_dX(nx))               ; m_dX = 0.0d0
      allocate(m_Zh(nx))               ; m_Zh = zhin
      allocate(m_dZh(nx))              ; m_dZh = 0.0d0
      allocate(m_U(nx))                ; m_U = 0.0d0
      allocate(m_Hwav(nx))             ; m_Hwav = 0.0d0
      allocate(m_Ewav(nx))             ; m_Ewav = 0.0d0
      allocate(m_Uwav(nx))             ; m_Uwav = 0.0d0
      allocate(m_tau(nx))              ; m_tau = 0.0d0
      allocate(m_Qb(nx))               ; m_Qb = 0.0d0
      allocate(m_Swg(nx))              ; m_Swg = 0.0d0
      allocate(m_Sbf(nx))              ; m_Sbf = 0.0d0
      allocate(m_Swc(nx))              ; m_Swc = 0.0d0
      allocate(m_Sbrk(nx))             ; m_Sbrk = 0.0d0
      allocate(m_Cz(nx))               ; m_Cz = 0.0d0
      ! hydrodynamics temporary allocatable arrays
      allocate(tmp_uhydroL(nvar,nx))   ; tmp_uhydroL = 0.0d0
      allocate(tmp_uhydroR(nvar,nx))   ; tmp_uhydroR = 0.0d0
      allocate(tmp_phi(nvar,nx))       ; tmp_phi = 0.0d0
      allocate(tmp_FL(nvar,nx))        ; tmp_FL = 0.0d0
      allocate(tmp_FR(nvar,nx))        ; tmp_FR = 0.0d0
      allocate(tmp_P(nvar,nx))         ; tmp_P = 0.0d0
      allocate(tmp_SRC(nvar,nx))       ; tmp_SRC = 0.0d0
      allocate(tmp_eigval(nvar,nx))    ; tmp_eigval = 0.0d0
      allocate(tmp_aL(nx))             ; tmp_aL = 0.0d0
      allocate(tmp_aR(nx))             ; tmp_aR = 0.0d0
      allocate(tmp_U(nx))              ; tmp_U = 0.0d0
      allocate(tmp_Cg(nx))             ; tmp_Cg = 0.0d0
      allocate(tmp_Nmax(nx))           ; tmp_Nmax = 0.0d0
      allocate(tmp_Qb(101))            ; tmp_Qb = 0.0d0
      ! user-defined allocatable arrays
      allocate(m_params%cD0(npft))     ; m_params%cD0 = 0.0d0
      allocate(m_params%ScD(npft))     ; m_params%ScD = 0.0d0
      allocate(m_params%alphaA(npft))  ; m_params%alphaA = 0.0d0
      allocate(m_params%betaA(npft))   ; m_params%betaA = 0.0d0
      allocate(m_params%alphaD(npft))  ; m_params%alphaD = 0.0d0
      allocate(m_params%betaD(npft))   ; m_params%betaD = 0.0d0
      allocate(m_forcings%pft(nx))     ; m_forcings%pft = -1
      allocate(m_forcings%Bag(nx))     ; m_forcings%Bag = 0.0d0
      allocate(m_forcings%Esed(nx))    ; m_forcings%Esed = 0.0d0
      allocate(m_forcings%Dsed(nx))    ; m_forcings%Dsed = 0.0d0
      
      m_uhydro = 0.0
      do ii = 1, nx, 1
         if (ii==1) then
            m_dX(ii) = 0.5 * (m_X(ii) + m_X(ii+1))
         else if (ii==nx) then
            m_dX(ii) = 0.5 * (m_X(ii-1) + m_X(ii))
         else
            m_dX(ii) = 0.5 * (m_X(ii+1) - m_X(ii-1))
         end if
         m_uhydro(1,ii) = max(-m_Zh(ii), 0.0)
         if (m_uhydro(1,ii)>0) then
            m_uhydro(5,ii) = 35.0d0
         end if
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

      ! deallocate rungekutta4 arrays
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
      ! deallocate hydrodynamics state arrays
      deallocate(m_uhydro)
      deallocate(m_X)
      deallocate(m_dX)
      deallocate(m_Zh)
      deallocate(m_dZh)
      deallocate(m_U)
      deallocate(m_Hwav)
      deallocate(m_Ewav)
      deallocate(m_Uwav)
      deallocate(m_tau)
      deallocate(m_Qb)
      deallocate(m_Swg)
      deallocate(m_Sbf)
      deallocate(m_Swc)
      deallocate(m_Sbrk)
      deallocate(m_Cz)
      ! deallocate hydrodynamics temporary arrays
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
      deallocate(tmp_Qb)
      ! deallocate user-defined arrays
      deallocate(m_params%cD0)
      deallocate(m_params%ScD)
      deallocate(m_params%alphaA)
      deallocate(m_params%betaA)
      deallocate(m_params%alphaD)
      deallocate(m_params%betaD)
      deallocate(m_forcings%Esed)
      deallocate(m_forcings%Dsed)
      deallocate(m_forcings%Bag)
      deallocate(m_forcings%pft)
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
      !f2py intent(in) :: d50, Cz0, Kdf, cbc, fr
      !f2py intent(in) :: alphaA, betaA, alphaD, betaD
      !f2py intent(in) :: cD0, ScD
      !f2py intent(hide), depend(alphaA) :: n = len(alphaA)
      real(kind=8) :: d50, Cz0, Kdf, cbc, fr
      real(kind=8), dimension(n) :: alphaA, betaA
      real(kind=8), dimension(n) :: alphaD, betaD
      real(kind=8), dimension(n) :: cD0, ScD
      integer :: n

      m_params%d50 = d50
      m_params%Cz0 = Cz0
      m_params%Kdf = Kdf
      m_params%cbc = cbc
      m_params%fr = fr
      m_params%alphaA = alphaA
      m_params%betaA = betaA
      m_params%alphaD = alphaD
      m_params%betaD = betaD
      m_params%cD0 = cD0
      m_params%ScD = ScD
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
      !f2py intent(in) :: zh, pft, Bag, Esed, Dsed
      !f2py intent(in) :: Twav, U10, h0, U0, Hwav0
      !f2py intent(in) :: Css0, Cj0
      !f2py intent(hide), depend(Zh) :: n = len(Zh)
      real(kind=8), dimension(n) :: zh, Bag
      real(kind=8), dimension(n) :: Esed, Dsed
      integer, dimension(n) :: pft
      real(kind=8) :: Twav, U10, h0, U0, Hwav0
      real(kind=8) :: Css0, Cj0
      integer :: n
      integer :: ii

      m_forcings%Esed = Esed
      m_forcings%Dsed = Dsed
      m_forcings%pft = pft
      m_forcings%Bag = Bag
      ! a typical wave period is 2s but increase greatly with the increase
      ! of wave speed (https://en.wikipedia.org/wiki/Wind_wave)
      m_forcings%Twav = Twav
      m_forcings%U10 = U10
      m_forcings%h0 = h0
      m_forcings%U0 = U0
      m_forcings%Hwav0 = Hwav0
      m_forcings%Css0 = Css0
      m_forcings%Cj0 = Cj0

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

      call UpdateGroundRoughness(m_forcings, m_uhydro(1,:), & 
                                 m_params, m_Cz)
      !call UpdateWaveNumber(m_forcings, m_uhydro(1,:), m_kwav)
      call UpdateWaveNumber2(m_forcings, m_uhydro(1,:), m_kwav)
      !call UpdateWaveBrkProb(m_uhydro(1,:), m_Hwav, m_params, m_Qb)
      call UpdateWaveBrkProb2(m_uhydro(1,:), m_Hwav, m_params, &
                              tmp_Qb, m_Qb)
      call UpdateWaveGeneration(m_forcings, m_uhydro(1,:), m_kwav, &
                                m_Ewav, m_Swg)
      call UpdateWaveBtmFriction(m_forcings, m_uhydro(1,:), m_Hwav, &
                                 m_kwav, m_Ewav, m_Qb, m_params, m_Sbf)
      call UpdateWaveWhiteCapping(m_forcings, m_Ewav, m_Swc)
      call UpdateWaveDepthBrking(m_forcings, m_uhydro(1,:), m_Hwav, &
                                 m_kwav, m_Ewav, m_Qb, m_params, m_Sbrk)

      ! boundary conditions
      m_uhydro(1,1) = h0
      m_uhydro(2,1) = U0 * h0
      m_uhydro(3,1) = 0.125*Roul*G*(Hwav0**2)/(2.0*PI/Twav)
      m_uhydro(4,1) = Css0
      m_uhydro(5,1) = Cj0
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Post-run model update.
   !          At depths smaller than 0.1 m, the flow velocity and sediment 
   !          transport are taken to zero (linearly interpolated to zero)??
   !
   !------------------------------------------------------------------------------
   subroutine ModelCallback(uhydro, n, m)
      implicit none
      !f2py intent(in) :: uhydro
      !f2py intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro
      integer :: n, m
      ! local variables
      real(kind=8) :: sigma, Hwav, h
      integer :: ii
      
      m_uhydro = uhydro
      sigma = 2.0*PI/m_forcings%Twav
      do ii = 1, m, 1
         if (m_uhydro(1,ii)<=0) then
            m_uhydro(2:n,ii) = 0.0d0
         end if
      end do
      m_U = m_uhydro(2,:) / max(0.1,m_uhydro(1,:))
      m_Ewav = sigma * m_uhydro(3,:)
      m_Hwav = sqrt(8.0*m_Ewav/G/Roul)
      do ii = 1, m, 1
         h = m_uhydro(1,ii)
         Hwav = m_Hwav(ii)
         if (h<=0) then
            m_Uwav(ii) = 0.0d0
         else
            m_Uwav(ii) = min(PI*Hwav/m_forcings%Twav/sinh(Karman*h), 20.0)
         end if
      end do
      call UpdateShearStress(m_forcings, m_uhydro(1,:), m_U, m_Hwav, &
                             m_Uwav, m_params, m_tau)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge convection flux. 
   !          Significant wave height is limited by fr times water depth.
   !
   !------------------------------------------------------------------------------
   subroutine CalcEdgeConvectionFlux(uhydro, fluxes, n, m)
      implicit none
      !f2py intent(in) :: uhydro
      !f2py intent(out) :: fluxes
      !f2py intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, fluxes
      integer :: n, m
      ! local variables
      real(kind=8) :: sigma
      integer :: indx

      sigma = 2.0*PI/m_forcings%Twav
      tmp_U = uhydro(2,:) / max(0.1,uhydro(1,:))
      tmp_Nmax = 0.125*Roul*G*(m_params%fr*uhydro(1,:))**2/sigma
      indx = count(uhydro(1,:)>0)
      tmp_Cg(1:indx) = 0.5*sigma*(1.0+2.0*m_kwav(1:indx)*uhydro(1,1:indx)/ &
         sinh(2.0*m_kwav(1:indx)*uhydro(1,1:indx)))/m_kwav(1:indx)
      tmp_Cg(indx+1:m) = 0.0d0
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
   !          ZGEEV('N', 'N', nvar, jacobi, nvar, b, DUMMY, 1, DUMMY, 1, &
   !                WORK, 2*nvar, WORK, err)
   !          complex(kind=8) :: jacobi(nvar,nvar), b(nvar)
   !          complex(kind=8) :: DUMMY(1,1), WORK(2*nvar)
   !          integer :: err
   !
   !------------------------------------------------------------------------------
   subroutine CalcEdgeMaxGradient(uhydro, gradient, n, m)
      implicit none
      !f2py intent(in) :: uhydro
      !f2py intent(out) :: gradient
      !f2py intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro
      real(kind=8), dimension(m) :: gradient
      integer :: n, m
      ! local variables
      real(kind=8) :: sigma
      integer :: indx

      sigma = 2.0*PI/m_forcings%Twav
      indx = count(uhydro(1,:)>0)
      tmp_Cg(1:indx) = 0.5*sigma*(1.0+2.0*m_kwav(1:indx)*uhydro(1,1:indx)/ &
         sinh(2.0*m_kwav(1:indx)*uhydro(1,1:indx)))/m_kwav(1:indx)
      tmp_Cg(indx+1:m) = 0.0d0
      tmp_U = uhydro(2,:) / max(0.1,uhydro(1,:))
      tmp_eigval(1,:) = 3.0*tmp_U + sqrt(5.0*(tmp_U**2)+G*uhydro(1,:))
      tmp_eigval(2,:) = 3.0*tmp_U - sqrt(5.0*(tmp_U**2)+G*uhydro(1,:))
      tmp_eigval(3,:) = tmp_Cg
      tmp_eigval(4,:) = uhydro(2,:)
      tmp_eigval(5,:) = uhydro(2,:)
      gradient = maxval(abs(tmp_eigval), dim=1)
   end subroutine

   subroutine CalcCellDiffusionFlux(uhydro, fluxes, n, m)
      implicit none
      !f2py intent(in) :: uhydro
      !f2py intent(out) :: fluxes
      !f2py intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, fluxes
      integer :: n, m

      fluxes = 0.0d0
      fluxes(4,2:m-1) = 0.5*m_params%Kdf*(uhydro(1,2:m-1)+uhydro(1,3:m))* &
         (uhydro(4,3:m)-uhydro(4,2:m-1))/m_dX(2:m-1)
      fluxes(5,2:m-1) = 0.5*m_params%Kdf*(uhydro(1,2:m-1)+uhydro(1,3:m))* &
         (uhydro(5,3:m)-uhydro(5,2:m-1))/m_dX(2:m-1)
   end subroutine

   subroutine CalcCellStateSources(uhydro, sources, n, m)
      implicit none
      !f2py intent(in) :: uhydro
      !f2py intent(out) :: sources
      !f2py intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, sources
      integer :: n, m
      real(kind=8) :: sigma

      sigma = 2.0*PI/m_forcings%Twav
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
   subroutine TAIHydroEquations(uhydro, duhydro, n, m)
      implicit none
      !f2py intent(in) :: uhydro
      !f2py intent(out) :: duhydro
      !f2py intent(hide), depend(uhydro) :: n = shape(uhydro,0)
      !f2py intent(hide), depend(uhydro) :: m = shape(uhydro,1)
      real(kind=8), dimension(n,m) :: uhydro, duhydro
      integer :: n, m
      ! local variables
      real(kind=8) :: Fminus(n), Fplus(n)
      real(kind=8) :: ap, am, dx
      integer :: ii

      ! calculate slope limiter
      call FVSKT_Superbee(uhydro, tmp_phi)
      ! calculate cell edge variable values
      call FVSKT_celledge(uhydro, tmp_phi, tmp_uhydroL, tmp_uhydroR) 
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
      do ii = 1, m, 1
         dx = m_dX(ii)
         ap = max(0.0, tmp_aR(ii), tmp_aL(ii+1))
         am = max(0.0, tmp_aR(ii-1), tmp_aL(ii))
         if (ii==0) then
            ! seaward boundary condition
            duhydro(:,ii) = 0.0d0
         else if (ii==m) then
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

   !------------------------------------------------------------------------------
   !
   ! Purpose: The 4th-order time step variable Runge-Kutta-Fehlberg method 
   !
   !------------------------------------------------------------------------------
   subroutine RK4Fehlberg(invars, mode, tol, curstep, outvars, &
                          ocurstep, nextstep, outerr, n, m)
      implicit none
      !f2py intent(callback) :: odeFunc
      external :: odeFunc
      !f2py real(kind=8), dimension(xn,xm), intent(in) :: x
      !f2py real(kind=8), dimension(xn,xm), intent(out) :: y
      !f2py integer, intent(hide), depend(x) :: xn = shape(x,0)
      !f2py integer, intent(hide), depend(x) :: xm = shape(x,1)
      !f2py call odeFunc(x,y,xn,xm)
      !f2py intent(in) :: invars, mode, tol, curstep
      !f2py intent(out) :: outvars, ocurstep, nextstep, outerr
      !f2py intent(hide), depend(invars) :: n = shape(invars,0)
      !f2py intent(hide), depend(invars) :: m = shape(invars,1)
      real(kind=8), dimension(n,m) :: invars, outvars
      real(kind=8), dimension(n) :: tol
      integer :: mode, n, m, outerr
      real(kind=8) :: curstep, ocurstep, nextstep
      ! local variables
      real(kind=8), dimension(n) :: dy, rdy, dyn
      real(kind=8), dimension(n) :: rel_tol
      real(kind=8), dimension(n) :: abs_rate
      real(kind=8), dimension(n) :: rel_rate
      real(kind=8) :: step, rate, delta
      logical  :: isLargeErr, isConstrainBroken
      integer  :: iter, ii

      isLargeErr = .True.
      isConstrainBroken = .False.
      outerr = 0
      step = curstep
      iter = 1
      rel_tol = TOL_REL
      call odeFunc(invars, rk4_K1)
      do while (isLargeErr .or. isConstrainBroken)
         if (iter>MAXITER) then
            outerr = 1
            return
         end if
         ocurstep = step
         rk4_interim = invars + step*0.25*rk4_K1
         call odeFunc(rk4_interim, rk4_K2)
         rk4_interim = invars + step*(0.09375*rk4_K1+0.28125*rk4_K2)
         call odeFunc(rk4_interim, rk4_K3)
         rk4_interim = invars + step*(0.87938*rk4_K1-3.27720*rk4_K2+ &
            3.32089*rk4_K3)
         call odeFunc(rk4_interim, rk4_K4)
         rk4_interim = invars + step*(2.03241*rk4_K1-8.0*rk4_K2+ &
            7.17349*rk4_K3-0.20590*rk4_K4)
         call odeFunc(rk4_interim, rk4_K5)
         rk4_nxt4th = invars + step*(0.11574*rk4_K1+0.54893*rk4_K3+ &
            0.53533*rk4_K4-0.2*rk4_K5)
         if (mode==fixed_mode) then
            nextstep = step
            outvars = rk4_nxt4th
            return
         end if
         rk4_interim = invars + step*(-0.29630*rk4_K1+2.0*rk4_K2- &
            1.38168*rk4_K3+0.45297*rk4_K4-0.275*rk4_K5)
         call odeFunc(rk4_interim, rk4_K6)
         rk4_nxt5th = invars + step*(0.11852*rk4_K1+0.51899*rk4_K3+ &
            0.50613*rk4_K4-0.18*rk4_K5+0.03636*rk4_K6)
         rk4_rerr = (rk4_nxt4th - rk4_nxt5th) / (rk4_nxt4th + INFTSML)
         call Norm(rk4_rerr, 1, rdy)
         call Norm(rk4_nxt4th-rk4_nxt5th, 1, dy)
         call Minimum(rk4_nxt4th, 1, dyn)
         ! check whether solution is converged
         isLargeErr = .False.
         isConstrainBroken = .False.
         do ii = 1, n, 1
            if (dy(ii)>tol(ii) .and. rdy(ii)>rel_tol(ii)) then
               isLargeErr = .True.
            end if
            if (dyn(ii)<-100*tol(ii)) then
               isConstrainBroken = .True.
            end if
         end do
         ! update time step
         if (isConstrainBroken) then
            step = 0.5*step
         else
            abs_rate = tol / (dy + INFTSML)
            rel_rate = rel_tol / (rdy + INFTSML)
            rate = max(minval(abs_rate), minval(rel_rate))
            delta = 0.84*rate**0.25
            if (delta<=0.1) then
               step = 0.1*step
            else if (delta>=4.0) then
               step = 4.0*step
            else
               step = delta*step
            end if
         end if
         iter = iter + 1
      end do
      nextstep = step
      outvars = rk4_nxt4th
   end subroutine

end module TAIHydroMOD

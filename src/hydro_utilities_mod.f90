module hydro_utilities_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module contains hydrodynamic utility subroutines in the model
!
!---------------------------------------------------------------------------------
   use data_buffer_mod,    only : par_d50, par_Cz0, par_Kdf, par_cbc, par_fr
   use data_buffer_mod,    only : par_alphaA, par_betaA, par_alphaD, par_betaD
   use data_buffer_mod,    only : par_cD0, par_ScD, par_cwc
   use data_buffer_mod,    only : rk4_K1, rk4_K2, rk4_K3, rk4_K4, rk4_K5
   use data_buffer_mod,    only : rk4_K6, rk4_nxt4th, rk4_nxt5th
   use data_buffer_mod,    only : rk4_interim, rk4_rerr

   implicit none
   public
   ! model control constants
   integer, parameter :: MAXITER = 100
   integer, parameter :: adaptive_mode = 101, fixed_mode = 102
   real(kind=8), parameter :: TOL_REL = 1.d-6
   real(kind=8), parameter :: INFTSML = 1.d-30
   real(kind=8), parameter :: INFNT = 1.d+30
   ! physical constants 
   real(kind=8), parameter :: PI = 3.14159265d+0
   real(kind=8), parameter :: e = 2.71828183d+0
   real(kind=8), parameter :: Roul = 1028.0
   real(kind=8), parameter :: Roua = 1.225
   real(kind=8), parameter :: Karman = 0.41
   real(kind=8), parameter :: G = 9.8

contains
   subroutine Norm(matrix, dir, values)
      implicit none
      real(kind=8), intent(in)  :: matrix(:,:)
      integer, intent(in)   :: dir
      real(kind=8), intent(out) :: values(:)
      integer :: ii, nn

      nn = size(matrix,dir)
      if (dir==1) then
         do ii = 1, nn, 1
            values(ii) = max( abs(minval(matrix(ii,:))), &
               abs(maxval(matrix(ii,:))) )
         end do
      else if (dir==2) then
         do ii = 1, nn, 1
            values(ii) = max( abs(minval(matrix(:,ii))), &
               abs(maxval(matrix(:,ii))) )
         end do
      end if
   end subroutine

   subroutine Minimum(matrix, dir, values)
      implicit none
      real(kind=8), intent(in) :: matrix(:,:)
      integer, intent(in) :: dir
      real(kind=8), intent(out) :: values(:)
      integer :: ii, nn

      nn = size(matrix,dir)
      if (dir==1) then
         do ii = 1, nn, 1
            values(ii) = minval(matrix(ii,:))
         end do
      else if (dir==2) then
         do ii = 1, nn, 1
            values(ii) = minval(matrix(:,ii))
         end do
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Solve non-linear equation using Brent's method.
   ! https://rosettacode.org/wiki/Roots_of_a_function#Brent.27s_Method
   !
   !------------------------------------------------------------------------------
   subroutine NonLRBrents(NonLREQ, coefs, xbounds, xtol, ytol, root, err)
      implicit none
      external :: NonLREQ
      real(kind=8), intent(in) :: coefs(:)
      real(kind=8), intent(in) :: xbounds(2)
      real(kind=8), intent(in) :: xtol
      real(kind=8), intent(in) :: ytol
      real(kind=8), intent(out) :: root
      integer, intent(out) :: err
      real(kind=8) :: xa, xb, ya, yb
      real(kind=8) :: xc, yc, xd, ys
      integer :: iter
      logical :: mflag

      err = 0
      xa = xbounds(1)
      xb = xbounds(2)
      call NonLREQ(xa, coefs, ya)
      call NonLREQ(xb, coefs, yb)
      if (abs(ya)<ytol .or. abs(yb)<ytol) then
         ! root at the boundary
         if (abs(ya)<=abs(yb)) then
            root = xa
         else
            root = xb
         end if
      else if (ya*yb>0) then
         ! root not available but set one
         if (abs(ya)<=abs(yb)) then
            root = xa
         else
            root = xb
         end if
         err = 1
      else
         if (abs(ya)<abs(yb)) then
            call swap(xa, xb)
            call swap(ya, yb)
         end if
         xc = xa
         yc = ya
         xd = 0.0 ! xd only used if mflag=false, not on first iteration
         mflag = .True.
         do iter = 1, MAXITER, 1
            ! stop if converged on root or error is less than tolerance
            if (abs(xb-xa)<xtol) then
               exit
            end if
            if (ya /= yc .and. yb /= yc) then
               ! use inverse qudratic interpolation
               root = xa*yb*yc/((ya-yb)*(ya-yc)) + &
                  xb*ya*yc/((yb-ya)*(yb-yc)) + xc*ya*yb/((yc-ya)*(yc-yb))
            else
               ! secant method
               root = xb - yb*(xb-xa)/(yb-ya)
            end if
            ! check to see whether we can use the faster converging quadratic
            ! && secant methods or if we need to use bisection
            if ( ( root<(3.0*xa+xb)*0.25 .or. root>xb ) .or. &
                  ( mflag .and. abs(root-xb)>=abs(xb-xc)*0.5 ) .or. &
                  ( (.NOT. mflag) .and. abs(root-xb)>=abs(xc-xd)*0.5 ) .or. &
                  ( mflag .and. abs(xb-xc)<xtol ) .or. &
                  ( (.NOT. mflag) .and. abs(xc-xd)<xtol ) ) then
               ! bisection method
               root = 0.5 * (xa + xb)
               mflag = .True.
            else
               mflag = .False.
            end if
            call NonLREQ(root, coefs, ys)
            xd = xc  ! first time xd is being used
            xc = xb  ! set xc equal to upper bound
            yc = yb
            ! update the bounds
            if (ya*ys<0) then
               xb = root
               yb = ys
            else
               xa = root
               ya = ys
            end if
            if (abs(ya)<abs(yb)) then
               call swap(xa, xb)
               call swap(ya, yb)
            end if
         end do
      end if
   end subroutine

   subroutine swap(a, b)
      implicit none
      real(kind=8), intent(inout) :: a
      real(kind=8), intent(inout) :: b
      real(kind=8) :: c

      c = a
      a = b
      b = c
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Slope limiter for the 1D Finite Volume Semi-discrete Kurganov and
   !          Tadmor (KT) central scheme.
   !
   !------------------------------------------------------------------------------
   subroutine FVSKT_Superbee(uhydro, phi, n, m)
      implicit none
      real(kind=8), intent(in) :: uhydro(n,m)
      real(kind=8), intent(out) :: phi(n,m)
      integer, intent(in) :: n, m
      real(kind=8) :: rr(m)
      real(kind=8) :: r0, r1
      integer :: ii, jj

      do ii = 1, n, 1
         if (ii==1 .or. ii==n) then
            phi(ii,:) = 0.0d0
         else
            do jj = 1, m, 1
               r0 = uhydro(ii,jj) - uhydro(ii-1,jj)
               r1 = uhydro(ii+1,jj) - uhydro(ii,jj)
               if (abs(r0)>INFTSML .and. abs(r1)>INFTSML) then
                  rr(jj) = r0 / r1
               else if (abs(r0)<=INFTSML) then
                  rr(jj) = 0.0
               else
                  rr(jj) = 1.0 / INFTSML
               end if
            end do
            phi(ii,:) = max(0.0,min(2.0*rr,1.0),min(rr,2.0))
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge state variables.
   !
   !------------------------------------------------------------------------------
   subroutine FVSKT_celledge(uhydro, phi, uhydroL, uhydroR, n, m)
      implicit none
      real(kind=8), intent(in) :: uhydro(n,m)
      real(kind=8), intent(in) :: phi(n,m)
      real(kind=8), intent(out) :: uhydroL(n,m)
      real(kind=8), intent(out) :: uhydroR(n,m)
      integer, intent(in) :: n, m
      real(kind=8) :: duhydro(m)
      integer :: ii

      do ii = 1, n, 1
         if (ii==1 .or. ii==n) then
            uhydroL(ii,:) = uhydro(ii,:)
            uhydroR(ii,:) = uhydro(ii,:)
         else
            duhydro = 0.5*phi(ii,:)*(uhydro(ii+1,:)-uhydro(ii,:))
            uhydroL(ii,:) = uhydro(ii,:) - duhydro
            uhydroR(ii,:) = uhydro(ii,:) + duhydro
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Performs a binary search of a sorted one-dimensional array for a
   !          specified element.
   !
   !------------------------------------------------------------------------------
   subroutine BinarySearch(array, obj, idx)
      implicit none
      real(kind=8), intent(in) :: array(:)
      real(kind=8), intent(in) :: obj
      integer, intent(out) :: idx
      integer :: middle, first, last
      logical :: ascend

      first = 1
      last = size(array)
      ascend = (array(first)<array(last))
      do while (last>first)
         middle = (first+last)/2
         if (array(middle)==obj) then
            last = middle
            exit
         else if (array(middle)<obj) then
            if (ascend) then
               first = middle + 1
            else
               last = middle
            end if
         else
            if (ascend) then
               last = middle
            else
               first = middle + 1
            end if
         end if
      end do
      idx = last - 1
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
   subroutine UpdateGroundRoughness(pft, Bag, h, Cz)
      implicit none
      integer, intent(in) :: pft(:)
      real(kind=8), intent(in) :: Bag(:)
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(out) :: Cz(:)     ! inverse of square of Cz
      real(kind=8), parameter :: Cb = 2.5d-3    ! bed drag coefficient (unitless)
      real(kind=8) :: cD0, ScD, alphaA, alphaD
      real(kind=8) :: betaA, betaD
      real(kind=8) :: asb, dsb, cD
      integer :: nx, ii

      nx = size(h)
      do ii = 1, nx, 1
         cD0 = par_cD0(pft(ii)+1)
         ScD = par_ScD(pft(ii)+1)
         alphaA = par_alphaA(pft(ii)+1)
         betaA = par_betaA(pft(ii)+1)
         alphaD = par_alphaD(pft(ii)+1)
         betaD = par_betaD(pft(ii)+1)
         asb = alphaA * Bag(ii)**betaA
         dsb = alphaD * Bag(ii)**betaD
         cD = cD0 + ScD * Bag(ii)
         !Cz(ii) = par_Cz0*sqrt(2.0/(cD*asb*h(ii)+2.0*(1.0-asb*dsb)*Cb))
         Cz(ii) = (0.5*cD*asb*h(ii)+(1.0-asb*dsb)*Cb)/par_Cz0**2
      end do
   end subroutine

   subroutine UpdateShearStress(Twav, h, U, Hwav, Uwav, tau)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: U(:)
      real(kind=8), intent(in) :: Hwav(:)
      real(kind=8), intent(in) :: Uwav(:)
      real(kind=8), intent(out) :: tau(:)
      real(kind=8) :: fcurr, fwave
      real(kind=8) :: tau_curr, tau_wave
      integer :: ii, nx

      nx = size(h)
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            ! bottom shear stress by currents
            fcurr = 0.24/(log(4.8*max(1.,h(ii)/par_d50)))**2
            tau_curr = 0.125*Roul*fcurr*U(ii)**2
            ! bottom shear stress by wave
            fwave = 1.39*(6.0*Twav/PI/par_d50)**(-0.52)
            tau_wave = 0.5*fwave*Roul*Uwav(ii)**1.48
            ! combined shear stress
            tau(ii) = tau_curr*(1.0+1.2*(tau_wave/(tau_curr+tau_wave))**3.2)
            tau(ii) = max(tau(ii), tau_curr+tau_wave)
         else
            tau(ii) = 0.0d0
         end if
      end do
   end subroutine

   subroutine CalcWaveReductionByVeg(x, dx, Bag, xref, fwave)
      implicit none
      real(kind=8), intent(in) :: x(:)
      real(kind=8), intent(in) :: dx(:)
      real(kind=8), intent(in) :: Bag(:)
      real(kind=8), intent(in) :: xref
      real(kind=8), intent(out) :: fwave(:)
      real(kind=8) :: Latt, xLatt
      integer :: ii, nx

      nx = size(x)
      do ii = 1, nx, 1
         if (x(ii)>xref) then
            if (fwave(ii-1)<TOL_REL) then
               fwave(ii) = 0.0
            else
               if (Bag(ii)>TOL_REL) then
                  xLatt = 1.0/(0.03*Bag(ii)/2.0)
                  fwave(ii) = fwave(ii-1) - 0.5*xLatt/dx(ii)
                  fwave(ii) = max(0.0, fwave(ii))
               else
                  fwave(ii) = fwave(ii-1)
               end if
            end if
         else
            fwave(ii) = 1.0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave number.
   !          UpdateWaveNumber2() is much more efficient but less accurate.
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

   subroutine UpdateWaveNumber(Twav, h, kwav)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: h
      real(kind=8), intent(out) :: kwav
      real(kind=8) :: coefs(2), xbounds(2)
      real(kind=8) :: sigma, xtol, ytol
      character(len=128) :: msg
      integer :: err

      sigma = 2.0*PI/Twav     ! wave frequency (dispersion)
      xtol = 1d-4
      ytol = 1d-10
      if (h>TOL_REL) then
         xbounds = (/2.51d-2, 6.2832d0/)
         coefs = (/Twav, h/)
         call NonLRBrents(WaveNumberEQ, coefs, xbounds, xtol, &
                           ytol, kwav, err)
         !if (err==1) then
         !   write(msg, "(F8.4, F8.4)") Twav, h
         !   print *, "Wave number isn't available: " // trim(msg)
         !end if
      else
         kwav = INFNT
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave depth breaking possibility.
   !          UpdateWaveBrkProb2() is much more efficient but less accurate.
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

   subroutine UpdateWaveBrkProb(h, Hwav, Qb)
      implicit none
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: Hwav(:)
      real(kind=8), intent(out) :: Qb(:)
      real(kind=8) :: xbounds(2), coefs(2)
      real(kind=8) :: xtol, ytol, Hrms, Hmax
      character(len=128) :: msg
      integer :: ii, nx, err

      nx = size(h)
      xtol = 1d-4
      ytol = 1d-10
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            Hmax = par_fr * h(ii)
            Hrms = Hwav(ii)
            xbounds = (/1d-10, 1.0-1d-10/)
            coefs = (/Hmax, Hrms/)
            call NonLRBrents(BreakProbEQ, coefs, xbounds, xtol, &
                             ytol, Qb(ii), err)
            if (err==1) then
               write(msg, "(E10.2E2, E10.2E2)") Hmax, Hrms
               print *, "Qb isn't available: " // trim(msg)
            end if
         else
            Qb(ii) = 1.0d0
         end if
      end do
   end subroutine

   subroutine UpdateWaveBrkProb2(h, Hwav, Qb)
      implicit none
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: Hwav(:)
      real(kind=8), intent(out) :: Qb(:)
      real(kind=8), parameter :: rawQb(101) = (/0.4,0.4637,0.5005,0.5260, &
         0.5461,0.5631,0.5780,0.5914,0.6035,0.6147,0.6252,0.6350,0.6442, &
         0.6530,0.6614,0.6694,0.6770,0.6844,0.6915,0.6984,0.7050,0.7115, &
         0.7177,0.7238,0.7298,0.7355,0.7412,0.7467,0.7521,0.7573,0.7625, &
         0.7676,0.7725,0.7774,0.7822,0.7869,0.7915,0.7960,0.8005,0.8049, &
         0.8092,0.8135,0.8177,0.8218,0.8259,0.8299,0.8339,0.8378,0.8417, &
         0.8455,0.8493,0.8531,0.8568,0.8604,0.8640,0.8676,0.8711,0.8746, &
         0.8781,0.8815,0.8849,0.8883,0.8916,0.8949,0.8981,0.9014,0.9046, &
         0.9078,0.9109,0.9140,0.9171,0.9202,0.9232,0.9262,0.9292,0.9322, &
         0.9352,0.9381,0.9410,0.9439,0.9467,0.9496,0.9524,0.9552,0.9580, &
         0.9607,0.9635,0.9662,0.9689,0.9716,0.9742,0.9769,0.9795,0.9821, &
         0.9847,0.9873,0.9899,0.9924,0.9950,0.9975,1.0/)
      real(kind=8) :: fHrms
      integer :: ii, nx, idx

      nx = size(h)
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            fHrms = Hwav(ii) / (par_fr*h(ii))
            if (fHrms<=rawQb(1)) then
               Qb(ii) = 0.0d0
            else if (fHrms>=rawQb(101)) then
               Qb(ii) = 1.0d0
            else
               call BinarySearch(rawQb, fHrms, idx)
               Qb(ii) = 0.01*DBLE(idx) - 0.005
            end if
         else
            Qb(ii) = 1.0d0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate stead-state wave regime.
   !
   !------------------------------------------------------------------------------
   subroutine SgnftWaveHeightEQ(Ewav, coefs, fval)
      implicit none
      real(kind=8), intent(in) :: Ewav
      real(kind=8), intent(in) :: coefs(4)
      real(kind=8), intent(out) :: fval
      real(kind=8), parameter :: rawQb(101) = (/0.4,0.4637,0.5005,0.5260, &
         0.5461,0.5631,0.5780,0.5914,0.6035,0.6147,0.6252,0.6350,0.6442, &
         0.6530,0.6614,0.6694,0.6770,0.6844,0.6915,0.6984,0.7050,0.7115, &
         0.7177,0.7238,0.7298,0.7355,0.7412,0.7467,0.7521,0.7573,0.7625, &
         0.7676,0.7725,0.7774,0.7822,0.7869,0.7915,0.7960,0.8005,0.8049, &
         0.8092,0.8135,0.8177,0.8218,0.8259,0.8299,0.8339,0.8378,0.8417, &
         0.8455,0.8493,0.8531,0.8568,0.8604,0.8640,0.8676,0.8711,0.8746, &
         0.8781,0.8815,0.8849,0.8883,0.8916,0.8949,0.8981,0.9014,0.9046, &
         0.9078,0.9109,0.9140,0.9171,0.9202,0.9232,0.9262,0.9292,0.9322, &
         0.9352,0.9381,0.9410,0.9439,0.9467,0.9496,0.9524,0.9552,0.9580, &
         0.9607,0.9635,0.9662,0.9689,0.9716,0.9742,0.9769,0.9795,0.9821, &
         0.9847,0.9873,0.9899,0.9924,0.9950,0.9975,1.0/) 
      real(kind=8), parameter :: Cd = 1.3d-3    ! drag coefficient for U10
      real(kind=8), parameter :: gammaPM = 0.055
      real(kind=8), parameter :: m = 2.0
      real(kind=8) :: Hrms, Hmax, fHrms, Qb
      real(kind=8) :: Twav, U10, kwav, h
      real(kind=8) :: sigma, alpha, beta, Cf
      real(kind=8) :: gammaR
      real(kind=8) :: Swg, Sbf, Swc, Sbrk
      integer :: idx

      Twav = coefs(1)
      U10 = coefs(2)
      h = coefs(3)
      kwav = coefs(4)
      Hmax = par_fr * h
      Hrms = sqrt(8.0*Ewav/G/Roul) 
      fHrms = Hrms / Hmax
      if (fHrms<=rawQb(1)) then
         Qb = 0.0d0
      else if (fHrms>=rawQb(101)) then
         Qb = 1.0d0
      else
         call BinarySearch(rawQb, fHrms, idx)
         Qb = 0.01*DBLE(idx) - 0.005
      end if
      sigma = 2.0*PI/Twav
      if (h>TOL_REL) then
         ! wave generation
         alpha = 80.0*sigma*(Roua*Cd*U10/Roul/G/kwav)**2
         beta = 5.0*Roua/Roul/Twav*max(0.0,U10*kwav/sigma-0.9)
         Swg = alpha + beta * Ewav
         ! wave reduction by bottom friction
         Cf = par_cbc*Hrms*sigma/sinh(kwav*h)
         Sbf = (1-Qb)*2.0*Cf*kwav*Ewav/sinh(2.0*kwav*h)
         ! wave reduction by white capping
         gammaR = Ewav*(sigma**4)/G**2
         Swc = par_cwc*sigma*Ewav*(gammaR/gammaPM)**m
         ! wave reduction by breaking
         Sbrk = 2.0*alpha/Twav*Qb*((Hmax/(1d-12+Hrms))**2)*Ewav
         ! wave energy balance equation
         fval = Swg - Sbf - Swc - Sbrk
      else
         fval = 0.0
      end if
   end subroutine

   subroutine UpdateSgnftWaveHeight(Twav, U10, h, kwav, Ewav)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: U10
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: kwav(:)
      real(kind=8), intent(out) :: Ewav(:)      ! wave energy
      real(kind=8) :: xbounds(2), coefs(4)
      real(kind=8) :: xtol, ytol
      character(len=128) :: msg
      integer :: ii, nx, err

      nx = size(h)
      xtol = 1d-6
      ytol = 1d-15
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            xbounds = (/0.0d0, 4.9d5/)
            coefs = (/Twav, U10, h(ii), kwav(ii)/)
            call NonLRBrents(SgnftWaveHeightEQ, coefs, xbounds, xtol, &
                             ytol, Ewav(ii), err)
            if (err==1) then
               write(msg, "(E10.2E2, E10.2E2, E10.2E2, E10.2E2)") &
                  Twav, U10, h(ii), kwav(ii)
               print *, "Ewav isn't available: " // trim(msg)
            end if
         else
            Ewav(ii) = 0d0
         end if 
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave dynamics using an empirical method (Tambroni & 
   !          Seminara, 2012, JGR).
   !
   !------------------------------------------------------------------------------
   subroutine UpdateSgnftWaveHeight2(U10, xfetch, h, Hwav, Twav)
      implicit none
      real(kind=8), intent(in) :: U10
      real(kind=8), intent(in) :: xfetch
      real(kind=8), intent(in) :: h
      real(kind=8), intent(out) :: Hwav   ! units: meter
      real(kind=8), intent(out) :: Twav   ! units: second
      real(kind=8) :: A1, A2, B1, B2

      if (h>TOL_REL) then
         A1 = 0.493*(h*G/U10**2)**0.75
         B1 = 3.13d-3*(xfetch*G/U10**2)**0.57
         A2 = 0.331*(h*G/U10**2)**1.01
         B2 = 5.21d-4*(xfetch*G/U10**2)**0.73
         Hwav = 0.17*U10**2/G*(tanh(A1)*tanh(B1/tanh(A1)))**0.87
         Twav = 7.518*U10/G*(tanh(A2)*tanh(B2/tanh(A2)))**0.37
      else
         Hwav = 0.0
         Twav = 2.0
      end if
   end subroutine   

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate wave sources and sinks.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateWaveGeneration(Twav, U10, h, kwav, Ewav, Swg)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: U10
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: kwav(:)
      real(kind=8), intent(in) :: Ewav(:)
      real(kind=8), intent(out) :: Swg(:)
      real(kind=8), parameter :: Cd = 1.3d-3    ! drag coefficient for U10
      real(kind=8) :: sigma, alpha, beta
      integer :: ii, nx

      nx = size(h)
      sigma = 2.0*PI/Twav
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            alpha = 80.0*sigma*(Roua*Cd*U10/Roul/G/kwav(ii))**2
            beta = 5.0*Roua/Roul/Twav*max(0.0,U10*kwav(ii)/sigma-0.9)
            Swg(ii) = alpha + beta * Ewav(ii)
         else
            Swg(ii) = 0.0d0
         end if
      end do
   end subroutine

   subroutine UpdateWaveBtmFriction(Twav, h, Hwav, kwav, Ewav, &
                                    Qb, Sbf)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: Hwav(:)
      real(kind=8), intent(in) :: kwav(:)
      real(kind=8), intent(in) :: Ewav(:)
      real(kind=8), intent(in) :: Qb(:)
      real(kind=8), intent(out) :: Sbf(:)
      real(kind=8) :: Cf
      integer :: ii, nx

      nx = size(h)
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            Cf = 2.0*par_cbc*PI*Hwav(ii)/Twav/sinh(kwav(ii)*h(ii))
            Sbf(ii) = (1-Qb(ii))*2.0*Cf*kwav(ii)*Ewav(ii)/ &
               sinh(2.0*kwav(ii)*h(ii))
         else
            Sbf(ii) = 0.0d0
         end if
      end do
   end subroutine

   subroutine UpdateWaveWhiteCapping(Twav, Ewav, Swc)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: Ewav(:)
      real(kind=8), intent(out) :: Swc(:)
      real(kind=8), parameter :: gammaPM = 0.055
      real(kind=8), parameter :: m = 2.0d0
      real(kind=8) :: sigma

      sigma = 2.0*PI/Twav
      Swc = par_cwc*sigma*Ewav*(Ewav*(sigma**4)/G**2/gammaPM)**m
   end subroutine

   subroutine UpdateWaveDepthBrking(Twav, U10, h, Hwav, kwav, Ewav, &
                                    Qb, Sbrk)
      implicit none
      real(kind=8), intent(in) :: Twav
      real(kind=8), intent(in) :: U10
      real(kind=8), intent(in) :: h(:)
      real(kind=8), intent(in) :: Hwav(:)
      real(kind=8), intent(in) :: kwav(:)
      real(kind=8), intent(in) :: Ewav(:)
      real(kind=8), intent(in) :: Qb(:)
      real(kind=8), intent(out) :: Sbrk(:)
      real(kind=8), parameter :: Cd = 1.3d-3    ! drag coefficient for U10
      real(kind=8) :: sigma, Hmax, alpha
      integer :: ii, nx

      nx = size(h)
      sigma = 2.0*PI/Twav
      do ii = 1, nx, 1
         if (h(ii)>TOL_REL) then
            Hmax = par_fr * h(ii)
            alpha = 80.0*sigma*(Roua*Cd*U10/Roul/G/kwav(ii))**2
            Sbrk(ii) = 2.0*alpha/Twav*Qb(ii)*((Hmax/Hwav(ii))**2)*Ewav(ii)
         else
            Sbrk(ii) = 0.0d0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: The 4th-order time step variable Runge-Kutta-Fehlberg method 
   !
   !------------------------------------------------------------------------------
   subroutine RK4Fehlberg(odeFunc, invars, mode, tol, dyncheck, outvars, &
                          curstep, nextstep, outerr)
      implicit none
      external :: odeFunc
      real(kind=8), intent(in) :: invars(:,:)
      integer, intent(in) :: mode
      real(kind=8), intent(in) :: tol(:)
      logical, intent(in) :: dyncheck(:)
      real(kind=8), intent(out) :: outvars(:,:)
      real(kind=8), intent(inout) :: curstep
      real(kind=8), intent(out) :: nextstep
      integer, intent(out) :: outerr
      ! local variables
      real(kind=8), dimension(size(tol)) :: dy, rdy, dyn
      real(kind=8), dimension(size(tol)) :: rel_tol
      real(kind=8), dimension(size(tol)) :: abs_rate
      real(kind=8), dimension(size(tol)) :: rel_rate
      real(kind=8) :: step, rate, delta
      logical  :: isLargeErr, isConstrainBroken
      integer  :: iter, ii, n, m

      n = size(invars,1)
      m = size(invars,2)
      isLargeErr = .True.
      isConstrainBroken = .False.
      outerr = 0
      step = curstep
      iter = 1
      rel_tol = TOL_REL
      call odeFunc(invars, rk4_K1, n, m)
      do while (isLargeErr .or. isConstrainBroken)
         if (iter>MAXITER) then
            outerr = 1
            return
         end if
         curstep = step
         rk4_interim = invars + step*0.25*rk4_K1
         call odeFunc(rk4_interim, rk4_K2, n, m)
         rk4_interim = invars + step*(0.09375*rk4_K1+0.28125*rk4_K2)
         call odeFunc(rk4_interim, rk4_K3, n, m)
         rk4_interim = invars + step*(0.87938*rk4_K1-3.27720*rk4_K2+ &
            3.32089*rk4_K3)
         call odeFunc(rk4_interim, rk4_K4, n, m)
         rk4_interim = invars + step*(2.03241*rk4_K1-8.0*rk4_K2+ &
            7.17349*rk4_K3-0.20590*rk4_K4)
         call odeFunc(rk4_interim, rk4_K5, n, m)
         rk4_nxt4th = invars + step*(0.11574*rk4_K1+0.54893*rk4_K3+ &
            0.53533*rk4_K4-0.2*rk4_K5)
         if (mode==fixed_mode) then
            nextstep = step
            outvars = rk4_nxt4th
            return
         end if
         rk4_interim = invars + step*(-0.29630*rk4_K1+2.0*rk4_K2- &
            1.38168*rk4_K3+0.45297*rk4_K4-0.275*rk4_K5)
         call odeFunc(rk4_interim, rk4_K6, n, m)
         rk4_nxt5th = invars + step*(0.11852*rk4_K1+0.51899*rk4_K3+ &
            0.50613*rk4_K4-0.18*rk4_K5+0.03636*rk4_K6)
         rk4_rerr = (rk4_nxt4th - rk4_nxt5th) / (rk4_nxt4th + INFTSML)
         call Norm(rk4_rerr, 2, rdy)
         call Norm(rk4_nxt4th-rk4_nxt5th, 2, dy)
         call Minimum(rk4_nxt4th, 2, dyn)
         ! check whether solution is converged
         isLargeErr = .False.
         isConstrainBroken = .False.
         do ii = 1, m, 1
            if (dy(ii)>tol(ii) .and. rdy(ii)>rel_tol(ii)) then
               isLargeErr = .True.
            end if
            if (dyn(ii)<-100*tol(ii) .and. dyncheck(ii)) then
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

end module hydro_utilities_mod

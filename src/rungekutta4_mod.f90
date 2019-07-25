module RungeKutta4
!---------------------------------------------------------------------------------
! Purpose:
!
! This module implements 4th-order Runge-Kutta-Fehlberg method
!
!---------------------------------------------------------------------------------
   implicit none
   integer, parameter :: MAXITER = 100
   integer, parameter :: adaptive_mode = 101, fixed_mode = 102
   real(kind=8), parameter :: TOL_REL = 1.d-6
   real(kind=8), parameter :: INFTSML = 1.d-30
   ! allocatable intermediate arrays 
   real(kind=8), allocatable, dimension(:,:) :: K1
   real(kind=8), allocatable, dimension(:,:) :: K2
   real(kind=8), allocatable, dimension(:,:) :: K3
   real(kind=8), allocatable, dimension(:,:) :: K4
   real(kind=8), allocatable, dimension(:,:) :: K5
   real(kind=8), allocatable, dimension(:,:) :: K6
   real(kind=8), allocatable, dimension(:,:) :: nxt4th
   real(kind=8), allocatable, dimension(:,:) :: nxt5th
   real(kind=8), allocatable, dimension(:,:) :: interim
   real(kind=8), allocatable, dimension(:,:) :: rerr
   ! intermediate scalers

contains
   subroutine InitRKDataBuffer(nvar, ncell)
      implicit none
      integer, intent(in) :: nvar
      integer, intent(in) :: ncell

      allocate(K1(nvar,ncell))
      allocate(K2(nvar,ncell))
      allocate(K3(nvar,ncell))
      allocate(K4(nvar,ncell))
      allocate(K5(nvar,ncell))
      allocate(K6(nvar,ncell))
      allocate(nxt4th(nvar,ncell))
      allocate(nxt5th(nvar,ncell))
      allocate(interim(nvar,ncell))
      allocate(rerr(nvar,ncell))
   end subroutine

   subroutine DestructRKDataBuffer()
      implicit none

      deallocate(K1)
      deallocate(K2)
      deallocate(K3)
      deallocate(K4)
      deallocate(K5)
      deallocate(K6)
      deallocate(nxt4th)
      deallocate(nxt5th)
      deallocate(interim)
      deallocate(rerr)
   end subroutine

   subroutine RK4Fehlberg(odeFunc, invars, outvars, mode, tol, &
                          curstep, nextstep, outerr)
      implicit none
      external :: odeFunc
      real(kind=8), intent(in) :: invars(:,:)    
      real(kind=8), intent(inout) :: outvars(:,:)
      integer, intent(in)  :: mode
      real(kind=8), intent(in) :: tol(:)
      real(kind=8), intent(inout) :: curstep(1)
      real(kind=8), intent(inout) :: nextstep(1)
      integer, intent(inout)  :: outerr(1)
      real(kind=8), dimension(size(invars,1)) :: dy, rdy, dyn
      real(kind=8), dimension(size(invars,1)) :: rel_tol
      real(kind=8), dimension(size(invars,1)) :: abs_rate
      real(kind=8), dimension(size(invars,1)) :: rel_rate
      real(kind=8) :: step, rate, delta
      logical  :: isLargeErr, isConstrainBroken
      integer  :: iter, ii, nx, ny

      nx = size(invars,1)
      ny = size(invars,2)
      isLargeErr = .True.
      isConstrainBroken = .False.
      outerr(1) = 0
      step = curstep(1)
      iter = 1
      rel_tol = TOL_REL
      call odeFunc(invars, K1)
      do while (isLargeErr .or. isConstrainBroken)
         if (iter>MAXITER) then
            outerr(1) = 1
            return
         end if
         curstep(1) = step
         interim = invars + step*0.25*K1
         call odeFunc(interim, K2)
         interim = invars + step*(0.09375*K1+0.28125*K2)
         call odeFunc(interim, K3)
         interim = invars + step*(0.87938*K1-3.27720*K2+3.32089*K3)
         call odeFunc(interim, K4)
         interim = invars + step*(2.03241*K1-8.0*K2+7.17349*K3-0.20590*K4)
         call odeFunc(interim, K5)
         nxt4th = invars + step*(0.11574*K1+0.54893*K3+0.53533*K4-0.2*K5)
         if (mode==fixed_mode) then
            nextstep(1) = step
            outvars = nxt4th
            return
         end if
         interim = invars + step*(-0.29630*K1+2.0*K2-1.38168*K3+0.45297*K4-0.275*K5)
         call odeFunc(interim, K6)
         nxt5th = invars + step*(0.11852*K1+0.51899*K3+0.50613*K4-0.18*K5+0.03636*K6)
         rerr = (nxt4th - nxt5th) / (nxt4th + INFTSML)
         call Norm(rerr, 1, rdy)
         call Norm(nxt4th-nxt5th, 1, dy)
         call Minimum(nxt4th, 1, dyn)
         ! check whether solution is converged
         isLargeErr = .False.
         isConstrainBroken = .False.
         do ii = 1, nx, 1
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
      nextstep(1) = step
      outvars = nxt4th
   end subroutine

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
   subroutine NonLRBrents(NonLREQ, coefs, xbounds, tol, root)
      implicit none
      external :: NonLREQ
      real(kind=8), intent(in) :: coefs(:)
      real(kind=8), intent(in) :: xbounds(2) 
      real(kind=8), intent(in) :: tol
      real(kind=8), intent(out) :: root
      real(kind=8) :: xa, xb, ya, yb
      real(kind=8) :: xc, yc, xd, ys
      integer :: iter
      logical :: mflag

      xa = xbounds(1)
      xb = xbounds(2)
      call NonLREQ(xa, coefs, ya)
      call NonLREQ(xb, coefs, yb)
      if (abs(ya)<1d-3 .or. abs(yb)<1d-3) then
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
         print *, "Root is not within xbounds"
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
            if (abs(xb-xa)<tol) then
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
                  ( mflag .and. abs(xb-xc)<tol ) .or. &
                  ( (.NOT. mflag) .and. abs(xc-xd)<tol ) ) then
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
   !			  Tadmor (KT) central scheme.
   !
   !------------------------------------------------------------------------------
   subroutine FVSKT_Superbee(uhydro, phi)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(out) :: phi(:,:)
      real(kind=8) :: rr(size(uhydro,1))
      integer :: ii, NX

      NX = size(uhydro,2)
      do ii = 1, NX, 1
         if (ii==1 .or. ii==NX) then
            phi(:,ii) = 0.0d0
         else
            rr = (uhydro(:,ii)-uhydro(:,ii-1))/(uhydro(:,ii+1)-uhydro(:,ii))
            phi(:,ii) = max(0.0,min(2.0*rr,1.0),min(rr,2.0))            
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate cell edge state variables.
   !
   !------------------------------------------------------------------------------
   subroutine FVSKT_celledge(uhydro, phi, uhydroL, uhydroR)
      implicit none
      real(kind=8), intent(in) :: uhydro(:,:)
      real(kind=8), intent(in) :: phi(:,:)
      real(kind=8), intent(out) :: uhydroL(:,:)
      real(kind=8), intent(out) :: uhydroR(:,:)
      real(kind=8) :: duhydro(size(uhydro,1))
      integer :: ii, NX

      NX = size(uhydro,2)
      do ii = 1, NX, 1
         if (ii==1 .or. ii==NX) then
            uhydroL(:,ii) = uhydro(:,ii)
            uhydroR(:,ii) = uhydro(:,ii)
         else
            duhydro = 0.5*phi(:,ii)*(uhydro(:,ii+1)-uhydro(:,ii))
            uhydroL(:,ii) = uhydro(:,ii) - duhydro 
            uhydroR(:,ii) = uhydro(:,ii) + duhydro 
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

end module RungeKutta4

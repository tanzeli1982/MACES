module hydro_utilities_mod

   implicit none
   public
   integer, parameter :: MAXITER = 100
   integer, parameter :: adaptive_mode = 101, fixed_mode = 102
   real(kind=8), parameter :: TOL_REL = 1.d-6
   real(kind=8), parameter :: INFTSML = 1.d-30
   ! intermediate scalers

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
   !          Tadmor (KT) central scheme.
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

end module hydro_utilities_mod

module RungeKutta4
!---------------------------------------------------------------------------------
! Purpose:
!
! This module implements 4th-order Runge-Kutta-Fehlberg method
!
!---------------------------------------------------------------------------------
   implicit none
   integer, parameter  :: MAXITER = 100
   integer, parameter  :: adaptive_mode = 101, fixed_mode = 102
   real(kind=8), parameter :: TOL_REL = 1.d-6
   real(kind=8), parameter :: INFTSML = 1.d-30
   ! allocatable intermediate arrays 
   real(kind=8), allocatable, dimension(:,:) :: K1
   real(kind=8), allocatable, dimension(:,:) :: K2
   real(kind=8), allocatable, dimension(:,:) :: K3
   real(kind=8), allocatable, dimension(:,:) :: K4
   real(kind=8), allocatable, dimension(:,:) :: K5
   real(kind=8), allocatable, dimension(:,:) :: K6
   real(kind=8), allocatable, dimension(:,:) :: nxt4th(:,:)
   real(kind=8), allocatable, dimension(:,:) :: nxt5th(:,:)
   real(kind=8), allocatable, dimension(:,:) :: interim(:,:)
   real(kind=8), allocatable, dimension(:,:) :: rerr(:,:)
   ! intermediate scalers

contains
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
            print *, "Runge-Kutta iteration is more than MAXITER"
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

   !subroutine odeFunc(matrix, values)
   !   implicit none
   !   real(kind=8), intent(in) :: matrix(:,:)
   !   real(kind=8), intent(inout) :: values(:,:)
!
!      values = 0.0
!   end subroutine

end module RungeKutta4

module data_type_mod
!----------------------------------------------------------------------------
! user self-defined types
!----------------------------------------------------------------------------
   public
   type ModelParams
      real(kind=8) :: d50                 ! sediment median diameter (m)
      real(kind=8) :: Cz0                 ! the Chézy friction coefficient (m^0.5/s)
      real(kind=8), pointer :: cD0(:)     ! bulk plant drag coefficient baseline (unitless)
      real(kind=8), pointer :: ScD(:)     ! the slope between cD and Bag (unknown)
      real(kind=8), pointer :: alphaA(:)  ! empirical coefficient for projected plant area
      real(kind=8), pointer :: betaA(:)   ! empirical coefficient for projected plant area
      real(kind=8), pointer :: alphaD(:)  ! empirical coefficient for stem diameter
      real(kind=8), pointer :: betaD(:)   ! empirical coefficient for stem diameter
      real(kind=8), pointer :: Kdf        ! sediment dispersion coefficient (m^2/s)
      real(kind=8), pointer :: cbc        ! wave dissipation coefficient (unknown)
      real(kind=8), pointer :: fr         ! maximum allowed wave height fraction
   end type

   type ForcingData
      real(kind=8) :: Twav                ! wave period (s)
      real(kind=8) :: U10                 ! wind speed (m/s)
      real(kind=8) :: h0                  ! seaward water level (msl)
      real(kind=8) :: U0                  ! seaward water flow velocity (+/- m/s)
      real(kind=8) :: Hwav0               ! seaward significant wave height (m)
      real(kind=8) :: Css0                ! seaward sediment conc (kg/m3)
      real(kind=8) :: Cj0                 ! seaward salinity (PSU)
      real(kind=8), pointer :: Esed(:)    ! sediment suspension (kg/m2/s)
      real(kind=8), pointer :: Dsed(:)    ! sediment deposition (kg/m2/s)
      real(kind=8), pointer :: Bag(:)     ! aboveground biomass (kg/m2)
      integer, pointer :: pft(:)          ! platform pft
   end type

   type RungeKuttaCache
      real(kind=8), pointer :: K1(:,:)
      real(kind=8), pointer :: K2(:,:)
      real(kind=8), pointer :: K3(:,:)
      real(kind=8), pointer :: K4(:,:)
      real(kind=8), pointer :: K5(:,:)
      real(kind=8), pointer :: K6(:,:)
      real(kind=8), pointer :: nxt4th(:,:)
      real(kind=8), pointer :: nxt5th(:,:)
      real(kind=8), pointer :: interim(:,:)
      real(kind=8), pointer :: rerr(:,:)
   end type
end module data_type_mod

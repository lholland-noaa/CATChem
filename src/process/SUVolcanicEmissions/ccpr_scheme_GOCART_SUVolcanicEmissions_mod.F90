!>
!! \file
!! \brief CCPr Scheme for Volcanic Emissions
!!
!!
!! Reference: Benchmarking GOCART-2G in the Goddard Earth Observing System (GEOS)
!! Allison B. Collow, Peter R. Colarco, Arlindo M. da Silva, Virginie Buchard,
!! Huisheng Bian, M Chin, Sampa Das, Ravi Govindaraju, Dongchul Kim, and Valentina Aquila,
!! Geosci. Model Development, 17, 14431468, 2024
!! https://doi.org/10.5194/gmd-17-1443-2024
!!
!! \author Lacey Holland
!! \date 07/2024
!!!>
module CCPr_Scheme_GOCART_SUVolcanicEmissions_Mod

   implicit none

   private

   public :: CCPr_Scheme_GOCART_SUVolcanicEmissions

contains

   !> \brief Brief description of the subroutine
   !!
   !! \param MetState     Meteorological Variables
   !! \param DiagState    Diagnostic Variables
   !! \param SUVolcanicEmissions  SUVolcanicEmissions Variables
   !! \param RC           Success or Failure
   !!
   !! Note that other state types may be required, e.g. one specific to the process group.
   !!!>

     ! Need to change this so that it includes only what is necessary
     !  May need to do something to connect iPoint, jPoint to vlat, vlon
     ! should only run this where we know there is a volcano????
   subroutine CCPr_Scheme_GOCART_SUVolcanicEmissions(nVolc,              &
      vStart,            &
      vEnd,            &
      vSO2,           &
      vElev,             &
      vCloud,           &
      iPoint,            &
      jPoint,           &
      nhms,      &
      SO2EMVN,              &
      SO2,              &
      nSO2,             &
      SU_emis,         &
      km, &
      cdt,          &
      grav,            &
      hghte,             &
      delp,             &
      area,        &
      vLat,         &
      vLon,         &
      RC)

      !subroutine SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, vElev, vCloud, iPoint, &
      !                             jPoint, nhms, SO2EMVN, SO2EMVE, SO2, nSO2, SU_emis, km, cdt, grav,&
      !                             hghte, delp, area, vLat, vLon, rc)


      USE GOCART2G_Process, only: SUVolcanicEmissions

      IMPLICIT NONE

      ! Arguments
      INTEGER, intent(in)              :: km                ! number of vertical levels
      INTEGER, intent(in)              :: cdt               ! model timestep

      INTEGER, intent(in)              :: vStart      ! Emissions Start time
      INTEGER, intent(in)              :: vEnd        ! Emissions end time


      REAL, allocatable, intent(in), DIMENSION(:) :: delp   ! Pressure Thickness [Pa]
      REAL, allocatable, intent(in), DIMENSION(:) :: hghte  ! Height [m]
      REAL, intent(in)                        :: g0


      ! Output
      integer, intent(out) :: RC                      ! Success or Failure

      ! Local Variables
      real, pointer :: GOCART_HGHTE(:,:,:)
      real, pointer :: GOCART_DELP(:,:,:)

      character(len=256) :: errMsg
      character(len=256) :: thisLoc

      ! Initialize
      errMsg = ''
      thisLoc = ' -> at CCPr_Scheme_GOCART_SUVolcanicEmissions &
             & (in CCPr_Scheme_GOCART_SUVolcanicEmissions_mod.F90)'
      RC = 0

      ! transform data for GOCART SUVolcanicEmissions call

      !  Need to re-write for individual variables
      ! put into the src/core directory as a module
      call PrepMetVarsForGOCARTSUV(km,     &
         tmpu,            &
         hghte,           &
         GOCART_tmpu,     &
         GOCART_HGHTE)

      !------------------
      ! Begin Scheme Code
      !------------------


      call SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, vElev, vCloud, iPoint, &
           jPoint, nhms, SO2EMVN, SO2EMVE, SO2, SUVolcanicEmissions%nSpeciesVolcanic, &
           SU_emis, km, cdt, g0, gocart_hghte, gocart_delp, area, vLat, vLon, rc)


      if (associated(GOCART_DELP)) nullify(GOCART_DELP)
      if (associated(GOCART_TMPU)) nullify(GOCART_TMPU)
     

! End GOCART Code


   end subroutine CCPr_Scheme_GOCART_SUVolcanicEmissions

   !>
   !! \brief PrepMetVarsForGOCART - Prep the meteorological variables for GOCART DryDeposition scheme
   !!
   !! \param [INOUT] metstate
   !! \param [INOUT] tmpu
   !! \param [INOUT] rhoa
   !! \param [INOUT] hghte
   !! \param [INOUT] oro
   !! \param [INOUT] ustar
   !! \param [INOUT] pblh
   !! \param [INOUT] shflux
   !! \param [INOUT] z0h
   !! \param [INOUT] u10m
   !! \param [INOUT] v10m
   !! \param [INOUT] fraclake
   !! \param [INOUT] gwettop
   !! \param [OUT] rc
   !!
   !! \ingroup core_modules
   !!!>
   subroutine PrepMetVarsForGOCARTSUV(km,        &
      delp,            &
      hghte,           &
      GOCART_delp,     &
      GOCART_HGHTE)


      IMPLICIT NONE

      ! INPUTS
      INTEGER, intent(in)                     :: km     ! number of vertical levels
      REAL,  intent(in), DIMENSION(:), target :: delp   ! Temperature [K]
      REAL,  intent(in), DIMENSION(:), target :: hghte  ! Height [m]

      ! INPUT/OUTPUTS
      REAL, intent(inout), pointer :: GOCART_DELP(:,:,:)   !< temperature [K]
      REAL, intent(inout), pointer, DIMENSION(:,:,:) :: GOCART_HGHTE  !< geometric height [m]

      ! OUTPUTS - Add error handling back in late
      !INTEGER :: rc !< Return code

      ! Error handling
      !character(len=255) :: thisloc

      allocate(GOCART_DELP(1, 1, km))
      allocate(GOCART_HGHTE(1, 1, km))

      GOCART_DELP(1,1,:) = delp ! temperature [K]
      GOCART_HGHTE = reshape(hghte, (/1, 1, km/))    ! top of layer geopotential height [m]


   end subroutine PrepMetVarsForGOCART


end module CCPr_Scheme_GOCART_SUVolcanicEmissions


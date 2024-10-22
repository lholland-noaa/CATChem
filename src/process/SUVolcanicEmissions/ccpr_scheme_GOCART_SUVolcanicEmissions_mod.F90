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
   !! \param EmisState    Emission State variables
   !! \param SUVolcanicEmissions  SUVolcanicEmissions Variables
   !! \param RC           Success or Failure
   !!
   !! Note that other state types may be required, e.g. one specific to the process group.
   !!!>

     ! Need to change this so that it includes only what is necessary
     !  May need to do something to connect iPoint, jPoint to vlat, vlon
     ! should only run this where we know there is a volcano????
   subroutine CCPr_Scheme_GOCART_SUVolcanicEmissions( MetState, DiagState, EmisState, &
      SUVolcanicEmissions, &
      RC)

!   call CCPr_Scheme_GOCART_SUVolcanicEmissions( MetState%NLEVS,   &
!           MetState%ZMID, MetState%DELP,    &
!           g0,               &
!           RC)

      ! Below is cut/paste from GOCART
      !subroutine SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, vElev, vCloud, iPoint, &
      !                             jPoint, nhms, SO2EMVN, SO2EMVE, SO2, nSO2, SU_emis, km, cdt, grav,&
      !                             hghte, delp, area, vLat, vLon, rc)


      USE GOCART2G_Process, only: SUVolcanicEmissions, ReadPointEmissions

      IMPLICIT NONE

      ! Arguments
      INTEGER, intent(in)              :: km                ! number of vertical levels
      INTEGER, intent(in)              :: cdt               ! model timestep [sec]

      INTEGER, intent(inout)           :: vStart      ! Emissions Start time
      INTEGER, intent(inout)           :: vEnd        ! Emissions end time
      INTEGER, intent(inout)           :: nVolc       ! number of volcanic sources
      INTEGER, intent(inout)           :: rc          ! error code - is this inout or out???
      INTEGER, intent(in)              :: iPoint, jPoint ! sub-domain - we only run this at the place/time of eruption??
      INTEGER, intent(in)              :: YMD
      INTEGER, intent(in)              :: HMS

      REAL, allocatable, intent(in), DIMENSION(:) :: delp   ! Pressure Thickness [Pa]
      REAL, allocatable, intent(in), DIMENSION(:) :: hghte  ! top of layer geopotential Height [m]
      REAL, intent(in)                        :: g0

      REAL, intent(inout)                     :: vSO2   ! volcanic emissions from file [kg]
 !!!!  Below can be figured out from ChemSpeciesState%nSpeciesSUVolcanicIndex???
 !     REAL, intent(inout)                     :: nSO2      ! index of SO2 relative to other sulfate tracers
      REAL, intent(inout)                     :: SO2       ! SO2 [kg kg-1]
      REAL, intent(inout)                     :: SU_emis   ! SU emissions, kg/m2/s
      REAL, intent(inout)                     :: vCloud    ! top elevation of emissions [m]
      REAL, intent(inout)                     :: vElev     ! bottom elevation of emissions [m]
      REAL, intent(inout)                     :: SO2EMVN   ! non-explosive volcanic emissions [kg m-2 s-1]
      REAL, intent(inout)                     :: SO2EMVE   ! explosive volcanic emissions [kg m-2 s-1]

      REAL, parameter :: fMassSulfur = 32.     !  gram molecular weights of species
      REAL, parameter :: fMassSO2 = 64.     !  gram molecular weights of species
      REAL, parameter :: fMassSO4 = 96.     !  gram molecular weights of species

      ! Output
      integer, intent(out) :: RC                      ! Success or Failure, success is 0

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
      ! put into the src/core directory as a module???
      call PrepMetVarsForGOCARTSUV(km,     &
         tmpu,            &
         hghte,           &
         GOCART_tmpu,     &
         GOCART_HGHTE)

      !------------------
      ! Begin Scheme Code
      !------------------
!!!!!!!!!!!!!!!!!! NEED TO EDIT THESE LINES !!!!!!!!!!!!!!!!!!
  integer,          optional, intent(in ) :: nymd    ! numeric year, month, day???
  integer,          optional, intent(in ) :: nhms    ! numeric hour, minute, second???
!! Above is from github.com/GEOS-ESM/MAPL/blob/main/base/StringTemplate.F90
!! The files only have hourly temporal resolution
!!!!!
    if(MetState%YMD /= nymd) then
       MetState%YMD = nymd
!      Get pointwise SO2 and altitude of volcanoes from a daily file data base
       if(index(self%volcano_srcfilen,'volcanic_') /= 0) then
          ! Need to replace StrTemplate with something else that doesn't rely on MAPL
          call StrTemplate(fname, self%volcano_srcfilen, xid='unknown', &
                            nymd=nymd, nhms=120000 )
          call ReadPointEmissions (YMD, fname, nVolc, vLat, vLon, &
                                   vElev, vCloud, vSO2, vStart, &
                                   vEnd, label='volcano', __RC__)
!!! I don't know if we need the line below, if we are just reading in??
          vSO2 = vSO2 * fMassSO2 / fMassSulfur
!         Special possible case
!!! I don't know why this is here, so I commented it out for now
          if(self%volcano_srcfilen(1:9) == '/dev/null') nVolc = 0  ! option for no volcanoes??
       end if
    end if
!!!!!!!!!!!!!!!!!!!!! NEED TO EDIT ABOVE !!!!!!!!!!!!!!!!!!!!!
!!!! Most of the volcanic stuff, is going directly into the call below

      call SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, vElev, vCloud, iPoint, &
           jPoint, nhms, SO2EMVN, SO2EMVE, SO2, SUVolcanicEmissions%nSpeciesSUVolcanic, &
           SU_emis, km, cdt, g0, gocart_hghte, gocart_delp, area, vLat, vLon, rc)


      if (associated(GOCART_DELP)) nullify(GOCART_DELP)
      if (associated(GOCART_TMPU)) nullify(GOCART_TMPU)
     

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

   ! Need to fix below subroutine to convert one variable at a time.
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
      GOCART_HGHTE(1,1,:) = hghte    ! top of layer geopotential height [m]


   end subroutine PrepMetVarsForGOCARTSUV


end module CCPr_Scheme_GOCART_SUVolcanicEmissions


      module casa_types

      use casavariable
      use phenvariable

      implicit none

      private

      ! Les 19 Jan 2011
!      TYPE (casa_biome)    , save   :: casabiome  !
!      TYPE (casa_pool)     , save   :: casapool   !
!      TYPE (casa_flux)     , save   :: casaflux   !
!      TYPE (casa_met)      , save   :: casamet    !
!      TYPE (casa_balance)  , save   :: casabal    !
!      TYPE (phen_variable) , save   :: phen       !
      TYPE (casa_biome)       :: casabiome  !
      TYPE (casa_pool)        :: casapool   !
      TYPE (casa_flux)        :: casaflux   !
      TYPE (casa_met)         :: casamet    !
      TYPE (casa_balance)     :: casabal    !
      TYPE (phen_variable)    :: phen       !
     !TYPE (casafiles_type)   :: casafile   !

      ! Save these so only have to do the allocation once.
      save casabiome, casapool, &
           casaflux, casamet, casabal, phen!, casafile
      public casabiome, casapool, &
           casaflux, casamet, casabal, phen!, casafile

      end module casa_types
   

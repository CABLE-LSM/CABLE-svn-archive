MODULE Totcnppools_casa_mod

#define  UM_BUILD YES

#ifndef UM_BUILD

IMPLICIT NONE

CONTAINS

  SUBROUTINE totcnppools(kloop,veg,casamet,casapool,bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
                               bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb,bmpsoilocc,bmarea)
  ! this subroutine is temporary, and its needs to be modified for multiple tiles within a cell
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER,                     INTENT(IN)  :: kloop
  TYPE (veg_parameter_type),   INTENT(IN)  :: veg  ! vegetation parameters
  TYPE(casa_pool),             INTENT(IN)  :: casapool
  TYPE(casa_met),              INTENT(IN)  :: casamet
  real,      dimension(5,mvtype,mplant)    :: bmcplant,  bmnplant,  bmpplant
  real,      dimension(5,mvtype,mlitter)   :: bmclitter, bmnlitter, bmplitter
  real,      dimension(5,mvtype,msoil)     :: bmcsoil,   bmnsoil,   bmpsoil
  real,      dimension(5,mvtype)           :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
  real,      dimension(mvtype)             :: bmarea
  ! local variables
  INTEGER  npt,nvt


      bmcplant(kloop,:,:)  = 0.0;  bmnplant(kloop,:,:)  = 0.0; bmpplant(kloop,:,:)  = 0.0
      bmclitter(kloop,:,:) = 0.0;  bmnlitter(kloop,:,:) = 0.0; bmplitter(kloop,:,:) = 0.0
      bmcsoil(kloop,:,:)   = 0.0;  bmnsoil(kloop,:,:)   = 0.0; bmpsoil(kloop,:,:)   = 0.0
      bmnsoilmin(kloop,:)  = 0.0;  bmpsoillab(kloop,:)  = 0.0; bmpsoilsorb(kloop,:) = 0.0;  bmpsoilocc(kloop,:) = 0.0

      bmarea(:) = 0.0

      do npt=1,mp
         nvt=veg%iveg(npt)
         bmcplant(kloop,nvt,:) = bmcplant(kloop,nvt,:)   + casapool%cplant(npt,:) * casamet%areacell(npt)
         bmnplant(kloop,nvt,:) = bmnplant(kloop,nvt,:)   + casapool%nplant(npt,:) * casamet%areacell(npt)
         bmpplant(kloop,nvt,:) = bmpplant(kloop,nvt,:)   + casapool%pplant(npt,:) * casamet%areacell(npt)

         bmclitter(kloop,nvt,:) = bmclitter(kloop,nvt,:) + casapool%clitter(npt,:) * casamet%areacell(npt)
         bmnlitter(kloop,nvt,:) = bmnlitter(kloop,nvt,:) + casapool%nlitter(npt,:) * casamet%areacell(npt)
         bmplitter(kloop,nvt,:) = bmplitter(kloop,nvt,:) + casapool%plitter(npt,:) * casamet%areacell(npt)

         bmcsoil(kloop,nvt,:) = bmcsoil(kloop,nvt,:)     + casapool%csoil(npt,:) * casamet%areacell(npt)
         bmnsoil(kloop,nvt,:) = bmnsoil(kloop,nvt,:)     + casapool%nsoil(npt,:) * casamet%areacell(npt)
         bmpsoil(kloop,nvt,:) = bmpsoil(kloop,nvt,:)     + casapool%psoil(npt,:) * casamet%areacell(npt)

         bmnsoilmin(kloop,nvt)  = bmnsoilmin(kloop,nvt)   + casapool%nsoilmin(npt) * casamet%areacell(npt)
         bmpsoillab(kloop,nvt)  = bmpsoillab(kloop,nvt)   + casapool%psoillab(npt) * casamet%areacell(npt)
         bmpsoilsorb(kloop,nvt) = bmpsoilsorb(kloop,nvt)  + casapool%psoilsorb(npt) * casamet%areacell(npt)
         bmpsoilocc(kloop,nvt)  = bmpsoilocc(kloop,nvt)   + casapool%psoilocc(npt) * casamet%areacell(npt)
         bmarea(nvt)  = bmarea(nvt) + casamet%areacell(npt)
      enddo

      do nvt=1,mvtype
         bmcplant(kloop,nvt,:) = bmcplant(kloop,nvt,:)/bmarea(nvt)
         bmnplant(kloop,nvt,:) = bmnplant(kloop,nvt,:)/bmarea(nvt)
         bmpplant(kloop,nvt,:) = bmpplant(kloop,nvt,:)/bmarea(nvt)

         bmclitter(kloop,nvt,:) = bmclitter(kloop,nvt,:)/bmarea(nvt)
         bmnlitter(kloop,nvt,:) = bmnlitter(kloop,nvt,:)/bmarea(nvt)
         bmplitter(kloop,nvt,:) = bmplitter(kloop,nvt,:)/bmarea(nvt)

         bmcsoil(kloop,nvt,:) = bmcsoil(kloop,nvt,:)/bmarea(nvt)
         bmnsoil(kloop,nvt,:) = bmnsoil(kloop,nvt,:)/bmarea(nvt)
         bmpsoil(kloop,nvt,:) = bmpsoil(kloop,nvt,:)/bmarea(nvt)

         bmnsoilmin(kloop,nvt)  = bmnsoilmin(kloop,nvt)/bmarea(nvt)
         bmpsoillab(kloop,nvt)  = bmpsoillab(kloop,nvt)/bmarea(nvt)
         bmpsoilsorb(kloop,nvt) = bmpsoilsorb(kloop,nvt)/bmarea(nvt)
         bmpsoilocc(kloop,nvt)  = bmpsoilocc(kloop,nvt)/bmarea(nvt)
      enddo

  END SUBROUTINE totcnppools
#endif

END MODULE Totcnppools_casa_mod

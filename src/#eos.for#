#include "defs.h"

!!!   Modified to take energy, entropy, temperature, and pressure from both eos.dat and eos-gluon.dat, then interpolate using fugacity factor given by fug(T0, Time, Teq)
!!!   Additionally modified to use uniform spacing in temperature rather than energy density - Andrew
      
      Subroutine InputRegulatedEOS
      Implicit none

      double precision :: e0, p0, cs2, e0gluon, p0gluon, cs2gluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon      

      character(len=1000) :: find_data_file

      ! read EOS data from binary file
      open(5, file=find_data_file('eos.dat'), status='old',
     &     access='stream')
      read(5) EOST0, EOSTend, PEOSdata, SEOSdata, EEOSdata
      close(5)

      open(6, file=find_data_file('eos-gluon.dat'), status='old',
     &     access='stream')
      read(6) EOSgluonT0, EOSgluonTend, PEOSgluondata, SEOSgluondata,
     &     EEOSgluondata
      close(6)
     
      ! save temp step size
      EOSdT = (EOSTend - EOST0)/(EOSnT - 1)
      EOSgluondT = (EOSgluonTend - EOSgluonT0)/(EOSgluonnT - 1)

      ! Prepare for extrapolating EOS beyond end of table.
      ! Assume a modified power-law ansatz for pressure in the high
      ! energy density limit:
      !
      !   p(e) = e/3 + a*e^b
      !
      ! For b < 1, this converges to the ideal gas limit as e -> inf and
      ! the speed of sound c_s^2 = dp/de -> 1/3 as e -> inf.
      ! The coefficients (a, b) may be determined from the pressure and
      ! its derivative (i.e. the speed of sound) at some reference
      ! energy density e0.  Use the second-to-last table entry as the
      ! reference point and estimate the derivative using a finite
      ! difference.
      e0 = EEOSdata(EOSnT - 1)
      p0 = PEOSdata(EOSnT - 1)
      cs2 = EOScs2data(EOSnT)
      e0gluon = EEOSgluondata(EOSgluonnT - 1)
      p0gluon = PEOSgluondata(EOSgluonnT - 1)
      cs2gluon = EOScs2gluondata(EOSgluonnT)

      ! Use variable 'epow' for the exponent 'b' and 'escale' for the
      ! scale factor 'a'.
      epow = e0*(1 - 3*cs2)/(e0 - 3*p0)
      escale = (p0 - e0/3) * e0**(-epow)
      epowgluon = e0gluon*(1 - 3*cs2gluon)/(e0gluon - 3*p0gluon)
      escalegluon = (p0gluon - e0gluon/3) * e0gluon**(-epowgluon)
      
      ! Having extrapolated p(e), the entropy density follows from the
      ! differential equation
      !
      !   (e + p(e))/s(e) = 1/(ds(e)/de)
      !
      ! (Note: both sides of this equation are the temperature.)
      ! For the chosen p(e) ansatz this diffeq can be analytically
      ! integrated for s(e).  These coefficients (spow, sscale) are used
      ! in SEOSL7, below.)
      spow = 3/(4*(1 - epow))
      sscale = SEOSdata(EOSnT - 1) * (e0**epow/(e0 + p0))**spow
      spowgluon = 3/(4*(1 - epowgluon))
      sscalegluon = SEOSgluondata(EOSgluonnT - 1)
     & * (e0gluon**epowgluon/(e0gluon + p0gluon))**spowgluon

      end


C====EOS from table===================================================
      Double Precision Function fug(T0, Time, Teq)
      double precision :: T0, Time, Teq

      if (Teq .EQ. 0) then
         fug = 1
      else
         fug = 1 - exp((T0 - Time)/Teq)
      end if
          
      return
      end

      Double Precision Function PEOSL7qcd(e) ! for lattice P(e)
      Implicit none

      double precision :: e, TEOSL7qcd

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e .lt. EEOSdata(1)) then
         PEOSL7qcd = PEOSdata(1)*e/EEOSdata(1)
      else if (e.lt.EEOSdata(EOSnT)) then
         call interpCubic(PEOSdata, EOSnT, EOST0, EOSdT,
     &    TEOSL7qcd(e), PEOSL7qcd)
      else
        ! extrapolate, see notes above
        PEOSL7qcd = e/3 + escale*e**epow
      endif
    
      return
      end

      Double Precision Function PEOSL7gluon(e) ! for lattice P(e)
      Implicit none

      double precision :: e, TEOSL7gluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e.lt.EEOSgluondata(1)) then
        PEOSL7gluon = PEOSgluondata(1)*e/EEOSgluondata(1)
      else if (e.lt.EEOSgluondata(EOSgluonnT)) then
         call interpCubic(PEOSgluondata, EOSgluonnT, EOSgluonT0,
     &    EOSgluondT, TEOSL7gluon(e), PEOSL7gluon)
      else
        ! extrapolate, see notes above
         PEOSL7gluon = e/3 + escalegluon*e**epowgluon
      endif     

      return
      end
      
      Double Precision Function PEOSL7(e, Tprop) ! for lattice P(e)
      Implicit none

      double precision :: e, Tprop

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

      double precision fugacity, fug, PEOSL7qcd, PEOSL7gluon
      
      fugacity = fug(T0, Tprop, Teq)
      
      e = abs(e)

      PEOSL7 = e/3d0 !fugacity * PEOSL7qcd(e) + (1 - fugacity) * PEOSL7gluon(e)
      
      return
      end      

      Double Precision Function SEOSL7qcd(e)  ! for lattice S(e)
      Implicit none

      double precision :: e
      double precision :: PEOSL7qcd, TEOSL7qcd
     
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSdata(1)) then
        SEOSL7qcd = SEOSdata(1)*e/EEOSdata(1)
      else if (e.lt.EEOSdata(EOSnT)) then
         call interpCubic(SEOSdata, EOSnT, EOST0, EOSdT, TEOSL7qcd(e),
     &    SEOSL7qcd)
      else
        ! extrapolate, see notes above
        SEOSL7qcd = sscale * ((e + PEOSL7qcd(e))/e**epow)**spow
      endif
      
      return
      end

      Double Precision Function SEOSL7gluon(e)  ! for lattice S(e)
      Implicit none

      double precision :: e
      double precision :: PEOSL7gluon, TEOSL7gluon
     
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSgluondata(1)) then
        SEOSL7gluon = SEOSgluondata(1)*e/EEOSgluondata(1)
      else if (e.lt.EEOSgluondata(EOSgluonnT)) then
         call interpCubic(SEOSgluondata, EOSgluonnT, EOSgluonT0,
     &       EOSgluondT, TEOSL7gluon(e), SEOSL7gluon)
      else
        ! extrapolate, see notes above
         SEOSL7gluon = sscalegluon
     & * ((e + PEOSL7gluon(e))/e**epowgluon)**spowgluon
      endif
      
      return
      end
      
      Double Precision Function SEOSL7(e, Tpin)  ! for lattice S(e)
      Implicit none

      double precision :: e, Tprop
      double precision,optional :: Tpin
      double precision :: PEOSL7, TEOSL7
     
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon


      double precision Time, Teq
      common /Time/ Time, Teq

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

      double precision fugacity, fug, SEOSL7qcd, SEOSL7gluon

      if (present(Tpin)) then
         Tprop = Tpin
      else
         Tprop = T0
      endif
      
      fugacity = fug(T0, Tprop, Teq)

      e = abs(e)

      SEOSL7 = (4d0/3d0)*e/TEOSL7(e)!fugacity * SEOSL7qcd(e) + (1 - fugacity) * SEOSL7gluon(e)
      
      return
      end

      Double Precision Function energyqcd(T) ! Returns energy given temperature -Andrew
      implicit None
      double precision :: T

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT                          ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT                         ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     &     EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
     
      call interpCubic(EEOSdata, EOSnT, EOST0, EOSdT, T, energyqcd)
      
      return
      end
      
      Double Precision Function energygluon(T) ! Returns energy given temperature -Andrew
      implicit None
      double precision :: T

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT

      call interpCubic(EEOSgluondata, EOSgluonnT, EOSgluonT0,
     & EOSgluondT, T, energygluon)

      return
      end
      
      Double Precision Function TEOSL7qcd(e)  ! for lattice T(e)
      Implicit none

      double precision :: e

      double precision :: PEOSL7qcd, SEOSL7qcd, energyqcd

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSdata(1)) then
        TEOSL7qcd = EOST0*e/EEOSdata(1)
      else if (e.lt.EEOSdata(EOSnT)) then
        TEOSL7qcd = energyqcd(EOST0)  ! Segmentation fault if function not called first?? - Andrew
        call invertFunction_binary(energyqcd, EOST0,
     &        EOSTend, 1D-16, 1D-6, e, TEOSL7qcd)
        
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        TEOSL7qcd = (e + PEOSL7qcd(e))/SEOSL7qcd(e)
      endif     
      
      return
      end

      Double Precision Function TEOSL7gluon(e)  ! for lattice T(e)
      Implicit none

      double precision :: e

      double precision :: PEOSL7gluon, SEOSL7gluon, energygluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e.lt.EEOSgluondata(1)) then
        TEOSL7gluon = EOSgluonT0*e/EEOSgluondata(1)
      else if (e.lt.EEOSgluondata(EOSgluonnT)) then
        TEOSL7gluon = energygluon(EOSgluonT0) ! Segmentation fault if function not called first?? - Andrew
        call invertFunction_binary(energygluon, EOSgluonT0,
     &        EOSgluonTend, 1D-16, 1D-6, e, TEOSL7gluon)
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        TEOSL7gluon = (e + PEOSL7gluon(e))/SEOSL7gluon(e)
      endif     
      
      return
      end
      
      Double Precision Function TEOSL7(e, Tprop)  ! for lattice T(e)
      Implicit none

      double precision :: e, Tprop

      double precision :: PEOSL7, SEOSL7

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

      double precision fugacity, fug, TEOSL7qcd, TEOSL7gluon
      
      fugacity = fug(T0, Tprop, Teq)
      
      e = abs(e)

      TEOSL7 = ((30d0/M_PI**2)*(M_HBARC**3)/(2*(3*3-1)+3.5*3*3)*e)
     &  **(1d0/4d0)!fugacity * TEOSL7qcd(e) + (1 - fugacity) * TEOSL7gluon(e)      
      
      return
      end


      Double Precision Function EOScs2qcd(e)  ! for lattice T(e)
      Implicit none

      double precision :: e, PEOSL7qcd, TEOSL7qcd

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSdata(1)) then
        EOScs2qcd = (PEOSdata(2)-PEOSdata(1))/(EEOSdata(2)-EEOSdata(1))
      else if (e.lt.EEOSdata(EOSnT)) then
        call interpCubic(EOScs2data, EOSnT, EOST0, EOSdT, TEOSL7qcd(e),
     &    EOScs2qcd)
        
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        EOScs2qcd = (PEOSL7qcd(e)-PEOSL7qcd(e*0.99))/(e-e*0.99)
      endif     
      
      return
      end


      Double Precision Function EOScs2gluon(e)  ! for lattice T(e)
      Implicit none

      double precision :: e, PEOSL7gluon, TEOSL7gluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSgluondata(1)) then
         EOScs2gluon = (PEOSgluondata(2)-PEOSgluondata(1))/
     &    (EEOSgluondata(2)-EEOSgluondata(1))
      else if (e.lt.EEOSgluondata(EOSnT)) then
         call interpCubic(EOScs2gluondata, EOSgluonnT, EOSgluonT0,
     &    EOSgluondT, TEOSL7gluon(e), EOScs2gluon)
        
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        EOScs2gluon = (PEOSL7gluon(e)-PEOSL7gluon(e*0.99))/(e-e*0.99)
      endif     
      
      return
      end
      

      Double Precision Function EOScs2(e, Tprop)  ! compute cs2 from lattice data
      Implicit none

      double precision :: e, Tprop

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

      double precision fugacity, fug, EOScs2qcd, EOScs2gluon
      
      fugacity = fug(T0, Tprop, Teq)
      
      e = abs(e)

      EOScs2 = 1d0/3d0!fugacity * EOScs2qcd(e) + (1 - fugacity) * EOScs2gluon(e)      
      
      return
      end


      Double Precision Function EOScstilde2qcd(e)  ! for lattice T(e)
      Implicit none

      double precision :: e, PEOSL7qcd, TEOSL7qcd

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSdata(1)) then
        EOScstilde2qcd = PEOSdata(1)/EEOSdata(1)
      else if (e.lt.EEOSdata(EOSnT)) then
        call interpCubic(EOScstilde2data, EOSnT, EOST0, EOSdT,
     &    TEOSL7qcd(e), EOScstilde2qcd)
        
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        EOScstilde2qcd = PEOSL7qcd(e)/e
      endif     
      
      return
      end


      Double Precision Function EOScstilde2gluon(e)  ! for lattice T(e)
      Implicit none

      double precision :: e, PEOSL7gluon, TEOSL7gluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EEOSgluondata(1)) then
         EOScstilde2gluon = PEOSgluondata(1)/EEOSgluondata(1)
      else if (e.lt.EEOSgluondata(EOSnT)) then
         call interpCubic(EOScstilde2gluondata, EOSgluonnT, EOSgluonT0,
     &    EOSgluondT, TEOSL7gluon(e), EOScstilde2gluon)
        
      else
        ! use extrapolated pressure and entropy density
        ! T = (e + p)/s
        EOScstilde2gluon = PEOSL7gluon(e)/e
      endif     
      
      return
      end

      

      Double Precision Function EOScstilde2(e, Tprop)  ! compute cs2 from lattice data
      Implicit none

      double precision :: e, Tprop

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    EEOSdata(EOSDATALENGTH),
     &                    EOScs2data(EOSDATALENGTH),
     &                    EOScstilde2data(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    EEOSgluondata(EOSGLUONDATALENGTH),
     &                    EOScs2gluondata(EOSGLUONDATALENGTH),
     &                    EOScstilde2gluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, EEOSdata, EOScs2data,
     &     EOScstilde2data, PEOSgluondata, SEOSgluondata, EEOSgluondata,
     &     EOScs2gluondata, EOScstilde2gluondata

      double precision :: EOST0, EOSe0         ! min (first) temperature
      double precision :: EOSdT, EOSde         ! temperature step
      double precision :: EOSTend, EOSEend       ! max (last) temperature
      Integer :: EOSnT = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluonT0, EOSgluone0        ! min (first) temperature
      double precision :: EOSgluondT, EOSgluonde         ! temperature step
      double precision :: EOSgluonTend, EOSgluonEend       ! max (last) temperature
      Integer :: EOSgluonnT = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOST0, EOSdT, EOSTend,
     & EOSgluonT0, EOSgluondT, EOSgluonTend, EOSnT, EOSgluonnT
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

      double precision fugacity, fug, EOScstilde2qcd, EOScstilde2gluon
      
      fugacity = fug(T0, Tprop, Teq)
      
      e = abs(e)

      EOScstilde2 = 1d0/3d0!fugacity * EOScstilde2qcd(e) +
   !  & (1 - fugacity) * EOScstilde2gluon(e)      
      
      return
      end
      

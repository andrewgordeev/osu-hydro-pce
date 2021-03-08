#include "defs.h"

!!!   Modified to take energy, entropy, temperature, and pressure from both eos.dat and eos-gluon.dat, then interpolate using fugacity factor given by fug(T0, Time, Teq) - Andrew
      
      Subroutine InputRegulatedEOS
      Implicit none

      double precision :: e0, p0, cs2, e0gluon, p0gluon, cs2gluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      

      character(len=1000) :: find_data_file

      ! read EOS data from binary file
      open(5, file=find_data_file('eos.dat'), status='old',
     &     access='stream')
      read(5) EOSe0, EOSEend, PEOSdata, SEOSdata, TEOSdata
      close(5)

      open(6, file=find_data_file('eos-gluon.dat'), status='old',
     &     access='stream')
      read(6) EOSgluone0, EOSgluonEend, PEOSgluondata, SEOSgluondata,
     &     TEOSgluondata
      close(6)
     
      ! save energy density step size
      EOSde = (EOSEend - EOSe0)/(EOSne - 1)
      EOSgluonde = (EOSgluonEend - EOSgluone0)/(EOSgluonne - 1)

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
      e0 = EOSEend - EOSde
      p0 = PEOSdata(EOSne - 1)
      cs2 = (PEOSdata(EOSne) - PEOSdata(EOSne - 2)) / (2*EOSde)
      e0gluon = EOSgluonEend - EOSgluonde
      p0gluon = PEOSgluondata(EOSgluonne - 1)
      cs2gluon = (PEOSgluondata(EOSgluonne) -
     & PEOSgluondata(EOSgluonne - 2)) / (2*EOSgluonde)

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
      sscale = SEOSdata(EOSne - 1) * (e0**epow/(e0 + p0))**spow
      spowgluon = 3/(4*(1 - epowgluon))
      sscalegluon = SEOSgluondata(EOSgluonne - 1)
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

      double precision :: e

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &     TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e.lt.EOSe0) then
        PEOSL7qcd = PEOSdata(1)*e/EOSe0
      else if (e.lt.EOSEend) then
        call interpCubic(PEOSdata, EOSne, EOSe0, EOSde, e, PEOSL7qcd)
      else
        ! extrapolate, see notes above
        PEOSL7qcd = e/3 + escale*e**epow
      endif
    
      return
      end

      Double Precision Function PEOSL7gluon(e) ! for lattice P(e)
      Implicit none

      double precision :: e

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &     TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e.lt.EOSgluone0) then
        PEOSL7gluon = PEOSgluondata(1)*e/EOSgluone0
      else if (e.lt.EOSgluonEend) then
         call interpCubic(PEOSgluondata, EOSgluonne, EOSgluone0,
     &    EOSgluonde, e, PEOSL7gluon)
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
     &     TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
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

      PEOSL7 = fugacity * PEOSL7qcd(e) + (1 - fugacity) * PEOSL7gluon(e)
      
      return
      end      

      Double Precision Function SEOSL7qcd(e)  ! for lattice S(e)
      Implicit none

      double precision :: e
      double precision :: PEOSL7qcd
     
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EOSe0) then
        SEOSL7qcd = SEOSdata(1)*e/EOSe0
      else if (e.lt.EOSEend) then
        call interpCubic(SEOSdata, EOSne, EOSe0, EOSde, e, SEOSL7qcd)
      else
        ! extrapolate, see notes above
        SEOSL7qcd = sscale * ((e + PEOSL7qcd(e))/e**epow)**spow
      endif
      
      return
      end

      Double Precision Function SEOSL7gluon(e)  ! for lattice S(e)
      Implicit none

      double precision :: e
      double precision :: PEOSL7gluon
     
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq

      e = abs(e)

      if (e.lt.EOSgluone0) then
        SEOSL7gluon = SEOSgluondata(1)*e/EOSgluone0
      else if (e.lt.EOSgluonEend) then
         call interpCubic(SEOSgluondata, EOSgluonne, EOSgluone0,
     &       EOSgluonde, e, SEOSL7gluon)
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
      double precision :: PEOSL7
     
      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &                    TEOSdata(EOSDATALENGTH)
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
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

      SEOSL7 = fugacity * SEOSL7qcd(e) + (1 - fugacity) * SEOSL7gluon(e)
      
      return
      end

      Double Precision Function TEOSL7qcd(e)  ! for lattice T(e)
      Implicit none

      double precision :: e

      double precision :: PEOSL7qcd, SEOSL7qcd

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &     TEOSdata(EOSDATALENGTH)      
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e.lt.EOSe0) then
        TEOSL7qcd = TEOSdata(1)*e/EOSe0
      else if (e.lt.EOSEend) then
        call interpCubic(TEOSdata, EOSne, EOSe0, EOSde, e, TEOSL7qcd)
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

      double precision :: PEOSL7gluon, SEOSL7gluon

      double precision :: PEOSdata(EOSDATALENGTH),
     &                    SEOSdata(EOSDATALENGTH),
     &     TEOSdata(EOSDATALENGTH)      
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
      double precision :: escale, epow, sscale, spow
      double precision :: escalegluon, epowgluon, sscalegluon, spowgluon
      common /EOScoeffs/ escale, epow, sscale, spow, escalegluon,
     & epowgluon, sscalegluon, spowgluon

      double precision Time, Teq
      common /Time/ Time, Teq
      
      e = abs(e)

      if (e.lt.EOSgluone0) then
        TEOSL7gluon = TEOSgluondata(1)*e/EOSgluone0
      else if (e.lt.EOSgluonEend) then
         call interpCubic(TEOSgluondata, EOSgluonne, EOSgluone0,
     &       EOSgluonde, e, TEOSL7gluon)
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
     &     TEOSdata(EOSDATALENGTH)      
      double precision :: PEOSgluondata(EOSGLUONDATALENGTH),
     &                    SEOSgluondata(EOSGLUONDATALENGTH),
     &                    TEOSgluondata(EOSGLUONDATALENGTH)
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata, PEOSgluondata,
     & SEOSgluondata, TEOSgluondata

      double precision :: EOSe0         ! min (first) energy density
      double precision :: EOSde         ! energy density step
      double precision :: EOSEend       ! max (last) energy density
      Integer :: EOSne = EOSDATALENGTH ! number of table rows
      double precision :: EOSgluone0         ! min (first) energy density
      double precision :: EOSgluonde         ! energy density step
      double precision :: EOSgluonEend       ! max (last) energy density
      Integer :: EOSgluonne = EOSGLUONDATALENGTH  ! number of table rows
      common /EOSdatastructure/ EOSe0, EOSde, EOSEend,
     & EOSgluone0, EOSgluonde, EOSgluonEend, EOSne, EOSgluonne
      
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

      TEOSL7 = fugacity * TEOSL7qcd(e) + (1 - fugacity) * TEOSL7gluon(e)      
      
      return
      end

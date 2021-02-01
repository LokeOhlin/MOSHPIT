!!#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine rate_eq(neq, time, y, ydot, rpar, ipar)
      implicit NONE
#include "param.h"
#include "cool.h"
#include "non_eq.h"
#include "cons.h"
c#include "current_photo.h"
      REAL G_dust

      integer neq
c
      REAL time, y(nspec), ydot(nspec), rpar(nrpar)
      integer ipar(nipar)
c
      integer itemp, i
      REAL abundances(nabn), rates(nrates), 
     $     rates_chem(nrates_chem),
     $     y_in(nspec)
c
      REAL abe     , abhp    , abHI    , abh2    ,
     $     abcp    , abco    , abo
c
      REAL ylam    , temp    ,
     $     yn      , dl      , divv    , yntot   , gamma , 
     $     energy  , dtemp   , ncrinv  , h2var0  , h2var1, 
     $     h2var2  , ykdh    , ykdh2   , ykdis   , NH    ,
     $     ylam_chem, f_ex   , rho     , cs2     , L_jeans,
     $     f_tot   , RH2     , Rdust   , ekn
c
      REAL fshield_H2, fshield_CO, fshield_dust, chi_mean, 
     $     fshield_dust_H2
c
#if CHEMISTRYNETWORK == 5
      REAL k0, k1, gamma_chx, gamma_co, beta
#endif
c
      REAL phi     , h_gr    , hgrvar1 , hgrvar2
c
      REAL cr1, cr2, cr10, cr11, xr_ion, xr_heat, phih, hvphih, h2phih,
     $     h2vphih, Grad, uvphih, uvgeom, G0conv
c
      REAL ch1     , ch2     , ch3     , ch4     , ch5  ,
     $     ch6     , ch11    , ch13    , ch25    ,
     $     ch34    , ch35    , ch89
c
c
      REAL dtch1   , dtch2   , dtch3   , dtch4   , dtch5  ,
     $     dtch6   , dtch11  , dtch13  , dtch25  ,
     $     dtch34  , dtch35  , dtch89
c
      REAL cdex, rdex, fpheat
c
      REAL pi
c needs more digits I'd say
      parameter (pi = 3.14159265358979323846)
c used in photoelectric heating, found in bakes & tielens 1994
      parameter (G0conv = 1.6e-3)
c
      REAL heating_per_ion, pi_heat, rec_cool, tcool
      integer pi_equilib
c     LO: added parameter for fpump to calculate # of absorbed photons by H2
c     should probably be linked to the one defined in Chemistry_data.F90
      REAL fpump
      parameter (fpump = 6.94)
c     LO: dust energy density absorption
      REAL DustEabs 
c
c Unpack parameters
c

c      print*,'rate_eq=',y,ydot,rpar,ipar
      yn   = rpar(1)
      dl   = rpar(2)
      divv = rpar(3)
      fshield_H2   = rpar(4)
      fshield_CO   = rpar(5)
      fshield_dust = rpar(6)
      chi_mean     = rpar(7)
      fshield_dust_H2 = rpar(8)
      phih    = rpar(9)
      hvphih  = rpar(10)
      h2phih  = rpar(11)
      h2vphih = rpar(12)

! for UV pumping and photo electric heating
      uvphih  = rpar(13)
      uvgeom  = rpar(14)
! dust absorption
      DustEabs = rpar(15)
      
      pi_equilib = ipar(2)
      if (pi_equilib .eq. 1 .and. phih .gt. 0.0) then
        heating_per_ion = hvphih / phih
      endif

! DR: if no_chem (i.e. no reactions but still cooling), set all heating terms from chemical processes to zero.
      if (no_chem .eq. 1 .or. pi_equilib .eq. 1) then
        phih    = 0.0
        hvphih  = 0.0
        h2phih  = 0.0
        h2vphih = 0.0
        uvphih  = 0.0
      endif
c
c Chemical species:
c
      do I = 1, nspec
        y_in(I) = y(I)
      enddo
c
c      print*,'rate_eq:',y,y_in
      call validate_input(y, ipar)
      if (ipar(1) .eq. 1) then
c Return with zeroed ydot array if input abundances invalid
        do I = 1, nspec
          ydot(I) = 0e0
        enddo
        return
      endif
c
c Tracked, non-equilibrium species
c
      abhp = max(0e0, y(ihp))
      abh2 = max(0e0, y(ih2))
      abco  = 0e0
#if CHEMISTRYNETWORK == 5
      abco = max(0e0, y(ico))
#endif
c
c Non-equilibrium species computed using conservation laws
c
      abHI = max(1e0 - 2e0 * abh2 - abhp, 0e0)
#if CHEMISTRYNETWORK == 5
      abcp = max(0e0, abundc - abco)
      abe  = abhp + abcp + abundsi
#else
      abcp = abundc
      abe  = abhp + abundc + abundsi
#endif
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute temperature
c
c [NB We ignore corrections to yntot & gamma arising from low abundance 
c     molecular species (e.g. H-, H2+) and metals]
c
      yntot = (1e0 + abhe - abh2 + abe) * yn
      ekn   = energy / (yntot * kboltz)
c For time being, assume gamma = 5/3 -- good approx. unless we have lots of hot H2
      gamma = 1.6667 !****5e0 / 3e0
      temp = (2e0 / 3e0) * ekn
c
c     call lookup_temp_gamma(abh2, ekn, temp, gamma)
c
c Bounds checking
c
c Temp:
c
      if (temp .le. 0e0) then
        itemp   = 1
        dtemp   = 0e0
        temp    = tmin
        ipar(1) = 1
      elseif (temp .le. tmin) then
        itemp = 1
        dtemp = 0e0
      elseif (temp .ge. tmax) then
        itemp = nmd
        dtemp = 0e0
      else
        itemp = int(log10(temp) / dtlog) + 1
        if (itemp .le. 0 .or. itemp .gt. nmd) then
          print*, 'Fatal error in rate_eq.F', itemp, temp
          ABORT(1)
        endif
        dtemp = temp - temptab(itemp)
      endif
c
c Read in coefficients
c
      dtch1  = dtchtab(1, itemp)
      dtch2  = dtchtab(2, itemp)
      dtch3  = dtchtab(3, itemp)
      dtch4  = dtchtab(4, itemp)
      dtch5  = dtchtab(5, itemp)
      dtch6  = dtchtab(6, itemp)
      dtch11 = dtchtab(11, itemp)
      dtch13 = dtchtab(13, itemp)
      dtch25 = dtchtab(25, itemp)
      dtch34 = dtchtab(34, itemp)
      dtch35 = dtchtab(35, itemp)
      dtch89 = dtchtab(89, itemp)
c
      ch1  = chtab(1, itemp) + dtemp * dtch1
      ch2  = chtab(2, itemp) + dtemp * dtch2
      ch3  = chtab(3, itemp) + dtemp * dtch3
      ch4  = chtab(4, itemp) + dtemp * dtch4
      ch5  = chtab(5, itemp) + dtemp * dtch5
      ch6  = chtab(6, itemp) + dtemp * dtch6
      ch11 = chtab(11, itemp) + dtemp * dtch11
#ifdef TESTATOMONLY
      ch13 = 2.6e-13
#else 
      ch13 = chtab(13, itemp) + dtemp * dtch13
#endif
      ch25 = chtab(25, itemp) + dtemp * dtch25
      ch34 = chtab(34, itemp) + dtemp * dtch34
      ch35 = chtab(35, itemp) + dtemp * dtch35
      ch89 = chtab(89, itemp) + dtemp * dtch89
c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
      ncrinv   = (2e0 * abh2 * (ch6 - ch5) + ch5)
      h2var0   = 1e0 / ( 1e0 + yn * ncrinv)
      h2var1   = ch3**h2var0
      h2var2   = ch4**h2var0
      ykdh     = ch1 * h2var1
      ykdh2    = ch2 * h2var2
c
c Photoelectric effect:
c first term is from ISRF, second term is from radiation source
c conversion is done by taking all the photons arriving in the zone and then
c dividing by the total number of photons that would be arriving in a unit Habing field
c from all sides of the zone
      Grad = (UVgeom/ G0conv)
      G_dust  = G0 * fshield_dust + Grad
c
c Photodissociation -- diffuse field
c
      ykdis = G0 * phtab(1) * fshield_H2 * fshield_dust_H2

#if CHEMISTRYNETWORK == 5
      gamma_co = G0 * phtab(10) * fshield_dust * fshield_CO
      gamma_chx =  2.94e-10 * G_dust
#endif
c                                                                                                                                                                
c X-ray ionization and heating rates
c
      call calc_xray(NH_ext, abe, xr_ion, xr_heat)
      xr_ion  = xr_ion  * xray_scaling
      xr_heat = xr_heat * xray_scaling
c
c Rates for recombination on grain surfaces. These rates must be 
c computed here because of their dependence on the electron number
c density. 
c
      if (abe .eq. 0e0) then
c If the fractional ionization is zero, then there won't be any recombination,
c so the value we use for phi doesn't matter too much -- 1d20 is simply an
c arbitrary large number
c 
        phi = 1e20
      else
        phi = G_dust * sqrt(temp) / (yn * abe)    
      endif
c
c HII:
c
c Fit below only valid for phi > 100, rate is close to constant below this
c We therefore limit value of phi to be >= 100 to avoid over-estimating rate
c in gas with very low phi
c
      if (phi .lt. 1d2) then
        phi = 1d2
      endif

      hgrvar1  = 8.074e-6 * phi**1.378e0
      hgrvar2  = (1e0 + ch34 * phi**ch35)
      h_gr     = 1.225e-13 * dust_to_gas_ratio /
     $           (1e0 + hgrvar1 * hgrvar2)

c Exponentially reduce rate at high temperature to account for fact that we 
c expect hot gas to be largely dust free. This also prevents us from getting
c an anomalously large grain recombination cooling rate when T is large
c and G_dust is very small or zero
      h_gr = h_gr * dexp(-temp / 3d4)
c
c      print*,'rate_eq: yn=',yn

      call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust,
     $                   DustEabs, 0)
c
c Cosmic ray ionization rates: precomputed in cheminmo.
c
c HI:
      cr1  = crtab(1)  
      cr2  = crtab(2)   ! H2 -> H2+ + e-
      cr10 = crtab(10)  ! H2 -> H + H+ + e
      cr11 = crtab(11)  ! H2 -> H + H
c
c Skip rate calculation if chemistry switched off
c
      if (no_chem .eq. 1 .or. pi_equilib .eq. 1) then
        ydot(ihp) = 0d0
        ydot(ih2) = 0d0
#if CHEMISTRYNETWORK == 5
        ydot(ico) = 0d0
#endif
        ylam_chem = 0d0
      else
#ifdef TESTATOMONLY
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          -  ch13 * abhp  * abe) * yn     ! Grain H recomb.
     $          + phih  * abHI ! ionisation rate
      ydot(ih2) = 0
#else
c
c Protons
c
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          -  ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $          -  h_gr * abhp  * 1e0) * yn     ! Grain H recomb.
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          + xr_ion* abHI
     $          + cr10  * abh2
     $          + phih  * abHI ! ionisation rate
c
c Molecular hydrogen
c
      ydot(ih2) =   (-ch25  * abe  * abh2          ! Collisional dissoc.: H2 + e- 
     $          -   ykdh  * abHI * abh2          ! Collisional dissoc.: H2 + H 
     $          -   ykdh2 * abh2 * abh2          ! Collisional dissoc.: H2 + H2
     $          +   Rdust * abHI * 1e0 ) * yn    ! H2 formation on dust
     $          -   ykdis * abh2                 ! Photodissoc:         H2 + PH
     $      - 0.5e0 * cr2 * abh2   ! Fudge factor of two gives correct HI production
     $          -   cr10  * abh2
     $          -   cr11  * abh2
     $          -   h2phih* abh2 	 ! dissociation from ionising radiation
     $          -   uvphih* abh2 	 ! dissociation from stellar UV
#endif

#if CHEMISTRYNETWORK == 5
#ifdef TESTATOMONLY
      ydot(ico) = 0
#else
c
c CO -- from NL97
c
      k0 = 5e-16
      k1 = 5e-10
      abo = max(0.0, abundo - abco)
      if (abo .eq. 0.0) then
        beta = 0.0
      elseif (gamma_chx .lt. 1e-30) then
        beta = 1.0
      else 
        beta = k1 * abo / (k1 * abo + 
     $       gamma_chx / yn)
      endif
c the h2phih dissociation rate is a bit iffy it comes from ionising radiation and how that 
c interacts with CO is not modeled
c the uvphih is fine though, the factor comes from the comparison the dissociation rates 
c of H2 and CO at UV
      ydot(ico) = k0 * abcp * abh2 * beta * yn 
     $          - ch89 * abco * (abh2 + abhI + abhp) * yn  ! CO collisional dissociation
     $          - gamma_co * abco
     $          - h2phih*3.861003861 * abco
     $          - uvphih*3.861003861 * abco

#endif
c 3.8.. factor from ratio of CO to H2 dissociation rate in 1 draine field, 2.0e-10/5.18e-11, Roellig et al 2007 numbers

c Collisional dissociation rate is actually for the reaction CO + H -> C + OH
c
c We assume (a) rate for dissociation by H2, H+ similar; (b) in conditions where this reaction important, OH immediately
c dissociated and C immediately collisionally ionized. These are likely very crude approximations but ensure that we don't
c end up with large CO fractions in very hot (T >> 1e4 K) gas
c
#endif
c
c H2 collisional dissociation cooling:
c
#ifdef TESTATOMONLY
      rates_chem(1) = 0
#else
      rates_chem(1) = 4.48e0 * eV * (ykdh  * abHI * abh2 + 
     $                               ykdh2 * abh2 * abh2 +
     $                               ch25  * abe  * abh2 ) *yn**2.e0
#endif
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c CB: added rate from point source uvphih, which is 11.2 to 13.6 eV energy bin
#ifdef TESTATOMONLY
      rates_chem(2) = 0
#else

      rates_chem(2) = - 4e-1 * eV * abH2 * yn * ( ykdis + uvphih ) 
     $                - h2vphih * abh2 * yn
#endif
c
c Heating by UV pumping of H2 -- effective for n > n_crit
c
c We assume 2 eV per UV photon (Burton, Hollenbach & Tielens 1990) (Appendix A),
c and adopt a pumping rate 8.5 times larger than the dissociation rate
c (Draine & Bertoldi 1996).
c
#ifdef TESTATOMONLY
      rates_chem(3) = 0
#else
      cdex = (1.4 * exp(-1.81e4 / (temp + 1200)) * abH2
     $     +  1.0 * exp(-1e3 / temp) * abHI) * 1e-12 * sqrt(temp) * yn
      rdex = 2e-7
      fpheat = cdex / (cdex + rdex)

c CB: added rate from point source uvphih
      rates_chem(3) = - 2e0 * eV * abH2 * yn * ( ykdis + uvphih )*6.94e0
     $              * fpheat
#endif
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c We assume that a fraction f_ex of the energy goes into rotational &
c vibrational excitation of the H2, and a further fraction f_kin goes
c directly into kinetic energy (i.e. heat). The remaining energy goes
c directly into the dust 
c
c Dust -- 4.48eV 
c
      if (RH2 .gt. 0e0) then
c        call compute_h2_heating_fraction(temp, yn, abh2, f_ex)
c        f_ex = f_ex * h2_form_ex
        f_ex = 0e0
      else
        f_ex = 0e0
      endif
      f_tot = h2_form_kin + f_ex
c
#ifdef TESTATOMONLY
      rates_chem(4) = 0
#else
      rates_chem(4) = - 4.48e0 * eV * f_tot * RH2
#endif
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(5) = 13.6e0 * eV * ch11 * abe * abHI *yn**2.e0
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(6) = kboltz * temp * abhp * (ch13 * abe+ h_gr)*yn**2.e0
c
c CB:photoelectric heating from Bakes & Tielens 1994 eq. 42, 43
c    and Wolfire 2003, eq. 20
c

c photoelectric heating efficiency
c this is done in cool_func

c temp is dust temperature, rate is in /cm^2/s, G0conv is in /cm^2/s
c -1, so if the radiation field is attenuated it will be cut off
c radiatian was attenuated in the raytracing step so don't use G_dust
C      Grad = (UVgeom/G0conv-1.0)*sqrt(temp)/(abe*PAH*yn)
C      eps  = 4.9e-2/(1.0+4.e-3*Grad**0.73) + 3.7e-2*(temp/1e4)**0.7/(1.0+2.e-4*Grad)
C
C      rates_chem(7) = -1e-24*eps* (UVgeom/G0conv-1.0) *yn * 
C     $               dust_to_gas_ratio
c sign has to be negative, see ydot(itmp)
c      rates_chem(7) = 0.0

c
      ylam_chem = rates_chem(1)  + rates_chem(2)  + rates_chem(3)  +
     $            rates_chem(4)  + rates_chem(5)  + rates_chem(6)  
c     $            rates_chem(7)
      
      endif  ! if no_chem = 1, we skipped the ydot stuff above
c
c
      pi_heat  = 0.0
      rec_cool = 0.0
      if (pi_equilib .eq. 1) then
c Photoion. heating rate is heating per ion event times recombination rate
        pi_heat = heating_per_ion * (ch13 * abhp * abe + h_gr * abhp)
     $          * yn**2
c Also need to include recomb. cooling
        rec_cool = kboltz * temp * abhp * (ch13 * abe+ h_gr)*yn**2.e0
      endif
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0e0
      enddo
c
      abundances(1)  = abh2
#if CHEMISTRYNETWORK == 5
c      abundances(2)  = max(0e0, abundo - abco)
      abundances(2)  = max(0e0, abundo - abco) * max(0e0, 1.0-abhp)
      abundances(5)  = abco
      abundances(7)  = abcp
#else
c      abundances(2)  = abundo
      abundances(2)  = abundo * max(0e0, 1.0-abhp)
      abundances(7)  = abundc
#endif
      abundances(9)  = abundsi
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
c Approximate OII abundance
      abundances(17) = abundo * abhp
c
c Compute local Jeans length (for use in cool_func)
c
      rho = (1e0 + 4e0 * abhe) * mh * yn
      cs2 = (gamma * (gamma - 1e0) * energy / rho)
      L_jeans  = sqrt(pi * cs2 /
     $           (guniv * rho))
c
      call cool_func(temp, yn, L_jeans, divv, G_dust, abundances, ylam,
     $               rates, 0)
c

#ifdef TESTATOMONLY
      ylam=0.e0
#endif
c CB: hvphih term is photoionisation heating
      ydot(itmp) = -(ylam + ylam_chem) + (xr_heat + hvphih) * abHI * yn
     $           + pi_heat - rec_cool
c
c LO: amount of photons absorbed by H0 and H2
      ydot(iion) = phih*abHI*yn
      ydot(idis) = (h2phih+uvphih*(1+fpump))*abH2*yn
#ifdef COOLING_RATE_SANITY_CHECK
      if (ydot(itmp) .lt. 0.0) then
        tcool = abs(y(itmp) / ydot(itmp))
        if (tcool .lt. 1e4 / yn .and. yn .lt. 1e2) then
c Suspiciously short cooling time for this density
          print*, "WARNING: strange cooling time!"
          print*, "rates = ", rates
          print*, "rates_chem = ", rates_chem
          print*, "xrays = ", xr_heat * abhI * yn
          print*, "photoion = ", hvphih * abHI * yn
          print*, "PI (eqb) = " pi_heat
          print*, "REC (eqb) = ", rec_cool
          print*, "ylam = ", ylam
          print*, "ylam_chem = ", ylam_chem
          print*, "y = ", y
          print*, "rpar = ", rpar
        endif
      endif
#endif
      call validate_output(y,ydot,rpar)

c      print*,'rate_eq=',y,ydot,rpar,ipar
c
      return      
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              R A T E _ E Q                \\\\\\\\\\
c
c=======================================================================
!!#endif /* CHEMISTRYNETWORK == 4 || 5 */

c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////              C H E M I N M O              \\\\\\\\\\
c
c=======================================================================
c
       subroutine cheminmo
c
c    written by: S. Glover, 2004-2007
c    - based in part on code originally in coolinmo.F, 
c      written by G. Suttner, M. D. Smith & G. Pavlovski 
c    
c    PURPOSE: Precompute and tabulate chemical rate coefficients
c
c    NB. The temperature table, temptab, is initialized in coolinmo, 
c    which MUST be called before cheminmo.
c
c    REFERENCES:
c
c      DM56   -- Dalgarno & McDowell, 1956, Proc. Phys. Soc. London A, 69, 615
c      S59    -- Shavitt, 1959, J. Chem. Phys., 31, 1359
c      JGC67  -- Jacobs, Giedt & Cohen, 1967, J. Chem. Phys., 47, 54
c      L67    -- Lotz, 1967, ApJS, 14, 207
c      SA67   -- Schulz & Asundi, 1967, Phys. Rev., 158, 25
c      SFF67  -- Schmeltekopf, Fehsenfeld & Ferguson, 1967, ApJ, 118, L155
c      FA69   -- Fairbairn, 1969, Proc. Roy. Soc. London A, 312, 207
c      MAP70  -- Moseley, Aberth & Peterson, 1970, Phys. Rev. Lett., 24, 435
c      AP73   -- Aldrovandi & Pequignot, 1973, A&A, 25, 137 
c      MB73   -- MacGregor & Berry, 1973, J. Phys. B, 6, 181
c      AAT75  -- Azatyan, Aleksandrov & Troshin, 1975, Kinet. Catal., 16, 306
c      WK75   -- Walkauskas & Kaufman, 1975, Symp. Int. Combust. Proc., 15, 691
c      RP76   -- Ramaker & Peek, 1976, Phys. Rev. A, 13, 58 
c      SLA76  -- Slack, 1976, J. Chem. Phys., 64, 228
c      MR77   -- Mitchell & Le Roy, 1977, J. Chem. Phys., 67, 1042
c      FP78   -- Frommhold & Pickett, 1978, Chem. Phys., 28, 441
c      L78    -- Langer, 1978, ApJ, 225, 860
c      P78    -- Poulaert et al, 1978, J. Phys. B, 11, L671
c      CW79   -- Cohen & Westberg, 1979, J. Phys. Chem., 83, 46
c      HM79   -- Hollenbach and McKee, 1979, ApJS, 41, 555
c      W79    -- Wishart, 1979, MNRAS, 187, 59P
c      FAS80  -- Field, Adams & Smith, 1980, MNRAS, 192, 1
c      PH80   -- Prasad & Huntress, 1980, ApJS, 43, 1
c      G82    -- Gerlich, 1982, in 'Symp. on Atomic & Surface Phys.', Kluwer
c      GMV82  -- Gargaud, McCarroll & Valiron, 1982, A&A, 106, 197
c      H82    -- Huq et al, 1982, J. Phys. B, 15, 951
c      SvS82  -- Shull & van Steenburg, 1982, ApJS, 48, 95
c      AAS83  -- Alge, Adams & Smith, 1983, J. Phys. B., 16. 1433
c      CW83   -- Cohen & Westberg, 1983, J. Phys. Chem. Ref. Data, 12, 531
c      LS83   -- Lepp & Shull, 1983, ApJ, 270, 578
c      MD83   -- Mitchell & Deveau, 1983, ApJ, 266, 646
c      PSS83  -- Palla, Salpeter & Stahler, 1983, ApJ, 271, 632
c      ASM84  -- Adams, Smith & Millar, 1984, MNRAS, 211, 857
c      B84    -- Barlow, 1984, PhD thesis, Univ. Colorado
c      WAR84  -- Warnatz, 1984, in `Combustion Chemistry', (Springer-Verlag:NY), 197
c      H85    -- Herbst, 1985, ApJ, 291, 226
c      MR86   -- Murrell & Rodriguez 1986, J. Mol. Struct. Theochem., 139, 267
c      MS86   -- Mac Low & Shull, 1986, ApJ, 302, 585
c      TH86   -- Tsang & Hampson, 1986, J. Phys. Chem. Ref. Data, 15, 1087
c      DL87   -- Dalgarno & Lepp, 1987, in 'Astrochemistry' (IAU Symposium 120)
c      DRCM87 -- Dove, Rusk, Cribb & Martin, 1987, ApJ, 318, 379
c      J87    -- Janev et al, 1987, 'Elementary Processes in 
c                Hydrogen-Helium Plasmas', Springer
c      Or87   -- Orel, 1987, J. Chem. Phys., 87, 314
c      NR87   -- Natarajan & Roth, 1987, Combust. Flame, 70, 267
c      PDAN89 -- Petuchowski, Dwek, Allen & Nuth, 1989, ApJ, 342, 406
c      Z89    -- Zygelman et al, 1989, Phys. Rev. A, 40, 2340
c      DD90   -- Dalgarno, Du & You, 1990, ApJ, 349, 675
c      M90    -- Mitchell, 1990, Phys. Rep., 186, 215
c      DDH91  -- Dean, Davidson & Hanson, 1991, J. Phys. Chem., 95, 183 
c      L91    -- Launay et al, 1991, A&A, 252, 842
c      T91    -- Takagi et al, 1991, J. Phys. B, 24, 711
c      BAU92  -- Baulch et al, 1992, J. Phys. Chem. Ref. Data, 21, 411
c      F92    -- Ferland et al, 1992, ApJ, 387, 95
c      OLD92  -- Oldenborg et al 1992, J. Phys. Chem., 96, 8426
c      SMT92  -- Sidhu, Miller & Tennyson, 1992, A&A, 255, 453
c      HGS93  -- Harding, Guadagnini & Schatz, 1993, J. Phys. Chem., 97, 5472
c      K93    -- Kimura, Lane, Dalgarno & Dixson, 1993, ApJ, 405, 801
c      KD93   -- Kimura et al, 1993, ApJ, 417, 812
c      MLTS94 -- Mielke et. al, 1994, J. Phys. Chem., 98, 8000
c      PH94   -- Peart & Hayton, 1994, J. Phys. B, 27, 2551
c      S94    -- Schneider et al, 1994, ApJ, 424, 983
c      G95    -- Guberman, 1995, J. Phys. Chem, 102, 1699
c      LJB95  -- Linder, Janev & Botero, 1995, in 'Atomic & Molecular
c                Processes in Fusion Edge Plasmas', Plenum Press, 397
c      N95    -- Nahar, 1995, ApJS, 101, 423
c      SSS95  -- Stromholm et al, 1995, Phys. Rev. A, 52, R4320 
c      K96    -- Kimura et al, 1996, ApJ, 473, 1114
c      KF96   -- Kingdon & Ferland, 1996, ApJS, 106, 205
c      MSM96  -- Martin, Schwarz & Mandy, 1996, ApJ, 461, 265
c      N96    -- Nahar, 1996, ApJS, 106, 213 (errata for N95)
c      A97    -- Abel et al, 1997, New Astron, 2, 181
c      AS97   -- Andreazza & Singh, 1997, MNRAS, 287, 287 
c      NP97   -- Nahar & Pradhan, 1997, ApJS, 111, 339
c      V97    -- Voronov, 1997, ADNDT, 65, 1
c      GP98   -- Galli & Palla, 1998, A&A, 335, 403
c      HS98   -- Hummer & Storey, 1998, MNRAS, 297, 1073
c      L98    -- Larson et al, 1998, ApJ, 505, 459
c      MKM98  -- Martin, Keogh & Mandy, 1998, ApJ, 499, 793 
c      MMCV98 -- Mazzotta, Mazzitelli, Colafrancesco & Vittorio, A&AS, 133, 403
c      NAJ98  -- Naji et al, 1998, J. Phys. B, 31, 4887 
c      R98    -- Rosen et al, 1998, Phys. Rev. A, 57, 4462
c      SH98   -- Stancil et al, 1998, ApJ, 502, 1006
c      SLD98  -- Stancil, Lepp & Dalgarno, 1998, ApJ, 509, 1
c      CDG99  -- Croft, Dickinson & Gadea, 1999, MNRAS, 304, 327
c      N99    -- Nahar, 1999, ApJS, 120, 131
c      SSBO99 -- Singh et al, 1999, MNRAS, 303, 235
c      SS99   -- Stancil et al, 1999, A&AS, 140, 225
c      ST99   -- Stibbe & Tennyson, 1999, ApJ 513, L147
c      JBS00  -- Jensen, Bilodeau, Safvan et al, 2000, ApJ, 543, 764
c      N00    -- Nahar, 2000, ApJS, 126, 537
c      RDS00  -- Rosen, Derkatch, Semaniak et al, 2000, Far. Disc., 115, 295
c      STA00  -- Stancil et al, 2000, A&AS, 142, 107
c      G01    -- Glover, 2001, PhD Thesis
c      WD01   -- Weingartner & Draine, 2001, ApJ, 563, 842
c      ABN02  -- Abel, Bryan & Norman, 2002, Science, 295, 93
c      GP02   -- Galli & Palla, 2002, Planetary & Space Science, 50, 1197
c      LSD02  -- Lepp, Stancil & Dalgarno, 2002, J. Phys. B, 35, 57
c      S02    -- Savin, 2002, ApJ, 566, 599
c      TT02a  -- Trevisan & Tennyson, 2002a, Plasma Phys. Cont. Fus., 44, 1263
c      MPS03  -- Mielke et al, 2003, Phys. Rev. Lett., 91, 063201
c      MAC04  -- McCall et al, 2004, Phys. Rev. A, 70, 052716
c      SHC04  -- Smith, Herbst & Chang, 2004, MNRAS, 350, 323
c      SKHS04 -- Savin, Krstic, Haiman & Stancil, 2004, ApJ, 606, L167
c                [See also Erratum -- ApJ, 607, L147]  
c      WFP04  -- Walmsley, Flower & Pineau des Forets, 2004, A&A, 418, 1035
c      Z04    -- Zhao et al, 2004, ApJ, 615, 1063
c      GEP05  -- Geppert et al, 2005, J. Phys.: Conf. Ser., 4, 26
c      BvH06  -- Barinovs & van Hemert, 2006, ApJ, 636, 923
c      CGK06  -- Carty et al, 2006, J. Phys. Chem. A, 110, 3101
c      G07    -- Glover, 2007, in 'First Stars III' eds. A. Heger, B. O'Shea & T. Abel
c
      implicit NONE
#include "param.h"
#include "cool.h"
#include "shield_data.h"
c
      integer k, itemp, ifix, I, J
c
      REAL temp    , temp2    , tfix     ,
     $     t4log   , yncrh    , yncrh2   , dtemp , tinv  , 
     $     te      , temp_1ev , lnte     , temp3 , u
c
      REAL h2_low_n, h2_high_n, alpha    , f1       , f2
c
      REAL heii_rec_di, heii_rec_gs, heii_rec_A, heii_rec_B 
c
      REAL si2p_rec_rr, si2p_rec_dr
c
      REAL NCO, NH2, dN
c
c Temperature in Kelvin such that kb * temp_1ev = 1eV
      parameter (temp_1ev = eV / kboltz)
c
      do 200 itemp = 1, nmd
         temp  = temptab(itemp)
         temp2 = temp * 1e-2
         temp3 = temp / 3e2
         tinv  = 1.e0 / temp
         t4log = log10(temp) - 4e0
c Temp. in units of eV
         te    = temp / temp_1ev
         lnte  = log(te)
c
c 
c Collision-induced dissociation (CID) of H2 by HI
c
c The best available rates are those of MSM96. However, they're unwieldy 
c for use in numerical simulations, and it would be fiddly to interpolate
c from them to the H2-H2 rates (discussed below) as the H2 fraction 
c increases. However, MSM96 demonstrate that a reasonable approximation to
c their rate can be obtained by using the low density rate quoted by MS86,
c together with the high density rate and interpolation formula from LS83,
c _provided that the critical density is reduced by a factor of 10_.
c [For more details, see MSM96, section 4.2.1, esp. figure 13]
c
c This approximation should be good to within a factor of 2 or so over
c the temperature and density range of interest.
c
c NB. For cosmological sims. where three-body H2 formation is an important
c process, we need to make sure that the ratio of the three-body rate to
c the high density limit of the collisional dissociation rate is consistent 
c with the known equilibrium constant [see e.g. FH07]
c
c (ch1) -- k_hH -- high density limit
c
         if (iflag_3bh2a .eq. 1) then
c Inverse of ABN02 three-body rate
           if (temp .le. 3e2) then
             chtab(1,itemp) = 1.1e-9 * temp**0.135e0 
     $                      * exp(-5.2e4 / temp)
           else
             chtab(1,itemp) = (3.7e-8 / temp**0.485e0)
     $                      * exp(-5.2e4 / temp)
           endif
         elseif (iflag_3bh2a .eq. 2) then
c Inverse of PSS83 three-body rate
           chtab(1,itemp) = (5.24e-7 / temp**0.485e0) 
     $                    * exp(-5.2e4 / temp)
         elseif (iflag_3bh2a .eq. 3) then
c Inverse of CW83 rate
           chtab(1,itemp) = 8.4e-11 * temp**0.515e0 * exp(-5.2e4 /temp)
         elseif (iflag_3bh2a .eq. 4) then
c Rate from JGC67, used by FH07 to construct their three-body rate
           chtab(1,itemp) = (1.38e-4 / temp**1.025e0)
     $                    * exp(-5.2e4 / temp)
         elseif (iflag_3bh2a .eq. 5) then
c High density rate from MSM96
           if (temp .le. 4.5e4) then
           chtab(1,itemp) = 1e1**(-178.4239e0 - 68.42243e0 * log10(temp) 
     $                    + 43.20243e0 * log10(temp)**2e0
     $                    - 4.633167e0 * log10(temp)**3e0 
     $                    + 69.70086e0 * log10(1e0 + 40870.38e0 / temp)
     $                    - (23705.7 / temp))
           else
c Fit from MSM96 valid only up to 45000K; at higher T, assume constant rate
           tfix = 4.5e4
           chtab(1,itemp) = 1e1**(-178.4239e0 - 68.42243e0 * log10(tfix) 
     $                    + 43.20243e0 * log10(tfix)**2e0
     $                    - 4.633167e0 * log10(tfix)**3e0 
     $                    + 69.70086e0 * log10(1e0 + 40870.38e0 / tfix)
     $                    - (23705.7 / tfix))
           endif
         endif
c
c NB: We no longer use the LS83 rate, as it's a very poor fit low T, which
c can cause real problems in cosmological sims.
c
c         chtab(1,itemp) = 3.52d-9 * dexp(-4.39d4 / temp)
c
c (ch3) -- ratio of low density limit to high density limit
c
         alpha = 1e0 + 5.48e0 * temp_1eV / temp
         f1 = sqrt(temp) * exp(-alpha)
c
         alpha = 1e0 + 5.48e0 * temp_1eV / 4.5e3
         f2 = sqrt(4.5e3) * exp(-alpha)
c
         h2_low_n  = 1.2e-16 * f1 / f2
         h2_high_n = chtab(1,itemp)
c
c Avoid division-by-zero. Note that this doesn't introduce any significant
c inaccuracy, since the test is only true when the dissociation rate is 
c tiny.
c
         if (h2_high_n .lt. 1e-100) then
           chtab(3,itemp) = 1e-100
         else
           chtab(3,itemp) = max(h2_low_n / h2_high_n, 1e-100)
         endif
c 
c Collision-induced dissociation (CID) of H2 by H2
c
c Low n limit from MKM98; high n limit from detailed balance applied
c to three-body rate; interpolation formula from SK87.
c We could really do with better molecular data here, but this seems to be
c the best option at present.
c
c (ch2) -- k_hH2 -- high density limit
c
         if (iflag_3bh2b .eq. 1) then
c From PSS83
           chtab(2,itemp) = 6.5e-8 * exp(-5.2e4 / temp) / temp**0.485e0
         elseif (iflag_3bh2b .eq. 2) then
c From CW83
           chtab(2,itemp) = 2.67e-9 * exp(-5.2e4 / temp) / temp**0.085e0
         elseif (iflag_3bh2b .eq. 3) then
c From FH07
           chtab(2,itemp) = 1.71e-5 * exp(-5.2e4 / temp) / temp**1.025e0
         elseif (iflag_3bh2b .eq. 4) then
           chtab(2,itemp) = chtab(1,itemp) / 8e0
         endif
c
c (ch4) -- k_lH2 / k_hH2 (1.18d-10 / 1.3d-9)
c
         h2_high_n = chtab(2,itemp)
         h2_low_n  = 5.996e-30 * temp**4.1881e0 *exp(-54657.4e0 / temp)/
     $               (1e0 + 6.761e-06 * temp)**5.6881e0
c
         if (h2_high_n .lt. 1e-100) then
           chtab(4,itemp) = 1e-100
         else
           chtab(4,itemp) = max(h2_low_n / h2_high_n, 1e-100)
         endif
c
c Inverse critical densities for H2-H and H2-H2 collisional dissociation
c Taken from LS83 (H2-H rate; but reduced by factor of 10, in light of
c discussion above), SK87 (H2-H2 rate). 
c
c (ch5) -- 1 / ncr
c
c (ch6) -- 1 / ncrh2
c
         yncrh          = 1e1**( 3e0 - 0.416e0 * t4log
     $                          - 0.327e0 * t4log * t4log )
         yncrh2         = 1e1**( 4.845e0 - 1.3e0 * t4log
     $                          + 1.62e0 * t4log * t4log )
c
         chtab(5,itemp) = 1e0 / yncrh
c
         chtab(6,itemp) = 1e0 / yncrh2
c
c (ch7) -- from HM79
c
c  deff is a free parameter, representing the efficiency of 
c  H_2 reformation on dust. Remainder of HM79 expression is calculated
c  when we actually use the rate. [We can't pre-calculate it as it depends
c  on tdust, which may vary]
c
         chtab(7,itemp) = dust_to_gas_ratio * deff * 3e-17 *temp2**0.5e0
c
c He+ + H-  -> He + H
c
c Rate from SLD98, determined from cross-section of PH94
c
        if (temp .lt. 2e4) then
          chtab(8, itemp) = 2.32e-7 * temp3**(-0.52e0) * 
     $                      exp(temp / 2.24e4)
        else
          chtab(8, itemp) = 2.32e-7 * temp3**(-0.52e0) * 
     $                      exp(2e4 / 2.24e4)
        endif
c
c He+ + H2  -> He + H + H+
c
c Rate from B84. WFP04 cite a rate that is roughly a factor of two
c larger at all temperatures.
c
        if (temp .ge. 3e0) then
          chtab(9, itemp) = 3.7e-14 * exp(-3.5e1 / temp)
        else
          chtab(9, itemp) = 3.7e-14 * exp(-3.5e1 / 3e0)
        endif
c
c (ch10) -- see HeII recombination below
c
c Collisional ionization of HI (H + e => H+ + 2e) 
c From A97; based on data from J87
c
         chtab(11,itemp) = exp(-32.71396786e0 
     $                     + 13.5365560e0  * lnte 
     $                     - 5.73932875e0  * (lnte**2e0) 
     $                     + 1.56315498e0  * (lnte**3e0) 
     $                     - 0.28770560e0  * (lnte**4e0) 
     $                     + 3.48255977e-2 * (lnte**5e0)
     $                     - 2.63197617e-3 * (lnte**6e0) 
     $                     + 1.11954395e-4 * (lnte**7e0) 
     $                     - 2.03914985e-6 * (lnte**8e0))

c
c Recombination of HII  (H+ + e => H)
c From F92.
c
c Case A:
c
         chtab(12,itemp) = 1.269e-13 * (315614e0 * tinv)**1.503e0 /
     $                     ((1e0 + (604625e0 * tinv)**0.47e0)**1.923e0) 
c
c Case B:
c
         chtab(13,itemp) = 2.753e-14 * (315614e0 * tinv)**1.500e0 / 
     $                     ((1e0 + (115188e0 * tinv)**0.407e0)**2.242e0)

c      
c Collisional ionization of HeI
c From A97; based on data from J87
c
         chtab(14,itemp) = exp(-44.09864886e0 
     $                        + 23.91596563e0 * lnte 
     $                        - 10.7532302e0  * (lnte**2e0)
     $                        + 3.05803875e0  * (lnte**3e0) 
     $                        - 0.56851189e0  * (lnte**4e0) 
     $                        + 6.79539123e-2 * (lnte**5e0) 
     $                        - 5.00905610e-3 * (lnte**6e0) 
     $                        + 2.06723616e-4 * (lnte**7e0) 
     $                        - 3.64916141e-6 * (lnte**8e0))
c
c Collisional ionization of HeII
c From A97; based on data from J87
c
        chtab(15,itemp) = exp(-68.71040990e0 
     $                        + 43.93347633e0 * lnte
     $                        - 18.4806699e0  * (lnte**2e0)
     $                        + 4.70162649e0  * (lnte**3e0)
     $                        - 0.76924663e0  * (lnte**4e0)
     $                        + 8.113042e-2   * (lnte**5e0)
     $                        - 5.32402063e-3 * (lnte**6e0)
     $                        + 1.97570531e-4 * (lnte**7e0)
     $                        - 3.16558106e-6 * (lnte**8e0))
c
c Recombination of He II
c Dielectronic: (from AP73)
c
        heii_rec_di = 1.9e-3 * exp(-473421e0 / temp) 
     $              * (1e0 + 0.3e0 * exp(-94684e0 / temp))
     $              / temp**1.5e0
c
c Dielectronic recombination of HeII at high T can be a significant source
c of cooling, so we need to store this part of the rate, as well as the total
c rate below.
c
        chtab(10, itemp) = heii_rec_di
c
c Case A radiative: (from HS98)
c
        heii_rec_A  = 1e-11 * (12.72e0 - 1.615e0*log10(temp) 
     $              - 0.3162e0*(log10(temp)**2e0) 
     $              + 4.930e-2*(log10(temp)**3e0)) / sqrt(temp)
c
c Case B radiative: (from HS98)
c
        heii_rec_B  = 1e-11 * (11.19e0 - 1.676e0*log10(temp) 
     $              - 0.2852e0*(log10(temp)**2e0) 
     $              + 4.433e-2*(log10(temp)**3e0)) / sqrt(temp)
c
c Rec. direct to ground state (i.e. case A - case B)
c
        heii_rec_gs = heii_rec_A - heii_rec_B
c
c Optically thin rate
c
        chtab(16,itemp) = heii_rec_di + heii_rec_A
c
c Optically thick rate
c
        chtab(17,itemp) = 0.677e0 * heii_rec_gs + heii_rec_B 
     $                  + heii_rec_di
c
c HI ionization due to HeII recombination emission.
c (NB This assumes gas is optically thick to ionizing radiation)
c
        chtab(18,itemp) = 0.96e0 * chtab(17,itemp) + 0.677e0 * 
     $                    heii_rec_gs
c
c (NB For the derivation of these formulae, and a detailed discussion,
c     see Osterbrock, 1989, 'Astrophysics of Gaseous Nebulae and AGN')
c
c Recombination of HeIII -- appropriately scaled versions of HII rates
c Case A:
c
        chtab(19,itemp) = 2e0 * 1.269e-13 
     $                  * (4e0 * 315614e0 * tinv)**1.503e0 
     $                  / (1e0 + (4e0*604625e0 * tinv)**0.47e0)**1.923e0
c
c Case B:
c
        chtab(20,itemp) = 2e0 * 2.753e-14
     $                 * (4e0 * 315614e0 * tinv)**1.5e0
     $                 / (1e0 + (4e0*115188e0 * tinv)**0.407e0)**2.242e0
c
c We don't include HI ionization by HeIII recombination emission, since
c if we have significant HeIII, it's hard to see how we'll have any HI 
c left. 
c
c We also don't include HeI ionization by HeIII recombination emission,
c which means that we might not get the ionization balance quite right
c in situations where this is significant
c

c
c Formation of H- by photo-attachment (H + e => H-)
c Fit by A97, data from W79
c
        if (temp .le. 6000) then
          chtab(21,itemp) = 10e0**(-17.845e0 
     $                      + 0.762e0  *  log10(temp) 
     $                      + 0.1523e0 *  log10(temp)**2e0 
     $                      - 3.274e-2 *  log10(temp)**3e0 ) 
        else
          chtab(21,itemp) = 10e0**(-16.4199e0 
     $                      + 0.1998e0  * log10(temp)**2e0
     $                      - 5.447e-3  * log10(temp)**4e0
     $                      + 4.0415e-5 * log10(temp)**6e0 ) 
        endif
c
c Formation of H_2 by associative detachment of H- (H + H- => H_2 + e)
c
        if (iflag_ad .eq. 1) then
c Fit by GP98, data from L91
          if (temp .le. 3e2) then
            chtab(22,itemp) = 1.5e-9
          else
            chtab(22,itemp) = 4.0e-9 / temp**0.17e0
          endif
        elseif (iflag_ad .eq. 2) then
c Largest theoretical value -- from Daniel Savin
          chtab(22,itemp) = 5.0e-9
        elseif (iflag_ad .eq. 3) then
c Rate from PH80, based on experimental measurement of SFF67.
          chtab(22,itemp) = 1.3e-9
        elseif (iflag_ad .eq. 4) then
c Half measured value -- from Daniel Savin; true value should lie 
c somewhere between this and value given by iflag_ad = 2
          chtab(22,itemp) = 0.65e-9
        else
          print*, "Unknown value specified for iflag_ad:", iflag_ad
          ABORT(426)
        endif
c
c H_2+ formation from HI (H + H+ => H_2+)
c Fit by GP98, data from RP76 
c
        chtab(23,itemp) = 10e0**(-19.38e0 
     $                    - 1.523e0 * log10(temp) 
     $                    + 1.118e0 * log10(temp)**2e0 
     $                    - 0.1269e0 *log10(temp)**3e0 )
c
c No point in tabulating rate for H2+ + H => H2 + H+, as it is
c temperature independent. For reference, the rate is 6.4d-10
c

c
c Destruction of H_2 by charge transfer (H_2 + H+ => H_2+ + H)
c For T < 3e4 K, we use a rate that is a fit to the calculations 
c of SKHS04. At higher temperature, the SKHS04 fit is not valid 
c and can potentially cause numerical problems. At these temperatures,
c we therefore assume a rate that is a constant prefactor * the
c exp(-1.83eV/T) exponential dependence, with the prefactor chosen
c such that the rate remains continuous at 3e4 K.
c
        if (temp .le. 3e4) then
          chtab(24,itemp) = (- 3.3232183e-7
     $                       + 3.3735382e-7  * log(temp)
     $                       - 1.4491368e-7  * log(temp)**2e0
     $                       + 3.4172805e-8  * log(temp)**3e0
     $                       - 4.7813720e-9  * log(temp)**4e0
     $                       + 3.9731542e-10 * log(temp)**5e0
     $                       - 1.8171411e-11 * log(temp)**6e0
     $                       + 3.5311932e-13 * log(temp)**7e0 )
     $                       * exp(-2.1237150e4 / temp)
       else
         chtab(24,itemp) = 3.7289e-10 * exp(-2.1237150e4 / temp)
       endif
c
c CID of H_2 by electrons (H_2 + e => 2H + e)  
c
c v=0 rate from TT02a
c
        chtab(25,itemp) = 4.4886e-9 * temp**0.109127e0 
     $                    * exp(-1.01858e5 / temp)
c
c Destruction of H- by electrons (H- + e => H + 2e)
c Fit by A97, data from J87
c 
        chtab(26,itemp) = exp(-18.01849334e0 
     $                        + 2.3608522e0   * lnte 
     $                        - 0.28274430e0  * (lnte**2e0)
     $                        + 1.62331664e-2 * (lnte**3e0)
     $                        - 3.36501203e-2 * (lnte**4e0)
     $                        + 1.17832978e-2 * (lnte**5e0)
     $                        - 1.65619470e-3 * (lnte**6e0)
     $                        + 1.06827520e-4 * (lnte**7e0)
     $                        - 2.63128581e-6 * (lnte**8e0))
c
c H- destruction by HI (H- + H => 2H + e)
c Fit by A97, data from J87
c
        if (te .le. 1e-1) then
          chtab(27,itemp) = 2.5634e-9 * Te**1.78186e0
        else
          chtab(27,itemp) = exp(-20.37260896e0 
     $                          + 1.13944933e0 * lnte 
     $                          - 0.14210135e0 * (lnte**2e0)
     $                          + 8.4644554e-3 * (lnte**3e0)
     $                          - 1.4327641e-3 * (lnte**4e0)
     $                          + 2.0122503e-4 * (lnte**5e0)
     $                          + 8.6639632e-5 * (lnte**6e0)
     $                          - 2.5850097e-5 * (lnte**7e0)
     $                          + 2.4555012e-6 * (lnte**8e0)
     $                          - 8.0683825e-8 * (lnte**9e0))
        endif
c
c Mutual neutralization (H- +  H+ => 2H)
c
        if (iflag_mn .eq. 1) then
c From GP98, based on data from MAP70.  
          chtab(28,itemp) = (5.7e-6 / sqrt(temp)) + 6.3e-8
     $                     - 9.2e-11 * sqrt(temp) + 4.4e-13 * temp
        elseif (iflag_mn .eq. 2) then
c From A97, based on DL87
          chtab(28,itemp) = 7e-7 / sqrt(temp)
        elseif (iflag_mn .eq. 3) then
c From PH80
          chtab(28,itemp) = 4e-6 / sqrt(temp)
        elseif (iflag_mn .eq. 4) then
c From CDG99
          chtab(28,itemp) = 2.4e-6 * (1e0 + temp / 2e4) / sqrt(temp)
        else
          print*, "Unknown value specified for iflag_mn:", iflag_mn
          ABORT(427)
        endif
c
c H_2+ formation from H- (H+ + H- => H_2+ + e)
c Fit by GP98 to data by P78
c
        if (temp .le. 8e3) then
          chtab(29,itemp) = 6.9e-9 / temp**0.35e0
        else
          chtab(29,itemp) = 9.6e-7 / temp**0.9e0
        endif
c
c Dissociative recombination of H_2+ (H2+ + e => 2H)
c Fit by A97 to data from S94. Bogus below 100K.
c
c The LTE rate given by S94, assuming T_e == T_n, is
c significantly smaller: k ~ 1.2e-8 / T^0.64 at T > 100K.
c
c However, the OSU2005 rate compilation cites a measurement
c by McGowan et al, 1979 Phys. Rev. Lett, 42. 373 which is
c much larger -- their rate is k = 2.53e-7 (T/300)**(-0.5), 
c which is 25x larger than our adopted value at 300K. This
c is for H2+ with significant vibrational, rotational 
c excitation, but suggests that our adopted rate may be too 
c low.
c
        if (temp .le. 6.17e2) then
          chtab(30,itemp) = 1e-8
        else
          chtab(30,itemp) = 1.32e-6 / temp**0.76e0
        endif
c
c H- electron detachment by He (H- + He -> H + He + e)
c From H82.
c 
        chtab(31,itemp) = 4.1e-17 * temp**2e0 * exp(-19870e0 / temp)
c
c Charge transfer, helium to hydrogen: (He+ + H -> He + H+)
c Fit by SLD98 to data by Z89.
c
c We assume that CT with DI has same rate as with HI
c
        chtab(32,itemp) = 1.2e-15 * (temp/3e2)**0.25e0  
c
c Charge transfer, hydrogen to helium: (He + H+ -> He+ + H)
c Above 10^4 K, fit by GP98 to data from K93.
c Below 10^4 K, fit by G01 to same data source. 
c
c We assume that CT with DI has same rate as with HI
c
        if (temp .le. 1e4) then
          chtab(33,itemp) = 1.26e-9 * exp(-1.275e5 / temp) /
     $                      temp**0.75e0
        else
          chtab(33,itemp) = 4e-37 * temp**4.74e0
        endif
c
c Grain assisted recombination. Rate from WD01.
c WD01 quote an accuracy of 20% for 10 < T < 10^3 K, 10^2 < phi < 10^6
c [where phi = G0 * sqrt(T) / ne]
c
c The WD01 fitting formula is:
c
c  alpha_g = 1e-14 * C0 * (1 + C1 * phi^C2 * (1 + C3 * T^C4 * phi * 
c            (- C5 - C6 ln T)))^-1
c
c We tabulate (C3*T^C4) and (- C5 - C6 ln T) for each ion of interest
c (currently HII, HeII, CII, SiII). The other portions of the fit need to be 
c computed on the fly in the chemistry routine, as they depend on the 
c electron density.
c
c The grain recombination rate for OII can be obtained from that for
c HII by a simple rescaling and so we need not tabulate it here.
c
c HII:
c
        chtab(34,itemp) = 5.087e2 * temp**1.586e-2
c
        chtab(35,itemp) = - 0.4723e0 - 1.102e-5 * log(temp)
c
c CII:
c
        chtab(36,itemp) = 4.331e2 * temp**4.845e-2
c
        chtab(37,itemp) = - 0.8120e0 - 1.333e-4 * log(temp)
c
c SiII:
c
        chtab(38,itemp) = 4.375e4 * temp**1.635e-6
c
        chtab(39,itemp) = - 0.8964e0 - 7.538e-5 * log(temp)
c
c HeII:
c
        chtab(135,itemp) = 5.115e3 * temp**3.903e-7
c
        chtab(136,itemp) = - 0.4956e0 - 5.494e-7 * log(temp)
c
c CII gas-phase recombination (radiative plus dielectronic) [case A rate]
c
c Based on data from NP97, fitting formulae from UMIST database.
c 
c Strictly speaking, the UMIST fit is not valid below 10K, but in 
c practice it is well behaved, so extrapolating it down to 1K 
c shouldn't cause too much of a problem. The fit also starts to 
c become inaccurate above 10^6 K, but this is well above the
c temperature range that we're interested in.
c
        if (temp .le. 7.95e3) then
          chtab(40,itemp) = 4.67e-12 * temp3**(-0.6e0)
        elseif (temp .le. 2.114e4) then
          chtab(40,itemp) = 1.23e-17 * temp3**2.49e0   *
     $                      exp(2.18456e4 / temp)
        else 
          chtab(40,itemp) = 9.62e-8 * temp3**(-1.37e0) *
     $                      exp(-1.157862e5 / temp)
        endif
c
c SiII gas-phase recombination (radiative plus dielectronic) [case A rate]
c
c Rate from N00; broken power-law fit by SCOG. Fit is accurate to within 
c about 5% for T < 10^4K; for 10^4 < T < 7.9 x 10^4K, fit accuracy 
c degrades to 25%, while at higher temperatures, rate is accurate only 
c to within a factor of 2. N00 estimates accuracy of calculated rates to 
c be 10-30%, so these fits should be adequate for the temperature range that
c we're interested in. [NB N00's tabulated data starts at 10K;
c below this, our rate is extrapolated, but it remains well-behaved]
c
c [NB The UMIST rate is 25-30% lower than the N00 rate at T < 1000K,
c  with the discrepancy increasing significantly at higher T -- it is
c  almost 50 times smaller at temps. of order 10^5K -- presumably, this
c  is because it does not include dielectronic recombination]
c
        if (temp .le. 2e3) then
          chtab(41,itemp) = 7.50e-12 / temp3**0.55e0
        elseif (temp .le. 1e4) then
          chtab(41,itemp) = 4.86e-12 / temp3**0.32e0
        elseif (temp .le. 7.9e4) then
          chtab(41,itemp) = 9.08e-14 * temp3**0.818e0
        else
          chtab(41,itemp) = 5.89e-9 / temp3**1.17e0
        endif
c 
c OII gas-phase recombination (radiative plus dielectronic) [case A rate]
c
c Rate from N99; broken power-law fit by SCOG (although form
c of high temperature portion is based on a fit from SvS82). 
c Fit is accurate to within about 10% for T < 10^7 K; at higher temperatures, 
c accuracy degrades, but at that point there's essentially no OII anyway, so
c we don't really care 
c
        if (temp .le. 4e2) then
          chtab(42,itemp) = 1.30e-10 / temp**0.64e0
        else
          chtab(42,itemp) = (1.41e-10 / temp**0.66e0) +
     $                      7.4e-4 * exp(-1.75e5 / temp) *
     $                      (1e0 + 6.2e-2 * exp(-1.45e5 / temp)) /
     $                      temp**1.5e0
        endif
c
c Collisional ionization: all rates from V97.
c
c CI:
c
        u = 11.26e0 / te
c
        if (u .gt. 1e2) then
          chtab(43,itemp) = 0e0
        else
          chtab(43,itemp) = (6.85e-8 / (0.193e0 + u)) * u**0.25e0 
     $                      * exp(-u)
        endif
c
c SiI:
c
        u = 8.2e0 / te
c
        if (u .gt. 1e2) then
          chtab(44,itemp) = 0e0
        else
          chtab(44,itemp) = (1.88e-7 * (1e0 + sqrt(u)) / (0.376e0 + u))
     $                      * u**0.25e0 * exp(-u)
        endif
c
c OI:
c
        u = 13.6e0 / te
c
        if (u .gt. 1e2) then
          chtab(45,itemp) = 0e0
        else
          chtab(45,itemp) = (3.59e-8 / (0.073e0 + u)) * u**0.34e0 
     $                      * exp(-u)
        endif
c
c Charge transfer:
c
c OII + HI  -> OI  + HII  -- from SS99; total for all three J states of OI
c
        chtab(46,itemp) = 2.08e-9  * (temp / 1e4)**0.405e0 +
     $                    1.11e-11 * (temp / 1e4)**(-0.458e0)
c
c OI  + HII -> OII + HI  -- from SS99.
c
c Assumes LTE level populations; note that depopulating the upper levels 
c relative to LTE will alter the rate coefficient by 50%, at most -- the
c rate is significantly altered only if the upper levels can be 
c overpopulated at low T (i.e. at temperatures where the 223K endothermicity
c of the CT reaction is a significant barrier).
c
        chtab(47,itemp) = (1.26e-9  * (temp / 1e4)**0.517e0  +
     $                     4.25e-10 * (temp / 1e4)**6.69e-3) * 
     $                     exp(-227e0 / temp)
c
c CI  + HII -> CII + HI -- from SH98. 
c
c Reverse reaction is given by ch58 below.
c
        if (temp .le. 1e5) then
          chtab(48,itemp) = 2.80e-15 * (temp / 1e4)**0.213e0
        else
          chtab(48,itemp) = 2.80e-15 * (temp / 1e4)**0.213e0 +
     $                      1.28e-15 * (temp / 1e4)**3.08e0 
     $                      * exp(-2.19e6 / temp)
        endif
c
c SiI + HII -> SiII + HI
c
c Produces SiII ions in the Si+(4P) and Si+(2D) states, which lie
c significantly above the ground state. Since these states are only
c populated at high gas density and high T, the rate of the reverse 
c process will generally be negligible.
c
c Quoted rate is my fit to data from K96, valid for 200 < T < 2e5 K.
c Fit is accurate to within 10% over this temperature range. 
c
        if (temp .le. 1e4) then
          chtab(49,itemp) = 5.88e-13 * temp**0.848e0
        else
          chtab(49,itemp) = 1.45e-13 * temp
        endif
c
c O + H  -> OH + gamma -- rate from UMIST99
c
        chtab(50,itemp) = 9.90e-19 / temp3**0.38e0
c
c H + OH -> H2O + gamma -- rate from UMIST99, based on FAS80
c
        chtab(51,itemp) = (5.26e-18 / temp3**5.22e0) * exp(-9e1 / temp)
c
c O + H2 -> OH + H -- rate from UMIST99, based on NR87
c
        chtab(52,itemp) = 3.14e-13 * temp3**2.7e0 * exp(-3.15e3 / temp)
c
c H + OH  -> O + H2 -- rate from UMIST99, based on TH86
c
        chtab(53,itemp) = 6.99e-14 * temp3**2.8e0 * exp(-1.95e3 / temp)
c
c H2 + OH -> H2O + H -- rate from UMIST99, based on OLD92
c
        chtab(54,itemp) = 2.05e-12 * temp3**1.52e0 * exp(-1.736e3 /temp)
c
c H + H2O -> OH + H2 -- rate from UMIST99, based on CW79
c
        chtab(55,itemp) = 1.59e-11 * temp3**1.2e0 * exp(-9.61e3 / temp)
c
c Dissociative recombination of HOC+:  HOC+ + e -> CO + H
c
c Rate from UMIST99, who assume that CO + H is the only outcome.
c GEP05 mention that it may also produce OH + C, but give no rate
c for HOC+ DR (they're concerned only with HCO+ DR).
c
c cf. HCO+ dissociative recombination (ch87).
c
        chtab(56, itemp) = 1.1e-7 / temp3
c
c CO formation by radiative association:   C + O -> CO + gamma
c
c Rate coefficient from UMIST99, based on SSBO99
c
        if (temp .le. 3e2) then
          chtab(57,itemp) = 2.1e-19
        else
          chtab(57,itemp) = 3.09e-17 * temp3**0.33e0 * 
     $                      exp(-1.629e3 / temp)
        endif
c
c C+ + H -> C + H+
c
c Rate coefficient from SH98, with typo corrected. Claimed to be
c accurate to within 50%.
c
        chtab(58,itemp) = 6.08e-14 * (temp / 1e4)**1.96e0 * 
     $                    exp(-1.7e5 / temp)
c
c CH2+ formation by radiative association: C+ + H2 -> CH2+ + gamma
c
c Rate from H85
c
        chtab(59, itemp) = 4.0e-16 / temp3**0.2e0
c
c CH+ formation from endothermic reaction between C+ and H2: C+ + H2 -> CH+ + H
c
c Rate from ASM84.
c
        chtab(60, itemp) = 1.0e-10 * exp(-4.64e3 / temp)
c
c CH+ dissociative recombination: CH+ + e- -> C + H
c 
c Rate from OSU2005; based on T91
c
        chtab(61, itemp) = 7.0e-8 / temp3**0.5e0
c
c CH2+ destruction by reaction with H: CH2+ + H -> CH+ + H2
c
c Rate from UMIST99
c
        chtab(62, itemp) = 1.0e-9  * exp(-7.08e3 / temp) 
c
c CH2+ dissociative recombination
c
c From UMIST06, based on L98.
c
c XXX: Now that T dependence is the same, we could merge these rates
c
        chtab(63, itemp) = 1.60e-7 / temp3**0.60e0 ! CH2+ + e- -> CH + H 
        chtab(64, itemp) = 4.03e-7 / temp3**0.60e0 ! CH2+ + e- -> C + H + H
        chtab(65, itemp) = 7.68e-8 / temp3**0.60e0 ! CH2+ + e- -> C + H2
c
c CH3+ destruction by atomic H: CH3+ + H -> CH2+ + H2
c
c Rate from UMIST99
c
        chtab(66, itemp) = 7e-10 * exp(-1.056e4 / temp)
c
c CH3+ dissociative recombination:
c
c Rates for ch67, ch68 are from UMIST99, based on M90; rate for ch69 is 
c from UMIST99 with no other source cited
c 
        chtab(67, itemp) = 7.75e-8 / sqrt(temp3)  ! CH3+ + e -> CH2 + H
        chtab(68, itemp) = 1.95e-7 / sqrt(temp3)  ! CH3+ + e -> CH + H2
        chtab(69, itemp) = 2.00e-7 / temp3**0.4e0  ! CH3+ + e -> CH + 2H
c
c CH2 formation from CH:  CH + H2 -> CH2 + H 
c
c From UMIST99, based on data from the NIST chemical kinetics database
c (exact source or sources unclear)
c
        chtab(70, itemp) = 5.46e-10 * exp(-1.943e3 / temp)
c
c H3O+ dissociative recombination, total rate; branching ratios in rate_eq
c
c Rate from UMIST06, based on JBS00.
c
        chtab(71, itemp) = 4.318e-7 / sqrt(temp3)
c
c OH+ dissociative recombination:  OH+ + e- -> O + H
c
c From OSU2005, based on G95.
c
        chtab(72, itemp) = 6.3e-9 * temp3**(-0.48e0)
c
c H2O+ dissociative recombination, total rate; branching ratios in rate_eq
c
c Rates from UMIST06, based on RDS00
c
        chtab(73, itemp) = 4.3e-7 / sqrt(temp3)
c
c SiH+ dissociative recombination: SiH+ + e -> Si + H
c
c Rate from UMIST99
c
        chtab(74, itemp) = 2e-7 / sqrt(temp3)
c
c SiH2+ dissociative recombination, total rate; branching ratios in rate_eq
c
c Currently not used.
c
        chtab(75, itemp) = 5e-7 / sqrt(temp3)
c
c ch76 -- formerly collisional dissoc. of OH by HI; currently unused
c
        chtab(76, itemp) = 0e0
c
c CO formation from CH:  CH + O -> CO + H
c
c From NIST kinetics database, based on MR86 at high T 
c and BAU92 at low T
c
        if (temp .le. 2e3) then
          chtab(77, itemp) = 6.6e-11
        else
          chtab(77, itemp) = 1.02e-10 * exp(-914e0 / temp)
        endif
c
c CH formation from C and H2:  C + H2 -> CH + H
c
c Rate from DDH91
c
        chtab(78, itemp) = 6.64e-10 * exp(-1.17e4 / temp)
c
c C2 formation by RA of carbon: C + C -> C2 + gamma
c
c Rate from AS97
c
        chtab(79, itemp) = 4.36e-18 * temp3**0.35e0 * 
     $                     exp(-1.613e2 / temp)
c
c CH destruction by atomic H: CH + H -> C + H2
c
c Rate from HGS93
c
      chtab(80, itemp) = 1.31e-10 * exp(-8e1 / temp)
c
c O2 formation from OH:   O + OH -> O2 + H
c
c Low T rate from CGK06, based on measurements between 39K and 142K. We extend
c it to a higher temperature of 261K so that it joins up smoothly with the high
c T rate.
c
c High T rate is from UMIST99, based on data from the NIST chemical kinetics
c database (the exact sources are unclear). 
c
        if (temp .le. 2.61e2) then
          chtab(81, itemp) = 3.5e-11 
        else
          chtab(81, itemp) = 1.77e-11 * exp(1.78e2 / temp)
        endif
c
c O2 formation by RA of oxygen:  O + O -> O2 + gamma
c
c Rate from PH80
c
        chtab(82, itemp) = 4.9e-20 * temp3**1.58e0
c
c CO formation from C and O2:  C + O2 -> CO + O
c
c Rate from UMIST06. Low T portion from SHC04. High T portion based
c on NIST kinetics database (exact sources unclear).
c
        if (temp .le. 2.95e2) then
          chtab(83, itemp) = 4.7e-11 * temp3**(-0.34e0)
        else
          chtab(83, itemp) = 2.48e-12 * temp3**1.54e0 * 
     $                       exp(6.13e2 / temp)
        endif
c
c O2+ dissociative recombination:  O2+ + e- -> O + O
c
c Rate from UMIST99, based on AAS83.
c
        chtab(84, itemp) = 1.95e-7 / temp3**0.7e0
c
c O2 destruction by atomic hydrogen: O2 + H -> OH + O
c
c Rate from UMIST99, based on NIST chemical kinetics database
c (exact sources unclear)
c
        chtab(85, itemp) = 2.61e-10 * exp(-8.156e3 / temp)
c
c O2 destruction by H2:  O2 + H2 -> OH + OH
c
c Rate from AAT75.
c
        chtab(86, itemp) = 3.16e-10 * exp(-2.189e4 / temp)
c
c HCO+ dissociative recombination: HCO+ + e -> products.
c 
c Rate from GEP05. Branching ratios are given in rate_eq().
c 
        chtab(87, itemp) = 3.0e-7 / temp3**0.64e0
c
c Mutual neutralization, generic rate.
c Taken from UMIST99.
c Used for:
c
c   1) H+   + C- -> C + H
c   2) H+   + O- -> O  + H
c   3) Mg+  + H- -> Mg + H  (currently unused)
c
c Mutual neutralzation of H- with +ve ions such as O+, C+,
c HCO+, H3O+ etc. assumed negligible compared to MN with 
c H+ or He+ & compared to DR of the +ve ions in question.
c
        chtab(88, itemp) = 2.3e-7 / sqrt(temp3)
c
c CO destruction by atomic hydrogen:  CO + H -> C + OH
c 
c Rate from UMIST99
c
        chtab(89, itemp) = 1.1e-10 * sqrt(temp3) * 
     $                     exp(-7.77e4 / temp)
c
c CO formation from C2 and O:  C2 + O -> CO + C
c
c UMIST99 lists fit for rate between 10 < T < 300K (based on
c MD83) and a single value at 8000K (based on FA69). The high T
c fit below is desgined to smoothly interpolate between these
c two rates and may be highly inaccurate at intermediate 
c temperatures
c
        if (temp .le. 3e2) then
          chtab(90, itemp) = 5e-11 * sqrt(temp3)
        else
          chtab(90, itemp) = 5e-11 * temp3**0.757e0
        endif
c
c HD+ formation by radiative association: H + D+ -> HD+ + gamma
c                                         D + H+ -> HD+ + gamma
c
c SLD98 give same rate for both processes -- functional form is same
c as their rate for the H + H+ analog, but scaling for each rate
c is slightly greater than 50% of the H + H+ rate. They cite RP76,
c FP78
c
c NB SLD98 don't quote a temperature range of validity for their fit,
c but clearly it is not going to work all the way down to 1K, thanks 
c to the exponential term; comparison of their rate for the H + H+
c reaction with that from GP98 suggests that the fit is good for 
c T > 7K, and for lower T, we assume a constant value. Note that 
c in cosmological sims., these temperatures are << T_CMB, so this 
c modification shouldn't have any effect.
c
c WFP04 have a similar rate, but don't have the small exponential 
c term; at T >> 20K, the rates are essentially identical
c
        if (temp .gt. 7e0) then
          chtab(91,itemp) = 3.9e-19 * (temp / 3e2)**1.8e0 * 
     $                      exp(2e1 / temp)
        else
          chtab(91,itemp) = 3.9e-19 * (7e0 / 3e2)**1.8e0 * 
     $                      exp(2e1 / 7e0)
        endif
c
c Dissociative charge transfer of CO:
c  
c  CO + He+ -> C+ + O + He
c
c Rate from PDAN89
c
        chtab(92, itemp) = 1.4e-9 / sqrt(temp3)
        if (chtab(92,itemp) .gt. 1.6e-9) then
          chtab(92,itemp) = 1.6e-9
        endif
c
c  CO + He+ -> C + O+ + He
c
c Rate from PDAN89
c
        chtab(93, itemp) = 1.4e-16 / sqrt(temp3)
c
c Charge transfer: C + He+ -> C+ + He
c
c Rate is my fit to data from KD93 and is valid for 10 < T < 10000K.
c The fit is accurate to within 20% over this range. Note that UMIST99 
c cites the same source, but fits only the restricted range 10 < T < 300K. 
c
        if (temp .le. 2e2) then
          chtab(94, itemp) = 8.58e-17 * temp**0.757e0
        elseif (temp .le. 2e3) then
          chtab(94, itemp) = 3.25e-17 * temp**0.968e0
        else
          chtab(94, itemp) = 2.77e-19 * temp**1.597e0
        endif
c
c Charge transfer: Si++ + H -> Si+ + H+
c
c Fit from KF96, quoting GMV82; valid for T < 10^6 K, assumed constant
c above this.
c
        if (temp .le. 1e6) then
          chtab(95, itemp) = 1.23e-9 * (temp / 1e4)**0.24e0 *
     $                     (1e0 + 3.17e0 * exp(4.18e-7 * temp))
        else
          ifix = int(1 + (6e0 / dtlog))
          chtab(95, itemp) = chtab(95, ifix)
        endif
c
c Charge transfer: Si+ + H+ -> Si++ + H
c
c Fit from KF96, quoting GMV82; valid for T < 10^6 K, assumed constant
c above this.
c
        if (temp .le. 1e6) then
          chtab(96, itemp) = 4.10e-10 * (temp / 1e4)**0.24e0 *
     $                     (1e0 + 3.17e0 * exp(4.18e-7 * temp))
     $                     * exp(-3.178e0 / te)
        else
          ifix = int(1 + (6e0 / dtlog))
          chtab(96, itemp) = chtab(96, ifix)
        endif
c
c Charge transfer: O + He+ -> O+ + He
c
c Rate coefficient from Z04, who quote a fitting accuracy of 6%. 
c Note that this is the rate for radiative charge transfer 
c (i.e. it is actually O + He+ -> O+ + He + PHOTON), which is the 
c dominant mechanism at all temperatures.
c
c Fit is valid for 10 < T < 1e6 K; at higher T, I assume that the 
c rate is constant, to avoid blow-up
c
        if (temp .le. 1e6) then
          chtab(97, itemp) = 4.991e-15 * (temp / 1e4)**0.3794e0 *
     $                       exp(-temp / 1.121e6) + 2.780e-15 *
     $                       (temp / 1e4)**(-0.2163e0) * 
     $                      exp(temp / 8.158e5)
        else
          ifix = int(1 + (6e0 / dtlog))
          chtab(97, itemp) = chtab(97, ifix)
        endif
c
c Deuterium chemistry.
c
c N.B. D+ recombination, D collisional ionization, D cosmic ray
c ionization rates are all assumed to be the same as for normal hydrogen
c
c Deuterium charge transfer:
c
c D + H+ -> D+ + H
c
c From S02. Published fit breaks down at v.low temperature and at
c temperatures greater than about 10^6 K, and gives negative values,
c which are clearly wrong. Since the range of validity of the fit is
c 2.7 < T < 2e5 K, and since we don't expect this process to be 
c important at temperatures outside of this range (since we probably
c won't reach T < T_CMB, and at high T everything will be ionized),
c the modifications below help to keep everything sane (and positive)
c without introducing any real error
c
        if (temp .le. 2.7e0) then
          chtab(98,itemp) = 2e-10 * temp**0.402e0 * exp(-37.1e0 / temp)
        elseif (temp .le. 2e5) then
          chtab(98,itemp) = 2e-10 * temp**0.402e0 * exp(-37.1e0 / temp)
     $                    - 3.31e-17 * temp**1.48e0
        else
          chtab(98,itemp) = 3.44e-10 * temp**0.35e0
        endif
c
c D+ + H -> D + H+
c
c Also from S02. This fit does not suffer from the same problems as the
c previous one, and so does not need to be modified.
c
        chtab(99,itemp) = 2.06e-10 * temp**0.396e0 * exp(-33e0 / temp)
     $                  + 2.03e-9 * temp**(-0.332e0)
c
c HD formation:
c
c D+ + H2 -> HD + H+
c
c From G82, as quoted by GP02.
c
        chtab(100,itemp) = (0.417e0 + 0.846e0 * log10(temp) 
     $                    - 0.137e0 * log10(temp)**2e0) * 1e-9
c
c  HD + H+ -> H2 + D+
c
c From G82, as quoted by GP02.
c
        chtab(101,itemp) = 1.1e-9 * exp(-4.88e2 / temp)
c
c  D + H2 -> H + HD
c
c Rate is my fit to data in MPS03
c
        if (temp .le. 2e3) then
          chtab(102,itemp) = 1e1**(-5.64737e1 + 5.88886e0 * log10(temp)
     $                            + 7.19692e0  * log10(temp)**2e0 
     $                            + 2.25069e0  * log10(temp)**3e0
     $                            - 2.16903e0  * log10(temp)**4e0
     $                            + 3.17887e-1 * log10(temp)**5e0) 
        else
          chtab(102,itemp) = 3.17e-10 * exp(-5.207e3 / temp)
        endif
c
c  HD + H -> H2 + D
c
c Rate from S59, as quoted by GP02. The fit in GP02 is only good down
c to about 200K, and so I have modified it at lower T to stop it from
c blowing up. 
c 
        if (temp .le. 2e2) then
          chtab(103,itemp) = 5.25e-11 * exp(-4.43e3 / temp) 
        else
          chtab(103,itemp) = 5.25e-11 * exp((-4.43e3 / temp) +
     $                      (1.739e5 / temp**2e0))
        endif
c
c CH+ formation by RA: C+ + H -> CH+ + gamma
c
c Rate from BvH06. Fit accurate to within 2% for 5 < T < 1000K
c
        chtab(104, itemp) = 4.46e-16 * exp(-4.93e0 * 
     $                      temp**(-2e0/3e0)) / sqrt(temp)
c
c H3+ formation: H2+ + H2 -> H3+ + H
c
c Rate from SLD98, based on cross-section from LJB95 
c
        chtab(105, itemp) = 2.24e-9 * temp3**0.042e0 * 
     $                     exp(-temp / 4.66e4)
c
c H3+ destruction: H3+ DR  (total)
c
c Fit from UMIST06, based on data from MAC04.
c Branching fractions are in rate_eq_mol.F

        chtab(106, itemp) = 6.7e-8 / temp3**0.52e0
c
c H3+ destruction:  H3+ + H -> H2+ + H2
c
c Rate from GP98, fit to SMT92
c
        chtab(107, itemp) = 7.7e-9 * exp(-1.756e4 / temp)
c
c Three-body H2 formation: H + H + H -> H2 + H
c
c NB. Highly uncertain -- the PSS83 & ABN02 rates, both of which are
c in common usage, disagree by more than an order of magnitude, and
c the CW83 rate doesn't agree with either of them
c
        if (iflag_3bh2a .eq. 1) then
c From ABN02, based on Or87 at T < 300K, simple power-law fit at T > 300K
          if (temp .le. 3e2) then
            chtab(108, itemp) = 1.3e-32 / temp3**0.38e0
          else
            chtab(108, itemp) = 1.3e-32 / temp3
          endif
        elseif (iflag_3bh2a .eq. 2) then
c From PSS83, based on JGC67
          chtab(108, itemp) = 5.5e-29 / temp
        elseif (iflag_3bh2a .eq. 3) then
c From CW83
          chtab(108, itemp) = 8.8e-33
        elseif (iflag_3bh2a .eq. 4) then
c From FH07
          chtab(108, itemp) = 1.44e-26 / temp**1.54e0
        elseif (iflag_3bh2a .eq. 5) then
c From G07
          chtab(108, itemp) = 7.7e-31 / temp**0.464e0
        else
          print*, "Unknown value specified for iflag_3bh2a:", 
     $            iflag_3bh2a
          ABORT(428)
        endif
c
c Three-body H2 formation: H + H + H2 -> H2 + H2
c
        if (iflag_3bh2b .eq. 1) then
c From PSS83, based on JGC67
          chtab(109, itemp) = 6.875e-30 / temp
        elseif (iflag_3bh2b .eq. 2) then
c From CW83; note that BAU92 give a very similar rate, but with
c a slightly smaller numerical coefficient: 2.7d-31.
          chtab(109, itemp) = 2.8e-31 / temp**0.6e0
        elseif (iflag_3bh2b .eq. 3) then
c From FH07 -- one-eighth of their H + H + H rate
          chtab(109, itemp) = 1.8e-27 / temp**1.54e0
        elseif (iflag_3bh2b .eq. 4) then
c Assume 1/8th of ch108, following PSS83
          chtab(109,itemp) = chtab(108,itemp) / 8e0
        else
          print*, "Unknown value specified for iflag_3bh2b:", 
     $            iflag_3bh2b
          ABORT(429)
        endif
c
c Three-body H2 formation: H + H + He -> H2 + He 
c
c Rate is from WK75, based on measurements between 77K and 298K.
c The value quoted by MR77 at 297K is approx. 18% lower, but this
c is within the quoted errors on the latter measurement.
c
        chtab(110, itemp) = 6.9e-32 * temp**(-0.4e0)
c
c THREE-BODY CHEMISTRY OF HEAVY MOLECULES. (CO, OH, H2O etc.)
c Unless otherwise noted, we assume that these rates are independent
c of the identity of the third body.
c
c Three-body CO formation: C + O + M -> CO + M
c
c Rate is from PDAN89, for M=He. We assume that the rates for M=H,
c M=H2 are similar. 
c
c Quoted fit valid for 2000 < T < 10000; at T < 2000, we assume that
c we can simply drop the exponential T dependence.
c
c NB: Three-body formation of HCO by: CO + H + M -> HCO + M can be 
c ignored, as almost all HCO resulting is destroyed by HCO + H -> CO + H2.
c Equilibrium abundance of HCO/CO is ~< 2e-8 at n <= 1e16 cm^-3.
c 

        if (temp .le. 2e3) then
          chtab(111, itemp) = 2.14e-29 * temp3**(-3.08e0) * 
     $                        exp(1.057e0)
        else
          chtab(111, itemp) = 2.14e-29 * temp3**(-3.08e0) * 
     $                        exp(2.114e3 / temp)
        endif
c
c Three-body CO+ formation: C+ + O  + M -> CO+ + M
c                           C  + O+ + M -> CO+ + M
c
c Rate is from PDAN89, for M=He. We assume that the rates for M=H,
c M=H2 are similar. Rate adjusted at T < 2000K as noted for ch111.
c Same rate is quoted for both reactions, so we only tabulate it 
c once
c
        if (temp .le. 2e3) then
          chtab(112, itemp) = 2.14e-27 * temp3**(-3.08e0) * 
     $                        exp(1.057e0)
        else
          chtab(112, itemp) = 2.14e-27 * temp3**(-3.08e0) * 
     $                        exp(2.114e3 / temp)
        endif
c
c Three-body OH formation: H + O + M -> OH + M
c
c Rate from NIST chemical kinetics database, based on TH86
c
        chtab(113, itemp) = 4.33e-32 / temp3
c
c Three-body H2O formation: H + OH + M -> H2O + M
c
c Rate from BAU92. (NB Assumes a bath gas of Argon)
c
        chtab(114, itemp) = 2.56e-31 * temp3**(-2.0e0)
c
c Three-body formation of C2: C + C + M -> C2 + M
c 
c Rate at T > 5000K is the fit given in UMIST99 to data from FA69
c and SLA76. At lower temperatures, this fitting formula blows up,
c so we assume instead that the T^-1.6 temperature dependence 
c reported by SLA76 continues to hold below 5000K. If this assumption
c is incorrect, then this rate may be significantly inaccurate at low
c temperatures
c
        if (temp .le. 5e3) then
          chtab(115, itemp) = 5.99e-33 * (temp / 5e3)**(-1.6e0)
        else
          chtab(115, itemp) = 5.99e-33 * (temp / 5e3)**(-0.64e0) * 
     $                      exp(5.255e3 / temp)
        endif
c
c Three-body formation of O2: O + O + M -> O2 + M
c
c Rate from NIST chemical kinetics database, based on WAR84.
c Note that this process is generally unimportant in comparison
c to ch113, unless the HI fraction of the gas is very small,
c but we include it here for completeness.
c
        chtab(116, itemp) = 9.2e-34 / temp3
c
c HCO+ formation: O + CH -> HCO+ + e- 
c
c Rate from UMIST99, based on MB73
c
        chtab(117, itemp) = 2e-11 * temp3**0.44e0
c
c H2O formation from OH-OH collisions: OH + OH -> H2O + O
c
c Rate from SHC04.
c
        chtab(118, itemp) = 1.65e-12 * temp3**1.14e0 * 
     $                        exp(-5e1 / temp)
c
c CO+ dissociative recombination:  CO+ + e -> C + O
c
c Rate from OSU2005, based on R98
c
        chtab(119, itemp) = 2.75e-7 * temp3**(-0.55e0)
c
c C2+ formation from CH:  C+ + CH -> C2+ + H
c
c Rate from UMIST06. OSU2005 quotes a rate that is ~ 10x larger at all T
c Currently, we don't account for C2+ so we don't include this, but we now
c have the rate handy should we choose to do so in the future.
c
        chtab(120, itemp) = 3.8e-10 / sqrt(temp3)
c
c SiIII gas-phase recombination.
c
c Radiative rate: my fit to data of N95 & N96. Accurate to within 10%
c
        si2p_rec_rr = 1.75e-12 / (temp / 1e4)**0.6346e0
c
c Dielectronic rate: rate from MMCV98.
c
        si2p_rec_dr = (2.2552e-11 * exp(-2.76e0 / te) + 
     $                 5.6058e-9 * exp(-1.013e1 / te)) / te**1.5e0

        chtab(121, itemp) = si2p_rec_rr + si2p_rec_dr
c
c SiH+ formation by radiative association: Si+ + H -> SiH+ + gamma
c
c Rate from STA00
c
        chtab(122, itemp) = 1.17e-17 * temp3**(-0.14e0)
c
c Mg+ recombination -- currently unused
c
        chtab(123, itemp) = 2.8e-12 / temp3**0.86e0
c
c CO+ formation by radiative association: C+ + O -> CO+ + gamma
c
c Rate from DD90.
c
        if (temp .le. 3e2) then
          chtab(124, itemp) = 2.5e-18
        else
          chtab(124, itemp) = 3.14e-18 * temp3**(-0.15e0) * 
     $                        exp(6.8e1 / temp)
        endif
c
c H2+ + H-  -> H3+ + e-
c
c From LSD02, based on the cross-sections of NAJ98.
c Not included in main molecular chem. network
c
        if (temp .le. 1e5) then
          chtab(125, itemp) = 2.7e-10 * temp3**(-0.485e0) * 
     $                        exp(temp / 3.12e4)
        else
          chtab(125, itemp) = 2.7e-10 * (1e5 / 3e2)**(-0.485e0) * 
     $                        exp(1e5 / 3.12e4)
        endif
c
c Dissociative attachment: H2 + e- -> H + H-
c
c Rate from GP98, based on data from SA67.
c Assumes that the H2 is in its vibrational ground state and so
c underestimates rate at high densities; however, since the 
c electron fraction will generally be very small at high n, this
c is probably not a major problem
c
c Not included in main molecular chem. network
c
        chtab(126, itemp) = 2.7e-8 * temp**(-1.27e0) * 
     $                      exp(-4.3e4 / temp)
c
c Neutralization of H2+ with H-:  H2+ + H- -> H2 + H
c                                          -> H + H + H
c
c SLD98 quote the same rate for both processes, taken from DL87
c Probably unimportant in comparison to DR
c
        chtab(127, itemp) = 1.4e-7 / sqrt(temp3)
c
c HD+ formation from H2+ and D: H2+ + D -> HD+ + H
c
c Rate from SLD98, based on cross-section of LJB95
c
        chtab(128, itemp) = 1.07e-9 * temp3**0.062e0 * 
     $                      exp(-temp / 4.14e4)
c
c HD+ dissociative recombination:  HD+ + e -> H + D
c 
c Rate from GP98, fit to cross-section of SSS95. Note that SLD98 
c provide an alternative fit to the same cross-section which 
c differs by ~ 50-70% at low T, ~ 15% at high T. WFP04 quote a 
c rate with a shallower T dependence that differs from the GP98
c rate by ~20% at 200K, but is almost the same at 3000K.
c
        chtab(129, itemp) = 7.2e-8 / sqrt(temp)
c
c H2+ formation from HD+ and H: HD+ + H -> H2+ + D
c
c Rate from DM56, but scaled by D reduced mass
c
        chtab(130, itemp) = 1.0e-9 * exp(-1.54e2 / temp)
c
c Collisional dissociation of OH by H: OH + H -> O + H + H
c
c Rate from UMIST06.
c
        chtab(131, itemp) = 6e-9 * exp(-5.09e4 / temp)
c
c Collisional dissociation of OH by H2: OH + H2 -> O + H + H2
c
c Rate from UMIST06.
c
        chtab(132, itemp) = 6e-9 * exp(-5.09e4 / temp)
c
c Collisional dissociation of H2O by H: H2O + H -> OH + H + H
c
c Rate from UMIST06.
c
        chtab(133, itemp) = 5.8e-9 * exp(-5.29e4 / temp)
c
c Collisional dissociation of H2O by H2: H2O + H2 -> OH + H + H2
c
c Rate from UMIST06.
c
        chtab(134, itemp) = 5.8e-9 * exp(-5.29e4 / temp)
c
c ch135--ch136: used for HeII grain-surface recombination rate; see above
c
c H2 collisional dissociation by He: H2 + He -> H + H + He
c
c Rate from DRCM87. ch137, ch138 & ch139 are low n limit, 
c high n limit and critical density resp. Data valid for
c 2000 < T < 10000K, but fit behaves OK at lower T; at higher
c T, we just use the 10000K value.
c
c
c NB. We now use rate for LTE limit determined by detailed balance
c applied to reverse reaction
c
        chtab(138,itemp) = 6.6e-10 * temp**0.115e0 * exp(-5.2e4 / temp)
c
        if (temp .le. 1e4) then
          chtab(137, itemp) = 1e1**(3.801e0 * log10(temp) - 
     $                        27.029 - 2.9487e4 / temp)
c          chtab(138, itemp) = 1d1**(-1.75 * log10(temp) - 2.729 - 
c     $                         2.3474d4 / temp) 
          chtab(139, itemp) = 1e1**(5.0792e0 * (1e0 - 1.23e-5 * 
     $                        (temp - 2e3)))
        else
          tfix = 1e4
          chtab(137, itemp) = 1e1**(3.801e0 * log10(tfix) - 
     $                        27.029e0 - 2.9487e4 / tfix)
c          chtab(138, itemp) = 1d1**(-1.75 * log10(tfix) - 2.729 - 
c     $                         2.3474d4 / tfix)
          chtab(139, itemp) = 1e1**(5.0792e0 * (1e0 - 1.23e-5 * 
     $                        (tfix - 2e3)))
        endif
c
        chtab(139, itemp) = 1e0 / chtab(139, itemp)
c
c CH3 formation from CH2: CH2 + H2 -> CH3 + H
c
        chtab(140, itemp) = 5.18e-11 * temp3**0.17e0 * 
     $                      exp(-6.4e3 / temp)
c
c CH3 destruction: CH3 + H -> CH2 + H2
c
        chtab(141, itemp) = 1.0e-10 * exp(-7.6e3 / temp)
c
c CH4 formation from CH3: CH3 + H2 -> CH4 + H
c
        chtab(142, itemp) = 6.86e-14 * temp3**2.74e0 *exp(-4.74e3/temp)
c
c CH4 destruction: CH4 + H -> CH3 + H2
c
        chtab(143, itemp) = 5.82e-13 * temp3**3.e0 
     $                    * exp(-4.045e3 / temp)
c
c HI collisional ionization by HI: H + H -> H + H+ + e
c
c Rate from GVZ87, as quoted by LCS91
c
        chtab(144, itemp) = 1.2e-17 * temp**1.2e0 * exp(-1.578e5 / temp)
c
c Surface chem:
c
c Accretion rate of CO (assuming standard dust properties); rates for other
c species then follow from simple rescaling (mco/mx)^1/2
c 
        chtab(145, itemp) = 1.1e-17 * sqrt(temp)
c
c Gas phase recombination of representative trace metal M (used in NL99 network)
c
        chtab(146, itemp) = 3.8e-10 / temp**(0.65e0)
c
c Formation of CO from OHx [NL99]: C + OHx -> CO + H
c
        chtab(147, itemp) = 5.8e-12 * sqrt(temp)
c
c Collisional ionization of H by H: H + H -> H+ + e- + H
c (Dominates when fractional ionization v. small or zero)
c Rate from Omukai, 2001; assumes only ground-state populated
c
        chtab(148, itemp) = 1.2e-17 * temp**1.2e0 * exp(-1.5782e5/temp)
c
c High-temperature collisional dissociation of H2+ into H + H+.
c Approximate rate (factor of a few uncertainty at least), based
c on McClure, 1965, Phys. Rev. A, 140, 770. Not really important
c for the H2+ chemistry itself, but useful for stabilizing the
c solver at high T
c
        chtab(149, itemp) = 1e-12 * sqrt(temp) * exp(-3.1e4 / temp)
c
c Formation of NO:  N + OH -> NO + H
c
c Rate from Smith, Herbst & Chang (2004)
c
        chtab(150, itemp) = 7.5e-11 / temp3**0.18e0
c
c Formation of N2 from NO: NO + N -> N2 + O
c
c Rate from NIST database
c
        chtab(151, itemp) = 3.75e-11 * exp(-26e0 / temp)
c
c Dissociative recombination of N2H+:  N2H+ + e -> N2 + H
c
c Rate from Geppert et al (2004).
c
c NB: This rate is actually for the total DR rate: N2H+ + e -> products
c Geppert et al find a branching ratio of 36% to N2 + H, 64% to NH + H
c Because our v. simple N chem. doesn't include NH, we assume that all
c of it is converted back to N2 by reactions with atomic N. In reality,
c some will be photodissociated, and so will react with C to form CN...
c
        chtab(152, itemp) = 1e-6 / temp3**0.51e0
c
 200  continue
c
c Approximate temperature derivatives
c
      do 300 itemp = 1,nmd-1
         dtemp           = 1e0 / ( temptab(itemp + 1) - temptab(itemp) )
         do 251 k = 1, nchtab
            dtchtab(k,itemp) = ( chtab(k, itemp+1)
     $                         - chtab(k, itemp  ) ) * dtemp
c
 251     continue
 300  continue
c
c Manually add on the final tabulated value at itemp = ntabtemp
c
      do 311 k = 1, nchtab
         dtchtab(k,nmd) = dtchtab(k,nmd-1)
 311  continue
c
c
c Primary cosmic ray ionization rates - calculated as scale factor * HI
c cosmic ray ion. rate. For cr3, cr4 we use scaling factors from 
c Liszt (2003, A&A, 398, 621), based on data from the UMIST database. 
c For cr5 & cr7, we use L67, L78. For cr6, we use the scaling relative to 
c the HI rate given by UMIST.
c
c [NB Does not include photoionizations caused by CR-induced photons 
c (i.e photons produced by excitation of HI or H2 by collisions with
c energetic primary and secondary electrons)]
c
c  HI:
c
      crtab(1) = cosmic_ray_ion_rate
c
c  H2: 
c
      crtab(2) = 2e0 * cosmic_ray_ion_rate   ! H2 -> H2+ + e-
      crtab(10) = 3.7e-2 * cosmic_ray_ion_rate ! H2 -> H + H+ + e
      crtab(11) = 2.2e-1 * cosmic_ray_ion_rate ! H2 -> H + H
      crtab(12) = 6.5e-4 * cosmic_ray_ion_rate ! H2 -> H+ + H-
c
c  CI:
c
      crtab(3) = 3.83e0 * cosmic_ray_ion_rate
c
c  OI:
c
      crtab(4) = 5.67e0 * cosmic_ray_ion_rate
c
c  SiI:  
c
      crtab(5) = 6.5e0 * cosmic_ray_ion_rate
c
c  HeI: 
c
      crtab(6) = 1.09e0 * cosmic_ray_ion_rate
c
c  Si+:
c
      crtab(7) = 2.5e0 * cosmic_ray_ion_rate
c
c  CO: [NB. Direct ionization by cosmic rays] 
c
      crtab(9) = 6.52e0 * cosmic_ray_ion_rate
c
c Cosmic-ray induced photoprocesses
c
c C -> C+ + e-
      crphot(1)  = 2.8e3 * cosmic_ray_ion_rate
c CH -> C + H
      crphot(2)  = 4e3 * cosmic_ray_ion_rate
c CH+ -> C+ + H
      crphot(3)  = 9.6e2 * cosmic_ray_ion_rate
c CH2 -> CH2+ + e
      crphot(4)  = 2.7e3 * cosmic_ray_ion_rate
c CH2 -> CH + H
      crphot(5)  = 2.7e3 * cosmic_ray_ion_rate
c C2 -> C + C
      crphot(6)  = 1.3e3 * cosmic_ray_ion_rate
c OH -> O + H
      crphot(7)  = 2.8e3 * cosmic_ray_ion_rate
c H2O -> OH + H
      crphot(8)  = 5.3e3 * cosmic_ray_ion_rate
c O2 -> O + O
      crphot(9)  = 4.1e3 * cosmic_ray_ion_rate
c O2 -> O2+ + e
      crphot(10) = 6.4e2 * cosmic_ray_ion_rate
c CO -> C + O
      crphot(11) = 0.21e0 * cosmic_ray_ion_rate
c NO -> N + O
c Note: this is actually the sum of two processes, direct photodissoc.
c and photoionization to NO+, which we assume then dissociatively 
c recombines to yield N + O...
      crphot(12) = 5.2e3 * cosmic_ray_ion_rate
c
c Setup tables for CO shielding
c
      do I = 1, imax_COss
        NCO = log10(NCO_shield_min) + 0.05e0 * (I - 1)        
        do J = num_CO_self_shield, 1, -1
          if (NCO .ge. CO_shield_columns(J)) then
            goto 312
          endif
        enddo
 312    continue
        dN = (NCO - CO_shield_columns(J)) / 
     $        (CO_shield_columns(J+1) - CO_shield_columns(J))
        CO_self_shielding(I) = dN * CO_shield_factors(J+1) 
     $                       + (1e0 - dN) * CO_shield_factors(J)
      enddo
c
      do I = 1, imax_COH2
        NH2 = log10(NH2_shield_min) + 0.05e0 * (I - 1)        
        do J = num_CO_H2_shield, 1, -1
          if (NH2 .ge. H2_shield_columns(J)) then
            goto 313
          endif
        enddo
 313    dN = (NH2 - H2_shield_columns(J)) / 
     $        (H2_shield_columns(J+1) - H2_shield_columns(J))
        CO_H2_shielding(I) = dN * H2_shield_factors(J+1) 
     $                       + (1e0 - dN) * H2_shield_factors(J)

      enddo
c
c Initialize optically-thin photochemical rates (e.g. CI ionization).
c This is done in a separate subroutine to allow us to easily 
c substitute new rates if the incident spectrum changes
c
      call photoinit
c
c Initialize constant rate coefficients -- i.e. ones which are not 
c temperature or density dependent
c
      call const_rates
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             C H E M I N M O               \\\\\\\\\\
c
c=======================================================================

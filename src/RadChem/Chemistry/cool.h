!
!  Originally based on work by G. Suttner (Univ. Wuerzburg, 1995) and  
!  M. D. Smith (Armagh Observatory, 2000-2001).
!  Completely rewritten and substantially extended by 
!  S. Glover (AMNH, 2002-2005, AIP 2006-2007)
!

#define UNUSED_PARAM(x)  x = x
#define ABORT(x) stop

#define DTCOOL_SCALE_FACTOR 1d-1

      REAL tiny_number
      parameter (tiny_number = 1e-16)
!
! He:H ratio by number (=> ratio by mass is 4*abhe)
!
      REAL abhe
      parameter(abhe = 0.1e0)

! Map symbolic indices to numbers

	integer ih2, ihp, itmp, iHatom, icp, ico, iion, idis
      integer ichx, iohx, ihcop, ihep, imp
      integer iHeatom, iCatom, iOatom, iMatom

      parameter(ih2 = 1)
      parameter(ihp = 2)
#if CHEMISTRYNETWORK == 4
      parameter(iHatom = 3)
      parameter(itmp   = 3)
	    parameter(iion   = 4)
	    parameter(idis   = 5) 
#endif
#if CHEMISTRYNETWORK == 5
      parameter(ico    = 3)
      parameter(iHatom = 4)
      parameter(icp    = 5)
      parameter(itmp   = 4)
	    parameter(iion   = 5)
	    parameter(idis   = 6) 
#endif
#if CHEMISTRYNETWORK == 15
      parameter(icp   = 3)
      parameter(ichx  = 4)
      parameter(iohx  = 5)
      parameter(ico   = 6)
      parameter(ihcop = 7)
      parameter(ihep  = 8)
      parameter(imp   = 9)
      parameter (iHatom  = 10)
      parameter (iHeatom = 11)
      parameter (iCatom  = 12)
      parameter (iOatom  = 13)
      parameter (iMatom  = 14)
      parameter(itmp  = 10)
		  parameter(iion  = 11)
	    parameter(idis  = 12) 
#endif
! Number of entries in cooling table
      integer nmd
      parameter(nmd = 10000)

! Number of cooling / heating rates computed in cooling fn.
      integer nrates
      parameter(nrates = 30)

! Number of cooling / heating rates computed in chemistry routines
      integer nrates_chem
      parameter(nrates_chem = 7)

! Total number of cooling / heating rates
      integer nrates_tot 
      parameter(nrates_tot = nrates + nrates_chem)

! Number of abundances passed to cooling function
! (Note that this is not necessarily the same as the number
! that we actually track as field variables)
      integer nabn
      parameter(nabn = 17)

! Boltzmann constant
      REAL kboltz
      parameter (kboltz = 1.38066e-16)

! One electron volt, in ergs
      REAL eV
      parameter (eV = 1.60219e-12)

! Number of different quantities stored in cooling look-up table
      integer ncltab
      parameter (ncltab = 86) 

! Number of different quantities stored in chemistry look-up table
      integer nchtab
      parameter (nchtab = 152) 

! Number of cosmic ray ionization rates tabulated
      integer ncrtab
      parameter(ncrtab = 12)

! Number of cosmic ray induced photoionizations/photodissociations tabulated
      integer ncrphot
      parameter(ncrphot = 12)

! Number of photochemical rates tabulated
      integer nphtab
      parameter(nphtab = 54)

! Number of constant rate coefficient initialized in const_rates
      integer nconst
      parameter(nconst = 91)

! These variables are initialized in cheminmo
      REAL chtab(nchtab,nmd), dtchtab(nchtab,nmd)
      REAL crtab(ncrtab), crphot(ncrphot)

! These variables are initialized in photoinit
      REAL phtab(nphtab), f_rsc

! This is initialized in const_rates
      REAL cst(nconst)

! These variables are initialized in coolinmo
      REAL temptab(nmd)
      REAL cltab(ncltab,nmd), dtcltab(ncltab,nmd)
      REAL dtlog,tmax,tmin
!
! CO rotational cooling
       integer nTco
       parameter (nTco = 1991)

       integer ncdco
       parameter (ncdco = 46)

       REAL co_temptab(nTco), co_colntab(ncdco)

       REAL co_L0(nTco), dTco_L0(nTco)
       REAL co_lte(ncdco,nTco), co_n05(ncdco,nTco), co_alp(ncdco,nTco)
       REAL dTco_lte(ncdco,nTco), dTco_n05(ncdco,nTco)
       REAL dTco_alp(ncdco,nTco)

! CO vibrational cooling
       integer nTco_vib
       parameter (nTco_vib = 3901)

       integer ncdco_vib
       parameter (ncdco_vib = 61)

       REAL co_vib_temptab(nTco_vib), co_vib_colntab(ncdco_vib)
       REAL co_vib_LTE_final(ncdco_vib, nTco_vib)
       REAL dTco_vib_LTE(ncdco_vib, nTco_vib)

       common /co_data/ co_temptab, co_colntab, co_L0, dTco_L0,            &
     &                  co_lte, co_n05, co_alp, dTco_lte, dTco_n05,        &
     &                  dTco_alp, co_vib_temptab, co_vib_colntab,          &
     &                  co_vib_LTE_final, dTco_vib_LTE

! H2O rotational cooling
       integer nTh2o
       parameter (nTh2o = 3991)

       integer ncdh2o
       parameter (ncdh2o = 91)

       REAL h2o_temptab(nTh2o), h2o_colntab(ncdh2o)

       REAL h2o_L0_ortho(nTh2o), dTh2o_L0_ortho(nTh2o)
       REAL  h2o_L0_para(nTh2o),  dTh2o_L0_para(nTh2o)

       REAL h2o_LTE_ortho(ncdh2o,nTh2o),                                  &
     &      h2o_n05_ortho(ncdh2o,nTh2o),                                  &
     &      h2o_alp_ortho(ncdh2o,nTh2o),                                  &
     &      h2o_LTE_para(ncdh2o,nTh2o),                                   &
     &      h2o_n05_para(ncdh2o,nTh2o),                                   &
     &      h2o_alp_para(ncdh2o,nTh2o)

       REAL dTh2o_LTE_ortho(ncdh2o,nTh2o),                                &
     &      dTh2o_n05_ortho(ncdh2o,nTh2o),                                &
     &      dTh2o_alp_ortho(ncdh2o,nTh2o),                                &
     &      dTh2o_LTE_para(ncdh2o,nTh2o),                                 &
     &      dTh2o_n05_para(ncdh2o,nTh2o),                                 &
     &      dTh2o_alp_para(ncdh2o,nTh2o)

! H2O vibrational cooling
       integer nTh2o_vib
       parameter (nTh2o_vib = 3901)

       integer ncdh2o_vib
       parameter (ncdh2o_vib = 61)

       REAL h2o_vib_temptab(nTh2o_vib), h2o_vib_colntab(ncdh2o_vib)
       REAL h2o_vib_LTE_final(ncdh2o_vib, nTh2o_vib)
       REAL dTh2o_vib_LTE(ncdh2o_vib, nTh2o_vib)

       common /h2o_data/ h2o_temptab, h2o_colntab, h2o_L0_ortho,          &
     &                   dTh2o_L0_ortho, h2o_L0_para, dTh2o_L0_para,      &
     &                   h2o_LTE_ortho, h2o_n05_ortho, h2o_alp_ortho,     &
     &                   h2o_LTE_para, h2o_n05_para, h2o_alp_para,        &
     &                   dTh2o_LTE_ortho, dTh2o_n05_ortho,                &
     &                   dTh2o_alp_ortho, dTh2o_LTE_para,                 &
     &                   dTh2o_n05_para, dTh2o_alp_para,                  &
     &                   h2o_vib_temptab, h2o_vib_colntab,                &
     &                   h2o_vib_LTE_final, dTh2o_vib_LTE

!
! These variables are initialized during problem setup
! 
      REAL deff, abundc, abundo, abundsi, abundD, abundM
      REAL abundN, G0
      REAL phi_pah, tdust, dust_to_gas_ratio
      REAL AV_conversion_factor, cosmic_ray_ion_rate, redshift
      REAL AV_ext, NH_ext, h2_form_ex, h2_form_kin, xray_scaling
      REAL Z_atom

      integer iphoto, iflag_mn, iflag_ad, iflag_atom
      integer iflag_3bh2a, iflag_3bh2b, iflag_h3pra
      integer no_chem, isrf_option, use_photo_eqb_approx

      common /coolr/ temptab, cltab, chtab, dtcltab, dtchtab,             &
     &               crtab, phtab, cst, dtlog, tdust, tmax, tmin,         &
     &               deff, abundc, abundo, abundsi, abundD,               &
     &               abundM, abundN, G0, f_rsc, phi_pah,                  &
     &               dust_to_gas_ratio, AV_conversion_factor,             &
     &               cosmic_ray_ion_rate, redshift, NH_ext,               &
     &               h2_form_ex, h2_form_kin, xray_scaling,               &
     &               AV_ext, Z_atom

      common /cooli/ iphoto, iflag_mn, iflag_ad, iflag_atom               &
     &,              iflag_3bh2a, iflag_3bh2b, iflag_h3pra                &
     &,              isrf_option, no_chem, use_photo_eqb_approx

#ifdef TREE_RAD
#define NPIX  12*NSIDE*NSIDE
#else
#define NPIX  1
#endif

!!
!! Written by S. Glover, AMNH, 2004-2005, AIP, 2006-2007
!!

      integer nchem_network 

#if CHEMISTRYNETWORK == 4 
      parameter(nchem_network = 4) !!CHEMISTRYNETWORK4
#endif

#if CHEMISTRYNETWORK == 5 
      parameter (nchem_network = 5) !!CHEMISTRYNETWORK5
#endif

!! Set up quantities (such as the absolute tolerances) that are used in 
!! multiple places in the non-equilibrium chemistry code. Note that most 
!! DVODE-specific setup should go in evolve_abundances.F -- nrpar & nipar
!! are exceptions, as they are used elsewhere, so it is useful to define
!! them here
!! 
      integer nrpar, nipar
      parameter (nrpar=15)
      parameter (nipar=3)

	integer num_non_eq_species
#if CHEMISTRYNETWORK == 4
      parameter (num_non_eq_species = 2)
#endif

#if CHEMISTRYNETWORK == 5
      parameter (num_non_eq_species = 3)
#endif

      integer nspec
      parameter (nspec = num_non_eq_species+3)
!!      REAL non_eq_abundances(num_non_eq_species)

      REAL ATOL(nspec), rtol(nspec)
      common /tolerance/ ATOL, rtol
!!!!      common /abundances/ non_eq_abundances

!! Amount by which abundances are allowed to stray over their theoretical
!! maximum before triggering an error in rate_eq -- set to a blanket value
!! of 1d-4 for the time being...
 
      REAL eps_max
      parameter (eps_max = 1e-4)

!!      REAL atol(nspec)
#define  RTOL_FIX  1e-4
#define  ATOL_H2   1e-3
#define  ATOL_HP   1e-8
#define  ATOL_CP   1e-12
#define  ATOL_HEP  1e-10
#define  ATOL_CO   1e-10
#define  ATOL_HCOP 1e-14
#define  ATOL_CHX  1e-10
#define  ATOL_OHX  1e-10
#define  ATOL_MP   1e-10
#define  ATOL_TMP  0e0
#define  ATOL_ION  1e0
#define  ATOL_DIS  1e0
!!      common /tolerance/ atol

!! DR: shutting off rate calculation if temperature very high (and overshoot risk)
!! this should be tested so more, so for now use unrealistically high values
!! housenumbers: 1e6 for ihp and 2e4 for ico
#define TMPIHP 1e99
#define TMPICO 1e99


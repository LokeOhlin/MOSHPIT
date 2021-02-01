c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////        I N I T _ C H E M I S T R Y        \\\\\\\\\\
c
c=======================================================================
c
c     Written by S. Glover (ITA 2009)
c
c     Calls chemistry setup routines, and sets default values for a 
c     number of parameters stored in the common blocks in cool.h
c
      subroutine init_chemistry
      implicit none
#include "param.h"
#include "cool.h"
c
      deff    = 1e0
      abundD  = 2.6e-5
      abundmg = 0e0
      f_rsc   = 1e0
      phi_pah = 0.5e0
      AV_conversion_factor = 5.348e-22
      h2_form_ex  = 0e0
      h2_form_kin = 0e0
c Metallicity to use for high T atomic cooling (in units of solar).
c At the moment, hard-wire this to 1.0 for simplicity. However, in
c the longer term, we should make this an adjustable parameter
      Z_atom = 1d0   
c
      iphoto   = 0
      iflag_mn = 1
      iflag_ad = 1
      iflag_atom = 3
      iflag_3bh2a = 1
      iflag_3bh2b = 1
      iflag_h3pra = 1
      idma_mass_option = 0
      no_chem = 0
c 
c All other values in the coolr and cooli common blocks in cool.h are either
c set within coolinmo & cheminmo (and the routines that they call), or must 
c be set from within FLASH itself
c
      call coolinmo
      call cheminmo
      call init_tolerances
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////        I N I T _ C H E M I S T R Y        \\\\\\\\\\
c
c=======================================================================

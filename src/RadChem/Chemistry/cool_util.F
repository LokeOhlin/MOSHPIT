c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                C O _ C O O L              \\\\\\\\\\
c
c=======================================================================
c
      subroutine co_cool(temp, N_co_eff, N_co_eff_vib, co_rot_L0, 
     $                   co_rot_lte, co_rot_alpha, co_rot_n05, 
     $                   co_vib_lte_out)
c
c     written by: Simon Glover, AMNH, 2004-2005, AIP, 2006
c
c  PURPOSE: Compute CO cooling rates based on tabulated Neufeld & 
c           Kaufman (1993) and Neufeld, Lepp & Melnick (1995) 
c           cooling functions [see coolinmo.F]
c
c  INPUT VARIABLES: temp, N_co_eff
c
c  OUTPUT VARIABLES: co_rot_L0, co_rot_lte, co_rot_alpha, co_rot_n05,
c                    co_vib_lte
c
c--------------------------------------------------------------------
c
      implicit none
#include "param.h"
#include "cool.h"
c
      REAL       temp, dtemp_co, dN_co, N_co_eff, N_co_eff_vib
      integer    itemp_co, iN_co
c
      REAL       co_rot_L0, co_rot_lte, co_rot_n05,
     $           co_rot_alpha, co_rot_lte_1, co_rot_lte_2,
     $           co_rot_alp_1, co_rot_alp_2, co_rot_n05_1,
     $           co_rot_n05_2, co_vib_lte_1, co_vib_lte_2,
     $           co_vib_lte_out
c
c CO rotational cooling
c
c No CO cooling below 5K, as no good data
c
      if (temp .lt. co_temptab(1)) then
        co_rot_L0    = 0.0
        co_rot_lte   = 0.0
        co_rot_n05   = 1.0
        co_rot_alpha = 1.0
        co_vib_lte_out = 0.0
        return
      elseif (temp .eq. co_temptab(1)) then
        itemp_co = 1
        dtemp_co = 0e0
      elseif (temp .ge. co_temptab(nTco)) then
        itemp_co = nTco
        dtemp_co = 0e0
      else
        itemp_co = int(temp) - 4    ! Table currently starts at 5K
        dtemp_co = temp - int(temp)
      endif
c
c For column densities that do not lie within the region covered by the
c NK93 or NLM95 data, we use the smallest or largest of the tabulated 
c values, as appropriate.
c
      if (N_co_eff .le. co_colntab(1)) then
        iN_co = 1
        dN_co = 0e0
      elseif (N_co_eff .ge. co_colntab(ncdco)) then
        iN_co = ncdco
        dN_co = 0e0
      else
        iN_co = int((10. * N_co_eff) - 144)
        dN_co = (N_co_eff - co_colntab(iN_co)) / 0.1e0
      endif
c
      co_rot_L0 = co_L0(itemp_co) + dtemp_co * dTco_L0(itemp_co)
c
      co_rot_lte_1 =   co_lte(iN_co,itemp_co) + dtemp_co *
     $               dTco_lte(iN_co,itemp_co)
      co_rot_alp_1 =   co_alp(iN_co,itemp_co) + dtemp_co *
     $               dTco_alp(iN_co,itemp_co)
      co_rot_n05_1 =   co_n05(iN_co,itemp_co) + dtemp_co *
     $               dTco_n05(iN_co,itemp_co)
      if (iN_co .eq. ncdco) then
        co_rot_lte   = co_rot_lte_1
        co_rot_alpha = co_rot_alp_1
        co_rot_n05   = co_rot_n05_1
      else
        co_rot_lte_2 =   co_lte(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_lte(iN_co+1,itemp_co)
        co_rot_alp_2 =   co_alp(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_alp(iN_co+1,itemp_co)
        co_rot_n05_2 =   co_n05(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_n05(iN_co+1,itemp_co)
c
        co_rot_lte = co_rot_lte_1 + (co_rot_lte_2 - co_rot_lte_1) *
     $               dN_co
        co_rot_alpha = co_rot_alp_1 + (co_rot_alp_2 - co_rot_alp_1) *
     $                 dN_co
        co_rot_n05 = co_rot_n05_1 + (co_rot_n05_2 - co_rot_n05_1) *
     $               dN_co
      endif
c
c Do final conversion to correct units:
c
      co_rot_L0  = 10e0**(-co_rot_L0)
      co_rot_lte = 10e0**(-co_rot_lte)
      co_rot_n05 = 10e0**(co_rot_n05)
c
c CO vibrational cooling
c
      if (temp .le. co_vib_temptab(1)) then
        itemp_co = 1
        dtemp_co = 0e0
      elseif (temp .ge. co_vib_temptab(nTco_vib)) then
        itemp_co = nTco_vib
        dtemp_co = 0e0
      else
        itemp_co = int(temp) - 99    ! Table starts at 100K
        dtemp_co = temp - int(temp)
      endif
c
      if (N_co_eff_vib .le. co_vib_colntab(1)) then
        iN_co = 1
        dN_co = 0e0
      elseif (N_co_eff_vib .ge. co_vib_colntab(ncdco_vib)) then
        iN_co = ncdco_vib
        dN_co = 0e0
      else
        iN_co = int((10 * N_co_eff) - 129)
        dN_co = (N_co_eff_vib - co_vib_colntab(iN_co)) / 0.1e0
      endif

      co_vib_lte_1 =   co_vib_lte_final(iN_co,itemp_co) + dtemp_co *
     $               dTco_vib_lte(iN_co,itemp_co)

      if (iN_co .eq. ncdco_vib) then
        co_vib_lte_out = co_vib_lte_1
      else
        co_vib_lte_2 =   co_vib_lte_final(iN_co+1,itemp_co) + dtemp_co *
     $                 dTco_vib_lte(iN_co+1,itemp_co)
        co_vib_lte_out = co_vib_lte_1 + (co_vib_lte_2 - co_vib_lte_1) *
     $               dN_co
      endif
c
      co_vib_lte_out = 10e0**(-co_vib_lte_out) * exp(-3.08e3 / temp)
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               C O _ C O O L               \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////          H 2 O _ R O T _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
      subroutine h2o_rot_cool(temp, N_h2o_eff, h2o_rot_L0, h2o_rot_lte, 
     $                        h2o_rot_alpha, h2o_rot_n05, iop_flag)
c
c     written by: Simon Glover, AMNH, 2004-2005, AIP 2006
c
c  PURPOSE: Compute H2O rotational cooling rates based on tabulated 
c           Neufeld & Kaufman (1993) and Neufeld, Lepp & Melnick (1995) 
c           cooling functions [see coolinmo.F]
c
c  INPUT VARIABLES: temp, N_h2o_eff, iop_flag (0 == para, 1 == ortho)
c
c  OUTPUT VARIABLES: h2o_rot_L0, h2o_rot_lte, h2o_rot_alpha, h2o_rot_n05
c
c--------------------------------------------------------------------
c
      implicit none
#include "param.h"
#include "cool.h"
c
      REAL     temp, dtemp_h2o, dN_h2o, N_h2o_eff
      integer  itemp_h2o, iN_h2o, iop_flag
c
      REAL     h2o_rot_L0, h2o_rot_lte, h2o_rot_n05,
     $         h2o_rot_alpha, h2o_rot_lte_1, h2o_rot_lte_2,
     $         h2o_rot_alp_1, h2o_rot_alp_2, h2o_rot_n05_1,
     $         h2o_rot_n05_2
c
      if (temp .le. h2o_temptab(1)) then
        itemp_h2o = 1
        dtemp_h2o = 0e0
      elseif (temp .ge. h2o_temptab(nTh2o)) then
        itemp_h2o = nTh2o
        dtemp_h2o = 0e0
      else
        itemp_h2o = int(temp) - 9    ! XXX: Table currently starts at 10K
        dtemp_h2o = temp - int(temp)
      endif
c
c For column densities that do not lie within the region covered by the
c NK93 or NLM95 data, we use the smallest or largest of the tabulated 
c values, as appropriate.
c
      if (N_h2o_eff .le. h2o_colntab(1)) then
        iN_h2o = 1
        dN_h2o = 0e0
      elseif (N_h2o_eff .ge. h2o_colntab(ncdh2o)) then
        iN_h2o = ncdh2o
        dN_h2o = 0e0
      else
        iN_h2o = int((10 * N_h2o_eff) - 144)
        dN_h2o = (N_h2o_eff - h2o_colntab(iN_h2o)) / 0.1e0
      endif
c
      if (iop_flag .eq. 0) then
        h2o_rot_L0 = h2o_L0_para(itemp_h2o) + dtemp_h2o * 
     $             dTh2o_L0_para(itemp_h2o)
        h2o_rot_lte_1 = h2o_lte_para(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_lte_para(iN_h2o,itemp_h2o)
        h2o_rot_alp_1 = h2o_alp_para(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_alp_para(iN_h2o,itemp_h2o)
        h2o_rot_n05_1 = h2o_n05_para(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_n05_para(iN_h2o,itemp_h2o)
        if (iN_h2o .eq. ncdh2o) then
          h2o_rot_lte   = h2o_rot_lte_1
          h2o_rot_alpha = h2o_rot_alp_1
          h2o_rot_n05   = h2o_rot_n05_1
        else
          h2o_rot_lte_2 = h2o_lte_para(iN_h2o+1,itemp_h2o) + 
     $      dtemp_h2o * dTh2o_lte_para(iN_h2o+1,itemp_h2o)
          h2o_rot_alp_2 = h2o_alp_para(iN_h2o+1,itemp_h2o) + 
     $      dtemp_h2o * dTh2o_alp_para(iN_h2o+1,itemp_h2o)
          h2o_rot_n05_2 = h2o_n05_para(iN_h2o+1,itemp_h2o) + 
     $      dtemp_h2o * dTh2o_n05_para(iN_h2o+1,itemp_h2o)
c
          h2o_rot_lte = h2o_rot_lte_1 + (h2o_rot_lte_2 - 
     $                  h2o_rot_lte_1) * dN_h2o
          h2o_rot_alpha = h2o_rot_alp_1 + (h2o_rot_alp_2 - 
     $                    h2o_rot_alp_1) * dN_h2o
          h2o_rot_n05 = h2o_rot_n05_1 + (h2o_rot_n05_2 - 
     $                  h2o_rot_n05_1) * dN_h2o
        endif
      else
        h2o_rot_L0 = h2o_L0_ortho(itemp_h2o) + dtemp_h2o * 
     $             dTh2o_L0_ortho(itemp_h2o)
        h2o_rot_lte_1 = h2o_lte_ortho(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_lte_ortho(iN_h2o,itemp_h2o)
        h2o_rot_alp_1 = h2o_alp_ortho(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_alp_ortho(iN_h2o,itemp_h2o)
        h2o_rot_n05_1 = h2o_n05_ortho(iN_h2o,itemp_h2o) + 
     $    dtemp_h2o * dTh2o_n05_ortho(iN_h2o,itemp_h2o)
        if (iN_h2o .eq. ncdh2o) then
          h2o_rot_lte   = h2o_rot_lte_1
          h2o_rot_alpha = h2o_rot_alp_1
          h2o_rot_n05   = h2o_rot_n05_1
        else 
          h2o_rot_lte_2 = h2o_lte_ortho(iN_h2o+1,itemp_h2o) + 
     $        dtemp_h2o * dTh2o_lte_ortho(iN_h2o+1,itemp_h2o)
          h2o_rot_alp_2 = h2o_alp_ortho(iN_h2o+1,itemp_h2o) + 
     $        dtemp_h2o * dTh2o_alp_ortho(iN_h2o+1,itemp_h2o)
          h2o_rot_n05_2 = h2o_n05_ortho(iN_h2o+1,itemp_h2o) + 
     $        dtemp_h2o * dTh2o_n05_ortho(iN_h2o+1,itemp_h2o)
c
          h2o_rot_lte = h2o_rot_lte_1 + (h2o_rot_lte_2 - 
     $                  h2o_rot_lte_1) * dN_h2o
          h2o_rot_alpha = h2o_rot_alp_1 + (h2o_rot_alp_2 - 
     $                    h2o_rot_alp_1) * dN_h2o
          h2o_rot_n05 = h2o_rot_n05_1 + (h2o_rot_n05_2 - 
     $                  h2o_rot_n05_1) * dN_h2o
        endif
      endif
c
c Do final conversion to correct units:
c
      h2o_rot_L0  = 10e0**(-h2o_rot_L0)
      h2o_rot_lte = 10e0**(-h2o_rot_lte)
      h2o_rot_n05 = 10e0**(h2o_rot_n05)
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////          H 2 O _ R O T _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////          H 2 O _ V I B _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
      subroutine h2o_vib_cool(temp, N_h2o_eff, h2o_vib_lte_out)
c
c     written by: Simon Glover, AMNH, 2004-2005, AIP 2006
c
c  PURPOSE: Compute H2O vibrational cooling rates based on tabulated 
c           Neufeld & Kaufman (1993) and Neufeld, Lepp & Melnick (1995) 
c           cooling functions [see coolinmo.F]
c
c  INPUT VARIABLES: temp, N_h2o_eff
c
c  OUTPUT VARIABLES: h2o_vib_lte_out
c
c--------------------------------------------------------------------
c
      implicit none
#include "param.h"
#include "cool.h"
c
      REAL     temp, dtemp_h2o, dN_h2o, N_h2o_eff
      integer  itemp_h2o, iN_h2o
c
      REAL     h2o_vib_lte_1, h2o_vib_lte_2,
     $         h2o_vib_lte_out
c
c H2O vibrational cooling
c
      if (temp .le. h2o_vib_temptab(1)) then
        itemp_h2o = 1
        dtemp_h2o = 0e0
      elseif (temp .ge. h2o_vib_temptab(nTh2o_vib)) then
        itemp_h2o = nTh2o_vib
        dtemp_h2o = 0e0
      else
        itemp_h2o = int(temp) - 99    ! Table starts at 100K
        dtemp_h2o = temp - int(temp)
      endif
c
      if (N_h2o_eff .le. h2o_vib_colntab(1)) then
        iN_h2o = 1
        dN_h2o = 0e0
      elseif (N_h2o_eff .ge. h2o_vib_colntab(ncdh2o_vib)) then
        iN_h2o = ncdh2o_vib
        dN_h2o = 0e0
      else
        iN_h2o = int((10 * N_h2o_eff) - 129)
        dN_h2o = (N_h2o_eff - h2o_vib_colntab(iN_h2o)) / 0.1e0
      endif

      h2o_vib_lte_1 =   h2o_vib_lte_final(iN_h2o,itemp_h2o) + 
     $               dtemp_h2o * dTh2o_vib_lte(iN_h2o,itemp_h2o)

      if (iN_h2o .eq. ncdh2o_vib) then
        h2o_vib_lte_out = h2o_vib_lte_1
      else
        h2o_vib_lte_2 =   h2o_vib_lte_final(iN_h2o+1,itemp_h2o) + 
     $                 dtemp_h2o * dTh2o_vib_lte(iN_h2o+1,itemp_h2o)
        h2o_vib_lte_out = h2o_vib_lte_1 + (h2o_vib_lte_2 - 
     $                    h2o_vib_lte_1) * dN_h2o
      endif
c
      h2o_vib_lte_out = 10e0**(-h2o_vib_lte_out) * exp(-2.325e3 / temp)
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////          H 2 O _ V I B _ C O O L          \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
c
c=======================================================================
c
      subroutine three_level_pops(r01, r02, r12, r10, r20, r21, 
     $                            n0, n1, n2)
      implicit none
#include "param.h"
#include "cool.h"
      REAL r01, r02, r12, r10, r20, r21
      REAL n0 , n1 , n2
      REAL a1 , a2 , a3 , b1 , b2 , b3
c
c If excitation rates are negligibly small, then we assume that all
c of the atoms are in level 0:
c
      if (r01 .eq. 0e0 .and. r02 .eq. 0e0) then
        n0 = 1e0
        n1 = 0e0
        n2 = 0e0
        return
      endif

      a1 = r01 + r02
      a2 = -r10
      a3 = -r20
      b1 = r01
      b2 = -(r10 + r12)
      b3 = r21
c
      n2 = -a1 * (a1 * b2 - b1 * a2) / ((a1 - a2) * 
     $     (a1 * b3 - b1 * a3) - (a1 - a3) * 
     $     (a1 * b2 - b1 * a2))
c
      n1 = (a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2
c
      n0 = 1e0 - n1 - n2
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////          C O M P U T E _ S T I M          \\\\\\\\\\
c
c=======================================================================
c
      subroutine compute_stim(a10, e10, rad_temp, b10)
      implicit none
#include "param.h"
#include "cool.h"
      REAL a10, e10, rad_temp, b10
      REAL x

      x = e10 / (kboltz * rad_temp)
      if (x .lt. 5e0) then
        b10 = a10 / (exp(x) - 1e0)
      else
        b10 = 0e0
      endif
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////          C O M P U T E _ S T I M          \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////        L O A D _ H 2 _ T A B L E          \\\\\\\\\\
c
c=======================================================================
c
      subroutine load_H2_table
#if CHEMISTRYNETWORK == 1
c At the moment, load this only when using primordial chemistry
      implicit none
#include "param.h"
#include "cool.h"
#include "mol_data.h"
      integer I, J

      open(12, file='H2-cooling-ratios.dat', status='old')
      do I = 1, nh2op
        do J = 1, nh2op
          read(12,*) h2_opac_temp(I), h2_opac_column(J), h2_opac(I,J)
        enddo
      enddo
      close (12, status='keep')
#endif /* CHEMISTRYNETWORK */
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////         L O A D _ H 2 _ T A B L E         \\\\\\\\\\
c
c=======================================================================

c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////    C O M P U T E _ H 2 _ O P A C I T Y    \\\\\\\\\\
c
c=======================================================================
c
      subroutine compute_h2_opacity(temp, N_H2_eff, opac)
      implicit none
#include "param.h"
#include "cool.h"
#include "mol_data.h"
      REAL temp, N_H2_eff, opac
      REAL column_min, column_max, logN, logT, diff, dN, dT
      REAL opac_tmp(2)
      integer ii, jj

      column_min = 10e0**(h2_opac_column(1))
      column_max = 10e0**(h2_opac_column(nh2op))

      if (N_H2_eff .le. column_min) then
        opac = 1e0
        return
      elseif (N_H2_eff .ge. column_max) then
        opac = 0e0
        return
      else
        logN = log10(N_H2_eff)
        jj   = 1e0 + int(10e0 * (logN - 17.0))
        diff = h2_opac_column(jj+1) - h2_opac_column(jj)
        dN   = (logN - h2_opac_column(jj)) / diff
c
        logT = log10(temp)
        if (logT .le. h2_opac_temp(1)) then
          ii = 1
          dT = 0e0
        elseif (logT .ge. h2_opac_temp(nh2op)) then  
          ii = nh2op
          dT = 0e0
        else
          ii   = 1 + int((logT - 1.5) / 0.03)
          diff = h2_opac_temp(ii+1) - h2_opac_temp(ii)
          dT   = (logT - h2_opac_temp(ii)) / diff
        endif
c
        if (dT .gt. 0e0) then
          opac_tmp(1) = h2_opac(ii,jj) + dT * (h2_opac(ii+1,jj) 
     $                - h2_opac(ii,jj))
          opac_tmp(2) = h2_opac(ii,jj+1) + dT * (h2_opac(ii+1,jj+1) 
     $                - h2_opac(ii,jj+1))
        else
          opac_tmp(1) = h2_opac(ii,jj)
          opac_tmp(2) = h2_opac(ii,jj+1)
        endif
c
        opac = opac_tmp(1) + dN * (opac_tmp(2) - opac_tmp(1))
        opac = 1e1**opac
      endif
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////    C O M P U T E _ H 2 _ O P A C I T Y    \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////        COMPUTE H2 HEATING FRACTION        \\\\\\\\\\
c
c=======================================================================
c
      subroutine compute_h2_heating_fraction(T, n, xh2, frac)
c Currently disabled
#if 0
      implicit none
#include "param.h"
#include "cool.h"
#include "h2heat.h"
      REAL T, n, xh2, frac
      REAL logT, logn, dT, dn, nscale, dh2
      REAL f1, f2, f3, f4, fr1, fr2, frac1, frac2
      integer idx1, idx2, idx3, offset

      logT = log10(T)
      logn = log10(n)

      if (logT .lt. 1e0) then
        idx1 = 1
        dT = 0e0
      elseif (logT .gt. 3.5e0) then
        idx1 = 11
        dT = 0e0
      else
        idx1 = 1 + int(4e0 * (logT - 1e0))
        dT = 4e0 * (logT - 1e0) - (idx1 - 1)
      endif

      if (logn .lt. 1e0) then
        idx2   = 1
        dn     = 0e0
        nscale = (n / 10.e0)
      elseif (logn .gt. 8e0) then
        idx2   = 29
        dn     = 0e0
        nscale = 1e0
      else
        idx2 = 1 + int(4e0 * (logn - 1e0))
        dn = 4e0 * (logn - 1e0) - (idx2 - 1)
        nscale = 1e0
      endif

      if (xh2 .eq. 0.5e0) then
        idx3 = 6
        dh2  = 0e0
      else
        idx3 = 1 + int(xh2 / 1e-1)
        dh2  = xh2 / 1e-1 - (idx3 - 1)
      endif
c
c Compute frac for idx3
c 
      offset = 319 * (idx3 - 1) + 29 * (idx1 - 1) + idx2

      if (dT .eq. 0e0) then
        if (dn .eq. 0e0) then
          frac1 = nscale * h2_rv_heat(offset)
        else
          f1 = h2_rv_heat(offset)
          f2 = h2_rv_heat(offset + 1)
          frac1 = f1 + dn * (f2 - f1)
        endif
      else
        if (dn .eq. 0e0) then
          f1 = nscale * h2_rv_heat(offset)
          f2 = nscale * h2_rv_heat(offset + 29)
          frac1 = f1 + dT * (f2 - f1)
        else
          f1 = h2_rv_heat(offset)
          f2 = h2_rv_heat(offset + 1)
          f3 = h2_rv_heat(offset + 29)
          f4 = h2_rv_heat(offset + 30)
          fr1 = f1 + dT * (f2 - f1)
          fr2 = f3 + dT * (f4 - f3)
          frac1 = fr1 + dn * (fr2 - fr1)
        endif
      endif
c
c Compute frac for idx3+1 (if required)
c 
      if (dh2 .eq. 0e0) then
        frac = frac1
      else
        offset = 319 * idx3 + 29 * (idx1 - 1) + idx2

        if (dT .eq. 0e0) then
          if (dn .eq. 0e0) then
            frac2 = nscale * h2_rv_heat(offset)
          else
            f1 = h2_rv_heat(offset)
            f2 = h2_rv_heat(offset + 1)
            frac2 = f1 + dn * (f2 - f1)
          endif
        else
          if (dn .eq. 0e0) then
            f1 = nscale * h2_rv_heat(offset)
            f2 = nscale * h2_rv_heat(offset + 29)
            frac2 = f1 + dT * (f2 - f1)
          else
            f1 = h2_rv_heat(offset)
            f2 = h2_rv_heat(offset + 1)
            f3 = h2_rv_heat(offset + 29)
            f4 = h2_rv_heat(offset + 30)
            fr1 = f1 + dT * (f2 - f1)
            fr2 = f3 + dT * (f4 - f3)
            frac2 = fr1 + dn * (fr2 - fr1)
          endif
        endif

        frac = frac1 + dh2 * (frac2 - frac1)
      endif
c
#endif /* 0 */
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////        COMPUTE H2 HEATING FRACTION        \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////        C A L C _ D U S T _ T E M P        \\\\\\\\\\
c
c=======================================================================
      subroutine calc_dust_temp(yn, chi_isrf, temp, abHI, RH2, Rdust,
     $                          Eabs, nocalc)
      implicit none
#include "param.h"
#include "cool.h"
      REAL yn, chi_isrf, temp, abHI, RH2, Rdust, Eabs
      REAL isrf_heat, tdust_new, dust_cool, gas_dust, h2_form_dust
      REAL temp2, tdust2, fa, stick, dfa, dtrdust
      REAL dust_h2_heat, F_nr, dF_nr, Rdust_rat
      REAL cl14, ch7, dtemp
      REAL field_strength
      integer K, itemp, nocalc

      integer ndust_iter
      parameter (ndust_iter = 20000)

c      print*,'dust_to_gas_ratio=',dust_to_gas_ratio
c      print*,'tdust=',tdust

      if (dust_to_gas_ratio .gt. 0e0) then        
        if (temp .le. tmin) then
          itemp = 1
          dtemp = 0e0
        elseif (temp .ge. tmax) then
          itemp = nmd
          dtemp = 0e0
        else
          itemp = int(log10(temp) / dtlog) + 1
          dtemp = temp - temptab(itemp)
        endif
c
        ch7  = chtab(7, itemp) + dtemp * dtchtab(7, itemp)
        cl14 = cltab(14,itemp) + dtemp * dtcltab(14,itemp)
c
        if (nocalc .eq. 1) then
c Only calculate new RH2, Rdust values, don't update tdust
          temp2  = temp * 1e-2
          tdust2 = tdust * 1e-2
          fa     = 1e0 / (1e0 + 1e4 * exp(-6e2 / tdust))
          stick  = 1e0 / (1e0 + 0.4e0 * (temp2 + tdust2)**0.5e0 +
     $                    0.2e0 * temp2 + 0.08e0 * temp2**2e0)
          Rdust  = ch7 * fa * stick
          RH2 = Rdust * abHI * yn**2e0
          return
        endif
c
c We solve here for the equilibrium dust temperature:
c
c Heating from interstellar radiation field. Depending on the value of the
c isrf_option parameter, we either use the ISRF from Mathis et al (1983)
c [isrf_option = 1], or the ISRF from Black et al (1994) [isrf_option = 2].
c The dust opacities used here are from Ossenkopf & Henning (1994; thick ice, 
c low n), extended to short wavelengths using values from Mathis et al (1983).
c
c We assume that the visual/IR/sub-mm background scales linearly with the UV background.
c In the future, we should make these independently configurable.
c 
        field_strength = G0 / 1.7e0
        if (isrf_option .eq. 1) then
          isrf_heat  = 5.8e-24 * chi_isrf * yn * field_strength
        else
          isrf_heat  = 5.73e-24 * chi_isrf * yn * field_strength
        endif 
c
c       LO: add rate of energy density from fervent to the isrf_heat
c           term
        isrf_heat = isrf_heat + Eabs
c
c Iterate using Newton's method
c
        tdust_new = tdust

c        print*,'cool_util: tdust',tdust,temp,yn
        do K = 1, ndust_iter
          tdust = tdust_new
c
c Cooling from dust (assuming optically thin emission) -- computed using same
c opacities as for heating (i.e. OH94, MMP83). Accurate to within a factor two
c over the temperature range 5 < Td < 100 K (==> max 10% error in Td)
c
          dust_cool = 4.68e-31 * tdust**6e0 * yn
c
c Gas-grain energy transfer -- Hollenbach & McKee (1989)
c
#ifdef NO_GAS_DUST_HEAT
          gas_dust  = 0e0
#else
          gas_dust  = cl14 * (temp - tdust) * yn**2.e0
#endif
c
c Heating from H2 formation on grain surface
c
          h2_form_dust = 1e0 - h2_form_ex - h2_form_kin
c
c H2 formation rate on dust -- we calculate this here, as this allows us to
c vary tdust without needing to call cheminmo again
c 
          temp2  = temp * 1e-2
          tdust2 = tdust * 1e-2
          fa     = 1e0 / (1e0 + 1e4 * exp(-6e2 / tdust))
          stick  = 1e0 / (1e0 + 0.4e0 * (temp2 + tdust2)**0.5e0 +
     $                    0.2e0 * temp2 + 0.08e0 * temp2**2e0)
          dfa    = - fa * (1 - fa) * (6e2 / tdust**2e0)
          dtrdust = ch7 * stick * (dfa - 2e-2 * stick * fa / 
     $            sqrt(tdust + temp))
          Rdust_rat = (dfa / fa - 2e-2 * stick / sqrt(tdust + temp))
#ifdef PDRtestChem
          Rdust = 3e-18*temp**0.5
          RH2   = Rdust * abHI * yn**2e0
          Rdust_rat = 0d0
#else
          Rdust  = ch7 * fa * stick
          RH2 = Rdust * abHI * yn**2.e0
#endif
c
#ifdef NO_DUST_H2_HEAT
          dust_h2_heat = 0e0
#else
          dust_h2_heat = (RH2 / dust_to_gas_ratio) * h2_form_dust 
     $                 * 4.48e0 * eV
#endif
c
          if (Rdust .eq. 0e0) then
c This can happen if we've set deff = 0
            F_nr = isrf_heat - dust_cool + gas_dust
#ifdef NO_GAS_DUST_HEAT
            dF_nr = -(6e0 * dust_cool / tdust)
#else
            dF_nr = -(6e0 * dust_cool / tdust) - cl14 * yn**2e0
#endif
          else
            F_nr  = isrf_heat - dust_cool + gas_dust + dust_h2_heat
#ifdef NO_GAS_DUST_HEAT
            dF_nr = -(6e0 * dust_cool / tdust)
     $            + dust_h2_heat * Rdust_rat
#else
            dF_nr = -(6e0 * dust_cool / tdust) - cl14 * yn**2e0 
     $            + dust_h2_heat * Rdust_rat
#endif
          endif
#if 0
c          if (temp .lt. 1.1d0) then
            print*, "Dust debug -- K, F_nr, dF_nr", K, F_nr, dF_nr
            print*, "Dust debug -- isrf_heat", isrf_heat
            print*, "Dust debug -- dust_cool", dust_cool
            print*, "Dust debug -- gas_dust", gas_dust
            print*, "Dust debug -- dust_h2_heat", dust_h2_heat
            print*, "Dust debug -- temp, tdust, yn, abhI,tdust_new", 
     $       temp,tdust, yn, abHI, tdust-F_nr/dF_nr
            print*, "Dust debug -- ch7, cl14, Rdust", ch7, cl14, Rdust
c          endif
#endif

          tdust_new = tdust - F_nr / dF_nr
          if (abs(tdust_new - tdust) .lt. 1e-3) then
            goto 101
          endif
        enddo
c      
        print*, "Dust temp iteration failed! ", tdust, tdust_new, 
     $          tdust-tdust_new
c        print*, "Particle :", id_current
        print*, "density, temp,", yn, temp
        print*, isrf_heat, dust_cool, gas_dust, dust_h2_heat
        print*, "Dust debug -- K, F_nr, dF_nr", K, F_nr, dF_nr
        print*, "Dust debug -- isrf_heat", isrf_heat
        print*, "Dust debug -- dust_cool", dust_cool
        print*, "Dust debug -- gas_dust", gas_dust
        print*, "Dust debug -- dust_h2_heat", dust_h2_heat
        print*, "Dust debug -- temp, tdust, yn, abhI,tdust_new", 
     $       temp,tdust, yn, abHI, tdust-F_nr/dF_nr
        print*, "Dust debug -- ch7, cl14, Rdust", ch7, cl14, Rdust

        ABORT(10)
      else
        RH2   = 0e0
        Rdust = 0e0
        tdust = temp
      endif
c     
 101  continue
c
      if (tdust .eq. 0e0) then
        print*, "Zero dust temp!"
        print*, isrf_heat, dust_cool, gas_dust, dust_h2_heat
      endif

      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////        C A L C _ D U S T _ T E M P        \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            C A L C _ X R A Y              \\\\\\\\\\
c
c=======================================================================
c
c     Computes X-ray heating & ionization rates as a function of H nuclei
c     number density and electron abundance. Based on Appendix A of 
c     Wolfire et al, 1995.
c
      subroutine calc_xray(NH, abe, ion, heat)
      implicit none
#include "param.h"
#include "cool.h"
      REAL NH, abe, ion, heat
      REAL p1, p2, f4, f5, f6

c We set a floor on NH of 1e18 cm^-2, as the fit in Wolfire et al is
c intended for column densities above this value, and at this v. low
c column density, we should, strictly speaking, also be worrying about
c Lyman continuum photons
c
c We set an upper limit on NH of 1e23 cm^-2, to avoid potential 
c numerical problems for larger NH, and because at this value,
c X-rays are ineffective compared to cosmic rays in any case
c
      p1 = log10(NH / 1e18)
      if (p1 .lt. 0e0) then
        p1 = 0e0
      elseif (p1 .gt. 5e0) then
        p1 = 5e0
      endif

      if (abe .lt. 1e-4) then
        p2 = -4.
      elseif (abe .gt. 0.1) then
        p2 = -1.
      else
        p2 = log10(abe)
      endif

      f4 = 1.06  + 4.08e-2 * p2 + 6.51e-3 * p2**2e0
      f5 = 1.90  + 0.678   * p2 + 0.113   * p2**2e0
      f6 = 0.990 - 2.74e-3 * p2 + 1.13e-3 * p2**2e0
      
      ion = f4 * (-15.6 - 1.10 * p1 + 9.13e-2 * p1**2e0) 
     $      + f5 * 0.87 * exp(-((p1 - 0.41) / 0.84)**2e0)

      heat = f6 * (-26.5 - 0.920 * p1 + 5.89e-2 * p1**2e0)
     $     + f6 * 0.96 * exp(-((p1 - 0.38) / 0.87)**2e0)

      ion = 1e1**(ion)
      heat = 1e1**(heat)

      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D    S U B R O U T I N E       //////////
c    //////////            C A L C _ X R A Y              \\\\\\\\\\
c
c=======================================================================


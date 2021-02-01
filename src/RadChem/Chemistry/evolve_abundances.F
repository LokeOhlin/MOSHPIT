c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////     E V O L V E _ A B U N D A N C E S     \\\\\\\\\\
c
c=======================================================================
c
c     Written by S. Glover (AMNH, 2004-2005, AIP 2006-2007)
c
      subroutine evolve_abundances(time, dl, yn, divv, energy, 
     $                             current_redshift, non_eq_abundances,
     $                             fshield_H2, fshield_CO,
     $                             AV_mean, chi_mean,
     $                             dust_temp, phih, hvphih,
     $                             h2phih, h2vphih, uvphih, uvgeom, 
     $                             tphoto, NionH0,NdisH2,dustEabs)
      implicit NONE
#include "param.h"
#include "cool.h"
#include "non_eq.h"
c
      REAL    time, dl, yn, divv, energy, current_redshift
      REAL    non_eq_abundances(num_non_eq_species)
      REAL    dust_temp, phih, hvphih, h2phih, h2vphih, uvphih,
     $        uvgeom, tphoto,
     $        NionH0,NdisH2
      REAL    abe, abhep, abHI
      REAL    yntot, ekn, temp, gamma
      REAL    t_start, t_end, dt_est, deltat, t_in, t_out
      REAL    deltat_old
      REAL    tdust_old
      REAL    dustEabs
      REAL    fshield_dust, fshield_H2, fshield_CO, AV_mean,
     $        chi_mean, fshield_dust_H2
      REAL    RH2, Rdust
      REAL    tion, tion_H, tion_H2, trec, tcool
      REAL    delta_e, heating_per_ion
      integer new_step, ntry
c
      integer nstep
      parameter (nstep = 5000000)
c
      REAL safety(nspec)
      integer used_for_stepsize_control(nspec)
#if CHEMISTRYNETWORK == 4
      DATA safety /1e-1, 1e-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 5
      DATA safety /1e-1, 1e-1, 1e-1, DTCOOL_SCALE_FACTOR, 1e-1, 1e-1/
      DATA used_for_stepsize_control /1, 1, 1, 1, 1, 1/
#endif
c
      integer I, J
c
c DVODE control variables -- NB atol declared in non_eq.h
c
      integer itask, iopt, itol, mf, istate, iflag_dvi
c
c Other DVODE variables
c
      integer lrw, liw
      parameter (lrw   = 2000)
      parameter (liw   = 2000)
      REAL    rpar(nrpar), rwork(lrw)
      integer ipar(nipar), iwork(liw)
      REAL    y(nspec), ydot(nspec), y_old(nspec)
c
      external jac, rate_eq
c
c Sanity check for our input time
c
      if (time .le. 0e0) then
#ifdef DEBUG_EVOLVE
        if (time .lt. -1e0) then
          print*, "Warning: negative timestep in evolve_abundances",
     $             time
        endif
#endif
        return
      endif        
c
c Set redshift
c
      redshift = current_redshift
c
c Initialize abundance vector, treating the internal energy as 
c an additional abundance. 
c
c N.B. We assume that abundances are stored consecutively, with
c no gaps, and that itmp == nspec (i.e. that the energy is the last
c value stored)
c
 
      do I = 1, nspec-3
        y(I) = non_eq_abundances(I)
      enddo
      y(itmp)  = energy
      y(iion)  = NionH0
      y(idis)  = NdisH2

c
c Setup initial parameters for DVODE integrator. These don't vary over 
c the course of the timestep.
c
c Optional input supplied -- see iwork, rwork init. below:
c
      iopt = 1
c
c Evolve ODEs from t -> tout, allowing overshooting & interpolation
c [See DVODE documentation for more details].
c
      itask = 1
c
c Use same relative tolerance for each species, different absolute tolerances
c (since magnitudes may vary by many orders of magnitude). 
c
      itol = 2
c
c Set relative tolerance -- for simplicity we keep this fixed at 0.1% for now  
c
      rtol(1:nspec) = RTOL_FIX !!1e-4
c
c ODE integration method: 22 => BDF with an internally generated Jacobian
c
      mf = 22
c
c Start DVODE in start-up mode (which forces it to initialize all of its
c internal arrays). Subsequent calls within the loop below use istate = 2
c
      istate = 1
c
c Optional input -- increase maximum number of steps that DVODE can take 
c internally without returning an to 20000, as 500 doesn't appear to be 
c enough for portions of our parameter space
c
      do I = 1, lrw      
        rwork(I) = 0e0
        iwork(I) = 0
      enddo
      do I = 5, 10      
        rwork(I) = 0e0
        iwork(I) = 0
      enddo
c      
      iwork(6) = 200000
c
      fshield_dust    = exp(-2.5e0 * AV_mean)
      fshield_dust_H2 = exp(-3.5e0 * AV_mean)
c
c Parameters to pass to rate_eq
c
      rpar(1) = yn
      rpar(2) = dl
      rpar(3) = divv
      rpar(4) = fshield_H2
      rpar(5) = fshield_CO
      rpar(6) = fshield_dust  ! Scaling used for CO, PE effect etc.
      rpar(7) = chi_mean
      rpar(8) = fshield_dust_H2
! from radiation 
      rpar(9)  = phih     ! ionisation rate
      rpar(10) = hvphih   ! ionisation heating rate
      rpar(11) = h2phih   ! direct h2 dissociation rate
      rpar(12) = h2vphih  ! direct h2 dissociation heating rate
      rpar(13) = uvphih   ! FUV dissociation rate
      rpar(14) = uvgeom   ! FUV geometric dissociation rate
      rpar(15) = DustEabs ! Dust energy density absorption
c Compute photoion timescale
      ipar(2) = 0
      tphoto  = 1e18
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
      tion_H  = 1e18
      tion_H2 = 1e18
      tion    = 1e18
c
      if (phih .gt. 0.0) then
        tion_H = 1.0 / phih
      endif
c
      if (h2phih .gt. 0.0) then
        tion_H2 = 1.0 / h2phih
      endif
c
      tion = min(tion_H, tion_H2)
c
      if (use_photo_eqb_approx .eq. 1 .and. tion .lt. 1e-4 * time) then
        trec  = 5e12 / yn
        tcool = 1e10 / yn
        if (tion .lt. 1e-4 * trec) then
          if (phih .gt. 0.0) then
            heating_per_ion = hvphih / phih
          else
            heating_per_ion = 0.0
          endif
c Assume we get similar amounts of heating per H, H2 photoion
          delta_e = heating_per_ion * (1.0 - y(ihp)) * yn
          if (delta_e .le. 0.01 * y(itmp) .or. 
     $           tion .le. 0.01 * tcool) then
c Can safely assume this cell is in photoionization equilibrium
            ipar(2) = 1
            y(itmp) = y(itmp) + delta_e   ! Heating due to initial photoion.
            y(ih2)  = 0.0
            y(ihp)  = 1.0
            y(ico)  = 0.0
          endif
        endif
      endif

c Compute fastest timescale on which radiation can change H abundance
      if (uvphih .gt. 0.0) then
        tphoto = min(tion, 1 / uvphih)
      else
        tphoto = tion
      endif

#endif /* CHEMISTRYNETWORK == 4 || 5 */
c
      ipar(1) = 0
      ipar(3) = 0
c
      t_start  = 0e0
      t_end    = time
      new_step = 1
      ntry     = 1
c
      do I = 1, nstep
c
c If we're starting a new step, then figure out a reasonable timestep to
c pass to DVODE -- we need to do this here to avoid DVODE producing -ve
c abundances or -ve temperatures, which can occur if it starts with too
c large a timestep.
c
c If we're not starting a new step, then dt has already been set. 
c
        if (new_step .eq. 1) then
c
c First, save a copy of the current abundances
c
          do J = 1, nspec
            y_old(J) = y(J)
          enddo
          tdust_old = tdust
c
c Now compute ydot. On the first step, we must use rate_eq directly,
c but on subsequent steps, we can just call dvindy [if we're using
c the dvode integrator]
c
          if (I .eq. 1) then
c            print*,'Calling rate_eq in step ',I
c            print*,nspec,t_start,y,ydot,rpar,ipar
            call rate_eq(nspec, t_start, y, ydot, rpar, ipar)
          else
c            print*,'Calling dvindy in step ',I
c            print*,t_start,rwork(21),nspec,ydot,iflag_dvi
            call dvindy(t_start, 1, rwork(21), nspec, ydot, iflag_dvi)
          endif
c
c Init. dt_est to an unrealistic value
c
          dt_est = 1e18

          do J = 1, nspec
            if (used_for_stepsize_control(J) .eq. 1) then
              if (ydot(J)   .ne. 0e0 .and.
     $            safety(J) .ne. 0e0 .and.
     $            y(J)      .gt. ATOL(J)) then
                dt_est = min(dt_est, abs(safety(J) * y(J) / ydot(J)))
              endif
            endif
          enddo
c
c Don't evolve past the end of the global timestep
c
          deltat = min(dt_est, t_end - t_start)
          deltat_old=deltat
        endif
c
        t_in  = t_start
        t_out = t_start + deltat
        if (deltat .lt. 0e0) then
          print*, 'Negative timestep!'
          ABORT(1)
        endif
        if (deltat .eq. 0e0) then 
          if (t_end .eq. t_start) then
            return
          else
            print*, 'Zero timestep!'
            ABORT(2)
          endif
        endif
        call dvode(rate_eq, nspec, y, t_in, t_out, itol,  
     $             itask, istate, iopt, rwork, lrw, iwork, liw, jac, 
     $             mf, rpar, ipar)
c Ensure that dvode has produced sane output values -- if not, we
c need to shrink the timestep. Note that if ipar(1) is non-zero, 
c or if istate != 2 (indicating an integrator error), we're going
c to shrink the timestep anyway, so we can skip the check.
c
        if (ipar(1) .eq. 0 .and. istate .eq. 2) then
          call check_abundances(y, ipar)
        endif

        if ((ipar(1) .ne. 0) .or. (istate .ne. 2)) then 
#ifdef DEBUG_EVOLVE
         if (istate .ne. 2) then
           print*, "Error: istate = ", istate
         endif
#endif
          do J = 1, nspec
            y(J) = y_old(J)
          enddo
          tdust    = tdust_old
          deltat       = deltat / 2e0
          if (deltat .lt. 1e3 .and. ntry .ge. 64) then
            print*, "Excessively small timestep! ", deltat, dt_est, 
     $              t_start, t_end
            print*, "rpar ", rpar
            print*, "y ", y
            print*, "istate", istate
	      print*, "ydot ", ydot
            print*, "ATOL ", ATOL(1), ATOL(2), ATOL(3), ATOL(4)
#ifdef TREE_RAD
            print*, "Tree rad totals: ", column_density_projection
#endif
#ifdef TREE_RAD_H2
            print*, "Tree rad h2: ", column_density_projection_h2
#endif
#ifdef TREE_RAD_CO
            print*, "Tree rad co: ", column_density_projection_co
#endif
            ABORT(100)
          endif
          new_step = 0
          ipar(1)  = 0
          istate   = 1
          ntry     = ntry + 1
        else
          t_start = t_start + deltat
          new_step = 1
          ntry     = 1
        endif
c

c        print*,'t_start, dt, t_end', t_start,deltat,t_end

        if (t_start .ge. t_end) then
          goto 100
        endif
      enddo
c
      write(*,'(A)') 'EVOLVE_ABUNDANCES: nstep exceeded,' 
      print*, deltat
      print*, rpar
      print*, non_eq_abundances
      print*, y
      print*, ydot
      ABORT(0)
      return

100   continue

      do I = 1, nspec-3
        non_eq_abundances(I) = y(I)
      enddo
      energy = y(itmp)
      NionH0 = y(iion)
      NdisH2 = y(idis)

c Calculate final dust temp.
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
      abHI  = max(0e0, 1e0 - 2e0 * y(ih2) - y(ihp))
#endif
      abe   = y(ihp) !+ abhep ! this appears to be a bug since abhep is never set  DS + SG 2014
      yntot = (1e0 + abhe - y(ih2) + abe) * yn
      temp  = 2e0 * energy / (3e0 * yntot * kboltz)
      call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust,
     $                    DustEabs, 0)
      dust_temp = tdust
c


      return      
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////     E V O L V E _ A B U N D A N C E S     \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      C H E C K _ A B U N D A N C E S    \\\\\\\\\\
c
c=======================================================================
c
      subroutine check_abundances(y, ipar)
      implicit NONE
#include "param.h"
#include "cool.h"
#include "non_eq.h"
      REAL y(nspec), max_abundance(nspec)
      REAL current_carbon, current_oxygen
      REAL diff
      integer ipar(nipar)
      integer J
c
      if (nchem_network .eq. 4) then
        max_abundance(ih2)   = 0.5e0
        max_abundance(ihp)   = 1.0e0
      endif
c
      if (nchem_network .eq. 5) then
        max_abundance(ih2)   = 0.5e0
        max_abundance(ihp)   = 1.0e0
        max_abundance(ico)   = min(abundc, abundo)
      endif
c
c Check for NaNs
c
      do J = 1, nspec
        if (y(J) .ne. y(J)) then
          ipar(1) = 1
          return
        endif
      enddo
c
c Check for negative abundances. Small negative values (which have
c absolute values less than the corresponding absolute tolerance)
c are accepted, but are clipped to zero. For large negative values,
c we repeat the step with a smaller dt.
c
      do J = 1, nspec
        if (y(J) .lt. 0e0) then
          if (y(J) .lt. -ATOL(J)) then
            ipar(1) = 1
            return
          else
            y(J) = 0e0
          endif
        endif
      enddo
c
c Don't allow ion abundances to exceed corresponding total metal
c abundances. Again, small errors are silently corrected, large
c errors lead to us repeating with a smaller timestep
c
      do J = 1, nspec
         if ((J .ne. itmp .and. J .ne. iion .and. J .ne. idis) .and.
     $        y(J).gt. max_abundance(J)) then
          diff = y(J) - max_abundance(J)
          if (diff .lt. 1e-4 * y(J)) then
            y(J) = max_abundance(J)
          else
#ifdef DEBUG_EVOLVE
            print*, "Excess abundance", y(J), max_abundance(J)
#endif
            ipar(1) = 1
            return
          endif
        endif
      enddo
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////      C H E C K _ A B U N D A N C E S    \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////       I N I T _ T O L E R A N C E S       \\\\\\\\\\
c
c=======================================================================
c
      subroutine init_tolerances
      implicit NONE
#include "param.h"
#include "non_eq.h"
#include "cool.h"
c
c Set-up absolute tolerances
c
      if (nchem_network .eq. 4) then
        ATOL(ihp) = ATOL_HP
        ATOL(ih2) = ATOL_H2
      endif
c
      if (nchem_network .eq. 5) then
        ATOL(ihp) = ATOL_HP
        ATOL(ih2) = ATOL_H2
        ATOL(ico) = ATOL_CO
      endif
c
      ATOL(itmp)  = ATOL_TMP
      ATOL(iion)  = ATOL_ION
      ATOL(idis)  = ATOL_DIS
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////       I N I T _ T O L E R A N C E S       \\\\\\\\\\
c
c=======================================================================

!!****
!!
!! NAME
!!
!!  Chemistry_FortranInit
!!
!!
!! SYNOPSIS
!!
!!  Chemistry_FortranInit()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once first. 
!!
!!***

subroutine Chemistry_FortranInit() bind(C, name="Chemistry_FortranInit")
implicit none

INTERFACE
  subroutine getRealChemistryPar(nam, val) bind(c)
    use iso_c_binding, only: c_char, c_double
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: nam
    real(kind = c_double), intent(out) :: val
  end subroutine
  subroutine getIntegerChemistryPar(nam, val) bind(c)
    use iso_c_binding, only: c_char, c_int
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: nam
    integer(kind = c_int), intent(out) :: val
  end subroutine
END INTERFACE

#include "cool.h"

!  integer,intent(IN) :: myPe
  character(len = 15) :: numToStr
  character(len = 80) :: strBuff


  strBuff = "ch_deff"
  call getRealChemistryPar(strBuff, deff)
  strBuff = "ch_abundC"
  call getRealChemistryPar(strBuff, abundc)
  strBuff = "ch_abundO"
  call getRealChemistryPar(strBuff, abundo)
  strBuff = "ch_abundSi"
  call getRealChemistryPar(strBuff, abundsi)
  strBuff = "ch_abundD"
  call getRealChemistryPar(strBuff, abundD)
  strBuff = "ch_abundM"
  call getRealChemistryPar(strBuff, abundM)
  strBuff = "ch_abundN"
  call getRealChemistryPar(strBuff, abundN)
  strBuff = "ch_G0"
  call getRealChemistryPar(strBuff, G0)
  strBuff = "ch_f_rsc"
  call getRealChemistryPar(strBuff, f_rsc)
  strBuff = "ch_phi_pah"
  call getRealChemistryPar(strBuff, phi_pah)
  strBuff = "ch_tdust"
  call getRealChemistryPar(strBuff, tdust)
  strBuff = "ch_dust_to_gas_ratio"
  call getRealChemistryPar(strBuff, dust_to_gas_ratio)
  strBuff = "ch_AV_conversion_factor"
  call getRealChemistryPar(strBuff, AV_conversion_factor)
  strBuff = "ch_cosmic_ray_ion_rate"
  call getRealChemistryPar(strBuff, cosmic_ray_ion_rate)
  strBuff = "ch_NH_ext"
  call getRealChemistryPar(strBuff, NH_ext)
  strBuff = "ch_h2_form_ex"
  call getRealChemistryPar(strBuff, h2_form_ex)
  strBuff = "ch_h2_form_kin"
  call getRealChemistryPar(strBuff, h2_form_kin)
  strBuff = "ch_xray_scaling"
  call getRealChemistryPar(strBuff, xray_scaling)
  strBuff = "ch_Z_atom"
  call getRealChemistryPar(strBuff, Z_atom)

! Set AV_ext based on specificed NH_ext, conversion factor and dust:gas ratio
  AV_ext = AV_conversion_factor * NH_ext * dust_to_gas_ratio

  strBuff = "ch_iphoto"
  call getIntegerChemistryPar(strBuff, iphoto)
  strBuff = "ch_iflag_mn"
  call getIntegerChemistryPar(strBuff, iflag_mn)
  strBuff = "ch_iflag_ad"
  call getIntegerChemistryPar(strBuff, iflag_ad)
  strBuff = "ch_iflag_atom"
  call getIntegerChemistryPar(strBuff, iflag_atom)
  strBuff = "ch_iflag_3bh2a"
  call getIntegerChemistryPar(strBuff, iflag_3bh2a)
  strBuff = "ch_iflag_3bh2b"
  call getIntegerChemistryPar(strBuff, iflag_3bh2b)
  strBuff = "ch_iflag_h3pra"
  call getIntegerChemistryPar(strBuff, iflag_h3pra)
  strBuff = "ch_isrf_option"
  call getIntegerChemistryPar(strBuff, isrf_option)
  strBuff = "ch_no_chem"
  call getIntegerChemistryPar(strBuff, no_chem)
  strBuff = "ch_use_photo_eqb_approx"
  call getIntegerChemistryPar(strBuff, use_photo_eqb_approx)


!
  call coolinmo
  call cheminmo
  call init_tolerances

 
  return
end subroutine Chemistry_FortranInit


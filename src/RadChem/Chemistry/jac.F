      SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)

      IMPLICIT NONE
#include "param.h"
      integer NEQ, ML, MU, NROWPD, IPAR
      REAL T, Y(NEQ), PD(NROWPD,NEQ), RPAR
      return
      end

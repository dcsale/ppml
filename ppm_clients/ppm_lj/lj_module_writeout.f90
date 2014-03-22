#define __REAL    1
#define __INTEGER 2
#define __RANK    3

  MODULE lj_module_writeout

    INTERFACE lj_writeout
       MODULE PROCEDURE lj_writeout_real
       MODULE PROCEDURE lj_writeout_integer
       MODULE PROCEDURE lj_writeout_rank
    END INTERFACE

  CONTAINS
        
#define __TYPE  __REAL
#include "lj_writeout.f90"
#undef  __TYPE
#define __TYPE  __INTEGER
#include "lj_writeout.f90"
#undef  __TYPE
#define __TYPE  __RANK
#include "lj_writeout.f90"
#undef  __TYPE

#undef __REAL
#undef __INTEGER
#undef __RANK

END MODULE lj_module_writeout

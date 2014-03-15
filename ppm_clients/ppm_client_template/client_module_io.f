      !-------------------------------------------------------------------------
      !  Module       :                   client_module_io
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all subroutines needed for 
      !                 file I/O.
      !                
      !  Remarks      : 
      !
      !  References   :
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE client_module_io
         !----------------------------------------------------------------------
         !  Define interface to client_write_output
         !----------------------------------------------------------------------
         INTERFACE client_write_output
             MODULE PROCEDURE client_write_output
         END INTERFACE

         !----------------------------------------------------------------------
         !  Define interface to client_write_diag
         !----------------------------------------------------------------------
         INTERFACE client_write_diag
             MODULE PROCEDURE client_write_diag
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS
#include "client_write_output.f"
#include "client_write_diag.f"

      END MODULE client_module_io

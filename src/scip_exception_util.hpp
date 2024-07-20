#pragma once

#include <iostream>
#include "scip_exception.hpp"

/** calls a SCIP function and prints an error message if anything went wrong
 *  this is required in place of SCIP_CALL_EXC in destructors as they are declared noexcept
 *
 */
#define SCIP_CALL_NOEXC(x) { \
   SCIP_RETCODE retcode; \
   if ((retcode = (x)) != SCIP_OKAY) { \
      char _msg[SCIP_MSG_MAX]; \
      if (SCIPgetErrorString(retcode, _msg, SCIP_MSG_MAX) == NULL) { \
         std::cerr << "SCIP_CALL_NOEXC unknown SCIP retcode " << (int) retcode << std::endl; \
      } else { \
         std::cerr << "SCIP_CALL_NOEXC " << _msg << std::endl; \
      } \
   } \
}

/* Umbrella header pulled in by bindgen.
   Add or remove solver headers as needed. */

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

/* ODE solvers */
#include <cvode/cvode.h>
#include <cvodes/cvodes.h>

/* DAE solvers */
#include <ida/ida.h>
#include <idas/idas.h>

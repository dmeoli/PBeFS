/*--------------------------------------------------------------------------*/
/*------------------------------ File core.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file that includes the core of PBeFS, i.e., what does not depend
 * on the kind of search: the helpers and the parameters of the parallel
 * searches [see SearchUtils.h], and the observers of the searches [see
 * SearchObserver.h and TraceRecorder.h]. The workers run on the executors
 * of ADS and exchange their messages through its mailboxes, those of the
 * processes of MPI included.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ DEFINITIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __PBeFS_core
 #define __PBeFS_core
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- INCLUDES --------------------------------*/
/*--------------------------------------------------------------------------*/

#include "PBeFS/core/SearchObserver.h"

#include "PBeFS/core/SearchUtils.h"

#include "PBeFS/core/TraceRecorder.h"

/*--------------------------------------------------------------------------*/

#endif  /* __PBeFS_core */

/*--------------------------------------------------------------------------*/
/*---------------------------- End File core.h -----------------------------*/
/*--------------------------------------------------------------------------*/

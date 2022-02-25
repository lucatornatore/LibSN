#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libimf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#ifndef LIBSN_VARS_H
#define LIBSN_VARS_H

#ifndef ALIGN
#define ALIGN 8
#endif

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

#define IA  0
#define II  1
#define AGB 2

#define DEFAULT_MEM_AMOUNT (1.5 * 1024*1024)
typedef long long unsigned int mem_t;
#define MEM_T_SIZE sizeof(mem_t)

extern char *options_string;

extern char  *Memory, *Memory_roof, *Memory_top;
extern size_t Memory_amount;

extern char *LSN_SnNAMES[LSN_N_SnSPECIES];

extern int   lsn_use[LSN_N_SnSPECIES]; 

extern double lsn_EjError;

extern int ThisTask;
extern int Task;

#define dprintf(LEVEL, TASK, ...) do{if( ((LEVEL) <= lsn_verbose_level) && (ThisTask == (TASK))) fprintf(stdout, __VA_ARGS__);} while(1 == 0)
#define dprintf_err(LEVEL, TASK, ...) do{if( ((LEVEL) <= lsn_verbose_level) && (ThisTask == (TASK))) fprintf(stderr, __VA_ARGS__);} while(1 == 0)
#define VDBG  4   // verbose level for debug
#define VDIAG 2   // verbose level for diagnostics
#define VMSG  1   // verbose level for flow messages
#define VXX   0   // essential messages
#define VERR  VXX // non letal errors
#define VXERR -1  // letal errors



#define DONT_WARN_UNUSED(x) (void)(x)   // __attribute__((unused)) works for gcc, but this way is more portable


typedef struct {
  int   IMFi;
  void *params;
} imf_plus_param_type;

/* /------------------------------------------------\  
 * |                                                |
 * | LIFETIME function                              |
 * |                                                |
 * \------------------------------------------------/ */

#define PADOVANI_MATTEUCCI_1993 0
#define MAEDER_MEYNET_1989 1

extern int lsn_LifeTime;
extern const char *lsn_lifetime_name[];


/* /------------------------------------------------\  
 * |                                                |
 * | YIELDS declaration                             |
 * |                                                |
 * \------------------------------------------------/ */

/* How many yield sets have been defined */
extern int Nset_ofYields[LSN_N_SnSPECIES];

/*  */
extern int extrap_sup[LSN_N_SnSPECIES];
extern int extrap_inf[LSN_N_SnSPECIES];

/* The length of Z-dimension for each set */
extern int      *Zbins_dim[LSN_N_SnSPECIES];
/* The values of Z bins for each set */
extern double  **Zbins[LSN_N_SnSPECIES];

/* The length of M-dimension for each set for each Zbin */
extern int     **Mbins_dim[LSN_N_SnSPECIES];
extern int      *Mbins_maxdim[LSN_N_SnSPECIES];
extern int      *Mbins_irr[LSN_N_SnSPECIES];

/* The values of M bins for each set for each Zbins*/
extern double ***Mbins[LSN_N_SnSPECIES];

/* The actual tables */
extern double ***Y[LSN_N_SnSPECIES], **Ej[LSN_N_SnSPECIES];

/* non proc flags */
extern int      *NonProcOn[LSN_N_SnSPECIES];

/* The filename base (if more than one set is present, the
 * corresponding files must be named $filename.000,
 * $filename.001 and so on.
 */
extern char     *string_space;
#define TOT_STR_SPACE        1000
#define MET_FNAME_START         0
#define IA_FNAME_START        100
#define II_FNAME_START        400
#define AGB_FNAME_START       700
extern char     *Yields_filename[LSN_N_SnSPECIES];

#define lsn_MAX(x, y) ( ((x) > (y))? (x) : (y) )
#define lsn_MIN(x, y) ( ((x) < (y))? (x) : (y) )


/* dimension of the integration space */
#define lsn_gsl_intspace_dim 10000
/* specify what function from gsl qag should be used */
#define lsn_gsl_intkey 6

#ifdef USE_QAGS
#define LSN_INT( F, INF, SUP, ABS, REL, LIMIT, KEY, WSPACE, RES, ERR) gsl_integration_qags((F), (INF), (SUP), (ABS), (REL), (LIMIT), (WSPACE), (RES), (ERR) )
#else
#define LSN_INT( F, INF, SUP, ABS, REL, LIMIT, KEY, WSPACE, RES, ERR) gsl_integration_qag((F), (INF), (SUP), (ABS), (REL), (LIMIT), (KEY), (WSPACE), (RES), (ERR) )
#endif


extern gsl_function lsn_F;
extern gsl_integration_workspace *lsn_w;

extern lsn_SDtype SD;

extern double *lsn_initial_metallicity;
extern double lsn_initial_abundance;

extern char *metals_file;

extern char lsn_prefix[10];
extern char lsn_err_prefix[20];
extern char lsn_warn_prefix[20];

extern char *lsn_error_msgs[];

extern int STELLAR_EVOLUTION_IMPOSSIBLE;



/* /------------------------------------------------\  
 * | U T I L I T I E S                              |
 * |                                                |
 * \------------------------------------------------/ */

int   get_parameters(char*);
int   read_yields   (int);
int   read_elements (char*);


void *manage_memory(size_t);
void *_allocate    (size_t, int);
void *_callocate   (size_t, int);
void *_reallocate  (char*, size_t, int);
void *_release     (char*);

#define allocate( ARG ) _allocate( (ARG), ALIGN)
#define callocate( ARG ) _callocate( (ARG), ALIGN)
#define release( ARG ) _release((char*)(ARG))
#define reallocate( ARG1, ARG2) _reallocate((char*)(ARG1), (ARG2), ALIGN)


/*  ===========================================================================  */



#define SET_ERR_MSG( ERRCODE) do { lsn_error = (ERRCODE); sprintf(lsn_error_message, "%s ", lsn_error_msgs[(ERRCODE)]);} while(1 == 0)

#define ADD_ERR_MSG( ERRCODE ) do { lsn_error = (ERRCODE); if(strlen(lsn_error_message) > 0) sprintf(&lsn_error_message[strlen(lsn_error_message)], "\n%s ", lsn_error_msgs[(ERRCODE)]); else sprintf(lsn_error_message, "%s ", lsn_error_msgs[(ERRCODE)]);} while(1 == 0)

#define ADD_toERR_MSG( STRING ) do { sprintf(&lsn_error_message[strlen(lsn_error_message)], "%s ", STRING); } while(1 == 0)


#define SET_ERROR( ERRCODE ) lsn_error = (ERRCODE); lsn_error_message[0]='\0'

#define NO_METALS_DEFINED 1

#define NO_YIELDS_DEFINED 2

#endif

#include "libsn_errors.h"

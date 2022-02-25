#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libimf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#ifndef LIBSN_H
#define LIBSN_H

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

#define LSN_N_SnSPECIES 3

#define lsn_SnII 2
#define lsn_SnIa 1
#define lsn_AGB  4
#define lsn_ALL  (lsn_SnII + lsn_SnIa + lsn_AGB)

#define PADOVANI_MATTEUCCI_1993 0
#define MAEDER_MEYNET_1989 1

#define BYTIME 0
#define BYMASS 1

#define LSN_OK 0

extern int lsn_verbose_level;

/* /------------------------------------------------\  
 * |                                                |
 * | LIFETIME function                              |
 * |                                                |
 * \------------------------------------------------/ */


extern       int   lsn_LifeTime;
extern const char *lsn_lifetime_name[];


/* /------------------------------------------------\  
 * |                                                |
 * | INTEGRATION stuffs                             |
 * |                                                |
 * \------------------------------------------------/ */


/* set by default ; it means that calculations are performed by mass */
extern int lsn_UseMassIntegration;   

typedef struct
/* 
 * this structure carries all the informations
 * between different routines
 */
{
  int Zbin, Zdim, Mdim, Mdim1, Yset, Mirr, Mmaxdim, extrap, extrap_inf;
  double Zstar;
  double *ZArray, *MArray, *MArray1;
  double *Y;
} lsn_SDtype;

typedef double (EXTERNALF)(double, void *);

/* /------------------------------------------------\  
 * |                                                |
 * | METALS SPECIES declaration                     |
 * |                                                |
 * \------------------------------------------------/ */


extern int lsn_NMet;                       /* total number of metals */
extern int lsn_Hyd, lsn_Hel, lsn_FillEl;   /* position of significative metals  */
extern int lsn_Extr;                       /*  */ 


extern char **lsn_metal_names;


/* /------------------------------------------------\  
 * |                                                |
 * | YIELDS declaration                             |
 * |                                                |
 * \------------------------------------------------/ */

extern double lsn_Mup, lsn_iMi, lsn_iMU;
extern double lsn_mean_lifetime, lsn_iinf_lifetime, lsn_isup_lifetime;

extern int    lsn_IaMode;
extern double lsn_BinFrac, lsn_BinFrac_normalized, lsn_k_alpha, lsn_MBm, lsn_MBM;

extern char   lsn_error_message[2000];

extern int    lsn_error;

/*  ===========================================================================  */
/*  PROTOTYPES  */


/* /------------------------------------------------\  
 * | APIs to access tables meta-data and data       |
 * |                                                |
 * \------------------------------------------------/ */


int      NUM_of_ZBINS (int, int);              // returns the number of metallicity bins for a given set
double  *ZBINS        (int, int);              // returns the metallicity bins for a given set
int     *MBINS_dim    (int, int);              // returns the start of Mbins dim for a given set
int      NUM_of_MBINS (int, int, int);         // returns the number of mass bins for a given zbin in a set
double **MBINS        (int, int);              // returns the start of mass bins for a given set
double  *MBINS_ZBIN   (int, int, int);         // returns the mass bins for a given zbin in a set
double  *YIELDS_all   (int, int, int);         // returns the start of yields of an element for a in a set
double  *YIELDS       (int, int, int, int);    // returns the start of yields of an element for a given zbins in a set



/* /------------------------------------------------\  
 * | Calculus                                       |
 * |                                                |
 * | > metals_integration                           |
 * | > num_time_integration                         |
 * \------------------------------------------------/ */



int    lsn_initialize_libsn(char *, char *);                   // initialize the library

void   lsn_finalize_libsn(void);                               // initialize the library

void   lsn_get_libsn_options(char **);                         // print the compilation options in the given char pointer

void   lsn_set_initial_metallicity(int, ...);                  // set initial metallicity of star

int    lsn_get_metal_position( char * );

int    lsn_get_metals(int, int, int, double,               // main integration routine
		      double, double, int, int*, double*);

double lsn_get_sn_number(int, int, double, double);        // main integration routine

void   lsn_dump_yields(int, int);                          // dump internal yields tables



double lsn_dm_dt(double, double);

double lsn_lifetime(double);

double lsn_dying_mass(double);

double lsn_nRSnII(double, void*);  

double lsn_mRSnII(double, void*);

double lsn_zmRSnII(double, void*);

double lsn_ztRSnII(double, void*);

double lsn_ejectaSnII(double, void*);

extern double (* lsn_nRSnIa)(double, void*);

extern double (* lsn_DTDfunction)(double, void*);

int    make_time_interval_bins_by_events(int, int, double, double, double,
					 double *, double *, double *);




#endif

#include "libsn.h"
#include "libsn_vars.h"


char  *Memory = NULL;
char  *Memory_roof = NULL;
char  *Memory_top = NULL;
size_t Memory_amount = DEFAULT_MEM_AMOUNT;

char *LSN_SnNAMES[LSN_N_SnSPECIES] = {"SnIa", "SnII", "LIMS"};

int lsn_verbose_level;
int ThisTask = 0;

/* /------------------------------------------------\  
 * |                                                |
 * | LIFETIME function                              |
 * |                                                |
 * \------------------------------------------------/ */

#define PADOVANI_MATTEUCCI_1993 0
#define MAEDER_MEYNET_1989 1

int lsn_LifeTime;
const char *lsn_lifetime_name[] = {"Padovani & Matteucci (1993)",
				   "Maeder & Meynet (1989)",
				   "Other"};


/* /------------------------------------------------\  
 * |                                                |
 * | INTEGRATION stuffs                             |
 * |                                                |
 * \------------------------------------------------/ */

/* set by default ; it means that calculations are performed by mass */
int lsn_UseMassIntegration;   

double *lsn_initial_metallicity;
double  lsn_initial_abundance;

/* /------------------------------------------------\  
 * |                                                |
 * | METALS SPECIES declaration                     |
 * |                                                |
 * \------------------------------------------------/ */


int lsn_NMet;                       /* total number of metals */
int lsn_Hyd, lsn_Hel, lsn_FillEl;   /* position of significative metals  */

char **lsn_metal_names;


/* /------------------------------------------------\  
 * |                                                |
 * | YIELDS declaration                             |
 * |                                                |
 * \------------------------------------------------/ */

int lsn_use[LSN_N_SnSPECIES] = {0};

double lsn_EjError;

/* How many yield sets have been defined */
int Nset_ofYields[LSN_N_SnSPECIES];

/*  */
int extrap_sup[LSN_N_SnSPECIES];
int extrap_inf[LSN_N_SnSPECIES];

/* The length of Z-dimension for each set */
int *Zbins_dim[LSN_N_SnSPECIES];
int *Mbins_maxdim[LSN_N_SnSPECIES];
int *Mbins_irr[LSN_N_SnSPECIES];
/* The values of Z bins for each set */
double **Zbins[LSN_N_SnSPECIES];

/* The length of M-dimension for each set */
int **Mbins_dim[LSN_N_SnSPECIES];
/* The values of M bins for each set */
double ***Mbins[LSN_N_SnSPECIES];

/* The actual tables */
double ***Y[LSN_N_SnSPECIES], **Ej[LSN_N_SnSPECIES];

int *NonProcOn[LSN_N_SnSPECIES];

/* The filename base (if more than one set is present, the
 * corresponding files must be named $filename.000,
 * $filename.001 and so on.
 */
char *string_space;
char *Yields_filename[LSN_N_SnSPECIES];


int    lsn_IaMode = -1;
double lsn_BinFrac, lsn_BinFrac_normalized, lsn_k_alpha;
double lsn_MBm, lsn_MBM, lsn_Mup;
double lsn_iMi, lsn_iMU;
double lsn_mean_lifetime, lsn_iinf_lifetime, lsn_isup_lifetime;



/* ===================================================================================  */

gsl_function lsn_F;
gsl_integration_workspace *lsn_w;

lsn_SDtype SD;

char *metals_file;

char lsn_error_message[2000];
char lsn_prefix[10] = "[lsn] ";
char lsn_err_prefix[20]  = "[lsn] >> err:: " ;
char lsn_warn_prefix[20] = "[lsn] >> warn:: ";

int  lsn_error;

int  STELLAR_EVOLUTION_IMPOSSIBLE = 0;

/* ===================================================================================  */


char *lsn_error_msgs[ERR_LAST] =
  {"[lsn] everything went on smoothly",
   "[lsn] unable to open paramfile",
   "[lsn] some fundamental parameters are missed in paramfile",
   "[lsn] not enough memory to allocate",
   
   "[lsn] some error while reading the IMF file",
   "[lsn] some error in IMF",
   "[lsn] it was impossible to allocate gsl integration space",
   "[lsn] some error while initializing the library",

   "[lsn] some severe error in integration routines",
   
   "[lsn] letal error occurred: it is impossible to proceed",        // letal errors

   "[lsn] some parameters are missed in paramfile",
   
   "[lsn] some error in reading the elements file",

   "[lsn] unable to open one (or more) Yields file(s)",
   "[lsn] some error while reading Yields files",
   "[lsn] some severe error in the Ej",
   "[lsn] some error in the Yields' mass bins",
   "[lsn] stellar evolution calculation not possible",
   "[lsn] stellar evolution calculation not possible: no metals are define",
   "[lsn] stellar evolution calculation not possible: no yields are defined",
   "[lsn] stellar evolution calculation not possible: neither metals nor yields are defined",
   "[lsn] stellar evolution calculation not possible: no metals are defined and so no yields have been loaded",

   "[lsn] severe error in stellar evolution: impossible to calculate", // severe SEv errors

   "[lsn] some mild error in the Ej",
   "[lsn] gsl reported some error in integration routine",
   "[lsn] some index was out-of-bounds",
   "[lsn] bad limits: sup was smaller than inf"
  };


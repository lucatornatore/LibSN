
#define ERR_PARAM_FILE_OPEN                    1
#define ERR_PARAM_FILE_NEEDED_PARAMS_MISS      ( ERR_PARAM_FILE_OPEN        + 1 )
#define ERR_MEMORY                             ( ERR_PARAM_FILE_NEEDED_PARAMS_MISS + 1 )
#define ERR_IMF_FILE                           ( ERR_MEMORY                 + 1 )
#define ERR_IMF                                ( ERR_IMF_FILE               + 1 )

#define ERR_GSL_ALLOC                          ( ERR_IMF                    + 1 )

#define ERR_INIT                               ( ERR_GSL_ALLOC              + 1 )

#define ERR_INTEGRATION                        ( ERR_INIT + 1 )
  
#define LETAL_ERROR                            ( ERR_INTEGRATION + 1 )

#define ERR_PARAM_FILE_PARAMS_MISS             ( LETAL_ERROR     + 1 )

#define ERR_METALS_FILE                        ( ERR_PARAM_FILE_PARAMS_MISS + 1 )
#define ERR_YIELDS_FILE_OPEN                   ( ERR_METALS_FILE            + 1 )
#define ERR_YIELDS_FILE                        ( ERR_YIELDS_FILE_OPEN       + 1 )
#define ERR_YIELDS_Ej_SEVERE                   ( ERR_YIELDS_FILE            + 1 )
#define ERR_YIELDS_MASSBINS                    ( ERR_YIELDS_Ej_SEVERE       + 1 )

#define ERR_SEv_IMPOSSIBLE                     ( ERR_YIELDS_MASSBINS         + 1 )
#define ERR_SEv_IMPOSSIBLE_noMETALS            ( ERR_SEv_IMPOSSIBLE          + 1 )
#define ERR_SEv_IMPOSSIBLE_noYIELDS            ( ERR_SEv_IMPOSSIBLE_noMETALS + 1 )
#define ERR_SEv_IMPOSSIBLE_noMETALS_noYIELDS   ( ERR_SEv_IMPOSSIBLE_noYIELDS + 1 )
#define ERR_SEv_IMPOSSIBLE_noMETALS_noYIELDS_loaded   ( ERR_SEv_IMPOSSIBLE_noMETALS_noYIELDS + 1 )

#define SEv_SEVERE_ERROR                       ( ERR_SEv_IMPOSSIBLE_noMETALS_noYIELDS_loaded + 1 )

#define ERR_YIELDS_Ej                          ( SEv_SEVERE_ERROR           + 1 )

#define ERR_GSL_INTEGRATION                    ( ERR_YIELDS_Ej     + 1)

#define ERR_BOUNDS                             ( ERR_GSL_INTEGRATION   + 1 )

#define ERR_LIMITS_INFSUP                      ( ERR_BOUNDS        + 1 )

#define ERR_LAST                               ( ERR_LIMITS_INFSUP + 1 )

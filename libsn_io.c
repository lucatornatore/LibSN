#define _GNU_SOURCE  
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "libsn.h"
#include "libsn_vars.h"


#define LL 2000
#define INF_LIMIT_for_EJ_DEFAULT 1e-6
#define INF_LIMIT_for_EJError_DEFAULT 1e-3

int         *_Zbins_dim, *_Mbins_dim, *_Mbins_maxdim, *_NonProcOn, *_Mbins_irr;
double      *_Zbins, **_Mbins, **Yields, *Ejecta;
void        *DataSpace;


int read_yields_file(FILE *, int);

// this macro purge end of lines in MAC- or Windows- style
#define PURGE_LINE(line) do {						\
    int L = strlen((line));						\
    int j;								\
    for( j = L-1; (j > 0) && ( (line)[j] != '\r' ); j-- );		\
    if( (line)[j] == '\r' )						\
      {									\
	(line)[j] = '\n';						\
	(line)[j+1] = '\0';						\
      }									\
  }									\
  while( 1 == 0);


int get_parameters(char *name)
{
  FILE *file;
  int n;
  char param[LL], line[LL], value[LL], *r;

  SET_ERROR(LSN_OK);
  
  /* default values */
  lsn_NMet          = 0;

  lsn_LifeTime      = PADOVANI_MATTEUCCI_1993;

  
  Yields_filename[IA] = NULL;
  Yields_filename[II] = NULL;
  Yields_filename[AGB]= NULL;

  Nset_ofYields[IA]  = 0;
  Nset_ofYields[II]  = 0;
  Nset_ofYields[AGB] = 0;

  extrap_sup[IA]    = 0;
  extrap_sup[II]    = 0;
  extrap_sup[AGB]   = 0;
  
  extrap_inf[IA]    = 1;
  extrap_inf[II]    = 1;
  extrap_inf[AGB]   = 1;

  lsn_iMU           = 100;
  lsn_Mup           = 8;
  lsn_iMi           = 0.8;
  lsn_MBM           = 16;
  lsn_MBm           = 3.0;
  lsn_BinFrac       = 0.1;
  
  lsn_IaMode        = 1;


  lsn_EjError       = INF_LIMIT_for_EJError_DEFAULT;
  
  lsn_verbose_level = VXX;
  
  if(name == 0x0)
    return -1;
  
  if((file = fopen(name, "r")) == 0x0)
    {
      SET_ERR_MSG(ERR_PARAM_FILE_OPEN);
      ADD_toERR_MSG(name);
      return -ERR_PARAM_FILE_OPEN;
    }

  n = 0;
  do
    {
      r = fgets(line, LL, file);
      n++;
      if(r != 0x0 && strchr("# \t\n", line[0]) == NULL)
	{
	  sscanf(line, "%s %s", &param[0], &value[0]);

	  if(!strcasecmp(param, "NMet"))
	    lsn_NMet = atoi(value);

	  else if(!strcasecmp(param, "Ia_Nset_ofYields"))
	    Nset_ofYields[IA] = atoi(value);

	  else if(!strcasecmp(param, "II_Nset_ofYields"))
	    Nset_ofYields[II] = atoi(value);

	  else if(!strcasecmp(param, "AGB_Nset_ofYields"))
	    Nset_ofYields[AGB] = atoi(value);

	  else if(!strcasecmp(param, "Ia_Yields_filename"))
	    {	      
	      Yields_filename[IA] = (char*)&string_space[IA_FNAME_START];
	      strcpy(Yields_filename[IA], value);
	    }

	  else if(!strcasecmp(param, "II_Yields_filename"))
	    {
	      Yields_filename[II] = (char*)&string_space[II_FNAME_START];
	      strcpy(Yields_filename[II], value);
	    }

	  else if(!strcasecmp(param, "AGB_Yields_filename"))
	    {
	      Yields_filename[AGB] = (char*)&string_space[AGB_FNAME_START];
	      strcpy(Yields_filename[AGB], value);

	    }

          else if(!strcasecmp(param, "Ia_extrap_sup"))
            extrap_sup[IA]= atoi(value);

          else if(!strcasecmp(param, "II_extrap_sup"))
            extrap_sup[II]= atoi(value);

          else if(!strcasecmp(param, "AGB_extrap_sup"))
            extrap_sup[AGB]= atoi(value);

	  else if(!strcasecmp(param, "Ia_extrap_inf"))
            extrap_inf[IA] = atoi(value);

          else if(!strcasecmp(param, "II_extrap_inf"))
            extrap_inf[II] = atoi(value);

          else if(!strcasecmp(param, "AGB_extrap_inf"))
            extrap_inf[AGB] = atoi(value);

	  else if(!strcasecmp(param, "MetalsFile"))
	    {
	      metals_file = (char*)&string_space[MET_FNAME_START];
	      strcpy(metals_file, value);
	    }

	  else if(!strcasecmp(param, "Mup"))
	    lsn_Mup = atof(value);

	  else if(!strcasecmp(param, "MB_min"))
	    lsn_MBm = atof(value);

	  else if(!strcasecmp(param, "MB_max"))
	    lsn_MBM = atof(value);

	  else if(!strcasecmp(param, "iMU"))
	    lsn_iMU = atof(value);

	  else if(!strcasecmp(param, "iMi"))
	    lsn_iMi = atof(value);

	  else if(!strcasecmp(param, "BinFrac"))
	    lsn_BinFrac = atof(value);
	  /* */

	  else if(!strcasecmp(param, "LifeT"))
	    lsn_LifeTime = atoi(value);

          else if(!strcasecmp(param, "IaMode"))
	    lsn_IaMode = atoi(value);

          else if(!strcasecmp(param, "Verbose"))
	    lsn_verbose_level = atoi(value);

	  else if(!strcasecmp(param, "Memory"))
	    Memory_amount = atoi(value);

	  else if(!strcasecmp(param, "EjError"))
	    lsn_EjError  = atof(value);
	  
	  else
	    dprintf(VXX, 0, "%s\t[w] parameter \"%s\" at line %i not known (ignored)\n", lsn_prefix, param, n);
	}
    }
  while( (r != 0x0) && !feof(file) );

  fclose(file);
  
  lsn_BinFrac_normalized = lsn_BinFrac;

  lsn_use[IA]  = (Yields_filename[IA] != 0x0);
  lsn_use[II]  = (Yields_filename[II] != 0x0);
  lsn_use[AGB] = (Yields_filename[AGB] != 0x0);
  
  if( !((lsn_use[IA] || lsn_use[II]) || lsn_use[AGB]) ||
      lsn_NMet == 0 ||
      metals_file == 0x0)
    SET_ERR_MSG(ERR_PARAM_FILE_PARAMS_MISS);

  return lsn_error;
}


/*

  This function read in the table used by metal pollution

*/

//#define DISCARD_COMMENTS_and_EMPTY_LINES  while((next_line < num_of_lines) && (strchr("%# \n", (line = lines_p[next_line++])[0]) != NULL) )

#define DISCARD_COMMENTS_and_EMPTY_LINES  if((strchr("%# \t\n", line[0]) != NULL) && (next_line < num_of_lines)) while( (next_line < num_of_lines) && (strchr("%# \n", (line = lines_p[next_line++])[0]) != NULL) )

#define GET_LINE line = lines_p[next_line++]

#define SOMETHING_WRONG ((1 << 30) + 1)



int get_Mbins(char **lines_p, int current_line, double *myMbins, int *ret_nbins, int num_of_lines)

// This function determines wheter there is a valid declaration of mass bins
//   in the text pointed by the "lines_p" pointer starting from current_line.
// If so, it returns back the number of bins and the bins themselves in myMbins
// num_of_lines states the maximum number of lines to be scanned
  
// RETURN VALUES
//   returned int value is the next_line
//   *ret_nbins:
//    integer value > 0: the number of mass bins
//    integer value < 0: some problem arose in the mass bins specifications
//    integer value = 0: no valid mass bins found in the region
//   *myMbins
//    the loaded mass bins
{

  long int  found = 0, nbins;
  int       next_line = current_line+1;
  char     *line = lines_p[current_line];

  /* if(strchr("#%", line[0]) == NULL) */
  /*   { */
  /*     *ret_nbins = 0; */
  /*     return current_line + 1; */
  /*   } */
  
  while ( !found && (next_line < num_of_lines) && (strchr("#%", line[0]) != NULL))
    {
      
      char *ppos = &line[1];
      char *pos  = strcasestr(&line[1], "Mbins");
	  
      // line must contain "Mbins"
      // and there must be anything else btw '#' and 'Mbins' than spaces or tabs
      if (pos != NULL )

	for( ; ppos < pos; ppos++)
	  if(strchr(" \t", *ppos) == NULL)
	    break;

      if( ppos == pos )
	{
	  // ok, we found a valid mass bins declaration
	  // let's see whether it follows a valid mass
	  // bins definition: i.e., the first two non-empty
	  // non-comment lines we found must be made of
	  // a single integer and a line of floats

	  // discard any comment or empty lines
	  DISCARD_COMMENTS_and_EMPTY_LINES ;
	  if(next_line == num_of_lines)
	    break;
	      
	  // check whether this line is made by a single integer
	  errno = 0;
	  nbins = strtol(line, &ppos, 10);
	  if( ( nbins > 0 ) &&
	      ( ppos > &line[0] ) &&
	      (*ppos != '\0' ) &&
	      ( errno == 0 ) )
	    {
	      // ok, the line begin with a valid integer

	      while (strchr(" \t", *ppos) != NULL)
		ppos++;
	      if(strchr("#\n", *ppos) != NULL)
		{
		  // ok, the line is made by a single valid integer plus, possibly, some comment

		  if( nbins > 1)
		    {
		      
		      // check whether the mass bins are also present
		      GET_LINE;
		      DISCARD_COMMENTS_and_EMPTY_LINES ;		      
		      if(next_line == num_of_lines)
			break;

		      // in case a problem arise in the following, the correct message will be passed through
		      found = -nbins;  		  
		      
		      int     acc = 0;
		      double *local_Mbins = alloca( nbins * sizeof(double) );
		      char    copy[LL];
		      sprintf(copy, "%s", line);
		      for(int k = 0; k < nbins; k++)
			acc += sscanf(copy, "%lg%[^n]s", &local_Mbins[k], &copy[0]);
		      if( acc / 2 != nbins )
			break;

		      found = nbins;
		      
		      if(myMbins != NULL )
			for(int k = 0; k < nbins; k++)
			  myMbins[k] = local_Mbins[k];
		      
		    }
		  else
		    {
		      found = nbins;
		      if(myMbins != NULL)
			myMbins[0] = 0;
		    }
		  GET_LINE;
		}
	    }
	  else
	    next_line--;
	}
      else
	GET_LINE;
    }


  *ret_nbins = found;
  
  return next_line;

}


int read_yields_file(FILE *file, int mode)
  
{
  int   num_of_lines = 0;
  char  s[LL], name[20];
  int   next_line;
  int   FillEl_present = -1.0;
  int   bookmark;
  int   status;
  char *lines, **lines_p, *line;
  long double *myEjecta;
  
  /*
   * read in the yields
   */

  do
    {
      if( fgets(s, LL, file) != NULL )
	num_of_lines++;
    }
  while(!feof(file));  
  
  lines   = (char* )malloc(num_of_lines * LL * sizeof(char));
  lines_p = (char**)malloc((num_of_lines + 1) * sizeof(char*));

  for(int i = 0; i < num_of_lines; i++)
    lines_p[i] = &lines[i * LL];
  lines_p[num_of_lines] = 0x0;

  rewind(file);
  {
    int   i = 0;
    do
      if( fgets(lines_p[i], LL, file) != NULL)
	{
	  /* { */
	  /*   int L = strlen((lines_p[i])); */
	  /*   int j; */
	  /*   for( j = L-1; (j >= 0) && ( lines_p[i][j] != '\r' ); j-- ); */
	  /*   if( lines_p[i][j] == '\r' ) */
	  /*     lines_p[i][j] = '\n'; */
	  /*   if( j == L-2 ) */
	  /*     lines_p[i][L-1]='\0'; */
	  /* } */
	
	  PURGE_LINE(lines_p[i]);
	  i++;
	}
    while(!feof(file));
  }

  // find a NonProcOn
  int np = 0;  
  for(np = 0; (np < num_of_lines) &&
	( (strcasestr(lines_p[np], "non proc on") == NULL) &&
	  (strcasestr(lines_p[np], "nonprocon"  ) == NULL) &&
	  (strcasestr(lines_p[np], "nonproc on" ) == NULL) ); np++);
  if( np < num_of_lines )
    *_NonProcOn = 1;
  else
    *_NonProcOn = 0;

  // discard comments and empty lines
  next_line = 0;
  line      = lines_p[0];
  DISCARD_COMMENTS_and_EMPTY_LINES ;
  
  // read in how many Z bins are in the file
  // that is supposed to be the first integer number found in the file
  
  *_Zbins_dim = atoi(line);

  if(!mode)
    {
      // does not allocate space for metal bins and read/store them
      // Zbins must have been allocated outside this routine
      // _Zbins = (double*)calloc(*_Zbins_dim, sizeof(double));
      
      if(*_Zbins_dim > 1)
	{
	  GET_LINE;
	  // if more than one metallicity bin is present, then they must be
	  // in the next non-comment non-empty lines
	  DISCARD_COMMENTS_and_EMPTY_LINES ;
	  for(int j = 0; j < *_Zbins_dim; j++)
	    sscanf(line, "%lg%[^n]s", &_Zbins[j], &line[0]);
	}
      else
	_Zbins[0] = 0;
    }

  bookmark = next_line+1;

  if(mode)
    {
      /* find if there are individual mass array for each Z bin */
      
      
      int nMbins = 0;
      status = 0;
      *_Mbins_maxdim = 0;
      
      GET_LINE;
      do
	{
	  if( strcasestr(&line[1], "Mbins") != NULL )
	    {	  
	      int local_n;
	      next_line = get_Mbins(lines_p, next_line-1, 0x0, &local_n, num_of_lines);
	      
	      if( local_n > 0)
		{
		  // account for one more array of masses
		  nMbins++;
		  
		  // store what is the largest of all the mass arrays
		  if( local_n > *_Mbins_maxdim)
		    *_Mbins_maxdim = local_n;
		}
	      else if(local_n < 0)
		{
		  status = 1;
		  break;
		}
	    }
	  
	  GET_LINE;
	}
      while( next_line < num_of_lines);
      
      if(status || (*_Mbins_maxdim == 0))
	return -1;

      free(lines_p);
      free(lines);

      // if not just one mass array is present in the file, then there must be
      // one mass array for each metallicity bin
      if(nMbins > 1 && nMbins != *_Zbins_dim)
	{
	  dprintf(VXERR, 0, "%s %d mass bins expected in the file but %d have been found \n",
		  lsn_err_prefix, *_Zbins_dim, nMbins);
	  return SOMETHING_WRONG;
	}

      
      // the function has been called with the purpose to know how many Z bins
      // there are in the file and how many M bins
      return nMbins;
    }

  
  // now, read in the yields

  next_line = bookmark;


  // if only one mass array is present, all the formal massa array pointers
  // will point to the same location; otherwise, each one will have its own
  // memory storage
  for(int j = 1; j < *_Zbins_dim; j++)
    if(*_Mbins_irr)
      _Mbins[j] = _Mbins[j-1] + *_Mbins_maxdim;
    else
      _Mbins[j] = _Mbins[0];

  // now reserve space for all the yields
  for(int j = 1; j < lsn_NMet; j++)
    Yields[j] = Yields[j-1] + *_Zbins_dim * *_Mbins_maxdim;
  
  myEjecta = (long double*)calloc(*_Zbins_dim * *_Mbins_maxdim, sizeof(long double));

  /* actually read yields. they are organized in subsequent blocks, one for each
     metal bin. each block is a table, whose rows refer to a single element and
     columns to the mass array */
  
  int j = 0;
  double *temp_store = (double*)calloc(*_Mbins_maxdim, sizeof(double));
  
  while(j < *_Zbins_dim)
    {

      // get to the next Z block, checking whether it begins with a
      // definition of mass bins
      // note: if the case (Mbins == 1 && j > 0) it won't be found a valid one,
      //       so the only effect will be to skip comment and empty lines

      int local_n;
      next_line = get_Mbins(lines_p, next_line-1, _Mbins[j], &local_n, num_of_lines);
      if( local_n > 0)
	// found a valid mass array
	_Mbins_dim[j] = local_n;

      // the mass array is equal for all Z arrays, and it was specified
      // only for the first Z bin.
      // copy it in the mass array of the current Z bin
      if(*_Mbins_irr == 0 && j > 0)
        {
          _Mbins_dim[j] = *_Mbins_maxdim;
          for(int k = 0; k < _Mbins_dim[j]; k++)
            _Mbins[j][k] = _Mbins[0][k];
        }

      // get_Mbins has changed the current line, then we must synchronize our line
      line = lines_p[next_line-1];
      // discard any other comment or empty line
      DISCARD_COMMENTS_and_EMPTY_LINES ;

      // now load the yields for this Z bin

      do
	{
	  // find the element name
	  sscanf(line, "%s %[^\n]s", name, &line[0]);        

	  int accumulate_ejecta = 0;
	  if( strcasecmp("ej", name) != 0 )
	    accumulate_ejecta = 1;
	  else
	    FillEl_present = 1;

	  memset(temp_store, 0, *_Mbins_maxdim * sizeof(double));
          for(int i = 0; i < _Mbins_dim[j]; i++)
            {
              sscanf(line, "%lg%[^\n]s", &temp_store[i], &line[0]);
	      if( accumulate_ejecta )
		// build up the total production array, in case it is not explicitly present
		// if the filling element "ej" is present, that is actually a vane effort
		// since later on the actual filling element will be used
		myEjecta[j * *_Mbins_maxdim + i] += temp_store[i];
            }

	  int k;
	  for(k = 0; k < lsn_NMet; k++)
	    if(strcasecmp(name, lsn_metal_names[k]) == 0)
	      break;
	  
	  // this is an element you want to use (as specified in metals.dat)
	  if(k < lsn_NMet)
            for(int i = 0; i < _Mbins_dim[j]; i++)
              Yields[k][*_Mbins_maxdim*j + i] = temp_store[i];

	  GET_LINE;
	}
      while((next_line < num_of_lines) && (strchr("%# \n", line[0]) == NULL));
      
      j++;      
    }

  free(temp_store);

  for(int i = 0; i < *_Zbins_dim * *_Mbins_maxdim; i++)
    Ejecta[i] = (double)myEjecta[i];
  
  free(myEjecta);
  free(lines_p);
  free(lines);
  
  return FillEl_present;
}


int read_yields(int type)

/*
  RETURN VALUES

  LSN_OK           means everything went smoothly
  < 0              letal errors
  > 0 && != LSN_OK some error occurred (typically related to Ej)
 */

{
  char  buff[300];
  FILE *file;

  if(!lsn_use[type])
    return -1;

  if(Yields_filename[type] == NULL)
    return -1;

  SET_ERROR(LSN_OK);
  
  int Nset = Nset_ofYields[type];
  
  NonProcOn[type]     = (int* )     allocate( Nset * sizeof(int) );

  Mbins_irr[type]     = (int* )     allocate( Nset * sizeof(int) );
  
  Mbins_maxdim [type] = (int* )     allocate( Nset * sizeof(int) );
  
  Zbins_dim[type]     = (int* )     allocate( Nset * sizeof(int) );

  Zbins[type]         = (double**)  allocate( Nset * sizeof(double*) );
  
  Mbins_dim[type]     = (int**)     allocate( Nset * sizeof(int*) );

  Mbins[type]         = (double***) allocate( Nset * sizeof(double**) );
  
  Y[type]             = (double***) allocate( Nset * sizeof(double**) );

  Ej[type]            = (double** ) allocate( Nset * sizeof(double*) );
  
  int i;
  for(i = strlen(Yields_filename[type])-1; i >= 0; i--)
    if( Yields_filename[type][i] == '/' )
      break;
  if(i >= 0)
    for(unsigned j = 0; j < strlen(Yields_filename[type])-i-1; j++)
      Yields_filename[type][j] = Yields_filename[type][i + j + 1];
  Yields_filename[type][strlen(Yields_filename[type]) - i -1] = '\0';
  
  for(int set = 0; set < Nset; set++)
    {
      sprintf(buff, "%s", Yields_filename[type]);
      
      if(Nset > 1)
	sprintf(&buff[strlen(Yields_filename[type])], ".%03d", set);
      
      if((file = fopen(buff, "r")) == NULL)
	{
	  SET_ERR_MSG(ERR_YIELDS_FILE_OPEN);
	  dprintf_err(VXERR, 0, "%s I can't open %s data input file: <%s>\n", lsn_err_prefix, LSN_SnNAMES[type], Yields_filename[type]);
	  return -ERR_YIELDS_FILE_OPEN;
	}

      _NonProcOn      = &NonProcOn[type][set];
      _Zbins_dim      = &Zbins_dim[type][set];      
      _Mbins_maxdim   = &Mbins_maxdim[type][set];
      _Mbins_irr      = &Mbins_irr[type][set];

      i               = read_yields_file(file, 1);
      if(i < 0)
	{
	  if(i == SOMETHING_WRONG)
	    return -ERR_YIELDS_FILE;
	}
      
      Mbins_irr[type][set] = i > 1;

      Zbins[type][set]     = (double*)allocate( Zbins_dim[type][set] * sizeof(double));
      _Zbins               = Zbins[type][set];
      
      Mbins_dim[type][set] = (int*)allocate( Zbins_dim[type][set] * sizeof(int));
      _Mbins_dim           = Mbins_dim[type][set];
      
      Mbins[type][set]     = (double**)allocate( Zbins_dim[type][set] * sizeof(double*));
      if(Mbins_irr[type][set])
	for(int aa = 0; aa < Zbins_dim[type][set]; aa++)
	  Mbins[type][set][aa] = (double*)allocate( Mbins_maxdim[type][set] * sizeof(double));

      else
	Mbins[type][set][0] = (double*)allocate( Mbins_maxdim[type][set] * sizeof(double));

      _Mbins               = Mbins[type][set];

      Y[type][set]         = (double**)allocate( lsn_NMet * sizeof(double*));
      
      for(int aa = 0; aa < lsn_NMet; aa++)
	Y[type][set][aa]   = (double*)callocate( Zbins_dim[type][set] *
						 Mbins_maxdim[type][set] *
						 sizeof(double));
      Yields               = Y[type][set];
      
      Ej[type][set]        = (double*)callocate( Zbins_dim[type][set] *
						 Mbins_maxdim[type][set] * sizeof(double));
      Ejecta               = Ej[type][set];


      rewind(file);
      i =  read_yields_file(file, 0);
      
      fclose(file);

      int FillEl_present = (i > 0);
            
      if(lsn_FillEl >= 0)
	{
	  if(FillEl_present)
	    
	    // if FillEl is already present, use the data read from the table
	    memcpy((void*)Ej[type][set], (void*)Y[type][set][lsn_FillEl],
		   (size_t)(Zbins_dim[type][set] * Mbins_maxdim[type][set] * sizeof(double)));

	  // now check that the game is a zero-sum game

	  int minor_flaws = 0;
	  
	  for(i = 0; i < Zbins_dim[type][set]; i++)
	    
	    for(int j = 0; j < Mbins_dim[type][set][i]; j++)
	      {
		
		// set the FillEl in the table to zero
		Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j] = 0;

		// now sum-up all the produced metals
		// -> this amount to the fact that the listed metals account for all the mass production
		long double accumulate = 0;
		for(int k = 0; k < lsn_NMet; k++)
		  if(k != lsn_FillEl)
		    accumulate += Y[type][set][k][i*Mbins_maxdim[type][set] + j];
		Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j] = accumulate;

		// now compare the sum of all the products with the content of Ej
		// - in case it was NOT present in the file, the two numbers must be equal. if the relative
		//   difference is tiny, let's put the FillEl to zero.
		long double compare;
		compare = fabs(Ej[type][set][i*Mbins_maxdim[type][set] + j] -
			       Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j]);
		
		if(Ej[type][set][i*Mbins_maxdim[type][set] + j] > 0)
		   compare /= Ej[type][set][i*Mbins_maxdim[type][set] + j];
		
		if( compare < INF_LIMIT_for_EJ_DEFAULT )
		  Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j] = 0;
	    
		else
		  // the difference is NOT negligible (so, the FillEl has been read from the file and it accounts for
		  // some elements that are not explicitly listed in the table.
		  // Then account for the FillEl.
		  {
		    Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j] =
		      Ej[type][set][i*Mbins_maxdim[type][set] + j] - Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j];
		
		    if(Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j] < 0)
		      // something is still wrong, probably only FP round-offs
		      {
			if( fabs(Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j]) > lsn_EjError )
			  {
			    char buffer[300];
			    sprintf(buffer, "\n%s in yields set %d for %s : the sum of yields is larger (%g vs %g) "
				    "than supposed total ejection for Zbin %d Mass bin %d\n", lsn_err_prefix,
				    set, LSN_SnNAMES[type],
				    Ej[type][set][i*Mbins_maxdim[type][set] + j] - Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j],
				    Ej[type][set][i*Mbins_maxdim[type][set] + j],
				    i, j);
			    SET_ERR_MSG(ERR_YIELDS_Ej_SEVERE);
			    ADD_toERR_MSG(buffer);
			    return lsn_error;
			  }
			else
			  minor_flaws++;
			Y[type][set][lsn_FillEl][i*Mbins_maxdim[type][set] + j] = 0;
		      }
		  }
	      }

	  if(minor_flaws)
	    {
	      
	      char buffer[300];
	      sprintf(buffer, "\n%s in  %d times the sum of yields for a given mass was "
		      "slightly larger than the Fill El. Due to small difference, I think it's due to limited precision\n",
		      lsn_warn_prefix, minor_flaws);

	      int lsn_error_save = lsn_error;
	      ADD_ERR_MSG(ERR_YIELDS_Ej);
	      ADD_toERR_MSG(buffer);

	      if(lsn_error_save < SEv_SEVERE_ERROR)
	      	lsn_error = lsn_error_save;
	    }		    
	  
	}


    }

  
  return lsn_error;
}



int read_elements(char *name)
{
  FILE *file;
  char line[100];

  int lsn_NMet_save   = lsn_NMet;
  int allocated_space = lsn_NMet * 2;
  
  if((name == NULL) || (file = fopen(name, "r")) == NULL)
    {
      //dprintf_err(VERR, 0, "%s unable to open metals' file %s\n", lsn_err_prefix, name);
      return -1;
    }

  lsn_metal_names = (char**)allocate(allocated_space * sizeof(char*));  
  lsn_Hyd = lsn_Hel = lsn_FillEl = -1;

  int n = 0;
  int r = 0;
  do
    {
      r = fscanf(file, "%[^\n]\n", &line[0]);
      PURGE_LINE(line);

      if( (r != EOF) && (strchr("# \t\n", line[0]) == NULL) )
	{
	  if(sscanf(line, "%[^\n]s", &line[0]) == 1)
	    {
	      if(n > lsn_NMet - 1)
		{
		  lsn_NMet++;
		  if(lsn_NMet > allocated_space)
		    {
		      allocated_space = lsn_NMet * 2;
		      lsn_metal_names = (char**)reallocate(lsn_metal_names,
							   (allocated_space * 2) * sizeof(char*));
		    }
		}

	      lsn_metal_names[n] = (char*)allocate(strlen(line)+1);
	      strcpy(lsn_metal_names[n], line);
	      
	      if(strcasecmp(line, "H")==0)
		lsn_Hyd = n;
	      else if(strcasecmp(line, "He")==0)
		lsn_Hel = n;
	      else if(strcasecmp(line, "Ej")==0)
		lsn_FillEl = n;
	      n++;
	    }
	}
    }
  while(r != EOF);

  // ensures that all the names' initials are upper case
  for(int i = 0; i < lsn_NMet; i++)
    lsn_metal_names[i][0] = toupper(lsn_metal_names[i][0]);
  
  if(lsn_NMet > lsn_NMet_save)
    dprintf_err(VMSG, 0, "%smore metals (%d) found in file %s than specified in param file (%d): data have been accordingly enlarged\n",
	    lsn_prefix, lsn_NMet, name, lsn_NMet_save); 
  
  if(lsn_FillEl == -1)
    dprintf_err(VMSG, 0, "%sFillEl is not explicitly accounted !\n", lsn_prefix);

  fclose(file);
  return n;
}



void lsn_dump_yields(int type, int set)
{
  type >>= 1;
  
  if(type < IA || type > AGB)
    return;

  if( Y[type] == NULL ||
      Y[type][set] == NULL ||
      Zbins_dim[type] == NULL ||
      Mbins_dim[type] == NULL ||
      Mbins_dim[type][set] == NULL )
    return;
  
  FILE *file;
  char buff[200];
  sprintf(buff, "%s_dump", Yields_filename[type]);
  int pos = strlen(buff);

  if(Nset_ofYields[type] > 1)
    sprintf(&buff[pos], ".%03d", set);
  
  file = fopen(buff, "w");
  
  fprintf(file, "# Zbins\n%d", Zbins_dim[type][set]);
  if(Zbins_dim[type][set] > 1)
    {
      fprintf(file,"\n");
      for(int Z = 0; Z < Zbins_dim[type][set]; Z++)	
	fprintf(file, "%12.7e ", Zbins[type][set][Z]);
    }
  fprintf(file, "\n#\n");
  
  for(int Z = 0; Z < Zbins_dim[type][set]; Z++)
    {
      fprintf(file, "# Z = %6.4e\n#\n", Zbins[type][set][Z]);
      fprintf(file, "# Mbins\n%d", Mbins_dim[type][set][Z]);
      if(Mbins_dim[type][set][Z] > 1)
	{
	  fprintf(file, "\n");
	  for(int M = 0; M < Mbins_dim[type][set][Z]; M++)
	    fprintf(file, "%12.7e ", Mbins[type][set][Z][M]);
	}
      fprintf(file, "\n#\n");
      for(int E = 0; E < lsn_NMet; E++)
	{
	  fprintf(file, "%-3s\t", lsn_metal_names[E]);
	  for(int M = 0; M < Mbins_dim[type][set][Z]; M++)
	    fprintf(file, "%12.7e\t", Y[type][set][E][Mbins_maxdim[type][set]*Z + M]);
	  fprintf(file, "\n");
	}
      fprintf(file, "#\n");
    }
  
  
  fclose(file);
  
  return;
}

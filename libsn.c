#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "libsn.h"
#include "libsn_vars.h"


/* /------------------------------------------------\  
 * | HIGH LEVEL API                                 |
 * |                                                |
 * | > initialize_libsn                             |
 * | > get_metals                                   |
 * | > get_sn_nuber                                 |
 * | > make_time_interval_bins_by_events            |
 * \------------------------------------------------/ */

int INLINE_FUNC getindex(double*, int, int, double*, int*);

double nRSnIa_DTD(double, void *);

double DTD_SD(double, void*);                                            /* Iamode = 1 */
double DTD_SD_oldstyle(double, void*);                                   /* Iamode = 5 */
double DTD_Mannucci_DellaValle_Panagia_2006(double, void*);              /* Iamode = 2 */
double DTD_gaussian(double, void*);                                      /* Iamode = 3 */
double DTD_powerlaw(double, void*);                                      /* Iamode = 4 */

double nRSnIa_SD_Greggio_and_Renzini(double, void*);


double dm_dt_PADOVANI_MATTEUCCI_1993(double, double);
double dm_dt_MAEDER_MEYNET_1989(double, double);
double dm_dt_generic(double, double);
double (*_dm_dt)(double, double);

double lifetime_PADOVANI_MATTEUCCI_1993(double);
double lifetime_MAEDER_MEYNET_1989(double);
double lifetime_generic(double);
double (*_lifetime)(double);


double dying_mass_PADOVANI_MATTEUCCI_1993(double);
double dying_mass_MAEDER_MEYNET_1989(double);
double dying_mass_generic(double);
double (*_dying_mass)(double);


double INLINE_FUNC sec_dist(double, double);
double initialize_Ia(int, double, double);


double (* lsn_nRSnIa)(double, void*);
double (* lsn_DTDfunction)(double, void*);


// -----------------------------------------------------------------
// -----------------------------------------------------------------
// define APIs to access tables meta-data and data
//



int NUM_of_ZBINS(int T, int SET)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( SET < Nset_ofYields[(T)] ) return Zbins_dim[T][SET]; else { SET_ERROR(ERR_BOUNDS); return 0; };
}

double *ZBINS(int T, int SET)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( SET < Nset_ofYields[T] ) return Zbins[T][SET]; else { SET_ERROR(ERR_BOUNDS); return NULL; };
}

int *MBINS_dim(int T, int SET)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( SET < Nset_ofYields[T] ) return Mbins_dim[T][SET]; else { SET_ERROR(ERR_BOUNDS); return NULL; };
}

int NUM_of_MBINS(int T, int SET, int ZBIN)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( (SET < Nset_ofYields[T]) && (ZBIN < Zbins_dim[T][SET])) return Mbins_dim[T][SET][ZBIN]; else { SET_ERROR(ERR_BOUNDS); return 0; };
}

double **MBINS(int T, int SET)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( SET < Nset_ofYields[T] ) return Mbins[T][SET]; else { SET_ERROR(ERR_BOUNDS); return NULL; };
}

double *MBINS_ZBIN(int T, int SET, int ZBIN)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( (SET < Nset_ofYields[T]) && (ZBIN < Zbins_dim[T][SET])) return Mbins[T][SET][ZBIN]; else { SET_ERROR(ERR_BOUNDS); return NULL; };
}

double *YIELDS_all(int T, int SET, int EL)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( (SET < Nset_ofYields[T]) && (EL < lsn_NMet)) return Y[T][SET][EL];
  else { SET_ERROR(ERR_BOUNDS); return NULL; };
}

double *YIELDS(int T, int SET, int ZBIN, int EL)
{
  SET_ERROR(LSN_OK);
  T = T>>1;
  if( (SET < Nset_ofYields[T]) && (ZBIN < Zbins_dim[T][SET]) && (EL < lsn_NMet))
    return &Y[T][SET][EL][ZBIN * Mbins_maxdim[T][SET] * lsn_NMet];
  else { SET_ERROR(ERR_BOUNDS); return NULL; };
}



// -----------------------------------------------------------------
// -----------------------------------------------------------------



int lsn_initialize_libsn(char *IMFfilename, char* PARAMfilename)
{
  SET_ERROR(LSN_OK);
  
  string_space = (char*)aligned_alloc(ALIGN, TOT_STR_SPACE);

  int res;


  // [ 1 ] read the parameter file and manage errors
  
  if(PARAMfilename == 0x0)
    res = get_parameters("snparam.dat");
  else
    res = get_parameters(PARAMfilename);

  if( (res <= LETAL_ERROR) && (lsn_error != LSN_OK) )
    {
      if(lsn_error > 0)
	lsn_error = -lsn_error;
      return lsn_error;
    }

  if ( Memory_amount > 0 )  
    if( (manage_memory(Memory_amount) == NULL ) )
      {
	SET_ERR_MSG(ERR_MEMORY); // an error in memory allocation is a letal error
	return -ERR_MEMORY;
      }

  // set-up lifetime
  if(lsn_LifeTime == PADOVANI_MATTEUCCI_1993)
    {
      _lifetime   = lifetime_PADOVANI_MATTEUCCI_1993;
      _dm_dt      = dm_dt_PADOVANI_MATTEUCCI_1993;
      _dying_mass = dying_mass_PADOVANI_MATTEUCCI_1993;
    }
  else if(lsn_LifeTime == MAEDER_MEYNET_1989)
    {
      _lifetime   = lifetime_MAEDER_MEYNET_1989;
      _dm_dt      = dm_dt_MAEDER_MEYNET_1989;
      _dying_mass = dying_mass_MAEDER_MEYNET_1989;
    }
  else
    {
      _lifetime   = lifetime_generic;
      _dm_dt      = dm_dt_generic;
      _dying_mass = dying_mass_generic;
    }
  

  // [ 2 ] get what elements to track
  
  if( read_elements("metals.dat") <= 0)
    {
      ADD_ERR_MSG(ERR_METALS_FILE);
      STELLAR_EVOLUTION_IMPOSSIBLE |= NO_METALS_DEFINED;
    }

  
  // [ 3 ] now read yields files
  
  lsn_initial_metallicity = (double*)calloc(lsn_NMet, sizeof(double));
  lsn_initial_abundance   = 0;
 
  if( !STELLAR_EVOLUTION_IMPOSSIBLE )
    {
      int not_ok = 0, mild_errors = 0;

      if(lsn_use[II])
	// II yields file(s) has(ve) been defined in paramfile
	{
	  int rep = read_yields(II); 
	  if( (rep < SEv_SEVERE_ERROR) && (rep != LSN_OK) )
	    {
	      not_ok++;
	      lsn_use[II] = 0;
	    }
	  else
	    {
	      if( rep != LSN_OK )
		mild_errors++;
	      dprintf(VMSG, 0, "%s using SnII ", lsn_prefix);
	      for(int i = 0; i < Nset_ofYields[II]; i++)
		{
		  if(NonProcOn[II][i]) {dprintf(VMSG, 0, "[ NonProcOn  ] ");}
		  else { dprintf(VMSG, 0, "[ NonProcOff ] ");}
		}
	      dprintf(VMSG, 0, "\n");
	    }
	}
    
      if(lsn_use[IA])
	// Ia yields file(s) has(ve) been defined in paramfile
	{
	  int rep = read_yields(IA); 
	  if( (rep < SEv_SEVERE_ERROR) && (rep != LSN_OK) )
	    {
	      not_ok++;
	      lsn_use[IA] = 0;
	    }
	  else
	    {
	      if( rep != LSN_OK )
		mild_errors++;
	      dprintf(VMSG, 0, "%s using SnIa ", lsn_prefix);
	      for(int i = 0; i < Nset_ofYields[IA]; i++)
		if(NonProcOn[IA][i]) { dprintf(VMSG, 0, "[ NonProcOn  ] "); }
		else { dprintf(VMSG, 0, "[ NonProcOff ] "); }
	      dprintf(VMSG, 0, "\n");
	    }

	}
	  
    
      if(lsn_use[AGB])
	// AGB yields file(s) has(ve) been defined in paramfile
	{
	  int rep = read_yields(AGB); 
	  if( (rep < SEv_SEVERE_ERROR) && (rep != LSN_OK) )
	    {
	      not_ok++;
	      lsn_use[AGB] = 0;
	    }
	  else
	    {
	      if( rep != LSN_OK )
		mild_errors++;	      
	      dprintf(VMSG, 0, "%s using AGB  ", lsn_prefix);
	      for(int i = 0; i < Nset_ofYields[AGB]; i++)
		if(NonProcOn[AGB][i]) { dprintf(VMSG, 0, "[ NonProcOn  ] ");}
		else { dprintf(VMSG, 0, "[ NonProcOff ] ");}
	      dprintf(VMSG, 0, "\n");
	    }
	}


      if(not_ok || mild_errors)
	{
	  if(mild_errors)
	    dprintf_err(VERR, 0, "%s\n", lsn_error_message);
	  
	  // some error occurred in one or more yields file
	  ADD_ERR_MSG(ERR_YIELDS_FILE);
	  if(not_ok == 3)
	    {
	      // errors occurred in all yields files
	      STELLAR_EVOLUTION_IMPOSSIBLE |= NO_YIELDS_DEFINED;

	      dprintf_err(VXERR, 0, "%s\n", lsn_error_message);
	    }
	}
    }
  else
    {
      SET_ERR_MSG(ERR_SEv_IMPOSSIBLE_noMETALS_noYIELDS_loaded);
      STELLAR_EVOLUTION_IMPOSSIBLE |= NO_YIELDS_DEFINED;
    }
  
  if( (lsn_w = gsl_integration_workspace_alloc(lsn_gsl_intspace_dim)) == NULL )
    {
      ADD_ERR_MSG(ERR_GSL_ALLOC);     // gsl allocation error is a letal error
      return lsn_error;
    }

  if( limf_initialize_imfs(IMFfilename, 0) < 0)
    {
      sprintf(lsn_error_message, " %s %s\n", lsn_prefix, limf_error_message); // imf error is a letal error
      return ERR_IMF;
    }
  
  //IMFp = &IMFs[0];
  lsn_mean_lifetime = lsn_lifetime(lsn_Mup);
    
  dprintf(VMSG, 0, "\n%s stellar mean lifetime:   mean (%5.4g Msun) = %g Gyrs\n\n", lsn_prefix, 
	  lsn_Mup, lsn_mean_lifetime);
  fflush(stdout);

  lsn_iinf_lifetime = lsn_lifetime(lsn_iMU);
  lsn_isup_lifetime = lsn_lifetime(lsn_iMi);

  for(int i = 0; i < limf_NIMFs; i++)
    {
      double A = limf_IntegrateIMF_byMass(IMFs[i].Mm, IMFs[i].MU, &IMFs[i], INC_BH);
      if( fabs(A-1)/A > 0.0001)
        {
          dprintf(VMSG, 0, "%s renormalizing.. (%g) ", lsn_prefix, fabs(A-1)/A);
          A = limf_renormalize_IMF(&IMFs[i], IMFs[i].Mm, IMFs[i].MU, INC_BH);
          limf_write_IMF_info(&IMFs[i], stdout);
        }
      
      IMFs[i].inf_lifetime = lsn_lifetime(IMFs[i].MU);
      if(IMFs[i].inf_lifetime > lsn_iinf_lifetime)
        {
          dprintf(VMSG, 0, "%s inf integration time raised from %g to %g, accordingly to limits in IMF %d\n",
		  lsn_prefix, lsn_iinf_lifetime, IMFs[i].inf_lifetime, i);
          lsn_iinf_lifetime = IMFs[i].inf_lifetime;
        }
      
      IMFs[i].sup_lifetime = lsn_lifetime(IMFs[i].Mm);
      if(IMFs[i].sup_lifetime < lsn_isup_lifetime)
        {
          dprintf(VMSG, 0, "%s sup integration time reduced from %g to %g, accordingly to limits in IMF %d\n",
		  lsn_prefix, lsn_isup_lifetime, IMFs[i].sup_lifetime, i);
          lsn_isup_lifetime = IMFs[i].sup_lifetime;
        }

      dprintf(VMSG, 0, "\n"
	      "%s   IMF %3d      inf (%5.4g Msun) = %g Gyrs\n"
	      "%s                sup (%5.4g Msun) = %g Gyrs\n",
	      lsn_prefix, 
	      i, IMFs[i].MU, IMFs[i].inf_lifetime,
	      lsn_prefix, 
	      IMFs[i].Mm, IMFs[i].sup_lifetime);
      
      fflush(stdout);
    }

  
  switch(lsn_IaMode)
    {
    case 0:
      lsn_nRSnIa = nRSnIa_SD_Greggio_and_Renzini;
      break;
      
    case 1:
      lsn_nRSnIa = nRSnIa_DTD;
      lsn_DTDfunction = DTD_SD;
      break;
      
    case 2:
      lsn_nRSnIa = nRSnIa_DTD;
      lsn_DTDfunction = DTD_Mannucci_DellaValle_Panagia_2006;
      break;
      
    case 3:
      lsn_nRSnIa = nRSnIa_DTD;
      lsn_DTDfunction = DTD_gaussian;
      break;
      
    case 4:
      lsn_nRSnIa = nRSnIa_DTD;
      lsn_DTDfunction = DTD_powerlaw;
      break;
      
    case 5:
      lsn_nRSnIa = nRSnIa_DTD;
      lsn_DTDfunction = DTD_SD_oldstyle;
      break;
      
    }
  
  if(lsn_IaMode > 0)
    for(int i = 0; i < limf_NIMFs; i++)
      initialize_Ia(i, lsn_mean_lifetime, lsn_isup_lifetime);
    
  // an error here is either < LETAL_ERROR - and that means it is impossible to continue -
  // or it is > LETAL_ERROR and < SEv_SEVERE_ERROR, which means that stellar evolution
  // calculations are impossible or partially impossible (either metals or yields are not
  // define or there has been some problem with one or more yields tables)
  
  if( lsn_error < LETAL_ERROR )
    return -lsn_error;
  else
    return lsn_error;
  
}




void lsn_finalize_libsn(void)
{
  SET_ERROR(LSN_OK);
  gsl_integration_workspace_free(lsn_w);

  release(Memory);

  free(string_space);
  
  return; 
}




void lsn_get_libsn_options(char **ostring)
{
  if ( ostring != NULL )
    {
      *ostring = (char*) malloc ( strlen(options_string) + 1 );
      sprintf ( *ostring, "%s", options_string );
    }
  
  return;
}


void lsn_set_initial_metallicity(int N, ...)
{
  va_list arglist;
  
  SET_ERROR(LSN_OK);
  
  assert (N == lsn_NMet);
  
  va_start(arglist, N);

  for(int i = 0; i < N; i++)
    lsn_initial_metallicity[i] = va_arg(arglist, double);
  
  va_end(arglist);

  long double sum = 0;
  for(int i = 0; i < lsn_NMet; i++)
    sum += lsn_initial_metallicity[i];

  if(lsn_Hyd >= 0 )
    sum -= lsn_initial_metallicity[lsn_Hyd];

  if(lsn_Hel >= 0 )
    sum -= lsn_initial_metallicity[lsn_Hel];

  lsn_initial_abundance = (double)(sum / (1.0 - sum));
  return;
}


int lsn_get_metal_position( char *metal_name )
{
  int pos;
  for ( pos = 0; pos < lsn_NMet; pos++ )
    if  ( strcmp( metal_name, lsn_metal_names[pos]) == 0 )
      break;

  return ( pos < lsn_NMet ? pos : -1 );
}


int lsn_get_metals(int type,
		   int mode,
                   int Yset,
		   double Zstar,
                   double inf,
		   double sup,
		   int IMFi,
                   int *check_elements,
                   double *mymetals)

     /*
      * type     =  1, 2, 4 for Ia, II and intermediate mass stars respectively, bitwise combined
      * mode     =  1, 0 for mass or time integration respectively
      * Yset     =  the yields set
      * Zstar    =  the metallicity of the stellar population
      * inf, sup =  the integration range (should be interpreted as masses or
                      time accordingly to the valu of mode)
      * check_elements =  an array of integer whose values set an element on or
                          off; calculations will be performed for the elements
                          set to 1
      * mymetals = the array in which the results are stored
      */

{
  int        gsl_status;
  int        calc;
  double     res, err, num;
  double     myinf, mysup;
#ifdef DEBUG
  unsigned int int_errors = 0;
#endif
  //lsn_SDtype SD.

  if(STELLAR_EVOLUTION_IMPOSSIBLE)
    {
      SET_ERR_MSG(ERR_SEv_IMPOSSIBLE + STELLAR_EVOLUTION_IMPOSSIBLE);
      return ERR_SEv_IMPOSSIBLE + STELLAR_EVOLUTION_IMPOSSIBLE;
    }

  sup = (sup > model->isup_lifetime ? model->isup_lifetime : sup);
  inf = (inf < model->iinf_lifetime ? model->iinf_lifetime : inf);

  if ( sup < inf )
    {
      for ( int i = 0; i < lsn_NMet; i++ )
	mymetals[i] = 0;
      SET_ERR_MSG(ERR_LIMITS_INFSUP);
      return ERR_LIMITS_INFSUP;
    }

  SET_ERROR(LSN_OK);
  
  SD.Yset = Yset;
  SD.Zstar = Zstar;

  calc = 1;


  //double *temp = (double*)calloc(lsn_NMet, sizeof(double));
  double temp[lsn_NMet];
  
  //lsn_F.params = (void*)&SD;
  
  if( (type & lsn_SnIa) &&
      lsn_use[IA] &&
      SD.Yset < Nset_ofYields[IA])
    /* make calculations for Ia = 1*/
    {
      SD.Zdim = Zbins_dim[IA][SD.Yset];
      for(SD.Zbin = Zbins_dim[IA][SD.Yset]-1; SD.Zbin > 0 && Zstar < Zbins[IA][SD.Yset][SD.Zbin]; SD.Zbin--)
	;
      SD.ZArray  = Zbins       [IA][SD.Yset];
      SD.MArray  = Mbins       [IA][SD.Yset][SD.Zbin];
      SD.Mdim    = Mbins_dim   [IA][SD.Yset][SD.Zbin];
      SD.Mmaxdim = Mbins_maxdim[IA][SD.Yset];
      if( (SD.Mirr = Mbins_irr[IA][SD.Yset]) )
	{
	  if( (SD.Zbin < Zbins_dim[IA][SD.Yset]-1) && (SD.Zbin >= 0))
	    {
	      SD.MArray1 = Mbins    [IA][SD.Yset][SD.Zbin+1];
	      SD.Mdim1   = Mbins_dim[IA][SD.Yset][SD.Zbin+1];
	    }
	}

      lsn_F.function = lsn_nRSnIa;
      lsn_F.params   = &IMFs[IMFi];

      if(mode)
	{
	  //mass integration
	  myinf = lsn_lifetime(sup);
	  mysup = lsn_lifetime(inf);
	}
      else
	{
	  //time integration
	  myinf = inf;
	  mysup = sup;
	}

      SD.extrap     = extrap_sup[IA];
      SD.extrap_inf = extrap_inf[IA];

      num = 0;
      if( (gsl_status = LSN_INT(&lsn_F, myinf, mysup,
				1e-6, 1e-4,
				lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
				&num, &err)) )
	{
	  //free(temp);
	  sprintf(lsn_error_message, "%s :: gsl integration error %d in Sn (nimf)"
		  " [%.6g - %.6g] : %.6g %.6g (line %d in file %s, at func %s)\n",
		  lsn_error_msgs[ERR_GSL_INTEGRATION],
		  gsl_status, myinf, mysup, num, err, __LINE__, __FILE__, __func__);
	  SET_ERROR( ERR_GSL_INTEGRATION );
	  return ERR_GSL_INTEGRATION;
	}

      if ( num > 0 )
	for ( int i = 0; i < lsn_NMet; i++ )
	  if ( check_elements[0] == -1 || check_elements[i] > 0 )
	    mymetals[i] += num * Y[IA][SD.Yset][i][0];
#ifdef DEBUG
      if ( num < 0 )
	{
	  int_errors++;
	  SET_ERROR( ERR_INTEGRATION );
	}
#endif
      
    }


  if ( (type & lsn_SnII) &&
       lsn_use[II] &&
       SD.Yset < Nset_ofYields[II])
    /* make calculations for II = 2 */
    {
      myinf = inf;
      mysup = sup;
      calc = 1;

      if ( mode )
	//mass integration
	{
	  lsn_F.function = &lsn_zmRSnII;
	  if ( mysup < lsn_Mup )
	    calc = 0;
	  else if ( myinf < lsn_Mup )
	    myinf = lsn_Mup;
	}
      else
	// time integration
	{
	  lsn_F.function = &lsn_ztRSnII;
	  if ( myinf > lsn_mean_lifetime )
	    calc = 0;
	  else if ( mysup > lsn_mean_lifetime )
	    mysup = lsn_mean_lifetime;	  
	}

      if ( calc )
	{
	  SD.Zdim = Zbins_dim[II][SD.Yset];
	  for ( SD.Zbin = Zbins_dim[II][SD.Yset]-1; SD.Zbin > 0 && Zstar < Zbins[II][SD.Yset][SD.Zbin]; SD.Zbin-- )
	    ;

	  SD.Mmaxdim = Mbins_maxdim[II][SD.Yset];
	  SD.Mdim    = Mbins_dim   [II][SD.Yset][SD.Zbin];
	  SD.ZArray  = Zbins       [II][SD.Yset];
	  SD.MArray  = Mbins       [II][SD.Yset][SD.Zbin];
	  if ( (SD.Mirr = Mbins_irr[II][SD.Yset]) )
	    {
	      if ( (SD.Zbin < Zbins_dim[II][SD.Yset]-1) && (SD.Zbin >= 0) )
		{
		  SD.MArray1 = Mbins    [II][SD.Yset][SD.Zbin+1];
		  SD.Mdim1   = Mbins_dim[II][SD.Yset][SD.Zbin+1];
		}
	    }
	  SD.extrap     = extrap_sup[II];
	  SD.extrap_inf = extrap_inf[II];

	  imf_plus_param_type params = { IMFi, NULL };
	  lsn_F.params = &params;

	  for (int i = 0; i < lsn_NMet; i++ )
	    if ( check_elements[0] == -1 || check_elements[i] )
	      {            
		SD.Y = Y[II][SD.Yset][i];
                
		if( (gsl_status = LSN_INT(&lsn_F, myinf, mysup,
					  1e-6, 1e-4,
					  lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
					  &temp[i], &err)) )
		  {
		    //free ( temp );
		    sprintf(lsn_error_message, "%s :: gsl integration error %i in Sn (nimf)"
			    " [%.6g - %.6g] : %.6g %.6g (line %d in file %s, at func %s)\n",
			    lsn_error_msgs[ERR_GSL_INTEGRATION],
			    gsl_status, myinf, mysup, num, err, __LINE__, __FILE__, __func__);
		    SET_ERROR(ERR_GSL_INTEGRATION);
			 
		    return ERR_GSL_INTEGRATION;
		  }

	      }

#ifdef DEBUG
	  {
	    for ( int i = 0; i < lsn_NMet; i++ )
	      if ( temp[i] < 0 )
		{
		  int_errors++;
		  SET_ERROR( ERR_INTEGRATION );
		}
	  }
#endif
	  
	  if ( NonProcOn[II][SD.Yset] && lsn_initial_abundance > 0 )
	    if ( lsn_Hyd || lsn_Hel >= 0 )
	      {
		double fact = 0;
		    
		if ( lsn_Hyd > 0 )
		  fact += temp[lsn_Hyd];

		if ( lsn_Hel > 0 )
		  fact += temp[lsn_Hel];

		fact *= (1 - lsn_initial_abundance) / (1 + lsn_initial_abundance);

		for ( int i = 0; i < lsn_NMet; i++ )
		  if( (i != lsn_Hyd) && (i != lsn_Hel) )
		    temp[i] += lsn_initial_metallicity[i] * fact;
		    
		temp[lsn_Hyd] /= (1 + lsn_initial_abundance);
		temp[lsn_Hel] /= (1 + lsn_initial_abundance);
	      }
	      
	  for ( int i = 0; i < lsn_NMet; i++ )
	    mymetals[i] += temp[i];

	}
    }


  if( (type & lsn_AGB) &&
      lsn_use[AGB] &&
      SD.Yset < Nset_ofYields[AGB])
    /* make calculations for AGB = 4 */
    {
      myinf = inf;
      mysup = sup;
      calc = 1;

      if(mode)
	// mass integration
	{
	  lsn_F.function = &lsn_zmRSnII;
	  if(myinf >= lsn_Mup)
	    calc = 0;
	  else if(mysup > lsn_Mup)
	    mysup = lsn_Mup;	  
	}
      else
	// time integration
	{
	  lsn_F.function = &lsn_ztRSnII;
	  if(mysup <= lsn_mean_lifetime)
	    calc = 0;
	  else if(myinf < lsn_mean_lifetime)
	    myinf = lsn_mean_lifetime;	  
	}
      
      if(calc)
	{            
      
	  SD.Zdim = Zbins_dim[AGB][SD.Yset];
	  for ( SD.Zbin = Zbins_dim[AGB][SD.Yset] - 1; SD.Zbin > 0 && Zstar < Zbins[AGB][SD.Yset][SD.Zbin]; SD.Zbin-- )
	    ;

	  SD.Mdim    = Mbins_dim   [AGB][SD.Yset][SD.Zbin];
	  SD.Mmaxdim = Mbins_maxdim[AGB][SD.Yset];
	  SD.ZArray  = Zbins       [AGB][SD.Yset];
	  SD.MArray  = Mbins       [AGB][SD.Yset][SD.Zbin];
	  if( (SD.Mirr = Mbins_irr[AGB][SD.Yset]) )
	    {
	      if( (SD.Zbin < Zbins_dim[AGB][SD.Yset]-1) && (SD.Zbin >= 0))
		{
		  SD.MArray1 = Mbins    [AGB][SD.Yset][SD.Zbin+1];
		  SD.Mdim1   = Mbins_dim[AGB][SD.Yset][SD.Zbin+1];
		}
	    }
	  SD.extrap     = extrap_sup[AGB];
	  SD.extrap_inf = extrap_inf[AGB];

	  imf_plus_param_type params = { IMFi, NULL };
	  lsn_F.params = &params;

	  for ( int i = 0; i < lsn_NMet; i++ )
	    if(check_elements[0] == -1 || check_elements[i])
	      {            
		SD.Y = Y[AGB][SD.Yset][i];

		if((gsl_status = LSN_INT(&lsn_F, myinf, mysup,
					 1e-6, 1e-4,
					 lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
					 &res, &err)))
		  {
		    //free(temp);
		    sprintf(lsn_error_message, "%s :: gsl integration error %i in Sn (nimf)"
			    " [%.6g - %.6g] : %.6g %.6g (line %d in file %s, at func %s)\n",
			    lsn_error_msgs[ERR_GSL_INTEGRATION],
			    gsl_status, myinf, mysup, num, err, __LINE__, __FILE__, __func__);
		    SET_ERROR(ERR_GSL_INTEGRATION);
		    return ERR_GSL_INTEGRATION;
		  }
		temp[i] = (1 - lsn_BinFrac_normalized) * res;
	      }

#ifdef DEBUG
	  {
	    for ( int i = 0; i < lsn_NMet; i++ )
	      if ( temp[i] < 0 )
		{
		  int_errors++;
		  SET_ERROR( ERR_INTEGRATION );
		}
	  }
#endif	  

	  if(NonProcOn[AGB][SD.Yset] && lsn_initial_abundance > 0)
	    if(lsn_Hyd || lsn_Hel >= 0)
	      {
		double fact = 0;
		    
		if(lsn_Hyd > 0)
		  fact += temp[lsn_Hyd];

		if(lsn_Hel > 0)
		  fact += temp[lsn_Hel];

		fact *= (1 - lsn_initial_abundance) / (1 + lsn_initial_abundance);

		for ( int i = 0; i < lsn_NMet; i++ )
		  if( (i != lsn_Hyd) && (i != lsn_Hel) )
		    temp[i] += lsn_initial_metallicity[i] * fact;
		    
		temp[lsn_Hyd] /= (1 + lsn_initial_abundance);
		temp[lsn_Hel] /= (1 + lsn_initial_abundance);
	      }
	      
	  for ( int i = 0; i < lsn_NMet; i++ )
	    mymetals[i] += temp[i];

	}
    }
  
  //free(temp);

#ifdef DEBUG
  if ( int_errors > 0 )
    return -int_errors;
  else
#endif
  return LSN_OK;
}




double lsn_get_sn_number(int type, int IMFi, double tinf, double tsup)
/*
  note: That should substitute num_time_integration()
    the main changes are:
  (1) type is no mopre 1, 2 and 3 for SnIa, SnII and AGB respectively but
      1, 2 and 4, bitwise combined (use defined symbols lsn_SnIa, lsn_SnII,
      lsn_AGB respectively).
  (2) you can compute at the same time the number of more than one type
      of events in whatever time interval (type = lsn_ALL returns the result
      for SnIa, SnII and AGB stars summed up).
 */
{
  int gsl_status;
  double mytinf, mytsup;
  double res, partial, err;

  SET_ERROR(LSN_OK);

  tsup = (tsup > model->isup_lifetime ? model->isup_lifetime : tsup);
  tinf = (tinf < model->iinf_lifetime ? model->iinf_lifetime : tinf);

  
  if(tsup <= tinf)
    {
      SET_ERR_MSG(ERR_LIMITS_INFSUP);
      return ERR_LIMITS_INFSUP;
    }

  res = 0;
  imf_plus_param_type params = { IMFi, NULL };
  
  mytinf = tinf;
  if( (type & lsn_SnII) && (tinf <= lsn_mean_lifetime) )
    {
      mytinf = tinf;
      if(tsup < lsn_mean_lifetime)
        {
          mytsup = tsup;
          tsup = 0;
        }
      else
        {
          mytsup = lsn_mean_lifetime;
          tinf = lsn_mean_lifetime;
        }
      
      lsn_F.function = &lsn_nRSnII;
      lsn_F.params = &params;

      if((gsl_status = LSN_INT(&lsn_F, mytinf, mytsup,
			       1e-6, 1e-4,
			       lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
			       &partial, &err)))
        {
	  sprintf(lsn_error_message, "%s :: gsl integration error %i in Sn II (nimf)"
		  " [%.6g - %.6g] : %.6g %.6g (line %d in file %s, at func %s)\n",
		  lsn_error_msgs[ERR_GSL_INTEGRATION],
		  gsl_status, mytinf, mytsup, partial, err, __LINE__, __FILE__, __func__);
	  SET_ERROR(ERR_GSL_INTEGRATION);
	  return ERR_GSL_INTEGRATION;
        }
      res = partial;
    }
  
  if(tsup > lsn_mean_lifetime)
    {
      if(tinf < lsn_mean_lifetime)
        tinf = lsn_mean_lifetime;
      
      if( type & lsn_SnIa)
        {
          lsn_F.function = lsn_nRSnIa;
	  lsn_F.params = &params;
          
          if((gsl_status = LSN_INT(&lsn_F, tinf, tsup,
				   1e-6, 1e-4,
				   lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
				   &partial, &err)))
            {
	      sprintf(lsn_error_message, "%s :: gsl integration error %i in Sn Ia (nimf)"
		      " [%.6g - %.6g] : %.6g %.6g (line %d in file %s, at func %s)\n",
		      lsn_error_msgs[ERR_GSL_INTEGRATION],
		      gsl_status, tinf, tsup, partial, err, __LINE__, __FILE__, __func__);
	      SET_ERROR(ERR_GSL_INTEGRATION);
	      return ERR_GSL_INTEGRATION;	      
            }
          res += partial;
        }
      
      if( type & lsn_AGB)
        {
          lsn_F.function = lsn_nRSnII;
	  lsn_F.params = &params;
	  
          if((gsl_status = LSN_INT(&lsn_F, tinf, tsup,
				   1e-6, 1e-4,
				   lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
				   &partial, &err)))
            {
	      sprintf(lsn_error_message, "%s :: gsl integration error %i in LIMS (nimf)"
		      " [%.6g - %.6g] : %.6g %.6g (line %d in file %s, at func %s)\n",
		      lsn_error_msgs[ERR_GSL_INTEGRATION],
		      gsl_status, tinf, tsup, partial, err, __LINE__, __FILE__, __func__);
	      SET_ERROR(ERR_GSL_INTEGRATION);
	      return ERR_GSL_INTEGRATION;
            }
          res += partial;
        }
      
    }

  return res;
}


int make_time_interval_bins_by_events(int mode, int IMFi, double time, double end_time, double Step_Prec,
                                      double *times, double *delta_times, double *bin_frac)
/*
  
 */
{
  double delta_time;

  double Left, Right, All, timed_frac;

  int i, limit;

  SET_ERROR(LSN_OK);
  
  All = lsn_get_sn_number(mode, IMFi, time, end_time);
  
  if(All == 0)      
    return 0;

  //IMF_Type *IMFp = &IMFs[IMFi];
  i = 0;
  /* initial guess for delta_time */
  if(time > 0)
    modf(log10(time) + 9, &delta_time);
  else
    modf(log10((time + end_time) * Step_Prec * Step_Prec) + 9, &delta_time);
  
  delta_time = pow(10, delta_time) / 1e9;

  /* cycle to calculate time steps */
  while(time < end_time)
    {
      timed_frac = All;
      Left = Right = limit = 0;
      
      while((timed_frac / All < Step_Prec * 0.9 || timed_frac / All > Step_Prec * 1.1) && !limit)
	{
	  if(time + delta_time > end_time)
	    delta_time = end_time - time;

	  timed_frac =
            lsn_get_sn_number(mode, IMFi, time, time+delta_time);

	  if(timed_frac / All < Step_Prec * 0.9)
	    {
              limit = (time + delta_time >= end_time *0.99);
	      if(!limit)
		{
		  Left = lsn_MAX(Left, delta_time);
		  if(Right == 0)
		    delta_time *= 2;
		  else
		    delta_time = (Left + Right) / 2;
		}
	    }
	  if(timed_frac / All > Step_Prec * 1.1)
	    {
	      if(Right == 0)
		Right = delta_time;
	      else
		Right = lsn_MIN(Right, delta_time);
	      if(Left == 0)
		delta_time /= 2;
	      else
		delta_time = (Left + Right) / 2;
	    }
	}
      
      if(timed_frac / All < Step_Prec / 2 && limit)
        {
          delta_times[i - 1] += delta_time;
          if(bin_frac != 0x0)
            bin_frac[i-1] += timed_frac / All;
        }
      else
	{
	  times[i] = time;
	  delta_times[i] = delta_time;
          if(bin_frac != 0x0)
            bin_frac[i] = timed_frac / All;
	  i++;
	}
      time += delta_time;
    }
  return i;
}


/* /------------------------------------------------\  
 * | LIFETIME RELATED FUNCTIONS                     |
 * |                                                |
 * | > dm_dt                                        |
 * | > lifetime                                     |
 * | > dying_mass                                   |
 * | > sec_dist                                     |
 * \------------------------------------------------/ */


double dm_dt_PADOVANI_MATTEUCCI_1993(double m, double t)
{
  if(t >= 0.039765318659064693)
    return -m / t * (1.338 - 0.1116 * (9 + log10(t)));
  else
    return -0.540540540540540540 * m / (t - 0.003);
  return 0;  // never reaches this point
}


double dm_dt_MAEDER_MEYNET_1989(double m, double t)
{
  if(m <= 1.3)
    return -m / t / 0.6545;
  if(m > 1.3 &&
     m <= 3)
    return  -m / t / 3.7;
  if(m > 3 &&
     m <= 7)
    return  -m / t / 2.51;
  if(m > 7 && m <= 15)
    return  -m / t / 1.78;
  if(m > 15 && m <= 53.054)
    return  -m / t / 0.86;
  if(m > 53.054)
    return -0.54054054054 * m / (t - 0.003);
  return 0;  // never reaches this point
}


double dm_dt_generic(double m, double t)
{
  if(t > 0.0302233)
    return -0.37037037037 * m / (t - 0.012);
  else
    return -0.54054054054 * m / (t - 0.003);
  return 0;  // never reaches this point
}


double dm_dt(double m, double t)
{
  return _dm_dt(m, t);
}


double lifetime_PADOVANI_MATTEUCCI_1993(double mass)
/* padovani & matteucci (1993) approach  */
/* move to gibson (1997) one             */
/*                                       */
/* mass is intended in solar units, life */
/* time in Gyr                           */
{
  if(mass < 0.6)
    // note: 0.5549921 is the value below which the argument of the sqrt here below is negative)
    return 160;
  if(mass <= 6.6)
    // original formula is
    //    pow(10, ((1.338 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) / 0.1116 ) - 9);    
    return pow(10, ((1.338 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) * 8.96057347670251) - 9);
  else
    return 1.2 * pow(mass, -1.85) + 0.003;
}


double lifetime_MAEDER_MEYNET_1989(double mass)
/* use the lifetime by
 * maeder & meynet (1989)
 * mass is intended in solar units, life
 * time in Gyr                         
*/

{
  if(mass <= 1.3)
    return pow(10, -0.6545 * log10(mass) + 1);
  if(mass > 1.3 && mass <= 3)
    return pow(10, -3.7 * log10(mass) + 1.35);
  if(mass > 3 && mass <= 7)
    return pow(10, -2.51 * log10(mass) + 0.77);
  if(mass > 7 && mass <= 15)
    return pow(10, -1.78 * log10(mass) + 0.17);
  if(mass > 15 && mass <= 53.054)
    return pow(10, -0.86 * log10(mass) - 0.94);
  if(mass > 53.054)
    return 1.2 * pow(mass, -1.85) + 0.003;
  return 0;  // never reaches this point
}


double lifetime_generic(double mass)
{
  if(mass <= 8)
    return 5 * pow(mass, -2.7) + 0.012;
  else
    return 1.2 * pow(mass, -1.85) + 0.003;
  return 0;  // never reaches this point
}


double lsn_lifetime(double m)
{
  return _lifetime(m);
}


double dying_mass_PADOVANI_MATTEUCCI_1993(double time)
{
  /* if(time < IMFp->inf_lifetime) */
  /*   return IMFp->MU; */
  /* if(time > IMFp->sup_lifetime) */
  /*   return IMFp->Mm; */

  double res = 0;
  int eval_1 = (time > 0.039765318659064811)*2;
  int eval_2 = (time > 0.003);

  switch( eval_1 + eval_2)
    {
    case 3: { double core = 1.338 - 0.1116 * (9 + log10(time));
	res = pow(10, 7.764 - (1.79 - core*core) * 4.48028673835125);}
      break;
    case 1: res = pow((time - 0.003) *0.833333333333, -0.540540540540540);
      break;
    default:
      res = 0;
    }

  return res;

      
  /* if(time > 0.039765318659064811) */
  /*   { */
  /*     // the inverse is: */
  /*     //  */
  /*     double core = 1.338 - 0.1116 * (9 + log10(time)); */
  /*     return pow(10, 7.764 - (1.79 - core*core) * 4.48028673835125); */

  /*     // the following formula is a fit (given by mathemathica) */
      
  /*     /\* double core = log(time); *\/ */
  /*     /\* double core2 = 1.0 / pow(time, 0.3336); *\/ */
  /*     /\* core *= core * 0.0242336; *\/ */
  /*     /\* core += 0.559282; *\/ */
  /*     /\* double core3 = exp(core) * core2; *\/ */
  /*     /\* return core3; *\/ */
  /*   } */
  /* else */
  /*   // note: */
  /*   // 0.54054054054054 = 1/1.85 */
  /*   // 0.83333333333333 = 1/1.2 */
  /*   return pow((time - 0.003) *0.833333333333, -0.540540540540540); */
  
}

double dying_mass_MAEDER_MEYNET_1989(double time)
{
  /* if(time < IMFp->inf_lifetime) */
  /*   return IMFp->MU; */
  /* if(time > IMFp->sup_lifetime) */
  /*   return IMFp->Mm; */

  if(time >= 8.4221714076)
    return pow(10, (1 - log10(time)) / 0.6545);
  if(time < 8.4221714076 &&
     time >= 0.38428316376)
    return pow(10, (1.35 - log10(time)) / 3.7);
  if(time < 0.38428316376 &&
     time >= 0.044545508363)
    return pow(10, (0.77 - log10(time)) / 2.51);
  if(time < 0.044545508363 &&
	 time >= 0.01192772338)
    return pow(10, (0.17 - log10(time)) / 1.78);
  if(time < 0.01192772338 &&
     time >= 0.0037734864318)
    return pow(10, -(0.94 + log10(time)) / 0.86);
  if(time < 0.0037734864318)
    return pow((time - 0.003) / 1.2, -0.54054);
  return 0;  // never reaches this point    
}

double dying_mass_generic(double time)
{
  /* if(time < IMFp->inf_lifetime) */
  /*   return IMFp->MU; */
  /* if(time > IMFp->sup_lifetime) */
  /*   return IMFp->Mm; */

  if(time > 0.0302233)
    return pow((time - 0.012) / 5, -0.37037037037);
  else
    return pow((time - 0.003) / 1.2, -0.54054054054);
  return 0;  // never reaches this point  
}

double lsn_dying_mass(double time)
     /* 
      * calculates mass dying at some time;
      * time is time_in_Gyr                
      */
{
  return _dying_mass(time);  
}


double INLINE_FUNC sec_dist(double gamma, double mu)
     /* 
      * retunrs the secondary distribution function
      * for Sn type Ia
      */
{
  return pow(2, 1 + gamma) * (1 + gamma) *pow(mu, gamma);
}





/* /------------------------------------------------\  
 * | SN II - LOW LEVEL FUNCTIONS                    |
 * |                                                |
 * | > nRSnII                                       |
 * | > mRSnII                                       |
 * | > zmRSnII                                      |
 * | > ztRSnII                                      |
 * | > ejectaSnII                                   |
 * \------------------------------------------------/ */



double lsn_nRSnII(double t, void *p)
     /* 
      * calculates Sn type II explosion rate in #/Gyr
      * at some time t, for a SSP of 1Msun.  
      * t is in Gyr 
      */
{

  if((t >= lsn_iinf_lifetime)&&
     (t <= lsn_isup_lifetime))
    {
      LIMF_FUNCTION *phi = IMFs[((imf_plus_param_type*)p)->IMFi].IMFfunc_byNum;
      //void *params = ((imf_plus_param_type*)p)->params;
      void *p_for_imf = (void*)&IMFs[((imf_plus_param_type*)p)->IMFi];
      double m     = _dying_mass(t);
      double fact  = phi(m, p_for_imf)*(-dm_dt(m, t));
      return fact;
    }
  else
    return 0;
}


double lsn_mRSnII(double t, void *p)
     /* 
      * calculates Sn type II explosion rate in Msun/Gyr
      * at some time t, for a SSP of 1Msun.
      * t is in Gyr 
      */

{
  if((t >= lsn_iinf_lifetime)&&
     (t <= lsn_isup_lifetime) )
    {
      LIMF_FUNCTION *phi = IMFs[((imf_plus_param_type*)p)->IMFi].IMFfunc_byNum;
      //void *params = ((imf_plus_param_type*)p)->params;
      void *p_for_imf = (void*)&IMFs[((imf_plus_param_type*)p)->IMFi];
      double m     = _dying_mass(t);
      double fact  = phi(m, p_for_imf)*(-dm_dt(m, t));
      return m * fact;
    }
  else
    return 0;
}



// ---------------------------------------------------------------------------

double lsn_zmRSnII(double m, void *p)
     /* 
      *	calculates Sn type II explosion rate in Msun/Gyr
      * at some mass m, for a SSP of 1Msun.      
      *
      * NOTE: THIS function accounts for yields that may have
      *       different mass bins for different Z bins
      */
  
{
  int multiply_by_phi = 1;

  IMF_Type *IMFp = &IMFs[((imf_plus_param_type*)p)->IMFi];
  
  if(m < 0)
    // when called by ztRSnII
    {
      m *= -1;
      multiply_by_phi = 0;      
    }

  {
    int BHi;
    for(BHi = 0; (BHi < IMFp->N_notBH_ranges); BHi++)
      if(m <= IMFp->notBH_ranges.sup[BHi] &&
	 m >=  IMFp->notBH_ranges.inf[BHi])
	break;
    
    if(BHi == IMFp->N_notBH_ranges)
      return 0;
  }

  if((m <= IMFp -> MU) &&	/* Mm Msun < m < MUP Msun */
     (m >= IMFp -> Mm ))
    {
      int    ir, ir1, zone = 0;
      double y = 0, y1 = 0, t, t1, u;
      
      // find the mass bin
      for(ir = SD.Mdim - 1; m < SD.MArray[ir] && ir > 0; ir--)
	;

#define Z_LARGE_   4
#define Z_SMALL_   5
#define M_SMALL_   1
#define M_LARGE_1_ 2
#define M_SMALL_1_ 3
      
#define Z_LARGE    (1 << Z_LARGE_)         // 16
#define Z_SMALL    (1 << Z_SMALL_)         // 32
#define M_LARGE    1
#define M_SMALL    (1 << M_SMALL_)         // 2
#define M_NOT_EQU 64
#define M_LARGE_1  (1 << M_LARGE_1_)       // 4
#define M_SMALL_1  (1 << M_SMALL_1_)       // 8
#define Z_SINGLE 128      
      
      if(SD.Zdim > 1)
	{
	  zone |= ((SD.Zstar >= SD.ZArray[SD.Zdim - 1]) << Z_LARGE_);
          zone |= ((SD.Zstar < SD.ZArray[0]           ) << Z_SMALL_);                     

	  zone |=  (m >= SD.MArray[SD.Mdim - 1]);   
	  zone |= ((m  < SD.MArray[0]          ) << M_SMALL_);

          if(SD.Mirr == 1)
            {
              zone |= M_NOT_EQU;
              zone |= ((m >= SD.MArray1[SD.Mdim - 1]) << M_LARGE_1_);
              zone |= ((m <  SD.MArray1[0]          ) << M_SMALL_1_);                
            }
	}
      else
	{
	  zone = Z_SINGLE;
	  zone |=  (m >= SD.MArray[SD.Mdim - 1]);
	  zone |= ((m  < SD.MArray[0]          ) << M_SMALL_);
	}


      // Z value beyond limits
      if(zone & (Z_LARGE + Z_SMALL))
        /* either Z > Zmax or Z < Zmin
         * in both cases we have to deal only with either the first or the last mass array
         * we do not extrapolate on the metallicity, yields are kept flat beyond the table limits
         */
        {
          if(!(zone & 15))
            {
              t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
              y = (1 - t) * SD.Y[SD.Zbin*SD.Mmaxdim + ir] + t * SD.Y[SD.Zbin*SD.Mmaxdim + ir + 1];
            }
          else
            {
              y = SD.Y[SD.Zbin*SD.Mmaxdim + ir];
              if( (zone & 1) && SD.extrap)
                y *= m / SD.MArray[SD.Mdim-1];
	      if( (zone & 2) && SD.extrap_inf)
                y *= m / SD.MArray[0];
            }
        }
      // Z value within limits
      else if(zone < Z_SINGLE)
        {
          // factor for Z interpolation
          u = (SD.Zstar - SD.ZArray[SD.Zbin]) / (SD.ZArray[SD.Zbin + 1] - SD.ZArray[SD.Zbin]);

          // interpolate on the first mass array
          if(!(zone & (M_LARGE + M_SMALL)))
            {
              // mass value within limits
              t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
              y = (1 - t) * SD.Y[SD.Zbin*SD.Mmaxdim + ir] + t * SD.Y[SD.Zbin*SD.Mmaxdim + ir + 1];
            }
          else
            {
              // mass value beyond limits
              y = SD.Y[SD.Zbin*SD.Mmaxdim + ir];

              if( (zone & M_LARGE) && SD.extrap)
                y *= m / SD.MArray[SD.Mdim-1];
	      
	      if( (zone & M_SMALL) && SD.extrap_inf)
                y *= m / SD.MArray[0];
            }

          // interpolate on the second mass array, if that is needed
          if(zone & M_NOT_EQU)
            {
	      for(ir1 = SD.Mdim1 - 1; m < SD.MArray1[ir1] && ir1 > 0; ir1--)
		;

              if(!(zone & (M_LARGE_1 + M_SMALL_1) ) )
                {
                  t1 = (m - SD.MArray1[ir1]) / (SD.MArray1[ir1 + 1] - SD.MArray1[ir1]);
		  y1 = (1 - t1) * SD.Y[(SD.Zbin+1)*SD.Mmaxdim + ir1] + t1 * SD.Y[(SD.Zbin+1)*SD.Mmaxdim + ir1 + 1];                  
                }
              else
                {
                  // mass value beyond limits
                  y1 = SD.Y[(SD.Zbin+1)*SD.Mmaxdim + ir1];

                  if( (zone & M_LARGE_1)  && SD.extrap )
		    y1 *= m / SD.MArray1[SD.Mdim1-1];

                  if( (zone & M_SMALL_1) && SD.extrap_inf)
                    y1 *= m / SD.MArray1[0];
                }
            }
          // mass array is the same
          else
            {
              if( !(zone & (M_LARGE + M_SMALL)) )
                // mass value within limits
                y1 = (1 - t) * SD.Y[(SD.Zbin+1)*SD.Mmaxdim + ir] + t * SD.Y[(SD.Zbin+1)*SD.Mmaxdim + ir + 1];
              else
                {
                  // mass value beyond limits
                  y1 = SD.Y[(SD.Zbin+1)*SD.Mmaxdim + ir];
		  
                  if( (zone & M_LARGE) && SD.extrap)
                    y1 *= m / SD.MArray[SD.Mdim-1];
		  
                  if( (zone & M_SMALL) && SD.extrap_inf)
                    y1 *= m / SD.MArray[0];
                }
            }

          // interpolate on metallicity
          y = (1 - u) * y + u * y1;
          
        }
      else
        {
          if(zone & M_LARGE)          /* M > M_max */
            {
              y = SD.Y[SD.Zbin*SD.Mmaxdim + SD.Mdim - 1];
              if(SD.extrap)
                y *= (m / SD.MArray[SD.Mdim-1]);
            }
          else if(zone & M_SMALL)     /* M < M_min */
	    {
	      y = SD.Y[SD.Zbin*SD.Mmaxdim];
	      if(SD.extrap_inf)
		y *= (m / SD.MArray[0]);
	    }
          else
            {
              t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
              y = (1 - t) * SD.Y[SD.Zbin*SD.Mmaxdim + ir] + t * SD.Y[SD.Zbin*SD.Mmaxdim + ir + 1];
            }
        }
      
      if(multiply_by_phi)
	{
	  LIMF_FUNCTION *phi = IMFp->IMFfunc_byNum;
	  void *params       = (void*)IMFp;
	  return y * phi(m, params);
	}
      else
        return y;
    }
  else
    return 0;
}



double lsn_ejectaSnII(double m, void *p)
  /* 
   * calculates Sn type II restored mass in Msun
   * for a star of mass m
   */
{

  int BHi;
  IMF_Type *IMFp = &IMFs[((imf_plus_param_type*)p)->IMFi];
  
  for(BHi = 0; (BHi < IMFp->N_notBH_ranges); BHi++)
    if(m <= IMFp->notBH_ranges.sup[BHi] &&
       m >=  IMFp->notBH_ranges.inf[BHi])
      break;
  
  if(BHi == IMFp->N_notBH_ranges)
    return 0;

  if( m >= IMFp -> Mm && m <= IMFp -> MU )
    {
      int ir;
      double mr;
      if(getindex((double *) &SD.MArray[0], 0, SD.Mdim - 1, &m, &ir) < 0)
	{
	  mr = SD.Y[SD.Zbin * SD.Mmaxdim + ir];
	  if(m < SD.MArray[0])
	    mr = SD.Y[SD.Zbin * SD.Mmaxdim] * (m / SD.MArray[0]);
	}
      else
	{
	  double t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  mr       = (1 - t) * SD.Y[SD.Zbin * SD.Mmaxdim + ir] + t * SD.Y[SD.Zbin * SD.Mmaxdim + ir + 1];
	}

      LIMF_FUNCTION *phi = IMFp->IMFfunc_byNum;
      void *params       = (void*)IMFp;
      
      return mr * phi(m, params);
    }
  else
    
    return 0;
}


double lsn_ztRSnII(double time, void *p)
     /* 
      * return the production rate in mass for elements at time t
      */
{

  IMF_Type *IMFp = &IMFs[((imf_plus_param_type*)p)->IMFi];
  
  if((time >= IMFp -> inf_lifetime) &&	/* Mm Msun < m < MUP Msun */
     (time <= IMFp -> sup_lifetime))
    {
      double m = _dying_mass(time);

      LIMF_FUNCTION *phi = IMFp->IMFfunc_byNum;
      void *params       = (void*)IMFp;
      
      double fact = phi(m, params)*(-dm_dt(m, time));

      return lsn_zmRSnII(-m, params) * fact;
    }
  else
    return 0;
}




/* /------------------------------------------------\  
 * | SN Ia - LOW LEVEL FUNCTIONS                    |
 * |                                                |
 * | > nRSnIa                                       |
 * \------------------------------------------------/ */


/* : -------------------------------------------- :
   :                                              :
   : General model as formultaed in Greggio 2005  :
   : -------------------------------------------- :
 */

double initialize_Ia(int IMFi, double tau_i, double tau_x)
/*
  returns the number of SnIa events over tau_i - tau_x
 */
{
  double A, res, err;
  int gsl_status;;

  lsn_F.function = lsn_DTDfunction;
  if((gsl_status = LSN_INT(&lsn_F, tau_i, tau_x,
                                       1e-6, 1e-4,
                                       lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
                                       &A, &err)))
    {
      dprintf_err(VXERR, 0, "%s  >> gsl integration error %i in Sn (nimf)"
		  " [%.6g - %.6g] : %.6g %.6g\n", lsn_err_prefix, 
		  gsl_status, tau_i, tau_x, A, err);
      exit(1);
    }
  lsn_DTDfunction(-A, 0x0);                                                  /* normalize DTD itself to 1 */

  lsn_F.function = lsn_DTDfunction;
  LSN_INT(&lsn_F, tau_i, tau_x,
                      1e-6, 1e-4,
                      lsn_gsl_intspace_dim, lsn_gsl_intkey, lsn_w,
                      &A, &err);
  
  lsn_k_alpha = limf_IntegrateIMF_byNum(IMFs[IMFi].Mm, IMFs[IMFi].MU, &IMFs[IMFi], INC_BH);

  lsn_BinFrac /= lsn_k_alpha;                                            /* the lsn_BinFrac is defined in the paramfile as
                                                                            the number of Ia events per solar mass of formed stars.
                                                                            this way it becomes, correctly, as the number of snIa
                                                                            events per number of stars formed.
                                                                         */

  res = lsn_get_sn_number(lsn_SnIa, IMFi, tau_i, tau_x);

  return res;
}

double nRSnIa_DTD(double t, void *params)
{
  DONT_WARN_UNUSED(params);
  
  if(t > lsn_isup_lifetime ||
     t < lsn_mean_lifetime)
    return 0;
  
  return lsn_k_alpha * lsn_BinFrac * lsn_DTDfunction(t, 0x0);
}




/* DTD for SD scenario as in Matteucci&Recchi 01, as stated in Bonaparte et al. 2013 */
double DTD_SD(double t, void *params)
{
  static double normalization = 1.0;
  static double gamma = 0.5;
  static double s = 2.35;
  static double Mm = 3.0;
  static double MM = 16.0;
  DONT_WARN_UNUSED(params);
  
  if(t < 0 )
    {
      normalization = -1.0/t;
      return normalization;
    }
  if(params != 0x0)
    {
      gamma = *(double*)params;
      s     = *((double*)params + 1);
      Mm    = *((double*)params + 2);
      MM    = *((double*)params + 3);
    }

  double M2      = _dying_mass(t);
  double s_gamma = -(s + gamma);
  double gamma_1 = (gamma + 1);

  double Mb = lsn_MAX(2 * M2, Mm);
  double MB = M2 + 0.5*MM;

  double f;
  
  f = normalization * pow(2, gamma_1) * gamma_1 * pow(M2, gamma) * (pow(Mb, s_gamma) - pow(MB, s_gamma)) / (s_gamma);
  
  return f * dm_dt(M2, t);
}


/* DTD for SD scenario as in Matteucci&Recchi 01, as stated in Bonaparte et al. 2013 */
double DTD_SD_oldstyle(double t, void *params)
{
  static double normalization = 1.0;
  static double gamma = 2.0;
  static double s = 2.35;
  static double Mm = 3.0;
  static double MM = 16.0;
  DONT_WARN_UNUSED(params);
  
  if(t < 0 )
    {
      normalization = -1.0/t;
      return normalization;
    }

  double M2      = _dying_mass(t);
  double s_gamma = -(s + gamma);
  double gamma_1 = (gamma + 1);

  double Mb = lsn_MAX(2 * M2, Mm);
  double MB = M2 + 0.5*MM;

  double f;
  
  f = normalization * pow(2, gamma_1) * gamma_1 * pow(M2, gamma) * (pow(Mb, s_gamma) - pow(MB, s_gamma)) / (s_gamma);
  
  return f * dm_dt(M2, t);
}


/* DTD by Mannucci, Della Valle and Panagia, 2006 */
double DTD_Mannucci_DellaValle_Panagia_2006(double t, void *params)
{
  static double normalization = 1.0;
  double exponent;
  double exponent_;
  DONT_WARN_UNUSED(params);
  
  if(t < 0)
    {
      normalization = -1.0/t;
      return normalization;
    }
  
  exponent_ = log10(t);
  
  if(t >= 0.0851)
    exponent = -0.71 - 0.9 * (exponent_ + 0.3) * (exponent_ + 0.3);
  else
    exponent = 1.4 - 50 * (exponent_ + 1.3) * (exponent_ + 1.3);

  return (normalization * pow(10, exponent));
}


/* Strolger et al 04, as reported by Yates et al 13 */
double DTD_gaussian(double t, void*params)
{
#define PI 3.141592654
  static double normalization = 1.0;
  const double tau_c = 1.0;  // Gyr
  const double sigma_tau2 = (0.2 * tau_c) * (0.2 * tau_c);
  const double factor = 1.0 / (sqrt(2 * PI * sigma_tau2));
  DONT_WARN_UNUSED(params);
  
  if(t < 0 )
    {
      normalization = -1.0/t;
      return normalization;
    }
  return normalization * factor * exp(-(t - tau_c)*(t - tau_c) / (2*sigma_tau2));
}

/* Maoz, Mannucci and Brandt 2012 */
double DTD_powerlaw(double t, void*params)
{
  static double normalization = 1.0;
  const double exponent = -1.12;
  DONT_WARN_UNUSED(params);
  
  if(t < 0 )
    {
      normalization = -1.0/t;
      return normalization;
    }
  
  return normalization * pow(t, exponent);
}



/* : ------------------------------ :
   :                                :
   : SD model by Greggio & Renzini  :
   : -------------------------------:
 */
double inner_integrand(double m, void *p)
{
  double m2 = m * m;
  LIMF_FUNCTION *phi = IMFs[((imf_plus_param_type*)p)->IMFi].IMFfunc_byNum;
  void *params       = ((imf_plus_param_type*)p)->params;
  DONT_WARN_UNUSED(params);
  
  return (phi(m, params) / (m2 * m2));
}

double nRSnIa_SD_Greggio_and_Renzini(double t, void *params)
  /* 
   * calculates Sn type Ia explosion rate in #/yr
   * at some time t, for a SSP of 1Msun.             
   * in snIa_mass return the total mass in SnIa.     
   * in Rejecta return the rate of ejecta in Msun/yr 
   * t is in Gyr 
   */
{
  double result;
  double m2, Mb_inf, Mb_sup;
  DONT_WARN_UNUSED(params);

  IMF_Type *IMFp = &IMFs[((imf_plus_param_type*)params)->IMFi];
  
  m2 = _dying_mass(t);

  if((m2 < lsn_iMi) || (m2 > lsn_Mup))
    return 0;

  Mb_inf = lsn_MAX(2 * m2, lsn_MBm);
  Mb_sup = 0.5 * lsn_MBM + m2;


  if(IMFp->type == power_law)
    {
      int sbin1, sbin2;
      
      if(IMFp->NSlopes == 1)
	sbin1 = sbin2 = 0;
      else
	{
	  int IMFi = ((imf_plus_param_type*)params)->IMFi;
	  sbin1 = get_IMF_SlopeBin(IMFp, Mb_sup);
	  sbin2 = get_IMF_SlopeBin(IMFp, Mb_inf);
	}

      if(sbin1 == sbin2)
	{
	  double fact2 = 1.0 / (3 + IMFp->Slopes.slopes[sbin1]);
	  result = ((fact2 * pow(Mb_inf, -(3 + IMFp->Slopes.slopes[sbin1]))) -
		    (fact2 * pow(Mb_sup, -(3 + IMFp->Slopes.slopes[sbin1])))) * IMFp -> A[0];
	}
      else
	{
	  double fact1 = Mb_inf;
	  for(int i = sbin1, result = 0; i <= sbin2; i++)
	    {
	      if(i > sbin1)
		Mb_sup = IMFp->Slopes.masses[i];

	      if(i < sbin2)
		Mb_sup = IMFp->Slopes.masses[i + 1];
	      else
		Mb_sup = fact1;

	      double fact2 = 1.0 / (3 + IMFp -> Slopes.slopes[i]);
	      result += ( (fact2 * pow(Mb_inf, -(3 + IMFp -> Slopes.slopes[i]))) -
			  (fact2 * pow(Mb_sup, -(3 + IMFp -> Slopes.slopes[i]))) ) * IMFp -> A[i];
	    }
	}
    }
  else
    {
      gsl_function localF;
      localF.function = &inner_integrand;
      localF.params = params;

      int    gsl_status;
      size_t n_eval;
      double err;
      if((gsl_status = gsl_integration_qng(&localF, Mb_inf, Mb_sup, 1, 1e-5, &result, &err, &n_eval)))
	{
	  dprintf_err(VXERR, 0, "%s gsl integration error %i in Sn (k) [%.6g - %.6g]: %.6g %.6g\n", lsn_err_prefix, 
		 gsl_status, Mb_inf, Mb_sup, result, err);
	  fflush(stdout);
	}
      result *= IMFp -> A[0]; /* if NSlopes > 1 this should appear inside the inner_integrand function */
    }

  double fact = lsn_BinFrac * 24 * m2 * m2 * (-dm_dt(m2, t));
  return ( fact * result );
}


int INLINE_FUNC getindex(double *table, int lower, int greatest, double *tval, int *index)
{
  double vh, maxval, minval;

  int found, half, upi, downi, res;

  downi = 0;

  maxval = table[greatest];
  minval = table[lower];

  upi = greatest;
  found = 0;

  if(*tval >= maxval)
    {
      *index = greatest;
      /**tval = maxval;*/
      res = -1;
    }
  else if(*tval <= minval)
    {
      *index = lower;
      /**tval = minval;*/
      res = -2;
    }
  else
    {
      while(found != 1)
	{
	  /* non-interpolated binary guess          
	   * half = (upi + downi) / 2;
	   */

	  /* linearly interpolate median index.
	   * this work if values are equally linearly distributed as is our case.
	   */

	  /*if((half = downi + (int)((*tval - table[downi]) * (upi - downi) / (table[upi] - table[downi]))) != downi) */
	  if((half = (upi + downi) / 2) != downi)
	    {
	      vh = table[half];

	      if(*tval < vh)
		upi = half;
	      else if(*tval > vh)
		downi = half;
	      else
		{
		  downi = half;
		  found = 1;
		}
	      if((upi - downi) == 1)
		found = 1;
	    }
	  else
	    found = 1;
	}

      res = 0;
      *index = downi;
    }
  return (res);
}


void *manage_memory(size_t bytes)
// initialize the memory used by the library
{
  if( Memory == NULL )
    {
      if( (Memory = (char*)aligned_alloc(ALIGN, bytes)) == NULL )
	return NULL;
      Memory_top    = Memory;
      Memory_amount = bytes;
      Memory_roof   = Memory + bytes;
    }

  else
    free(Memory);
  
  return Memory;
}

 void *_allocate(size_t bytes, int alignment)
 {
   char *ptr;
   
   // --------------------------------------
   // keep memory aligned
   
   // add space to store the amount of allocated bytes;
   bytes += MEM_T_SIZE;
   while( bytes % alignment )
     bytes++;
   
   while( (Memory_top - Memory) % alignment )
     Memory_top++;
   
   //
   // --------------------------------------
   
   if( Memory_top + bytes > Memory_roof )
     // not wnough memory available
     ptr = NULL;
   
   else
     {
       *(mem_t*)Memory_top = bytes - MEM_T_SIZE;
       ptr = (void*)(Memory_top + MEM_T_SIZE);
       Memory_top += bytes;
     }

   return (void*)ptr;
 }

 void *_callocate(size_t bytes, int alignment)
 {   
   void *ptr = _allocate(bytes, alignment);

   if(ptr == NULL)
     return NULL;

   for(size_t i = 0; i < bytes; i++)
     *(char*)(ptr + i) = 0;
   
   return ptr;
 }

 
 void *_reallocate(char *ptr, size_t bytes, int alignment)
 {
   mem_t old_bytes = *(mem_t*)(ptr - MEM_T_SIZE);

   if( Memory_roof == ptr + old_bytes )
     // this was anyway the last block, so we can just let it grow
     {
       while( bytes % alignment )
	 bytes++;
       *(mem_t*)(ptr - MEM_T_SIZE) = bytes;

       Memory_top = ptr + bytes;

       return ptr;
     }
   
   while( (Memory_top - Memory) % alignment)
     Memory_top++;
   
   if( (Memory_top + bytes + MEM_T_SIZE > Memory_roof) ||
       (bytes <= old_bytes) )
     return NULL;

   *(mem_t*)Memory_top = bytes;
   Memory_top += MEM_T_SIZE;
   char *myptr = Memory_top;

   memcpy(Memory_top, ptr, old_bytes);

   Memory_top += bytes;

   return (void*)myptr;
 }
 
 

void *_release(char *ptr)
{
  if( (ptr > Memory_roof) || (ptr - MEM_T_SIZE < Memory) )
    return NULL;
      
  else
    Memory_top = ptr - MEM_T_SIZE;
  
  return (void*)Memory_top;
}

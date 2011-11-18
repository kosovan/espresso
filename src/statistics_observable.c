#include "statistics_observable.h"
#include "particle_data.h"
#include "parser.h"
#include "integrate.h"
#include "lb.h"
#include "pressure.h"

observable** observables = 0;
int n_observables = 0; 


static int convert_types_to_ids(IntList * type_list, IntList * id_list); 
  
int tclcommand_parse_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, profile_data** pdata_);
int tclcommand_parse_radial_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, radial_profile_data** pdata);
int sf_print_usage(Tcl_Interp* interp);


int tclcommand_observable_print(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  char buffer[TCL_DOUBLE_SPACE];
  double* values=malloc(obs->n*sizeof(double));
  (*obs->fun)(obs->args, values, obs->n);
  if (argc==0) {
    for (int i = 0; i<obs->n; i++) {
      Tcl_PrintDouble(interp, values[i], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
    }
  } else if (argc>0 && ARG0_IS_S("formatted")) {
    tclcommand_observable_print_formatted(interp, argc-1, argv+1, change, obs, values);
  } else {
    Tcl_AppendResult(interp, "Unknown argument to observable print: ", argv[0], "\n", (char *)NULL );
    return TCL_ERROR;
  }
  free(values);
  return TCL_OK;
}

int tclcommand_observable_print_parameters(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  int i;
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 3];
  double* values=malloc(obs->n*sizeof(double));
  (*obs->fun)(obs->args, values, obs->n);

  if (argc==0) {
    sprintf(buffer," %d ",obs->n);
    Tcl_AppendResult(interp, obs->fun_name, buffer, (char *)NULL );
    for (i = 0; i<obs->n; i++) {
      Tcl_PrintDouble(interp, values[i], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
    }
  } else { 
    Tcl_AppendResult(interp, "Unknown argument to observable print_parameters: ", argv[0], "\n", (char *)NULL );
    return TCL_ERROR;
  }
  free(values);
  return TCL_OK;
}

int tclcommand_observable_print_formatted(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs, double* values) {
  if (0) {
#ifdef LB
  } else if (obs->fun == (&observable_lb_velocity_profile)) {
    return tclcommand_observable_print_profile_formatted(interp, argc, argv, change, obs, values, 3, 0);
#endif
  } else if (obs->fun == (&observable_density_profile)) {
    return tclcommand_observable_print_profile_formatted(interp, argc, argv, change, obs, values, 1, 1);
#ifdef LB
  } else if (obs->fun == (&observable_lb_radial_velocity_profile)) {
    return tclcommand_observable_print_radial_profile_formatted(interp, argc, argv, change, obs, values, 3, 0);
#endif
  } else if (obs->fun == (&observable_radial_density_profile)) {
    return tclcommand_observable_print_radial_profile_formatted(interp, argc, argv, change, obs, values, 1, 1);
  } else if (obs->fun == (&observable_radial_flux_density_profile)) {
    return tclcommand_observable_print_radial_profile_formatted(interp, argc, argv, change, obs, values, 3, 1);
  } else if (obs->fun == (&observable_flux_density_profile)) {
    return tclcommand_observable_print_profile_formatted(interp, argc, argv, change, obs, values, 3, 1);
  } else { 
    Tcl_AppendResult(interp, "Observable can not be printed formatted\n", (char *)NULL );
    return TCL_ERROR;
  }

}

int tclcommand_observable_print_profile_formatted(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs, double* values, int groupsize, int shifted) {
  profile_data* pdata=(profile_data*) obs->args;
  char buffer[TCL_DOUBLE_SPACE];
  double data;
  int linear_index;
  double x_offset, y_offset, z_offset, x_incr, y_incr, z_incr;
  if (shifted) {
    x_incr = (pdata->maxx-pdata->minx)/(pdata->xbins);
    y_incr = (pdata->maxy-pdata->miny)/(pdata->ybins);
    z_incr = (pdata->maxz-pdata->minz)/(pdata->zbins);
    x_offset = pdata->minx + 0.5*x_incr;
    y_offset = pdata->miny + 0.5*y_incr;
    z_offset = pdata->minz + 0.5*z_incr;
  } else {
    x_incr = (pdata->maxx-pdata->minx)/(pdata->xbins-1);
    y_incr = (pdata->maxy-pdata->miny)/(pdata->ybins-1);
    z_incr = (pdata->maxz-pdata->minz)/(pdata->zbins-1);
    x_offset = pdata->minx;
    y_offset = pdata->miny;
    z_offset = pdata->minz;
  }
  for (int i = 0; i < pdata->xbins; i++) 
    for (int j = 0; j < pdata->ybins; j++) 
      for (int k = 0; k < pdata->zbins; k++) {
        if (pdata->xbins>1) {
          Tcl_PrintDouble(interp, x_offset + i*x_incr , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        } else {
          Tcl_PrintDouble(interp, 0 , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }

        if (pdata->ybins>1) {
          Tcl_PrintDouble(interp,y_offset + j*y_incr, buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        } else {
          Tcl_PrintDouble(interp, 0 , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }
        if (pdata->zbins>1) {
          Tcl_PrintDouble(interp,z_offset + k*z_incr, buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        } else {
          Tcl_PrintDouble(interp, 0 , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }
        linear_index = 0;
        if (pdata->xbins > 1)
          linear_index += i*pdata->ybins*pdata->zbins;
        if (pdata->ybins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;

        for (int l = 0; l<groupsize; l++) {
          data=values[groupsize*linear_index+l];
          Tcl_PrintDouble(interp, data , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }
        Tcl_AppendResult(interp, "\n", (char *)NULL );
      }
  return TCL_OK;
//  Tcl_AppendResult(interp, "\n", (char *)NULL );
}

int tclcommand_observable_print_radial_profile_formatted(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs, double* values, int groupsize, int shifted) {
  radial_profile_data* pdata=(radial_profile_data*) obs->args;
  char buffer[TCL_DOUBLE_SPACE];
  double data;
  int linear_index;
  double r_offset, phi_offset, z_offset, r_incr, phi_incr, z_incr;
  if (shifted) {
    r_incr = (pdata->maxr-pdata->minr)/(pdata->rbins);
    phi_incr = (pdata->maxphi-pdata->minphi)/(pdata->phibins);
    z_incr = (pdata->maxz-pdata->minz)/(pdata->zbins);
    r_offset = pdata->minr + 0.5*r_incr;
    phi_offset = pdata->minphi + 0.5*phi_incr;
    z_offset = pdata->minz + 0.5*z_incr;
  } else {
    r_incr = (pdata->maxr-pdata->minr)/(pdata->rbins-1);
    phi_incr = (pdata->maxphi-pdata->minphi)/(pdata->phibins-1);
    z_incr = (pdata->maxz-pdata->minz)/(pdata->zbins-1);
    r_offset = pdata->minr;
    phi_offset = pdata->minphi;
    z_offset = pdata->minz;
  }
  for (int i = 0; i < pdata->rbins; i++) 
    for (int j = 0; j < pdata->phibins; j++) 
      for (int k = 0; k < pdata->zbins; k++) {
        if (pdata->rbins>1) {
          Tcl_PrintDouble(interp, r_offset + i*r_incr , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        } else {
          Tcl_PrintDouble(interp, 0 , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }

        if (pdata->phibins>1) {
          Tcl_PrintDouble(interp, phi_offset + j*phi_incr , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        } else {
          Tcl_PrintDouble(interp, 0 , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }
        if (pdata->zbins>1) {
          Tcl_PrintDouble(interp, z_offset + k*z_incr , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        } else {
          Tcl_PrintDouble(interp, 0 , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }
        linear_index = 0;
        if (pdata->rbins > 1)
          linear_index += i*pdata->phibins*pdata->zbins;
        if (pdata->phibins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;
        
        for (int j = 0; j<groupsize; j++) {
          data=values[groupsize*linear_index+j];
          Tcl_PrintDouble(interp, data , buffer);
          Tcl_AppendResult(interp, buffer, " ", (char *)NULL );
        }
        Tcl_AppendResult(interp, "\n", (char *)NULL );
      }
  return TCL_OK;
}

int tclcommand_observable_particle_velocities(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp;
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
    return TCL_ERROR;
  obs->fun=&observable_particle_velocities;
  obs->fun_name=strdup("particle_velocities");
  obs->args=ids;
  obs->n=3*ids->n;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_com_velocity(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp, blocksize;
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
    return TCL_ERROR;
  argc-=temp+1;
  argv+=temp+1;
  for ( int i = 0; i < argc; i++) {
    printf("%s\n", argv[i]);
  }
  if (argc>0 && ARG0_IS_S("blocked")) {
    if (argc >= 2 && ARG1_IS_I(blocksize) && (ids->n % blocksize ==0 )) {
      obs->fun=&observable_blocked_com_velocity;
      obs->fun_name=strdup("blocked_com_velocity");
      obs->args=ids;
      obs->n=3*ids->n/blocksize;
      *change=3+temp;
      printf("found %d ids and a blocksize of %d, that makes %d dimensions\n", ids->n, blocksize, obs->n);
      return TCL_OK;
    } else {
      Tcl_AppendResult(interp, "com_velocity blocked expected integer argument that fits the number of particles\n", (char *)NULL );
      return TCL_ERROR;
    }
  } else /* if nonblocked com is to be taken */ {
    obs->fun=&observable_com_velocity;
    obs->fun_name=strdup("com_velocity");
    obs->args=ids;
    obs->n=3;
    *change=1+temp;
    return TCL_OK;
  }
}

int tclcommand_observable_com_position(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp, blocksize;
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
    return TCL_ERROR;
  argc-=temp+1;
  argv+=temp+1;
  for ( int i = 0; i < argc; i++) {
    printf("%s\n", argv[i]);
  }
  if (argc>0 && ARG0_IS_S("blocked")) {
    if (argc >= 2 && ARG1_IS_I(blocksize) && (ids->n % blocksize ==0 )) {
      obs->fun=&observable_blocked_com_position;
      obs->fun_name=strdup("blocked_com_position");
      obs->args=ids;
      obs->n=3*ids->n/blocksize;
      *change=3+temp;
      printf("found %d ids and a blocksize of %d, that makes %d dimensions\n", ids->n, blocksize, obs->n);
      return TCL_OK;
    } else {
      Tcl_AppendResult(interp, "com_velocity blocked expected integer argument that fits the number of particles\n", (char *)NULL );
      return TCL_ERROR;
    }
  } else /* if nonblocked com is to be taken */ {
    obs->fun=&observable_com_position;
    obs->fun_name=strdup("com_position");
    obs->args=ids;
    obs->n=3;
    *change=1+temp;
    return TCL_OK;
  }
}

int tclcommand_observable_particle_positions(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList* ids;
  int temp;
  if (parse_id_list(interp, argc-1, argv+1, &temp, &ids) != TCL_OK ) 
     return TCL_ERROR;
  obs->fun = &observable_particle_positions;
  obs->fun_name=strdup("particle_positions");
  obs->args=(void*)ids;
  obs->n=3*ids->n;
  *change=1+temp;
  return TCL_OK;
}


int tclcommand_observable_stress_tensor(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  obs->fun = &observable_stress_tensor;
  obs->fun_name=strdup("stress_tensor");
  obs->args=(void*)NULL;
  obs->n=9;
  *change=1;
  return TCL_OK;
}


int tclcommand_observable_stress_tensor_acf_obs(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  obs->fun = &observable_stress_tensor_acf_obs;
  obs->fun_name=strdup("stress_tensor_acf_obs");
  obs->args=(void*)NULL;
  obs->n=6;
  *change=1;
  return TCL_OK;
}


int tclcommand_observable_density_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
  profile_data* pdata;
  obs->fun = &observable_density_profile;
  obs->fun_name=strdup("density_profile");
  if (! tclcommand_parse_profile(interp, argc-1, argv+1, &temp, &obs->n, &pdata) == TCL_OK ) 
    return TCL_ERROR;
  if (pdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in radial_profile: particle ids/types not specified\n" , (char *)NULL);
    return TCL_ERROR;
  }
  obs->args=(void*)pdata;
  obs->n=pdata->xbins*pdata->ybins*pdata->zbins;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_lb_velocity_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
#ifndef LB
  return TCL_ERROR;
#else
  profile_data* pdata;
  obs->fun = &observable_lb_velocity_profile;
  obs->fun_name=strdup("lb_velocity_profile");
  if (! tclcommand_parse_profile(interp, argc-1, argv+1, &temp, &obs->n, &pdata) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)pdata;
  obs->n=3*pdata->xbins*pdata->ybins*pdata->zbins;
  *change=1+temp;
  return TCL_OK;
#endif
}


int tclcommand_observable_radial_density_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
  radial_profile_data* rpdata;
  obs->fun = &observable_radial_density_profile;
  obs->fun_name=strdup("radial_dnesity_profile");
  if (! tclcommand_parse_radial_profile(interp, argc-1, argv+1, &temp, &obs->n, &rpdata) == TCL_OK ) 
     return TCL_ERROR;
  if (rpdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in radial_profile: particle ids/types not specified\n" , (char *)NULL);
    return TCL_ERROR;
  }
  obs->args=(void*)rpdata;
  obs->n=rpdata->rbins*rpdata->phibins*rpdata->zbins;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_radial_flux_density_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
  radial_profile_data* rpdata;
  obs->fun = &observable_radial_flux_density_profile;
  obs->fun_name=strdup("radial_flux_dnesity_profile");
  if (! tclcommand_parse_radial_profile(interp, argc-1, argv+1, &temp, &obs->n, &rpdata) == TCL_OK ) 
     return TCL_ERROR;
  if (rpdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in radial_profile: particle ids/types not specified\n" , (char *)NULL);
    return TCL_ERROR;
  }
  obs->args=(void*)rpdata;
  obs->n=3*rpdata->rbins*rpdata->phibins*rpdata->zbins;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_flux_density_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
  profile_data* pdata;
  obs->fun = &observable_flux_density_profile;
  obs->fun_name=strdup("flux_dnesity_profile");
  if (! tclcommand_parse_profile(interp, argc-1, argv+1, &temp, &obs->n, &pdata) == TCL_OK ) 
     return TCL_ERROR;
  if (pdata->id_list==0) {
    Tcl_AppendResult(interp, "Error in radial_profile: particle ids/types not specified\n" , (char *)NULL);
    return TCL_ERROR;
  }
  obs->args=(void*)pdata;
  obs->n=3*pdata->xbins*pdata->ybins*pdata->zbins;
  *change=1+temp;
  return TCL_OK;
}

int tclcommand_observable_lb_radial_velocity_profile(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
  int temp;
#ifndef LB
  return TCL_ERROR;
#else

  radial_profile_data* rpdata;
  obs->fun = &observable_lb_radial_velocity_profile;
  obs->fun_name=strdup("lb_radial_velocity_profile");
  if (! tclcommand_parse_radial_profile(interp, argc-1, argv+1, &temp, &obs->n, &rpdata) == TCL_OK ) 
     return TCL_ERROR;
  obs->args=(void*)rpdata;
  obs->n=3*rpdata->rbins*rpdata->phibins*rpdata->zbins;
  *change=1+temp;
  return TCL_OK;
#endif
}

int tclcommand_observable_particle_currents(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs){
#ifdef ELECTROSTATICS
  int temp;
  IntList* ids;
  obs->fun = &observable_particle_currents;
  obs->fun_name=strdup("particle_currents");
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)ids;
  obs->n=3*ids->n;
  *change=1+temp;
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "Feature ELECTROSTATICS needed for observable particle_currents", (char *)NULL);
  return TCL_ERROR;
#endif
  
}

int tclcommand_observable_currents(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
#ifdef ELECTROSTATICS
  int temp;
  IntList* ids;
  obs->fun = &observable_currents;
  obs->fun_name=strdup("currents");
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)ids;
  obs->n=3;
  *change=1+temp;
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "Feature ELECTROSTATICS needed for observable currents", (char *)NULL);
  return TCL_ERROR;
#endif
} 

int tclcommand_observable_dipole_moment(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
#ifdef ELECTROSTATICS
  int temp;
  IntList* ids;
  obs->fun = &observable_dipole_moment;
  obs->fun_name=strdup("dipole_moment");
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids) == TCL_OK ) 
    return TCL_ERROR;
  obs->args=(void*)ids;
  obs->n=3;
  *change=1+temp;
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "Feature ELECTROSTATICS needed for observable currents", (char *)NULL);
  return TCL_ERROR;
#endif
}

int tclcommand_observable_structure_factor(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  int order;
  int* order_p;

  Tcl_AppendResult(interp, "Structure Factor not available yet!!", (char *)NULL);
  return TCL_ERROR;
  if (argc > 1 && ARG1_IS_I(order)) {
    obs->fun = &observable_structure_factor;
    obs->fun_name=strdup("structure_factor");
    order_p=malloc(sizeof(int));
    *order_p=order;
    obs->args=(void*) order_p;
    int order2,i,j,k,l,n ; 
    order2=order*order;
    l=0;
    // lets counter the number of entries for the DSF
    for(i=-order; i<=order; i++) 
      for(j=-order; j<=order; j++) 
        for(k=-order; k<=order; k++) {
          n = i*i + j*j + k*k;
          if ((n<=order2) && (n>=1)) 
            l=l+2;
  }
    obs->n=l;
    *change=2;
    return TCL_OK;
  } else { 
    sf_print_usage(interp);
    return TCL_ERROR; 
  }
}

int tclcommand_observable_interacts_with(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList *ids1, *ids2;
  int temp;
  double cutoff;
  obs->fun = &observable_interacts_with;
  obs->fun_name=strdup("interacts_with");
  ids1=(IntList*)malloc(sizeof(IntList));
  ids2=(IntList*)malloc(sizeof(IntList));
  iw_params* iw_params_p=(iw_params*) malloc(sizeof(iw_params));
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids1) == TCL_OK ) {
    free(ids1);
    free(ids2);
    free(iw_params_p);
    return TCL_ERROR;
  }
  iw_params_p=(iw_params*)malloc(sizeof(iw_params));
  iw_params_p->ids1=ids1;
  *change=1+temp;
  if (! parse_id_list(interp, argc-3, argv+3, &temp, &ids2) == TCL_OK ) {
    free(ids1);
    free(ids2);
    free(iw_params_p);
    return TCL_ERROR;
  }
  *change+=temp;
  iw_params_p->ids2=ids2;
  if ( argc < 5 || !ARG_IS_D(5,cutoff)) {
    Tcl_AppendResult(interp, "aUsage: analyze correlation ... interacts_with id_list1 id_list2 cutoff", (char *)NULL);
    free(ids1);
    free(ids2);
    free(iw_params_p);
  return TCL_ERROR;
  } 
  *change+=1;
  iw_params_p->cutoff=cutoff;
  obs->args=(void*)iw_params_p;
  obs->n=ids1->n; // number of ids from the 1st argument
  return TCL_OK;
}


int tclcommand_observable_nearest_neighbour_conditional(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList *ids1, *ids2, *partners, *conditions;
  char usage_msg[]="Usage: analyze correlation ... nearest_neighbour_conditional id_list1 id_list2 cutoff chain_length";
  int temp, i;
  double cutoff;
  int chain_length;
  obs->fun = &observable_nearest_neighbour_conditional;
  obs->fun_name=strdup("nearest_neighbour_conditional");
  nn_cond_params* nn_cond_params_p=(nn_cond_params*) malloc(sizeof(nn_cond_params));
  ids1=(IntList*)malloc(sizeof(IntList));
  ids2=(IntList*)malloc(sizeof(IntList));
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids1) == TCL_OK ) {
    free(ids1);
    free(ids2);
    free(nn_cond_params_p);
    return TCL_ERROR;
  }
  nn_cond_params_p=(nn_cond_params*)malloc(sizeof(nn_cond_params));
  nn_cond_params_p->ids1=ids1;
  *change=1+temp;
  if (! parse_id_list(interp, argc-3, argv+3, &temp, &ids2) == TCL_OK ) {
    free(ids1);
    free(ids2);
    free(nn_cond_params_p);
    return TCL_ERROR;
  }
  *change+=temp;
  nn_cond_params_p->ids2=ids2;
  if ( argc < 6 || !ARG_IS_D(5,cutoff)) {
    Tcl_AppendResult(interp, usage_msg, (char *)NULL);
    free(ids1);
    free(ids2);
    free(nn_cond_params_p);
  return TCL_ERROR;
  }
  *change+=1;
  nn_cond_params_p->cutoff=cutoff;
  nn_cond_params_p->cutoff2=cutoff*cutoff;
  if ( argc < 7 || !ARG_IS_I(6,chain_length)) {
    Tcl_AppendResult(interp, usage_msg, (char *)NULL);
    free(ids1);
    free(ids2);
    free(nn_cond_params_p); 
    return TCL_ERROR;
  }
  *change+=1; 
  // allocate space for condition arguments
  nn_cond_params_p->conditions=(IntList*)malloc(sizeof(IntList));
  init_intlist(nn_cond_params_p->conditions); 
  alloc_intlist(nn_cond_params_p->conditions,ids1->n); 
  
  nn_cond_params_p->prev_partners=(IntList*)malloc(sizeof(IntList));
  init_intlist(nn_cond_params_p->prev_partners); 
  alloc_intlist(nn_cond_params_p->prev_partners,ids1->n);
  for (i=0; i< ids1->n; i++) {
    nn_cond_params_p->prev_partners->e[i]=-1;
    nn_cond_params_p->conditions->e[i]=-1;
  }
  nn_cond_params_p->prev_partners->n=i;
  nn_cond_params_p->prev_partners->n=i;
  
  nn_cond_params_p->chain_length=chain_length;
  obs->args=(void*)nn_cond_params_p;
  obs->n=1+2*(ids1->n); // number of ids from the 1st argument
  return TCL_OK;
}

int tclcommand_observable_particle_positions_conditional(Tcl_Interp* interp, int argc, char** argv, int* change, observable* obs) {
  IntList *ids1, *ids2, *partners, *conditions;
  char usage_msg[]="Usage: analyze correlation ... particle_positions_conditional id_list1 id_list2 cutoff";
  int temp, i;
  double cutoff;
  obs->fun = &observable_particle_positions_conditional;
  obs->fun_name=strdup("particle_positions_conditional");
  nn_cond_params* nn_cond_params_p=(nn_cond_params*) malloc(sizeof(nn_cond_params));
  ids1=(IntList*)malloc(sizeof(IntList));
  ids2=(IntList*)malloc(sizeof(IntList));
  if (! parse_id_list(interp, argc-1, argv+1, &temp, &ids1) == TCL_OK ) {
    free(ids1);
    free(ids2);
    free(nn_cond_params_p);
    return TCL_ERROR;
  }
  //printf("n: %d\n", ids1->n); fflush(stdout); //exit(1);
  nn_cond_params_p=(nn_cond_params*)malloc(sizeof(nn_cond_params));
  nn_cond_params_p->ids1=ids1;
  *change=1+temp;
  if (! parse_id_list(interp, argc-3, argv+3, &temp, &ids2) == TCL_OK ) {
    free(ids1);
    free(ids2);
    free(nn_cond_params_p);
    return TCL_ERROR;
  }
  *change+=temp;
  nn_cond_params_p->ids2=ids2;
  if ( argc < 6 || !ARG_IS_D(5,cutoff)) {
    Tcl_AppendResult(interp, usage_msg, (char *)NULL);
    free(ids1);
    free(ids2);
    free(nn_cond_params_p);
  return TCL_ERROR;
  }
  *change+=1;
  nn_cond_params_p->cutoff=cutoff;
  nn_cond_params_p->cutoff2=cutoff*cutoff;
  // allocate space for condition arguments
  nn_cond_params_p->conditions=(IntList*)malloc(sizeof(IntList));
  init_intlist(nn_cond_params_p->conditions); 
  alloc_intlist(nn_cond_params_p->conditions,ids1->n); 
  
  nn_cond_params_p->prev_partners=(IntList*)malloc(sizeof(IntList));
  init_intlist(nn_cond_params_p->prev_partners); 
  alloc_intlist(nn_cond_params_p->prev_partners,ids1->n);
  for (i=0; i< ids1->n; i++) {
    nn_cond_params_p->prev_partners->e[i]=-1;
    nn_cond_params_p->conditions->e[i]=-1;
  }
  nn_cond_params_p->prev_partners->n=i;
  nn_cond_params_p->prev_partners->n=i;
  
  obs->args=(void*)nn_cond_params_p;
  obs->n=4*(ids1->n); // 1 condition + 3 dimensions per ID
  return TCL_OK;
}

//  if (ARG0_IS_S("textfile")) {
//    // We still can only handle full files
//    if ( argc>1 ) {
//      fds= malloc(sizeof(file_data_source));
//      error = file_data_source_init(fds, argv[1], 0);
//      *change=2;
//      if (!error) {
//        *A_args=(void*) fds;
//        *dim_A = fds->n_columns;
//        *A_fun = (void*)&file_data_source_readline;
//        return TCL_OK;
//      } else {
//        Tcl_AppendResult(interp, "Error reading file ", argv[1] ,"\n", (char *)NULL);
//        Tcl_AppendResult(interp, file_data_source_init_errors[error] ,"\n", (char *)NULL);
//        return TCL_ERROR;
//      }
//    } else {
//      Tcl_AppendResult(interp, "Error in parse_observable textfile: no filename given" , (char *)NULL);
//      return TCL_ERROR;
//    }
//    return TCL_OK;
//  } 
//  if (ARG0_IS_S("tclinput")) {
//    if (argc>1 && ARG1_IS_I(temp)) {
//      *dim_A = temp;
//      *A_fun = &tcl_input;
//      *A_args=malloc(sizeof(tcl_input_data));
//      *change =2;
//      return TCL_OK;
//    } else {
//      Tcl_AppendResult(interp, "\nError in parse_observable tclinfo. You must pass the dimension of the observable." , (char *)NULL);
//      return TCL_ERROR;
//    }
//  }else {
//    return observable_usage(interp);
//  }
//  return 0 ;














#define REGISTER_OBSERVABLE(name,parser,id) \
  if (ARG_IS_S(2,#name)) { \
    observables[id]=malloc(sizeof(observable)); \
    if (parser(interp, argc-2, argv+2, &temp, observables[n_observables]) == TCL_OK) { \
      n_observables++; \
      argc-=1+temp; \
      argv+=1+temp; \
      sprintf(buffer,"%d",id); \
      Tcl_AppendResult(interp,buffer,(char *)NULL);\
      return TCL_OK; \
    } else { \
      free(observables[n_observables]);\
      Tcl_AppendResult(interp, "\nError parsing observable ", #name, "\n", (char *)NULL); \
      return TCL_ERROR; \
    } \
  } 



int tclcommand_observable(ClientData data, Tcl_Interp *interp, int argc, char **argv){
//  file_data_source* fds;
  char buffer[TCL_INTEGER_SPACE];
  int n;
  int id;
  int temp;
  int no;
  

  if (argc<2) {
    Tcl_AppendResult(interp, "Usage!!!\n", (char *)NULL);
    return TCL_ERROR;
  }

  if (argc > 1 && ARG_IS_S(1, "n_observables")) {
	  sprintf(buffer, "%d", n_observables);
    Tcl_AppendResult(interp, buffer, (char *)NULL );
    return TCL_OK;
  }

  
//  if (argc > 1 && ARG1_IS_I(no)) {
// }
  if (argc > 2 && ARG1_IS_S("new") ) {

    // find the next free observable id
    for (id=0;id<n_observables;id++) {
      if ( observables+id == 0 ) 
	break;
    }
    if (id==n_observables) 
      observables=(observable**) realloc(observables, (n_observables+1)*sizeof(observable*)); 

    REGISTER_OBSERVABLE(particle_velocities, tclcommand_observable_particle_velocities,id);
    REGISTER_OBSERVABLE(com_velocity, tclcommand_observable_com_velocity,id);
    REGISTER_OBSERVABLE(com_position, tclcommand_observable_com_position,id);
    REGISTER_OBSERVABLE(particle_positions, tclcommand_observable_particle_positions,id);
    REGISTER_OBSERVABLE(stress_tensor, tclcommand_observable_stress_tensor,id);
    REGISTER_OBSERVABLE(stress_tensor_acf_obs, tclcommand_observable_stress_tensor_acf_obs,id);
    REGISTER_OBSERVABLE(particle_currents, tclcommand_observable_particle_currents,id);
    REGISTER_OBSERVABLE(currents, tclcommand_observable_currents,id);
    REGISTER_OBSERVABLE(dipole_moment, tclcommand_observable_dipole_moment,id);
    REGISTER_OBSERVABLE(structure_factor, tclcommand_observable_structure_factor,id);
    REGISTER_OBSERVABLE(interacts_with, tclcommand_observable_interacts_with,id);
    REGISTER_OBSERVABLE(nearest_neighbour_conditional, tclcommand_observable_nearest_neighbour_conditional,id);
   REGISTER_OBSERVABLE(particle_positions_conditional, tclcommand_observable_particle_positions_conditional,id);
  //  REGISTER_OBSERVABLE(obs_nothing, tclcommand_observable_obs_nothing,id);
  //  REGISTER_OBSERVABLE(flux_profile, tclcommand_observable_flux_profile,id);
    REGISTER_OBSERVABLE(density_profile, tclcommand_observable_density_profile,id);
    REGISTER_OBSERVABLE(lb_velocity_profile, tclcommand_observable_lb_velocity_profile,id);
    REGISTER_OBSERVABLE(radial_density_profile, tclcommand_observable_radial_density_profile,id);
    REGISTER_OBSERVABLE(radial_flux_density_profile, tclcommand_observable_radial_flux_density_profile,id);
    REGISTER_OBSERVABLE(flux_density_profile, tclcommand_observable_flux_density_profile,id);
    REGISTER_OBSERVABLE(lb_radial_velocity_profile, tclcommand_observable_lb_radial_velocity_profile,id);
    Tcl_AppendResult(interp, "Unknown observable ", argv[2] ,"\n", (char *)NULL);
    return TCL_ERROR;
  }
  
  if (ARG1_IS_I(n)) {
    if (n>=n_observables || observables+n == NULL ) {
      sprintf(buffer,"%d \n",n);
      Tcl_AppendResult(interp, "Observable with id ", buffer, (char *)NULL);
      Tcl_AppendResult(interp, "is not defined\n", (char *)NULL);
      return TCL_ERROR;
    }
    if (argc > 2 && ARG_IS_S(2,"print")) {
      return tclcommand_observable_print(interp, argc-3, argv+3, &temp, observables[n]);
    }
    if (argc > 2 && ARG_IS_S(2,"print_parameters")) {
      return tclcommand_observable_print_parameters(interp, argc-3, argv+3, &temp, observables[n]);
    }
  }
  Tcl_AppendResult(interp, "Unknown observable ", argv[1] ,"\n", (char *)NULL);
  return TCL_ERROR;
}


int parse_id_list(Tcl_Interp* interp, int argc, char** argv, int* change, IntList** ids ) {
  int i,ret;
//  char** temp_argv; int temp_argc;
//  int temp;
  IntList* input=malloc(sizeof(IntList));
  IntList* output=malloc(sizeof(IntList));
  init_intlist(input);
  alloc_intlist(input,1);
  init_intlist(output);
  alloc_intlist(output,1);
  if (argc<1) {
    Tcl_AppendResult(interp, "Error parsing id list\n", (char *)NULL);
    return TCL_ERROR;
  }


  if (ARG0_IS_S("id")) {
    if (!parse_int_list(interp, argv[1],input)) {
      Tcl_AppendResult(interp, "Error parsing id list\n", (char *)NULL);
      return TCL_ERROR;
    } 
    *ids=input;
    for (i=0; i<input->n; i++) {
      if (input->e[i] >= n_total_particles) {
        Tcl_AppendResult(interp, "Error parsing ID list. Given particle ID exceeds the number of existing particles\n", (char *)NULL);
        return TCL_ERROR;
      }
    }
    *change=2;
    return TCL_OK;

  } else if ( ARG0_IS_S("types") ) {
    if (!parse_int_list(interp, argv[1],input)) {
      Tcl_AppendResult(interp, "Error parsing types list\n", (char *)NULL);
      return TCL_ERROR;
    } 
    if( (ret=convert_types_to_ids(input,output))<=0){
        Tcl_AppendResult(interp, "Error parsing types list. No particle with given types found.\n", (char *)NULL);
        return TCL_ERROR;
    } else { 
      *ids=output;
    }
    *change=2;
    return TCL_OK;
  } else if ( ARG0_IS_S("all") ) {
    if( (ret=convert_types_to_ids(NULL,output))<=0){
        Tcl_AppendResult(interp, "Error parsing keyword all. No particle found.\n", (char *)NULL);
        return TCL_ERROR;
    } else { 
      *ids=output;
    }
    *change=1;
    return TCL_OK;
  }

  Tcl_AppendResult(interp, "unknown keyword given to observable: ", argv[0] , (char *)NULL);
  return TCL_ERROR;
}


static int convert_types_to_ids(IntList * type_list, IntList * id_list){ 
      int i,j,n_ids=0,flag;
      sortPartCfg();
      for ( i = 0; i<n_total_particles; i++ ) {
         if(type_list==NULL) { 
		/* in this case we select all particles */
               flag=1 ;
         } else {  
                flag=0;
                for ( j = 0; j<type_list->n ; j++ ) {
                    if(partCfg[i].p.type == type_list->e[j])  flag=1;
	        }
         }
	 if(flag==1){
              realloc_intlist(id_list, id_list->n=n_ids+1);
	      id_list->e[n_ids] = i;
	      n_ids++;
	 }
      }
      return n_ids;
}

int observable_particle_velocities(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    A[3*i + 0] = partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

#ifdef ELECTROSTATICS
int observable_particle_currents(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    A[3*i + 0] = charge * partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = charge * partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}
int observable_currents(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].m.v[0]/time_step;
    j[1] += charge * partCfg[ids->e[i]].m.v[1]/time_step;
    j[2] += charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}
int observable_dipole_moment(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_total_particles)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].r.p[0];
    j[1] += charge * partCfg[ids->e[i]].r.p[1];
    j[2] += charge * partCfg[ids->e[i]].r.p[2];
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}
#endif

int observable_com_velocity(void* idlist, double* A, unsigned int n_A) {
/* TODO: this does not work with MASS ... */
  unsigned int i;
  double v_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    v_com[0] += partCfg[ids->e[i]].m.v[0]/time_step;
    v_com[1] += partCfg[ids->e[i]].m.v[1]/time_step;
    v_com[2] += partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=v_com[0]/ids->n;
  A[1]=v_com[1]/ids->n;
  A[2]=v_com[2]/ids->n;
  return 0;
}

int observable_blocked_com_velocity(void* idlist, double* A, unsigned int n_A) {
/* TODO: this does not work with MASS ... */
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  n_blocks=n_A/3; 
  blocksize=ids->n/n_blocks;
  for ( block = 0; block < n_blocks; block++ ) {
    for ( i = 0; i < blocksize; i++ ) {
      id = ids->e[block*blocksize+i];
      if (ids->e[i] >= n_total_particles)
        return 1;
      A[3*block+0] +=  partCfg[id].m.v[0]/time_step;
      A[3*block+1] +=  partCfg[id].m.v[1]/time_step;
      A[3*block+2] +=  partCfg[id].m.v[2]/time_step;
    }
    A[3*block+0] /=  blocksize;
    A[3*block+1] /=  blocksize;
    A[3*block+2] /=  blocksize;
  }
  return 0;
}

int observable_blocked_com_position(void* idlist, double* A, unsigned int n_A) {
/* TODO: this does not work with MASS ... */
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  n_blocks=n_A/3; 
  blocksize=ids->n/n_blocks;
  for ( block = 0; block < n_blocks; block++ ) {
    for ( i = 0; i < blocksize; i++ ) {
      id = ids->e[block*blocksize+i];
      if (ids->e[i] >= n_total_particles)
        return 1;
      A[3*block+0] +=  partCfg[id].r.p[0];
      A[3*block+1] +=  partCfg[id].r.p[1];
      A[3*block+2] +=  partCfg[id].r.p[2];
    }
    A[3*block+0] /=  blocksize;
    A[3*block+1] /=  blocksize;
    A[3*block+2] /=  blocksize;
  }
  return 0;
}

int observable_com_position(void* idlist, double* A, unsigned int n_A) {
/* TODO: this does not work with MASS ... */
  unsigned int i;
  double p_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
    p_com[0] += partCfg[ids->e[i]].r.p[0];
    p_com[1] += partCfg[ids->e[i]].r.p[1];
    p_com[2] += partCfg[ids->e[i]].r.p[2];
  }
  A[0]=p_com[0]/ids->n;
  A[1]=p_com[1]/ids->n;
  A[2]=p_com[2]/ids->n;
  return 0;
}

int observable_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int binx, biny, binz;
  double ppos[3];
  int img[3];
  IntList* ids;
  profile_data* pdata;
  sortPartCfg();
  pdata=(profile_data*) pdata_;
  ids=pdata->id_list;
  double bin_volume=(pdata->maxx-pdata->minx)*(pdata->maxy-pdata->miny)*(pdata->maxz-pdata->minz)/pdata->xbins/pdata->ybins/pdata->zbins;
    
  for ( i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    binx= (int) floor( pdata->xbins*  (ppos[0]-pdata->minx)/(pdata->maxx-pdata->minx));
    biny= (int) floor( pdata->ybins*  (ppos[1]-pdata->miny)/(pdata->maxy-pdata->miny));
    binz= (int) floor( pdata->zbins*  (ppos[2]-pdata->minz)/(pdata->maxz-pdata->minz));
    if (binx>=0 && binx < pdata->xbins && biny>=0 && biny < pdata->ybins && binz>=0 && binz < pdata->zbins) {
      A[binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz] += 1./bin_volume;
    } 
  }
  return 0;
}

#ifdef LB
int observable_lb_velocity_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i, j, k;
  unsigned int maxi, maxj, maxk;
  double xoffset, yoffset, zoffset;
  double x_incr, y_incr, z_incr;
  double p[3], v[3];
  profile_data* pdata;
  pdata=(profile_data*) pdata_;
  int linear_index;

    
  for ( i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  double normalization_factor = 1.;
  if ( pdata->xbins == 1 ) {
    maxi = (int) floor(box_l[0]/lbpar.agrid);
    normalization_factor/=maxi;
    xoffset=0;
    x_incr=lbpar.agrid;
  } else {
    maxi = pdata->xbins;
    xoffset=pdata->minx;
    x_incr=(pdata->maxx-pdata->minx)/(pdata->xbins-1);
  }
  if ( pdata->ybins == 1 ) {
    maxj = (int) floor(box_l[1]/lbpar.agrid);
    normalization_factor/=maxj;
    yoffset=0;
    y_incr=lbpar.agrid;
  } else {
    maxj = pdata->ybins;
    yoffset=pdata->miny;
    y_incr=(pdata->maxy-pdata->miny)/(pdata->ybins-1);
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) floor(box_l[2]/lbpar.agrid);
    normalization_factor/=maxk;
    zoffset=0;
    z_incr=lbpar.agrid;
  } else {
    maxk = pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }

  for ( i = 0; i < maxi; i++ ) {
    for ( j = 0; j < maxj; j++ ) {
      for ( k = 0; k < maxk; k++ ) {
        p[0]=xoffset + i*x_incr;
        p[1]=yoffset + j*y_incr;
        p[2]=zoffset + k*z_incr;
        if (lb_lbfluid_get_interpolated_velocity(p, v)!=0)
          return 1;
        linear_index = 0;
        if (pdata->xbins > 1)
          linear_index += i*pdata->ybins*pdata->zbins;
        if (pdata->ybins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;

        A[3*linear_index+0]+=v[0];
        A[3*linear_index+1]+=v[1];
        A[3*linear_index+2]+=v[2];
      }
    }
  }
  
  for ( i = 0; i<n_A; i++ ) {
    A[i]*=normalization_factor;
  }

  
  return 0;
}
#endif


#ifdef LB
int observable_lb_radial_velocity_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i, j, k;
  unsigned int maxi, maxj, maxk;
  double roffset, phioffset, zoffset;
  double r, phi, z;
  double r_incr, phi_incr, z_incr;
  double p[3], v[3];
  double v_r, v_phi, v_z;
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  int linear_index;

    
  for ( i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  double normalization_factor = 1.;
  if ( pdata->rbins == 1 ) {
    return 1;
  } else {
    maxi = pdata->rbins;
    roffset=pdata->minr;
    r_incr=(pdata->maxr-pdata->minr)/(pdata->rbins-1);
  }
  if ( pdata->phibins == 1 ) {
    maxj = (int)floor( 2*3.1415*pdata->maxr/lbpar.agrid ) ; 
    normalization_factor/=maxj;
    phioffset=0;
    phi_incr=2*3.1415/maxj;
  } else {
    maxj = pdata->phibins;
    phioffset=pdata->minphi;
    phi_incr=(pdata->maxphi-pdata->minphi)/(pdata->phibins-1);
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) floor(box_l[2]/lbpar.agrid);
    normalization_factor/=maxk;
    zoffset=-pdata->center[2];
    z_incr=lbpar.agrid;
  } else {
    maxk = pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }

  for ( i = 0; i < maxi; i++ ) {
    for ( j = 0; j < maxj; j++ ) {
      for ( k = 0; k < maxk; k++ ) {
        r = roffset + i*r_incr;
        phi = phioffset + j*phi_incr;
        z = zoffset + k*z_incr;
        p[0]=r*cos(phi)+pdata->center[0];
        p[1]=r*sin(phi)+pdata->center[1];
        p[2]=z+pdata->center[2];
        if (lb_lbfluid_get_interpolated_velocity(p, v)!=0)
          return 1;
        linear_index = 0;
        if (pdata->rbins > 1)
          linear_index += i*pdata->phibins*pdata->zbins;
        if (pdata->phibins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;
        if (r>0) {
          v_r = 1/r*((p[0]-pdata->center[0])*v[0] + (p[1]-pdata->center[1])*v[1]); 
          v_phi = 1/r/r*((p[0]-pdata->center[0])*v[1]-(p[1]-pdata->center[1])*v[0]);
        } else {
          v_r = 0;
          v_phi = 0;
        }
        v_z = v[2];

        A[3*linear_index+0]+=v_r;
        A[3*linear_index+1]+=v_phi;
        A[3*linear_index+2]+=v_z;
      }
    }
  }
  
  for ( i = 0; i<n_A; i++ ) {
    A[i]*=normalization_factor;
  }

  
  return 0;
}
#endif
void transform_to_cylinder_coordinates(double x, double y, double z_, double* r, double* phi, double* z) {
  *z =  z_;
  *r =  sqrt(x*x+y*y);
  *phi = atan2(y,x);
}

int observable_radial_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int binr, binphi, binz;
  double ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  sortPartCfg();
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
    
  for ( i = 0; i< n_A; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    transform_to_cylinder_coordinates(ppos[0]-pdata->center[0], ppos[1]-pdata->center[1], ppos[2]-pdata->center[2], &r, &phi, &z);
    //printf("%f %f %f %f %f %f\n", ppos[0], ppos[1], ppos[2], r*cos(phi)+pdata->center[0], r*sin(phi)+pdata->center[1], z+pdata->center[2]);
    binr  =(int)floor((r-pdata->minr)/rbinsize);
    binphi=(int)floor((phi-pdata->minphi)/phibinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);

    if (binr>=0 && binr < pdata->rbins && binphi>=0 && binphi < pdata->phibins && binz>=0 && binz < pdata->zbins) {
      bin_volume=PI*((pdata->minr+(binr+1)*rbinsize)*(pdata->minr+(binr+1)*rbinsize) - (pdata->minr+(binr)*rbinsize)*(pdata->minr+(binr)*rbinsize)) *zbinsize * phibinsize/2/PI;
      A[binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz] += 1./bin_volume;
    } 
  }
  return 0;
}

int observable_radial_flux_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int binr, binphi, binz;
  double ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  sortPartCfg();
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_r, v_phi, v_z;
    
  for ( i = 0; i< n_A; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    v[0]=partCfg[ids->e[i]].m.v[0]/time_step;
    v[1]=partCfg[ids->e[i]].m.v[1]/time_step;
    v[2]=partCfg[ids->e[i]].m.v[2]/time_step;
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    transform_to_cylinder_coordinates(ppos[0]-pdata->center[0], ppos[1]-pdata->center[1], ppos[2]-pdata->center[2], &r, &phi, &z);
    binr  =(int)floor((r-pdata->minr)/rbinsize);
    binphi=(int)floor((phi-pdata->minphi)/phibinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);

    if (binr>=0 && binr < pdata->rbins && binphi>=0 && binphi < pdata->phibins && binz>=0 && binz < pdata->zbins) {
      bin_volume=PI*((pdata->minr+(binr+1)*rbinsize)*(pdata->minr+(binr+1)*rbinsize) - (pdata->minr+(binr)*rbinsize)*(pdata->minr+(binr)*rbinsize)) *zbinsize * phibinsize/2/PI;
      v_r = 1/r*((ppos[0]-pdata->center[0])*v[0] + (ppos[1]-pdata->center[1])*v[1]); 
      v_phi = 1/r/r*((ppos[0]-pdata->center[0])*v[1]-(ppos[1]-pdata->center[1])*v[0]);
      v_z = v[2];
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 0] += v_r/bin_volume;
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 1] += v_phi/bin_volume;
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 2] += v_z/bin_volume;
    } 
  }
  return 0;
}

int observable_flux_density_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i;
  int binx, biny, binz;
  double ppos[3];
  double x, y, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  sortPartCfg();
  profile_data* pdata;
  pdata=(profile_data*) pdata_;
  ids=pdata->id_list;
  double xbinsize=(pdata->maxx - pdata->minx)/pdata->xbins;
  double ybinsize=(pdata->maxy - pdata->miny)/pdata->ybins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_x, v_y, v_z;
    
  for ( i = 0; i< n_A; i++ ) {
    A[i]=0;
  }
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
/* We use folded coordinates here */
    v[0]=partCfg[ids->e[i]].m.v[0]*time_step;
    v[1]=partCfg[ids->e[i]].m.v[1]*time_step;
    v[2]=partCfg[ids->e[i]].m.v[2]*time_step;
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    x=ppos[0];
    y=ppos[1];
    z=ppos[2];
    binx  =(int)floor((x-pdata->minx)/xbinsize);
    biny  =(int)floor((y-pdata->miny)/ybinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);


    if (binx>=0 && binx < pdata->xbins && biny>=0 && biny < pdata->ybins && binz>=0 && binz < pdata->zbins) {
      bin_volume=xbinsize*ybinsize*zbinsize;
      v_x=v[0];
      v_y=v[1];
      v_z=v[2];
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 0] += v_x/bin_volume;
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 1] += v_y/bin_volume;
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 2] += v_z/bin_volume;
    } 
  }
  return 0;
}
int observable_particle_positions(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  IntList* ids;
  sortPartCfg();
  ids=(IntList*) idlist;
  for ( i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_total_particles)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].r.p[0];
      A[3*i + 1] = partCfg[ids->e[i]].r.p[1];
      A[3*i + 2] = partCfg[ids->e[i]].r.p[2];
  }
  return 0;
}


int observable_particle_positions_conditional(void* params_p, double* A, unsigned int n_A) {
  nn_cond_params* params = (nn_cond_params*)params_p;
  unsigned int i, j;
  double pos1[3]; // position of the traced particle
  double pos2[3]; // position of the possible interaction partner
  double dist[3]; // minimum image vector between pos1 and pos2
  double dist2, dx, dy, dz; // position of the possible interaction partner
  double cut2 = params->cutoff2;
  int partner = -1; // by default ther is no partner
  IntList* ids1=params->ids1;
  IntList* ids2=params->ids2;
  sortPartCfg();

  for ( i = 0; i<ids1->n; i++ ) {
    // This checking should not be necessary if all functions are implemented properly
    if (ids1->e[i] >= n_total_particles)
      return 1; 
    pos1[0] = partCfg[ids1->e[i]].r.p[0]; 
    pos1[1] = partCfg[ids1->e[i]].r.p[1]; 
    pos1[2] = partCfg[ids1->e[i]].r.p[2];
    // first check for the last known partner
    if ( params->prev_partners->e[i] != -1) {
      j = params->prev_partners->e[i];
      pos2[0] = partCfg[j].r.p[0]; 
      pos2[1] = partCfg[j].r.p[1]; 
      pos2[2] = partCfg[j].r.p[2]; 
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if (dist2 < cut2) { 
	partner=ids2->e[j];
      }
    }
    // if not, then check all other candidates
    for ( j = 0; (j<ids2->n) && (partner == -1 ); j++ ) {
      if (ids2->e[j] >= n_total_particles)
        return 1; 
      pos2[0] = partCfg[ids2->e[j]].r.p[0]; 
      pos2[1] = partCfg[ids2->e[j]].r.p[1]; 
      pos2[2] = partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if (dist2 < cut2) { 
	partner=ids2->e[j];
      }
    }
    // check for a detachment event 
    if ( partner == -1 && params->prev_partners->e[i] != -1) {
      params->conditions->e[i] += 1;
      params->conditions->e[i] *= -1;
    }
    if ( partner > -1 && params->prev_partners->e[i] == -1) {
      params->conditions->e[i] *= -1;
    }
    A[4*i+0]=(double)params->conditions->e[i];
    A[4*i+1]=pos1[0];
    A[4*i+2]=pos1[1];
    A[4*i+3]=pos1[2];
    params->prev_partners->e[i]=partner;
  }
  return 0;
}

int observable_stress_tensor(void* params_p, double* A, unsigned int n_A) {
  unsigned int i;
  sortPartCfg();
  observable_compute_stress_tensor(1,A,n_A);
  return 0;
}


int observable_stress_tensor_acf_obs(void* params_p, double* A, unsigned int n_A) {
  unsigned int i;
  double stress_tensor[9];
  sortPartCfg();
  observable_compute_stress_tensor(1,&stress_tensor,9);
  A[0]=stress_tensor[1];
  A[1]=stress_tensor[5];
  A[2]=stress_tensor[6];
  A[3]=stress_tensor[0]-stress_tensor[4];
  A[4]=stress_tensor[0]-stress_tensor[8];
  A[5]=stress_tensor[4]-stress_tensor[8];
  return 0;
}

int observable_structure_factor(void* params_p, double* A, unsigned int n_A) {
  int i,j,k,l,p;
  int order, order2, n;
  double twoPI_L, C_sum, S_sum, qr;
  observable_sf_params* params;
  params = (observable_sf_params*)params_p;
  order = params->order;
  order2=order*order;
  twoPI_L = 2*PI/box_l[0];
  
  sortPartCfg();

    for(p=0; p<n_A; p++) {
       A[p]   = 0.0;
    }

    l=0;
    //printf("n_A: %d, dim_sf: %d\n",n_A, params.dim_sf); fflush(stdout);
    for(i=-order; i<=order; i++) {
      for(j=-order; j<=order; j++) {
        for(k=-order; k<=order; k++) {
	  n = i*i + j*j + k*k;
	  if ((n<=order2) && (n>=1)) {
	    C_sum = S_sum = 0.0;
            //printf("l: %d, n: %d %d %d\n",l,i,j,k); fflush(stdout);
	    for(p=0; p<n_total_particles; p++) {
	      qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
	      C_sum+= partCfg[p].p.scattering_length * cos(qr);
	      S_sum-= partCfg[p].p.scattering_length * sin(qr);
	    }
            A[l]   =C_sum;
            A[l+1] =S_sum;
            l=l+2;
	  }
	}
      }
    }
    //printf("finished calculating sf\n"); fflush(stdout);
    return 0;
}

int observable_interacts_with (void* params_p, double* A, unsigned int n_A) {
  iw_params *params=(iw_params*)params_p;
  IntList* ids1;
  IntList* ids2;
  int i,j;
//  double dx,dy,dz;
  double dist2;
  double cutoff2=params->cutoff*params->cutoff;
  double pos1[3], pos2[3], dist[3];
  ids1=params->ids1;
  ids2=params->ids2;
  sortPartCfg();
  for ( i = 0; i<ids1->n; i++ ) {
    if (ids1->e[i] >= n_total_particles)
      return 1;
    pos1[0]=partCfg[ids1->e[i]].r.p[0];
    pos1[1]=partCfg[ids1->e[i]].r.p[1];
    pos1[2]=partCfg[ids1->e[i]].r.p[2];
    for ( j = 0; j<ids2->n; j++ ) {
      if (ids2->e[j] >= n_total_particles)
        return 1;
      A[i] = 0;
      pos2[0]=partCfg[ids2->e[j]].r.p[0];
      pos2[1]=partCfg[ids2->e[j]].r.p[1];
      pos2[2]=partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if(dist2<cutoff2) {
        A[i] = 1;
	break;
	// interaction found for i, go for next
      }
    }
  }
  return 0;
}


int observable_nearest_neighbour_conditional (void* params_p, double* A, unsigned int n_A) {
  nn_cond_params *params=(nn_cond_params*)params_p;
  /* assume that each of the dim_A entries maps to the following sub-entries:
     A[i] -> int condition
     A[i] -> int position;
     */
  IntList* ids1;
  IntList* ids2;
  int i,j;
  double dist2;
  double mindist2;
  double cutoff2=params->cutoff2;
  int partner=-1;
  double pos1[3], pos2[3], dist[3];
  ids1=params->ids1;
  ids2=params->ids2;
  sortPartCfg();
  A[n_A-1]=params->chain_length; 
  if (n_A < 2*(ids1->n)+1) {
    // this should never happen
    printf("Error: n_A should be %d but is %d\n",2*(ids1->n)+1,n_A); fflush(stdout);
    return 1;
  }
  for ( i = 0; i<ids1->n; i++ ) {
    // just a safety check
    if (ids1->e[i] >= n_total_particles)
      return 1;
    pos1[0]=partCfg[ids1->e[i]].r.p[0];
    pos1[1]=partCfg[ids1->e[i]].r.p[1];
    pos1[2]=partCfg[ids1->e[i]].r.p[2];
    mindist2=cutoff2; partner=-1;
    for ( j = 0; j<ids2->n; j++ ) {
      if (ids2->e[j] >= n_total_particles)
        return 1;
      pos2[0]=partCfg[ids2->e[j]].r.p[0];
      pos2[1]=partCfg[ids2->e[j]].r.p[1];
      pos2[2]=partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if ( mindist2 > dist2 ) { 
	mindist2=dist2; 
	partner=ids2->e[j];
      }
    }
    // check for a detachment event 
    if ( partner == -1 && params->prev_partners->e[i] != -1) {
      params->conditions->e[i] += 1;
      params->conditions->e[i] *= -1;
    }
    if ( partner > -1 && params->prev_partners->e[i] == -1) {
      params->conditions->e[i] *= -1;
    }
    A[2*i+0]=params->conditions->e[i];
    A[2*i+1]=partner;
    params->prev_partners->e[i]=partner;
  }
  return 0;
}



//int file_data_source_init(file_data_source* self, char* filename, IntList* columns) {
//  int counter=1;
//  char* token;
//  if (filename==0)
//    return 1;
//  self->f = fopen(filename, "r");
//  if (! self->f )
//    return 2;
//  fgets(self->last_line, MAXLINELENGTH, self->f);
//  while (self->last_line && self->last_line[0] == 35) {
//    fgets(self->last_line, MAXLINELENGTH, self->f);
//  }
//  if (!self->last_line)
//    return 3;
//// Now lets count the tokens in the first line
//  token=strtok(self->last_line, " \t\n");
//  while (token) {
////    printf("reading token **%s**\n", token);
//    token=strtok(NULL, " \t\n");
//    counter++;
//  }
//  self->n_columns = counter;
//  rewind(self->f);
//  self->data_left=1;
////  printf("I found out that your file has %d columns\n", self->n_columns);
//  if (columns !=0)
//    /// Here we would like to check if we can handle the desired columns, but this has to be implemented!
//    return -1;
//  return 0;
//}

int file_data_source_readline(void* xargs, double* A, int dim_A) {
//  file_data_source* self = xargs;
//  int counter=0;
//  char* token;
//  char* temp;
//
//  temp=fgets(self->last_line, MAXLINELENGTH, self->f);
//  while (temp!= NULL && self->last_line && self->last_line[0] == 35) {
//    temp=fgets(self->last_line, MAXLINELENGTH, self->f);
//  }
//  if (!self->last_line || temp==NULL) {
////    printf("nothing left\n");
//    self->data_left=0;
//    return 3;
//  }
//  token=strtok(self->last_line, " \t\n");
//  while (token) {
////    printf("reading token: ");
//    A[counter]=atof(token);
////    printf("%f ", A[counter]);
//    token=strtok(NULL, " \t\n");
//    counter++;
//    if (counter >= dim_A) {
////      printf("urgs\n");
//      return 4;
//    }
//  }
////  printf("\n");
//  return 0;
//}
//
//int tcl_input(void* data, double* A, unsigned int n_A) {
//  tcl_input_data* input_data = (tcl_input_data*) data;
//  int i, tmp_argc;
//  const char  **tmp_argv;
//  Tcl_SplitList(input_data->interp, input_data->argv[0], &tmp_argc, &tmp_argv);
//  // function prototype from man page:
//  // int Tcl_SplitList(interp, list, argcPtr, argvPtr)
//  if (tmp_argc < n_A) {
//    Tcl_AppendResult(input_data->interp, "Not enough arguments passed to analyze correlation update", (char *)NULL);
//    return 1;
//  }
//  for (i = 0; i < n_A; i++) {
//    if (Tcl_GetDouble(input_data->interp, tmp_argv[i], &A[i]) != TCL_OK) {
//      Tcl_AppendResult(input_data->interp, "error parsing argument ", input_data->argv[i],"\n", (char *)NULL);
//      return 1;
//    }
//  }
  return 0;
}


int tclcommand_parse_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, profile_data** pdata_) {
  int temp;
  *change=0;
  profile_data* pdata=(profile_data*)malloc(sizeof(profile_data));
  *pdata_ = pdata;
  pdata->id_list=0;
  pdata->minx=0;
  pdata->maxx=box_l[0];
  pdata->xbins=1;
  pdata->miny=0;
  pdata->maxy=box_l[1];
  pdata->ybins=1;
  pdata->minz=0;
  pdata->maxz=box_l[2];
  pdata->zbins=1;
  while (argc>0) {
    if (ARG0_IS_S("id") || ARG0_IS_S("type") || ARG0_IS_S("all")) {
      if (!parse_id_list(interp, argc, argv, &temp, &pdata->id_list )==TCL_OK) {
        Tcl_AppendResult(interp, "Error reading profile: Error parsing particle id information\n" , (char *)NULL);
        return TCL_ERROR;
      } else {
        *change+=temp;
        argc-=temp;
        argv+=temp;
      }
    } else if ( ARG0_IS_S("minx")){
      if (argc>1 && ARG1_IS_D(pdata->minx)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else  if ( ARG0_IS_S("maxx") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxx)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else  if (ARG0_IS_S("xbins")) {
      if (argc>1 && ARG1_IS_I(pdata->xbins)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read nbins\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("miny")){
      if (argc>1 && ARG1_IS_D(pdata->miny)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else  if ( ARG0_IS_S("maxy") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxy)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else  if (ARG0_IS_S("ybins")) {
      if (argc>1 && ARG1_IS_I(pdata->ybins)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read nbins\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("minz")){
      if (argc>1 && ARG1_IS_D(pdata->minz)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else  if ( ARG0_IS_S("maxz") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxz)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else  if (ARG0_IS_S("zbins")) {
      if (argc>1 && ARG1_IS_I(pdata->zbins)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read nbins\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: understand argument ", argv[0], "\n" , (char *)NULL);
      return TCL_ERROR;
    }
  }
  
  temp=0;
  if (pdata->xbins <= 0 || pdata->ybins <=0 || pdata->zbins <= 0) {
    Tcl_AppendResult(interp, "Error in profile: the bin number in each direction must be >=1\n" , (char *)NULL);
    temp=1;
  }
  if (temp)
    return TCL_ERROR;
  else
    return TCL_OK;
}


int tclcommand_parse_radial_profile(Tcl_Interp* interp, int argc, char** argv, int* change, int* dim_A, radial_profile_data** pdata_) {
  int temp;
  *change=0;
  radial_profile_data* pdata=(radial_profile_data*)malloc(sizeof(radial_profile_data));
  *pdata_ = pdata;
  pdata->id_list=0;
  pdata->maxr=sqrt(box_l[0]*box_l[0]/4.+ box_l[1]*box_l[1]/4.);
  pdata->minr=0;
  pdata->maxphi=PI;
  pdata->minphi=-PI;
  pdata->minz=-box_l[2]/2.;
  pdata->maxz=+box_l[2]/2.;
  pdata->center[0]=box_l[0]/2.;pdata->center[1]=box_l[1]/2.;pdata->center[2]=box_l[2]/2.;
  pdata->rbins=1;
  pdata->zbins=1;
  pdata->phibins=1;
  pdata->axis[0]=0.;
  pdata->axis[1]=0.;
  pdata->axis[2]=1.;
  if (argc < 1) {
    Tcl_AppendResult(interp, "Usage radial_profile id $ids center $x $y $z maxr $r_max nbins $n\n" , (char *)NULL);
    return TCL_ERROR;
  }
  while (argc>0) {
    if (ARG0_IS_S("id") || ARG0_IS_S("type") || ARG0_IS_S("all")) {
      if (!parse_id_list(interp, argc, argv, &temp, &pdata->id_list )==TCL_OK) {
        Tcl_AppendResult(interp, "Error reading profile: Error parsing particle id information\n" , (char *)NULL);
        return TCL_ERROR;
      } else {
        *change+=temp;
        argc-=temp;
        argv+=temp;
      }
    } else if ( ARG0_IS_S("center")){
      if (argc>3 && ARG1_IS_D(pdata->center[0]) && ARG_IS_D(2,pdata->center[1]) && ARG_IS_D(3,pdata->center[2])) {
        argc-=4;
        argv+=4;
        *change+=4;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read center\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("axis")){
        Tcl_AppendResult(interp, "Using arbitrary axes does not work yet!\n" , (char *)NULL);
        return TCL_ERROR;
      if (argc>3 && ARG1_IS_D(pdata->axis[0]) && ARG_IS_D(2,pdata->axis[1]) && ARG_IS_D(3,pdata->axis[2])) {
        argc-=4;
        argv+=4;
        *change+=4;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read center\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if  ( ARG0_IS_S("maxr") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxr)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read maxr\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if  ( ARG0_IS_S("minr") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxr)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read maxr\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("minz")){
      if (argc>1 && ARG1_IS_D(pdata->minz)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("maxz") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxz)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("minphi")){
      if (argc>1 && ARG1_IS_D(pdata->minphi)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read minz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if ( ARG0_IS_S("maxphi") ) {
      if (argc>1 && ARG1_IS_D(pdata->maxphi)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in profile: could not read maxz\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if (ARG0_IS_S("rbins")) {
      if (argc>1 && ARG1_IS_I(pdata->rbins)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read rbins\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if (ARG0_IS_S("zbins")) {
      if (argc>1 && ARG1_IS_I(pdata->zbins)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read rbins\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else if (ARG0_IS_S("phibins")) {
      if (argc>1 && ARG1_IS_I(pdata->phibins)) {
        argc-=2;
        argv+=2;
        *change+=2;
      } else {
        Tcl_AppendResult(interp, "Error in radial_profile: could not read rbins\n" , (char *)NULL);
        return TCL_ERROR;
      } 
    } else {
      Tcl_AppendResult(interp, "Error in radial_profile: understand argument ", argv[0], "\n" , (char *)NULL);
      return TCL_ERROR;
    }
  }
  
  temp=0;
//  if (pdata->center[0]>1e90) {
//    Tcl_AppendResult(interp, "Error in radial_profile: center not specified\n" , (char *)NULL);
//    temp=1;
//  }
//  if (pdata->maxr>1e90) {
//    Tcl_AppendResult(interp, "Error in radial_profile: maxr not specified\n" , (char *)NULL);
//    temp=1;
//  }
//  if (pdata->rbins<1) {
//    Tcl_AppendResult(interp, "Error in radial_profile: rbins not specified\n" , (char *)NULL);
//    temp=1;
//  }
  if (temp)
    return TCL_ERROR;
  else
    return TCL_OK;
}


int sf_print_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "\nusage: structure_factor order delta_t tau_max  tau_lin", (char *)NULL);
  return TCL_ERROR;
}

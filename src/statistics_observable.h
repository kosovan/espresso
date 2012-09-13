/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _STATISTICS_OBSERVABLE_H
#define _STATISTICS_OBSERVABLE_H

#include "config.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct s_observable;

struct s_observable {
  char* obs_name;
  void* container;
  int n;
  int (*update)    ( struct s_observable* obs );
  int (*calculate) ( struct s_observable* obs );
  double* last_value;
  int last_update;
};

typedef struct s_observable observable;

extern observable** observables;
extern int n_observables; 

int observable_calculate(observable* self);
int observable_update(observable* self);

/* Here we have the particular observables listed */
int observable_calc_particle_velocities(observable* self_);
int observable_calc_com_velocity(observable* self); 
int observable_calc_blocked_com_velocity(observable* self); 
/** Obtain the particle positions.
 * TODO: Folded or unfolded?
 */ 
int observable_calc_particle_positions(observable* self);
int observable_calc_particle_forces(observable* self);
int observable_calc_com_force(observable* self);
int observable_calc_blocked_com_force(observable* self);
int observable_calc_stress_tensor(observable* self);
int observable_calc_stress_tensor_acf_obs(observable* self);
int observable_calc_com_position(observable* self);
int observable_calc_blocked_com_position(observable* self);

#ifdef ELECTROSTATICS
int observable_calc_particle_currents(observable* self);
int observable_calc_currents(observable* self);
int observable_calc_dipole_moment(observable* self);
#endif

/* Here go the one-time simple derived observable */
/* block average */
int observable_init_block_average(observable* self, observable* reference_observable, unsigned int blocksize, unsigned int stride);
int observable_calc_block_average(observable* self);
typedef struct {
  observable* reference_observable;
  unsigned int blocksize;
  unsigned int stride;
  unsigned int n_sweeps;
  unsigned int n_blocks;
} block_container;

/* block sum */
int observable_init_block_sum(observable* self, observable* reference_observable, unsigned int blocksize, unsigned int stride);
int observable_calc_block_sum(observable* self);



/* average */
int observable_update_average(observable* self);
int observable_reset_average(observable* self);
typedef struct {
  observable* reference_observable;
  unsigned int n_sweeps;
} observable_average_container;

/* variance */
int observable_init_variance(observable* self, observable* reference_observable);
int observable_update_variance(observable* self);
int observable_calc_variance(observable* self);
int observable_reset_variance(observable* self);
typedef struct {
  observable* reference_observable;
  unsigned int n_sweeps;
  double* sum;
  double* sum_squares;
} observable_variance_container;

/* stddev, also uses variance_container */
int observable_init_stddev(observable* self, observable* reference_observable);
int observable_update_stddev(observable* self);
int observable_calc_stddev(observable* self);
int observable_reset_stddev(observable* self);

/** Calculate structure factor from positions and scattering length */
int observable_calc_structure_factor(observable* self);
typedef struct {
// FIXME finish the implementation of scattering length
  IntList* id_list;
  DoubleList *scattering_length; // Scattering lengths of particles
  int order;
  int dim_sf; // number of q vectors
  int *q_vals; // values of q vectors
  double *q_density; // number of q vectors per bin
  // entries for spherical averaging
} observable_sf_params;

/** See if particles from idList1 interact with any of the particles in idList2 
input parameters are passed via struct iw_params
*/
int observable_calc_interacts_with(observable* self);
typedef struct {
  double cutoff;
  IntList *ids1;
  IntList *ids2;
} iw_params;

/** If a particle changes its state from "bound" to "unbound", record the lag time between 
binding and unbiding.
If the state is bound, it positive, otherwise it is zero.
Mask of (-1) should be used to switch behaviour to measuring
unbound lifetimes.
*/
int observable_calc_interaction_lifetimes(observable* self);
int observable_update_interaction_lifetimes(observable* self);
typedef struct {
  double cutoff;
  IntList *ids1;
  IntList *ids2;
  IntList *old_states;  // tracker for states in the past
  DoubleList *times;  // times when the old_states have become positive
  int max_n; // last entry where we have written
  int mask; // 1 = binding, -1 = unbinding
} lft_params;

/** Do nothing */
int observable_calc_obs_nothing (observable* self);

int observable_calc_flux_density_profile(observable* self);
typedef struct { 
  IntList* id_list;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
  int xbins;
  int ybins;
  int zbins;
} profile_data;

int observable_calc_density_profile(observable* self);

int observable_calc_lb_velocity_profile(observable* self);

int observable_calc_radial_density_profile(observable* self);
int observable_calc_radial_flux_density_profile(observable* self);
int observable_calc_lb_radial_velocity_profile(observable* self);
typedef struct {
  IntList* id_list;
  double minr;
  double maxr;
  double minphi;
  double maxphi;
  double minz;
  double maxz;
  double center[3];
  double axis[3];
  int phibins;
  int rbins;
  int zbins;
} radial_profile_data;



#endif

include "myconfig.pxi"

from libcpp cimport bool

cdef extern from "reaction_ensemble.hpp":
    ctypedef struct single_reaction:
        int* educt_types
        int len_educt_types
        int* educt_coefficients
        int* product_types
        int len_product_types
        int* product_coefficients
        double equilibrium_constant
        int nu_bar


    ctypedef struct reaction_system:
        int nr_single_reactions
        single_reaction** reactions
        int* type_index
        int nr_different_types
        double* charges_of_types
        int water_type
        double standard_pressure_in_simulation_units
        double given_length_in_SI_units
        double given_length_in_simulation_units
        double temperature_reaction_ensemble
        double exclusion_radius

    cdef extern reaction_system current_reaction_system
    int do_reaction()
    int find_index_of_type(int type)
    int calculate_nu_bar(int* educt_coefficients, int len_educt_types,  int* product_coefficients, int len_product_types)
    int update_type_index(int* educt_types, int len_educt_types , int* product_types, int len_product_types)




#///////////////////////////////////////////// Wang-Landau algorithm

    ctypedef struct collective_variable:
        double CV_minimum
        double CV_maximum
        double delta_CV
        double (*determine_current_state_in_collective_variable_with_index) (int)
        int* corresponding_acid_types
        int nr_corresponding_acid_types
        int associated_type
        char* energy_boundaries_filename

    ctypedef struct wang_landau_system:
        int* histogram
        int len_histogram
        double* wang_landau_potential
        int nr_collective_variables
        collective_variable** collective_variables
        int* nr_subindices_of_collective_variable
        double wang_landau_parameter
        double initial_wang_landau_parameter
        int int_fill_value
        double double_fill_value
        int number_of_monte_carlo_moves_between_check_of_convergence
        double final_wang_landau_parameter
        int used_bins
        int monte_carlo_trial_moves
        int wang_landau_steps
        char* output_filename
        double* minimum_energies_at_flat_index
        double* maximum_energies_at_flat_index
        bool do_energy_reweighting
        int counter_ion_type
        int polymer_start_id
        int polymer_end_id
        bool fix_polymer
        bool do_not_sample_reaction_partition_function
        bool use_hybrid_monte_carlo

    cdef extern wang_landau_system current_wang_landau_system
    int initialize_wang_landau()
    int do_reaction_wang_landau()
    void free_wang_landau()
    int update_maximum_and_minimum_energies_at_current_state()
    void write_out_preliminary_energy_run_results(char* filename)


    int write_wang_landau_checkpoint(char* identifier)
    int load_wang_landau_checkpoint(char* identifier)


cdef class ReactionEnsemble:
    cdef _params





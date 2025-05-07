// --------------------------------------------------------------------------
//
// Copyright (C) 2009 - 2016 by the adaflo authors
//
// This file is part of the adaflo library.
//
// The adaflo library is free software; you can use it, redistribute it,
// and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.  The full text of the
// license can be found in the file LICENSE at the top level of the adaflo
// distribution.
//
// --------------------------------------------------------------------------

#ifndef __adaflo_parameters_h
#define __adaflo_parameters_h

#include <deal.II/base/parameter_handler.h>

#include <adaflo/time_stepping.h>

#include <fstream>
#include <iostream>

namespace adaflo
{
  using namespace dealii;

  struct FlowParameters
  {
  public:
    // allows to defer construction to a derived class. Handle with care, as
    // parameters are not default-initialized
    FlowParameters();

    FlowParameters(const std::string &parameter_filename);

    void
    add_parameters(ParameterHandler &prm);

    void
    parse_parameters(const std::string parameter_filename, ParameterHandler &prm);

    void
    post();

    void
    check_for_file(const std::string &parameter_filename, ParameterHandler &prm) const;

    unsigned int dimension                  = 2;
    unsigned int global_refinements         = 1;
    unsigned int adaptive_refinements       = 0;
    bool         use_anisotropic_refinement = false;
    bool         use_simplex_mesh           = false;

    enum PhysicalType
    {
      incompressible,
      incompressible_stationary,
      stokes
    } physical_type;

    std::map<std::string, double> get_beta_formulation_convective_term_momentum_balance =
      {{"conservative", 1.0}, {"convective", 0.0}, {"skew-symmetric", 0.5}};

    double       beta_convective_term_momentum_balance = 0.5;
    unsigned int velocity_degree                       = 2;
    bool         augmented_taylor_hood                 = false;
    double       viscosity                             = 1;
    double       density                               = 1;
    double       damping                               = 0.0;
    double       tau_grad_div                          = 0;
    unsigned int max_nl_iteration                      = 10;
    double       tol_nl_iteration                      = 1e-6;
    enum Linearization
    {
      coupled_implicit_newton,
      coupled_implicit_picard,
      coupled_velocity_semi_implicit,
      coupled_velocity_explicit,
      projection
    } linearization;

    unsigned int max_lin_iteration = 500;
    double       tol_lin_iteration = 1e-3;
    bool         rel_lin_iteration = true;
    enum PreconditionVelocity
    {
      u_ilu,
      u_ilu_scalar,
      u_amg_linear,
      u_amg
    } precondition_velocity;
    enum PreconditionPressure
    {
      p_mass_diag,
      p_mass_ilu
    } precondition_pressure;
    unsigned int iterations_before_inner_solvers = 50;

    std::string  output_filename       = "";
    unsigned int output_verbosity      = 2;
    double       output_frequency      = 1;
    unsigned int print_solution_fields = 0;
    bool         output_wall_times     = false;
    bool         output_memory         = false;

    TimeSteppingParameters::Scheme time_step_scheme =
      TimeSteppingParameters::Scheme::bdf_2;
    double start_time           = 0;
    double end_time             = 1;
    double time_step_size_start = 1e-2;
    double time_stepping_cfl    = 0.8;
    double time_stepping_coef2  = 10;
    double time_step_tolerance  = 1e-2;
    double time_step_size_max   = 1.0;
    double time_step_size_min   = 0.1;

    // Two-phase specific parameters
    double density_diff   = 0;
    double viscosity_diff = 0;

    double surface_tension  = 1;
    double gravity          = 0;
    double epsilon          = 1;
    double diffusion_length = 0.1; // only useful in Cahn-Hilliard
    double contact_angle    = 0;   // only useful in Cahn-Hilliard

    bool pressure_constraint = true;

    unsigned int concentration_subdivisions     = 2;
    unsigned int curvature_correction           = 0;     // only for level set
    bool         interpolate_grad_onto_pressure = false; // only for level set
    bool         surface_tension_from_heaviside = true;  // only for level set
    bool         approximate_projections        = false; // only for level set
    bool         ch_do_newton                   = true;  // only for Cahn-Hilliard
    bool         do_iteration                   = false; // only for Cahn-Hilliard
    unsigned int n_reinit_steps                 = 2;     // only for level set
    unsigned int n_initial_reinit_steps         = 0;     // only for level set
    bool         convection_stabilization       = false;

    enum ConstitutiveType
    {
      newtonian_compressible_stokes_hypothesis,
      newtonian_incompressible,
      user_defined
    } constitutive_type;

  private:
    // for parameter parsing only
    double      two_phase_density                            = -1;
    double      two_phase_viscosity                          = -1;
    std::string physical_type_str                            = "incompressible";
    std::string constitutive_type_str                        = "newtonian incompressible";
    std::string formulation_convective_term_momentum_balance = "skew-symmetric";
    std::string linearization_str                            = "coupled implicit Newton";
    std::string uprec                                        = "amg linear";
    std::string pprec                                        = "ilu";
    std::string time_step_scheme_str                         = "bdf2";

    // temporary values for boolean variables
    struct BooleanVariables
    {
      int use_anisotropic_refinement     = false;
      int use_simplex_mesh               = false;
      int augmented_taylor_hood          = false;
      int rel_lin_iteration              = true;
      int pressure_constraint            = true;
      int output_wall_times              = false;
      int output_memory                  = false;
      int interpolate_grad_onto_pressure = false; 
      int surface_tension_from_heaviside = true;  
      int approximate_projections        = false; 
      int ch_do_newton                   = true;  
      int do_iteration                   = false; 
      int convection_stabilization       = false;
    } to_bool;
  };
} // namespace adaflo

#endif

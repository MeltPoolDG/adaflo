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

#include <deal.II/base/mpi.h>

#include <adaflo/parameters.h>


using namespace dealii;


adaflo::FlowParameters::FlowParameters()
  : dimension(numbers::invalid_unsigned_int)
{
  // do nothing
}



adaflo::FlowParameters::FlowParameters(const std::string &parameter_filename)
{
  ParameterHandler prm;
  add_parameters(prm);
  check_for_file(parameter_filename, prm);
  parse_parameters(parameter_filename, prm);
}



void
adaflo::FlowParameters::check_for_file(const std::string &parameter_filename,
                                       ParameterHandler & /*prm*/) const
{
  std::ifstream parameter_file(parameter_filename.c_str());

  if (!parameter_file)
    {
      parameter_file.close();

      std::ostringstream message;
      message << "Input parameter file <" << parameter_filename
              << "> not found. Please make sure the file exists!" << std::endl;

      AssertThrow(false, ExcMessage(message.str().c_str()));
    }
}

void
adaflo::FlowParameters::add_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Navier-Stokes");
  prm.add_parameter("dimension",
                    dimension,
                    "Defines the dimension of the problem. Not essential "
                    "to the Navier-Stokes class, but useful in many "
                    "applications.");
  prm.add_parameter("global refinements",
                    global_refinements,
                    "Defines the number of initial global refinements. "
                    "Not used in the Navier-Stokes class, but useful in "
                    "many applications.");
  prm.add_parameter("anisotropic refinement",
                    to_bool.use_anisotropic_refinement,
                    "defines whether the mesh should be refined "
                    "anisotropically in normal direction to the interface, "
                    "0 means no anisotropy");
  prm.add_parameter("simplex mesh",
                    to_bool.use_simplex_mesh,
                    "defines whether a simplex mesh has been provided, "
                    "0 means mesh with only quadrilaterals (2D) and hexahedra "
                    "(3D) has been provided");
  prm.add_parameter("adaptive refinements",
                    adaptive_refinements,
                    "Defines the number of adaptive refinements. Not used "
                    "in the Navier-Stokes class, but useful in many "
                    "applications.");
  prm.add_parameter("velocity degree",
                    velocity_degree,
                    "Sets the degree for velocity. Pressure degree is "
                    "velocity degree minus one. Currently implemented for "
                    "orders 2 to 6");
  prm.add_parameter("augmented Taylor-Hood elements",
                    to_bool.augmented_taylor_hood,
                    "Option to choose the pressure space FE_Q_DG0(p_degree) "
                    "instead of the standard space FE_Q(p_degree). This "
                    "adds a constant discontinuous part to the pressure "
                    "basis and gives element-wise divergence-free solutions. "
                    "It produces solutions that are in general better but "
                    "also a bit more expensive to compute.");
  prm.add_parameter("viscosity",
                    viscosity,
                    "Defines the fluid dynamic viscosity");
  prm.add_parameter("density", density, "Defines the fluid density", Patterns::Double());
  prm.add_parameter("damping", damping, "Defines the fluid damping", Patterns::Double());
  prm.add_parameter("physical type",
                    physical_type_str,
                    "Sets the type of equations, Navier-Stokes or Stokes. For "
                    "Navier-Stokes, one can choose between a stationary and a "
                    "time-dependent variant. The time-dependent Navier-Stokes "
                    "equations are the default.",
                    Patterns::Selection(
                      "incompressible|incompressible stationary|stokes"));
  prm.add_parameter(
    "constitutive type",
    constitutive_type_str,
    "Sets the type of constitutive equations. The incompressible Newtonian "
    "fluid assumption is the default case. Alternatively, a compressible "
    "Newtonian fluid formulation exploiting the Stokes hypothesis or a "
    "user defined type can be chosen.",
    Patterns::Selection(
      "newtonian incompressible|newtonian compressible stokes hypothesis|user defined"));
  prm.add_parameter("formulation convective term momentum balance",
                    formulation_convective_term_momentum_balance,
                    "Sets the formulation of the convective term in the "
                    "momentum balance of the Navier-Stokes equations, i.e. "
                    "∇·(u x u) =(u·∇)u + βu(∇·u). The parameter β will be "
                    "set to 1 for the conservative form, to 0 for the "
                    "convective form and to 0.5 for the skew-symmetric "
                    "form (default formulation).",
                    Patterns::Selection("skew-symmetric|convective|conservative"));

  prm.enter_subsection("Solver");
  prm.add_parameter("NL max iterations",
                    max_nl_iteration,
                    "Defines the maximum number of nonlinear Newton "
                    "iterations.");
  prm.add_parameter("NL tolerance",
                    tol_nl_iteration,
                    "Defines the tolerance in the residual l2 norm in "
                    "the nonlinear Newton iteration.");
  prm.add_parameter(
    "linearization scheme",
    linearization_str,
    "Sets how to treat the coupled nonlinear Navier-Stokes "
    "system. The 'coupled' variants solve for the full block "
    "system, whereas 'projection' applies a fractional-step "
    "pressure correction method with the solution of a pressure "
    "Poisson matrix. "
    "The nonlinear convective term can be treated by a"
    "full Newton iteration, a Picard iteration (fixed-point "
    "like), a semi-implicit approach with the same term as in "
    "the fixed-point like iteration but velocity extrapolated "
    "from the old time, and an approach where the complete "
    "convective term is treated explicitly. "
    "For the projection scheme, only the semi-implicit "
    "velocity treatment is implemented because iterating "
    "out the nonlinearity makes no sense.",
    Patterns::Selection(
      "coupled implicit Newton|coupled implicit Picard|coupled velocity semi-implicit|coupled velocity explicit|projection"));
  prm.add_parameter("tau grad div",
                    tau_grad_div,
                    "Adds the term (div(v), tau div(u))"
                    "to the weak form the momentum equation, which is "
                    "consistent with the Navier-Stokes equations but "
                    "penalizes the divergence more. This term is usually "
                    "referred to as grad-div stabilization. It simplifies "
                    "the solution of linear systems if tau is on the order "
                    "of unity but not too large (as the added term is "
                    "singular).");

  prm.add_parameter("lin max iterations",
                    max_lin_iteration,
                    "Maximum number of linear iterations");
  prm.add_parameter("lin tolerance",
                    tol_lin_iteration,
                    "Tolerance for the linear solver");
  prm.add_parameter("lin relative tolerance",
                    to_bool.rel_lin_iteration,
                    "Sets whether the residual for the linear solver "
                    "should be measured relative to the nonlinear residual "
                    "(recommended option).");
  prm.add_parameter("lin velocity preconditioner",
                    uprec,
                    "Sets the preconditioner for approximating the inverse "
                    "of the velocity matrix in the Schur complement "
                    "preconditioner. 'amg linear' uses a matrix based on "
                    "subdividing FE_Q into several linear elements to "
                    "create a matrix hierarchy. This might decrease "
                    "interpolation quality, but AMG is typically much better "
                    "for linears, so it is recommended for more complex "
                    "problems with relatively large time steps or large "
                    "viscosities, otherwise ILU. The method 'ilu scalar' "
                    "is a simplified ILU that only constructs the ILU for "
                    "one velocity block and applies the same operator to "
                    "all components. It is cheaper to apply but approximates "
                    "somewhat worse.",
                    Patterns::Selection("ilu|ilu scalar|amg linear|amg"));
  prm.add_parameter("lin pressure mass preconditioner",
                    pprec,
                    "Sets whether the pressure mass matrix in the Schur "
                    "complement should be represented by the diagonal only "
                    "or by an ILU based on the full pressure mass matrix.",
                    Patterns::Selection("ilu|diagonal"));
  prm.add_parameter("lin its before inner solvers",
                    iterations_before_inner_solvers,
                    "The linear solver comes in two flavors. A simple "
                    "solver which uses only AMG V-cycles or ILUs as "
                    "preconditioner components in the Schur complement, "
                    "or a stronger solver with inner iterations. The "
                    "variant with inner solves is less efficient when "
                    "only a few iterations are needed, but much more "
                    "robust and more efficient for many iterations. This "
                    "option sets how many linear iterations with the cheap "
                    "preconditioners should be made before the stronger "
                    "version with more iterations starts.");
  prm.leave_subsection();
  prm.leave_subsection();

  prm.enter_subsection("Output options");
  prm.add_parameter("output filename",
                    output_filename,
                    "Sets the base name for the file output.");
  prm.add_parameter("output verbosity",
                    output_verbosity,
                    "Sets the amount of information from the "
                    "Navier-Stokes solver that is printed to screen. "
                    "0 means no output at all, and larger numbers mean an "
                    "increasing amount of output (maximum value: 3). "
                    "A value of 3 not only includes solver iterations "
                    "but also details on solution time and some memory "
                    "statistics.");
  prm.add_parameter("output frequency",
                    output_frequency,
                    "defines at with time interface the solution "
                    "should be written to file (in supported routines)");
  prm.add_parameter("output vtk files",
                    print_solution_fields,
                    "defines whether to output vtk files with the "
                    "whole solution field or just collected point data");
  prm.add_parameter("output wall times",
                    to_bool.output_wall_times,
                    "Defines whether to output wall times. 0 means no output.");
  prm.add_parameter("output memory",
                    to_bool.output_memory,
                    "Defines whether to output memory. 0 means no output.");
  prm.leave_subsection();

  prm.enter_subsection("Two phase");
  prm.add_parameter("density",
                    two_phase_density,
                    "Density of fluid 1 (negative region of level set function). "
                    "If given a positive value, overwrites density in "
                    "Navier-Stokes subsection.");
  prm.add_parameter("density difference",
                    density_diff,
                    "absolute difference in density compared to fluid 1");
  prm.add_parameter("viscosity",
                    two_phase_viscosity,
                    "Dynamic viscosity of fluid 1 (negative region of level "
                    "set function). If given a positive value, overwrites "
                    "density in Navier-Stokes subsection.");
  prm.add_parameter("viscosity difference",
                    viscosity_diff,
                    "absolute difference in viscosity compared to fluid 1");

  prm.add_parameter("surface tension", surface_tension, "surface tension coefficient");
  prm.add_parameter("epsilon",
                    epsilon,
                    "Width of diffuse interface, relative to mesh size "
                    "for Level-Set method, but absolute for Cahn-Hilliard.");
  prm.add_parameter("gravity", gravity, "Gravity.");
  prm.add_parameter("diffusion length",
                    diffusion_length,
                    "Diffusion length scale in Cahn-Hilliard. Its square "
                    "equals the mobility and inverse Peclet number.");
  prm.add_parameter("contact angle",
                    contact_angle,
                    "defines the contact angle at solid interfaces, "
                    "at boundaries with indicator 0 or 2");
  prm.add_parameter("pressure constraint",
                    to_bool.pressure_constraint,
                    "Fixes value of pressure in one point to zero");
  prm.add_parameter("concentration subdivisions",
                    concentration_subdivisions,
                    "Number of subdivision of Q1 elements in smaller elements "
                    "to generate higher accuracy in level set/phase field");
  prm.add_parameter("curvature correction",
                    curvature_correction,
                    "if 1, extend the curvature to the value "
                    "at the interface in normal direction");
  prm.add_parameter("grad pressure compatible",
                    to_bool.interpolate_grad_onto_pressure,
                    "if 1, the gradient in the surface tension force "
                    "is interpolated from the pressure gradient",
                    Patterns::Selection("0|1"));
  prm.add_parameter("localize surface tension",
                    to_bool.surface_tension_from_heaviside,
                    "if 1, the surface tension is computed from a gradient "
                    "that is localized around the interface (from a "
                    "reconstructed distance function), otherwise it is "
                    "computed from the tanh profile (i.e., nonzero "
                    "everywhere)",
                    Patterns::Selection("0|1"));
  prm.add_parameter("approximate projections",
                    to_bool.approximate_projections,
                    "if 0, the normal and curvature in the level set method "
                    "are computed by proper projection (full mass matrix "
                    "and little diffusion), otherwise with diagonal mass "
                    "matrix and time-dependent diffusion",
                    Patterns::Selection("0|1"));
  prm.add_parameter("Cahn-Hilliard do Newton",
                    to_bool.ch_do_newton,
                    "Sets whether a Newton iteration should be done on the "
                    "Cahn-Hilliard equation (if on that model). If 0 is "
                    "selected, use a convexity splitting as proposed by "
                    "Eyre.",
                    Patterns::Selection("0|1"));
  prm.add_parameter("full nonlinear iteration",
                    do_iteration,
                    "iterates between Navier-Stokes and concentration "
                    "if enabled",
                    Patterns::Selection("0|1"));
  prm.add_parameter("number reinit steps",
                    n_reinit_steps,
                    "number of iterations in reinitialization");
  prm.add_parameter("number initial reinit steps",
                    n_initial_reinit_steps,
                    "reinitialization steps before starting the time "
                    "loop (for bad initial profiles)");
  prm.add_parameter("convection stabilization",
                    convection_stabilization,
                    "add stabilization terms to advection equation if "
                    "set to 1 (typically not necessary)",
                    Patterns::Selection("0|1"));
  prm.leave_subsection();


  prm.enter_subsection("Time stepping");
  prm.add_parameter("start time", start_time, "Sets the start time for the simulation");
  prm.add_parameter("end time", end_time, "Sets the final time for the simulation");
  prm.add_parameter("step size",
                    time_step_size_start,
                    "Sets the step size for time stepping. For non-uniform "
                    "time stepping, this sets the size of the first time "
                    "step.");
  prm.add_parameter("CFL number",
                    time_stepping_cfl,
                    "Limits the time step size in terms of a condition "
                    "dt <= CFL * dx / |u|, where u is a characteristic velocity. "
                    "For two-phase flow, we typically take the velocity of "
                    "the bubble");
  prm.add_parameter("CFL number capillary",
                    time_stepping_coef2,
                    "Limits the time step size in terms of a condition "
                    "dt <= CFL_cap * sqrt(rho/sigma) * dx^1.5, i.e., it "
                    "represents a capillarity time step limit.");
  prm.add_parameter("tolerance",
                    time_step_tolerance,
                    "Sets the tolerance for time step selection in "
                    "non-uniform time stepping strategies.");
  prm.add_parameter("max step size",
                    time_step_size_max,
                    "Defines the maximum time step size in non-uniform "
                    "strategies.");
  prm.add_parameter("min step size",
                    time_step_size_min,
                    "Defines the minimum time step size in non-uniform "
                    "strategies.");
  prm.add_parameter("scheme",
                    time_step_scheme_str,
                    "Sets the time stepping scheme. Allowed options are "
                    "explicit_euler, implicit_euler, crank_nicolson "
                    "fractional0, fractional1, new_variant, and bdf_2.",
                    Patterns::Selection("explicit_euler|implicit_euler|"
                                        "crank_nicolson|bdf_2"));
  prm.leave_subsection();
}



void
adaflo::FlowParameters::parse_parameters(const std::string parameter_file,
                                         ParameterHandler &prm)
{
  try
    {
      if (parameter_file.substr(parameter_file.find_last_of(".") + 1) == "json")
        {
          std::ifstream file;
          file.open(parameter_file);
          prm.parse_input_from_json(file, true/*skip_undefined*/);
        }
      else if (parameter_file.substr(parameter_file.find_last_of(".") + 1) == "prm")
        prm.parse_input(parameter_file);
      else
        AssertThrow(false, ExcMessage("Parameterhandler cannot handle current file"));
    }
  catch (...)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        prm.print_parameters(std::cout, ParameterHandler::PRM);
      AssertThrow(false, ExcMessage("Invalid input parameter file."));
    }

  this->post();
}

void
adaflo::FlowParameters::post()
{
  AssertThrow(velocity_degree > 1, ExcNotImplemented());
  damping =
    -damping; // sign(damping) = minus: damping; sign(damping) = plus: acceleration
  if (physical_type_str == "stokes")
    physical_type = stokes;
  else if (physical_type_str == "incompressible")
    physical_type = incompressible;
  else if (physical_type_str == "incompressible stationary")
    physical_type = incompressible_stationary;
  else
    Assert(false, ExcNotImplemented());

  if (physical_type == stokes)
    density = 0;

  if (constitutive_type_str == "newtonian incompressible")
    constitutive_type = newtonian_incompressible;
  else if (constitutive_type_str == "newtonian compressible stokes hypothesis")
    constitutive_type = newtonian_compressible_stokes_hypothesis;
  else if (constitutive_type_str == "user defined")
    constitutive_type = user_defined;
  else
    Assert(false, ExcNotImplemented());
  beta_convective_term_momentum_balance =
    get_beta_formulation_convective_term_momentum_balance
      [formulation_convective_term_momentum_balance];

  if (linearization_str == "coupled implicit Newton")
    linearization = coupled_implicit_newton;
  else if (linearization_str == "coupled implicit Picard")
    linearization = coupled_implicit_picard;
  else if (linearization_str == "coupled velocity semi-implicit")
    linearization = coupled_velocity_semi_implicit;
  else if (linearization_str == "coupled velocity explicit")
    linearization = coupled_velocity_explicit;
  else if (linearization_str == "projection")
    linearization = projection;
  else
    Assert(false,
           ExcMessage(("Linearization " + linearization_str + " not available").c_str()));

  if (physical_type == incompressible_stationary)
    Assert(linearization == coupled_implicit_newton,
           ExcMessage("Only coupled implicit Newton linearization available for "
                      "stationary equation"));

  AssertThrow(tau_grad_div >= 0., ExcMessage("Invalid parameter value"));

  if (uprec == "ilu")
    precondition_velocity = u_ilu;
  else if (uprec == "ilu scalar")
    precondition_velocity = u_ilu_scalar;
  else if (uprec == "amg linear")
    precondition_velocity = u_amg_linear;
  else if (uprec == "amg")
    precondition_velocity = u_amg;
  else
    AssertThrow(false, ExcMessage("Invalid name"));

  if (pprec == "ilu")
    precondition_pressure = p_mass_ilu;
  else
    precondition_pressure = p_mass_diag;

  Assert(output_verbosity <= 3, ExcInternalError());

  if (print_solution_fields > 2)
    print_solution_fields = 1;

  if (two_phase_density > 0)
    density = two_phase_density;
  if (physical_type == stokes)
    density = density_diff = 0;

  if (two_phase_viscosity > 0)
    viscosity = two_phase_viscosity;

  AssertThrow(diffusion_length > 0, ExcMessage("Diffusion length must be positive"));
  AssertThrow(epsilon > 0, ExcMessage("Diffusion length must be positive"));

  // no adaptive time stepping in case the start step was large
  if (time_step_size_min > time_step_size_start)
    time_step_size_max = time_step_size_min = time_step_size_start;

  const std::string& schem = time_step_scheme_str;
  if (schem == "implicit_euler")
    time_step_scheme = TimeSteppingParameters::Scheme::implicit_euler;
  else if (schem == "explicit_euler")
    time_step_scheme = TimeSteppingParameters::Scheme::explicit_euler;
  else if (schem == "crank_nicolson")
    time_step_scheme = TimeSteppingParameters::Scheme::crank_nicolson;
  else if (schem == "bdf_2")
    time_step_scheme = TimeSteppingParameters::Scheme::bdf_2;
  else
    // parameter handler should make sure that we
    // never end up here
    AssertThrow(false, ExcInternalError());

  // convert boolean variables
  use_anisotropic_refinement     = to_bool.use_anisotropic_refinement    ;        
  use_simplex_mesh               = to_bool.use_simplex_mesh              ;
  augmented_taylor_hood          = to_bool.augmented_taylor_hood         ;
  rel_lin_iteration              = to_bool.rel_lin_iteration             ;
  pressure_constraint            = to_bool.pressure_constraint           ;
  output_wall_times              = to_bool.output_wall_times             ;
  output_memory                  = to_bool.output_memory                 ;
  interpolate_grad_onto_pressure = to_bool.interpolate_grad_onto_pressure;
  surface_tension_from_heaviside = to_bool.surface_tension_from_heaviside;
  approximate_projections        = to_bool.approximate_projections       ;
  ch_do_newton                   = to_bool.ch_do_newton                  ;
  do_iteration                   = to_bool.do_iteration                  ;
  convection_stabilization       = to_bool.convection_stabilization      ;
}

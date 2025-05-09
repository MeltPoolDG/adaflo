// --------------------------------------------------------------------------
//
// Copyright (C) 2008 - 2016 by the adaflo authors
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


#ifndef __adaflo_level_set_okz_h
#define __adaflo_level_set_okz_h

#include <adaflo/diagonal_preconditioner.h>
#include <adaflo/level_set_base.h>
#include <adaflo/level_set_okz_advance_concentration.h>
#include <adaflo/level_set_okz_compute_curvature.h>
#include <adaflo/level_set_okz_compute_normal.h>
#include <adaflo/level_set_okz_reinitialization.h>


namespace adaflo
{
  using namespace dealii;


  // forward declarations
  class BlockMatrixExtension;
  class BlockILUExtension;
  namespace AssemblyData
  {
    struct Data;
  }



  template <int dim>
  class LevelSetOKZSolver : public LevelSetBaseAlgorithm<dim>
  {
  public:
    LevelSetOKZSolver(const FlowParameters &parameters,
                      Triangulation<dim>   &triangulation);
    virtual ~LevelSetOKZSolver()
    {}

    virtual void
    initialize_data_structures() override;

    virtual void
    transform_distance_function(
      LinearAlgebra::distributed::Vector<double> &vector) const override;

  private:
    // compute the force term and variable density/viscosity for the
    // Navier--Stokes equations
    virtual void
    compute_force() override;

    // advection step
    virtual void
    advance_concentration() override;

    // computes normal direction vector by projection of level set gradient (not
    // scaled to have norm 1!)
    virtual void
    compute_normal(const bool fast_computation) override;

    // computes curvature by projecting the divergence of the normal vector
    // (scaled to norm 1 now)
    virtual void
    compute_curvature(const bool diffuse_large_values = false) override;

    // performs reinitialization
    virtual void
    reinitialize(const unsigned int stab_steps,
                 const unsigned int diff_steps                      = 0,
                 const bool diffuse_cells_with_large_curvature_only = false) override;

    virtual void
    compute_heaviside() override;

    // matrix-free worker operations for various operations
    template <int ls_degree, int velocity_degree>
    void
    local_compute_force(const MatrixFree<dim, double>                    &data,
                        LinearAlgebra::distributed::Vector<double>       &dst,
                        const LinearAlgebra::distributed::Vector<double> &src,
                        const std::pair<unsigned int, unsigned int>      &cell_range);

    void
    local_projection_matrix(
      const MatrixFree<dim, double>                                    &data,
      std::shared_ptr<Threads::ThreadLocalStorage<AssemblyData::Data>> &scratch,
      const unsigned int &,
      const std::pair<unsigned int, unsigned int> &cell_range);

    template <int ls_degree>
    void
    local_projection_matrix(
      const MatrixFree<dim, double>                                    &data,
      std::shared_ptr<Threads::ThreadLocalStorage<AssemblyData::Data>> &scratch,
      const std::pair<unsigned int, unsigned int>                      &cell_range);

    bool                           first_reinit_step;
    double                         global_max_velocity;
    DiagonalPreconditioner<double> preconditioner;

    std::shared_ptr<BlockMatrixExtension> projection_matrix;
    std::shared_ptr<BlockILUExtension>    ilu_projection_matrix;

    std::unique_ptr<LevelSetOKZSolverComputeNormal<dim>>        normal_operator;
    std::unique_ptr<LevelSetOKZSolverComputeCurvature<dim>>     curvatur_operator;
    std::unique_ptr<LevelSetOKZSolverReinitialization<dim>>     reinit_operator;
    std::unique_ptr<LevelSetOKZSolverAdvanceConcentration<dim>> advection_operator;
  };
} // namespace adaflo


#endif

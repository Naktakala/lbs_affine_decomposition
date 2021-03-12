#ifndef AFFINE_DECOMPOSITIONER_H
#define AFFINE_DECOMPOSITIONER_H

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"

namespace LinearBoltzmann
{

class AffineDecompositioner : public LinearBoltzmann::Solver
{
public:
  enum class ProjectionMode
  {
    Galerkin         = 1,
    Petrov_Galerkin  = 2
  };
  struct
  {
    unsigned int num_modes=0;
    ProjectionMode projection_mode =
      ProjectionMode::Galerkin;
  }options;
private:
  std::vector<std::vector<double>> U;
  std::set<int> unique_material_ids;
public:
  AffineDecompositioner();

  void Initialize() override;

  void ReadPODModes(const std::string& file_base_name,
                    unsigned int groupset_num=0);

private:
  void ReadSinglePODModeBinary(const std::string& file_name,
                               const LBSGroupset& groupset,
                               std::vector<double>& U_i);

public:
  void WriteTotalInteractionOperator(const std::string& file_base_name,
                                     const LBSGroupset& groupset);
};

}

#endif //AFFINE_DECOMPOSITIONER_H
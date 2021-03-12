#include "affine_decompositioner.h"

//###################################################################
/**Overrides LBS initialization routine.*/
void LinearBoltzmann::AffineDecompositioner::Initialize()
{
  //============================================= Calling parent class initialize
  Solver::Initialize();

  //============================================= Develop list of unique
  //                                              material ids
  for (auto& cell : grid->local_cells)
    unique_material_ids.insert(cell.material_id);
}
#include "ChiLua/chi_lua.h"

#include "../affine_decompositioner.h"

#include "ChiPhysics/chi_physics.h"

//###################################################################
/**Creates an AffineDecompositioner object and pushes it to the
 * physics stack. Returns a handle.*/
int chiAffineDecompositionerCreate(lua_State* L)
{
  auto affine_decompositioner = new LinearBoltzmann::AffineDecompositioner;

  auto& physics = ChiPhysics::GetInstance();
  physics.solver_stack.push_back(affine_decompositioner);

  lua_pushnumber(L,physics.solver_stack.size()-1);
  return 1;
}
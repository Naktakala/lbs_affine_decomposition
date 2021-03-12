#include "ChiLua/chi_lua.h"

#include "../affine_decompositioner.h"

#include "ChiPhysics/chi_physics.h"
#include "chi_log.h"

//###################################################################
/**Initializes an AffineDecompositioner
 *
\param Handle int Handle to an existing AffineDecompositioner object.*/
int chiAffineDecompositionerInitialize(lua_State* L)
{
  auto& physics = ChiPhysics::GetInstance();
  auto& chi_log = ChiLog::GetInstance();

  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int handle = lua_tonumber(L,1);

  chi_physics::Solver* solver;
  try {
    solver = physics.solver_stack.at(handle);
  }
  catch (const std::out_of_range& oor)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Invalid handle to solver.";
    exit(EXIT_FAILURE);
  }

  auto aff_decomp = dynamic_cast<LinearBoltzmann::AffineDecompositioner*>(solver);
  if (not aff_decomp)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Invalid handle to solver. The handle does not point"
                         " to a solver that is of type "
                         "LinearBoltzmann::AffineDecompositioner.";
    exit(EXIT_FAILURE);
  }

  aff_decomp->Initialize();

  return 0;
}

//###################################################################
/**Initializes an AffineDecompositioner
 *
\param Handle int Handle to an existing AffineDecompositioner object.
\param Property int What property to set.
\param PropertyValue varying Depending on the property.

 No return.*/
int chiAffineDecompositionSetProperty(lua_State* L)
{
  auto& physics = ChiPhysics::GetInstance();
  auto& chi_log = ChiLog::GetInstance();

  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int handle = lua_tonumber(L,1);

  chi_physics::Solver* solver;
  try {
    solver = physics.solver_stack.at(handle);
  }
  catch (const std::out_of_range& oor)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Invalid handle to solver.";
    exit(EXIT_FAILURE);
  }

  auto aff_decomp = dynamic_cast<LinearBoltzmann::AffineDecompositioner*>(solver);
  if (not aff_decomp)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Invalid handle to solver. The handle does not point"
                         " to a solver that is of type "
                         "LinearBoltzmann::AffineDecompositioner.";
    exit(EXIT_FAILURE);
  }

  //============================================= Handle what property
  LuaCheckNilValue(__FUNCTION__, L, 2);
  const std::string property = lua_tostring(L,2);

  LuaCheckNilValue(__FUNCTION__, L, 3);

  if (property == "num_modes")
  {
    aff_decomp->options.num_modes = lua_tonumber(L,3);
  }
  else if (property == "projection_mode")
  {
    const std::string mode = lua_tostring(L,3);

    if (mode == "Galerkin")
      aff_decomp->options.projection_mode =
        LinearBoltzmann::AffineDecompositioner::ProjectionMode::Galerkin;
    else if (mode == "Petrov-Galerkin")
      aff_decomp->options.projection_mode =
        LinearBoltzmann::AffineDecompositioner::ProjectionMode::Petrov_Galerkin;
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Unsupported projection mode " << mode
        << " used for property projection_mode.";
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Unsupported property " << property << " supplied"
                         " for second argument.";
    exit(EXIT_FAILURE);
  }

  return 0;
}

//###################################################################
/**Initializes an AffineDecompositioner
 *
\param Handle int Handle to an existing AffineDecompositioner object.
\param FileBaseNum string Basename of the POD modes. A postfix of the mode number
                          and ".data" will be added to each base name.
\param Groupsetnum int The groupset number for which this applied

 No return.*/
int chiAffineDecompositionReadPODModes(lua_State* L)
{
  auto& physics = ChiPhysics::GetInstance();
  auto& chi_log = ChiLog::GetInstance();

  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int handle = lua_tonumber(L,1);

  chi_physics::Solver* solver;
  try {
    solver = physics.solver_stack.at(handle);
  }
  catch (const std::out_of_range& oor)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Invalid handle to solver.";
    exit(EXIT_FAILURE);
  }

  auto aff_decomp = dynamic_cast<LinearBoltzmann::AffineDecompositioner*>(solver);
  if (not aff_decomp)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": Invalid handle to solver. The handle does not point"
                         " to a solver that is of type "
                         "LinearBoltzmann::AffineDecompositioner.";
    exit(EXIT_FAILURE);
  }

  //============================================= File base name
  LuaCheckNilValue(__FUNCTION__, L, 2);
  const std::string file_base_name = lua_tostring(L,2);

  //============================================= Groupset num
  LuaCheckNilValue(__FUNCTION__, L, 3);
  unsigned int groupset_num = lua_tonumber(L,3);

  LBSGroupset* groupset;
  try{
    groupset = &aff_decomp->group_sets.at(groupset_num);
  }
  catch (const std::out_of_range& oor)
  {
    chi_log.Log(LOG_ALLERROR)
      << __FUNCTION__ << ": groupset number.";
    exit(EXIT_FAILURE);
  }

  aff_decomp->ReadPODModes(file_base_name,groupset_num);

  return 0;
}
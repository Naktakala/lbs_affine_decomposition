#include "affine_decompositioner.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "chi_log.h"

//###################################################################
/**Reads a number of POD modes.*/
void LinearBoltzmann::AffineDecompositioner::
  ReadPODModes(const std::string &file_base_name,
               unsigned int groupset_num/*=0*/)
{
  U.resize(options.num_modes);
  const auto& groupset = group_sets[groupset_num];
  for (unsigned int i=0; i<options.num_modes; ++i)
  {
    std::string file_name = file_base_name + std::to_string(i) + ".data";
    ReadSinglePODModeBinary(file_name,groupset,U[i]);
  }

  WriteTotalInteractionOperator(file_base_name,groupset);
}

//###################################################################
/**Reads a single POD mode into the designated vector.*/
void LinearBoltzmann::AffineDecompositioner::
  ReadSinglePODModeBinary(const std::string &file_name,
                          const LBSGroupset& groupset,
                          std::vector<double> &U_i)
{
  auto& chi_log = ChiLog::GetInstance();
  U_i.clear();

  //============================================= Open file
  std::ifstream file(file_name,
                     std::ofstream::binary | //binary file
                     std::ofstream::in);     //no accidental writing

  //============================================= Check file is open
  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLWARNING)
      << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;
  auto fe = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(discretization);
  if (not fe)
  {
    file.close();
    chi_log.Log(LOG_ALLWARNING) << "Angular flux file reading cancelled "
                                   "because a spatial discretization has not "
                                   "been initialized.";
    return;
  }

  size_t num_local_nodes   = discretization->GetNumLocalDOFs(NODES_ONLY);
  size_t num_angles        = groupset.quadrature->abscissae.size();
  size_t num_groups        = groupset.groups.size();
  size_t num_local_dofs    = groupset.num_psi_unknowns_local;
  auto   dof_handler       = groupset.psi_uk_man;

  size_t file_num_local_nodes;
  size_t file_num_angles     ;
  size_t file_num_groups     ;
  size_t file_num_local_dofs ;


  //============================================= Read header
  chi_log.Log() << "Reading POD-mode file " << file_name;
  char header_bytes[320]; header_bytes[319] = '\0';
  file.read(header_bytes,319);

  file.read((char*)&file_num_local_nodes, sizeof(size_t));
  file.read((char*)&file_num_angles     , sizeof(size_t));
  file.read((char*)&file_num_groups     , sizeof(size_t));
  file.read((char*)&file_num_local_dofs , sizeof(size_t));

  //============================================= Check compatibility
  if (file_num_local_nodes != num_local_nodes or
      file_num_angles      != num_angles      or
      file_num_groups      != num_groups      or
      file_num_local_dofs  != num_local_dofs)
  {
    std::stringstream outstr;
    outstr << "num_local_nodes: " << file_num_local_nodes << "\n";
    outstr << "num_angles     : " << file_num_angles << "\n";
    outstr << "num_groups     : " << file_num_groups << "\n";
    outstr << "num_local_dofs : " << file_num_local_dofs << "\n";
    chi_log.Log(LOG_ALL)
      << "Incompatible DOF data found in file " << file_name << "\n"
      << outstr.str();
    file.close();
    return;
  }

  //============================================= Commit to reading the file
  U_i.reserve(file_num_local_dofs);
  std::set<uint64_t> cells_touched;
  for (size_t dof=0; dof < file_num_local_dofs; ++dof)
  {
    uint64_t     cell_global_id;
    unsigned int node;
    unsigned int angle_num;
    unsigned int group;
    double       psi_value;

    file.read((char*)&cell_global_id,sizeof(uint64_t));
    file.read((char*)&node          ,sizeof(unsigned int));
    file.read((char*)&angle_num     ,sizeof(unsigned int));
    file.read((char*)&group         ,sizeof(unsigned int));
    file.read((char*)&psi_value     ,sizeof(double));

    cells_touched.insert(cell_global_id);

    const auto& cell = grid->cells[cell_global_id];

    size_t imap = fe->MapDOFLocal(cell,node,dof_handler,angle_num,group);

    U_i[imap] = psi_value;
  }

  chi_log.Log(LOG_ALL) << "Done reading POD-mode file " << file_name
                       << ". Number of cells read: " << cells_touched.size();

  //============================================= Clean-up
  file.close();
}
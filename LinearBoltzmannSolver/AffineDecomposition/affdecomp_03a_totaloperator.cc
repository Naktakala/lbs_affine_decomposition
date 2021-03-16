#include "affine_decompositioner.h"

#include "ChiMath/dynamic_matrix.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"
#include "ChiMath/dynamic_vector.h"

#include "chi_log.h"

//###################################################################
/**Creates and writes the total interaction operator.*/
void LinearBoltzmann::AffineDecompositioner::
  WriteTotalInteractionOperator(const std::string &file_base_name,
                                const LBSGroupset &groupset)
{
  auto& chi_log = ChiLog::GetInstance();
  chi_log.Log() << "Creating total interaction operator.";

  //============================================= Get FE spatial discretization
  auto fe = std::dynamic_pointer_cast<SpatialDiscretization_FE>(discretization);
  if (not fe) throw std::logic_error(std::string(__FUNCTION__) +
                                     " error using spatial discretization.");

  //============================================= Define some types
  typedef chi_math::DynamicMatrix<double> Matrix;
  typedef chi_math::DynamicVector<double> Vector;
  typedef chi_math::finite_element::UnitIntegralData CellMatrices;
  typedef std::vector<Matrix> VecMats;
  typedef std::vector<VecMats> MatMats;
  auto dof_handler = groupset.psi_uk_man;

  //============================================= Lambda to extract DOFs
  auto Extractor_DOFs = [fe,dof_handler](
    const chi_mesh::Cell& cell, const CellMatrices& fe_matrices,
    const std::vector<double>& ref_U, size_t d, unsigned int g)
  {
    Vector U_local(fe_matrices.NumNodes(),0.0);

    for (unsigned n=0; n < fe_matrices.NumNodes(); ++n) //node
    {
      size_t dof_index = fe->MapDOFLocal(cell, n, dof_handler, d, g);
      U_local[n] = ref_U[dof_index];
    }

    return U_local;
  };

  //============================================= Open the file
  const std::string file_name = file_base_name + ".totalop";
  std::ofstream file(file_name, std::ofstream::out);

  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLERROR) << __FUNCTION__ << "Failed to open " << file_name;
    return;
  } //Puke
  file << "Data structure:\n"
       << "Material k Direction d Row i Column j Group g Value M_kdijg\n";

  //============================================= Start building
  const size_t r          = options.num_modes;
  const size_t num_angles = groupset.quadrature->abscissae.size();

  for (auto& k : unique_material_ids)
    for (size_t d=0; d<num_angles; ++d) //direction
      for (size_t i=0; i<r; ++i)
        for (size_t j=0; j<r; ++j)
          for (unsigned int g=0; g<groupset.groups.size(); ++g)
          {
            double M_kdijg = 0.0;

            for (const auto& cell : grid->local_cells)
            {
              if (cell.material_id != k) continue;

              const auto&   fe_matrices = fe->GetUnitIntegrals(cell);

              Matrix mass_C(fe_matrices.GetIntV_shapeI_shapeJ());
              Vector U_idc = Extractor_DOFs(cell,fe_matrices,U[i],d,g);
              Vector U_jdc = Extractor_DOFs(cell,fe_matrices,U[j],d,g);

              Vector aux_i = mass_C * U_idc;
              Vector aux_j = mass_C * U_jdc;

              if (options.projection_mode == ProjectionMode::Petrov_Galerkin)
                M_kdijg += aux_i.Dot(aux_j);
              else
                M_kdijg += U_idc.Dot(aux_j);

            }//for cell

            file << k << " " << d << " "
                 << i << " " << j << " "
                 << g << " " << M_kdijg << "\n";
          }//for g

  file.close();

  chi_log.Log() << "Done creating total interaction operator. "
                   "File saved to " << file_name;
}
#include "affine_decompositioner.h"

#include "ChiMath/dynamic_matrix.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"
#include "ChiMath/dynamic_vector.h"

#include "chi_log.h"

//###################################################################
/**Creates and writes the total interaction operator.*/
void LinearBoltzmann::AffineDecompositioner::
  WriteStreamingOperator(const std::string &file_base_name,
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
    const chi_mesh::Cell& cell, const size_t num_nodes,
    const std::vector<double>& ref_U, size_t d, unsigned int g)
  {
    Vector U_local(num_nodes,0.0);

    for (unsigned n=0; n < num_nodes; ++n) //node
    {
      size_t dof_index = fe->MapDOFLocal(cell, n, dof_handler, d, g);
      U_local[n] = ref_U[dof_index];
    }

    return U_local;
  };

  //============================================= Open the file
  const std::string file_name = file_base_name + ".streamop";
  std::ofstream file(file_name, std::ofstream::out);

  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLERROR) << __FUNCTION__ << "Failed to open " << file_name;
    return;
  } //Puke

  //============================================= Start building
  const size_t r          = options.num_modes;
  const size_t num_angles = groupset.quadrature->abscissae.size();

  for (unsigned int g=0; g<groupset.groups.size(); ++g)
    for (size_t i=0; i<r; ++i)
      for (size_t j=0; j<r; ++j)
        for (size_t d=0; d<num_angles; ++d) //direction
          for (const auto& cell : grid->local_cells)
          {
            const auto&   fe_matrices = fe->GetUnitIntegrals(cell);
            size_t num_nodes = fe_matrices.NumNodes();

            Matrix grad_C(num_nodes,num_nodes,0.0);
            const auto& G = fe_matrices.GetIntV_shapeI_gradshapeJ();
            const auto& omega = groupset.quadrature->omegas[d];

            for (unsigned ni=0; ni<num_nodes; ++ni)
              for (unsigned nj=0; nj<num_nodes; ++nj)
                grad_C[ni][nj] = omega.Dot(G[ni][nj]);

            Vector U_idc = Extractor_DOFs(cell,num_nodes,U[i],d,g);
            Vector U_jdc = Extractor_DOFs(cell,num_nodes,U[j],d,g);

            Vector aux_i = grad_C * U_idc;
            Vector aux_j = grad_C * U_jdc;

            size_t num_faces = cell.faces.size();

            std::vector<Vector> U_idc_upstream(num_faces,Vector(num_nodes,0.0));
            std::vector<Vector> U_jdc_upstream(num_faces,Vector(num_nodes,0.0));

            std::vector<Vector> aux_face_i(num_faces,Vector(num_nodes,0.0));
            std::vector<Vector> aux_face_j(num_faces,Vector(num_nodes,0.0));

            for (int f=0; f<num_faces; ++f)
              if (cell.faces[f].has_neighbor)
                if (cell.faces[f].normal.Dot(omega) < 0.0)
                {
                  const auto&  adj_cell        = grid->cells[cell.faces[f].neighbor_id];
                  const auto&  adj_fe_matrices = fe->GetUnitIntegrals(adj_cell);
                  const size_t adj_num_nodes   = adj_fe_matrices.NumNodes();

                  Matrix mass_F(fe_matrices.GetIntS_shapeI_shapeJ()[f]);

                  U_idc_upstream[f] = Extractor_DOFs(adj_cell,adj_num_nodes,U[i],d,g);
                  U_jdc_upstream[f] = Extractor_DOFs(adj_cell,adj_num_nodes,U[j],d,g);

                  aux_face_i[f] = mass_F * U_idc_upstream[f];
                  aux_face_j[f] = mass_F * U_jdc_upstream[f];
                }//if mu<0.0

            double G_dkij;
            if (options.projection_mode == ProjectionMode::Petrov_Galerkin)
            {
              G_dkij = aux_i.Dot(aux_j);
              for (int f=0; f<num_faces; ++f)
                if (cell.faces[f].has_neighbor)
                  if (cell.faces[f].normal.Dot(omega) < 0.0)
                    G_dkij += aux_face_i[f].Dot(aux_face_j[f]);
            }
            else
            {
              G_dkij = U_idc.Dot(aux_j);
              for (int f=0; f<num_faces; ++f)
                if (cell.faces[f].has_neighbor)
                  if (cell.faces[f].normal.Dot(omega) < 0.0)
                    G_dkij += U_idc_upstream[f].Dot(aux_face_j[f]);
            }

            file << d << " " << g << " "
                 << i << " " << j << " " << G_dkij << "\n";
          }//for cell

  file.close();

  chi_log.Log() << "Done creating total interaction operator. "
                   "File saved to " << file_name;
}
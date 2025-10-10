#include "linear_elastic_material.h"

namespace Solid
{
  template <int dim>
  LinearElasticMaterial<dim>::LinearElasticMaterial(double young,
                                                    double poisson,
                                                    double rho,
                                                    double eta,
                                                    std::string material_type,
                                                    dealii::Tensor<1, dim> initial_fiber,
                                                    double young1,
                                                    double young2,
                                                    double poisson12,
                                                    double poisson23,
                                                    double shear12
                                                    )
    : Material<dim>(rho), E(young), nu(poisson), eta(eta), material_type(material_type), initial_fiber(initial_fiber), E1(young1), E2(young2), nu12(poisson12), nu23(poisson23), G12(shear12)
  {
    this->lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    this->mu = E / (2 * (1 + nu));
  }

  template <int dim>
  dealii::SymmetricTensor<4, dim>
  LinearElasticMaterial<dim>::get_elasticity() const
  {
  //create elasticity tensor for any material
  dealii::SymmetricTensor<4, dim> elasticity;
  
  //Create isotropic material
  if (this->material_type == "Isotropic")
  {
    for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    elasticity[i][j][k][l] =
                      (i == k && j == l ? this->mu : 0.0) +
                      (i == l && j == k ? this->mu : 0.0) +
                      (i == j && k == l ? this->lambda : 0.0);
                  }
              }
          }
      }
  //Create planar isotropic material
  }
  else if(this->material_type == "PlanarIsotropic")
  {
    //find k before filling elasticity tensor
    // Pre-compute constant values outside the loop to avoid repeated calculations
    const double constk = 1 - 2*(this->E2*(1+this->nu23)*pow(this->nu12,2))/this->E1 - pow(this->nu23,2);

    // define elasticity components (only computed once instead of in every iteration)
    // Cmn --> m and n represent the resp. matrix entry
    const double C11 = this->E1*(1-pow(this->nu23,2))/constk;
    const double C22 = this->E2*(1-this->E2/this->E1*pow(this->nu12,2))/constk;
    const double C12 = this->E2*this->nu12*(1+this->nu23)/constk;
    const double C23 = this->E2*(this->E2/this->E1*pow(this->nu12,2)+this->nu23)/constk;
    const double C44 = (this->E1*(1-pow(this->nu23,2))-this->E2*this->nu12*(1+this->nu23))/(2*constk);
    
    // lambda function to convert tensor indices (i,j) to Voigt notation
    // returns values 1-6, with 0 indicating an invalid combination
    auto to_voigt = [](unsigned int i, unsigned int j) -> int {
        if (i == 0 && j == 0) return 1;
        if (i == 1 && j == 1) return 2;
        if (i == 2 && j == 2) return 3;
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 4;
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 5;
        if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 6;
        return 0;
    };
    
    for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            // Compute Voigt index m once per (i,j) pair instead of for every (i,j,k,l)
            const int m = to_voigt(i, j);
            
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    // Compute Voigt index n once per (k,l) pair
                    const int n = to_voigt(k, l);
                    
                    // Use pre-computed values and simplified conditional logic
                    if (m == 1 && n == 1)
                      {
                        elasticity[i][j][k][l] = C11;
                      }
                    else if ((m == 2 && n == 2) || (m == 3 && n == 3))
                      {
                        elasticity[i][j][k][l] = C22;
                      }
                    else if (m == 1 && (n == 2 || n == 3))
                      {
                        elasticity[i][j][k][l] = C12;
                        elasticity[k][l][i][j] = C12;  // Symmetry
                      }
                    else if (m == 2 && n == 3)
                      {
                        elasticity[i][j][k][l] = C23;
                        elasticity[k][l][i][j] = C23;  // Symmetry
                      }
                    else if (m == 4 && n == 4)
                      {
                        elasticity[i][j][k][l] = C44;
                      }
                    else if ((m == 5 && n == 5) || (m == 6 && n == 6))
                      {
                        elasticity[i][j][k][l] = this->G12;
                      }
                  }
              }
          }
      }
  }
  return elasticity;
  }

  template <int dim>
  dealii::SymmetricTensor<4, dim>
  LinearElasticMaterial<dim>::get_viscosity() const
  {
    dealii::SymmetricTensor<4, dim> viscosity;
    for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    viscosity[i][j][k][l] =
                      (i == k && j == l ? this->eta / 2 : 0.0) +
                      (i == l && j == k ? this->eta / 2 : 0.0);
                  }
              }
          }
      }
    return viscosity;
  }

  //TODO can look into default
  //std optional SymmetricTensor
  //set field as optional
  //Can also be used to rotate viscosity if anisotropic viscosity is required
  template <int dim>
  dealii::SymmetricTensor<4, dim>
  LinearElasticMaterial<dim>::rotate_tensor(dealii::Tensor<1, dim> current_fiber, dealii::SymmetricTensor<4, dim> elasticity) const
  {
    dealii::Tensor<1, dim> current_fiberxy, xaxis;
    current_fiberxy[0] = current_fiber[0];
    current_fiberxy[1] = current_fiber[1];
    xaxis[0] = 1;

    double theta = current_fiberxy.norm() == 0 ? 0 : dealii::Physics::VectorRelations::angle(current_fiberxy, xaxis);
    //double theta = fiberxy.norm() == 0 ? 0 : (fiber[1] > 0 ? dealii::Physics::VectorRelations::angle(fiberxy, xaxis) : -dealii::Physics::VectorRelations::angle(fiberxy, xaxis));
    theta = current_fiber[1] > 0 ? theta : -theta;
  
    dealii::Tensor<2, dim> R = dealii::unit_symmetric_tensor<dim>(), Rz = dealii::unit_symmetric_tensor<dim>();
  
    //Rotate about z axis, same process for both 2d and 3d 
    Rz[0][0] = cos(theta);
    Rz[0][1] = -sin(theta);
    Rz[1][0] = sin(theta);
    Rz[1][1] = cos(theta);

    //Rotation tensor creation depends on dimension of simulations
    if(dim == 2)
    {
      //2d case can be rotated about the z axis, since the whole simulation is within the xy plane
      R = Rz;
    }
    else if(dim == 3)
    {
      double phi = current_fiberxy.norm() == 0 ? M_PI/2 : dealii::Physics::VectorRelations::angle(current_fiber, current_fiberxy);
      phi = current_fiber[2] > 0 ? phi : -phi;

      dealii::Tensor<2, dim> Ry = dealii::unit_symmetric_tensor<dim>();

      Ry[0][0] = cos(phi);
      Ry[0][2] = -sin(phi);
      Ry[2][0] = sin(phi);
      Ry[2][2] = cos(phi);
    
      R=Rz*Ry;
    }
    //Create temporary asymmetric tensor for multiplications then converting to symmetric after
    dealii::Tensor<4, dim> temp, temp2;

    //Rotate about i and l
    temp = R*elasticity*transpose(R);
    //Flip so j and k are on outside
    for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
      {
        for (unsigned int k = 0; k < dim; k++)
        {
          for (unsigned int l = 0; l < dim; l++)
          {
            temp2[j][i][l][k] = temp[i][j][k][l];
          }
        }
      }
    }
    //Rotate about j and k
    temp2 = R*temp2*transpose(R);
    //Rotate back to create original temp tensor
    for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
      {
        for (unsigned int k = 0; k < dim; k++)
        {
          for (unsigned int l = 0; l < dim; l++)
          {
            elasticity[i][j][k][l] = temp2[j][i][l][k];
          }
        }
      }
    }

    return elasticity;
  }
  
  //Might move to solid_solver or utilities 
  template <int dim>
  dealii::Tensor<1, dim>
  LinearElasticMaterial<dim>::get_current_fiber_direction(dealii::Tensor<2, dim> deformation_gradient) const
  {
    dealii::Tensor<1, dim> current_fiber = (deformation_gradient + dealii::unit_symmetric_tensor<dim>())*this->initial_fiber;
    return current_fiber;
  }


  // explicit instantiation
  template class LinearElasticMaterial<2>;
  template class LinearElasticMaterial<3>;
} // namespace Solid

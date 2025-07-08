#ifndef LINEAR_ELASTIC_MATERIAL
#define LINEAR_ELASTIC_MATERIAL
#define _USE_MATH_DEFINES

#include <deal.II/base/symmetric_tensor.h>
//needed for vector relations for anisotropy
#include <deal.II/physics/vector_relations.h>
#include <deal.II/base/tensor.h>
#include <math.h>

#include "material.h"

namespace Solid
{
  /*! \brief Linear elastic material.
   */
  template <int dim>
  class LinearElasticMaterial : public Material<dim>
  {
  public:
    LinearElasticMaterial()
      : Material<dim>(), E(0.0), nu(0.0), lambda(0.0), mu(0.0), eta(0.0), material_type("Isotropic"), initial_fiber(), E1(0.0), E2(0.0), nu12(0.0), nu23(0.0), G12(0.0)
    {
    }
    /**
     * Constructor using Young's modulus and Poisson's ratio.
     */
    LinearElasticMaterial(double, double, double, double, std::string, dealii::Tensor<1,dim>, double, double, double, double, double);
    dealii::SymmetricTensor<4, dim> get_elasticity() const;
    dealii::SymmetricTensor<4, dim> get_viscosity() const;
    dealii::SymmetricTensor<4, dim> rotate_tensor(dealii::Tensor<1, dim> current_fiber, dealii::SymmetricTensor<4, dim> elasticity) const;
    //dealii::Tensor<1, dim> get_current_fiber_direction(dealii::Tensor<2, dim> deformation_gradient, dealii::Tensor<1, dim> i_fiber) const;
    dealii::Tensor<1, dim> get_current_fiber_direction(dealii::Tensor<2, dim> deformation_gradient) const;
    std::string material_type;


  protected:
    double E;      //!< Young's modulus
    double nu;     //!< Poisson's ratio
    double lambda; //!< First lame parameter
    double mu;     //!< Second lame parameter
    double eta;    //!< Viscosity
    //Anisotropic material properties
    dealii::Tensor<1, dim> initial_fiber;
    double E1;
    double E2;
    double nu12;
    double nu23;
    double G12;
  };
} // namespace Solid

#endif

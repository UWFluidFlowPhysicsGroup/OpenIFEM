#ifndef LINEAR_ELASTIC_MATERIAL
#define LINEAR_ELASTIC_MATERIAL

#include <deal.II/base/symmetric_tensor.h>
//needed for vector relations for anisotropy
#include <deal.II/physics/vector_relations.h>
#include <deal.II/base/tensor.h>

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
      : Material<dim>(), E(0.0), nu(0.0), lambda(0.0), mu(0.0), eta(0.0), material_type("Isotropic"), fiber()
    {
    }
    /**
     * Constructor using Young's modulus and Poisson's ratio.
     */
    LinearElasticMaterial(double, double, double, double, std::string, std::vector<double>);
    dealii::SymmetricTensor<4, dim> get_elasticity() const;
    dealii::SymmetricTensor<4, dim> get_viscosity() const;
    dealii::SymmetricTensor<4, dim> rotate_tensor(dealii::Tensor<2, dim> grad_u, dealii::SymmetricTensor<4, dim> elasticity) const;

  protected:
    double E;      //!< Young's modulus
    double nu;     //!< Poisson's ratio
    double lambda; //!< First lame parameter
    double mu;     //!< Second lame parameter
    double eta;    //!< Viscosity
    std::string material_type;
    std::vector<double> fiber;
  };
} // namespace Solid

#endif

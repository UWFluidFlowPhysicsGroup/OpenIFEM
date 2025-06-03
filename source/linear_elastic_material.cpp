#include "linear_elastic_material.h"

namespace Solid
{
  template <int dim>
  LinearElasticMaterial<dim>::LinearElasticMaterial(double young,
                                                    double poisson,
                                                    double rho,
                                                    double eta)
    : Material<dim>(rho), E(young), nu(poisson), eta(eta)
  {
    this->lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    this->mu = E / (2 * (1 + nu));
  }

  template <int dim>
  dealii::SymmetricTensor<4, dim>
  LinearElasticMaterial<dim>::get_elasticity() const
  {
  //Variables and principal matrix creation
  double E1 = 0, E2 = 0, G12 = 0, nu12 = 0, nu23 = 0;
  //TODO import variable values from parameters, hard coding for now
  E1 = 131e9; //Pa
  E2 = 10.3e9; //Pa
  nu12 = 0.22;
  nu23 = nu12;
  G12 = 6.9e9;
  
  dealii::Tensor<1, dim> fiber, fiberxy, xaxis;
  fiber[0] = 1;
  fiber[1] = 0;
  //define z axis only for 3D case (otherwise out of bounds)
  if (dim == 3){
    fiber[2] = 0;
  }

  fiberxy[0] = fiber[0];
  fiberxy[1] = fiber[1];
  xaxis[0] = 1;
    

  //create elasticity tensor in principal coordinates
  dealii::SymmetricTensor<4, dim> elasticity;
  
  /*
  TODO create JSON files that map the results of each if statement condition in 4d space, comprised of 1 and 0s (pseudo identity matrix)
  Then just run through each material type, multiply by that mapping matrix and add to elasticity matrix
  */
  //find k before filling elasticity tensor
  const double constk = 1-2*(E2*(1+nu23)*pow(nu12,2))/E1-pow(nu23,2);
  //TODO: defining G12 same as later on in the for loop, isotropic test case rn
  //G12 = (E1*(1-pow(nu23,2))-E2*nu12*(1+nu23))/(2*constk);

  //SymmetricTensor object automatically applies symmetries in ijkl=jikl=ijlk, but still need to manually input ijkl=klij
  int m=0, n=0;
  for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
          {
            for (unsigned int k = 0; k < dim; ++k)
              {
                for (unsigned int l = 0; l < dim; ++l)
                  {
                    //Create temporary indices m, n to write equivalent Voigt notation for ijkl, makes it easier to process and visualize each case
                    //TODO make another function to convert Einstein to Voigt? (input a, b and output c to avoid duplicate loops for m, n?)
                    m=0;
                    n=0;

                    if (i==0 && j==0) {
                      m=1;
                    } else if (i==1 && j==1) {
                      m=2;
                    } else if (i==2 && j==2) {
                      m=3;
                    } else if ((i==1 && j==2)||((i==2 && j==1))) {
                      m=4;
                    } else if ((i==0 && j==2)||((i==2 && j==0))) {
                      m=5;
                    } else if ((i==0 && j==1)||((i==1 && j==0))) {
                      m=6;
                    }

                    if (k==0 && l==0) {
                      n=1;
                    } else if (k==1 && l==1) {
                      n=2;
                    } else if (k==2 && l==2) {
                      n=3;
                    } else if ((k==1 && l==2)||((k==2 && l==1))) {
                      n=4;
                    } else if ((k==0 && l==2)||((k==2 && l==0))) {
                      n=5;
                    } else if ((k==0 && l==1)||((k==1 && l==0))) {
                      n=6;
                    }

                    if(m==1 && n==1)
                    {
                      elasticity[i][j][k][l] = E1*(1-pow(nu23,2))/constk;
                    }
                    else if((m==2 && n==2)||(m==3 && n==3))
                    {
                      elasticity[i][j][k][l] = E2*(1-E2/E1*pow(nu12,2))/constk;
                    }
                    else if(m==1 && (n==2 || n==3))
                    {
                      elasticity[i][j][k][l] = E2*nu12*(1+nu23)/constk;
                      //applying symmetry along main diagonal (not automatically done for SymmetricTensor)
                      elasticity[k][l][i][j] = elasticity[i][j][k][l];
                    }
                    else if(m==2 && n==3)
                    {
                      elasticity[i][j][k][l] = E2*(E2/E1*pow(nu12,2)+nu23)/constk;
                      //applying symmetry along main diagonal (not automatically done for SymmetricTensor)
                      elasticity[k][l][i][j] = elasticity[i][j][k][l];
                    }
                    else if(m==4 && n==4)
                    {
                      //C11 and C12 are already known based on order of for loops (all of i=1 is done first), but writing explicitly just to be safe
                      elasticity[i][j][k][l] = (E1*(1-pow(nu23,2))-E2*nu12*(1+nu23))/(2*constk);
                    }
                    else if((m==5 && n==5)||(m==6 && n==6))
                    {
                      elasticity[i][j][k][l] = G12;
                    }
                  }
              }
          }
      }

      double theta;
    //fiberxy.norm is 0 only when xy length is 0, ie fiberxy is 0.
    if(fiberxy.norm() == 0){
      //if fiberxy is 0, assume no theta rotation required (fiber is purely vertical)
      theta = 0;
    }else{
      //find theta from angle between x axis and xy projection
      theta = dealii::Physics::VectorRelations::angle(fiberxy, xaxis);
    }
    
    //TODO find how to prepopulate R, Rz and Ry from identity tensor using SymmetricTensor::unit_symmetric_tensor() function
    //Prepopulating not really necessary, defaults to 0 value, so identity only ensures the "1" along main diagonal is filled for 3d
    dealii::Tensor<2, dim> R;
    dealii::Tensor<2, dim> Rz;
    //temporary identity matrix function (could be replaced with if(dim == 3){R[2][2] = 1})
    for (int i = 0; i < dim; i++){
      R[i][i] = 1;
      Rz[i][i] = 1;
    }

    //Rotate about z axis, same process for both 2d and 3d 
    //TODO find way to condense this part?
    Rz[0][0] = cos(theta);
    Rz[0][1] = -sin(theta);
    Rz[1][0] = sin(theta);
    Rz[1][1] = cos(theta);

    //Rotation tensor creation depends on dimension of simulations
    if(dim == 2){
      //2d case can be rotated about the z axis, since the whole simulation is within the xy plane
      R = Rz;
      //display rotation tensor for debugging
    }else if(dim == 3){
      double phi;
      if(fiberxy.norm() == 0){
        //if fiberxy is 0, then vector is purely vertical, which is assumed to be 90 degrees
        const double pi = 3.14159265358979323846;
        phi = pi/2;
      }else{
        //3d case needs phi to account for components in the z direction, which is found using xy projection and full fiber direction
        phi = dealii::Physics::VectorRelations::angle(fiber, fiberxy);
      }
      //std::cout << "theta = " << 57.3*theta << "\nphi = " << 57.3*phi << "\n";

      //same issues with Rz, identity tensor and condense assigning sin and cos values
      dealii::Tensor<2, dim> Ry;
      for (int i = 0; i < dim; i++){
        Ry[i][i] = 1;
      }
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
    for (unsigned int i = 0; i < dim; i++){
      for (unsigned int j = 0; j < dim; j++){
        for (unsigned int k = 0; k < dim; k++){
          for (unsigned int l = 0; l < dim; l++){
            temp2[j][i][l][k] = temp[i][j][k][l];
          }
        }
      }
    }
    //Rotate about j and k
    temp2 = R*temp2*transpose(R);
    //Rotate back to create original temp tensor
    for (unsigned int i = 0; i < dim; i++){
      for (unsigned int j = 0; j < dim; j++){
        for (unsigned int k = 0; k < dim; k++){
          for (unsigned int l = 0; l < dim; l++){
            temp[i][j][k][l] = temp2[j][i][l][k];
          }
        }
      }
    }

    //for loops required to move generic tensor object to symmetric object for output
    for (unsigned int i = 0; i < dim; i++){
      for (unsigned int j = 0; j < dim; j++){
        for (unsigned int k = 0; k < dim; k++){
          for (unsigned int l = 0; l < dim; l++){
            //for cases where ijkl=jikl=ijlk will be overwritten with the last entry, assumes it is already symmetric
            elasticity[i][j][k][l] = temp[i][j][k][l];
          }
        }
      }
    }
    return elasticity;
    /*
    dealii::SymmetricTensor<4, dim> elasticity;
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
    return elasticity;*/
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

  template <int dim>
  dealii::SymmetricTensor<4, dim>
  LinearElasticMaterial<dim>::rotate_tensor(dealii::Tensor<2, dim> grad_u, dealii::Tensor<4, dim> elasticity) const
  {
    //TODO replace this with getting initial fiber coordinates from parameters/material file
    //creating array to get fiber coordinates before creating tensor
    dealii::Tensor<1, dim> fiber;
    fiber[0] = 1;
    fiber[1] = 0;
    //define z axis only for 3D case (otherwise out of bounds)
    if (dim == 3){
      fiber[2] = 0;
    }
    //Define deformation gradient F
    //TODO change to use unit_symmetric_tensor() later
    dealii::Tensor<2, dim> F = grad_u;
    for (int i = 0; i < dim; i++){
      F[i][i]++;
    }
    //dealii::Tensor<2, dim> F = grad_u + dealii::SymmetricTensor<2, dim>::unit_symmetric_tensor();

    //rotate fiber direction from initial to current fiber direction 
    fiber = F*fiber;
    //std::cout << fiber[0] << ", " << fiber[1] << "\n";
    //Create projection of fiber onto xy plane to find angles
    //fiberxy will be the same as input fiber for 2d case, creating this way to avoid duplicate code for theta
    dealii::Tensor<1, dim> fiberxy;
    fiberxy[0] = fiber[0];
    fiberxy[1] = fiber[1];
    //Create vector of x axis direction
    dealii::Tensor<1, dim> xaxis;
    xaxis[0] = 1;
    
    double theta;
    //fiberxy.norm is 0 only when xy length is 0, ie fiberxy is 0.
    if(fiberxy.norm() == 0){
      //if fiberxy is 0, assume no theta rotation required (fiber is purely vertical)
      theta = 0;
    }else{
      //find theta from angle between x axis and xy projection
      theta = dealii::Physics::VectorRelations::angle(fiberxy, xaxis);
    }
    
    //TODO find how to prepopulate R, Rz and Ry from identity tensor using SymmetricTensor::unit_symmetric_tensor() function
    //Prepopulating not really necessary, defaults to 0 value, so identity only ensures the "1" along main diagonal is filled for 3d
    dealii::Tensor<2, dim> R;
    dealii::Tensor<2, dim> Rz;
    //temporary identity matrix function (could be replaced with if(dim == 3){R[2][2] = 1})
    for (int i = 0; i < dim; i++){
      R[i][i] = 1;
      Rz[i][i] = 1;
    }

    //Rotate about z axis, same process for both 2d and 3d 
    //TODO find way to condense this part?
    Rz[0][0] = cos(theta);
    Rz[0][1] = -sin(theta);
    Rz[1][0] = sin(theta);
    Rz[1][1] = cos(theta);

    //Rotation tensor creation depends on dimension of simulations
    if(dim == 2){
      //2d case can be rotated about the z axis, since the whole simulation is within the xy plane
      R = Rz;
      //display rotation tensor for debugging
    }else if(dim == 3){
      double phi;
      if(fiberxy.norm() == 0){
        //if fiberxy is 0, then vector is purely vertical, which is assumed to be 90 degrees
        const double pi = 3.14159265358979323846;
        phi = pi/2;
      }else{
        //3d case needs phi to account for components in the z direction, which is found using xy projection and full fiber direction
        phi = dealii::Physics::VectorRelations::angle(fiber, fiberxy);
      }
      //std::cout << "theta = " << 57.3*theta << "\nphi = " << 57.3*phi << "\n";

      //same issues with Rz, identity tensor and condense assigning sin and cos values
      dealii::Tensor<2, dim> Ry;
      for (int i = 0; i < dim; i++){
        Ry[i][i] = 1;
      }
      Ry[0][0] = cos(phi);
      Ry[0][2] = -sin(phi);
      Ry[2][0] = sin(phi);
      Ry[2][2] = cos(phi);
    
      R=Rz*Ry;
    }

    //Create temporary asymmetric tensor for multiplications then converting to symmetric after
    dealii::Tensor<4, dim> temp, temp2;
    dealii::SymmetricTensor<4, dim> elasticityCartesian;

    //Rotate about i and l
    temp = R*elasticity*transpose(R);
    //Flip so j and k are on outside
    for (unsigned int i = 0; i < dim; i++){
      for (unsigned int j = 0; j < dim; j++){
        for (unsigned int k = 0; k < dim; k++){
          for (unsigned int l = 0; l < dim; l++){
            temp2[j][i][l][k] = temp[i][j][k][l];
          }
        }
      }
    }
    //Rotate about j and k
    temp2 = R*temp2*transpose(R);
    //Rotate back to create original temp tensor
    for (unsigned int i = 0; i < dim; i++){
      for (unsigned int j = 0; j < dim; j++){
        for (unsigned int k = 0; k < dim; k++){
          for (unsigned int l = 0; l < dim; l++){
            temp[i][j][k][l] = temp2[j][i][l][k];
          }
        }
      }
    }

    //for loops required to move generic tensor object to symmetric object for output
    for (unsigned int i = 0; i < dim; i++){
      for (unsigned int j = 0; j < dim; j++){
        for (unsigned int k = 0; k < dim; k++){
          for (unsigned int l = 0; l < dim; l++){
            //for cases where ijkl=jikl=ijlk will be overwritten with the last entry, assumes it is already symmetric
            elasticityCartesian[i][j][k][l] = temp[i][j][k][l];
          }
        }
      }
    }
    return elasticityCartesian;
  }

  // explicit instantiation
  template class LinearElasticMaterial<2>;
  template class LinearElasticMaterial<3>;
} // namespace Solid

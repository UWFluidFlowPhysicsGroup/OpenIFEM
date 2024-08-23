/* 
Seemed interesting to implement when starting to set up simulation program, but upon closer inspection is not that useful
Dimension Independent Programming (DIP) is supposed to assist with method generation that can handle both 2D and 3D objects
Upon reflection of the project scope, DIP is not needed since the extent of dimensionless programming will only include
importing the file and reading if it is 2d or 3d and extruding the 2d shape to 3d.

Keeping this code for future reference
*/

//import dealII libraries
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

//import OpenIFEM libraries
#include "parameters.h"
#include "utilities.h"


//import c++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>

using namespace dealii;

//Vars in unnamed namespace to avoid reading from other files
template <int dim>
class Sim{
public:
  //extern template class Solid::LinearElasticity<dim>;
  Sim();
  int loadMesh(std::string meshName);
  Triangulation<3> extrude();
  int refine(int refinement);
private:
  //TODO change to dimension independent programming using <dim> (step 4)
  Triangulation<dim> tria;
  DoFHandler<dim>    dof_handler;
  GridIn<dim> gridIn;
  
};

namespace {
const std::string squareMeshName = "squareMesh";
const std::string vocalFoldMeshName = "vocalFold2D";
const std::string meshPath = "meshes/";
std::string paramsPath("parameters.prm");

GridOut gridOut;
}

template <int dim>
Sim<dim>::Sim()
  : dof_handler(tria)
{}

//imports a mesh and outputs svg file in the XY plane
template <int dim>
int Sim<dim>::loadMesh(std::string meshName){
  //define 2D GridIn object to receive 2d mesh
  gridIn.attach_triangulation(tria);
  //identifies mesh to be imported from meshes folder
  std::ifstream f(meshPath + meshName + ".msh");
  //checks if desired mesh can be read
  if (!f){
    //Display error handler that file cannot be found
    std::cerr << "----------------------------------------------------"
              << "ERROR FINDING MESH FILE " << meshName
              << "----------------------------------------------------";
    //return to kill the class
    return -1;
  }
  //imports mesh from selected area
  gridIn.read_msh(f);
  return 1;
}

//takes input 2d mesh from before and extrudes to a 3d shape, exports shape to .geo file
template <int dim>
Triangulation<3> Sim<dim>::extrude(){
  Triangulation<3> tria3d;  
  //2d input, number of slices, height, output height, output triangulation
  GridGenerator::extrude_triangulation(tria, 7, 12.0, tria3d);
  std::ofstream out(meshPath + "vocalFold3D.msh");
  gridOut.write_msh(tria3d, out);
  return tria3d;
}

template <int dim>
int Sim<dim>::refine(int refinement){
  //refine_global is set to 1 subdivision because mesh is subdivided from previous loop 
  //one further step into refinement
  tria.refine_global(1);
  //output the refined mesh with a different name based on refinement level
  std::ofstream out(meshPath + "vocalFold3d" + std::to_string(refinement) + ".msh");
  gridOut.write_msh(tria, out);
  return 1;
}

int main(){
  //read parameters file to determine the dimensions present
  Parameters::AllParameters params(paramsPath);
  
  if (params.dimension == 2){
    Sim<2> DIPTest;
    DIPTest.loadMesh(vocalFoldMeshName);
    //TODO do something with the returned extruded 3d mesh, maybe create a Sim<3> object?
    DIPTest.extrude();
    for(int i = 1; i <= 2; i++){
      DIPTest.refine(i);
    }  
  } else if (params.dimension == 3){
    Sim<3> DIPTest;
    DIPTest.loadMesh(vocalFoldMeshName);
    //no extrude since it is already 3d
    for(int i = 1; i <= 2; i++){
      DIPTest.refine(i);
    }
  } else {
    std::cerr << "Cannot find dimension from parameters file" << std::endl
              << "Check if " << paramsPath << "exists";
    return 1;
  }
}

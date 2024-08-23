//import dealII libraries
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

//import OpenIFEM libraries
#include "linear_elasticity.h"
#include "parameters.h"
#include "utilities.h"


//import c++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>

extern template class Solid::LinearElasticity<2>;
extern template class Solid::LinearElasticity<3>;

using namespace dealii;

namespace {
const std::string squareMeshName = "squareMesh";
const std::string vocalFoldMeshName = "vocalFold2D";
const std::string cubeMeshName = "cubeMesh";
const std::string meshPath = "meshes/";
std::string paramsPath2d("parameters2d.prm");
std::string paramsPath3d("parameters3d.prm");

//define objects for both 2d and 3d mesh manipulation
Triangulation<2> tria2d;
Triangulation<3,3> tria3d;

GridIn<2> gridIn2d;
GridIn<3> gridIn3d;

GridOut gridOut;
}

//imports a mesh and outputs svg file in the XY plane
int loadMesh2d(std::string meshName){
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

  //TODO write if statement to check if 2d or 3d mesh being imported, maybe try catch as 3d and have 2d in the catch segment and merge with extrude class?
  //define 2D GridIn object to receive 2d mesh
  gridIn2d.attach_triangulation(tria2d);
  //imports mesh from selected area
  gridIn2d.read_msh(f);
  //prepares squareMesh.svg file
  std::ofstream out(meshName + ".svg");
  //writes refined mesh to svg in XY plane
  gridOut.write_svg(tria2d, out);
  return 1;
}

//imports a 3D mesh and outputs svg file in the XY plane
int loadMesh3d(std::string meshName){
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

  //TODO write if statement to check if 2d or 3d mesh being imported, maybe try catch as 3d and have 2d in the catch segment and merge with extrude class?
  //define 2D GridIn object to receive 2d mesh
  gridIn3d.attach_triangulation(tria3d);
  //imports mesh from selected area
  gridIn3d.read_msh(f);
  
  //prepares squareMesh.svg file
  //std::ofstream out(meshName + ".svg");
  return 1;
}

int importParams2d(std::string paramName){
    //import params from .prm file and assign to 2d square
    Parameters::AllParameters params(paramName);
    Solid::LinearElasticity<2> solid(tria2d, params);
    solid.run();
    
    return 0;
}

int importParams3d(std::string paramName){
    //import params from .prm file
    Parameters::AllParameters params(paramName);
    Solid::LinearElasticity<3> solid(tria3d, params);
    solid.run();
    
    return 0;
}

//takes input 2d mesh from before and extrudes to a 3d shape, exports shape to .geo file
int extrude(){  
  //2d input, number of slices, height, output height, output triangulation
  GridGenerator::extrude_triangulation(tria2d, 7, 12.0, tria3d);
  std::ofstream out(meshPath + "vocalFold3d.msh");
  gridOut.write_msh(tria3d, out);
  return 1;
}

int refine(int i){
  //refine_global is set to 1 subdivision because mesh is subdivided from previous loop 
  //one further step into refinement
  tria3d.refine_global(1);
  //output the refined mesh with a different name based on refinement levevl
  std::ofstream out(meshPath + "vocalFold3d" + std::to_string(i) + ".msh");
  gridOut.write_msh(tria3d, out);
  return 1;
}

//pulled from OpenIFEM Solid_beam_bending_linearelastic test
int generateMesh(){
  double L = 8.0, H = 1.0;
  std::ofstream out(meshPath + "solid_beam_mesh.msh");

  Parameters::AllParameters params("solid_beam_bending_linearelastic.prm");
  
  Triangulation<2> tria;
          dealii::GridGenerator::subdivided_hyper_rectangle(
            tria, {32, 4}, Point<2>(0, 0), Point<2>(L, H), true);
  Solid::LinearElasticity<2> solidBeam(tria, params);
  solidBeam.run();
  return 1;
}

int main(){

  loadMesh3d(cubeMeshName);
  //generateMesh();
  importParams3d(paramsPath3d);
  
  
  //loadMesh2d(squareMeshName);
  //extrude();
  //for(int i = 1; i <= 3; i++){
  //  refine(i);
  //}
  //importParams2d(paramsPath2d);
  
}

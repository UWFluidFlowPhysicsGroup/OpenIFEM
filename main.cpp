//import dealII libraries
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

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

namespace {
const std::string squareMeshName = "squareMesh";
const std::string vocalFoldMeshName = "vocalFold2D";
const std::string cubeMeshName = "cubeMesh";
const std::string meshPath = "meshes/";

//define objects for both 2d and 3d mesh manipulation
Triangulation<2> tria2d;
Triangulation<3> tria3d;

GridIn<2> gridIn2d;
GridIn<3> gridIn3d;

GridOut gridOut;
}

//imports a mesh and outputs svg file in the XY plane
int loadMesh(std::string meshName){
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

//takes input 2d mesh from before and extrudes to a 3d shape, exports shape to .geo file
int extrude(){  
  //2d input, number of slices, height, output height, output triangulation
  GridGenerator::extrude_triangulation(tria2d, 7, 12.0, tria3d);
  std::ofstream out(meshPath + "vocalFold3D.msh");
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

int main(){
  loadMesh(squareMeshName);
  extrude();
  for(int i = 1; i <= 3; i++){
    refine(i);
  }
}

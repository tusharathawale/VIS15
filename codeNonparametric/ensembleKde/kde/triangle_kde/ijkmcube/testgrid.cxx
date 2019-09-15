
#include <iostream>

#include "ijkmcube_datastruct.h"

using namespace std;

int main() {

  int DIM3 = 3;
  int axis_size[] = { 5, 6, 7 };
  IJKMCUBE::MC_SCALAR_GRID grid(DIM3, axis_size);

  int num_cube_vertices = grid.NumCubeVertices();
  cerr << "Num cube vertices = " << num_cube_vertices << endl;
  cerr << "Axis increment: ";
  for (int d = 0; d < DIM3; d++) {
    cerr << " " << grid.AxisIncrement(d);
  }
  cerr << endl;
  cerr << "Cube vertex increment: ";
  for (int i = 0; i < num_cube_vertices; i++) {
    cerr << " " << grid.CubeVertexIncrement(i);
  }
  cerr << endl;
  cerr << "Unit cube vertex coordinates: ";
  for (int i = 0; i < num_cube_vertices; i++) {
    cerr << " (";
    for (int d = 0; d < DIM3; d++) {
      cerr << grid.UnitCubeCoord(i,d);
      if (d+1 < DIM3) { cerr << ","; };
    }
    cerr << ")";
  }
  cerr << endl;
}

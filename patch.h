#include "Cell.h"
#include <stdio.h>
#include <list>

bool Delete(Cell& c) {return (c.MustDelete == true);}

using namespace std;

class patch
{
 public:
  patch();
  patch(int nx_, int ny_);
  ~patch() {};
  

  int nx, ny, nz;
  list<Cell*> CellList;

  list<patch*> neighbour;

};


patch::patch(){


}

patch::patch(int nx_, int ny_)
{
nx = nx_, ny = ny_;
  list<Cell*> CellList;
}


#include <stdio.h>
#include <list>
#include <iostream> 
#include <stdlib.h>

using namespace std;

class Pltpair
{
 public:
  Pltpair();
 ~Pltpair() {};

  Pltpair(Cell* &plt1_, Cell* &plt2_, int level_); 

  Cell* plt1;
  Cell* plt2;
  int level;
  bool MustDelete;

  double age;

};


Pltpair::Pltpair(){


}


Pltpair::Pltpair(Cell* &plt1_, Cell* &plt2_, int level_){
plt1 = plt1_;
plt2 = plt2_;
level = level_;
age = 0;
MustDelete = false;
}






#include <list>
#include <stdio.h>
#include <stdlib.h>
#include "parameters_cell.h"
#include "Vect.h"
#include <iostream> // for the use of 'cout'
#include <vector> // vector containers
#include <math.h>
#include <cstdlib> // standard library
#define _USE_MATH_DEFINES
#define SQ(x) ((x) * (x))
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

using namespace std;
class Cell
{
 public:
  Cell();
  Cell(Vect2D pos_, double R0_, int type_, double FasL_, int infectious, int inf_, double age_, int card_, int diab_, int bp_, int cd_, int can_, int mobile);


  double *node_x;
  double *node_y;

  list<Cell*> aggCells;

  double age;
  double time;
  double R0;
  int num_nodes;
  int type;
  Vect2D f;
  Vect2D v;
  Vect2D pos;
  double r;
  double phaseG01;
  double ERK;
  double Fasl;
  double incubation_time;
  double hospitalization_time;
  double release_time;
  double death_time;
  bool MustDelete;
  bool MustDivide;
  int inf;
  int infectious;
  int symp;
  int immune;
  int isolated;
  int newInf;
  int die;
  int indirect;
  int mobile;
  double death_prob;


  int card, diab, bp, cd, can;
  void infected(double, int);
  void Step(double dt, int tot_case, double Ca);
  void Step2(double dt);
  void Move(double dt, int t);
  void Move2(double dt, int t);

};

void Cell::Step(double dt, int tot_case, double Ca)
{

    double sm = (double) rand()/RAND_MAX;
	

    if (e_thr*(Ca/100) >= sm && immune == 0) {infected(dt, tot_case); indirect = 1; }



}



void Cell::Step2(double dt)
{
    time += dt;

    if (time >= incubation_time - presymp && time < incubation_time - presymp + dt) {infectious = 1;}

    if (time >= incubation_time && time < incubation_time + dt) {symp = 1; newInf = 1;}

    if (time >= hospitalization_time && time < hospitalization_time + dt) {isolated = 1; infectious = 0;
	if (age >=18 && age < 40) death_prob = 0.002;
	if (age >=40 && age < 49) death_prob = 0.004;
	if (age >=50 && age < 59) death_prob = 0.013;
	if (age >=60 && age < 69) death_prob = 0.036;
	if (age >=70 && age < 79) death_prob = 0.08;
	if (age >=80 && age < 90) death_prob = 0.148;

	if (cd == 1) death_prob = max(death_prob, 0.13);

	if (diab == 1) death_prob = max(death_prob, 0.092);
	
	if (bp == 1) death_prob = max(death_prob, 0.084);

	if (cd == 1) death_prob = max(death_prob, 0.08);

	if (can == 1) death_prob = max(death_prob, 0.076);
	


	double death_rate = (double) rand()/RAND_MAX;
	
	if (death_rate <= death_prob) {die = 1;
		  const gsl_rng_type * T;
		  gsl_rng * r;	

		  gsl_rng_env_setup();
		 
		  T = gsl_rng_default;
		  r = gsl_rng_alloc (T);

		  unsigned long int mySeed = rand();

		  gsl_rng_set(r, mySeed);

		  death_time = 13 + 24*3*gsl_ran_gaussian(r, 1.0);

		  release_time = 10000000;

	 }

    }

    

    if (time >= release_time) {inf = 0; symp = 0; isolated = 0; infectious = 0; immune = 1;}

    if (die == 1 && time >= death_time) MustDelete = true;
    
}



void Cell::infected(double dt, int tot_case)
{


  const gsl_rng_type * T;
  gsl_rng * r;	

  gsl_rng_env_setup();
 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  unsigned long int mySeed = rand();

  gsl_rng_set(r, mySeed);

  incubation_time = 24*gsl_ran_lognormal(r, 1.621, 0.418);

  double rel_perturb = gsl_ran_gaussian(r, 1.0);

  double hosp_perturb = gsl_ran_gaussian(r, 1.0);




  hospitalization_time = incubation_time + 24*3.9 + hosp_perturb*24*1;

  release_time = hospitalization_time + 24*13.0 + 24*7*rel_perturb;

  cout << incubation_time << endl;

  gsl_rng_free (r);
  tot_case += 1;

  inf = 1;

}




void Cell::Move2(double dt,int t)
{
  if (f.x == 0) v.x = 0;
  else v.x = (v.x + dt*(f.x)  );
  if (f.y == 0) v.y = 0;
  else v.y = (v.y + dt*(f.y) );


  if (pos.x > 0 && pos.x < 500 && pos.y > 0 && pos.y < 500) pos+=v*dt;

  if (pos.x <= 0) {pos.x += v.x*dt + 500; pos.y += v.y*dt;}

  if (pos.x >= 500) {pos.x += v.x*dt - 500; pos.y += v.y*dt;}

  if (pos.y <= 0) {pos.x += v.x*dt; pos.y += v.y*dt + 500;}

  if (pos.y >= 500) {pos.x += v.x*dt; pos.y += v.y*dt - 500;}


  for (int i = 0; i < num_nodes; i++){

   node_x[i] = pos.x + 3 * sin(2*M_PI*i/num_nodes);
   node_y[i] = pos.y + 3 * cos(2*M_PI*i/num_nodes);
  }
}

Cell::Cell()
{
  pos = 0;
  r = 0;

  double node_x[num_nodes];
  double node_y[num_nodes];
 for (int i=0; i<num_nodes; i++){
 node_x[i] = 0;
 node_y[i] = 0;
 }

}


Cell::Cell(Vect2D pos_, double r0_,int type_, double Fasl_, int infectious_, int inf_, double age_, int card_, int diab_, int bp_, int cd_, int can_, int mob_)
{

  age = age_;
  time = 0;
  pos=pos_; r=r0_;
  type = type_;
  death_prob = 0;
  v.x = 0, v.y = 0, f.x = 0, f.y = 0, 
  card = card_, diab = diab_, bp = bp_, cd = cd_, can = can_;
  Fasl = Fasl_; MustDelete = false;
  MustDivide = false;
  inf = inf_;
  infectious = infectious_;
  newInf = 0;
  symp = 0;
  die = 0;
  indirect = 0;
  mobile = mob_;
  death_time = 0;

  incubation_time = 0;

  hospitalization_time = 0;

  release_time = 0;
  
  immune = 0;
  isolated = 0;
  if (inf == 1){
  const gsl_rng_type * T;
  gsl_rng * ro;	

  gsl_rng_env_setup();
  //gsl_rng_default_seed = (int) time(0);

  unsigned long int mySeed = rand();

  T = gsl_rng_default;
  ro = gsl_rng_alloc (T);
  gsl_rng_set(ro, mySeed);

  //incubation_time = 24*gsl_ran_lognormal(ro, 1.621, 0.418);
  incubation_time = 180.4;

  double rel_perturb = gsl_ran_gaussian(ro, 1.0);

  double hosp_perturb = gsl_ran_gaussian(ro, 1.0);




  hospitalization_time = incubation_time + 24*3.9 + hosp_perturb*24*1;

  release_time = hospitalization_time + 24*13.0 + 24*7*rel_perturb;

  gsl_rng_free (ro);}

  double r1 = (double) rand()/RAND_MAX;


  num_nodes = 16;

  node_x = new double[num_nodes];
  node_y = new double[num_nodes];

  for (int n=0; n<num_nodes; n++){



   node_x[n] = pos.x + 3 * sin(2. * M_PI * (double) n / num_nodes);
   node_y[n] = pos.y + 3 * cos(2. * M_PI * (double) n / num_nodes);





}

}





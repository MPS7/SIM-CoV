

#include "patch.h"
#include "Pltpair.h"
#include "parameters.h"
#include <vector> // vector containers
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#include <iomanip>
#include <stdio.h>
#include <time.h>       /* time */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std; // permanently use the standard namespace

bool DeleteC(Cell& cc) {return (cc.MustDelete == true);}
list <Cell> CC;

int  NPatchesx, NPatchesy;
list <patch> PatchList;
patch** Patches;



#ifdef allcells
#undef allcells
#endif
#define allcells(c) for(list <Cell>::iterator c = CC.begin(); c != CC.end();c++)




int numbAll, numb3, numb4, numb_inf, numb_infectious, numb_symp, numb_isolated, numb_dailyInf, numb_death, numb_immune, numb_dir, numb_ind;
int tot_case = 1;

double cumV = 0;




const int t_num = 2500000; // number of time steps (running from 1 to t_num)
const int t_disk = 100;  // steps to save results
const int t_info = 100;  // time steps to write info




/// *****************
/// DECLARE VECTORS
/// *****************



double **factorT;
double **factorTnew;

double **factorCa;
double **factorCanew;
/// *****************
/// DECLARE FUNCTIONS
/// *****************

// The following functions are used in the simulation code.

void initialize(); // allocate memory and initialize variables
void setCells();



void ConcentrationT();
void intraCell(int, double);


void computeTransmission();
void newArrays();
void updatePatches();
void interpolate();
double getT(Vect2D pos);
double getCa(Vect2D pos);

double delta(double);


void countCells(int);
void computeMotion(double, int);
void computeForces(gsl_rng * r);

void compute_average();


void write_particle1_vtk(int, int, int);
void write_particleInf_vtk(int);
void write_particleSymp_vtk(int);
void write_particleHosp_vtk(int);
void write_particleImm_vtk(int);
void write_T_vtk(int);
void write_Ca_vtk(int);

void write_data(int);
void write_data2(int);
/// *************
/// MAIN FUNCTION
/// *************

// This is the main function, containing the simulation initialization and the simulation loop.

int main() {

  srand(time(0));
  const gsl_rng_type * T;
  gsl_rng * r;	

  gsl_rng_env_setup();
 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  unsigned long int mySeed = rand();

  gsl_rng_set(r, mySeed);

  /// ************
  /// PREPARATIONS
  /// ************



  initialize();
  newArrays();
  setCells();
  write_data2(0);
  /// Report derived parameters

  cout << "simulation of hybrid discrete-continuous system" << endl;
  cout << "===============================================" << endl;


  cout << endl;
  cout << "starting simulation" << endl;



  for(int t = 1; t <= t_num; ++t) {

    updatePatches();
    countCells(t);

    if(numb4 == 0) break;
    interpolate();

    intraCell(t, dt);
    computeTransmission();

 //   if (numb_inf = 0) CC.push_back(Cell(Vect2D(250, 250), 0.5, 4, 0, 0, 1, 25, 0, 0, 0, 0, 0, 1));


    computeForces(r);
    computeMotion(dt, t);

    ConcentrationT(); //coagulation fac solver

    if(t % t_disk == 0) { //write only each t_disk

       // write_particle1_vtk(t , 3, numb3);
      write_particle1_vtk(t , 4, numb4);
      write_particleInf_vtk(t);
      write_particleSymp_vtk(t);
      write_particleHosp_vtk(t);
      write_particleImm_vtk(t);

      // write_T_vtk(t);
       write_Ca_vtk(t);

    }

    /// Report end of time step

    if(t % t_info == 0) {
      compute_average();
      write_data(t);
      //write_data2(t);
      cout << "completed time step " << t*dt << " in [1, " << t_num*dt << "]" << endl;
    }
}
  /// Report successful end of simulation

  cout << "simulation complete" << endl;

  return 0;

} // end of main function



void countCells(int time){
numbAll = 0, numb3 = 0, numb4 = 0, numb_inf = 0, numb_symp = 0, numb_infectious = 0, numb_isolated = 0; numb_immune = 0, numb_dailyInf = 0, numb_ind = 0;


allcells(c){
	numbAll++;
	if (c->type == 3) numb3++;
	if (c->type == 4) numb4++;

	if (c->type == 4 && c->inf == 1) numb_inf++;
	if (c->type == 4 && c->symp == 1) numb_symp++;
	if (c->type == 4 && c->infectious == 1) numb_infectious++;
	if (c->type == 4 && c->immune == 1) numb_immune++;
	if (c->type == 4 && c->isolated == 1) numb_isolated++;

	if (c->type == 4 && c->newInf == 1) {tot_case++; numb_dailyInf += 1; c->newInf = 0;}

	if (c->type == 4 && c->indirect == 1) {numb_ind += 1; c->indirect = 0;}

	numb_death = nbre_nk - numb4;
}
}



void setCells(){

int n = 0;

    srand(time(0));
while (n < nbre_nk){

    double r1 = (double) rand()/RAND_MAX;

    double r2 = (double) rand()/RAND_MAX;

    double cposx = 0;
    cposx = 10 + r1*480;
    double cposy = 0;
    cposy = 10 + r2*480;

    int indk = 1;

    allcells(c){
        double d = sqrt((cposx - c->pos.x)*(cposx - c->pos.x)+(cposy - c->pos.y)*(cposy - c->pos.y));

        if (d < 2*c->r) indk = indk*0;
    }

        if (indk == 1) {
	    double age_f = (double) rand()/RAND_MAX;
	    double ak;
	    int card1, diab1, bp1, cd1, can1;

	    if (age_f <= 0.769) ak = rand() % 42 + 18;	
            if (age_f > 0.769 && age_f <= 0.915) ak = rand() % 10 + 60;
            if (age_f > 0.915 && age_f <= 0.956) ak = rand() % 5 + 70;
            if (age_f > 0.956 && age_f <= 1) ak = rand() % 15 + 75;

	    double rdm = (double) rand()/RAND_MAX;



	    if (rdm <= 0.2) card1 = 1; 
	    else card1 = 0;

	    rdm = (double) rand()/RAND_MAX;

	    if (rdm <= 0.08) diab1 = 1;
	    else diab1 = 0;

	    rdm = (double) rand()/RAND_MAX;

	    if (rdm <= 0.25) bp1 = 1;
	    else bp1 = 0;

	    rdm = (double) rand()/RAND_MAX;

	    if (rdm <= 0.2) cd1 = 1;
	    else cd1 = 0;

	    rdm = (double) rand()/RAND_MAX;

	    if (rdm <= 0.005) can1 = 1;
	    else can1 = 0;

if (n <= nbre_nk){

	    if (n == 10) CC.push_back(Cell(Vect2D(cposx, cposy), 0.5, 4, 0, 0, 1, ak, card1, diab1, bp1, cd1, can1, 1));


else CC.push_back(Cell(Vect2D(cposx, cposy), 0.5, 4, 0, 0, 0, ak, card1, diab1, bp1, cd1, can1, 1));}

else CC.push_back(Cell(Vect2D(cposx, cposy), 0.5, 4, 0, 0, 0, ak, card1, diab1, bp1, cd1, can1, 0));
            n++;
        }



}

}


void initialize() {


  int ignore; // ignore return value of system calls
  ignore = system("mkdir -p vtk_PDE"); // create folder if not existing
  ignore = system("rm -f vtk_PDE/*.*"); //
  ignore = system("rm -f vtk_particle/*.*"); //
  ignore = system("mkdir -p vtk_particle"); // create folder if not existing
  ignore = system("rm -f data.txt"); // delete file if existing
  ignore = system("rm -f data2.txt"); // delete file if existing



  factorT=new double*[Nx];
  factorTnew=new double*[Nx];


  factorCa=new double*[Nx];
  factorCanew=new double*[Nx];

  for(int X = 0; X < Nx; ++X) {

    factorT[X]= new double[Ny];
    factorTnew[X]= new double[Ny];


    factorCa[X]= new double[Ny];
    factorCanew[X]= new double[Ny];

  }





  for(int X = 0; X < Nx; ++X) {
    for(int Y = 0; Y < Ny; ++Y) {
      factorT[X][Y] = 0;
      factorTnew[X][Y] = 0;

      factorCa[X][Y] = 0;
      factorCanew[X][Y] = 0;
      }
}

  return;
}


void newArrays(){
  NPatchesx = 20;
  NPatchesy = 20;

  Patches = new patch*[NPatchesx];
  for (int i = 0; i < NPatchesx; i++)
      Patches[i] = new patch[NPatchesy];


  for (int i = 0; i < NPatchesx; i++)
 	 for (int j = 0; j < NPatchesy; j++)
		 PatchList.push_back(patch(i,j));



	for(list <patch>::iterator p = PatchList.begin(); p != PatchList.end();p++){
		for(list <patch>::iterator p2 = PatchList.begin(); p2 != PatchList.end();p2++){
			if(p2->nx == p->nx+1) p->neighbour.push_back(&(*p2));
			if(p2->ny == p->ny+1) p->neighbour.push_back(&(*p2));
			if(p2->nx == p->nx+1) p->neighbour.push_back(&(*p2));
			if(p2->nx == p->nx-1) p->neighbour.push_back(&(*p2));
		}
	}
}


void intraCell(int t, double dt){
allcells(c){


        if (t % 24/dt == 0) if (c->type == 4 && c->inf == 0 && c->immune == 0) c->Step(dt, 0, getCa(c->pos));
	if (c->type == 4 && c->inf == 1) c->Step2(dt);
}

}





void ConcentrationT() {

    for(int X = 1; X < Nx -1 ; ++X) {
    for(int Y = 1; Y < Ny -1; ++Y) {

    factorCanew[X][Y]=factorCa[X][Y]+dt*(Dc/(SQ(h)))*(factorCa[X-1][Y]+factorCa[X+1][Y]+factorCa[X][Y+1]+factorCa[X][Y-1]-4*factorCa[X][Y]) - dt*gamma2*factorCa[X][Y];

     }
   }

    for(int X = 1; X < Nx -1 ; ++X) {
    for(int Y = 1; Y < Ny -1; ++Y) {

	     factorCa[X][Y]=factorCanew[X][Y];

     }
   }

   return;
}

void computeMotion(double dt, int t){
allcells(c){

	if (c->type == 4 && c->isolated == 0 && c->mobile == 1){
		c->Move2(dt,t);
		}

}

CC.remove_if(DeleteC);
}



void computeForces(gsl_rng * r){

  allcells(c)   {c->f = 0;
    if (c->type == 4){
        double r1 = (double) rand()/RAND_MAX;
        double r2 = (double) rand()/RAND_MAX;

        c->f.x +=  (1.34*3600*gsl_ran_gaussian(r, 0.26) - c->v.x)/(dt*0.55) + ff*(- 0.5 + r1);
        c->f.y +=  (1.34*3600*gsl_ran_gaussian(r, 0.26) - c->v.y)/(dt*0.55) + ff*(- 0.5 + r2);
    }

  }
	
}




void computeTransmission(){

	for(list <patch>::iterator pat = PatchList.begin(); pat != PatchList.end();pat++){


		for(list <Cell*>::iterator ic = pat->CellList.begin(); ic != pat->CellList.end(); ic++){
		
		Cell* c=*ic;
		if (c-> infectious == 1){


			for(list <Cell*>::iterator icp = pat->CellList.begin(); icp != pat->CellList.end(); icp++){

				Cell* cp=*icp;
				if (cp != c && cp -> immune == 0 && cp-> inf == 0){
				double d = sqrt((cp->pos.x - c->pos.x)*(cp->pos.x - c->pos.x)+(cp->pos.y - c->pos.y)*(cp->pos.y - c->pos.y));

;
	    			double dref = 1;
	    			if (d <= dref){
				 double rend = (double) rand()/RAND_MAX;
	   			 if (rend < p_d) {cp->infected(dt, 0); numb_dir++;}

 }} 

			}



			for(list <patch*>::iterator nei = pat->neighbour.begin(); nei != pat->neighbour.end();nei++){
				for(list <Cell*>::iterator icp = (*nei)->CellList.begin(); icp != (*nei)->CellList.end(); icp++){

					Cell* cp=*icp;
					if (cp != c && cp -> immune == 0 && cp-> inf == 0){
					double d = sqrt((cp->pos.x - c->pos.x)*(cp->pos.x - c->pos.x)+(cp->pos.y - c->pos.y)*(cp->pos.y - c->pos.y));


		    			double dref = 1;
		    			if (d <= dref){
					 double rend = (double) rand()/RAND_MAX;
	   				 if (rend < p_d) {cp->infected(dt,0); numb_dir++;}
					 }}
				}
			}
		}}
}
}




void updatePatches() {


	for(list <patch>::iterator p = PatchList.begin(); p != PatchList.end();p++){
		p->CellList.clear();
		allcells(c){
			if (c->pos.x >= p->nx*11 - 5 && c->pos.x < p->nx*11 + 11 - 5 && c->pos.y >= p->ny*11 - 5 && c->pos.y < p->ny*11 + 11 - 5)
			p->CellList.push_back(&(*c));
		}

	}



}

void compute_average(){



long long cumul;

allcells(c){
		cumul += sqrt(SQ(c->v.x) + SQ(c->v.y))/numb4;
	}

cout << cumul << endl;

/*
cumV = 0;

    for(int X = 1; X < Nx -1 ; ++X) {
    for(int Y = 1; Y < Ny -1; ++Y) {
		cumV += factorCa[X][Y];

    }}

*/
 
}


void interpolate(){



allcells(c){


		int kx = (int) (c->pos.x/h);
		int ky = (int) (c->pos.y/h);

		if (kx > 1 && kx < Nx-2 && ky > 1 && ky < Ny-2){
			for (int i = kx-1; i <= kx+1; i++)
				for (int j = ky-1; j <= ky+1; j ++){
					double r = sqrt(SQ(c->pos.x - i*h) + SQ(c->pos.y - j*h));

					if (c->type == 4 && c->infectious == 1) factorCa[i][j] += delta(r)*beta2*dt;


        }}



	}



}


double getT(Vect2D pos){

	double fac;
	//double rep = 1/h;
        int kx = (int) (pos.x/h);
	int ky = (int) (pos.y/h);


	fac = (((ky+1)*h - pos.y)/h)*(((kx+1)*h-pos.x)*factorT[kx][ky] + (pos.x-kx*h)*factorT[kx+1][ky])/h + ((pos.y - ky*h)/h)*(((kx+1)*h - pos.x)*factorT[kx][ky+1]/h + (pos.x-kx*h)*factorT[kx+1][ky+1]/h);



	return fac;

}


double getCa(Vect2D pos){

	double fac;
	//double rep = 1/h;
        int kx = (int) (pos.x/h);
	int ky = (int) (pos.y/h);

	if (pos.x > 40 && pos.x < 400 && pos.y > 40 && pos.y < 400){
	fac = (((ky+1)*h - pos.y)/h)*(((kx+1)*h-pos.x)*factorCa[kx][ky] + (pos.x-kx*h)*factorCa[kx+1][ky])/h + ((pos.y - ky*h)/h)*(((kx+1)*h - pos.x)*factorCa[kx][ky+1]/h + (pos.x-kx*h)*factorCa[kx+1][ky+1]/h);

	}

	return fac;

}


/// **************************************
/// WRITE PDE & PARTICLE STATE TO VTK FILE
/// **************************************




void write_T_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_PDE/T_t" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "ConcentrationT_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET RECTILINEAR_GRID\n";
  output_file << "DIMENSIONS " << Nx << " " << Ny << " 1" << "\n";
  output_file << "X_COORDINATES " << Nx << " float\n";

  for(int X = 0; X < Nx; ++X) {
    output_file << X*h  << " ";
  }

  output_file << "\n";
  output_file << "Y_COORDINATES " << Ny << " float\n";

  for(int Y = 0; Y < Ny; ++Y) {
    output_file << Y*h  << " ";
  }

  output_file << "\n";
  output_file << "Z_COORDINATES " << 1 << " float\n";
  output_file << 0 << "\n";
  output_file << "POINT_DATA " << (Nx) * (Ny) << "\n";


  /// Write T

  output_file << "SCALARS T float\n";
  output_file << "LOOKUP_TABLE default\n";

  for(int Y = 0; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << fixed << setprecision(5) << factorT[X][Y] << "\n";
    }
  }
 /// Close file

  output_file.close();

  return;
}




void write_Ca_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_PDE/Ca_t" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "ConcentrationT_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET RECTILINEAR_GRID\n";
  output_file << "DIMENSIONS " << Nx << " " << Ny << " 1" << "\n";
  output_file << "X_COORDINATES " << Nx << " float\n";

  for(int X = 0; X < Nx; ++X) {
    output_file << X*h  << " ";
  }

  output_file << "\n";
  output_file << "Y_COORDINATES " << Ny << " float\n";

  for(int Y = 0; Y < Ny; ++Y) {
    output_file << Y*h  << " ";
  }

  output_file << "\n";
  output_file << "Z_COORDINATES " << 1 << " float\n";
  output_file << 0 << "\n";
  output_file << "POINT_DATA " << (Nx) * (Ny) << "\n";


  /// Write T

  output_file << "SCALARS Ca float\n";
  output_file << "LOOKUP_TABLE default\n";

  for(int Y = 0; Y < Ny; ++Y) {
    for(int X = 0; X < Nx; ++X) {
      output_file << fixed << setprecision(5) << factorCa[X][Y] << "\n";
    }
  }
 /// Close file

  output_file.close();

  return;
}



double delta(double r){
	if (r <= 10) return 0.25*(1+ cos(M_PI*r/10));
	else return 0;

}


void write_data(int time) {

  /// Create filename

  string output_filename("data.txt");
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.c_str(), fstream::app);

  /// Write data

  output_file << dt*time << " "; // 
  output_file << numb_inf << " "; // 
  output_file << numb_symp << " "; // 
  output_file << numb_isolated << " "; // 
  output_file << numb_death << " "; // 
  output_file << numb_immune << " "; // t
  output_file << tot_case << " "; // 
  output_file << cumV << " "; // t
  output_file << numb_dailyInf << " "; // *
  output_file << numb_dir << " "; // *
  output_file << numb_ind << " "; // *
  output_file << numb4 << "\n"; // *

  /// Close file

  output_file.close();

  return;
}



void write_data2(int time) {

  /// Create filename

  string output_filename("data2.txt");
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.c_str(), fstream::app);

  /// Write data
allcells(c){
  output_file << c->age << " "; // 
  output_file << c->card << " "; // 
  output_file << c->diab << " "; // 
  output_file << c->bp << " "; // 
  output_file << c->cd << " "; // 
  output_file << c->can << " "; // 
  output_file << "\n"; // 
}

  /// Close file

  output_file.close();

  return;
}


void write_particle1_vtk(int time, int type_, int numbs) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_particle/cell_t" << type_ << "_" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << 16*numbs << " float\n";
allcells(c){
if (c->type == type_){
  for(int n = 0; n <  c->num_nodes; ++n) {
   output_file << c->node_x[n] << " " << c->node_y[n] << " 0.0001\n";
  }}
}


  output_file << "POLYGONS " << numbs << " " << numbs*(16+1) << "\n";
for (int j=0; j<numbs; j++){
  output_file << 16 << " ";

  for(int n = 0; n < 16; ++n) {
    output_file << j*16 + n << " ";
  }
    output_file << "\n";
}
  /// Close file

  output_file.close();

  return;
}




void write_particleInf_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_particle/cell_t_Inf" << "_" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << 16*numb_inf << " float\n";
allcells(c){
if (c->inf == true){
  for(int n = 0; n <  c->num_nodes; ++n) {
   output_file << c->node_x[n] << " " << c->node_y[n] << " 0.001\n";
  }}
}


  output_file << "POLYGONS " << numb_inf << " " << numb_inf*(16+1) << "\n";
for (int j=0; j<numb_inf; j++){
  output_file << 16 << " ";

  for(int n = 0; n < 16; ++n) {
    output_file << j*16 + n << " ";
  }
    output_file << "\n";
}
  /// Close file

  output_file.close();

  return;
}


void write_particleSymp_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_particle/cell_t_symp" << "_" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << 16*numb_symp << " float\n";
allcells(c){
if (c->symp == true){
  for(int n = 0; n <  c->num_nodes; ++n) {
   output_file << c->node_x[n] << " " << c->node_y[n] << " 0.003\n";
  }}
}


  output_file << "POLYGONS " << numb_symp << " " << numb_symp*(16+1) << "\n";
for (int j=0; j<numb_symp; j++){
  output_file << 16 << " ";

  for(int n = 0; n < 16; ++n) {
    output_file << j*16 + n << " ";
  }
    output_file << "\n";
}
  /// Close file

  output_file.close();

  return;
}




void write_particleHosp_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_particle/cell_t_hosp" << "_" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << 16*numb_isolated << " float\n";
allcells(c){
if (c->isolated == true){
  for(int n = 0; n <  c->num_nodes; ++n) {
   output_file << c->node_x[n] << " " << c->node_y[n] << " 0.005\n";
  }}
}


  output_file << "POLYGONS " << numb_isolated  << " " << numb_isolated *(16+1) << "\n";
for (int j=0; j<numb_isolated; j++){
  output_file << 16 << " ";

  for(int n = 0; n < 16; ++n) {
    output_file << j*16 + n << " ";
  }
    output_file << "\n";
}
  /// Close file

  output_file.close();

  return;
}


void write_particleImm_vtk(int time) {

  /// Create filename

  stringstream output_filename;
  output_filename << "vtk_particle/cell_t_imm" << "_" << time << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << 16*numb_immune << " float\n";
allcells(c){
if (c->immune == true){
  for(int n = 0; n <  c->num_nodes; ++n) {
   output_file << c->node_x[n] << " " << c->node_y[n] << " 0.01\n";
  }}
}


  output_file << "POLYGONS " << numb_immune  << " " << numb_immune *(16+1) << "\n";
for (int j=0; j<numb_immune; j++){
  output_file << 16 << " ";

  for(int n = 0; n < 16; ++n) {
    output_file << j*16 + n << " ";
  }
    output_file << "\n";
}
  /// Close file

  output_file.close();

  return;
}






/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Fri Sep. 21 2012
 * 
 *      
 */


#ifndef __INITIAL_SETTING_H__
#define __INITIAL_SETTING_H__


#include <cstring>
#include <list>
#include <vector>

using namespace std;



struct Time {
  double start;
  double end;
  double dt;
};


struct Domain {
  bool use;
  double xLower;
  double xUpper;
  int xGridNumber;
  double yLower;
  double yUpper;
  int yGridNumber;
  double zLower;
  double zUpper;
  int zGridNumber;
};


struct Visualization {
  bool use;
  int precision;
  double timeinterval;
  string classname;
  string filename;
  Domain domain;
};


struct FieldSetting {
  string solver;
  string epsilon;
  string mu;
  string initial_distribution_ex;
  string initial_distribution_ey;
  string initial_distribution_ez;
  string initial_distribution_hx;
  string initial_distribution_hy;
  string initial_distribution_hz;
  Visualization visualization;
};


/*
struct ParticleGrid {
  string mode;
  double dx;
  double dy;
  double dz;
};
*/

/*
struct UserGrid {
  double dx;
  double dy;
  double dz;
};
*/



struct PhysicalData {
  string charge;
  string mass;
  string velocity_u;
  string velocity_v;
  string velocity_w;
};

// UserGrid must be coarse than ParticleGrid.
struct UserSpecificSetting {
  bool use;
  string mode;
  int group;
  string classname;
  double dx;
  //Domain domain;
  //ParticleGrid particleGrid;
  //UserGrid userGrid;

  PhysicalData physical_data;
};

// Direct initial particle setting.
struct InitialParticle {
  string x;
  string y;
  string z;
  
  PhysicalData physical_data;
};



//This defines a subregion we will be filling with particles. 
struct LevelFunction {
  string function;
  char range;
  bool equality;
};


//several level set functions will comprise the region we will fill with particles
struct Bunch {
  int group;
  bool use;
  string mode;
  Domain domain;
  vector<LevelFunction> functions;
  PhysicalData physical_data;
};


struct NormalRandomBunch {
  bool use;
  int group;
  size_t num_particles;
  Domain domain;
  double mu_x;
  double mu_y;
  double mu_z;
  double sigma_x;
  double sigma_y;
  double sigma_z;
  PhysicalData physical_data;
};

struct UniformRandomBunch {
  bool use;
  int group;
  size_t num_particles;
  Domain domain;
  vector<LevelFunction> functions;
  PhysicalData physical_data;
};


struct ParticleInitialDistribution {
  size_t buffer_size;
  int directSetting_group;
  list<InitialParticle> directSetting;
  UserSpecificSetting userSetting;
  vector<Bunch> bunchs;
  vector<NormalRandomBunch> normal_bunchs;
  vector<UniformRandomBunch> uniform_bunchs;
};


struct ParticleSetting {
  // gamma : Lorentz factor (relativistic factor)
  //double gamma;
  ParticleInitialDistribution distribution;
  Visualization visualization;
};



//This is a superstructure which encompasses everything to do with initial setting of the problem
//Things might need to be added to this. or else we can do without this and directly start with the Iniitial Solver Setting.
struct InitialSetting {
  int dim;
  size_t omp_num_threads;
  Time time;
  Domain domain;
  string controller;
  FieldSetting field;
  ParticleSetting particle;
  //Logger logger;
};



int loadInitialSetting(const char* filename, InitialSetting* is);



size_t getInitialHexagonalPacking3D(double localXXmin, double localYYmin, double localZZmin, double localXXmax, double localYYmax, double localZZmax, double dz, std::vector<double>& coord_x, std::vector<double>& coord_y, std::vector<double>& coord_z);

size_t getInitialRectangularPacking3D(double localXXmin, double localYYmin, double localZZmin, double localXXmax, double localYYmax, double localZZmax, double dz, std::vector<double>& coord_x, std::vector<double>& coord_y, std::vector<double>& coord_z);

size_t getInitialHexagonalPacking2D(double localXXmin, double localZZmin, double localXXmax, double localZZmax, double dx, std::vector<double>& coord_x, std::vector<double>& coord_z);

size_t getInitialRectangularPacking2D(double localXXmin, double localZZmin, double localXXmax, double localZZmax, double dx, std::vector<double>& coord_x, std::vector<double>& coord_z);
  









#endif

/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jun. 05 2012
 * 
 *      
 */



#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

#include <parser.h>

#include "time_controller.h"
#include "parallel_solver.h"
#include "serial_solver.h"
#include "resultview.h"
#include "em_error.h"
#include "initial_setting.h"
#include "constants.h"
#include "serial_mover.h"
#include "user_specific_initial_particle_setting.h"


#ifndef _MSC_VER  
#include <getopt.h>
#include <unistd.h>
#endif 


using namespace mparser;
using namespace std;



#ifndef _MSC_VER  
void parseParameters(int argc, char *argv[], double& startTime, double& endTime, double& dt, double& leftX, double& rightX, double& dx, double& leftY, double& rightY, double& dy, double& leftZ, double& rightZ, double& dz, int& vis_interval, double& vis_leftX, double& vis_rightX, double& vis_leftY, double& vis_rightY, double& vis_leftZ, double& vis_rightZ, int& errormode);
#endif 

int validateInitialSetting(InitialSetting* is);
int generateInitialParticles(ParticleInitialDistribution& pd, ParticleMover* pm);
bool isIncludedBunch(Parser& parser, double x, double y, double z, const Bunch& bunch);
bool isIncludedUniformBunch(Parser& parser, double x, double y, double z, const UniformRandomBunch& bunch);
double calculateLorentzFactor(double u, double v, double w);

int main(int argc, char *argv[]) {

/*  
  double startTime = 0;
  double endTime = 0.00000001;
  //double dt =      0.000000000001;
  double dt = 75; // 75%
  double leftX = 0;
  double rightX = 0.256;
  double leftY = 0;
  double rightY = 0.256;
  double leftZ = 0;
  double rightZ = 1.024;
  double dx = 0.008;
  double dy = 0.008;
  double dz = 0.008;
  int vis_interval = 1;
  double vis_leftX = leftX;
  double vis_rightX = rightX;
  double vis_leftY = leftY;
  double vis_rightY = rightY;
  double vis_leftZ = leftZ;
  double vis_rightZ = rightZ;
  int errormode = 0;

#ifndef _MSC_VER
  parseParameters(argc, argv, startTime, endTime, dt, leftX, rightX, dx, leftY, rightY, dy, leftZ, rightZ, dz, vis_interval, vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ, errormode);
#endif 
*/

  //cout << "OMP_NUM_THREAD : " << omp_get_max_threads() << endl;


  double startTime;
  double endTime;
  double dt;
  double leftX;
  double rightX;
  double leftY;
  double rightY;
  double leftZ;
  double rightZ;
  double dx;
  double dy;
  double dz;
  int field_vis_interval;
  int particle_vis_interval;
  double vis_leftX;
  double vis_rightX;
  double vis_leftY;
  double vis_rightY;
  double vis_leftZ;
  double vis_rightZ;
  int errormode = 0;

  int error;



  if(argv[1] == NULL) {
    cout << "Input initial setting file error." << endl;
    return EM_ERR_SYS_INITIAL_SETTING_GEN;
  }


  InitialSetting is;
  //cout << "is address : " << &is << endl;
  error = loadInitialSetting(argv[1], &is);
  if(error) return error;


  error = validateInitialSetting(&is);
  if(error) return error;

#ifdef _OPENMP
  omp_set_num_threads(static_cast<int>(is.omp_num_threads));
  cout << "OMP_NUM_THREAD : " << omp_get_max_threads() << endl;
#endif

  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : Converting data from InitialSetting.
  ////////////////////////////////////////////////////////////////////////////////////////
  startTime = is.time.start;
  endTime = is.time.end;

  leftX = is.domain.xLower;
  rightX = is.domain.xUpper;
  leftY = is.domain.yLower;
  rightY = is.domain.yUpper;
  leftZ = is.domain.zLower;
  rightZ = is.domain.zUpper;


  dx = (rightX - leftX) / is.domain.xGridNumber;
  dy = (rightY - leftY) / is.domain.yGridNumber;
  dz = (rightZ - leftZ) / is.domain.zGridNumber;
  
  dt = is.time.dt;
  if( (dt <= 0.0) || (dt >= 100.0) ) dt = 75.0;
  double ndt = dt / 100 / C0 / sqrt((1/dx/dx) + (1/dy/dy) + (1/dz/dz));

  field_vis_interval = static_cast<int>((is.field.visualization.timeinterval)/ndt);
  cout << "field_vis_interval : " << field_vis_interval << endl;

  particle_vis_interval = static_cast<int>((is.particle.visualization.timeinterval)/ndt);
  cout << "particle_vis_interval : " << particle_vis_interval << endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : Converting data from InitialSetting.
  ////////////////////////////////////////////////////////////////////////////////////////




  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : FieldSolver Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////
  FieldSolver* solver = 0;
  if(is.field.solver == "SerialFieldSolver") {
    solver = new SerialFieldSolver(ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "SerialTestSolver") {
    solver = new SerialTestSolver(ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "SerialNullFieldSolver") {
    solver = new SerialNullFieldSolver(ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "SerialTestSolverWithParticle1") {
    solver = new SerialTestSolverWithParticle1(ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "SerialTestSolverWithParticle2") {
    solver = new SerialTestSolverWithParticle2(ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "SerialTestSolverWithParticle3") {
    solver = new SerialTestSolverWithParticle3(ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else {
    // Error
    EM_ERROR(STR_ERR_SYS_OBJECT_CREATION_FIELDSOLVER);
    return EM_ERR_SYS_OBJECT_CREATION_FIELDSOLVER;
  }

  solver->initializeSolver();
  solver->setErrorMode(errormode);
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : FieldSolver Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////



/*
#ifndef _MSC_VER
  FieldViewer* fieldViewer = new VTKFieldViewer(solver, "output/EMField", vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ, VTKFieldViewer::VECTOR_SEPARATED);
#else
  FieldViewer* fieldViewer = new VTKFieldViewer(solver, "output\\EMField", vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ, VTKFieldViewer::VECTOR_SEPARATED);
#endif 
*/


  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : FieldViewer Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////
  FieldViewer* fieldViewer = 0;
  if(is.field.visualization.use) {
    vis_leftX = is.field.visualization.domain.xLower;
    vis_rightX = is.field.visualization.domain.xUpper;
    vis_leftY = is.field.visualization.domain.yLower;
    vis_rightY = is.field.visualization.domain.yUpper;
    vis_leftZ = is.field.visualization.domain.zLower;
    vis_rightZ = is.field.visualization.domain.zUpper;

    if(is.field.visualization.classname == "VTKFieldViewer") {
      fieldViewer = new VTKFieldViewer(solver, is.field.visualization.filename, vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ, is.field.visualization.precision, VTKFieldViewer::VECTOR_SEPARATED);
    } else {
      EM_ERROR(STR_ERR_SYS_OBJECT_CREATION_FIELDVIEWER);
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : FieldViewer Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////





  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : ParticleMover Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////

  ParticleMover* mover = new SerialParticleMover(dynamic_cast<SerialFieldSolver*>(solver), ndt);

  ////////////////////////////////////////////////////////////////////////////////////////
  // End : ParticleMover Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////

  



  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : Initial Particle Setting.
  ////////////////////////////////////////////////////////////////////////////////////////
  double startParticleInitialization = omp_get_wtime();
  error = generateInitialParticles(is.particle.distribution , mover);
  double endParticleInitialization = omp_get_wtime();
  if(error) {
    EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
    return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  }
  cout << "Initial Particle Generating Time is " << endParticleInitialization - startParticleInitialization << " seconds." << endl;

  mover->calculateInitialVelocity();

  double start = omp_get_wtime();
  mover->initElectricIntensityByParticles();
  double end = omp_get_wtime();
  cout << "Initializing Electric Field by Particles Time is " << end - start << " seconds." << endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : Initial Particle Setting.
  ////////////////////////////////////////////////////////////////////////////////////////



/*
#ifndef _MSC_VER
  ParticleViewer* particleViewer = new VTKParticleViewer(mover, "output/Particles", vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ);
#else
  ParticleViewer* particleViewer = new VTKParticleViewer(mover, "output\\Particles", vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ);
#endif
*/


  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : ParticleViewer Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////
  ParticleViewer* particleViewer = 0;
  if(is.particle.visualization.use) {
    if(is.particle.visualization.classname == "VTKParticleViewer") {
      // Although this creator designate visualization domain, this domain is meaningless. This domain setting is used for later development.
      particleViewer = new VTKParticleViewer(mover, is.particle.visualization.filename, leftX, rightX, leftY, rightY, leftZ, rightZ, is.particle.visualization.precision);
    } else {
      EM_ERROR(STR_ERR_SYS_OBJECT_CREATION_PARTICLEVIEWER);
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : ParticleViewer Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////




  TimeController controller(startTime, endTime, ndt, solver, mover, fieldViewer, particleViewer, field_vis_interval, particle_vis_interval);
  //TimeController controller(startTime, endTime, ndt, solver, mover, 0, particleViewer, vis_interval);
  //TimeController controller(startTime, endTime, ndt, solver, fieldViewer, vis_interval);
  //TimeController controller(startTime, endTime, dt, solver);


  map< int, list<Particle> >& particlesGroup = mover->getParticlesGroup();
  size_t totalnumparticle = 0;
  for(map< int, list<Particle> >::iterator it2 = particlesGroup.begin() ; it2 != particlesGroup.end() ; ++it2) {
    totalnumparticle += it2->second.size();
  }

  cout << "===============================================================================" << endl;
  cout << "leftX : " << leftX << ", rightX : " << rightX << ", dx : " << dx << ", Grid Size X : " << is.domain.xGridNumber << endl;
  cout << "leftY : " << leftY << ", rightY : " << rightY << ", dy : " << dy << ", Grid Size Y : " << is.domain.yGridNumber << endl;
  cout << "leftZ : " << leftZ << ", rightZ : " << rightZ << ", dz : " << dz << ", Grid Size Z : " << is.domain.zGridNumber << endl;
  cout << "dt : " << ndt << ", last time index : " << static_cast<int>(((endTime - startTime)/ndt) + ndt) << endl;
  cout << "The total number of particles : " << totalnumparticle << endl;
  cout << "===============================================================================" << endl;

  

  if(controller.solve()) {
    EM_ERROR("Controller solve() Error");

    if(fieldViewer != 0) delete fieldViewer;
    if(particleViewer != 0) delete particleViewer;

    delete solver;
    delete mover;

    return 0;
  }



  /*
  vector<double>* currentEx;
  vector<double>* currentEy;
  vector<double>* currentEz;
  vector<double>* currentHx;
  vector<double>* currentHy;
  vector<double>* currentHz;


  ofstream num_sol_ex("num_sol_ex.txt");
  ofstream num_sol_ey("num_sol_ey.txt");
  ofstream num_sol_ez("num_sol_ez.txt");
  ofstream num_sol_hx("num_sol_hx.txt");
  ofstream num_sol_hy("num_sol_hy.txt");
  ofstream num_sol_hz("num_sol_hz.txt");


  solver->getCurrentEx(&currentEx);
  solver->getCurrentEy(&currentEy);
  solver->getCurrentEz(&currentEz);
  solver->getCurrentHx(&currentHx);
  solver->getCurrentHy(&currentHy);
  solver->getCurrentHz(&currentHz);

  size_t size = currentEx->size();
  for(size_t i=0 ; i<size-1 ; i++) {
    num_sol_ex << (*currentEx)[i] << ",";
    num_sol_ey << (*currentEy)[i] << ",";
    num_sol_ez << (*currentEz)[i] << ",";
    num_sol_hx << (*currentHx)[i] << ",";
    num_sol_hy << (*currentHy)[i] << ",";
    num_sol_hz << (*currentHz)[i] << ",";
  }
  num_sol_ex << (*currentEx)[size-1] << endl;
  num_sol_ey << (*currentEy)[size-1] << endl;
  num_sol_ez << (*currentEz)[size-1] << endl;
  num_sol_hx << (*currentHx)[size-1] << endl;
  num_sol_hy << (*currentHy)[size-1] << endl;
  num_sol_hz << (*currentHz)[size-1] << endl;

  num_sol_ex.close();
  num_sol_ey.close();
  num_sol_ez.close();
  num_sol_hx.close();
  num_sol_hy.close();
  num_sol_hz.close();
  */

  if(fieldViewer != 0) delete fieldViewer;
  if(particleViewer != 0) delete particleViewer;

  delete solver;
  delete mover;

  return 0;
}





int validateInitialSetting(InitialSetting* is) {

  // Domain Check.
  if(is->domain.xLower >= is->domain.xUpper) {
    cout << "<x lower=\"" << is->domain.xLower << "\" upper=\"" << is->domain.xUpper << "\" gridnumber=\"" << is->domain.xGridNumber << "\" />" << endl;
    EM_ERROR(STR_ERR_COM_DOM_SCOPE_X);
    return EM_ERR_COM_DOM_SCOPE_X;
  }
  if(is->domain.yLower >= is->domain.yUpper) {
    cout << "<y lower=\"" << is->domain.yLower << "\" upper=\"" << is->domain.yUpper << "\" gridnumber=\"" << is->domain.yGridNumber << "\" />" << endl;
    EM_ERROR(STR_ERR_COM_DOM_SCOPE_Y);
    return EM_ERR_COM_DOM_SCOPE_Y;
  }
  if(is->domain.zLower >= is->domain.zUpper) {
    cout << "<z lower=\"" << is->domain.zLower << "\" upper=\"" << is->domain.zUpper << "\" gridnumber=\"" << is->domain.zGridNumber << "\" />" << endl;
    EM_ERROR(STR_ERR_COM_DOM_SCOPE_Z);
    return EM_ERR_COM_DOM_SCOPE_Z;
  }

  if(is->omp_num_threads < 0) {
    cout << "Wrong OMP_NUM_THREADS setting : " << is->omp_num_threads << " : Check <system> node."  << endl;
    EM_ERROR(STR_ERR_SYS_GEN);
    return EM_ERR_SYS_GEN;
  }


  return EM_SUCCESS;
}





int generateInitialParticle_DirectSetting(list<InitialParticle>& directSetting, ParticleMover* pm, int group) {

  if(!directSetting.size()) return EM_SUCCESS;

  list<Particle>* particles = pm->getParticles(group);
  if(particles == 0) {
    list<Particle> p;
    pm->getParticlesGroup().insert(pair< int, list<Particle> >(group, p));
    particles = pm->getParticles(group);
  }

  try {
    Parser parser;

    parser.add_variable("x");
    parser.add_variable("y");
    parser.add_variable("z");

    list<InitialParticle> p_list = directSetting;
    for(list<InitialParticle>::iterator iter = p_list.begin() ; iter!=p_list.end() ; iter++) {
      double x = parser.parse((*iter).x.c_str());
      double y = parser.parse((*iter).y.c_str());
      double z = parser.parse((*iter).z.c_str());

      parser.set_variable(0, x);
      parser.set_variable(1, y);
      parser.set_variable(2, z);

      double u = parser.parse((*iter).physical_data.velocity_u.c_str());
      double v = parser.parse((*iter).physical_data.velocity_v.c_str());
      double w = parser.parse((*iter).physical_data.velocity_w.c_str());
      double charge = parser.parse((*iter).physical_data.charge.c_str());
      double mass = parser.parse((*iter).physical_data.mass.c_str());
      
      mass *= calculateLorentzFactor(u,v,w);

      Particle p(x, y, z, u, v, w, charge, mass);
      particles->push_back(p);
    }
  } catch(...) {
    EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
    return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  }

  return EM_SUCCESS;

}

int generateInitialParticle_Analytic(vector<Bunch>& bunchs, ParticleMover* pm) {


  //double start,stop;
  list<Particle>* particles;
  int group;

  for(size_t k=0 ; k<bunchs.size() ; k++) {

    if(!bunchs[k].use) continue;
    
    group = bunchs[k].group;
    particles = pm->getParticles(group);
    if(particles == 0) {
      list<Particle> p;
      pm->getParticlesGroup().insert(pair< int, list<Particle> >(group, p));
      particles = pm->getParticles(group);
    }


    double range_x_min = bunchs[k].domain.xLower;
    double range_x_max = bunchs[k].domain.xUpper;
    double range_y_min = bunchs[k].domain.yLower;
    double range_y_max = bunchs[k].domain.yUpper;
    double range_z_min = bunchs[k].domain.zLower;
    double range_z_max = bunchs[k].domain.zUpper;



    //double dx = (range_x_max - range_x_min) / bunchs[k].domain.xGridNumber;
    //double dy = (range_y_max - range_y_min) / bunchs[k].domain.yGridNumber;
    double dz = (range_z_max - range_z_min) / bunchs[k].domain.zGridNumber;



    //size_t max_x_step = (size_t)((range_x_max - range_x_min)/dx + 1);
    //size_t max_y_step = (size_t)((range_y_max - range_y_min)/dy + 1);
    //size_t max_z_step = (size_t)((range_z_max - range_z_min)/dz + 1);
    //cout << "max_x_step = " << max_x_step <<endl;
    //cout << "max_y_step = " << max_y_step <<endl;
    //cout << "max_z_step = " << max_z_step <<endl;


    //char math_expression[255] = "x^2/3 + y^2/4 + z^2/5 - 1";
    //char math_expression[255] = "x^2 + y^2 - z";
    //char range = '<';
    //bool range_equality = false;



    Parser parser;
    
    parser.add_variable("x");
    parser.add_variable("y");
    parser.add_variable("z");


    /*
    vector<double> particle_x;
    vector<double> particle_y;
    vector<double> particle_z;
    vector<double> density;
    vector<double> mass;
    vector<double> velocity_u;
    vector<double> velocity_v;
    vector<double> velocity_w;
    vector<double> pressure;
    vector<double> thermal_energy;
    */

    //size_t num_particles = 0;



    vector<double> initPackingX;
    vector<double> initPackingY;
    vector<double> initPackingZ;
    size_t num_particle_packing = 0;


    if(bunchs[k].mode.compare("rect") == 0) {
      num_particle_packing = getInitialRectangularPacking3D(range_x_min, range_y_min, range_z_min, range_x_max, range_y_max, range_z_max, dz, initPackingX, initPackingY, initPackingZ);
    } else { // defalut : "hex"
      num_particle_packing = getInitialHexagonalPacking3D(range_x_min, range_y_min, range_z_min, range_x_max, range_y_max, range_z_max, dz, initPackingX, initPackingY, initPackingZ);
    }



    try {
      double x;
      double y;
      double z;

      double u;
      double v;
      double w;

      double mass;
      double charge;
      //double gamma; // Lorentz factor

      bool diagnose;
      
      //size_t total_num_particles = 0;

      for(size_t i=0 ; i<num_particle_packing ; ++i) {

        x = initPackingX[i];
        y = initPackingY[i];
        z = initPackingZ[i];
      
        parser.set_variable(0, x);
        parser.set_variable(1, y);
        parser.set_variable(2, z);

        diagnose = isIncludedBunch(parser, x,y,z, bunchs[k]);

        if(diagnose) {
          //total_num_particles++;
          //particle_x.push_back(x);
          //particle_y.push_back(y);
          //particle_z.push_back(z);
          u = parser.parse(bunchs[k].physical_data.velocity_u.c_str());
          //velocity_u.push_back(answer);
          v = parser.parse(bunchs[k].physical_data.velocity_v.c_str());
          //velocity_v.push_back(answer);
          w = parser.parse(bunchs[k].physical_data.velocity_w.c_str());
          //velocity_w.push_back(answer);

          charge = parser.parse(bunchs[k].physical_data.charge.c_str());
          //charge.push_back(answer);
          mass = parser.parse(bunchs[k].physical_data.mass.c_str());
          //mass.push_back(answer);

          mass *= calculateLorentzFactor(u,v,w);

          Particle p(x, y, z, u, v, w, charge, mass);
          particles->push_back(p);
        

          /*
          cout << "x=" << x <<", y=" << y << ", z=" << z << endl;
          cout << "density=" << density.at(num_particles-1) << ", mass=" << mass.at(num_particles-1) << endl;
          cout << "u=" << velocity_u.at(num_particles-1) <<", v=" << velocity_v.at(num_particles-1) << ", w=" << velocity_w.at(num_particles-1) << endl;
          cout << "pressure=" << pressure.at(num_particles-1) << ", thermal_energy=" << thermal_energy.at(num_particles-1) << endl;
          */
          //break;
        }



      } // for i
      
      //EM_LOG("Generated particles number by the analytic function : " + total_num_particles);
      //cout << "Generated particles number by the analytic function : " << total_num_particles << endl;
      

    } catch(Error err) {
      EM_ERROR(err.get_msg());
      EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
      //printf("\tError: Unknown error occured in parser\n");

      return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
    }


    //cout<<"\nSetting initial data by the given level functions took "<< (stop-start)/(double) CLOCKS_PER_SEC << "  seconds" << endl;

    //cout<<"The number of particles generated by the level functions is "<< particle_x.size() << endl;


    /*
    string fn_data_x = "generated_particles_x.txt";
    string fn_data_y = "generated_particles_y.txt";
    string fn_data_z = "generated_particles_z.txt";
    ofstream fs_data_x(fn_data_x.c_str());
    ofstream fs_data_y(fn_data_y.c_str());
    ofstream fs_data_z(fn_data_z.c_str());
    */
    /*
    string fn_initial_data = "generated_particles.txt";
    ofstream fs_data(fn_initial_data.c_str());

    for(int i=0; i<particle_x.size() ; i++) {
      //fs_data_x << particle_x[i] << ",";
      //fs_data_y << particle_y[i] << ",";
      //fs_data_z << particle_z[i] << ",";
      //fs_data_x << particle_x[i] << endl;
      //fs_data_y << particle_y[i] << endl;
      //fs_data_z << particle_z[i] << endl;
      //fs_data << particle_x[i] << "\t" << particle_y[i] << "\t" << particle_z[i] <<endl;
      fs_data << particle_x[i] << "," << particle_y[i] << "," << particle_z[i] <<endl;
    }
    //fs_data_x << endl; // endl include flushing.
    //fs_data_y << endl; // endl include flushing.
    //fs_data_z << endl; // endl include flushing.
    //fs_data_x.close();
    //fs_data_y.close();
    //fs_data_z.close();
    fs_data.close();
    */

  }

  return EM_SUCCESS;

}


int generateInitialParticle_NormalRandom(vector<NormalRandomBunch>& normal_bunchs, ParticleMover* pm) {


  //double start,stop;
  int group;
  list<Particle>* particles;

  for(size_t k=0 ; k<normal_bunchs.size() ; k++) {

    if(!normal_bunchs[k].use) continue;

    group = normal_bunchs[k].group;
    particles = pm->getParticles(group);
    if(particles == 0) {
      list<Particle> p;
      pm->getParticlesGroup().insert(pair< int, list<Particle> >(group, p));
      particles = pm->getParticles(group);
    }

    double range_x_min = normal_bunchs[k].domain.xLower;
    double range_x_max = normal_bunchs[k].domain.xUpper;
    double range_y_min = normal_bunchs[k].domain.yLower;
    double range_y_max = normal_bunchs[k].domain.yUpper;
    double range_z_min = normal_bunchs[k].domain.zLower;
    double range_z_max = normal_bunchs[k].domain.zUpper;

    double u,v,w,charge,mass;
    double mu_x = normal_bunchs[k].mu_x;
    double mu_y = normal_bunchs[k].mu_y;
    double mu_z = normal_bunchs[k].mu_z;
    double sigma_x = normal_bunchs[k].sigma_x;
    double sigma_y = normal_bunchs[k].sigma_y;
    double sigma_z = normal_bunchs[k].sigma_z;

    Parser parser;
    
    try {

      u = parser.parse(normal_bunchs[k].physical_data.velocity_u.c_str());
      v = parser.parse(normal_bunchs[k].physical_data.velocity_v.c_str());
      w = parser.parse(normal_bunchs[k].physical_data.velocity_w.c_str());
      charge = parser.parse(normal_bunchs[k].physical_data.charge.c_str());
      mass = parser.parse(normal_bunchs[k].physical_data.mass.c_str());

      mass *= calculateLorentzFactor(u,v,w);

    } catch(Error err) {
      EM_ERROR(err.get_msg());
      EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
      //printf("\tError: Unknown error occured in parser\n");

      return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
    }

    

    for(size_t i=0 ; i<normal_bunchs[k].num_particles ; ++i) {
      double r1; 
      double r2;
      double x,y,z;

      // r1, r2 : uniform distribution, (0..1]
      r1 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
      r2 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
      x = mu_x + sigma_x*(sqrt(-2*log(r1)) * cos(2*mparser::PI*r2));
      if((x < range_x_min) || (x > range_x_max)) continue;

      r1 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
      r2 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
      y = mu_y + sigma_y*(sqrt(-2*log(r1)) * cos(2*mparser::PI*r2));
      if((y < range_y_min) || (y > range_y_max)) continue;

      r1 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
      r2 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
      z = mu_z + sigma_z*(sqrt(-2*log(r1)) * cos(2*mparser::PI*r2));
      if((z < range_z_min) || (z > range_z_max)) continue;

      Particle p(x, y, z, u, v, w, charge, mass);
      particles->push_back(p);
    }


  }

  return EM_SUCCESS;

}


int generateInitialParticle_UniformRandom(vector<UniformRandomBunch>& uniform_bunchs, ParticleMover* pm) {


  //double start,stop;
  int group;
  list<Particle>* particles;



  for(size_t k=0 ; k<uniform_bunchs.size() ; k++) {

    if(!uniform_bunchs[k].use) continue;

    group = uniform_bunchs[k].group;
    particles = pm->getParticles(group);
    if(particles == 0) {
      list<Particle> p;
      pm->getParticlesGroup().insert(pair< int, list<Particle> >(group, p));
      particles = pm->getParticles(group);
    }

    double range_x_min = uniform_bunchs[k].domain.xLower;
    double range_x_max = uniform_bunchs[k].domain.xUpper;
    double range_y_min = uniform_bunchs[k].domain.yLower;
    double range_y_max = uniform_bunchs[k].domain.yUpper;
    double range_z_min = uniform_bunchs[k].domain.zLower;
    double range_z_max = uniform_bunchs[k].domain.zUpper;

    double u,v,w,charge,mass;


    Parser parser;

    parser.add_variable("x");
    parser.add_variable("y");
    parser.add_variable("z");
    
    try {

      u = parser.parse(uniform_bunchs[k].physical_data.velocity_u.c_str());
      v = parser.parse(uniform_bunchs[k].physical_data.velocity_v.c_str());
      w = parser.parse(uniform_bunchs[k].physical_data.velocity_w.c_str());
      charge = parser.parse(uniform_bunchs[k].physical_data.charge.c_str());
      mass = parser.parse(uniform_bunchs[k].physical_data.mass.c_str());

      mass *= calculateLorentzFactor(u,v,w);

    
      double x,y,z;
      bool diagnose;
      for(size_t i=0 ; i<uniform_bunchs[k].num_particles ; ++i) {
        // uniform distribution, (0..1]
        x = (std::rand() + 1.0) / (RAND_MAX + 1.0);
        y = (std::rand() + 1.0) / (RAND_MAX + 1.0);
        z = (std::rand() + 1.0) / (RAND_MAX + 1.0);

        x = (range_x_max - range_x_min)*x + range_x_min;
        y = (range_y_max - range_y_min)*y + range_y_min;
        z = (range_z_max - range_z_min)*z + range_z_min;


        parser.set_variable(0, x);
        parser.set_variable(1, y);
        parser.set_variable(2, z);

        diagnose = isIncludedUniformBunch(parser, x,y,z, uniform_bunchs[k]);

        if(diagnose) {
          Particle p(x, y, z, u, v, w, charge, mass);
          particles->push_back(p);
        }

      }


    } catch(Error err) {
      EM_ERROR(err.get_msg());
      EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
      //printf("\tError: Unknown error occured in parser\n");

      return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
    }


  } // for k



  return EM_SUCCESS;

}



int generateInitialParticle_UserSpecific(UserSpecificSetting& is, ParticleMover* pm) {

  if(!is.use) return EM_SUCCESS;

  int group = is.group;
  list<Particle>* particles = pm->getParticles(group);
  if(particles == 0) {
    list<Particle> p;
    pm->getParticlesGroup().insert(pair< int, list<Particle> >(group, p));
    particles = pm->getParticles(group);
  }

  UserSpecificInitialParticleSetting* particleSetting = 0;

  if(is.classname == "DefaultUserSpecificInitialParticleSetting") {
    particleSetting = new DefaultUserSpecificInitialParticleSetting;
  } else if(is.classname == "TestParticleSetting1") {
    particleSetting = new TestParticleSetting1;
  } else {
    // Error
    EM_ERROR(STR_ERR_SYS_OBJECT_CREATION_USERSPECPARTICLE);
    return EM_ERR_SYS_OBJECT_CREATION_USERSPECPARTICLE;
  }

  double range_x_min, range_y_min, range_z_min, range_x_max, range_y_max, range_z_max;
  size_t xgridsize, ygridsize, zgridsize;

  particleSetting->getDomain(range_x_min, range_y_min, range_z_min, range_x_max, range_y_max, range_z_max);
  particleSetting->getGridSize(xgridsize, ygridsize, zgridsize);

  std::vector<bool> lsf;
  particleSetting->fillGrid(lsf);

  double sdz = is.dx;

  double dx = (range_x_max-range_x_min)/xgridsize;
  double dy = (range_y_max-range_y_min)/ygridsize;
  double dz = (range_z_max-range_z_min)/zgridsize;


  vector<double> initPackingX;
  vector<double> initPackingY;
  vector<double> initPackingZ;
  size_t num_particle_packing = 0;


  if(is.mode.compare("rect") == 0) {
    num_particle_packing = getInitialRectangularPacking3D(range_x_min, range_y_min, range_z_min, range_x_max, range_y_max, range_z_max, sdz, initPackingX, initPackingY, initPackingZ);
  } else { // defalut : "hex"
    num_particle_packing = getInitialHexagonalPacking3D(range_x_min, range_y_min, range_z_min, range_x_max, range_y_max, range_z_max, sdz, initPackingX, initPackingY, initPackingZ);
  }


  Parser parser;
    
  parser.add_variable("x");
  parser.add_variable("y");
  parser.add_variable("z");


  for (size_t count = 0; count < num_particle_packing; ++count) {  
	  //Find the grid cell which the points x[count] y[count] z[count] belongs to
	  int i =  (int)floor( (initPackingX[count] - range_x_min)/dx ) ;
	  int j =  (int)floor( (initPackingY[count] - range_y_min)/dy ) ;
	  int k =  (int)floor( (initPackingZ[count] - range_z_min)/dz ) ;
	  double cellvol = dx * dy * dz;
	  double temp; // An accumulator variable for the bilinear interpollation.
    size_t index = i + j*xgridsize + k*xgridsize*ygridsize;
    temp = (lsf[index + 1 + xgridsize] ? ( fabs( initPackingX[count] - (range_x_min + i*dx) ) * fabs( initPackingY[count] - (range_y_min + j*dy) ) * fabs(initPackingZ[count] - (range_z_min + (k+1)*dz)) ) : 0.0) /cellvol;
    temp += (lsf[index] ? ( fabs( initPackingX[count] - (range_x_min + (i+1) * dx)) * fabs( initPackingY[count] -(range_y_min + (j+1) * dy ) ) * fabs(initPackingZ[count] - (range_z_min + (k+1)*dz)) ) : 0.0) / cellvol;
    temp += (lsf[index + xgridsize] ? ( fabs( initPackingX[count] - (range_x_min + (i+1) * dx))  * fabs( initPackingY[count] -(range_y_min + j*dy) )  * fabs(initPackingZ[count] - (range_z_min + (k+1)* dz)) ) : 0) / cellvol;
	  temp += (lsf[index + 1] ? ( fabs( initPackingX[count] - (range_x_min + i * dx)  ) * fabs( initPackingY[count] -(range_y_min + (j+1) * dy ) )  * fabs(initPackingZ[count] - (range_z_min + (k+1)* dz)) ) : 0) / cellvol;
	  temp += (lsf[index + 1 + xgridsize + xgridsize*ygridsize] ? ( fabs( initPackingX[count] - (range_x_min + i * dx) ) * fabs( initPackingY[count] -(range_y_min + j*dy) )  * fabs(initPackingZ[count] - (range_z_min + (k)* dz)) ) : 0) /cellvol;
	  temp += (lsf[index + xgridsize*ygridsize] ? ( fabs( initPackingX[count] - (range_x_min + (i+1) * dx))  * fabs( initPackingY[count] -(range_y_min + (j+1) * dy ) ) * fabs(initPackingZ[count] - (range_z_min + (k)* dz)) ) : 0) / cellvol;
	  temp += (lsf[index + xgridsize + xgridsize*ygridsize] ? ( fabs( initPackingX[count] - (range_x_min + (i+1) * dx))  * fabs( initPackingY[count] -(range_x_min + j*dy) ) * fabs(initPackingZ[count] - (range_z_min + (k)* dz)) ) : 0) / cellvol;
	  temp += (lsf[index + 1 + xgridsize*ygridsize] ? ( fabs( initPackingX[count] - (range_x_min + i * dx)  ) * fabs( initPackingY[count] -(range_y_min + (j+1) * dy ) )  * fabs(initPackingZ[count] - (range_z_min + (k)* dz)) ) : 0) / cellvol;
	  
    //Based on the value of temp select a particle
	  if ( temp > particleSetting->getTolerance() ) {
      double x = initPackingX[count];
      double y = initPackingY[count];
      double z = initPackingZ[count];
    
      parser.set_variable(0, x);
      parser.set_variable(1, y);
      parser.set_variable(2, z);
      double u = parser.parse(is.physical_data.velocity_u.c_str());
      double v = parser.parse(is.physical_data.velocity_v.c_str());
      double w = parser.parse(is.physical_data.velocity_w.c_str());
      double mass = parser.parse(is.physical_data.mass.c_str());
      double charge = parser.parse(is.physical_data.charge.c_str());

      mass *= calculateLorentzFactor(u,v,w);

      Particle p(x, y, z, u, v, w, charge, mass);
      particles->push_back(p);
	  }
	}// End the for-loop. 



  return EM_SUCCESS;
}



int generateInitialParticles(ParticleInitialDistribution& pd, ParticleMover* pm) {

  int error;

  error = generateInitialParticle_DirectSetting(pd.directSetting, pm, pd.directSetting_group);
  if(error) return error;

  error = generateInitialParticle_Analytic(pd.bunchs, pm);
  if(error) return error;

  error = generateInitialParticle_NormalRandom(pd.normal_bunchs, pm);
  if(error) return error;

  error = generateInitialParticle_UniformRandom(pd.uniform_bunchs, pm);
  if(error) return error;

  error = generateInitialParticle_UserSpecific(pd.userSetting, pm);
  if(error) return error;


  return EM_SUCCESS;
  
}



bool isIncludedBunch(Parser& parser, double x, double y, double z, const Bunch& bunch) {

  size_t functionsSize = bunch.functions.size();
  bool diagnose = false;
  double answer = 0.0;

  for(size_t m=0 ; m<functionsSize ; ++m) {
    diagnose = true;
    bool temp_diag = true;
    LevelFunction levelfunction = bunch.functions[m];
    /*
    cout << "\nlevel function = " << levelfunction.function << endl;
    cout << "x = " << x <<endl;
    cout << "y = " << y <<endl;
    cout << "z = " << z <<endl;
    //*/
    answer = parser.parse(levelfunction.function.c_str());
    //cout << "answer = " << answer <<endl;

    if(levelfunction.range == '<') {
      if(levelfunction.equality) {
        temp_diag = (answer <= 0) ? true : false;
      } else {
        temp_diag = (answer < 0) ? true : false;
      }
    } else {
      if(levelfunction.equality) {
        temp_diag = (answer >= 0) ? true : false;
      } else {
        temp_diag = (answer > 0) ? true : false;
      }
    }
    //cout << "diagnose = " << diagnose <<endl;
    diagnose = diagnose && temp_diag;
    if(!diagnose) break;
    
  } // for m

  return diagnose;

}



bool isIncludedUniformBunch(Parser& parser, double x, double y, double z, const UniformRandomBunch& bunch) {

  size_t functionsSize = bunch.functions.size();
  bool diagnose = false;
  double answer = 0.0;

  for(size_t m=0 ; m<functionsSize ; ++m) {
    diagnose = true;
    bool temp_diag = true;
    LevelFunction levelfunction = bunch.functions[m];
    /*
    cout << "\nlevel function = " << levelfunction.function << endl;
    cout << "x = " << x <<endl;
    cout << "y = " << y <<endl;
    cout << "z = " << z <<endl;
    //*/
    answer = parser.parse(levelfunction.function.c_str());
    //cout << "answer = " << answer <<endl;

    if(levelfunction.range == '<') {
      if(levelfunction.equality) {
        temp_diag = (answer <= 0) ? true : false;
      } else {
        temp_diag = (answer < 0) ? true : false;
      }
    } else {
      if(levelfunction.equality) {
        temp_diag = (answer >= 0) ? true : false;
      } else {
        temp_diag = (answer > 0) ? true : false;
      }
    }
    //cout << "diagnose = " << diagnose <<endl;
    diagnose = diagnose && temp_diag;
    if(!diagnose) break;
    
  } // for m

  return diagnose;

}



double calculateLorentzFactor(double u, double v, double w) {
  double beta_2 = (u*u + v*v + w*w) / C0 / C0;
  return 1.0 / sqrt(1.0 - beta_2);
}




#ifndef _MSC_VER
void parseParameters(int argc, char *argv[], double& startTime, double& endTime, double& dt, double& leftX, double& rightX, double& dx, double& leftY, double& rightY, double& dy, double& leftZ, double& rightZ, double& dz, int& vis_interval, double& vis_leftX, double& vis_rightX, double& vis_leftY, double& vis_rightY, double& vis_leftZ, double& vis_rightZ, int& errormode) {

  int param_opt;
  struct option  longopts[23];

  longopts[0].name = "start_time";
  longopts[0].has_arg = 2;
  longopts[0].flag = NULL;
  longopts[0].val = -2;
  
  longopts[1].name = "end_time";
  longopts[1].has_arg = 2;
  longopts[1].flag = NULL;
  longopts[1].val = -3;

  longopts[2].name = "dt";
  longopts[2].has_arg = 2;
  longopts[2].flag = NULL;
  longopts[2].val = -4;
  
  longopts[3].name = "left_x";
  longopts[3].has_arg = 2;
  longopts[3].flag = NULL;
  longopts[3].val = -5;
  
  longopts[4].name = "right_x";
  longopts[4].has_arg = 2;
  longopts[4].flag = NULL;
  longopts[4].val = -6;

  longopts[5].name = "dx";
  longopts[5].has_arg = 2;
  longopts[5].flag = NULL;
  longopts[5].val = -7;

  longopts[6].name = "left_y";
  longopts[6].has_arg = 2;
  longopts[6].flag = NULL;
  longopts[6].val = -8;

  longopts[7].name = "right_y";
  longopts[7].has_arg = 2;
  longopts[7].flag = NULL;
  longopts[7].val = -9;

  longopts[8].name = "dy";
  longopts[8].has_arg = 2;
  longopts[8].flag = NULL;
  longopts[8].val = -10;

  longopts[9].name = "left_z";
  longopts[9].has_arg = 2;
  longopts[9].flag = NULL;
  longopts[9].val = -11;

  longopts[10].name = "right_z";
  longopts[10].has_arg = 2;
  longopts[10].flag = NULL;
  longopts[10].val = -12;

  longopts[11].name = "dz";
  longopts[11].has_arg = 2;
  longopts[11].flag = NULL;
  longopts[11].val = -13;

  longopts[12].name = "vis_interval";
  longopts[12].has_arg = 2;
  longopts[12].flag = NULL;
  longopts[12].val = -17;

  longopts[13].name = "vis_left_x";
  longopts[13].has_arg = 2;
  longopts[13].flag = NULL;
  longopts[13].val = -18;
  
  longopts[14].name = "vis_right_x";
  longopts[14].has_arg = 2;
  longopts[14].flag = NULL;
  longopts[14].val = -19;

  longopts[15].name = "vis_left_y";
  longopts[15].has_arg = 2;
  longopts[15].flag = NULL;
  longopts[15].val = -20;

  longopts[16].name = "vis_right_y";
  longopts[16].has_arg = 2;
  longopts[16].flag = NULL;
  longopts[16].val = -21;

  longopts[17].name = "vis_left_z";
  longopts[17].has_arg = 2;
  longopts[17].flag = NULL;
  longopts[17].val = -22;

  longopts[18].name = "vis_right_z";
  longopts[18].has_arg = 2;
  longopts[18].flag = NULL;
  longopts[18].val = -23;

  longopts[19].name = "errormode";
  longopts[19].has_arg = 2;
  longopts[19].flag = NULL;
  longopts[19].val = -24;



  while(-1 != (param_opt = getopt_long(argc, argv, ":",longopts, NULL)))
  {
    //printf("optarg = \'%s\'\n", optarg);
    //printf("optind = \'%d\'\n", optind);
    //printf("param_opt= \'%d\'\n", param_opt);
    switch(param_opt)
    {
      case -2 : 
        startTime = atof(optarg);
        break;
         
      case -3 :
        endTime = atof(optarg);
        break;

      case -4 :
        dt = atof(optarg);
        break;

      case -5 :
        leftX = atof(optarg);
        break;

      case -6 :
        rightX = atof(optarg);
        break;

      case -7 :
        dx = atof(optarg);
        break;

      case -8 :
        leftY = atof(optarg);
        break;

      case -9 :
        rightY = atof(optarg);
        break;

      case -10 :
        dy = atof(optarg);
        break;

      case -11 :
        leftZ = atof(optarg);
        break;

      case -12 :
        rightZ = atof(optarg);
        break;

      case -13 :
        dz = atof(optarg);
        break;

      /*
      case -14 :
        procNumX = atoi(optarg);
        break;

      case -15 :
        procNumY = atoi(optarg);
        break;

      case -16 :
        procNumZ = atoi(optarg);
        break;
      */

      case -17 :
        vis_interval = atoi(optarg);
        break;

      case -18 :
        vis_leftX = atof(optarg);
        break;

      case -19 :
        vis_rightX = atof(optarg);
        break;

      case -20 :
        vis_leftY = atof(optarg);
        break;

      case -21 :
        vis_rightY = atof(optarg);
        break;

      case -22 :
        vis_leftZ = atof(optarg);
        break;

      case -23 :
        vis_rightZ = atof(optarg);
        break;

      case -24 :
        errormode = atoi(optarg);
        break;


    }
   
  }

  
}


#endif 

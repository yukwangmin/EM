/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jun. 05 2012
 * 
 *      
 */



#include <cstdlib>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include <omp.h>
#include <deque>


#include "time_controller.h"
#include "parallel_solver.h"
#include "serial_solver.h"
#include "resultview.h"
#include "em_error.h"
#include "initial_setting.h"
#include "parser.h"
#include "constants.h"
#include "parallel_mover.h"
#include "user_specific_initial_particle_setting.h"


#ifndef _MSC_VER  
#include <getopt.h>
#include <unistd.h>
#endif 


using namespace mparser;



typedef struct struct_MPIParticle {
  double x;
  double y;
  double z;
  double u;
  double v;
  double w;
  double charge;
  double mass;
} MPIParticle;


/*
#ifndef _MSC_VER  
void parseParameters(int argc, char *argv[], double& startTime, double& endTime, double& dt, double& leftX, double& rightX, double& dx, double& leftY, double& rightY, double& dy, double& leftZ, double& rightZ, double& dz, int& procNumX, int& procNumY, int& procNumZ, int& vis_interval, double& vis_leftX, double& vis_rightX, double& vis_leftY, double& vis_rightY, double& vis_leftZ, double& vis_rightZ, int& errormode);
#endif 
*/

#ifndef _MSC_VER  
void parseParameters(int argc, char *argv[], int& procNumX, int& procNumY, int& procNumZ, char** filename);
#endif 


int validateInitialSetting(InitialSetting* is);
//int generateInitialParticles(ParticleInitialDistribution& pd, deque<MPIParticle>& particles);
bool isIncludedBunch(Parser& parser, double x, double y, double z, const Bunch& bunch);
double calculateLorentzFactor(double u, double v, double w);
void buildMPIDerivedType_MPIParticle(MPIParticle* ctype, MPI_Datatype* mpitype);
int generateInitialParticles(ParticleMover* mover, int totalProcNum, int rank, InitialSetting& is);
void distributeParticles(int rank, int totalProcNum, deque<MPIParticle>& particleQueue, int group, ParticleMover* pm);


int main(int argc, char *argv[]) {
/*
  double startTime = 0;
  double endTime = 0.000000001;
  //double dt =      0.000000000001;
  double dt = 75; // 75%
  double leftX = 0;
  double rightX = 0.256;
  double leftY = 0;
  double rightY = 0.256;
  double leftZ = 0;
  double rightZ = 0.256;
  double dx = 0.008;
  double dy = 0.008;
  double dz = 0.008;
  int vis_interval = -1;
  double vis_leftX = leftX;
  double vis_rightX = rightX;
  double vis_leftY = leftY;
  double vis_rightY = rightY;
  double vis_leftZ = leftZ;
  double vis_rightZ = rightZ;
  int errormode = 0;



  int procNumX = 2;
  int procNumY = 2;
  int procNumZ = 1;

  int totalProcNum;
  int rank;


#ifndef _MSC_VER
  parseParameters(argc, argv, startTime, endTime, dt, leftX, rightX, dx, leftY, rightY, dy, leftZ, rightZ, dz, procNumX, procNumY, procNumZ, vis_interval, vis_leftX, vis_rightX, vis_leftY, vis_rightY, vis_leftZ, vis_rightZ, errormode);
#endif 
*/


  
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

  int procNumX = 1;
  int procNumY = 1;
  int procNumZ = 1;

  int totalProcNum;
  int rank;


  
  int error;




  char* filename = const_cast<char*>("setting_MPI.xml");
  
#ifndef _MSC_VER
  parseParameters(argc, argv, procNumX, procNumY, procNumZ, &filename);
#endif 

  //cout << "filename : " << filename << endl;

  InitialSetting is;
  
  error = loadInitialSetting(filename, &is);
  if(error) return error;

  error = validateInitialSetting(&is);
  if(error) return error;



#ifdef _OPENMP
  omp_set_num_threads(static_cast<int>(is.omp_num_threads));
  //cout << "OMP_NUM_THREAD : " << omp_get_max_threads() << endl;
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
  //cout << "field_vis_interval : " << field_vis_interval << endl;

  particle_vis_interval = static_cast<int>((is.particle.visualization.timeinterval)/ndt);
  //cout << "particle_vis_interval : " << particle_vis_interval << endl;
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : Converting data from InitialSetting.
  ////////////////////////////////////////////////////////////////////////////////////////





  int gridSizeX = is.domain.xGridNumber;
  int gridSizeY = is.domain.yGridNumber;
  int gridSizeZ = is.domain.zGridNumber;
  if( (gridSizeX % procNumX != 0) || (gridSizeY % procNumY != 0) || (gridSizeZ % procNumZ != 0) ) {
    cout << "Invalid grid number and processor number : " << endl;
    cout << "Grid size along X : " << gridSizeX << ", But the number of processors along X : " << procNumX << endl;
    cout << "Grid size along Y : " << gridSizeY << ", But the number of processors along Y : " << procNumY << endl;
    cout << "Grid size along Z : " << gridSizeZ << ", But the number of processors along Z : " << procNumZ << endl;
	  exit(0);
  }
  


	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&totalProcNum);
  if(totalProcNum != procNumX*procNumY*procNumZ) {
    cout << "Invalid Processor number : " << totalProcNum << ". But # of X : " << procNumX << ", # of Y : " << procNumY << ", # of Z : " << procNumX << "." << endl;
    MPI_Finalize();
	  exit(0);
  }
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#ifndef _MSC_VER
  cout << "[PID:" << getpid() << "] = [RANK:" << rank << "]" << endl;
#else
  cout << "[PID:" << _getpid() << "] = [RANK:" << rank << "]" << endl;
#endif

  //cout << "Rank is " << rank << " out of " << totalProcNum << " processors." << endl;




  

  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : Initial Particle Setting.
  ////////////////////////////////////////////////////////////////////////////////////////
  
  ParticleMover* mover = new MPIParticleMover;
  error = generateInitialParticles(mover, totalProcNum, rank, is);
  if(error) return error;

  cout << "[Rank:" << rank << "]Total number of particles : " << mover->getTotalNumberOfParticles() << endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : Initial Particle Setting.
  ////////////////////////////////////////////////////////////////////////////////////////




  ////////////////////////////////////////////////////////////////////////////////////////
  // Start : FieldSolver Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////
  FieldSolver* solver = 0;
  if(is.field.solver == "MPIFieldSolver") {
    solver = new MPIFieldSolver(rank, procNumX, procNumY, procNumZ, ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "MPITestSolver") {
    solver = new MPITestSolver(rank, procNumX, procNumY, procNumZ, ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "MPINullFieldSolver") {
    solver = new MPINullFieldSolver(rank, procNumX, procNumY, procNumZ, ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else if(is.field.solver == "MPITestSolverWithParticle1") {
    solver = new MPITestSolverWithParticle1(rank, procNumX, procNumY, procNumZ, ndt, leftX, rightX, is.domain.xGridNumber, leftY, rightY, is.domain.yGridNumber, leftZ, rightZ, is.domain.zGridNumber);
  } else {
    // Error
    EM_ERROR(STR_ERR_SYS_OBJECT_CREATION_FIELDSOLVER);
    return EM_ERR_SYS_OBJECT_CREATION_FIELDSOLVER;
  }

  solver->initializeSolver();
  solver->setErrorMode(errormode);


  //MPI_Finalize();
  //exit(0);
  ////////////////////////////////////////////////////////////////////////////////////////
  // End : FieldSolver Object Creation.
  ////////////////////////////////////////////////////////////////////////////////////////





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
  // Start : ParticleMover Object Setting.
  ////////////////////////////////////////////////////////////////////////////////////////

  mover->setFieldSolver(solver);
  mover->setdt(ndt);
  mover->calculateInitialVelocity();

  ////////////////////////////////////////////////////////////////////////////////////////
  // End : ParticleMover Object Setting.
  ////////////////////////////////////////////////////////////////////////////////////////

  







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




  MPITimeController controller(rank, startTime, endTime, ndt, solver, mover, fieldViewer, particleViewer, field_vis_interval, particle_vis_interval);
  //TimeController controller(startTime, endTime, ndt, solver, mover, 0, particleViewer, vis_interval);
  //TimeController controller(startTime, endTime, ndt, solver, fieldViewer, vis_interval);
  //TimeController controller(startTime, endTime, dt, solver);










  if(rank == 0) {
    cout << "===============================================================================" << endl;
    cout << "leftX : " << leftX << ", rightX : " << rightX << ", dx : " << dx << ", Grid Size X : " << gridSizeX << endl;
    cout << "leftY : " << leftY << ", rightY : " << rightY << ", dy : " << dy << ", Grid Size Y : " << gridSizeY << endl;
    cout << "leftZ : " << leftZ << ", rightZ : " << rightZ << ", dz : " << dz << ", Grid Size Z : " << gridSizeZ << endl;
    cout << "dt : " << ndt << ", last time index : " << static_cast<int>(((endTime - startTime)/ndt) + ndt) << endl;
    cout << "particle_vis_interval : " << particle_vis_interval << endl;
    cout << "===============================================================================" << endl;
  }


  //MPI_Finalize();
  //return 0;


  if(controller.solve()) {
    EM_ERROR("Controller solve() Error");

    if(fieldViewer != 0) delete fieldViewer;
    if(particleViewer != 0) delete particleViewer;

    delete solver;
    delete mover;

    return 0;
  }



  if(fieldViewer != 0) delete fieldViewer;
  if(particleViewer != 0) delete particleViewer;

  delete solver;
  delete mover;


  MPI_Finalize();
  return 0;


  
}




int validateInitialSetting(InitialSetting* is) {

  // Domain Check.
  if(is->domain.xLower >= is->domain.xUpper) {
    cout << STR_ERR_COM_DOM_SCOPE_X << endl;
    cout << "<x lower=\"" << is->domain.xLower << "\" upper=\"" << is->domain.xUpper << "\" gridnumber=\"" << is->domain.xGridNumber << "\" />" << endl;
    return EM_ERR_COM_DOM_SCOPE_X;
  }
  if(is->domain.yLower >= is->domain.yUpper) {
    cout << STR_ERR_COM_DOM_SCOPE_Y << endl;
    cout << "<y lower=\"" << is->domain.yLower << "\" upper=\"" << is->domain.yUpper << "\" gridnumber=\"" << is->domain.yGridNumber << "\" />" << endl;
    return EM_ERR_COM_DOM_SCOPE_Y;
  }
  if(is->domain.zLower >= is->domain.zUpper) {
    cout << STR_ERR_COM_DOM_SCOPE_Z << endl;
    cout << "<z lower=\"" << is->domain.zLower << "\" upper=\"" << is->domain.zUpper << "\" gridnumber=\"" << is->domain.zGridNumber << "\" />" << endl;
    return EM_ERR_COM_DOM_SCOPE_Z;
  }


  return EM_SUCCESS;
}







int generateInitialParticle_DirectSetting(int rank, list<InitialParticle>& directSetting, deque<MPIParticle>& particleQueue, int group) {
  if(rank == 0) {

    try {
      Parser parser;

      parser.add_variable("x");
      parser.add_variable("y");
      parser.add_variable("z");

      list<InitialParticle> particles = directSetting;
      for(list<InitialParticle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
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

        MPIParticle p;
        p.x = x;
        p.y = y;
        p.z = z;
        p.u = u;
        p.v = v;
        p.w = w;
        p.charge = charge;
        p.mass = mass;
        particleQueue.push_back(p);
      }
    } catch(...) {
      EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
      return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
    }

  }

  return EM_SUCCESS;



}

int generateInitialParticle_Analytic(int rank, int totalProcNum, vector<Bunch>& bunchs, deque<MPIParticle>& particleQueue, ParticleMover* pm, size_t& totalNumParticle) {

  //double start,stop;
  int group;

  for(size_t k=0 ; k<bunchs.size() ; k++) {

    if(!bunchs[k].use) continue;

    group = bunchs[k].group;


    if(rank == 0) {
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

        bool diagnose;
        
        size_t total_num_particles = 0;

        for(size_t i=0 ; i<num_particle_packing ; ++i) {

          x = initPackingX[i];
          y = initPackingY[i];
          z = initPackingZ[i];
        
          parser.set_variable(0, x);
          parser.set_variable(1, y);
          parser.set_variable(2, z);

          diagnose = isIncludedBunch(parser, x,y,z, bunchs[k]);

          if(diagnose) {
            total_num_particles++;
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

            MPIParticle p;
            p.x = x;
            p.y = y;
            p.z = z;
            p.u = u;
            p.v = v;
            p.w = w;
            p.charge = charge;
            p.mass = mass;
            particleQueue.push_back(p);
          

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

    } // if(rank == 0)

    totalNumParticle += particleQueue.size();
    distributeParticles(rank, totalProcNum, particleQueue, group, pm);


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

  } // for k

  return EM_SUCCESS;

}

int generateInitialParticle_NormalRandom(int rank, int totalProcNum, vector<NormalRandomBunch>& normal_bunchs, deque<MPIParticle>& particleQueue, ParticleMover* pm, size_t& totalNumParticle) {

  //double start,stop;
  int group;
  //cout << "normal_bunchs.size() : " << normal_bunchs.size() << endl;

  for(size_t k=0 ; k<normal_bunchs.size() ; k++) {

    if(!normal_bunchs[k].use) continue;

    group = normal_bunchs[k].group;
    //cout << "rank=" << rank << ", k=" << k << ", group number : " << group << endl;

    if(rank == 0) {

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

        MPIParticle p;
        p.x = x;
        p.y = y;
        p.z = z;
        p.u = u;
        p.v = v;
        p.w = w;
        p.charge = charge;
        p.mass = mass;
        particleQueue.push_back(p);
      }
    } //if(rank == 0)


    totalNumParticle += particleQueue.size();
    //cout << "rank=" << rank << ", k=" << k << ", Just before calling distributeParticles() : group : " << group << endl;

    //MPI_Barrier(MPI_COMM_WORLD);
    distributeParticles(rank, totalProcNum, particleQueue, group, pm);


  } // for k

  return EM_SUCCESS;

}


int generateInitialParticle_UserSpecific(int rank, int totalProcNum, UserSpecificSetting& is, deque<MPIParticle>& particleQueue) {

  if(!is.use) return EM_SUCCESS;

  if(rank == 0) {

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
        
        MPIParticle p;
        p.x = x;
        p.y = y;
        p.z = z;
        p.u = u;
        p.v = v;
        p.w = w;
        p.charge = charge;
        p.mass = mass;
        particleQueue.push_back(p);
	    }
	  }// End the for-loop. 

  } //if(rank == 0)

  return EM_SUCCESS;
}


/*
int generateInitialParticles(ParticleInitialDistribution& pd, deque<MPIParticle>& particles) {

  int error;


  error = generateInitialParticle_DirectSetting(rank, pd.directSetting, particles, pd.directSetting_group);
  if(error) return error;

  error = generateInitialParticle_Analytic(rank, pd.bunchs, particles);
  if(error) return error;
  
  error = generateInitialParticle_NormalRandom(rank, pd.normal_bunchs, particles);
  if(error) return error;

  //error = generateInitialParticle_UniformRandom(rank, pd.uniform_bunchs, particles);
  //if(error) return error;

  error = generateInitialParticle_UserSpecific(rank, pd.userSetting, particles);
  if(error) return error;



  return EM_SUCCESS;
  
}
*/





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






double calculateLorentzFactor(double u, double v, double w) {
  double beta_2 = (u*u + v*v + w*w) / C0 / C0;
  return 1.0 / sqrt(1.0 - beta_2);
}







void buildMPIDerivedType_MPIParticle(MPIParticle* ctype, MPI_Datatype* mpitype) {
  const int num_var = 8;

  int block_length[num_var];
  MPI_Aint mem_location[num_var];
  MPI_Datatype data_type[num_var];

  MPI_Aint baseAddr;

  MPI_Address(ctype, &baseAddr);


  MPI_Address(&(ctype->x), &mem_location[0]);
  MPI_Address(&(ctype->y), &mem_location[1]);
  MPI_Address(&(ctype->z), &mem_location[2]);
  MPI_Address(&(ctype->u), &mem_location[3]);
  MPI_Address(&(ctype->v), &mem_location[4]);
  MPI_Address(&(ctype->w), &mem_location[5]);
  MPI_Address(&(ctype->charge), &mem_location[6]);
  MPI_Address(&(ctype->mass), &mem_location[7]);


  for(int i=0 ; i<num_var ; ++i) {
    block_length[i] = 1;
    mem_location[i] -= baseAddr;
    data_type[i] = MPI_DOUBLE;
  }

  MPI_Type_struct(num_var, block_length, mem_location, data_type, mpitype);
  MPI_Type_commit(mpitype);

}


int generateInitialParticles(ParticleMover* mover, int totalProcNum, int rank, InitialSetting& is) {

  int error;
  size_t totalNumParticle = 0; // Only valid for rank 0.

  deque<MPIParticle> particleQueue;
  

  error = generateInitialParticle_DirectSetting(rank, is.particle.distribution.directSetting, particleQueue, is.particle.distribution.directSetting_group);
  if(error) {
    EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
    return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  }
  
  totalNumParticle += particleQueue.size();
  distributeParticles(rank, totalProcNum, particleQueue, is.particle.distribution.directSetting_group, mover);

  if(rank == 0) {
    if(particleQueue.size()) cout << "Error : particleQueue is not empty after particle distributing of direct setting particles: particleQueue.size() : " << particleQueue.size() << endl;
    
  }


  error = generateInitialParticle_Analytic(rank, totalProcNum, is.particle.distribution.bunchs, particleQueue, mover, totalNumParticle);
  if(error) {
    EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
    return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  }
  if(rank == 0) {
    if(particleQueue.size()) cout << "Error : particleQueue is not empty after particle distributing of analytic particles: particleQueue.size() : " << particleQueue.size() << endl;
  }
  
  error = generateInitialParticle_NormalRandom(rank, totalProcNum, is.particle.distribution.normal_bunchs, particleQueue, mover, totalNumParticle);
  if(error) {
    EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
    return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  }
  if(rank == 0) {
    if(particleQueue.size()) cout << "Error : particleQueue is not empty after particle distributing of normal random particles: particleQueue.size() : " << particleQueue.size() << endl;
  }

  //error = generateInitialParticle_UniformRandom(rank, totalProcNum, is.particle.distribution.uniform_bunchs, particleQueue, mover, totalNumParticle);
  //if(error) {
  //  EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
  //  return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  //}

  error = generateInitialParticle_UserSpecific(rank, totalProcNum, is.particle.distribution.userSetting, particleQueue);
  if(error) {
    EM_ERROR(STR_ERR_SYS_INITIAL_SETTING_PARTICLE);
    return EM_ERR_SYS_INITIAL_SETTING_PARTICLE;
  }

  totalNumParticle += particleQueue.size();
  distributeParticles(rank, totalProcNum, particleQueue, is.particle.distribution.userSetting.group, mover);

  if(rank == 0) {
    if(particleQueue.size()) cout << "Error : particleQueue is not empty after particle distributing of user specific particles: particleQueue.size() : " << particleQueue.size() << endl;
  }

  if(rank == 0) cout << "The total number of generated particles is " << totalNumParticle << endl;
  
  return EM_SUCCESS;

}


void distributeParticles(int rank, int totalProcNum, deque<MPIParticle>& particleQueue, int group, ParticleMover* pm) {

  //cout << "rank=" << rank << ", distributeParticles() : param group : " << group << endl;

  size_t numParticlesPerProc;
  size_t residue;

  MPIParticle* sendBuffer;
  MPIParticle* recvBuffer;

  MPIParticle cMPI_particle;
  MPI_Datatype MPI_Particle;
  buildMPIDerivedType_MPIParticle(&cMPI_particle, &MPI_Particle);
  size_t numScatter; // valid for only rank 0.

  //start = omp_get_wtime();
  if(rank == 0)  {
    size_t numParticles = particleQueue.size();
    numParticlesPerProc = numParticles/totalProcNum;
    residue = numParticles%totalProcNum;
    numScatter = numParticlesPerProc*totalProcNum;
  }

  MPI_Bcast(&numParticlesPerProc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&residue, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  if(numParticlesPerProc == 0 && residue == 0) return;


  if(rank == 0) {
    sendBuffer = new MPIParticle[numScatter]; //(MPIParticle*)(malloc(numScatter*sizeof(MPIParticle)));

    for(size_t i=0 ; i<numScatter ; ++i) {
      MPIParticle& p = particleQueue.front();
      sendBuffer[i].x = p.x;
      sendBuffer[i].y = p.y;
      sendBuffer[i].z = p.z;
      sendBuffer[i].u = p.u;
      sendBuffer[i].v = p.v;
      sendBuffer[i].w = p.w;
      sendBuffer[i].charge = p.charge;
      sendBuffer[i].mass = p.mass;

      particleQueue.pop_front();
    }
  }


  
  

  recvBuffer = new MPIParticle[numParticlesPerProc]; //(MPIParticle*)(malloc(numParticlesPerProc*sizeof(MPIParticle)));


  MPI_Scatter(sendBuffer, numParticlesPerProc, MPI_Particle, recvBuffer, numParticlesPerProc, MPI_Particle, 0, MPI_COMM_WORLD);

  list<Particle>* particles = pm->getParticles(group);
  if(particles == 0) {
    list<Particle> p;
    pm->getParticlesGroup().insert(pair< int, list<Particle> >(group, p));
    particles = pm->getParticles(group);
  }
  
  for(size_t i=0 ; i<numParticlesPerProc ; ++i) {
    MPIParticle& mpi = recvBuffer[i];
    Particle p(mpi.x, mpi.y, mpi.z, mpi.u, mpi.v, mpi.w, mpi.charge, mpi.mass);
    particles->push_back(p);
  }

  if(residue > 0) {
    if(rank == 0) {
      MPIParticle& mpi = particleQueue.front();
      Particle p(mpi.x, mpi.y, mpi.z, mpi.u, mpi.v, mpi.w, mpi.charge, mpi.mass);
      particles->push_back(p);
      particleQueue.pop_front();

      for(size_t i=1 ; i<residue ; ++i) {
        MPIParticle& p = particleQueue.front();
        MPI_Send(&p, 1, MPI_Particle, i,  0, MPI_COMM_WORLD);
        particleQueue.pop_front();
      }
      //cout << "Initial Particle Distributing Time is " << end - start << " seconds." << endl;
    } else {
      if(rank < static_cast<int>(residue)) {
        MPIParticle mpi;
        MPI_Recv(&mpi, 1, MPI_Particle, 0,  MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        Particle p(mpi.x, mpi.y, mpi.z, mpi.u, mpi.v, mpi.w, mpi.charge, mpi.mass);
        particles->push_back(p);
      }
    }
  }


  delete[] recvBuffer;//free(recvBuffer);
  if(rank == 0) delete[] sendBuffer;//free(sendBuffer);



}



#ifndef _MSC_VER
void parseParameters(int argc, char *argv[], int& procNumX, int& procNumY, int& procNumZ, char** filename) {

  int param_opt;
  struct option  longopts[4];

  //cout <<

  longopts[0].name = "procNumX";
  longopts[0].has_arg = 2;
  longopts[0].flag = NULL;
  longopts[0].val = -2;

  longopts[1].name = "procNumY";
  longopts[1].has_arg = 2;
  longopts[1].flag = NULL;
  longopts[1].val = -3;

  longopts[2].name = "procNumZ";
  longopts[2].has_arg = 2;
  longopts[2].flag = NULL;
  longopts[2].val = -4;

  longopts[3].name = "setting";
  longopts[3].has_arg = 2;
  longopts[3].flag = NULL;
  longopts[3].val = -5;
  

  while(-1 != (param_opt = getopt_long(argc, argv, ":",longopts, NULL)))
  {
    //printf("optarg = \'%s\'\n", optarg);
    //printf("optind = \'%d\'\n", optind);
    //printf("param_opt= \'%d\'\n", param_opt);
    switch(param_opt)
    {
      case -2 : 
        procNumX = atoi(optarg);
        break;
         
      case -3 :
        procNumY = atoi(optarg);
        break;

      case -4 :
        procNumZ = atoi(optarg);
        break;

      case -5 :
        *filename = optarg;
        break;
      
    }
   
  }

  
}


#endif 












/*
#ifndef _MSC_VER
void parseParameters(int argc, char *argv[], double& startTime, double& endTime, double& dt, double& leftX, double& rightX, double& dx, double& leftY, double& rightY, double& dy, double& leftZ, double& rightZ, double& dz, int& procNumX, int& procNumY, int& procNumZ, int& vis_interval, double& vis_leftX, double& vis_rightX, double& vis_leftY, double& vis_rightY, double& vis_leftZ, double& vis_rightZ, int& errormode) {

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

  longopts[12].name = "procNumX";
  longopts[12].has_arg = 2;
  longopts[12].flag = NULL;
  longopts[12].val = -14;

  longopts[13].name = "procNumY";
  longopts[13].has_arg = 2;
  longopts[13].flag = NULL;
  longopts[13].val = -15;

  longopts[14].name = "procNumZ";
  longopts[14].has_arg = 2;
  longopts[14].flag = NULL;
  longopts[14].val = -16;
  
  longopts[15].name = "vis_interval";
  longopts[15].has_arg = 2;
  longopts[15].flag = NULL;
  longopts[15].val = -17;
  
  longopts[16].name = "vis_left_x";
  longopts[16].has_arg = 2;
  longopts[16].flag = NULL;
  longopts[16].val = -18;
  
  longopts[17].name = "vis_right_x";
  longopts[17].has_arg = 2;
  longopts[17].flag = NULL;
  longopts[17].val = -19;

  longopts[18].name = "vis_left_y";
  longopts[18].has_arg = 2;
  longopts[18].flag = NULL;
  longopts[18].val = -20;

  longopts[19].name = "vis_right_y";
  longopts[19].has_arg = 2;
  longopts[19].flag = NULL;
  longopts[19].val = -21;

  longopts[20].name = "vis_left_z";
  longopts[20].has_arg = 2;
  longopts[20].flag = NULL;
  longopts[20].val = -22;

  longopts[21].name = "vis_right_z";
  longopts[21].has_arg = 2;
  longopts[21].flag = NULL;
  longopts[21].val = -23;
  
  longopts[22].name = "errormode";
  longopts[22].has_arg = 2;
  longopts[22].flag = NULL;
  longopts[22].val = -24;



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

      case -14 :
        procNumX = atoi(optarg);
        break;

      case -15 :
        procNumY = atoi(optarg);
        break;

      case -16 :
        procNumZ = atoi(optarg);
        break;
        
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
*/



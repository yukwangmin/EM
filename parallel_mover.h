/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Tue Aug. 07 2012
 * 
 *      
 */


#ifndef __PARALLEL_MOVER_H__
#define __PARALLEL_MOVER_H__


#include <vector>

#include "particle.h"
#include "parallel_solver.h"
#include "particle_mover.h"




class MPIParticleMover : public ParticleMover {
private:
  MPIFieldSolver* m_pFieldSolver;
  double m_dt;
  map< int, list<Particle> > m_mParticlesGroup;
  //list<Particle> m_lParticles;
  //size_t m_iTotalNumParticles;

  ////////////////////////////////////////////////////////////////////
  // MPIFieldSolver info.
  ////////////////////////////////////////////////////////////////////
  /*
  int m_iTotalNumProc;
  int m_iRank;
  int m_iProcNumX;
  int m_iProcNumY;
  int m_iProcNumZ;
  int m_iGridSizeX;
  int m_iGridSizeY;
  int m_iGridSizeZ;
  double m_dx;
  double m_dy;
  double m_dz;
  double m_fGlobalLowerX;
  double m_fGlobalUpperX;
  double m_fGlobalLowerY;
  double m_fGlobalUpperY;
  double m_fGlobalLowerZ;
  double m_fGlobalUpperZ;
  */
  ////////////////////////////////////////////////////////////////////
  // MPIFieldSolver info.
  ////////////////////////////////////////////////////////////////////


  


public:
  MPIParticleMover() {}
  MPIParticleMover(MPIFieldSolver* solver, double dt) : m_pFieldSolver(solver), m_dt(dt) {}
  virtual ~MPIParticleMover() {}

public:
  virtual void setFieldSolver(FieldSolver* solver) {m_pFieldSolver = dynamic_cast<MPIFieldSolver*>(solver);}
  virtual void setdt(double dt) {m_dt = dt;}
  virtual int addParticle(Particle& p, int group = 0);
  virtual size_t getNumberOfParticles(int group = 0);
  virtual size_t getTotalNumberOfParticles();
  virtual list<Particle>* getParticles(int group = 0);
  virtual map< int, list<Particle> >& getParticlesGroup() {return m_mParticlesGroup;}

  /**
   * Interpolate from EM grid in *solver to particle.
   * After that apply Lorentz force.
   * After that calculate current.
   * Then update current density to E in *solver
   */
  virtual int moveParticles(double time, size_t timestep);
  
  // Calculate initial half timestep velocity.
  virtual void calculateInitialVelocity();


private:
  

  /**
   * @return -1 if a particle passes more than 1 cell along a coordination direction. This means a violation of CFL.
   */
  int updateCurrentDensity(double time, double charge, double x, double y, double z, double new_x, double new_y, double new_z);

  void updateCurrentDensityInOneCell(double time, double charge, int cellX, int cellY, int cellZ, double x, double y, double z, double dx, double dy, double dz, double new_x, double new_y, double new_z, double ndx, double ndy, double ndz);

  
  /*
  int getElectricIntensityAt(double time, double x, double y, double z, double& ex, double& ey, double& ez);
  int getMagneticDensityAt(double time, double x, double y, double z, double& bx, double& by, double& bz);
  int findCellHavingPosition(double x, double y, double z, int& blockX, int& blockY, int& blockZ, int& cellX, int& cellY, int& cellZ, double& dx, double& dy, double& dz);
  */

  
};





#endif //__PARALLEL_MOVER_H__



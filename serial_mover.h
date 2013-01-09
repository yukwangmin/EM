/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Tue Aug. 07 2012
 * 
 *      
 */


#ifndef __SERIAL_MOVER_H__
#define __SERIAL_MOVER_H__


#include <vector>
#include <list>
#include <map>


#include "particle.h"
#include "serial_solver.h"
#include "particle_mover.h"



class SerialParticleMover : public ParticleMover {
private:
  SerialFieldSolver* m_pFieldSolver;
  double m_dt;
  map< int, list<Particle> > m_mParticlesGroup;
  //list<Particle> m_lParticles;
  //size_t m_iTotalNumParticles;


public:
  SerialParticleMover(SerialFieldSolver* solver, double dt) : m_pFieldSolver(solver), m_dt(dt) {}
  virtual ~SerialParticleMover() {}

public:
  virtual void setFieldSolver(FieldSolver* solver) {m_pFieldSolver = dynamic_cast<SerialFieldSolver*>(solver);}
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
  virtual void initElectricIntensityByParticles();


private:
  /**
   * @return -1 if a particle passes more than 1 cell along a coordination direction. This means a violation of CFL.
   */
  int updateCurrentDensity(double time, double charge, double x, double y, double z, double new_x, double new_y, double new_z);

  void updateCurrentDensityInOneCell(double time, double charge, int cellX, int cellY, int cellZ, double x, double y, double z, double dx, double dy, double dz, double new_x, double new_y, double new_z, double ndx, double ndy, double ndz);

};





#endif //__SERIAL_MOVER_H__



/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Tue Aug. 07 2012
 * 
 *      
 */


#ifndef __PARTICLE_MOVER_H__
#define __PARTICLE_MOVER_H__


#include <map>
#include <list>

#include "particle.h"
#include "solver.h"


using namespace std;



class ParticleMover {

public:
  virtual void setFieldSolver(FieldSolver* solver) = 0;
  virtual void setdt(double dt) = 0;
  virtual int addParticle(Particle& p, int group) = 0;
  virtual size_t getNumberOfParticles(int group) = 0;
  virtual size_t getTotalNumberOfParticles() = 0;
  virtual list<Particle>* getParticles(int group) = 0;
  virtual map< int, list<Particle> >& getParticlesGroup() = 0;
  virtual int moveParticles(double time, size_t timestep) = 0;

  // Calculate initial half timestep velocity.
  virtual void calculateInitialVelocity() = 0;
  virtual void initElectricIntensityByParticles() = 0;
  
};



#endif //__PARTICLE_MOVER_H__
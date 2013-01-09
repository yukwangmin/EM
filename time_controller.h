/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jun. 05 2012
 * 
 *      
 */



#ifndef __TIME_CONTROLLER_H__
#define __TIME_CONTROLLER_H__


#include <vector>


#include "solver.h"
#include "resultview.h"
#include "particle_mover.h"


using namespace std;


class TimeController {
protected:
  double m_fStartTime;
  double m_fEndTime;
  double m_dt;
  size_t m_iLastTimeIndex;
  int m_iWriteFieldInterval;
  int m_iWriteParticleInterval;
  //string m_sFilename;
  //string m_sExtension;

  FieldViewer* m_pFieldViewer;
  ParticleViewer* m_pParticleViewer;
  FieldSolver* m_pSolver;
  ParticleMover* m_pMover;



public:
  //*
  TimeController(double startTime, double endTime, double dt, FieldSolver* solver);
  TimeController(double startTime, double endTime, double dt, FieldSolver* solver, FieldViewer* viewer, int writeInterval = 1);
  TimeController(double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover);
  TimeController(double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover, FieldViewer* viewer, int writeInterval = 1);
  //*/
  TimeController(double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover, FieldViewer* fviewer, ParticleViewer* pviewer, int writeFieldInterval = 1, int writeParticleInterval = 1);
  


public:
  virtual int solve();




};



class MPITimeController : public TimeController {
private:
  int m_iRank;
  /*
  double m_fStartTime;
  double m_fEndTime;
  double m_dt;
  size_t m_iLastTimeIndex;
  
  Solver* m_pSolver;  
  */



public:
  MPITimeController(int rank, double startTime, double endTime, double dt, FieldSolver* solver);
  MPITimeController(int rank, double startTime, double endTime, double dt, FieldSolver* solver, FieldViewer* viewer, int saveInterval = 1);
  MPITimeController(int rank, double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover, FieldViewer* fviewer, ParticleViewer* pviewer, int writeFieldInterval = 1, int writeParticleInterval = 1);


public:
  virtual int solve();


};


#endif
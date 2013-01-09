/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jun. 05 2012
 * 
 *      
 */


#include <ctime>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <cmath>
#include <sstream>


#include "time_controller.h"
#include "em_error.h"


int getLastTimeIndex(double start, double end, double dt) {return static_cast<int>(((end - start)/dt) + dt);}


//*
TimeController::TimeController(double startTime, double endTime, double dt, FieldSolver* solver)
: m_fStartTime(startTime), m_fEndTime(endTime), m_dt(dt), m_pSolver(solver), m_pMover(0), m_pFieldViewer(0), m_iWriteFieldInterval(1)
{
  m_iLastTimeIndex = getLastTimeIndex(m_fStartTime, m_fEndTime, m_dt);
}


TimeController::TimeController(double startTime, double endTime, double dt, FieldSolver* solver, FieldViewer* viewer, int writeInterval)
: m_fStartTime(startTime), m_fEndTime(endTime), m_dt(dt), m_pSolver(solver), m_pMover(0), m_pFieldViewer(viewer), m_iWriteFieldInterval(writeInterval)
{
  m_iLastTimeIndex = static_cast<int>(((m_fEndTime - m_fStartTime)/dt) + m_dt);
}

TimeController::TimeController(double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover)
: m_fStartTime(startTime), m_fEndTime(endTime), m_dt(dt), m_pSolver(solver), m_pMover(mover), m_pFieldViewer(0), m_iWriteFieldInterval(1)
{
  m_iLastTimeIndex = static_cast<int>(((m_fEndTime - m_fStartTime)/dt) + m_dt);
}

TimeController::TimeController(double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover, FieldViewer* viewer, int writeInterval)
: m_fStartTime(startTime), m_fEndTime(endTime), m_dt(dt), m_pSolver(solver), m_pMover(mover), m_pFieldViewer(viewer), m_iWriteFieldInterval(writeInterval)
{
  m_iLastTimeIndex = static_cast<int>(((m_fEndTime - m_fStartTime)/dt) + m_dt);
}
//*/
TimeController::TimeController(double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover, FieldViewer* fviewer, ParticleViewer* pviewer, int writeFieldInterval, int writeParticleInterval)
: m_fStartTime(startTime), m_fEndTime(endTime), m_dt(dt), m_pSolver(solver), m_pMover(mover), m_pFieldViewer(fviewer), m_pParticleViewer(pviewer), m_iWriteFieldInterval(writeFieldInterval), m_iWriteParticleInterval(writeParticleInterval)
{
  m_iLastTimeIndex = getLastTimeIndex(m_fStartTime, m_fEndTime, m_dt);
}



int TimeController::solve() {
  double start, end;
  double total_start, total_end;

  total_start = omp_get_wtime();

  if(m_pFieldViewer != 0) {
    m_pFieldViewer->writeResult(0, 0);
  }

  if(m_pParticleViewer != 0) {
    m_pParticleViewer->writeResult(0, 0);
  }

  
  for(size_t i=1 ; i<=m_iLastTimeIndex ; i++) {
    //if(i == 4) {
    //  cout << endl;
    //}
    double time = m_dt*i;

    //clock_t start,stop;
    //start = clock();
    start = omp_get_wtime();
    m_pSolver->solve(time);
    //stop = clock();
    end = omp_get_wtime();
    //cout << "i : "<< i << " : Solving Time is " << (stop-start)/(double) CLOCKS_PER_SEC << " seconds." << endl;
    cout << "i : "<< i << " : Filed Solving Time is " << end-start << " seconds." << endl;

    start = omp_get_wtime();
    if(m_pMover != 0) {
      if(int err = m_pMover->moveParticles(time, i)) {
        EM_ERROR("Moving Particles Error");
        return err;
      }
    }
    end = omp_get_wtime();
    cout << "i : "<< i << " : Particle Moving Time is " << end-start << " seconds." << endl;



    
    if((m_pFieldViewer != 0) && (i % m_iWriteFieldInterval == 0)) {
      double wstart = omp_get_wtime();
      m_pFieldViewer->writeResult(time, static_cast<int>(i));
      double wend = omp_get_wtime();
      cout << "i : "<< i << " : FieldViewer Saving Time is " << wend-wstart << " seconds." << endl;
    }

    if((m_pParticleViewer != 0) && (i % m_iWriteParticleInterval == 0)) {
      double wstart = omp_get_wtime();
      m_pParticleViewer->writeResult(time, static_cast<int>(i));
      double wend = omp_get_wtime();
      cout << "i : "<< i << " : ParticleViewer Saving Time is " << wend-wstart << " seconds." << endl;
    }
    

    /*
    double Ex = 0.0;
    double Ey = 0.0;
    double Ez = 0.0;
    double Hx = 0.0;
    double Hy = 0.0;
    double Hz = 0.0;

    m_pSolver->evaluateError(i*m_dt, Ex, Ey, Ez, Hx, Hy, Hz, 1e-20);
    cout << "Timestep : " << i << " : Error : Ex : " << Ex << endl;
    cout << "Timestep : " << i << " : Error : Ey : " << Ey << endl;
    cout << "Timestep : " << i << " : Error : Ez : " << Ez << endl;
    cout << "Timestep : " << i << " : Error : Hx : " << Hx << endl;
    cout << "Timestep : " << i << " : Error : Hy : " << Hy << endl;
    cout << "Timestep : " << i << " : Error : Hz : " << Hz << endl;
    //*/

  }
  total_end = omp_get_wtime();
  cout << "Solving Time is " << total_end-total_start << " seconds." << endl;


  
  /*
  double Ex = 0.0;
  double Ey = 0.0;
  double Ez = 0.0;
  double Hx = 0.0;
  double Hy = 0.0;
  double Hz = 0.0;

  m_pSolver->evaluateError(m_iLastTimeIndex*m_dt, Ex, Ey, Ez, Hx, Hy, Hz, 1e-20);
  cout << "Error : Ex : " << Ex << endl;
  cout << "Error : Ey : " << Ey << endl;
  cout << "Error : Ez : " << Ez << endl;
  cout << "Error : Hx : " << Hx << endl;
  cout << "Error : Hy : " << Hy << endl;
  cout << "Error : Hz : " << Hz << endl;
  //*/

  return 0;

}





MPITimeController::MPITimeController(int rank, double startTime, double endTime, double dt, FieldSolver* solver)
: TimeController(startTime, endTime, dt, solver), m_iRank(rank)
{
  //m_iLastTimeIndex = static_cast<int>(((m_fEndTime - m_fStartTime)/dt) + m_dt);
}



MPITimeController::MPITimeController(int rank, double startTime, double endTime, double dt, FieldSolver* solver, FieldViewer* viewer, int writeInterval)
: TimeController(startTime, endTime, dt, solver, viewer, writeInterval), m_iRank(rank)
{
  //m_iLastTimeIndex = static_cast<int>(((m_fEndTime - m_fStartTime)/dt) + m_dt);
}

MPITimeController::MPITimeController(int rank, double startTime, double endTime, double dt, FieldSolver* solver, ParticleMover* mover, FieldViewer* fviewer, ParticleViewer* pviewer, int writeFieldInterval, int writeParticleInterval)
: TimeController(startTime, endTime, dt, solver, mover, fviewer, pviewer, writeFieldInterval, writeParticleInterval), m_iRank(rank)
{
  //m_iLastTimeIndex = getLastTimeIndex(m_fStartTime, m_fEndTime, m_dt);
}



int MPITimeController::solve() {
  double start, end;
  double total_start, total_end;

  total_start = omp_get_wtime();

  if(m_pFieldViewer != 0) {
    m_pFieldViewer->writeResult(0, 0, m_iRank);
  }

  if(m_pParticleViewer != 0) {
    m_pParticleViewer->writeResult(0, 0, m_iRank);
  }


  
  for(size_t i=1 ; i<=m_iLastTimeIndex ; i++) {
    double time = m_dt*i;

    //clock_t start,stop;
    //start = clock();
    start = omp_get_wtime();
    m_pSolver->solve(time);
    //stop = clock();
    end = omp_get_wtime();
    cout << i << "-th loop of processor No. " << m_iRank << " : Solving Time is   " << end-start << "   seconds." << endl;

    start = omp_get_wtime();
    if(m_pMover != 0) {
      if(int err = m_pMover->moveParticles(time, i)) {
        EM_ERROR("Moving Particles Error");
        return err;
      }
    }
    end = omp_get_wtime();
    cout << i << "-th loop of processor No. " << m_iRank << " : Particle Moving Time is " << end-start << " seconds." << endl;

    
    //*
    if((m_pFieldViewer != 0) && (i % m_iWriteFieldInterval == 0)) {
      double wstart = omp_get_wtime();
      
      m_pFieldViewer->writeResult(time, static_cast<int>(i), m_iRank);
      //fs.close();

      double wend = omp_get_wtime();
      cout << i << "-th loop of processor No. " << m_iRank << " : FieldViewer Saving Time is   " << wend-wstart << "   seconds." << endl;
    }
    //*/

    if((m_pParticleViewer != 0) && (i % m_iWriteParticleInterval == 0)) {
      double wstart = omp_get_wtime();
      m_pParticleViewer->writeResult(time, static_cast<int>(i), m_iRank);
      double wend = omp_get_wtime();
      cout << "i : "<< i << " : ParticleViewer Saving Time is " << wend-wstart << " seconds." << endl;
    }

    /*
    double Ex = 0.0;
    double Ey = 0.0;
    double Ez = 0.0;
    double Hx = 0.0;
    double Hy = 0.0;
    double Hz = 0.0;


    m_pSolver->evaluateError(i*m_dt, Ex, Ey, Ez, Hx, Hy, Hz, 1e-20);
    cout << "Timestep : " << i << " : Processor No. " << m_iRank << " : Error : Ex : " << Ex << endl;
    cout << "Timestep : " << i << " : Processor No. " << m_iRank << " : Error : Ey : " << Ey << endl;
    cout << "Timestep : " << i << " : Processor No. " << m_iRank << " : Error : Ez : " << Ez << endl;
    cout << "Timestep : " << i << " : Processor No. " << m_iRank << " : Error : Hx : " << Hx << endl;
    cout << "Timestep : " << i << " : Processor No. " << m_iRank << " : Error : Hy : " << Hy << endl;
    cout << "Timestep : " << i << " : Processor No. " << m_iRank << " : Error : Hz : " << Hz << endl;
    //*/
    
  }
  total_end = omp_get_wtime();
  cout << "Processor No. " << m_iRank << " : Solving Time is   " << total_end-total_start << "   seconds." << endl;





  /*
  double Ex = 0.0;
  double Ey = 0.0;
  double Ez = 0.0;
  double Hx = 0.0;
  double Hy = 0.0;
  double Hz = 0.0;


  m_pSolver->evaluateError(m_iLastTimeIndex*m_dt, Ex, Ey, Ez, Hx, Hy, Hz, 1e-20);
  cout << "Processor No. " << m_iRank << " : Error : Ex : " << Ex << endl;
  cout << "Processor No. " << m_iRank << " : Error : Ey : " << Ey << endl;
  cout << "Processor No. " << m_iRank << " : Error : Ez : " << Ez << endl;
  cout << "Processor No. " << m_iRank << " : Error : Hx : " << Hx << endl;
  cout << "Processor No. " << m_iRank << " : Error : Hy : " << Hy << endl;
  cout << "Processor No. " << m_iRank << " : Error : Hz : " << Hz << endl;
  //*/

  return 0;
  
}

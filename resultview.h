/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 13 2012
 * 
 *      
 */


#ifndef __RESULTVIEW_H__
#define __RESULTVIEW_H__


#include <ostream>

#include "solver.h"
#include "particle_mover.h"


class FieldViewer {

public:
  virtual void setFieldSolver(FieldSolver* solver) = 0;
  virtual FieldSolver* getFieldSolver() = 0;
  virtual void setMode(int mode) = 0;
  virtual int getMode() = 0;
  virtual void setOutputName(const string& name) = 0;
  virtual string getOutputName() = 0;
  virtual void setPrecision(size_t p) = 0;

  virtual int writeResult(double time, int timestep, int proc = -1) = 0;
};




class VTKFieldViewer : public FieldViewer {
public:
  enum SCHEME {SCALAR=0, VECTOR_SEPARATED=1, VECTOR_MERGED=2};

private:
  FieldSolver* m_pSolver;
  int m_iMode;
  size_t m_iPrecision;
  string m_sOutputName;

  double m_fLowerX; // for visualization. i.e., visualization boundary.
  double m_fUpperX; // for visualization. i.e., visualization boundary.
  double m_fLowerY; // for visualization. i.e., visualization boundary.
  double m_fUpperY; // for visualization. i.e., visualization boundary.
  double m_fLowerZ; // for visualization. i.e., visualization boundary.
  double m_fUpperZ; // for visualization. i.e., visualization boundary.

public:
  VTKFieldViewer(FieldSolver* solver, string outputname, int precision, int mode = 0);
  VTKFieldViewer(FieldSolver* solver, string outputname, double lowerX, double upperX, double lowerY, double upperY, double lowerZ, double upperZ, int precision, int mode = 0);


public:
  virtual void setFieldSolver(FieldSolver* solver) {m_pSolver = solver;}
  virtual FieldSolver* getFieldSolver() {return m_pSolver;}
  virtual void setMode(int mode) {m_iMode = mode;}
  virtual int getMode() {return m_iMode;}
  virtual void setOutputName(const string& name) {m_sOutputName = name;}
  virtual string getOutputName() {return m_sOutputName;}
  virtual void setPrecision(size_t p) {m_iPrecision = p;}

  virtual int writeResult(double time, int timestep, int proc = -1);
};






////////////////////////////////////////////////////////////////////////////////////////////





class ParticleViewer {


public:
  virtual void setParticleMover(ParticleMover* mover) = 0;
  virtual ParticleMover* getParticleSolver() = 0;
  //virtual void setMode(int mode) = 0;
  //virtual int getMode() = 0;
  virtual void setOutputName(const string& name) = 0;
  virtual string getOutputName() = 0;
  virtual void setPrecision(size_t p) = 0;

  virtual int writeResult(double time, int timestep, int proc = -1) = 0;
};




class VTKParticleViewer : public ParticleViewer {
public:
  //enum SCHEME {SCALAR=0, VECTOR_SEPARATED=1, VECTOR_MERGED=2};

private:
  ParticleMover* m_pMover;
  //int m_iMode;
  size_t m_iPrecision;
  string m_sOutputName;

  double m_fLowerX; // for visualization. i.e., visualization boundary.
  double m_fUpperX; // for visualization. i.e., visualization boundary.
  double m_fLowerY; // for visualization. i.e., visualization boundary.
  double m_fUpperY; // for visualization. i.e., visualization boundary.
  double m_fLowerZ; // for visualization. i.e., visualization boundary.
  double m_fUpperZ; // for visualization. i.e., visualization boundary.


public:
  //VTKParticleViewer(ParticleMover* mover, string outputname, int mode = 0);
  VTKParticleViewer(ParticleMover* mover, string outputname, double lowerX, double upperX, double lowerY, double upperY, double lowerZ, double upperZ, int precision);


public:
  virtual void setParticleMover(ParticleMover* mover) {m_pMover = mover;}
  virtual ParticleMover* getParticleSolver() {return m_pMover;}
  //virtual void setMode(int mode) {m_iMode = mode;}
  //virtual int getMode() {return m_iMode;}
  virtual void setOutputName(const string& name) {m_sOutputName = name;}
  virtual string getOutputName() {return m_sOutputName;}
  virtual void setPrecision(size_t p) {m_iPrecision = p;}

  virtual int writeResult(double time, int timestep, int proc = -1);

};








#endif

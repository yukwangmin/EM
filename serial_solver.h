/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 11 2012
 * 
 *      
 */


#ifndef __SERIAL_SOLVER_H__
#define __SERIAL_SOLVER_H__

#include "solver.h"


#define CELL_LOWER_X 0x00000001
#define CELL_UPPER_X 0x00000010
#define CELL_LOWER_Y 0x00000100
#define CELL_UPPER_Y 0x00001000
#define CELL_LOWER_Z 0x00010000
#define CELL_UPPER_Z 0x00100000




class SerialFieldSolver : public FieldSolver {
protected:
  size_t m_iDim;

  double m_dt;
  int m_iErrorMode;

  double m_fLowerX;
  double m_fUpperX;
  double m_dx;
  int m_iGridSizeX;
  double m_fLowerY;
  double m_fUpperY;
  double m_dy;
  int m_iGridSizeY;
  double m_fLowerZ;
  double m_fUpperZ;
  double m_dz;
  int m_iGridSizeZ;

  int m_iTotalDomainMemorySize;

  vector<double> *m_vCurrentEx;
  vector<double> *m_vCurrentEy;
  vector<double> *m_vCurrentEz;
  vector<double> *m_vCurrentHx;
  vector<double> *m_vCurrentHy;
  vector<double> *m_vCurrentHz;
  vector<double> *m_vNewEx;
  vector<double> *m_vNewEy;
  vector<double> *m_vNewEz;
  vector<double> *m_vNewHx;
  vector<double> *m_vNewHy;
  vector<double> *m_vNewHz;

  vector<double>* m_vInitialE;
  vector<double>* m_vInitialH;

//For Implementation
private:
  vector<double> *m_vMuX;
  vector<double> *m_vMuY;
  vector<double> *m_vMuZ;
  vector<double> *m_vEpsilonX;
  vector<double> *m_vEpsilonY;
  vector<double> *m_vEpsilonZ;
  vector<double> *m_vSigmaX;
  vector<double> *m_vSigmaY;
  vector<double> *m_vSigmaZ;



public:
  //SerialFieldSolver(double dt, double leftX, double rightX, double dx, double leftY, double rightY, double dy, double leftZ, double rightZ, double dz);
  SerialFieldSolver(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ);
  virtual ~SerialFieldSolver();

public: // interface.

  void setdt(double dt) {m_dt = dt;}
  double getdt() {return m_dt;}

  virtual void setErrorMode(int mode) {m_iErrorMode = mode;}
  virtual int getErrorMode() {return m_iErrorMode;}

  virtual double getLowerX() {return m_fLowerX;}
  virtual double getUpperX() {return m_fUpperX;}
  virtual double getdx() {return m_dx;}
  virtual int getGridSizeX() {return m_iGridSizeX;}

  virtual double getLowerY() {return m_fLowerY;}
  virtual double getUpperY() {return m_fUpperY;}
  virtual double getdy() {return m_dy;}
  virtual int getGridSizeY() {return m_iGridSizeY;}

  virtual double getLowerZ() {return m_fLowerZ;}
  virtual double getUpperZ() {return m_fUpperZ;}
  virtual double getdz() {return m_dz;}
  virtual int getGridSizeZ() {return m_iGridSizeZ;}


  /*!
   * This function must be called explicitly. 
   * This function initializes grid size and computational domain. Also this initializes memory saving data and load initial values on the data.
   */
  virtual void initializeSolver();

  virtual int solve(double time);


  virtual void getCurrentEx(vector<double>** result) {*result = m_vCurrentEx;}
  virtual void getCurrentEy(vector<double>** result) {*result = m_vCurrentEy;}
  virtual void getCurrentEz(vector<double>** result) {*result = m_vCurrentEz;}

  virtual void getCurrentHx(vector<double>** result) {*result = m_vCurrentHx;}
  virtual void getCurrentHy(vector<double>** result) {*result = m_vCurrentHy;}
  virtual void getCurrentHz(vector<double>** result) {*result = m_vCurrentHz;}

  // return H values at n by averaging (n-1/2 + n+1/2) / 2
  //virtual void getMagneticIntensityAtExactTime(vector<double>& hx, vector<double>& hy, vector<double>& hz);



  //BC
  //*
  virtual double getBoundaryXLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t) {return 0.0;}
  //*/


  virtual double getInitialEx(double x, double y, double z);
  virtual double getInitialEy(double x, double y, double z);
  virtual double getInitialEz(double x, double y, double z);

  virtual double getInitialHx(double x, double y, double z);
  virtual double getInitialHy(double x, double y, double z);
  virtual double getInitialHz(double x, double y, double z);

  virtual double getMu(double x, double y, double z, double time);
  virtual double getEpsilon(double x, double y, double z, double time);
  virtual double getSigma(double x, double y, double z, double time);

  virtual int updateCurrent(double time, int cellX, int cellY, int cellZ, double jx_y0z0, double jx_y1z0, double jx_y0z1, double jx_y1z1, double jy_x0z0, double jy_x1z0, double jy_x0z1, double jy_x1z1, double jz_x0y0, double jz_x1y0, double jz_x0y1, double jz_x1y1);



public: // Methods of this class
  /**
   * @return If the given position is not included in this computation domain, -1 will be returned.
   */
  virtual int getElectricIntensityAt(double time, size_t timestep, double x, double y, double z, double& ex, double& ey, double& ez);

  /**
   * @return If the given position is not included in this computation domain, -1 will be returned.
   */
  virtual int getMagneticDensityAt(double time, size_t timestep, double x, double y, double z, double& bx, double& by, double& bz);

  virtual int findCellHavingPosition(double x, double y, double z, int& cellX, int& cellY, int& cellZ, double& dx, double& dy, double& dz);

  


private:
  
  void updateE();
  void updateH();
  void loadBoundaryValues(double t);


  



public:
  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHz(double x, double y, double z, double t) {return 0.0;}

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {}
};






class SerialTestSolver : public SerialFieldSolver {
public:
  static const double mu;
  static const double epsilon;
  static const double frequency;
  static const double omega;

private:
  double m_fPIa;
  double m_fPIb;
  double m_fBeta;
  double m_fH;
  double m_fHsquare;

public:
  SerialTestSolver(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ);

public:

  virtual void initializeSolver();



  virtual double getBoundaryXLowerEx(double x, double y, double z, double t);
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t);
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t);
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t);
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t);
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t);
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t);
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t);
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t);
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t);
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t);
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t);
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t);
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t);
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t);
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t);
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t);
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t);

  
  virtual double getInitialEx(double x, double y, double z);
  virtual double getInitialEy(double x, double y, double z);
  virtual double getInitialEz(double x, double y, double z);

  virtual double getInitialHx(double x, double y, double z);
  virtual double getInitialHy(double x, double y, double z);
  virtual double getInitialHz(double x, double y, double z);



  virtual double getMu(double x, double y, double z, double time) {return mu;}
  virtual double getEpsilon(double x, double y, double z, double time) {return epsilon;}
  virtual double getSigma(double x, double y, double z, double time) {return 0.0;}

private:
  


public:
  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t);
  virtual double getExactEy(double x, double y, double z, double t);
  virtual double getExactEz(double x, double y, double z, double t);
  virtual double getExactHx(double x, double y, double z, double t);
  virtual double getExactHy(double x, double y, double z, double t);
  virtual double getExactHz(double x, double y, double z, double t);

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold);
};



class SerialNullFieldSolver : public SerialFieldSolver {
public:
  //static const double mu;
  //static const double epsilon;
  //static const double frequency;
  //static const double omega;

private:
  double m_fPIa;
  double m_fPIb;
  double m_fBeta;
  double m_fH;
  double m_fHsquare;

public:
  SerialNullFieldSolver(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
   : SerialFieldSolver(dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ) {}




  virtual double getBoundaryXLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t) {return 0.0;}

  
  virtual double getInitialEx(double x, double y, double z) {return 0.0;}
  virtual double getInitialEy(double x, double y, double z) {return 0.0;}
  virtual double getInitialEz(double x, double y, double z) {return 0.0;}

  virtual double getInitialHx(double x, double y, double z) {return 0.0;}
  virtual double getInitialHy(double x, double y, double z) {return 0.0;}
  virtual double getInitialHz(double x, double y, double z) {return 0.0;}



  virtual double getMu(double x, double y, double z, double time) {return 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;}
  virtual double getEpsilon(double x, double y, double z, double time) {return 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;}
  virtual double getSigma(double x, double y, double z, double time) {return 0.0;}

private:
  


public:
  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHz(double x, double y, double z, double t) {return 0.0;}

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {}
};





// Electric Static Only - Z direction
class SerialTestSolverWithParticle1 : public SerialFieldSolver {
public:
  //static const double mu;
  //static const double epsilon;
  //static const double frequency;
  //static const double omega;

private:
  double m_fPIa;
  double m_fPIb;
  double m_fBeta;
  double m_fH;
  double m_fHsquare;

public:
  SerialTestSolverWithParticle1(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
   : SerialFieldSolver(dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ) {}




  virtual double getBoundaryXLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t) {return 10.0;}

  
  virtual double getInitialEx(double x, double y, double z) {return 0.0;}
  virtual double getInitialEy(double x, double y, double z) {return 0.0;}
  virtual double getInitialEz(double x, double y, double z) {return 10.0;}

  virtual double getInitialHx(double x, double y, double z) {return 0.0;}
  virtual double getInitialHy(double x, double y, double z) {return 0.0;}
  virtual double getInitialHz(double x, double y, double z) {return 0.0;}



  virtual double getMu(double x, double y, double z, double time) {return 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;}
  virtual double getEpsilon(double x, double y, double z, double time) {return 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;}
  virtual double getSigma(double x, double y, double z, double time) {return 0.0;}

private:
  


public:
  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHz(double x, double y, double z, double t) {return 0.0;}

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {}
};







// Electric Static Only - X direction
class SerialTestSolverWithParticle2 : public SerialFieldSolver {
public:
  //static const double mu;
  //static const double epsilon;
  //static const double frequency;
  //static const double omega;

private:
  double m_fPIa;
  double m_fPIb;
  double m_fBeta;
  double m_fH;
  double m_fHsquare;

public:
  SerialTestSolverWithParticle2(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
   : SerialFieldSolver(dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ) {}




  virtual double getBoundaryXLowerEx(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t) {return 0.0;}

  
  virtual double getInitialEx(double x, double y, double z) {return 10.0;}
  virtual double getInitialEy(double x, double y, double z) {return 0.0;}
  virtual double getInitialEz(double x, double y, double z) {return 0.0;}

  virtual double getInitialHx(double x, double y, double z) {return 0.0;}
  virtual double getInitialHy(double x, double y, double z) {return 0.0;}
  virtual double getInitialHz(double x, double y, double z) {return 0.0;}



  virtual double getMu(double x, double y, double z, double time) {return 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;}
  virtual double getEpsilon(double x, double y, double z, double time) {return 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;}
  virtual double getSigma(double x, double y, double z, double time) {return 0.0;}

private:
  


public:
  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHz(double x, double y, double z, double t) {return 0.0;}

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {}
};







// Electric Static Only - Y direction
class SerialTestSolverWithParticle3 : public SerialFieldSolver {
public:
  //static const double mu;
  //static const double epsilon;
  //static const double frequency;
  //static const double omega;

private:
  double m_fPIa;
  double m_fPIb;
  double m_fBeta;
  double m_fH;
  double m_fHsquare;

public:
  SerialTestSolverWithParticle3(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
   : SerialFieldSolver(dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ) {}




  virtual double getBoundaryXLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t) {return 10.0;}
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t) {return 0.0;}

  
  virtual double getInitialEx(double x, double y, double z) {return 0.0;}
  virtual double getInitialEy(double x, double y, double z) {return 10.0;}
  virtual double getInitialEz(double x, double y, double z) {return 0.0;}

  virtual double getInitialHx(double x, double y, double z) {return 0.0;}
  virtual double getInitialHy(double x, double y, double z) {return 0.0;}
  virtual double getInitialHz(double x, double y, double z) {return 0.0;}



  virtual double getMu(double x, double y, double z, double time) {return 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;}
  virtual double getEpsilon(double x, double y, double z, double time) {return 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;}
  virtual double getSigma(double x, double y, double z, double time) {return 0.0;}

private:
  


public:
  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHz(double x, double y, double z, double t) {return 0.0;}

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {}
};





#endif
/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 11 2012
 * 
 *      
 */


#ifndef __PARALLEL_SOLVER_H__
#define __PARALLEL_SOLVER_H__

#include "solver.h"





struct BlockDomainField {
  double* ex;
  double* ey;
  double* ez;

  double* bx;
  double* by;
  double* bz;
};



class MPIFieldSolver : public FieldSolver {
protected:
  size_t m_iDim;
  int m_iTotalProcNum;
  int m_iProcNumX;
  int m_iProcNumY;
  int m_iProcNumZ;
  int m_iRank;

  int m_iLowerXNeighbor;
  int m_iUpperXNeighbor;
  int m_iLowerYNeighbor;
  int m_iUpperYNeighbor;
  int m_iLowerZNeighbor;
  int m_iUpperZNeighbor;


  double m_dt;
  int m_iErrorMode;


  double m_fGlobalLowerX;
  double m_fGlobalUpperX;
  double m_fGlobalLowerY;
  double m_fGlobalUpperY;
  double m_fGlobalLowerZ;
  double m_fGlobalUpperZ;

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


  /*
  vector<double> *m_vMuX;
  vector<double> *m_vMuY;
  vector<double> *m_vMuZ;
  vector<double> *m_vEpsilonX;
  vector<double> *m_vEpsilonY;
  vector<double> *m_vEpsilonZ;
  vector<double> *m_vSigmaX;
  vector<double> *m_vSigmaY;
  vector<double> *m_vSigmaZ;
  */
  double m_fMu;
  double m_fEpsilon;
  double m_fSigma;


//For Implementation
private:
  vector<BlockDomainField> m_vGolableFieldData;
  double* m_vSendBuffer;
  double* m_vRecvBuffer;

  vector<double> *m_vGlobalEx;
  vector<double> *m_vGlobalEy;
  vector<double> *m_vGlobalEz;
  vector<double> *m_vGlobalBx;
  vector<double> *m_vGlobalBy;
  vector<double> *m_vGlobalBz;

  int m_iGlobalGridSizeX;
  int m_iGlobalGridSizeY;
  int m_iGlobalGridSizeZ;


  /*
  double *m_aBuf_lower_x_Ey;
  double *m_aBuf_lower_x_Ez;
  double *m_aBuf_upper_x_Ey;
  double *m_aBuf_upper_x_Ez;
  double *m_aBuf_lower_y_Ex;
  double *m_aBuf_lower_y_Ez;
  double *m_aBuf_upper_y_Ex;
  double *m_aBuf_upper_y_Ez;
  double *m_aBuf_lower_z_Ex;
  double *m_aBuf_lower_z_Ey;
  double *m_aBuf_upper_z_Ex;
  double *m_aBuf_upper_z_Ey;

  int m_iBuf_size_x_Ey;
  int m_iBuf_size_x_Ez;
  int m_iBuf_size_y_Ex;
  int m_iBuf_size_y_Ez;
  int m_iBuf_size_z_Ex;
  int m_iBuf_size_z_Ey;

  double *m_aBuf_lower_x_Hx;
  double *m_aBuf_upper_x_Hx;
  double *m_aBuf_lower_y_Hy;
  double *m_aBuf_upper_y_Hy;
  double *m_aBuf_lower_z_Hz;
  double *m_aBuf_upper_z_Hz;

  int m_iBuf_size_x_Hx;
  int m_iBuf_size_y_Hy;
  int m_iBuf_size_z_Hz;
  */




public:
  MPIFieldSolver(int rank, int procNumX, int procNumY, int procNumZ, double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ);
  virtual ~MPIFieldSolver();

public:  // Solver class Overriding

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


  virtual double getInitialEx(double x, double y, double z) {return 0.0;}
  virtual double getInitialEy(double x, double y, double z) {return 0.0;}
  virtual double getInitialEz(double x, double y, double z) {return 0.0;}

  virtual double getInitialHx(double x, double y, double z) {return 0.0;}
  virtual double getInitialHy(double x, double y, double z) {return 0.0;}
  virtual double getInitialHz(double x, double y, double z) {return 0.0;}



  virtual double getMu(double x, double y, double z, double time);
  virtual double getEpsilon(double x, double y, double z, double time);
  virtual double getSigma(double x, double y, double z, double time);

  virtual int updateCurrent(double time, int cellX, int cellY, int cellZ, double jx_y0z0, double jx_y1z0, double jx_y0z1, double jx_y1z1, double jy_x0z0, double jy_x1z0, double jy_x0z1, double jy_x1z1, double jz_x0y0, double jz_x1y0, double jz_x0y1, double jz_x1y1);


public: // Methods of this class
  void getGlobalFieldData(size_t timestep);
  void updateGlobalElectricIntensityByCurrent(size_t timestep);

  /**
   * @return If the given position is not included in this computation domain, -1 will be returned.
   */
  virtual int getElectricIntensityAt(double time, size_t timestep, double x, double y, double z, double& ex, double& ey, double& ez);

  /**
   * @return If the given position is not included in this computation domain, -1 will be returned.
   */
  virtual int getMagneticDensityAt(double time, size_t timestep,  double x, double y, double z, double& bx, double& by, double& bz);

  virtual int findCellHavingPosition(double x, double y, double z, int& cellX, int& cellY, int& cellZ, double& dx, double& dy, double& dz);

  /*
  void getGlobalDomainInformation(double& lowerX, double& upperX, double& lowerY, double& upperY, double& lowerZ, double& upperZ) {lowerX = m_fGlobalLowerX; upperX = m_fGlobalUpperX; lowerY = m_fGlobalLowerY; upperY = m_fGlobalUpperY; lowerZ = m_fGlobalLowerZ; upperZ = m_fGlobalUpperZ;}
  void getProcInformation(int& rank, int& totalNumProc, int& numProcX, int& numProcY, int& numProcZ) {rank = m_iRank; totalNumProc = m_iTotalProcNum; numProcX = m_iProcNumX; numProcY = m_iProcNumY; numProcZ = m_iProcNumZ;}
  */
  double getGlobalLowerX() {return m_fGlobalLowerX;}
  double getGlobalUpperX() {return m_fGlobalUpperX;}
  double getGlobalLowerY() {return m_fGlobalLowerY;}
  double getGlobalUpperY() {return m_fGlobalUpperY;}
  double getGlobalLowerZ() {return m_fGlobalLowerZ;}
  double getGlobalUpperZ() {return m_fGlobalUpperZ;}
  
  void resetGlobalElectricIntensity();
  


private:
  
  void updateE();
  void updateH();

  void initGlobalFieldDomain();
  



public: // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactEz(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHx(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHy(double x, double y, double z, double t) {return 0.0;}
  virtual double getExactHz(double x, double y, double z, double t) {return 0.0;}

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {}



};




class MPITestSolver : public MPIFieldSolver {
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
  MPITestSolver(int rank, int procNumX, int procNumY, int procNumZ, double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ);



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


  //*
  virtual double getMu(double x, double y, double z, double time) {return mu;}
  virtual double getEpsilon(double x, double y, double z, double time) {return epsilon;}
  virtual double getSigma(double x, double y, double z, double time) {return 0.0;}
  //*/
  /*
  virtual double getMu() {return mu;}
  virtual double getEpsilon() {return epsilon;}
  virtual double getSigma() {return 0.0;}
  */


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





class MPINullFieldSolver : public MPIFieldSolver {
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
  MPINullFieldSolver(int rank, int procNumX, int procNumY, int procNumZ, double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
   : MPIFieldSolver(rank, procNumX, procNumY, procNumZ, dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ) {}



  /*
  virtual void initializeSolver() {
    m_fMu = getMu(0,0,0,0);
    m_fEpsilon = getEpsilon(0,0,0,0);
    m_fSigma = getSigma(0,0,0,0);
    
    MPIFieldSolver::initializeSolver();
  }
  */
  
  
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
class MPITestSolverWithParticle1 : public MPIFieldSolver {
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
  MPITestSolverWithParticle1(int rank, int procNumX, int procNumY, int procNumZ, double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
   : MPIFieldSolver(rank, procNumX, procNumY, procNumZ, dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ) {}



  /*
  virtual void initializeSolver() {
    m_fMu = getMu(0,0,0,0);
    m_fEpsilon = getEpsilon(0,0,0,0);
    m_fSigma = getSigma(0,0,0,0);
    
    MPIFieldSolver::initializeSolver();
  }
  */
  
  
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





#endif
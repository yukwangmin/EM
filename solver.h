/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jun. 05 2012
 * 
 *      
 */



#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <vector>


using namespace std;



class FieldSolver {


public:
  enum Face{LOWER_X = -1, UPPER_X = -2, LOWER_Y = -3, UPPER_Y = -4, LOWER_Z = -5, UPPER_Z = -6};

  virtual void setdt(double dt) = 0;
  virtual double getdt() = 0;

  virtual void setErrorMode(int mode) = 0;
  virtual int getErrorMode() = 0;



  virtual double getLowerX() = 0;
  virtual double getUpperX() = 0;
  virtual double getdx() = 0;
  virtual int getGridSizeX() = 0;

  virtual double getLowerY() = 0;
  virtual double getUpperY() = 0;
  virtual double getdy() = 0;
  virtual int getGridSizeY() = 0;

  virtual double getLowerZ() = 0;
  virtual double getUpperZ() = 0;
  virtual double getdz() = 0;
  virtual int getGridSizeZ() = 0;



  virtual void initializeSolver() = 0;

  virtual int solve(double time) = 0;



  virtual void getCurrentEx(vector<double>** result) = 0;
  virtual void getCurrentEy(vector<double>** result) = 0;
  virtual void getCurrentEz(vector<double>** result) = 0;

  virtual void getCurrentHx(vector<double>** result) = 0;
  virtual void getCurrentHy(vector<double>** result) = 0;
  virtual void getCurrentHz(vector<double>** result) = 0;

  // return H values at n by averaging (n-1/2 + n+1/2) / 2
  //virtual void getMagneticIntensityAtExactTime(vector<double>& hx, vector<double>& hy, vector<double>& hz) = 0;




  virtual double getBoundaryXLowerEx(double x, double y, double z, double t) = 0;
  virtual double getBoundaryXLowerEy(double x, double y, double z, double t) = 0;
  virtual double getBoundaryXLowerEz(double x, double y, double z, double t) = 0;
  virtual double getBoundaryXUpperEx(double x, double y, double z, double t) = 0;
  virtual double getBoundaryXUpperEy(double x, double y, double z, double t) = 0;
  virtual double getBoundaryXUpperEz(double x, double y, double z, double t) = 0;
  virtual double getBoundaryYLowerEx(double x, double y, double z, double t) = 0;
  virtual double getBoundaryYLowerEy(double x, double y, double z, double t) = 0;
  virtual double getBoundaryYLowerEz(double x, double y, double z, double t) = 0;
  virtual double getBoundaryYUpperEx(double x, double y, double z, double t) = 0;
  virtual double getBoundaryYUpperEy(double x, double y, double z, double t) = 0;
  virtual double getBoundaryYUpperEz(double x, double y, double z, double t) = 0;
  virtual double getBoundaryZLowerEx(double x, double y, double z, double t) = 0;
  virtual double getBoundaryZLowerEy(double x, double y, double z, double t) = 0;
  virtual double getBoundaryZLowerEz(double x, double y, double z, double t) = 0;
  virtual double getBoundaryZUpperEx(double x, double y, double z, double t) = 0;
  virtual double getBoundaryZUpperEy(double x, double y, double z, double t) = 0;
  virtual double getBoundaryZUpperEz(double x, double y, double z, double t) = 0;


  virtual double getInitialEx(double x, double y, double z) = 0;
  virtual double getInitialEy(double x, double y, double z) = 0;
  virtual double getInitialEz(double x, double y, double z) = 0;
  virtual double getInitialHx(double x, double y, double z) = 0;
  virtual double getInitialHy(double x, double y, double z) = 0;
  virtual double getInitialHz(double x, double y, double z) = 0;



  virtual double getMu(double x, double y, double z, double time) = 0;
  virtual double getEpsilon(double x, double y, double z, double time) = 0;
  virtual double getSigma(double x, double y, double z, double time) = 0;

  virtual int updateCurrent(double time, int cellX, int cellY, int cellZ, double jx_y0z0, double jx_y1z0, double jx_y0z1, double jx_y1z1, double jy_x0z0, double jy_x1z0, double jy_x0z1, double jy_x1z1, double jz_x0y0, double jz_x1y0, double jz_x0y1, double jz_x1y1) = 0;


  // Temporary for test.
  virtual double getExactEx(double x, double y, double z, double t) = 0;
  virtual double getExactEy(double x, double y, double z, double t) = 0;
  virtual double getExactEz(double x, double y, double z, double t) = 0;
  virtual double getExactHx(double x, double y, double z, double t) = 0;
  virtual double getExactHy(double x, double y, double z, double t) = 0;
  virtual double getExactHz(double x, double y, double z, double t) = 0;

  virtual void evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) = 0;

};








#endif
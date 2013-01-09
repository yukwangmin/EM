/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Thu Oct. 27 2012
 *
 */

#ifndef __USER_SPECIFIC_INITIAL_PARTICLE_SETTING_H__
#define __USER_SPECIFIC_INITIAL_PARTICLE_SETTING_H__


#include <stdio.h>
#include <vector>

/*
struct UserGrid
{
  //Spacing along each direction.
  size_t xGridSize;
  size_t yGridSize;
  size_t zGridSize;

  //Grid Bounding box corners.
  double xleft, xright;
  double yleft, yright;
  double zleft, zright;

  //For the binary interpolation.Usually set to some non zero number in ( 0, 1 ) This should be chosen empirically.
  double tolerance;


  std::vector<bool> level_set_function;

};
*/


class UserSpecificInitialParticleSetting {
public:
  virtual void getDomain(double& lowerX, double& lowerY, double& lowerZ, double& upperX, double& upperY, double& upperZ) = 0;
  virtual void getGridSize(size_t& x, size_t& y, size_t& z) = 0;
  virtual double getTolerance() = 0;

  /*!
   *
   * return Then number of vector size.
   *
   * Level set functions for geometry. 
   * 3 dimensional index.
   * size of level_set_function is (xGridSize + 1)*(yGridSize + 1)*(zGridSize + 1)
   * (i,j,k) index : i + (xGridSize+1)*j + (xGridSize+1)*(yGridSize+1)*k
   */
  virtual size_t fillGrid(std::vector<bool>& level_set) = 0;
};


class DefaultUserSpecificInitialParticleSetting : public UserSpecificInitialParticleSetting {
  virtual void getDomain(double& lowerX, double& lowerY, double& lowerZ, double& upperX, double& upperY, double& upperZ) {
    lowerX = 0.0;
    lowerY = 0.0;
    lowerZ = 0.0;

    upperX = 0.0;
    upperY = 0.0;
    upperZ = 0.0;
  }

  virtual void getGridSize(size_t& x, size_t& y, size_t& z) {x = 1; y = 1; z = 1;}
  virtual double getTolerance() {return 0.0;}

  virtual size_t fillGrid(std::vector<bool>& level_set) {return 0;}
};




class TestParticleSetting1 : public UserSpecificInitialParticleSetting {
  virtual void getDomain(double& lowerX, double& lowerY, double& lowerZ, double& upperX, double& upperY, double& upperZ) {
    lowerX = 0.064;
    lowerY = 0.064;
    lowerZ = 0.512;

    upperX = 0.192;
    upperY = 0.192;
    upperZ = 0.64;
  }

  virtual void getGridSize(size_t& x, size_t& y, size_t& z) {x = 10; y = 10; z = 10;}

  virtual double getTolerance() {return 0.5;}

  virtual size_t fillGrid(std::vector<bool>& level_set);
};






#endif //__USER_SPECIFIC_INITIAL_PARTICLE_SETTING_H__ 
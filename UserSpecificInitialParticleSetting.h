/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 11 2012
 * 
 *      
 */


#ifndef __USERSPECIFICINITIALPARTICLESETTING_H__
#define __USERSPECIFICINITIALPARTICLESETTING_H__


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


 
#endif


/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Thu Oct. 27 2012
 *
 */

#include <cmath>

#include "user_specific_initial_particle_setting.h"


size_t TestParticleSetting1::fillGrid(std::vector<bool>& level_set) {

  size_t xgridsize, ygridsize, zgridsize;
  double xleft, yleft, zleft, xright, yright, zright;
  double dx, dy, dz;

  getGridSize(xgridsize, ygridsize, zgridsize);
  getDomain(xleft, yleft, zleft, xright, yright, zright);

  dx = (xright-xleft)/xgridsize;
  dy = (yright-yleft)/ygridsize;
  dz = (zright-zleft)/zgridsize;

  level_set.clear();

  size_t size = (xgridsize+1)*(ygridsize+1)*(zgridsize+1);

  level_set.resize(size);

  //Now fill up the grid.
  size_t index;
  for (size_t k = 0; k < zgridsize; ++k) {
    for (size_t j = 0; j < ygridsize; ++j) {
	    for (size_t i = 0; i < xgridsize; ++i) {
        index = i + j*xgridsize + k*xgridsize*ygridsize;
        //A diamond shaped region.
        //double temp = fabs( xleft + i*dx - 0.128) + fabs( yleft + j*dy - 0.128) + fabs( zleft + k*dz - 0.128);
        if ( fabs( xleft + i*dx - 0.128) + fabs( yleft + j*dy - 0.128) + fabs( zleft + k*dz - 0.576) <= 0.128 ) {
	     	  level_set[index] = true;
	     	} else {
	     	  level_set[index] = false;
	     	}     
         
	    }
	  }
  }

  return size;
}
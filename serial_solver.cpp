/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 11 2012
 * 
 *      
 */



#include <iostream>
#include <omp.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>



#include "serial_solver.h"
#include "em_error.h"
#include "constants.h"


using namespace std;






////////////////////////////////////////////////////////////////////////////////
// Start : SerialFieldSolver
////////////////////////////////////////////////////////////////////////////////

/*
SerialFieldSolver::SerialFieldSolver(double dt, double leftX, double rightX, double dx, double leftY, double rightY, double dy, double leftZ, double rightZ, double dz) 
: m_dt(dt), m_fLowerX(leftX), m_fUpperX(rightX), m_dx(dx), m_fLowerY(leftY), m_fUpperY(rightY), m_dy(dy), m_fLowerZ(leftZ), m_fUpperZ(rightZ), m_dz(dz), m_iErrorMode(0)
{
}
*/


SerialFieldSolver::SerialFieldSolver(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
: m_dt(dt), m_fLowerX(leftX), m_fUpperX(rightX), m_iGridSizeX(gridNumX+1), m_fLowerY(leftY), m_fUpperY(rightY), m_iGridSizeY(gridNumY+1), m_fLowerZ(leftZ), m_fUpperZ(rightZ), m_iGridSizeZ(gridNumZ+1), m_iErrorMode(0)
{
  m_dx = (rightX - leftX) / gridNumX;
  m_dy = (rightY - leftY) / gridNumY;
  m_dz = (rightZ - leftZ) / gridNumZ;
}


SerialFieldSolver::~SerialFieldSolver() {
  delete m_vCurrentEx;
  delete m_vCurrentEy;
  delete m_vCurrentEz;
  delete m_vCurrentHx;
  delete m_vCurrentHy;
  delete m_vCurrentHz;
  delete m_vNewEx;
  delete m_vNewEy;
  delete m_vNewEz;
  delete m_vNewHx;
  delete m_vNewHy;
  delete m_vNewHz;

  /*
  delete m_vBoundaryLowerX;
  delete m_vBoundaryUpperX;
  delete m_vBoundaryLowerY;
  delete m_vBoundaryUpperY;
  delete m_vBoundaryLowerZ;
  delete m_vBoundaryUpperZ;
  */

  delete m_vMuX;
  delete m_vMuY;
  delete m_vMuZ;
  delete m_vEpsilonX;
  delete m_vEpsilonY;
  delete m_vEpsilonZ;
  delete m_vSigmaX;
  delete m_vSigmaY;
  delete m_vSigmaZ;
}


void SerialFieldSolver::initializeSolver() {
  /*
  double temp;
  if(m_fLowerX > m_fUpperX) {
    temp = m_fLowerX;
    m_fLowerX = m_fUpperX;
    m_fUpperX = temp;
  }
  if(m_fLowerY > m_fUpperY) {
    temp = m_fLowerY;
    m_fLowerY = m_fUpperY;
    m_fUpperY = temp;
  }
  if(m_fLowerZ > m_fUpperZ) {
    temp = m_fLowerZ;
    m_fLowerZ = m_fUpperZ;
    m_fUpperZ = temp;
  }
  */


  //m_iGridSizeX = static_cast<int>(((m_fUpperX - m_fLowerX)/m_dx) + m_dx) + 1;
  //m_iGridSizeY = static_cast<int>(((m_fUpperY - m_fLowerY)/m_dy) + m_dy) + 1;
  //m_iGridSizeZ = static_cast<int>(((m_fUpperZ - m_fLowerZ)/m_dz) + m_dz) + 1;
  m_iTotalDomainMemorySize = m_iGridSizeX*m_iGridSizeY*m_iGridSizeZ;

  m_vCurrentEx = new vector<double>(m_iTotalDomainMemorySize);
  m_vCurrentEy = new vector<double>(m_iTotalDomainMemorySize);
  m_vCurrentEz = new vector<double>(m_iTotalDomainMemorySize);
  m_vCurrentHx = new vector<double>(m_iTotalDomainMemorySize);
  m_vCurrentHy = new vector<double>(m_iTotalDomainMemorySize);
  m_vCurrentHz = new vector<double>(m_iTotalDomainMemorySize);
  m_vNewEx = new vector<double>(m_iTotalDomainMemorySize);
  m_vNewEy = new vector<double>(m_iTotalDomainMemorySize);
  m_vNewEz = new vector<double>(m_iTotalDomainMemorySize);
  m_vNewHx = new vector<double>(m_iTotalDomainMemorySize);
  m_vNewHy = new vector<double>(m_iTotalDomainMemorySize);
  m_vNewHz = new vector<double>(m_iTotalDomainMemorySize);

  /*
  m_vBoundaryLowerX = new vector<double>(m_iGridSizeY*m_iGridSizeZ);
  m_vBoundaryUpperX = new vector<double>(m_iGridSizeY*m_iGridSizeZ);
  m_vBoundaryLowerY = new vector<double>(m_iGridSizeX*m_iGridSizeZ);
  m_vBoundaryUpperY = new vector<double>(m_iGridSizeX*m_iGridSizeZ);
  m_vBoundaryLowerZ = new vector<double>(m_iGridSizeX*m_iGridSizeY);
  m_vBoundaryUpperZ = new vector<double>(m_iGridSizeX*m_iGridSizeY);
  */

  m_vMuX = new vector<double>(m_iTotalDomainMemorySize);
  m_vMuY = new vector<double>(m_iTotalDomainMemorySize);
  m_vMuZ = new vector<double>(m_iTotalDomainMemorySize);
  m_vEpsilonX = new vector<double>(m_iTotalDomainMemorySize);
  m_vEpsilonY = new vector<double>(m_iTotalDomainMemorySize);
  m_vEpsilonZ = new vector<double>(m_iTotalDomainMemorySize);
  m_vSigmaX = new vector<double>(m_iTotalDomainMemorySize);
  m_vSigmaY = new vector<double>(m_iTotalDomainMemorySize);
  m_vSigmaZ = new vector<double>(m_iTotalDomainMemorySize);

  size_t index;

  // Initialize mu, epsilon, sigma
  for(int i=0 ; i<m_iGridSizeZ ; i++) {
    for(int j=0 ; j<m_iGridSizeY ; j++) {
      for(int k=0 ; k<m_iGridSizeX ; k++) {
        index = k + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*i;
        (*m_vMuX)[index] = getMu(m_fLowerX + k*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + i*m_dz + 0.5*m_dz, 0);
        (*m_vMuY)[index] = getMu(m_fLowerX + k*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + i*m_dz + 0.5*m_dz, 0);
        (*m_vMuZ)[index] = getMu(m_fLowerX + k*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + i*m_dz, 0);
        (*m_vEpsilonX)[index] = getEpsilon(m_fLowerX + k*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + i*m_dz, 0);
        (*m_vEpsilonY)[index] = getEpsilon(m_fLowerX + k*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + i*m_dz, 0);
        (*m_vEpsilonZ)[index] = getEpsilon(m_fLowerX + k*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + i*m_dz + 0.5*m_dz, 0);
        (*m_vSigmaX)[index] = getSigma(m_fLowerX + k*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + i*m_dz, 0);
        (*m_vSigmaY)[index] = getSigma(m_fLowerX + k*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + i*m_dz, 0);
        (*m_vSigmaZ)[index] = getSigma(m_fLowerX + k*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + i*m_dz + 0.5*m_dz, 0);
      }
    }
  }



  ////////////////////////////////////////////////////////////////////////////////
  // Start : loading Initial values
  ////////////////////////////////////////////////////////////////////////////////

  // For E.
  for(int k=0 ; k<m_iGridSizeZ ; k++) {
    for(int j=0 ; j<m_iGridSizeY ; j++) {
      for(int i=0 ; i<m_iGridSizeX ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vCurrentEx)[index] = getInitialEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz);
        (*m_vCurrentEy)[index] = getInitialEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz);
        (*m_vCurrentEz)[index] = getInitialEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz);
      }
    }
  }
  
  // For H.
  for(int k=0 ; k<m_iGridSizeZ ; k++) {
    for(int j=0 ; j<m_iGridSizeY ; j++) {
      for(int i=0 ; i<m_iGridSizeX ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vCurrentHx)[index] = getInitialHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz);
        (*m_vCurrentHy)[index] = getInitialHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz);
        (*m_vCurrentHz)[index] = getInitialHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz);
      }
    }
  }

  // Calculating half time step of H.
  //#pragma omp parallel for private(i,j,k,index)
  int i,j,k;
  for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=0 ; j<m_iGridSizeY-1 ; j++) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = (*m_vCurrentHx)[index] + (0.5*m_dt/(*m_vMuX)[index]) * ( ((*m_vCurrentEy)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEy)[index])/m_dz - ((*m_vCurrentEz)[index+m_iGridSizeX] - (*m_vCurrentEz)[index])/m_dy );
        (*m_vNewHy)[index] = (*m_vCurrentHy)[index] + (0.5*m_dt/(*m_vMuY)[index]) * ( ((*m_vCurrentEz)[index+1] - (*m_vCurrentEz)[index])/m_dx - ((*m_vCurrentEx)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEx)[index])/m_dz );
        (*m_vNewHz)[index] = (*m_vCurrentHz)[index] + (0.5*m_dt/(*m_vMuZ)[index]) * ( ((*m_vCurrentEx)[index+m_iGridSizeX] - (*m_vCurrentEx)[index])/m_dy - ((*m_vCurrentEy)[index+1] - (*m_vCurrentEy)[index])/m_dx );
      }
    }
  }

  updateH();


  ////////////////////////////////////////////////////////////////////////////////
  // End : loading Initial values
  ////////////////////////////////////////////////////////////////////////////////


  loadBoundaryValues(0);




}


int SerialFieldSolver::solve(double time) {

  size_t index;

  /*
#ifndef _MSC_VER
  size_t i,j,k;
#else
  int i,j,k;
#endif 
  */

    /*
    {
      ofstream fsEx("timestep_0_Ex.txt");
      ofstream fsEy("timestep_0_Ey.txt");
      ofstream fsEz("timestep_0_Ez.txt");
  
      for(int k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(int j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(int i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
            fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEx->at(index) << endl;
            fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEy->at(index) << endl;
            fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEz->at(index) << endl;
          }
        }
      }
      fsEx.close();
      fsEy.close();
      fsEz.close();
    }
    //*/

  int i,j,k;


  //double v1, v2, v3, v4, v5, v6, v7, v8;
  // Solving E.
  #pragma omp parallel for private(i,j,k,index)
  for(k=1 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=1 ; j<m_iGridSizeY-1 ; j++) {
      for(i=1 ; i<m_iGridSizeX-1 ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEx)[index] = (*m_vCurrentEx)[index] + (m_dt/(*m_vEpsilonX)[index]) * ( ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-m_iGridSizeX])/m_dy - ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-m_iGridSizeX*m_iGridSizeY])/m_dz ) + (*m_vSigmaX)[index]*(*m_vCurrentEx)[index];
        (*m_vNewEy)[index] = (*m_vCurrentEy)[index] + (m_dt/(*m_vEpsilonY)[index]) * ( ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX*m_iGridSizeY])/m_dz - ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-1])/m_dx ) + (*m_vSigmaY)[index]*(*m_vCurrentEy)[index];
        (*m_vNewEz)[index] = (*m_vCurrentEz)[index] + (m_dt/(*m_vEpsilonZ)[index]) * ( ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-1])/m_dx - ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX])/m_dy ) + (*m_vSigmaZ)[index]*(*m_vCurrentEz)[index];
      }
    }
  }
  #pragma omp parallel for private(j,k,index)
  for(k=1 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=1 ; j<m_iGridSizeY-1 ; j++) {
      index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vNewEx)[index] = (*m_vCurrentEx)[index] + (m_dt/(*m_vEpsilonX)[index]) * ( ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-m_iGridSizeX])/m_dy - ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-m_iGridSizeX*m_iGridSizeY])/m_dz ) + (*m_vSigmaX)[index]*(*m_vCurrentEx)[index];
    }
  }
  #pragma omp parallel for private(i,k,index)
  for(k=1 ; k<m_iGridSizeZ-1 ; k++) {
    for(i=1 ; i<m_iGridSizeX-1 ; i++) {
      index = i + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vNewEy)[index] = (*m_vCurrentEy)[index] + (m_dt/(*m_vEpsilonY)[index]) * ( ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX*m_iGridSizeY])/m_dz - ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-1])/m_dx ) + (*m_vSigmaY)[index]*(*m_vCurrentEy)[index];
    }
  }
  #pragma omp parallel for private(i,j,index)
  for(j=1 ; j<m_iGridSizeY-1 ; j++) {
    for(i=1 ; i<m_iGridSizeX-1 ; i++) {
      index = i + m_iGridSizeX*j;
      (*m_vNewEz)[index] = (*m_vCurrentEz)[index] + (m_dt/(*m_vEpsilonZ)[index]) * ( ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-1])/m_dx - ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX])/m_dy ) + (*m_vSigmaZ)[index]*(*m_vCurrentEz)[index];
    }
  }

  updateE();
  loadBoundaryValues(time);

  
  // Solving H.
  #pragma omp parallel for private(i,j,k,index)
  //*
  for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=0 ; j<m_iGridSizeY-1 ; j++) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = (*m_vCurrentHx)[index] + (m_dt/(*m_vMuX)[index]) * ( ((*m_vCurrentEy)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEy)[index])/m_dz - ((*m_vCurrentEz)[index+m_iGridSizeX] - (*m_vCurrentEz)[index])/m_dy );
        (*m_vNewHy)[index] = (*m_vCurrentHy)[index] + (m_dt/(*m_vMuY)[index]) * ( ((*m_vCurrentEz)[index+1] - (*m_vCurrentEz)[index])/m_dx - ((*m_vCurrentEx)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEx)[index])/m_dz );
        (*m_vNewHz)[index] = (*m_vCurrentHz)[index] + (m_dt/(*m_vMuZ)[index]) * ( ((*m_vCurrentEx)[index+m_iGridSizeX] - (*m_vCurrentEx)[index])/m_dy - ((*m_vCurrentEy)[index+1] - (*m_vCurrentEy)[index])/m_dx );
      }
    }
  }


  updateH();
  //*/

  



  /*
  vector<double> vExactEx(m_iTotalDomainMemorySize);
  vector<double> vExactEy(m_iTotalDomainMemorySize);
  vector<double> vExactEz(m_iTotalDomainMemorySize);
  vector<double> vExactHx(m_iTotalDomainMemorySize);
  vector<double> vExactHy(m_iTotalDomainMemorySize);
  vector<double> vExactHz(m_iTotalDomainMemorySize);


  vector<double> vDiffEx(m_iTotalDomainMemorySize);
  vector<double> vDiffEy(m_iTotalDomainMemorySize);
  vector<double> vDiffEz(m_iTotalDomainMemorySize);
  vector<double> vDiffHx(m_iTotalDomainMemorySize);
  vector<double> vDiffHy(m_iTotalDomainMemorySize);
  vector<double> vDiffHz(m_iTotalDomainMemorySize);


  //double time = (timestep-1)*m_dt;
  //double time = timestep*m_dt;
  for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=0 ; j<m_iGridSizeY-1 ; j++) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;

        double exactEx = getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time);
        double exactEy = getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time);
        double exactEz = getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time);

        vDiffEx[index] = fabs( (*m_vCurrentEx)[index] - exactEx );
        vDiffEy[index] = fabs( (*m_vCurrentEy)[index] - exactEy );
        vDiffEz[index] = fabs( (*m_vCurrentEz)[index] - exactEz );

        double exactHx = getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt);
        double exactHy = getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt);
        double exactHz = getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time + 0.5*m_dt);

        vDiffHx[index] = fabs( (*m_vCurrentHx)[index] - exactHx );
        vDiffHy[index] = fabs( (*m_vCurrentHy)[index] - exactHy );
        vDiffHz[index] = fabs( (*m_vCurrentHz)[index] - exactHz );

      }
    }
  }
  //*/



  /*
  {
    ofstream fsEx("timestep_1_Ex.txt");
    ofstream fsEy("timestep_1_Ey.txt");
    ofstream fsEz("timestep_1_Ez.txt");

    for(int k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(int j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(int i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEx->at(index) << endl;
          fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEy->at(index) << endl;
          fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEz->at(index) << endl;
        }
      }
    }
    fsEx.close();
    fsEy.close();
    fsEz.close();
  }
  //*/



  return 0;
}


void SerialFieldSolver::updateE() {
  vector<double> *temp;

  temp = m_vCurrentEx;
  m_vCurrentEx = m_vNewEx;
  m_vNewEx = temp;
  
  temp = m_vCurrentEy;
  m_vCurrentEy = m_vNewEy;
  m_vNewEy = temp;
  
  temp = m_vCurrentEz;
  m_vCurrentEz = m_vNewEz;
  m_vNewEz = temp;
}




void SerialFieldSolver::updateH() {
  vector<double> *temp;

  temp = m_vCurrentHx;
  m_vCurrentHx = m_vNewHx;
  m_vNewHx = temp;
  
  temp = m_vCurrentHy;
  m_vCurrentHy = m_vNewHy;
  m_vNewHy = temp;
  
  temp = m_vCurrentHz;
  m_vCurrentHz = m_vNewHz;
  m_vNewHz = temp;
}



void SerialFieldSolver::loadBoundaryValues(double t) {

  ////////////////////////////////////////////////////////////////////////////////
  // Start : loading BC values of E
  ////////////////////////////////////////////////////////////////////////////////
  size_t index;
  size_t index2;
  //double t = timestep*m_dt;
  int i,j,k;

  // X : Ey
  #pragma omp parallel for private(j,k,index,index2)
  for(k=0 ; k<m_iGridSizeZ ; k++ ) {
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEy)[index] = getBoundaryXLowerEy(m_fLowerX, m_fLowerY + j*m_dy + 0.5*m_dx, m_fLowerZ + k*m_dz, t);

      index2 = (m_iGridSizeX - 1) + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEy)[index2] = getBoundaryXUpperEy(m_fUpperX, m_fLowerY + j*m_dy + 0.5*m_dx, m_fLowerZ + k*m_dz, t);
    }
  }
  // X  : Ez
  #pragma omp parallel for private(j,k,index,index2)
  for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEz)[index] = getBoundaryXLowerEz(m_fLowerX, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, t);

      index2 = (m_iGridSizeX - 1) +  m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEz)[index2] = getBoundaryXUpperEz(m_fUpperX, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, t);
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  // Y : Ex
  #pragma omp parallel for private(i,k,index,index2)
  for(k=0 ; k<m_iGridSizeZ ; k++ ) {
    for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
      index = i + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEx)[index] = getBoundaryYLowerEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY, m_fLowerZ + k*m_dz, t);

      index2 = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEx)[index2] = getBoundaryYUpperEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fUpperY, m_fLowerZ + k*m_dz, t);
    }
  }
  // Y : Ez
  #pragma omp parallel for private(i,k,index,index2)
  for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
    for(i=0 ; i<m_iGridSizeX ; i++ ) {
      index = i + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEz)[index] = getBoundaryYLowerEz(m_fLowerX + i*m_dx, m_fLowerY, m_fLowerZ + k*m_dz + 0.5*m_dz, t);

      index2 = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vCurrentEz)[index2] = getBoundaryYUpperEz(m_fLowerX + i*m_dx, m_fUpperY, m_fLowerZ + k*m_dz + 0.5*m_dz, t);
    }
  }

  ////////////////////////////////////////////////////////////////////////////

  // Z : Ex
  #pragma omp parallel for private(i,j,index,index2)
  for(j=0 ; j<m_iGridSizeY ; j++ ) {
    for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
      index = i + m_iGridSizeX*j;
      (*m_vCurrentEx)[index] = getBoundaryZLowerEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ, t);

      index2 = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
      (*m_vCurrentEx)[index2] = getBoundaryZUpperEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fUpperZ, t);
    }
  }
  // Z : Ey
  #pragma omp parallel for private(i,j,index,index2)
  for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
    for(i=0 ; i<m_iGridSizeX ; i++ ) {
      index = i + m_iGridSizeX*j;
      (*m_vCurrentEy)[index] = getBoundaryZLowerEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ, t);

      index2 = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
      (*m_vCurrentEy)[index2] = getBoundaryZUpperEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fUpperZ, t);
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  // End : loading BC values of E
  ////////////////////////////////////////////////////////////////////////////////

}


int SerialFieldSolver::findCellHavingPosition(double x, double y, double z, int& cellX, int& cellY, int& cellZ, double& dx, double& dy, double& dz) {
  if( (x < m_fLowerX) || (x > m_fUpperX) || (y < m_fLowerY) || (y > m_fUpperY) || (z < m_fLowerZ) || (z > m_fUpperZ) ) {
    return EM_ERR_COM_DOM_SCOPE;
  }

  double distanceX = x - m_fLowerX;
  double distanceY = y - m_fLowerY;
  double distanceZ = z - m_fLowerZ;

  cellX = static_cast<int>(floor(distanceX/m_dx));
  cellY = static_cast<int>(floor(distanceY/m_dy));
  cellZ = static_cast<int>(floor(distanceZ/m_dz));

  dx = distanceX/m_dx - cellX; // (distanceX - cellX*m_dx)/m_dx
  dy = distanceY/m_dy - cellY; // (distanceY - cellX*m_dy)/m_dy
  dz = distanceZ/m_dz - cellZ; // (distanceZ - cellX*m_dz)/m_dz

  if(cellX == m_iGridSizeX) {
    cellX--;
    dx = 1.0;
  }
  if(cellY == m_iGridSizeY) {
    cellY--;
    dy = 1.0;
  }
  if(cellZ == m_iGridSizeZ) {
    cellZ--;
    dz = 1.0;
  }

  return EM_SUCCESS;
}

int SerialFieldSolver::updateCurrent(double time, int cellX, int cellY, int cellZ, double jx_y0z0, double jx_y1z0, double jx_y0z1, double jx_y1z1, double jy_x0z0, double jy_x1z0, double jy_x0z1, double jy_x1z1, double jz_x0y0, double jz_x1y0, double jz_x0y1, double jz_x1y1) {
  int index = cellX + m_iGridSizeX*cellY + m_iGridSizeX*m_iGridSizeY*cellZ;

  m_vCurrentEx->at(index) -= m_dt*jx_y0z0/getEpsilon(cellX*m_dx + 0.5*m_dx, cellY*m_dy, cellZ*m_dz, time);
  m_vCurrentEx->at(index + m_iGridSizeX) -= m_dt*jx_y1z0/getEpsilon(cellX*m_dx + 0.5*m_dx, (cellY+1)*m_dy, cellZ*m_dz, time);
  m_vCurrentEx->at(index + m_iGridSizeX*m_iGridSizeY) -= m_dt*jx_y0z1/getEpsilon(cellX*m_dx + 0.5*m_dx, cellY*m_dy, (cellZ+1)*m_dz, time);
  m_vCurrentEx->at(index + m_iGridSizeX + m_iGridSizeX*m_iGridSizeY) -= m_dt*jx_y1z1/getEpsilon(cellX*m_dx + 0.5*m_dx, (cellY+1)*m_dy, (cellZ+1)*m_dz, time);

  m_vCurrentEy->at(index) -= m_dt*jy_x0z0/getEpsilon(cellX*m_dx, cellY*m_dy + 0.5*m_dy, cellZ*m_dz, time);
  m_vCurrentEy->at(index + 1) -= m_dt*jy_x1z0/getEpsilon((cellX+1)*m_dx, cellY*m_dy + 0.5*m_dy, cellZ*m_dz, time);
  m_vCurrentEy->at(index + m_iGridSizeX*m_iGridSizeY) -= m_dt*jy_x0z1/getEpsilon(cellX*m_dx, cellY*m_dy + 0.5*m_dy, (cellZ+1)*m_dz, time);
  m_vCurrentEy->at(index + 1 + m_iGridSizeX*m_iGridSizeY) -= m_dt*jy_x1z1/getEpsilon((cellX+1)*m_dx, cellY*m_dy + 0.5*m_dy, (cellZ+1)*m_dz, time);

  m_vCurrentEz->at(index) -= m_dt*jz_x0y0/getEpsilon(cellX*m_dx, cellY*m_dy, cellZ*m_dz + 0.5*m_dz, time);
  m_vCurrentEz->at(index + 1) -= m_dt*jz_x1y0/getEpsilon((cellX+1)*m_dx, cellY*m_dy, cellZ*m_dz + 0.5*m_dz, time);
  m_vCurrentEz->at(index + m_iGridSizeX) -= m_dt*jz_x0y1/getEpsilon(cellX*m_dx, (cellY+1)*m_dy, cellZ*m_dz + 0.5*m_dz, time);
  m_vCurrentEz->at(index + 1 + m_iGridSizeX) -= m_dt*jz_x1y1/getEpsilon((cellX+1)*m_dx, (cellY+1)*m_dy, cellZ*m_dz + 0.5*m_dz, time);


  /*
  cout << "index : " << index << endl;
  cout << "Ez at <index> : " << m_vCurrentEz->at(index) << endl;
  cout << "Ez at <index + 1> : " << m_vCurrentEz->at(index + 1) << endl;
  cout << "Ez at <index + m_iGridSizeX> : " << m_vCurrentEz->at(index + m_iGridSizeX) << endl;
  cout << "Ez at <index + 1 + m_iGridSizeX> : " << m_vCurrentEz->at(index + 1 + m_iGridSizeX) <<"\n"<<endl;
  */
  /*
  {
    size_t i,j,k,index;
    stringstream fnEx;
    stringstream fnEy;
    stringstream fnEz;
    fnEx << "after_updateCurrent_Ex.txt";
    fnEy << "after_updateCurrent_Ey.txt";
    fnEz << "after_updateCurrent_Ez.txt";
    FILE* fsEx = fopen(fnEx.str().c_str(), "w");
    FILE* fsEy = fopen(fnEy.str().c_str(), "w");
    FILE* fsEz = fopen(fnEz.str().c_str(), "w");
    
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fprintf(fsEx, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, m_vCurrentEx->at(index));
          fprintf(fsEy, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, m_vCurrentEy->at(index));
          fprintf(fsEz, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, m_vCurrentEz->at(index));
        }
      }
    }
    fclose(fsEx);
    fclose(fsEy);
    fclose(fsEz);
    
    //_exit(0);
  }
  //*/

  /*
  {
    size_t i,j,k,index;
    stringstream fnEx;
    stringstream fnEy;
    stringstream fnEz;
    fnEx << "after_updateCurrent_Ex.txt";
    fnEy << "after_updateCurrent_Ey.txt";
    fnEz << "after_updateCurrent_Ez.txt";
    ofstream fsEx(fnEx.str().c_str());
    ofstream fsEy(fnEy.str().c_str());
    ofstream fsEz(fnEz.str().c_str());
    fsEx.precision(60);
    fsEy.precision(60);
    fsEz.precision(60);
    
    for(k=0 ; k<m_iGridSizeZ ; k++) {
      for(j=0 ; j<m_iGridSizeY ; j++) {
        for(i=0 ; i<m_iGridSizeX ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEx)[index] << endl;
          fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEy)[index] << endl;
          fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEz)[index] << endl;
        }
      }
    }
    fsEx.close();
    fsEx.close();
    fsEx.close();
    
    _exit(0);
  }
  */


  return EM_SUCCESS;
}


double SerialFieldSolver::getMu(double x, double y, double z, double time) {
  return 0.0;
}

double SerialFieldSolver::getEpsilon(double x, double y, double z, double time) {
  return 0.0;
}

double SerialFieldSolver::getSigma(double x, double y, double z, double time) {
  return 0.0;
}

int SerialFieldSolver::getElectricIntensityAt(double time, size_t timestep, double x, double y, double z, double& ex, double& ey, double& ez) {

  int cellX;
  int cellY;
  int cellZ;
  int i;
  double dx;
  double dy;
  double dz;

  int error = findCellHavingPosition(x,y,z,cellX,cellY,cellZ, dx, dy, dz);
  if(error) {
    switch(error) {
      case EM_ERR_COM_DOM_SCOPE:
        EM_ERROR(STR_ERR_COM_DOM_SCOPE);
        break;

      default:
        EM_ERROR(STR_ERR_UNKNOWN);
        break;
    }
    return error;
  }

  double lxlylzEx;
  double uxlylzEx;
  double lxuylzEx;
  double uxuylzEx;
  double lxlyuzEx;
  double uxlyuzEx;
  double lxuyuzEx;
  double uxuyuzEx;

  double lxlylzEy;
  double uxlylzEy;
  double lxuylzEy;
  double uxuylzEy;
  double lxlyuzEy;
  double uxlyuzEy;
  double lxuyuzEy;
  double uxuyuzEy;

  double lxlylzEz;
  double uxlylzEz;
  double lxuylzEz;
  double uxuylzEz;
  double lxlyuzEz;
  double uxlyuzEz;
  double lxuyuzEz;
  double uxuyuzEz;


  i = cellX + cellY*m_iGridSizeX + cellZ*m_iGridSizeX*m_iGridSizeY;
  int x1 = i + 1;
  int y1 = i + m_iGridSizeX;
  int z1 = i + m_iGridSizeX*m_iGridSizeY;
  int x1y1 = x1 + m_iGridSizeX;
  int x1z1 = z1 + 1;
  int y1z1 = z1 + m_iGridSizeX;

  if(cellX == 0) {
    lxlylzEx = getBoundaryXLowerEx(m_fLowerX, m_fLowerY + m_dy*cellY, m_fLowerZ + m_dz*cellZ, time);
    lxuylzEx = getBoundaryXLowerEx(m_fLowerX, m_fLowerY + m_dy*(cellY+1), m_fLowerZ + m_dz*cellZ, time);
    lxlyuzEx = getBoundaryXLowerEx(m_fLowerX, m_fLowerY + m_dy*cellY, m_fLowerZ + m_dz*(cellZ+1), time);
    lxuyuzEx = getBoundaryXLowerEx(m_fLowerX, m_fLowerY + m_dy*(cellY+1), m_fLowerZ + m_dz*(cellZ+1), time);
  } else {
    lxlylzEx = 0.5*(m_vCurrentEx->at(i) + m_vCurrentEx->at(i-1));
    lxuylzEx = 0.5*(m_vCurrentEx->at(y1) + m_vCurrentEx->at(y1-1));
    lxlyuzEx = 0.5*(m_vCurrentEx->at(z1) + m_vCurrentEx->at(z1-1));
    lxuyuzEx = 0.5*(m_vCurrentEx->at(y1z1) + m_vCurrentEx->at(y1z1-1));
  }

  if(cellX == m_iGridSizeX - 1) {
    uxlylzEx = getBoundaryXUpperEx(m_fUpperX, m_fLowerY + m_dy*cellY, m_fLowerZ + m_dz*cellZ, time);
    uxuylzEx = getBoundaryXUpperEx(m_fUpperX, m_fLowerY + m_dy*(cellY+1), m_fLowerZ + m_dz*cellZ, time);
    uxlyuzEx = getBoundaryXUpperEx(m_fUpperX, m_fLowerY + m_dy*cellY, m_fLowerZ + m_dz*(cellZ+1), time);
    uxuyuzEx = getBoundaryXUpperEx(m_fUpperX, m_fLowerY + m_dy*(cellY+1), m_fLowerZ + m_dz*(cellZ+1), time);
  } else {
    uxlylzEx = 0.5*(m_vCurrentEx->at(i) + m_vCurrentEx->at(i+1));
    uxuylzEx = 0.5*(m_vCurrentEx->at(y1) + m_vCurrentEx->at(y1+1));
    uxlyuzEx = 0.5*(m_vCurrentEx->at(z1) + m_vCurrentEx->at(z1+1));
    uxuyuzEx = 0.5*(m_vCurrentEx->at(y1z1) + m_vCurrentEx->at(y1z1+1));
  }
  


  if(cellY == 0) {
    lxlylzEy = getBoundaryYLowerEy(m_fLowerX + m_dx*cellX, m_fLowerY, m_fLowerZ + m_dz*cellZ, time);
    uxlylzEy = getBoundaryYLowerEy(m_fLowerX + m_dx*(cellX+1), m_fLowerY, m_fLowerZ + m_dz*cellZ, time);
    lxlyuzEy = getBoundaryYLowerEy(m_fLowerX + m_dx*cellX, m_fLowerY, m_fLowerZ + m_dz*(cellZ+1), time);
    uxlyuzEy = getBoundaryYLowerEy(m_fLowerX + m_dx*(cellX+1), m_fLowerY, m_fLowerZ + m_dz*(cellZ+1), time);
  } else {
    lxlylzEy = 0.5*(m_vCurrentEy->at(i) + m_vCurrentEy->at(i-m_iGridSizeX));
    uxlylzEy = 0.5*(m_vCurrentEy->at(x1) + m_vCurrentEy->at(x1-m_iGridSizeX));
    lxlyuzEy = 0.5*(m_vCurrentEy->at(z1) + m_vCurrentEy->at(z1-m_iGridSizeX));
    uxlyuzEy = 0.5*(m_vCurrentEy->at(x1z1) + m_vCurrentEy->at(x1z1-m_iGridSizeX));
  }

  if(cellY == m_iGridSizeY - 1) {
    lxuylzEy = getBoundaryYUpperEy(m_fLowerX + m_dx*cellX, m_fUpperY, m_fLowerZ + m_dz*cellZ, time);
    uxuylzEy = getBoundaryYUpperEy(m_fLowerX + m_dx*(cellX+1), m_fUpperY, m_fLowerZ + m_dz*cellZ, time);
    lxuyuzEy = getBoundaryYUpperEy(m_fLowerX + m_dx*cellX, m_fUpperY, m_fLowerZ + m_dz*(cellZ+1), time);
    uxuyuzEy = getBoundaryYUpperEy(m_fLowerX + m_dx*(cellX+1), m_fUpperY, m_fLowerZ + m_dz*(cellZ+1), time);
  } else {
    lxuylzEy = 0.5*(m_vCurrentEy->at(i) + m_vCurrentEy->at(i+m_iGridSizeX));
    uxuylzEy = 0.5*(m_vCurrentEy->at(x1) + m_vCurrentEy->at(x1+m_iGridSizeX));
    lxuyuzEy = 0.5*(m_vCurrentEy->at(z1) + m_vCurrentEy->at(z1+m_iGridSizeX));
    uxuyuzEy = 0.5*(m_vCurrentEy->at(x1z1) + m_vCurrentEy->at(x1z1+m_iGridSizeX));
  }



  if(cellZ == 0) {
    lxlylzEz = getBoundaryZLowerEz(m_fLowerX + m_dx*cellX, m_fLowerY + m_dy*cellY, m_fLowerZ, time);
    uxlylzEz = getBoundaryZLowerEz(m_fLowerX + m_dx*(cellX+1), m_fLowerY + m_dy*cellY, m_fLowerZ, time);
    lxuylzEz = getBoundaryZLowerEz(m_fLowerX + m_dx*cellX, m_fLowerY + m_dy*(cellY+1), m_fLowerZ, time);
    uxuylzEz = getBoundaryZLowerEz(m_fLowerX + m_dx*(cellX+1), m_fLowerY + m_dy*(cellY+1), m_fLowerZ, time);
  } else {
    lxlylzEz = 0.5*(m_vCurrentEz->at(i) + m_vCurrentEz->at(i-m_iGridSizeX*m_iGridSizeY));
    uxlylzEz = 0.5*(m_vCurrentEz->at(x1) + m_vCurrentEz->at(x1-m_iGridSizeX*m_iGridSizeY));
    lxuylzEz = 0.5*(m_vCurrentEz->at(y1) + m_vCurrentEz->at(y1-m_iGridSizeX*m_iGridSizeY));
    uxuylzEz = 0.5*(m_vCurrentEz->at(x1y1) + m_vCurrentEz->at(x1y1-m_iGridSizeX*m_iGridSizeY));
  }

  if(cellZ == m_iGridSizeZ - 1) {
    lxlyuzEz = getBoundaryZUpperEz(m_fLowerX + m_dx*cellX, m_fLowerY + m_dy*cellY, m_fUpperZ, time);
    uxlyuzEz = getBoundaryZUpperEz(m_fLowerX + m_dx*(cellX+1), m_fLowerY + m_dy*cellY, m_fUpperZ, time);
    lxuyuzEz = getBoundaryZUpperEz(m_fLowerX + m_dx*cellX, m_fLowerY + m_dy*(cellY+1), m_fUpperZ, time);
    uxuyuzEz = getBoundaryZUpperEz(m_fLowerX + m_dx*(cellX+1), m_fLowerY + m_dy*(cellY+1), m_fUpperZ, time);
  } else {
    lxlyuzEz = 0.5*(m_vCurrentEz->at(i) + m_vCurrentEz->at(i+m_iGridSizeX*m_iGridSizeY));
    uxlyuzEz = 0.5*(m_vCurrentEz->at(x1) + m_vCurrentEz->at(x1+m_iGridSizeX*m_iGridSizeY));
    lxuyuzEz = 0.5*(m_vCurrentEz->at(y1) + m_vCurrentEz->at(y1+m_iGridSizeX*m_iGridSizeY));
    uxuyuzEz = 0.5*(m_vCurrentEz->at(x1y1) + m_vCurrentEz->at(x1y1+m_iGridSizeX*m_iGridSizeY));
  }
  /*
  cout << "lxlylzEz : " << lxlyuzEz << endl;
  cout << "uxlylzEz : " << uxlylzEz << endl;
  cout << "lxuylzEz : " << lxuylzEz << endl;
  cout << "uxuylzEz : " << uxuylzEz << endl;
  cout << "lxlyuzEz : " << lxlyuzEz << endl;
  cout << "uxlyuzEz : " << uxlyuzEz << endl;
  cout << "lxuyuzEz : " << lxuyuzEz << endl;
  cout << "uxuyuzEz : " << uxuyuzEz << endl;
  */
  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////


  ex = lxlylzEx*(1-dx)*(1-dy)*(1-dz) + uxlylzEx*dx*(1-dy)*(1-dz) + lxuylzEx*(1-dx)*dy*(1-dz) + uxuylzEx*dx*dy*(1-dz) + lxlyuzEx*(1-dx)*(1-dy)*dz + uxlyuzEx*dx*(1-dy)*dz + lxuyuzEx*(1-dx)*dy*dz + uxuyuzEx*dx*dy*dz;
  ey = lxlylzEy*(1-dx)*(1-dy)*(1-dz) + uxlylzEy*dx*(1-dy)*(1-dz) + lxuylzEy*(1-dx)*dy*(1-dz) + uxuylzEy*dx*dy*(1-dz) + lxlyuzEy*(1-dx)*(1-dy)*dz + uxlyuzEy*dx*(1-dy)*dz + lxuyuzEy*(1-dx)*dy*dz + uxuyuzEy*dx*dy*dz;
  ez = lxlylzEz*(1-dx)*(1-dy)*(1-dz) + uxlylzEz*dx*(1-dy)*(1-dz) + lxuylzEz*(1-dx)*dy*(1-dz) + uxuylzEz*dx*dy*(1-dz) + lxlyuzEz*(1-dx)*(1-dy)*dz + uxlyuzEz*dx*(1-dy)*dz + lxuyuzEz*(1-dx)*dy*dz + uxuyuzEz*dx*dy*dz;
  //ex = lxlylzEx*dx*dy*dz + uxlylzEx*(1-dx)*dy*dz + lxuylzEx*dx*(1-dy)*dz + uxuylzEx*(1-dx)*(1-dy)*dz + lxlyuzEx*dx*dy*(1-dz) + uxlyuzEx*(1-dx)*dy*(1-dz) + lxuyuzEx*dx*(1-dy)*(1-dz) + uxuyuzEx*(1-dx)*(1-dy)*(1-dz);
  //ey = lxlylzEy*dx*dy*dz + uxlylzEy*(1-dx)*dy*dz + lxuylzEy*dx*(1-dy)*dz + uxuylzEy*(1-dx)*(1-dy)*dz + lxlyuzEy*dx*dy*(1-dz) + uxlyuzEy*(1-dx)*dy*(1-dz) + lxuyuzEy*dx*(1-dy)*(1-dz) + uxuyuzEy*(1-dx)*(1-dy)*(1-dz);
  //ez = lxlylzEz*dx*dy*dz + uxlylzEz*(1-dx)*dy*dz + lxuylzEz*dx*(1-dy)*dz + uxuylzEz*(1-dx)*(1-dy)*dz + lxlyuzEz*dx*dy*(1-dz) + uxlyuzEz*(1-dx)*dy*(1-dz) + lxuyuzEz*dx*(1-dy)*(1-dz) + uxuyuzEz*(1-dx)*(1-dy)*(1-dz);


  return EM_SUCCESS;
}

int SerialFieldSolver::getMagneticDensityAt(double time, size_t timestep, double x, double y, double z, double& bx, double& by, double& bz) {
  
  int cellX;
  int cellY;
  int cellZ;
  int i;
  double dx;
  double dy;
  double dz;

  int positionIndex = findCellHavingPosition(x,y,z,cellX,cellY,cellZ, dx, dy, dz);
  if(positionIndex == -1) return -1;


  double lxlylzBx;
  double uxlylzBx;
  double lxuylzBx;
  double uxuylzBx;
  double lxlyuzBx;
  double uxlyuzBx;
  double lxuyuzBx;
  double uxuyuzBx;

  double lxlylzBy;
  double uxlylzBy;
  double lxuylzBy;
  double uxuylzBy;
  double lxlyuzBy;
  double uxlyuzBy;
  double lxuyuzBy;
  double uxuyuzBy;

  double lxlylzBz;
  double uxlylzBz;
  double lxuylzBz;
  double uxuylzBz;
  double lxlyuzBz;
  double uxlyuzBz;
  double lxuyuzBz;
  double uxuyuzBz;

  i = cellX + cellY*m_iGridSizeX + cellZ*m_iGridSizeX*m_iGridSizeY;

  int lxlylz;
  int uxlylz;
  int lxuylz;
  int uxuylz;
  int lxlyuz;
  int uxlyuz;
  int lxuyuz;
  int uxuyuz;

  /*
  int tmp_lxlylz;
  int tmp_lxuylz;
  int tmp_lxlyuz;
  int tmp_lxuyuz;
  int tmp_uxlylz;
  int tmp_uxuylz;
  int tmp_uxlyuz;
  int tmp_uxuyuz;
  //*/

  double rdx, rdy, rdz; // relative dx, dy, dz to adjusted new vertices.
  double mucoeff = 0.5*getMu(x,y,z,time);

  if( (x >= m_fLowerX + 0.5*m_dx) && (x <= m_fUpperX - 0.5*m_dx) && (y >= m_fLowerY + 0.5*m_dy) && (y <= m_fUpperY - 0.5*m_dy) && (z >= m_fLowerZ + 0.5*m_dz) && (z <= m_fUpperZ - 0.5*m_dz) ) {
    // Hx
    if( (dy < 0.5) && (dz < 0.5) ) {
      lxlylz = i - m_iGridSizeX - m_iGridSizeX*m_iGridSizeY;
      lxuylz = i - m_iGridSizeX*m_iGridSizeY;
      lxlyuz = i - m_iGridSizeX;
      lxuyuz = i;
      uxlylz = i - m_iGridSizeX - m_iGridSizeX*m_iGridSizeY + 1;
      uxuylz = i - m_iGridSizeX*m_iGridSizeY + 1;
      uxlyuz = i - m_iGridSizeX + 1;
      uxuyuz = i + 1;
      rdx = dx;
      rdy = dy + 0.5;
      rdz = dz + 0.5;
    } else if ( (dy < 0.5) && (dz >= 0.5) ) {
      lxlylz = i - m_iGridSizeX;
      lxuylz = i;
      lxlyuz = i - m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      lxuyuz = i + m_iGridSizeX*m_iGridSizeY;
      uxlylz = i - m_iGridSizeX + 1;
      uxuylz = i + 1;
      uxlyuz = i - m_iGridSizeX + m_iGridSizeX*m_iGridSizeY + 1;
      uxuyuz = i + m_iGridSizeX*m_iGridSizeY + 1;
      rdx = dx;
      rdy = dy + 0.5;
      rdz = dz - 0.5;
    } else if ( (dy >= 0.5) && (dz < 0.5) ) {
      lxlylz = i - m_iGridSizeX*m_iGridSizeY;
      lxuylz = i + m_iGridSizeX - m_iGridSizeX*m_iGridSizeY;
      lxlyuz = i;
      lxuyuz = i + m_iGridSizeX;
      uxlylz = i - m_iGridSizeX*m_iGridSizeY + 1;
      uxuylz = i + m_iGridSizeX - m_iGridSizeX*m_iGridSizeY + 1;
      uxlyuz = i + 1;
      uxuyuz = i + m_iGridSizeX + 1;
      rdx = dx;
      rdy = dy - 0.5;
      rdz = dz + 0.5;
    } else {
      lxlylz = i;
      lxuylz = i + m_iGridSizeX;
      lxlyuz = i + m_iGridSizeX*m_iGridSizeY;
      lxuyuz = i + m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      uxlylz = i + 1;
      uxuylz = i + m_iGridSizeX + 1;
      uxlyuz = i + m_iGridSizeX*m_iGridSizeY + 1;
      uxuyuz = i + m_iGridSizeX + m_iGridSizeX*m_iGridSizeY + 1;
      rdx = dx;
      rdy = dy - 0.5;
      rdz = dz - 0.5;
    }
    
    lxlylzBx = mucoeff*(m_vCurrentHx->at(lxlylz) + m_vNewHx->at(lxlylz));
    //cout << "lxlylz=" << lxlylz << ", lxlylzBx=" << lxlylzBx << endl;
    uxlylzBx = mucoeff*(m_vCurrentHx->at(uxlylz) + m_vNewHx->at(uxlylz));
    //cout << "uxlylz=" << uxlylz << ", uxlylzBx=" << uxlylzBx << endl;
    lxuylzBx = mucoeff*(m_vCurrentHx->at(lxuylz) + m_vNewHx->at(lxuylz));
    //cout << "lxuylz=" << lxuylz << ", lxuylzBx=" << lxuylzBx << endl;
    uxuylzBx = mucoeff*(m_vCurrentHx->at(uxuylz) + m_vNewHx->at(uxuylz));
    //cout << "uxuylz=" << uxuylz << ", uxuylzBx=" << uxuylzBx << endl;
    lxlyuzBx = mucoeff*(m_vCurrentHx->at(lxlyuz) + m_vNewHx->at(lxlyuz));
    //cout << "lxlyuz=" << lxlyuz << ", lxlyuzBx=" << lxlyuzBx << endl;
    uxlyuzBx = mucoeff*(m_vCurrentHx->at(uxlyuz) + m_vNewHx->at(uxlyuz));
    //cout << "uxlyuz=" << uxlyuz << ", uxlyuzBx=" << uxlyuzBx << endl;
    lxuyuzBx = mucoeff*(m_vCurrentHx->at(lxuyuz) + m_vNewHx->at(lxuyuz));
    //cout << "lxuyuz=" << lxuyuz << ", lxuyuzBx=" << lxuyuzBx << endl;
    uxuyuzBx = mucoeff*(m_vCurrentHx->at(uxuyuz) + m_vNewHx->at(uxuyuz));
    //cout << "uxuyuz=" << uxuyuz << ", uxuyuzBx=" << uxuyuzBx << endl;
    bx = lxlylzBx*(1-rdx)*(1-rdy)*(1-rdz) + uxlylzBx*rdx*(1-rdy)*(1-rdz) + lxuylzBx*(1-rdx)*rdy*(1-rdz) + uxuylzBx*rdx*rdy*(1-rdz) + lxlyuzBx*(1-rdx)*(1-rdy)*rdz + uxlyuzBx*rdx*(1-rdy)*rdz + lxuyuzBx*(1-rdx)*rdy*rdz + uxuyuzBx*rdx*rdy*rdz;

    /*
    tmp_lxlylz = lxlylz;
    tmp_lxuylz = lxuylz;
    tmp_lxlyuz = lxlyuz;
    tmp_lxuyuz = lxuyuz;
    tmp_uxlylz = uxlylz;
    tmp_uxuylz = uxuylz;
    tmp_uxlyuz = uxlyuz;
    tmp_uxuyuz = uxuyuz;
    //*/


    // Hy
    if( (dx < 0.5) && (dz < 0.5) ) {
      lxlylz = i - 1 - m_iGridSizeX*m_iGridSizeY;
      uxlylz = i - m_iGridSizeX*m_iGridSizeY;
      lxlyuz = i - 1;
      uxlyuz = i;
      lxuylz = i - 1 - m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      uxuylz = i - m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      lxuyuz = i - 1 + m_iGridSizeX;
      uxuyuz = i + m_iGridSizeX;
      rdx = dx + 0.5;
      rdy = dy;
      rdz = dz + 0.5;
    } else if ( (dx < 0.5) && (dz >= 0.5) ) {
      lxlylz = i - 1;
      uxlylz = i;
      lxlyuz = i - 1 + m_iGridSizeX*m_iGridSizeY;
      uxlyuz = i + m_iGridSizeX*m_iGridSizeY;
      lxuylz = i - 1 + m_iGridSizeX;
      uxuylz = i + m_iGridSizeX;
      lxuyuz = i - 1 + m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      uxuyuz = i + m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      rdx = dx + 0.5;
      rdy = dy;
      rdz = dz - 0.5;
    } else if ( (dx >= 0.5) && (dz < 0.5) ) {
      lxlylz = i - m_iGridSizeX*m_iGridSizeY;
      uxlylz = i + 1 - m_iGridSizeX*m_iGridSizeY;
      lxlyuz = i;
      uxlyuz = i + 1;
      lxuylz = i - m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      uxuylz = i + 1 - m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      lxuyuz = i + m_iGridSizeX;
      uxuyuz = i + 1 + m_iGridSizeX;
      rdx = dx - 0.5;
      rdy = dy;
      rdz = dz + 0.5;
    } else {
      lxlylz = i;
      uxlylz = i + 1;
      lxlyuz = i + m_iGridSizeX*m_iGridSizeY;
      uxlyuz = i + 1 + m_iGridSizeX*m_iGridSizeY;
      lxuylz = i + m_iGridSizeX;
      uxuylz = i + 1 + m_iGridSizeX;
      lxuyuz = i + m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      uxuyuz = i + 1 + m_iGridSizeX*m_iGridSizeY + m_iGridSizeX;
      rdx = dx - 0.5;
      rdy = dy;
      rdz = dz - 0.5;
    }
    
    lxlylzBy = mucoeff*(m_vCurrentHy->at(lxlylz) + m_vNewHy->at(lxlylz));
    uxlylzBy = mucoeff*(m_vCurrentHy->at(uxlylz) + m_vNewHy->at(uxlylz));
    lxuylzBy = mucoeff*(m_vCurrentHy->at(lxuylz) + m_vNewHy->at(lxuylz));
    uxuylzBy = mucoeff*(m_vCurrentHy->at(uxuylz) + m_vNewHy->at(uxuylz));
    lxlyuzBy = mucoeff*(m_vCurrentHy->at(lxlyuz) + m_vNewHy->at(lxlyuz));
    uxlyuzBy = mucoeff*(m_vCurrentHy->at(uxlyuz) + m_vNewHy->at(uxlyuz));
    lxuyuzBy = mucoeff*(m_vCurrentHy->at(lxuyuz) + m_vNewHy->at(lxuyuz));
    uxuyuzBy = mucoeff*(m_vCurrentHy->at(uxuyuz) + m_vNewHy->at(uxuyuz));
    by = lxlylzBy*(1-rdx)*(1-rdy)*(1-rdz) + uxlylzBy*rdx*(1-rdy)*(1-rdz) + lxuylzBy*(1-rdx)*rdy*(1-rdz) + uxuylzBy*rdx*rdy*(1-rdz) + lxlyuzBy*(1-rdx)*(1-rdy)*rdz + uxlyuzBy*rdx*(1-rdy)*rdz + lxuyuzBy*(1-rdx)*rdy*rdz + uxuyuzBy*rdx*rdy*rdz;

    // Hz
    if( (dx < 0.5) && (dy < 0.5) ) {
      lxlylz = i - m_iGridSizeX - 1;
      uxlylz = i - m_iGridSizeX;
      lxuylz = i - 1;
      uxuylz = i;
      lxlyuz = i - m_iGridSizeX - 1 + m_iGridSizeX*m_iGridSizeY;
      uxlyuz = i - m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      lxuyuz = i - 1 + m_iGridSizeX*m_iGridSizeY;
      uxuyuz = i + m_iGridSizeX*m_iGridSizeY;
      rdx = dx + 0.5;
      rdy = dy + 0.5;
      rdz = dz;
    } else if ( (dx < 0.5) && (dy >= 0.5) ) {
      lxlylz = i - 1;
      uxlylz = i;
      lxuylz = i - 1 + m_iGridSizeX;
      uxuylz = i + m_iGridSizeX;
      lxlyuz = i - 1 + m_iGridSizeX*m_iGridSizeY;
      uxlyuz = i + m_iGridSizeX*m_iGridSizeY;
      lxuyuz = i - 1 + m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      uxuyuz = i + m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      rdx = dx + 0.5;
      rdy = dy - 0.5;
      rdz = dz;
    } else if ( (dx >= 0.5) && (dy < 0.5) ) {
      lxlylz = i - m_iGridSizeX;
      uxlylz = i - m_iGridSizeX + 1;
      lxuylz = i;
      uxuylz = i + 1;
      lxlyuz = i - m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      uxlyuz = i - m_iGridSizeX + 1 + m_iGridSizeX*m_iGridSizeY;
      lxuyuz = i + m_iGridSizeX*m_iGridSizeY;
      uxuyuz = i + 1 + m_iGridSizeX*m_iGridSizeY;
      rdx = dx - 0.5;
      rdy = dy + 0.5;
      rdz = dz;
    } else {
      lxlylz = i;
      uxlylz = i + 1;
      lxuylz = i + m_iGridSizeX;
      uxuylz = i + m_iGridSizeX + 1;
      lxlyuz = i + m_iGridSizeX*m_iGridSizeY;
      uxlyuz = i + 1 + m_iGridSizeX*m_iGridSizeY;
      lxuyuz = i + m_iGridSizeX + m_iGridSizeX*m_iGridSizeY;
      uxuyuz = i + m_iGridSizeX + 1 + m_iGridSizeX*m_iGridSizeY;
      rdx = dx - 0.5;
      rdy = dy - 0.5;
      rdz = dz;
    }
    
    lxlylzBz = mucoeff*(m_vCurrentHz->at(lxlylz) + m_vNewHz->at(lxlylz));
    uxlylzBz = mucoeff*(m_vCurrentHz->at(uxlylz) + m_vNewHz->at(uxlylz));
    lxuylzBz = mucoeff*(m_vCurrentHz->at(lxuylz) + m_vNewHz->at(lxuylz));
    uxuylzBz = mucoeff*(m_vCurrentHz->at(uxuylz) + m_vNewHz->at(uxuylz));
    lxlyuzBz = mucoeff*(m_vCurrentHz->at(lxlyuz) + m_vNewHz->at(lxlyuz));
    uxlyuzBz = mucoeff*(m_vCurrentHz->at(uxlyuz) + m_vNewHz->at(uxlyuz));
    lxuyuzBz = mucoeff*(m_vCurrentHz->at(lxuyuz) + m_vNewHz->at(lxuyuz));
    uxuyuzBz = mucoeff*(m_vCurrentHz->at(uxuyuz) + m_vNewHz->at(uxuyuz));
    bz = lxlylzBz*(1-rdx)*(1-rdy)*(1-rdz) + uxlylzBz*rdx*(1-rdy)*(1-rdz) + lxuylzBz*(1-rdx)*rdy*(1-rdz) + uxuylzBz*rdx*rdy*(1-rdz) + lxlyuzBz*(1-rdx)*(1-rdy)*rdz + uxlyuzBz*rdx*(1-rdy)*rdz + lxuyuzBz*(1-rdx)*rdy*rdz + uxuyuzBz*rdx*rdy*rdz;
  } else { // in the case that poisiton is located at boundary side.
    // TODO : implementation.
    return -1;
  }







  

  /*
  if(timestep == 2) {
    int i,j,k,index;
    double mu = getMu(0,0,0,0);
    stringstream fnBx;
    stringstream fnBy;
    stringstream fnBz;
    fnBx << "serial_Bx.txt";
    fnBy << "serial_By.txt";
    fnBz << "serial_Bz.txt";
    ofstream fsBx(fnBx.str().c_str());
    ofstream fsBy(fnBy.str().c_str());
    ofstream fsBz(fnBz.str().c_str());
    fsBx.precision(60);
    fsBy.precision(60);
    fsBz.precision(60);


    lxlylzBx = mucoeff*(m_vCurrentHx->at(tmp_lxlylz) + m_vNewHx->at(tmp_lxlylz));
    uxlylzBx = mucoeff*(m_vCurrentHx->at(tmp_uxlylz) + m_vNewHx->at(tmp_uxlylz));
    lxuylzBx = mucoeff*(m_vCurrentHx->at(tmp_lxuylz) + m_vNewHx->at(tmp_lxuylz));
    uxuylzBx = mucoeff*(m_vCurrentHx->at(tmp_uxuylz) + m_vNewHx->at(tmp_uxuylz));
    lxlyuzBx = mucoeff*(m_vCurrentHx->at(tmp_lxlyuz) + m_vNewHx->at(tmp_lxlyuz));
    uxlyuzBx = mucoeff*(m_vCurrentHx->at(tmp_uxlyuz) + m_vNewHx->at(tmp_uxlyuz));
    lxuyuzBx = mucoeff*(m_vCurrentHx->at(tmp_lxuyuz) + m_vNewHx->at(tmp_lxuyuz));
    uxuyuzBx = mucoeff*(m_vCurrentHx->at(tmp_uxuyuz) + m_vNewHx->at(tmp_uxuyuz));


    fsBx << "timestep : " << timestep << ", cellX=" << cellX << ", cellY=" << cellY << ", cellZ=" << cellZ << endl;
    fsBx << "timestep : " << timestep << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << endl;
    fsBx << "timestep : " << timestep << ", rdx=" << rdx << ", rdy=" << rdy << ", rdz=" << rdz << endl;
    fsBx << "timestep : " << timestep << ", dx=" << dx << endl;
    fsBx << "timestep : " << timestep << ", dy=" << dy << endl;
    fsBx << "timestep : " << timestep << ", dz=" << dz << endl;
    fsBx << "timestep : " << timestep << ", bx=" << bx << endl;
    fsBx << "timestep : " << timestep << ", lxlylzBx=" << lxlylzBx << endl;
    fsBx << "timestep : " << timestep << ", uxlylzBx=" << uxlylzBx << endl;
    fsBx << "timestep : " << timestep << ", lxuylzBx=" << lxuylzBx << endl;
    fsBx << "timestep : " << timestep << ", uxuylzBx=" << uxuylzBx << endl;
    fsBx << "timestep : " << timestep << ", lxlyuzBx=" << lxlyuzBx << endl;
    fsBx << "timestep : " << timestep << ", uxlyuzBx=" << uxlyuzBx << endl;
    fsBx << "timestep : " << timestep << ", lxuyuzBx=" << lxuyuzBx << endl;
    fsBx << "timestep : " << timestep << ", uxuyuzBx=" << uxuyuzBx << endl;
    
    fsBx << "timestep : " << timestep << ", lxlylz=" << tmp_lxlylz << endl;
    fsBx << "timestep : " << timestep << ", uxlylz=" << tmp_uxlylz << endl;
    fsBx << "timestep : " << timestep << ", lxuylz=" << tmp_lxuylz << endl;
    fsBx << "timestep : " << timestep << ", uxuylz=" << tmp_uxuylz << endl;
    fsBx << "timestep : " << timestep << ", lxlyuz=" << tmp_lxlyuz << endl;
    fsBx << "timestep : " << timestep << ", uxlyuz=" << tmp_uxlyuz << endl;
    fsBx << "timestep : " << timestep << ", lxuyuz=" << tmp_lxuyuz << endl;
    fsBx << "timestep : " << timestep << ", uxuyuz=" << tmp_uxuyuz << endl;
    fsBx << "timestep : " << timestep << ", mu=" << mu << endl;
    fsBx << "timestep : " << timestep << ", Bx->at(8065)=" << mu*0.5*(m_vCurrentHx->at(8065) + m_vNewHx->at(8065)) << endl;
    fsBx << endl;
    
    
    fsBy << "timestep : " << timestep << ", cellX=" << cellX << ", cellY=" << cellY << ", cellZ=" << cellZ << endl;
    fsBy << "timestep : " << timestep << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << endl;
    fsBy << "timestep : " << timestep << ", rdx=" << rdx << ", rdy=" << rdy << ", rdz=" << rdz << endl;
    fsBy << "timestep : " << timestep << ", dx=" << dx << endl;
    fsBy << "timestep : " << timestep << ", dy=" << dy << endl;
    fsBy << "timestep : " << timestep << ", by=" << by << endl;
    fsBy << "timestep : " << timestep << ", lxlylzBy=" << lxlylzBy << endl;
    fsBy << "timestep : " << timestep << ", uxlylzBy=" << uxlylzBy << endl;
    fsBy << "timestep : " << timestep << ", lxuylzBy=" << lxuylzBy << endl;
    fsBy << "timestep : " << timestep << ", uxuylzBy=" << uxuylzBy << endl;
    fsBy << "timestep : " << timestep << ", lxlyuzBy=" << lxlyuzBy << endl;
    fsBy << "timestep : " << timestep << ", uxlyuzBy=" << uxlyuzBy << endl;
    fsBy << "timestep : " << timestep << ", lxuyuzBy=" << lxuyuzBy << endl;
    fsBy << "timestep : " << timestep << ", uxuyuzBy=" << uxuyuzBy << endl;
    
    fsBz << "timestep : " << timestep << ", cellX=" << cellX << ", cellY=" << cellY << ", cellZ=" << cellZ << endl;
    fsBz << "timestep : " << timestep << ", dx=" << dx << ", dy=" << dy << ", dz=" << dz << endl;
    fsBz << "timestep : " << timestep << ", rdx=" << rdx << ", rdy=" << rdy << ", rdz=" << rdz << endl;
    fsBz << "timestep : " << timestep << ", dx=" << dx << endl;
    fsBz << "timestep : " << timestep << ", dy=" << dy << endl;
    fsBz << "timestep : " << timestep << ", bz=" << bz << endl;
    fsBz << "timestep : " << timestep << ", lxlylzBz=" << lxlylzBz << endl;
    fsBz << "timestep : " << timestep << ", uxlylzBz=" << uxlylzBz << endl;
    fsBz << "timestep : " << timestep << ", lxuylzBz=" << lxuylzBz << endl;
    fsBz << "timestep : " << timestep << ", uxuylzBz=" << uxuylzBz << endl;
    fsBz << "timestep : " << timestep << ", lxlyuzBz=" << lxlyuzBz << endl;
    fsBz << "timestep : " << timestep << ", uxlyuzBz=" << uxlyuzBz << endl;
    fsBz << "timestep : " << timestep << ", lxuyuzBz=" << lxuyuzBz << endl;
    fsBz << "timestep : " << timestep << ", uxuyuzBz=" << uxuyuzBz << endl;



    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fsBx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << mu*0.5*(m_vCurrentHx->at(index) + m_vNewHx->at(index)) << endl;
          fsBy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << mu*0.5*(m_vCurrentHy->at(index) + m_vNewHy->at(index)) << endl;
          fsBz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << mu*0.5*(m_vCurrentHz->at(index) + m_vNewHz->at(index)) << endl;
        }
      }
    }
    fsBx.close();
    fsBy.close();
    fsBz.close();
    //_exit(0);
  }
  //*/



  /*
  if(timestep == 2) {
    int i,j,k,index;
    double mu = getMu(0,0,0,0);
    stringstream fnEx;
    stringstream fnEy;
    stringstream fnEz;
    stringstream fnBx;
    stringstream fnBy;
    stringstream fnBz;
    stringstream fn2Ex;
    stringstream fn2Ey;
    stringstream fn2Ez;
    stringstream fn2Bx;
    stringstream fn2By;
    stringstream fn2Bz;
    fnEx << "serial_current_Ex.txt";
    fnEy << "serial_current_Ey.txt";
    fnEz << "serial_current_Ez.txt";
    fnBx << "serial_current_Hx.txt";
    fnBy << "serial_current_Hy.txt";
    fnBz << "serial_current_Hz.txt";
    fn2Ex << "serial_new_Ex.txt";
    fn2Ey << "serial_new_Ey.txt";
    fn2Ez << "serial_new_Ez.txt";
    fn2Bx << "serial_new_Hx.txt";
    fn2By << "serial_new_Hy.txt";
    fn2Bz << "serial_new_Hz.txt";
    ofstream fsEx(fnEx.str().c_str());
    ofstream fsEy(fnEy.str().c_str());
    ofstream fsEz(fnEz.str().c_str());
    ofstream fsBx(fnBx.str().c_str());
    ofstream fsBy(fnBy.str().c_str());
    ofstream fsBz(fnBz.str().c_str());
    ofstream fs2Ex(fn2Ex.str().c_str());
    ofstream fs2Ey(fn2Ey.str().c_str());
    ofstream fs2Ez(fn2Ez.str().c_str());
    ofstream fs2Bx(fn2Bx.str().c_str());
    ofstream fs2By(fn2By.str().c_str());
    ofstream fs2Bz(fn2Bz.str().c_str());

    fsEx.precision(60);
    fsEy.precision(60);
    fsEz.precision(60);
    fsBx.precision(60);
    fsBy.precision(60);
    fsBz.precision(60);
    fs2Ex.precision(60);
    fs2Ey.precision(60);
    fs2Ez.precision(60);
    fs2Bx.precision(60);
    fs2By.precision(60);
    fs2Bz.precision(60);

    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEx->at(index) << endl;
          fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEy->at(index) << endl;
          fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEz->at(index) << endl;
          fsBx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentHx->at(index) << endl;
          fsBy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentHy->at(index) << endl;
          fsBz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentHz->at(index) << endl;
          fs2Ex << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewEx->at(index) << endl;
          fs2Ey << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewEy->at(index) << endl;
          fs2Ez << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewEz->at(index) << endl;
          fs2Bx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewHx->at(index) << endl;
          fs2By << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewHy->at(index) << endl;
          fs2Bz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewHz->at(index) << endl;
        }
      }
    }
    fsEx.close();
    fsEy.close();
    fsEz.close();
    fsBx.close();
    fsBy.close();
    fsBz.close();
    fs2Ex.close();
    fs2Ey.close();
    fs2Ez.close();
    fs2Bx.close();
    fs2By.close();
    fs2Bz.close();
    //_exit(0);
  }
  //*/




  return EM_SUCCESS;
}

double SerialFieldSolver::getInitialEx(double x, double y, double z) {return 0.0;}
double SerialFieldSolver::getInitialEy(double x, double y, double z) {return 0.0;}
double SerialFieldSolver::getInitialEz(double x, double y, double z) {return 0.0;}

double SerialFieldSolver::getInitialHx(double x, double y, double z) {return 0.0;}
double SerialFieldSolver::getInitialHy(double x, double y, double z) {return 0.0;}
double SerialFieldSolver::getInitialHz(double x, double y, double z) {return 0.0;}



////////////////////////////////////////////////////////////////////////////////
// End : SerialFieldSolver
////////////////////////////////////////////////////////////////////////////////







 

////////////////////////////////////////////////////////////////////////////////
// Start : SerialTestSolver
////////////////////////////////////////////////////////////////////////////////
//*
const double SerialTestSolver::mu = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
const double SerialTestSolver::epsilon = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
const double SerialTestSolver::frequency = 2.45e9;
const double SerialTestSolver::omega = 2*PI*SerialTestSolver::frequency;
//*/


SerialTestSolver::SerialTestSolver(double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
: SerialFieldSolver(dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ)
{
}


void SerialTestSolver::initializeSolver() {
  
  double a = m_fUpperX - m_fLowerX;
  double b = m_fUpperY - m_fLowerY;

  m_fPIa = PI/a;
  m_fPIb = PI/b;

  m_fBeta = sqrt(omega*omega*mu*epsilon - m_fPIa*m_fPIa - m_fPIb*m_fPIb);
  m_fHsquare = m_fPIa*m_fPIa + m_fPIb*m_fPIb;
  m_fH = sqrt(m_fHsquare);

  SerialFieldSolver::initializeSolver();
}



double SerialTestSolver::getBoundaryXLowerEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*m_fLowerX) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryXLowerEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*m_fLowerX) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryXLowerEz(double x, double y, double z, double t) {
  return sin(m_fPIa*m_fLowerX) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryXUpperEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*m_fUpperX) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryXUpperEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*m_fUpperX) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryXUpperEz(double x, double y, double z, double t) {
  return sin(m_fPIa*m_fUpperX) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryYLowerEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*m_fLowerY) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryYLowerEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*m_fLowerY) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryYLowerEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*m_fLowerY) * cos(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryYUpperEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*m_fUpperY) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryYUpperEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*m_fUpperY) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryYUpperEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*m_fUpperY) * cos(omega*t - m_fBeta*z);
}

double SerialTestSolver::getBoundaryZLowerEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*m_fLowerZ);
}

double SerialTestSolver::getBoundaryZLowerEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*m_fLowerZ);
}

double SerialTestSolver::getBoundaryZLowerEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*m_fLowerZ);
}

double SerialTestSolver::getBoundaryZUpperEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*m_fUpperZ);
}

double SerialTestSolver::getBoundaryZUpperEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*m_fUpperZ);
}

double SerialTestSolver::getBoundaryZUpperEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*m_fUpperZ);
}




double SerialTestSolver::getInitialEx(double x, double y, double z) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(-m_fBeta*z);
}

double SerialTestSolver::getInitialEy(double x, double y, double z) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(-m_fBeta*z);
}

double SerialTestSolver::getInitialEz(double x, double y, double z) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(-m_fBeta*z);
}

double SerialTestSolver::getInitialHx(double x, double y, double z) {
  //return -(omega*epsilon/m_fHsquare)*m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*0.5*m_dt - m_fBeta*z);
  return -(omega*epsilon/m_fHsquare)*m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(- m_fBeta*z);
}

double SerialTestSolver::getInitialHy(double x, double y, double z) {
  //return (omega*epsilon/m_fHsquare)*m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*0.5*m_dt - m_fBeta*z);
  return (omega*epsilon/m_fHsquare)*m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(- m_fBeta*z);
}

double SerialTestSolver::getInitialHz(double x, double y, double z) {return 0.0;}






double SerialTestSolver::getExactEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getExactEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getExactEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*z);
}

double SerialTestSolver::getExactHx(double x, double y, double z, double t) {
  return -(omega*epsilon/m_fHsquare)*m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getExactHy(double x, double y, double z, double t) {
  return (omega*epsilon/m_fHsquare)*m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double SerialTestSolver::getExactHz(double x, double y, double z, double t) {return 0.0;}


void SerialTestSolver::evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {


  size_t index;

  /*
#ifndef _MSC_VER
  size_t i,j,k;
#else
  int i,j,k;
#endif 
  */
  int i,j,k;

  double difference_ex = 0.0;
  double difference_ey = 0.0;
  double difference_ez = 0.0;
  double max_difference_ex = 0.0;
  double max_difference_ey = 0.0;
  double max_difference_ez = 0.0;

  double difference_hx = 0.0;
  double difference_hy = 0.0;
  double difference_hz = 0.0;
  double max_difference_hx = 0.0;
  double max_difference_hy = 0.0;
  double max_difference_hz = 0.0;

  int sumCountEx = 0;
  int sumCountEy = 0;
  int sumCountEz = 0;
  int sumCountHx = 0;
  int sumCountHy = 0;
  int sumCountHz = 0;


  switch(m_iErrorMode) {
    case -1: // Absolute L_1 norm.
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;


            difference_ex = fabs( (*m_vCurrentEx)[index] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) );
            difference_ey = fabs( (*m_vCurrentEy)[index] - getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time) );
            difference_ez = fabs( (*m_vCurrentEz)[index] - getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );

            max_difference_ex += difference_ex;
            max_difference_ey += difference_ey;
            max_difference_ez += difference_ez;

            difference_hx = fabs( (*m_vCurrentHx)[index] - getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_hy = fabs( (*m_vCurrentHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_hz = fabs( (*m_vCurrentHz)[index] - getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time + 0.5*m_dt) );

            max_difference_hx += difference_hx;
            max_difference_hy += difference_hy;
            max_difference_hz += difference_hz;

          }
        }
      }
      break;

    case -2: // Absolute L_2 norm
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;

            


            difference_ex = fabs( (*m_vCurrentEx)[index] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) );
            difference_ey = fabs( (*m_vCurrentEy)[index] - getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time) );
            difference_ez = fabs( (*m_vCurrentEz)[index] - getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );

            max_difference_ex += difference_ex*difference_ex;
            max_difference_ey += difference_ey*difference_ey;
            max_difference_ez += difference_ez*difference_ez;

            difference_hx = fabs( (*m_vCurrentHx)[index] - getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_hy = fabs( (*m_vCurrentHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_hz = fabs( (*m_vCurrentHz)[index] - getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time + 0.5*m_dt) );

            max_difference_hx += difference_hx*difference_hx;
            max_difference_hy += difference_hy*difference_hy;
            max_difference_hz += difference_hz*difference_hz;
            
          }
        }
      }
      max_difference_ex = sqrt(max_difference_ex);
      max_difference_ey = sqrt(max_difference_ey);
      max_difference_ez = sqrt(max_difference_ez);
      max_difference_hx = sqrt(max_difference_hx);
      max_difference_hy = sqrt(max_difference_hy);
      max_difference_hz = sqrt(max_difference_hz);
      break;

    case 2: // Relative L_2 norm
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;

            double exactEx = getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time);
            double exactEy = getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time);
            double exactEz = getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time);
            double exactHx = getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt);
            double exactHy = getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt);
            double exactHz = getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time + 0.5*m_dt);

            if(fabs(exactEx) > threshold) {
              difference_ex = fabs( ((*m_vCurrentEx)[index]) - exactEx ) / fabs(exactEx);
              max_difference_ex += difference_ex*difference_ex;
              sumCountEx++;
            }

            if(fabs(exactEy) > threshold) {
              difference_ey = fabs( ((*m_vCurrentEy)[index]) - exactEy ) / fabs(exactEy);
              max_difference_ey += difference_ey*difference_ey;
              sumCountEy++;
            }

            if(fabs(exactEz) > threshold) {
              difference_ez = fabs( ((*m_vCurrentEz)[index]) - exactEz ) / fabs(exactEz);
              max_difference_ez += difference_ez*difference_ez;
              sumCountEz++;
            }

            if(fabs(exactHx) > threshold) {
              difference_hx = fabs( ((*m_vCurrentHx)[index]) - exactHx ) / fabs(exactHx);
              max_difference_hx += difference_hx*difference_hx;
              sumCountHx++;
            }

            if(fabs(exactHy) > threshold) {
              difference_hy = fabs( ((*m_vCurrentHy)[index]) - exactHy ) / fabs(exactHy);
              max_difference_hy += difference_hy*difference_hy;
              sumCountHy++;
            }

            if(fabs(exactHz) > threshold) {
              difference_hz = fabs( ((*m_vCurrentHz)[index]) - exactHz ) / fabs(exactHz);
              max_difference_hz += difference_hz*difference_hz;
              sumCountHz++;
            }


            
          }
        }
      }
      max_difference_ex = sqrt(max_difference_ex);
      max_difference_ey = sqrt(max_difference_ey);
      max_difference_ez = sqrt(max_difference_ez);
      max_difference_hx = sqrt(max_difference_hx);
      max_difference_hy = sqrt(max_difference_hy);
      max_difference_hz = sqrt(max_difference_hz);
      break;

    case 0:
    default:
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;


            difference_ex = fabs( (*m_vCurrentEx)[index] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) );
            difference_ey = fabs( (*m_vCurrentEy)[index] - getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time) );
            difference_ez = fabs( (*m_vCurrentEz)[index] - getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );

            max_difference_ex = max_difference_ex  > difference_ex ? max_difference_ex : difference_ex;
            max_difference_ey = max_difference_ey  > difference_ey ? max_difference_ey : difference_ey;
            max_difference_ez = max_difference_ez  > difference_ez ? max_difference_ez : difference_ez;

            difference_hx = fabs( (*m_vCurrentHx)[index] - getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_hy = fabs( (*m_vCurrentHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_hz = fabs( (*m_vCurrentHz)[index] - getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time + 0.5*m_dt) );

            max_difference_hx = max_difference_hx  > difference_hx ? max_difference_hx : difference_hx;
            max_difference_hy = max_difference_hy  > difference_hy ? max_difference_hy : difference_hy;
            max_difference_hz = max_difference_hz  > difference_hz ? max_difference_hz : difference_hz;

          }
        }
      }
      break;

  }


  /*
  diff_sol_ex << endl;
  diff_sol_ey << endl;
  diff_sol_ez << endl;
  diff_sol_hx << endl;
  diff_sol_hy << endl;
  diff_sol_hz << endl;

  diff_sol_ex.close();
  diff_sol_ey.close();
  diff_sol_ez.close();
  diff_sol_hx.close();
  diff_sol_hy.close();
  diff_sol_hz.close();
  */

  //

  if(m_iErrorMode > 0) { // relative error mode. infinite error mode has only absolute error mode.
    Ex = (sumCountEx == 0) ? 0 : max_difference_ex/sumCountEx;
    Ey = (sumCountEy == 0) ? 0 : max_difference_ey/sumCountEy;
    Ez = (sumCountEz == 0) ? 0 : max_difference_ez/sumCountEz;

    Hx = (sumCountHx == 0) ? 0 : max_difference_hx/sumCountHx;
    Hy = (sumCountHy == 0) ? 0 : max_difference_hy/sumCountHy;
    Hz = (sumCountHz == 0) ? 0 : max_difference_hz/sumCountHz;
  } else if (m_iErrorMode < 0) { // absolute error mode.
    size_t countNum = (m_iGridSizeX-1)*(m_iGridSizeY-1)*(m_iGridSizeZ-1);
    Ex = max_difference_ex/countNum;
    Ey = max_difference_ey/countNum;
    Ez = max_difference_ez/countNum;

    Hx = max_difference_hx/countNum;
    Hy = max_difference_hy/countNum;
    Hz = max_difference_hz/countNum;
  } else { // infinite error mode.
    Ex = max_difference_ex;
    Ey = max_difference_ey;
    Ez = max_difference_ez;

    Hx = max_difference_hx;
    Hy = max_difference_hy;
    Hz = max_difference_hz;
  }
}


////////////////////////////////////////////////////////////////////////////////
// End : SerialTestSolver
////////////////////////////////////////////////////////////////////////////////




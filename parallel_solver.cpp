/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 11 2012
 * 
 *      
 */



#include <mpi.h>
#include <iostream>
#include <omp.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>



#include "parallel_solver.h"
#include "constants.h"
#include "em_error.h"


using namespace std;









////////////////////////////////////////////////////////////////////////////////
// Start : MPIFieldSolver
////////////////////////////////////////////////////////////////////////////////

MPIFieldSolver::MPIFieldSolver(int rank, int procNumX, int procNumY, int procNumZ, double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ) 
: m_iRank(rank), m_iProcNumX(procNumX), m_iProcNumY(procNumY), m_iProcNumZ(procNumZ), m_dt(dt), m_fGlobalLowerX(leftX), m_fGlobalUpperX(rightX), m_iGridSizeX((gridNumX/procNumX)+1), m_fGlobalLowerY(leftY), m_fGlobalUpperY(rightY), m_iGridSizeY((gridNumY/procNumY)+1), m_fGlobalLowerZ(leftZ), m_fGlobalUpperZ(rightZ), m_iGridSizeZ((gridNumZ/procNumZ)+1), m_iErrorMode(0)
{
  m_iTotalProcNum = m_iProcNumX*m_iProcNumY*m_iProcNumZ;
  m_dx = (rightX - leftX) / gridNumX;
  m_dy = (rightY - leftY) / gridNumY;
  m_dz = (rightZ - leftZ) / gridNumZ;

  m_iGlobalGridSizeX = gridNumX + 1;
  m_iGlobalGridSizeY = gridNumY + 1;
  m_iGlobalGridSizeZ = gridNumZ + 1;

  m_vGlobalEx = new vector<double>(m_iGlobalGridSizeX*m_iGlobalGridSizeY*m_iGlobalGridSizeZ);
  m_vGlobalEy = new vector<double>(m_iGlobalGridSizeX*m_iGlobalGridSizeY*m_iGlobalGridSizeZ);
  m_vGlobalEz = new vector<double>(m_iGlobalGridSizeX*m_iGlobalGridSizeY*m_iGlobalGridSizeZ);

  m_vGlobalBx = new vector<double>(m_iGlobalGridSizeX*m_iGlobalGridSizeY*m_iGlobalGridSizeZ);
  m_vGlobalBy = new vector<double>(m_iGlobalGridSizeX*m_iGlobalGridSizeY*m_iGlobalGridSizeZ);
  m_vGlobalBz = new vector<double>(m_iGlobalGridSizeX*m_iGlobalGridSizeY*m_iGlobalGridSizeZ);
  
  m_vSendBuffer = 0;
  m_vRecvBuffer = 0;

  
}

MPIFieldSolver::~MPIFieldSolver() {
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

  delete m_vGlobalEx;
  delete m_vGlobalEy;
  delete m_vGlobalEz;

  delete m_vGlobalBx;
  delete m_vGlobalBy;
  delete m_vGlobalBz;

  /*
  delete m_vBoundaryLowerX;
  delete m_vBoundaryUpperX;
  delete m_vBoundaryLowerY;
  delete m_vBoundaryUpperY;
  delete m_vBoundaryLowerZ;
  delete m_vBoundaryUpperZ;
  */

  /*
  delete m_vMuX;
  delete m_vMuY;
  delete m_vMuZ;
  delete m_vEpsilonX;
  delete m_vEpsilonY;
  delete m_vEpsilonZ;
  delete m_vSigmaX;
  delete m_vSigmaY;
  delete m_vSigmaZ;
  */


  /*
  if(m_aBuf_lower_x_Ey != 0) delete[] m_aBuf_lower_x_Ey;
  if(m_aBuf_lower_x_Ez != 0) delete[] m_aBuf_lower_x_Ez;
  if(m_aBuf_upper_x_Ey != 0) delete[] m_aBuf_upper_x_Ey;
  if(m_aBuf_upper_x_Ez != 0) delete[] m_aBuf_upper_x_Ez;

  if(m_aBuf_lower_y_Ex != 0) delete[] m_aBuf_lower_y_Ex;
  if(m_aBuf_lower_y_Ez != 0) delete[] m_aBuf_lower_y_Ez;
  if(m_aBuf_upper_y_Ex != 0) delete[] m_aBuf_upper_y_Ex;
  if(m_aBuf_upper_y_Ez != 0) delete[] m_aBuf_upper_y_Ez;

  if(m_aBuf_lower_z_Ex != 0) delete[] m_aBuf_lower_z_Ex;
  if(m_aBuf_lower_z_Ey != 0) delete[] m_aBuf_lower_z_Ey;
  if(m_aBuf_upper_z_Ex != 0) delete[] m_aBuf_upper_z_Ex;
  if(m_aBuf_upper_z_Ey != 0) delete[] m_aBuf_upper_z_Ey;


  if(m_aBuf_lower_x_Hx != 0) delete[] m_aBuf_lower_x_Hx;
  if(m_aBuf_upper_x_Hx != 0) delete[] m_aBuf_upper_x_Hx;
  if(m_aBuf_lower_y_Hy != 0) delete[] m_aBuf_lower_y_Hy;
  if(m_aBuf_upper_y_Hy != 0) delete[] m_aBuf_upper_y_Hy;
  if(m_aBuf_lower_z_Hz != 0) delete[] m_aBuf_lower_z_Hz;
  if(m_aBuf_upper_z_Hz != 0) delete[] m_aBuf_upper_z_Hz;
  */


}


void MPIFieldSolver::initializeSolver() {
  //cout << "MPISolver::initializeSolver() : start " << endl;
  
  m_fMu = getMu(0,0,0,0);
  m_fEpsilon = getEpsilon(0,0,0,0);
  m_fSigma = getSigma(0,0,0,0);

  double temp;
  if(m_fGlobalLowerX > m_fGlobalUpperX) {
    temp = m_fGlobalLowerX;
    m_fGlobalLowerX = m_fGlobalUpperX;
    m_fGlobalUpperX = temp;
  }
  if(m_fGlobalLowerY > m_fGlobalUpperY) {
    temp = m_fGlobalLowerY;
    m_fGlobalLowerY = m_fGlobalUpperY;
    m_fGlobalUpperY = temp;
  }
  if(m_fGlobalLowerZ > m_fGlobalUpperZ) {
    temp = m_fGlobalLowerZ;
    m_fGlobalLowerZ = m_fGlobalUpperZ;
    m_fGlobalUpperZ = temp;
  }

  int indexZ = m_iRank / (m_iProcNumX*m_iProcNumY);
  int temp2 = m_iRank % (m_iProcNumX*m_iProcNumY);
  int indexY = temp2 / m_iProcNumX;
  int indexX = temp2 % m_iProcNumX;

  if(indexX == 0) {
    m_iLowerXNeighbor = FieldSolver::LOWER_X;
  } else {
    m_iLowerXNeighbor = m_iRank - 1;
  }
  if(indexX == (m_iProcNumX-1)) {
    m_iUpperXNeighbor = FieldSolver::UPPER_X;
  } else {
    m_iUpperXNeighbor = m_iRank + 1;
  }

  if(indexY == 0) {
    m_iLowerYNeighbor = FieldSolver::LOWER_Y;
  } else {
    m_iLowerYNeighbor = m_iRank - m_iProcNumX;
  }
  if(indexY == (m_iProcNumY-1)) {
    m_iUpperYNeighbor = FieldSolver::UPPER_Y;
  } else {
    m_iUpperYNeighbor = m_iRank + m_iProcNumX;
  }

  if(indexZ == 0) {
    m_iLowerZNeighbor = FieldSolver::LOWER_Z;
  } else {
    m_iLowerZNeighbor = m_iRank - m_iProcNumX*m_iProcNumY;
  }
  if(indexZ == (m_iProcNumZ-1)) {
    m_iUpperZNeighbor = FieldSolver::UPPER_Z;
  } else {
    m_iUpperZNeighbor = m_iRank + m_iProcNumX*m_iProcNumY;
  }



  double cellLengthX = (m_fGlobalUpperX - m_fGlobalLowerX) / m_iProcNumX;
  double cellLengthY = (m_fGlobalUpperY - m_fGlobalLowerY) / m_iProcNumY;
  double cellLengthZ = (m_fGlobalUpperZ - m_fGlobalLowerZ) / m_iProcNumZ;

  m_fLowerX = m_fGlobalLowerX + cellLengthX*indexX;
  m_fUpperX = m_fGlobalLowerX + cellLengthX*(indexX+1);
  m_fLowerY = m_fGlobalLowerY + cellLengthY*indexY;
  m_fUpperY = m_fGlobalLowerY + cellLengthY*(indexY+1);
  m_fLowerZ = m_fGlobalLowerZ + cellLengthZ*indexZ;
  m_fUpperZ = m_fGlobalLowerZ + cellLengthZ*(indexZ+1);


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

  /*
  m_vMuX = new vector<double>(m_iTotalDomainMemorySize);
  m_vMuY = new vector<double>(m_iTotalDomainMemorySize);
  m_vMuZ = new vector<double>(m_iTotalDomainMemorySize);
  m_vEpsilonX = new vector<double>(m_iTotalDomainMemorySize);
  m_vEpsilonY = new vector<double>(m_iTotalDomainMemorySize);
  m_vEpsilonZ = new vector<double>(m_iTotalDomainMemorySize);
  m_vSigmaX = new vector<double>(m_iTotalDomainMemorySize);
  m_vSigmaY = new vector<double>(m_iTotalDomainMemorySize);
  m_vSigmaZ = new vector<double>(m_iTotalDomainMemorySize);
  */

  int index, index2;
  // Initialize mu, epsilon, sigma
  /*
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
  */
  
  initGlobalFieldDomain();

  

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
        (*m_vNewHx)[index] = (*m_vCurrentHx)[index] + (0.5*m_dt/m_fMu) * ( ((*m_vCurrentEy)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEy)[index])/m_dz - ((*m_vCurrentEz)[index+m_iGridSizeX] - (*m_vCurrentEz)[index])/m_dy );
        (*m_vNewHy)[index] = (*m_vCurrentHy)[index] + (0.5*m_dt/m_fMu) * ( ((*m_vCurrentEz)[index+1] - (*m_vCurrentEz)[index])/m_dx - ((*m_vCurrentEx)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEx)[index])/m_dz );
        (*m_vNewHz)[index] = (*m_vCurrentHz)[index] + (0.5*m_dt/m_fMu) * ( ((*m_vCurrentEx)[index+m_iGridSizeX] - (*m_vCurrentEx)[index])/m_dy - ((*m_vCurrentEy)[index+1] - (*m_vCurrentEy)[index])/m_dx );
      }
    }
  }



  MPI_Barrier(MPI_COMM_WORLD);

  //*
  ////////////////////////////////////////////////////////////////////////////////
  // Start : Communication of H
  ////////////////////////////////////////////////////////////////////////////////
  
  
  double *buf_lower_x_Hx = 0;
  double *buf_lower_x_Hy = 0;
  double *buf_lower_x_Hz = 0;
  double *buf_upper_x_Hx = 0;
  double *buf_upper_x_Hy = 0;
  double *buf_upper_x_Hz = 0;
  double *buf_lower_y_Hx = 0;
  double *buf_lower_y_Hy = 0;
  double *buf_lower_y_Hz = 0;
  double *buf_upper_y_Hx = 0;
  double *buf_upper_y_Hy = 0;
  double *buf_upper_y_Hz = 0;
  double *buf_lower_z_Hx = 0;
  double *buf_lower_z_Hy = 0;
  double *buf_lower_z_Hz = 0;
  double *buf_upper_z_Hx = 0;
  double *buf_upper_z_Hy = 0;
  double *buf_upper_z_Hz = 0;

  int buf_size_x_H = (m_iGridSizeY)*(m_iGridSizeZ);
  int buf_size_y_H = (m_iGridSizeX)*(m_iGridSizeZ);
  int buf_size_z_H = (m_iGridSizeX)*(m_iGridSizeY);
  
  
  MPI_Request requestXHx;
  MPI_Request requestXHy;
  MPI_Request requestXHz;
  MPI_Request requestYHx;
  MPI_Request requestYHy;
  MPI_Request requestYHz;
  MPI_Request requestZHx;
  MPI_Request requestZHy;
  MPI_Request requestZHz;
  MPI_Status statusXHx;
  MPI_Status statusXHy;
  MPI_Status statusXHz;
  MPI_Status statusYHx;
  MPI_Status statusYHy;
  MPI_Status statusYHz;
  MPI_Status statusZHx;
  MPI_Status statusZHy;
  MPI_Status statusZHz;
  
  // Start : Sending Hx, Hy, Hz to Lower(Left) Cell.
  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    buf_upper_x_Hx = new double[buf_size_x_H];
    buf_upper_x_Hy = new double[buf_size_x_H];
    buf_upper_x_Hz = new double[buf_size_x_H];
    MPI_Irecv(buf_upper_x_Hx, buf_size_x_H, MPI_DOUBLE, m_iUpperXNeighbor, 0, MPI_COMM_WORLD, &requestXHx);
    MPI_Irecv(buf_upper_x_Hy, buf_size_x_H, MPI_DOUBLE, m_iUpperXNeighbor, 1, MPI_COMM_WORLD, &requestXHy);
    MPI_Irecv(buf_upper_x_Hz, buf_size_x_H, MPI_DOUBLE, m_iUpperXNeighbor, 2, MPI_COMM_WORLD, &requestXHz);
  }


  if(m_iLowerXNeighbor >= 0) { // Not Lower(Left) Boundary along x-coordinates.
    buf_lower_x_Hx = new double[buf_size_x_H];
    buf_lower_x_Hy = new double[buf_size_x_H];
    buf_lower_x_Hz = new double[buf_size_x_H];

    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index = j + (m_iGridSizeY)*k;
        index2 = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_x_Hx[index] = (*m_vNewHx)[index2];
        buf_lower_x_Hy[index] = (*m_vNewHy)[index2];
        buf_lower_x_Hz[index] = (*m_vNewHz)[index2];
      }
    }

    MPI_Send(buf_lower_x_Hx, buf_size_x_H, MPI_DOUBLE, m_iLowerXNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_x_Hy, buf_size_x_H, MPI_DOUBLE, m_iLowerXNeighbor, 1, MPI_COMM_WORLD);
    MPI_Send(buf_lower_x_Hz, buf_size_x_H, MPI_DOUBLE, m_iLowerXNeighbor, 2, MPI_COMM_WORLD);
  }


  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    MPI_Wait(&requestXHx, &statusXHx);
    MPI_Wait(&requestXHy, &statusXHy);
    MPI_Wait(&requestXHz, &statusXHz);
    
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index2 = j + (m_iGridSizeY)*k;
        index = (m_iGridSizeX - 1) + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = buf_upper_x_Hx[index2];
        (*m_vNewHy)[index] = buf_upper_x_Hy[index2];
        (*m_vNewHz)[index] = buf_upper_x_Hz[index2];
      }
    }
  }



  MPI_Barrier(MPI_COMM_WORLD);



  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    buf_upper_y_Hx = new double[buf_size_y_H];
    buf_upper_y_Hy = new double[buf_size_y_H];
    buf_upper_y_Hz = new double[buf_size_y_H];
    MPI_Irecv(buf_upper_y_Hx, buf_size_y_H, MPI_DOUBLE, m_iUpperYNeighbor, 0, MPI_COMM_WORLD, &requestYHx);
    MPI_Irecv(buf_upper_y_Hy, buf_size_y_H, MPI_DOUBLE, m_iUpperYNeighbor, 1, MPI_COMM_WORLD, &requestYHy);
    MPI_Irecv(buf_upper_y_Hz, buf_size_y_H, MPI_DOUBLE, m_iUpperYNeighbor, 2, MPI_COMM_WORLD, &requestYHz);
  }


  if(m_iLowerYNeighbor >= 0) { // Not Lower(Left) Boundary along y-coordinates.
    buf_lower_y_Hx = new double[buf_size_y_H];
    buf_lower_y_Hy = new double[buf_size_y_H];
    buf_lower_y_Hz = new double[buf_size_y_H];

    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + (m_iGridSizeX)*k;
        index2 = i + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_y_Hx[index] = (*m_vNewHx)[index2];
        buf_lower_y_Hy[index] = (*m_vNewHy)[index2];
        buf_lower_y_Hz[index] = (*m_vNewHz)[index2];
      }
    }

    MPI_Send(buf_lower_y_Hx, buf_size_y_H, MPI_DOUBLE, m_iLowerYNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_y_Hy, buf_size_y_H, MPI_DOUBLE, m_iLowerYNeighbor, 1, MPI_COMM_WORLD);
    MPI_Send(buf_lower_y_Hz, buf_size_y_H, MPI_DOUBLE, m_iLowerYNeighbor, 2, MPI_COMM_WORLD);
  }


  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    MPI_Wait(&requestYHx, &statusYHx);
    MPI_Wait(&requestYHy, &statusYHy);
    MPI_Wait(&requestYHz, &statusYHz);
    
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + (m_iGridSizeX)*k;
        index = i + m_iGridSizeX*(m_iGridSizeY-1) + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = buf_upper_y_Hx[index2];
        (*m_vNewHy)[index] = buf_upper_y_Hy[index2];
        (*m_vNewHz)[index] = buf_upper_y_Hz[index2];
      }
    }
  }



  MPI_Barrier(MPI_COMM_WORLD);



  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    buf_upper_z_Hx = new double[buf_size_z_H];
    buf_upper_z_Hy = new double[buf_size_z_H];
    buf_upper_z_Hz = new double[buf_size_z_H];
    MPI_Irecv(buf_upper_z_Hx, buf_size_z_H, MPI_DOUBLE, m_iUpperZNeighbor, 0, MPI_COMM_WORLD, &requestZHx);
    MPI_Irecv(buf_upper_z_Hy, buf_size_z_H, MPI_DOUBLE, m_iUpperZNeighbor, 1, MPI_COMM_WORLD, &requestZHy);
    MPI_Irecv(buf_upper_z_Hz, buf_size_z_H, MPI_DOUBLE, m_iUpperZNeighbor, 2, MPI_COMM_WORLD, &requestZHz);
  }

  
  if(m_iLowerZNeighbor >= 0) { // Not Lower(Left) Boundary along z-coordinates.
    buf_lower_z_Hx = new double[buf_size_z_H];
    buf_lower_z_Hy = new double[buf_size_z_H];
    buf_lower_z_Hz = new double[buf_size_z_H];

    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + (m_iGridSizeX)*j;
        index2 = i + m_iGridSizeX*j;
        buf_lower_z_Hx[index] = (*m_vNewHx)[index2];
        buf_lower_z_Hy[index] = (*m_vNewHy)[index2];
        buf_lower_z_Hz[index] = (*m_vNewHz)[index2];
      }
    }

    MPI_Send(buf_lower_z_Hx, buf_size_z_H, MPI_DOUBLE, m_iLowerZNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_z_Hy, buf_size_z_H, MPI_DOUBLE, m_iLowerZNeighbor, 1, MPI_COMM_WORLD);
    MPI_Send(buf_lower_z_Hz, buf_size_z_H, MPI_DOUBLE, m_iLowerZNeighbor, 2, MPI_COMM_WORLD);
  }

  
  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    MPI_Wait(&requestZHx, &statusZHx);
    MPI_Wait(&requestZHy, &statusZHy);
    MPI_Wait(&requestZHz, &statusZHz);
    
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + (m_iGridSizeX)*j;
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        (*m_vNewHx)[index] = buf_upper_z_Hx[index2];
        (*m_vNewHy)[index] = buf_upper_z_Hy[index2];
        (*m_vNewHz)[index] = buf_upper_z_Hz[index2];
      }
    }
  }
  
  // End : Sending Hx, Hy, Hz to Lower(Left) Cell.




  if(buf_lower_x_Hx != 0) delete[] buf_lower_x_Hx;
  if(buf_lower_x_Hy != 0) delete[] buf_lower_x_Hy;
  if(buf_lower_x_Hz != 0) delete[] buf_lower_x_Hz;

  if(buf_upper_x_Hx != 0) delete[] buf_upper_x_Hx;
  if(buf_upper_x_Hy != 0) delete[] buf_upper_x_Hy;
  if(buf_upper_x_Hz != 0) delete[] buf_upper_x_Hz;

  if(buf_lower_y_Hx != 0) delete[] buf_lower_y_Hx;
  if(buf_lower_y_Hy != 0) delete[] buf_lower_y_Hy;
  if(buf_lower_y_Hz != 0) delete[] buf_lower_y_Hz;

  if(buf_upper_y_Hx != 0) delete[] buf_upper_y_Hx;
  if(buf_upper_y_Hy != 0) delete[] buf_upper_y_Hy;
  if(buf_upper_y_Hz != 0) delete[] buf_upper_y_Hz;

  if(buf_lower_z_Hx != 0) delete[] buf_lower_z_Hx;
  if(buf_lower_z_Hy != 0) delete[] buf_lower_z_Hy;
  if(buf_lower_z_Hz != 0) delete[] buf_lower_z_Hz;

  if(buf_upper_z_Hx != 0) delete[] buf_upper_z_Hx;
  if(buf_upper_z_Hy != 0) delete[] buf_upper_z_Hy;
  if(buf_upper_z_Hz != 0) delete[] buf_upper_z_Hz;
  ////////////////////////////////////////////////////////////////////////////////
  // End : Communication of H
  ////////////////////////////////////////////////////////////////////////////////
  //*/


  updateH();











  ////////////////////////////////////////////////////////////////////////////////
  // End : loading Initial values
  ////////////////////////////////////////////////////////////////////////////////

  //cout << "MPISolver::initializeSolver() : end " << endl;



  /*
  int m_iBuf_size_x_Ey = (m_iGridSizeY-1)*m_iGridSizeZ;
  int m_iBuf_size_x_Ez = m_iGridSizeY*(m_iGridSizeZ-1);
  int m_iBuf_size_y_Ex = (m_iGridSizeX-1)*m_iGridSizeZ;
  int m_iBuf_size_y_Ez = m_iGridSizeX*(m_iGridSizeZ-1);
  int m_iBuf_size_z_Ex = (m_iGridSizeX-1)*m_iGridSizeY;
  int m_iBuf_size_z_Ey = m_iGridSizeX*(m_iGridSizeY-1);

  if(m_iLowerXNeighbor >= 0) { // Not Lower(Left) Boundary along x-coordinates.
    m_aBuf_lower_x_Ey = new double[m_iBuf_size_x_Ey];
    m_aBuf_lower_x_Ez = new double[m_iBuf_size_x_Ez];
  } else {
    m_aBuf_lower_x_Ey = 0;
    m_aBuf_lower_x_Ez = 0;
  }
  if(m_iLowerYNeighbor >= 0) { // Not Lower(Left) Boundary along y-coordinates.
    m_aBuf_lower_y_Ex = new double[m_iBuf_size_y_Ex];
    m_aBuf_lower_y_Ez = new double[m_iBuf_size_y_Ez];
  } else {
    m_aBuf_lower_y_Ex = 0;
    m_aBuf_lower_y_Ez = 0;
  }
  if(m_iLowerZNeighbor >= 0) { // Not Lower(Left) Boundary along z-coordinates.
    m_aBuf_lower_z_Ex = new double[m_iBuf_size_z_Ex];
    m_aBuf_lower_z_Ey = new double[m_iBuf_size_z_Ey];
  } else {
    m_aBuf_lower_z_Ex = 0;
    m_aBuf_lower_z_Ey = 0;
  }

  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    m_aBuf_upper_x_Ey = new double[m_iBuf_size_x_Ey];
    m_aBuf_upper_x_Ez = new double[m_iBuf_size_x_Ez];
  } else {
    m_aBuf_upper_x_Ey = 0;
    m_aBuf_upper_x_Ez = 0;
  }
  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    m_aBuf_upper_y_Ex = new double[m_iBuf_size_y_Ex];
    m_aBuf_upper_y_Ez = new double[m_iBuf_size_y_Ez];
  } else {
    m_aBuf_upper_y_Ex = 0;
    m_aBuf_upper_y_Ez = 0;
  }
  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    m_aBuf_upper_z_Ex = new double[m_iBuf_size_z_Ex];
    m_aBuf_upper_z_Ey = new double[m_iBuf_size_z_Ey];
  } else {
    m_aBuf_upper_z_Ex = 0;
    m_aBuf_upper_z_Ey = 0;
  }



  int m_iBuf_size_x_Hx = (m_iGridSizeY-1)*(m_iGridSizeZ-1);
  int m_iBuf_size_y_Hy = (m_iGridSizeX-1)*(m_iGridSizeZ-1);
  int m_iBuf_size_z_Hz = (m_iGridSizeX-1)*(m_iGridSizeY-1);

  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    m_aBuf_upper_x_Hx = new double[m_iBuf_size_x_Hx];
  } else {
    m_aBuf_upper_x_Hx = 0;
  }
  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    m_aBuf_upper_y_Hy = new double[m_iBuf_size_y_Hy];
  } else {
    m_aBuf_upper_y_Hy = 0;
  }
  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    m_aBuf_upper_z_Hz = new double[m_iBuf_size_z_Hz];
  } else {
    m_aBuf_upper_z_Hz = 0;
  }

  if(m_iLowerXNeighbor >= 0) { // Not Lower(Left) Boundary along x-coordinates.
    m_aBuf_lower_x_Hx = new double[m_iBuf_size_x_Hx];
  } else {
    m_aBuf_lower_x_Hx = 0;
  }
  if(m_iLowerYNeighbor >= 0) { // Not Lower(Left) Boundary along y-coordinates.
    m_aBuf_lower_y_Hy = new double[m_iBuf_size_y_Hy];
  } else {
    m_aBuf_lower_y_Hy = 0;
  }
  if(m_iLowerZNeighbor >= 0) { // Not Lower(Left) Boundary along z-coordinates.
    m_aBuf_lower_z_Hz = new double[m_iBuf_size_z_Hz];
  } else {
    m_aBuf_lower_z_Hz = 0;
  }
  */




}


int MPIFieldSolver::solve(double time) {
  //cout << "MPISolver::solve() : start " << endl;
  //double start, end;

  //start = omp_get_wtime();

  //double time = timestep*m_dt;
  int index;
  int index2;
  //int index3;

  //*
  double *buf_lower_x_Ey = 0;
  double *buf_lower_x_Ez = 0;
  double *buf_upper_x_Ey = 0;
  double *buf_upper_x_Ez = 0;
  double *buf_lower_y_Ex = 0;
  double *buf_lower_y_Ez = 0;
  double *buf_upper_y_Ex = 0;
  double *buf_upper_y_Ez = 0;
  double *buf_lower_z_Ex = 0;
  double *buf_lower_z_Ey = 0;
  double *buf_upper_z_Ex = 0;
  double *buf_upper_z_Ey = 0;

  int buf_size_x_Ey = (m_iGridSizeY-1)*m_iGridSizeZ;
  int buf_size_x_Ez = m_iGridSizeY*(m_iGridSizeZ-1);
  int buf_size_y_Ex = (m_iGridSizeX-1)*m_iGridSizeZ;
  int buf_size_y_Ez = m_iGridSizeX*(m_iGridSizeZ-1);
  int buf_size_z_Ex = (m_iGridSizeX-1)*m_iGridSizeY;
  int buf_size_z_Ey = m_iGridSizeX*(m_iGridSizeY-1);

  double *buf_lower_x_Hx = 0;
  double *buf_lower_x_Hy = 0;
  double *buf_lower_x_Hz = 0;
  double *buf_upper_x_Hx = 0;
  double *buf_upper_x_Hy = 0;
  double *buf_upper_x_Hz = 0;
  double *buf_lower_y_Hx = 0;
  double *buf_lower_y_Hy = 0;
  double *buf_lower_y_Hz = 0;
  double *buf_upper_y_Hx = 0;
  double *buf_upper_y_Hy = 0;
  double *buf_upper_y_Hz = 0;
  double *buf_lower_z_Hx = 0;
  double *buf_lower_z_Hy = 0;
  double *buf_lower_z_Hz = 0;
  double *buf_upper_z_Hx = 0;
  double *buf_upper_z_Hy = 0;
  double *buf_upper_z_Hz = 0;

  int buf_size_x_H = (m_iGridSizeY)*(m_iGridSizeZ);
  int buf_size_y_H = (m_iGridSizeX)*(m_iGridSizeZ);
  int buf_size_z_H = (m_iGridSizeX)*(m_iGridSizeY);
  //*/

  //end = omp_get_wtime();
  //cout << "Initializing Time " << end-start << " seconds." << endl;

  
  int i,j,k;


  //cout << "MPISolver::solve() : 1 " << endl;


  //start = omp_get_wtime();
  // Solving E.
  //double value;
  #pragma omp parallel for private(i,j,k,index)
  for(k=1 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=1 ; j<m_iGridSizeY-1 ; j++) {
      for(i=1 ; i<m_iGridSizeX-1 ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEx)[index] = (*m_vCurrentEx)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-m_iGridSizeX])/m_dy - ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-m_iGridSizeX*m_iGridSizeY])/m_dz ) + m_fSigma*(*m_vCurrentEx)[index];
        (*m_vNewEy)[index] = (*m_vCurrentEy)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX*m_iGridSizeY])/m_dz - ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-1])/m_dx ) + m_fSigma*(*m_vCurrentEy)[index];
        (*m_vNewEz)[index] = (*m_vCurrentEz)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-1])/m_dx - ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX])/m_dy ) + m_fSigma*(*m_vCurrentEz)[index];
      }
    }
  }

  // When i=0 on above setting.
  #pragma omp parallel for private(j,k,index)
  for(k=1 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=1 ; j<m_iGridSizeY-1 ; j++) {
      index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vNewEx)[index] = (*m_vCurrentEx)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-m_iGridSizeX])/m_dy - ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-m_iGridSizeX*m_iGridSizeY])/m_dz ) + m_fSigma*(*m_vCurrentEx)[index];
    }
  }

  // When j=0 on above setting.
  #pragma omp parallel for private(i,k,index)
  for(k=1 ; k<m_iGridSizeZ-1 ; k++) {
    for(i=1 ; i<m_iGridSizeX-1 ; i++) {
      index = i + m_iGridSizeX*m_iGridSizeY*k;
      (*m_vNewEy)[index] = (*m_vCurrentEy)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX*m_iGridSizeY])/m_dz - ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-1])/m_dx ) + m_fSigma*(*m_vCurrentEy)[index];
    }
  }
  //end = omp_get_wtime();
  //cout << "Solving E Time : when j=0 " << end-start << " seconds." << endl;
  //start = omp_get_wtime();
  // When k=0 on above setting.
  #pragma omp parallel for private(i,j,index)
  for(j=1 ; j<m_iGridSizeY-1 ; j++) {
    for(i=1 ; i<m_iGridSizeX-1 ; i++) {
      index = i + m_iGridSizeX*j;
      (*m_vNewEz)[index] = (*m_vCurrentEz)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-1])/m_dx - ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX])/m_dy ) + m_fSigma*(*m_vCurrentEz)[index];
    }
  }
  //end = omp_get_wtime();
  //cout << "Solving E Time : when k=0 " << end-start << " seconds." << endl;

  //end = omp_get_wtime();
  //cout << "Solving E Time " << end-start << " seconds." << endl;



  ////////////////////////////////////////////////////////////////////////////////
  // Start : loading upper boundary values of E
  ////////////////////////////////////////////////////////////////////////////////
  //start = omp_get_wtime();

  if(m_iUpperXNeighbor != FieldSolver::UPPER_X) {
    // X : Ey
    #pragma omp parallel for private(j,k,index)
    for(k=1 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index = (m_iGridSizeX - 1) +  m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k; // Upper(Right) side of x-coordinates.
        (*m_vNewEy)[index] = (*m_vCurrentEy)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX*m_iGridSizeY])/m_dz - ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-1])/m_dx ) + m_fSigma*(*m_vCurrentEy)[index];
      }
    }
    // X : Ez
    #pragma omp parallel for private(j,k,index)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=1 ; j<m_iGridSizeY ; j++ ) {
        index = (m_iGridSizeX - 1) +  m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k; // Upper(Right) side of x-coordinates.
        (*m_vNewEz)[index] = (*m_vCurrentEz)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-1])/m_dx - ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX])/m_dy ) + m_fSigma*(*m_vCurrentEz)[index];
      }
    }
  }

  if(m_iUpperYNeighbor != FieldSolver::UPPER_Y) {
    // Y : Ex
    #pragma omp parallel for private(i,k,index)
    for(k=1 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + m_iGridSizeX*(m_iGridSizeY-1) + m_iGridSizeX*m_iGridSizeY*k; // Upper(Right) side of y-coordinates.
        (*m_vNewEx)[index] = (*m_vCurrentEx)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-m_iGridSizeX])/m_dy - ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-m_iGridSizeX*m_iGridSizeY])/m_dz ) + m_fSigma*(*m_vCurrentEx)[index];
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=1 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*(m_iGridSizeY-1) + m_iGridSizeX*m_iGridSizeY*k; // Upper(Right) side of y-coordinates.
        (*m_vNewEz)[index] = (*m_vCurrentEz)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-1])/m_dx - ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX])/m_dy ) + m_fSigma*(*m_vCurrentEz)[index];
      }
    }
  }

  if(m_iUpperZNeighbor != FieldSolver::UPPER_Z) {
    // Z : Ex
    #pragma omp parallel for private(i,j,index)
    for(j=1 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1); // Upper(Right) side of z-coordinates.
        (*m_vNewEx)[index] = (*m_vCurrentEx)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-m_iGridSizeX])/m_dy - ((*m_vCurrentHy)[index] - (*m_vCurrentHy)[index-m_iGridSizeX*m_iGridSizeY])/m_dz ) + m_fSigma*(*m_vCurrentEx)[index];
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=1 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1); // Upper(Right) side of z-coordinates.
        (*m_vNewEy)[index] = (*m_vCurrentEy)[index] + (m_dt/m_fEpsilon) * ( ((*m_vCurrentHx)[index] - (*m_vCurrentHx)[index-m_iGridSizeX*m_iGridSizeY])/m_dz - ((*m_vCurrentHz)[index] - (*m_vCurrentHz)[index-1])/m_dx ) + m_fSigma*(*m_vCurrentEy)[index];
      }
    }
  }

  //////////////////////////////////////////////////////////////////

  if(m_iUpperXNeighbor == FieldSolver::UPPER_X) {
    // X : Ey
    #pragma omp parallel for private(j,k,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index2 = (m_iGridSizeX - 1) + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEy)[index2] = getBoundaryXUpperEy(m_fUpperX, m_fLowerY + j*m_dy + 0.5*m_dx, m_fLowerZ + k*m_dz, time);
      }
    }
    // X  : Ez
    #pragma omp parallel for private(j,k,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index2 = (m_iGridSizeX - 1) +  m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEz)[index2] = getBoundaryXUpperEz(m_fUpperX, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time);
      }
    }
  }

  if(m_iUpperYNeighbor == FieldSolver::UPPER_Y) {
    // Y : Ex
    #pragma omp parallel for private(i,k,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index2 = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEx)[index2] = getBoundaryYUpperEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fUpperY, m_fLowerZ + k*m_dz, time);
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEz)[index2] = getBoundaryYUpperEz(m_fLowerX + i*m_dx, m_fUpperY, m_fLowerZ + k*m_dz + 0.5*m_dz, time);
      }
    }
  }

  if(m_iUpperZNeighbor == FieldSolver::UPPER_Z) {
    // Z : Ex
    #pragma omp parallel for private(i,j,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index2 = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        (*m_vNewEx)[index2] = getBoundaryZUpperEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fUpperZ, time);
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index2)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        (*m_vNewEy)[index2] = getBoundaryZUpperEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fUpperZ, time);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // End : loading upper boundary values of E
  ////////////////////////////////////////////////////////////////////////////////



  /*
  if(timestep == 1) {
    stringstream fnEx;
    stringstream fnEy;
    stringstream fnEz;
    fnEx << "before_datatest_Ex_" << m_iRank << ".txt";
    fnEy << "before_datatest_Ey_" << m_iRank << ".txt";
    fnEz << "before_datatest_Ez_" << m_iRank << ".txt";
    ofstream fsEx(fnEx.str().c_str());
    ofstream fsEy(fnEy.str().c_str());
    ofstream fsEz(fnEz.str().c_str());

    fsEx << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsEx << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsEx << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsEx << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsEx << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsEx << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsEy << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsEy << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsEy << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsEy << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsEy << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsEy << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsEz << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsEz << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsEz << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsEz << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsEz << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsEz << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;



    double difference_ex;
    double difference_ey;
    double difference_ez;
    double index;
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY ; j++) {
        for(i=0 ; i<m_iGridSizeX ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          difference_ex = fabs( (*m_vNewEx)[index] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) );
          difference_ey = fabs( (*m_vNewEy)[index] - getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time) );
          difference_ez = fabs( (*m_vNewEz)[index] - getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );
          fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_ex << endl;
          fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_ey << endl;
          fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_ez << endl;

          //fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEx)[index] << endl;
          //fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEy)[index] << endl;
          //fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEz)[index] << endl;
        }
      }
    }
    fsEx.close();
    fsEy.close();
    fsEz.close();
  }
  //*/


  /*
  if(timestep == 2) {
    stringstream fnHx;
    stringstream fnHy;
    stringstream fnHz;
    fnHx << "before_datatest_Hx_" << m_iRank << ".txt";
    fnHy << "before_datatest_Hy_" << m_iRank << ".txt";
    fnHz << "before_datatest_Hz_" << m_iRank << ".txt";
    ofstream fsHx(fnHx.str().c_str());
    ofstream fsHy(fnHy.str().c_str());
    ofstream fsHz(fnHz.str().c_str());

    fsHx << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsHx << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsHx << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsHx << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsHx << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsHx << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsHy << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsHy << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsHy << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsHy << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsHy << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsHy << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsHz << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsHz << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsHz << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsHz << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsHz << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsHz << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;



    double difference_hx;
    double difference_hy;
    double difference_hz;
    double index;
    for(k=0 ; k<m_iGridSizeZ ; k++) {
      for(j=0 ; j<m_iGridSizeY ; j++) {
        for(i=0 ; i<m_iGridSizeX ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          difference_hx = fabs( (*m_vCurrentHx)[index] - getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time - 0.5*m_dt) );
          difference_hy = fabs( (*m_vCurrentHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time - 0.5*m_dt) );
          difference_hz = fabs( (*m_vCurrentHz)[index] - getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time - 0.5*m_dt) );
          fsHx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_hx << endl;
          fsHy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_hy << endl;
          fsHz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_hz << endl;

          //fsHx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentHx)[index] << endl;
          //fsHy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentHy)[index] << endl;
          //fsHz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentHz)[index] << endl;
        }
      }
    }
    fsHx.close();
    fsHy.close();
    fsHz.close();
  }
  */






  MPI_Barrier(MPI_COMM_WORLD);


  ////////////////////////////////////////////////////////////////////////////////
  // Start : Communication of E
  ////////////////////////////////////////////////////////////////////////////////
  //start = omp_get_wtime();

  MPI_Request requestXEy;
  MPI_Request requestXEz;
  MPI_Request requestYEx;
  MPI_Request requestYEz;
  MPI_Request requestZEx;
  MPI_Request requestZEy;
  MPI_Status statusXEy;
  MPI_Status statusXEz;
  MPI_Status statusYEx;
  MPI_Status statusYEz;
  MPI_Status statusZEx;
  MPI_Status statusZEy;

  MPI_Request requestXHx;
  MPI_Request requestXHy;
  MPI_Request requestXHz;
  MPI_Request requestYHx;
  MPI_Request requestYHy;
  MPI_Request requestYHz;
  MPI_Request requestZHx;
  MPI_Request requestZHy;
  MPI_Request requestZHz;
  MPI_Status statusXHx;
  MPI_Status statusXHy;
  MPI_Status statusXHz;
  MPI_Status statusYHx;
  MPI_Status statusYHy;
  MPI_Status statusYHz;
  MPI_Status statusZHx;
  MPI_Status statusZHy;
  MPI_Status statusZHz;


  // Start : Sending Ex, Ey, Ez to Upper(Right) Cell.
  if(m_iLowerXNeighbor < 0) { // Lower(Left) Boundary along x-coordinates.
    // X : Ey
    #pragma omp parallel for private(j,k,index)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEy)[index] = getBoundaryXLowerEy(m_fLowerX, m_fLowerY + j*m_dy + 0.5*m_dx, m_fLowerZ + k*m_dz, time);
      }
    }
    // X  : Ez
    #pragma omp parallel for private(j,k,index)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEz)[index] = getBoundaryXLowerEz(m_fLowerX, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time);
      }
    }
  } else {
    buf_lower_x_Ey = new double[buf_size_x_Ey];
    buf_lower_x_Ez = new double[buf_size_x_Ez];
    MPI_Irecv(buf_lower_x_Ey, buf_size_x_Ey, MPI_DOUBLE, m_iLowerXNeighbor, 0, MPI_COMM_WORLD, &requestXEy);
    MPI_Irecv(buf_lower_x_Ez, buf_size_x_Ez, MPI_DOUBLE, m_iLowerXNeighbor, 1, MPI_COMM_WORLD, &requestXEz);
  }

  
  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    buf_upper_x_Ey = new double[buf_size_x_Ey];
    buf_upper_x_Ez = new double[buf_size_x_Ez];

    // X : Ey
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index = j + (m_iGridSizeY-1)*k;
        index2 = (m_iGridSizeX - 1) + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        buf_upper_x_Ey[index] = (*m_vNewEy)[index2];
      }
    }
    // X  : Ez
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index = j + m_iGridSizeY*k;
        index2 = (m_iGridSizeX - 1) +  m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        buf_upper_x_Ez[index] = (*m_vNewEz)[index2];
      }
    }

    MPI_Send(buf_upper_x_Ey, buf_size_x_Ey, MPI_DOUBLE, m_iUpperXNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_upper_x_Ez, buf_size_x_Ez, MPI_DOUBLE, m_iUpperXNeighbor, 1, MPI_COMM_WORLD);
  }

  
  if(m_iLowerXNeighbor >= 0) { // Not Lower(Left) Boundary along x-coordinates.
    MPI_Wait(&requestXEy, &statusXEy);
    MPI_Wait(&requestXEz, &statusXEz);
    
    // X : Ey
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index2 = j + (m_iGridSizeY-1)*k;
        index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEy)[index] = buf_lower_x_Ey[index2];
      }
    }
    // X  : Ez
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index2 = j + m_iGridSizeY*k;
        index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEz)[index] = buf_lower_x_Ez[index2];
      }
    }
  }





  MPI_Barrier(MPI_COMM_WORLD);






  if(m_iLowerYNeighbor < 0) { // Lower(Left) Boundary along y-coordinates.
    // Y : Ex
    #pragma omp parallel for private(i,k,index)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEx)[index] = getBoundaryYLowerEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY, m_fLowerZ + k*m_dz, time);
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEz)[index] = getBoundaryYLowerEz(m_fLowerX + i*m_dx, m_fLowerY, m_fLowerZ + k*m_dz + 0.5*m_dz, time);
      }
    }
  } else {
    buf_lower_y_Ex = new double[buf_size_y_Ex];
    buf_lower_y_Ez = new double[buf_size_y_Ez];
    MPI_Irecv(buf_lower_y_Ex, buf_size_y_Ex, MPI_DOUBLE, m_iLowerYNeighbor, 0, MPI_COMM_WORLD, &requestYEx);
    MPI_Irecv(buf_lower_y_Ez, buf_size_y_Ez, MPI_DOUBLE, m_iLowerYNeighbor, 1, MPI_COMM_WORLD, &requestYEz);
  }


  
  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    buf_upper_y_Ex = new double[buf_size_y_Ex];
    buf_upper_y_Ez = new double[buf_size_y_Ez];

    // Y : Ex
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + (m_iGridSizeX-1)*k;
        index2 = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
        buf_upper_y_Ex[index] = (*m_vNewEx)[index2];
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*k;
        index2 = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
        buf_upper_y_Ez[index] = (*m_vNewEz)[index2];
      }
    }

    MPI_Send(buf_upper_y_Ex, buf_size_y_Ex, MPI_DOUBLE, m_iUpperYNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_upper_y_Ez, buf_size_y_Ez, MPI_DOUBLE, m_iUpperYNeighbor, 1, MPI_COMM_WORLD);
  }


  
  if(m_iLowerYNeighbor >= 0) { // Not Lower(Left) Boundary along y-coordinates.
    MPI_Wait(&requestYEx, &statusYEx);
    MPI_Wait(&requestYEz, &statusYEz);

    // Y : Ex
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index2 = i + (m_iGridSizeX-1)*k;
        index = i + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEx)[index] = buf_lower_y_Ex[index2];
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + m_iGridSizeX*k;
        index = i + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewEz)[index] = buf_lower_y_Ez[index2];
      }
    }
  }




  MPI_Barrier(MPI_COMM_WORLD);





  if(m_iLowerZNeighbor < 0) { // Lower(Left) Boundary along z-coordinates.
    // Z : Ex
    #pragma omp parallel for private(i,j,index)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + m_iGridSizeX*j;
        (*m_vNewEx)[index] = getBoundaryZLowerEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ, time);
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*j;
        (*m_vNewEy)[index] = getBoundaryZLowerEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ, time);
      }
    }
  } else {
    buf_lower_z_Ex = new double[buf_size_z_Ex];
    buf_lower_z_Ey = new double[buf_size_z_Ey];
    MPI_Irecv(buf_lower_z_Ex, buf_size_z_Ex, MPI_DOUBLE, m_iLowerZNeighbor, 0, MPI_COMM_WORLD, &requestZEx);
    MPI_Irecv(buf_lower_z_Ey, buf_size_z_Ey, MPI_DOUBLE, m_iLowerZNeighbor, 1, MPI_COMM_WORLD, &requestZEy);
  }


  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    buf_upper_z_Ex = new double[buf_size_z_Ex];
    buf_upper_z_Ey = new double[buf_size_z_Ey];

    // Z : Ex
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + (m_iGridSizeX-1)*j;
        index2 = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        buf_upper_z_Ex[index] = (*m_vNewEx)[index2];
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*j;
        index2 = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        buf_upper_z_Ey[index] = (*m_vNewEy)[index2];
      }
    }

    MPI_Send(buf_upper_z_Ex, buf_size_z_Ex, MPI_DOUBLE, m_iUpperZNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_upper_z_Ey, buf_size_z_Ey, MPI_DOUBLE, m_iUpperZNeighbor, 1, MPI_COMM_WORLD);
  }


  if(m_iLowerZNeighbor >= 0) { // Not Lower(Left) Boundary along z-coordinates.
    MPI_Wait(&requestZEx, &statusZEx);
    MPI_Wait(&requestZEy, &statusZEy);
    
    // Z : Ex
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index2 = i + (m_iGridSizeX-1)*j;
        index = i + m_iGridSizeX*j;
        (*m_vNewEx)[index] = buf_lower_z_Ex[index2];
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + m_iGridSizeX*j;
        index = i + m_iGridSizeX*j;
        (*m_vNewEy)[index] = buf_lower_z_Ey[index2];
      }
    }
  }
  // End : Sending Ex, Ey, Ez to Upper(Right) Cell.



  //MPI_Barrier(MPI_COMM_WORLD);


  if(buf_lower_x_Ey != 0) delete[] buf_lower_x_Ey;
  if(buf_lower_x_Ez != 0) delete[] buf_lower_x_Ez;
  if(buf_upper_x_Ey != 0) delete[] buf_upper_x_Ey;
  if(buf_upper_x_Ez != 0) delete[] buf_upper_x_Ez;

  if(buf_lower_y_Ex != 0) delete[] buf_lower_y_Ex;
  if(buf_lower_y_Ez != 0) delete[] buf_lower_y_Ez;
  if(buf_upper_y_Ex != 0) delete[] buf_upper_y_Ex;
  if(buf_upper_y_Ez != 0) delete[] buf_upper_y_Ez;

  if(buf_lower_z_Ex != 0) delete[] buf_lower_z_Ex;
  if(buf_lower_z_Ey != 0) delete[] buf_lower_z_Ey;
  if(buf_upper_z_Ex != 0) delete[] buf_upper_z_Ex;
  if(buf_upper_z_Ey != 0) delete[] buf_upper_z_Ey;





    /*
    stringstream fn3Ex;
    stringstream fn3Ey;
    stringstream fn3Ez;
    fn3Ex << "right_after_m_vCurrent_Ex_" << m_iRank << ".txt";
    fn3Ey << "right_after_m_vCurrent_Ey_" << m_iRank << ".txt";
    fn3Ez << "right_after_m_vCurrent_Ez_" << m_iRank << ".txt";
    ofstream fs3Ex(fn3Ex.str().c_str());
    ofstream fs3Ey(fn3Ey.str().c_str());
    ofstream fs3Ez(fn3Ez.str().c_str());


    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
          fs3Ex << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEx->at(index) << endl;
          fs3Ey << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEy->at(index) << endl;
          fs3Ez << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEz->at(index) << endl;

          //fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEx)[index] << endl;
          //fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEy)[index] << endl;
          //fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEz)[index] << endl;
        }
      }
    }
    fs3Ex.close();
    fs3Ey.close();
    fs3Ez.close();





    stringstream fn5Ex;
    stringstream fn5Ey;
    stringstream fn5Ez;
    fn5Ex << "right_after_m_vNew_Ex_" << m_iRank << ".txt";
    fn5Ey << "right_after_m_vNew_Ey_" << m_iRank << ".txt";
    fn5Ez << "right_after_m_vNew_Ez_" << m_iRank << ".txt";
    ofstream fs5Ex(fn5Ex.str().c_str());
    ofstream fs5Ey(fn5Ey.str().c_str());
    ofstream fs5Ez(fn5Ez.str().c_str());


    fs5Ex << "m_iGridSizeX : " << m_iGridSizeX << endl;
    fs5Ex << "m_iGridSizeY : " << m_iGridSizeY << endl;
    fs5Ex << "m_iGridSizeZ : " << m_iGridSizeZ << endl;
    fs5Ex << "m_fEpsilon : " << m_fEpsilon << endl;
    fs5Ex << "m_fMu : " << m_fMu << endl;
    fs5Ex << "m_fSigma : " << m_fSigma << endl;
    fs5Ey << "m_iGridSizeX : " << m_iGridSizeX << endl;
    fs5Ey << "m_iGridSizeY : " << m_iGridSizeY << endl;
    fs5Ey << "m_iGridSizeZ : " << m_iGridSizeZ << endl;
    fs5Ey << "m_fEpsilon : " << m_fEpsilon << endl;
    fs5Ey << "m_fMu : " << m_fMu << endl;
    fs5Ey << "m_fSigma : " << m_fSigma << endl;
    fs5Ez << "m_iGridSizeX : " << m_iGridSizeX << endl;
    fs5Ez << "m_iGridSizeY : " << m_iGridSizeY << endl;
    fs5Ez << "m_iGridSizeZ : " << m_iGridSizeZ << endl;
    fs5Ez << "m_fEpsilon : " << m_fEpsilon << endl;
    fs5Ez << "m_fMu : " << m_fMu << endl;
    fs5Ez << "m_fSigma : " << m_fSigma << endl;
    
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
          fs5Ex << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewEx->at(index) << endl;
          fs5Ey << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewEy->at(index) << endl;
          fs5Ez << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vNewEz->at(index) << endl;

          //fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEx)[index] << endl;
          //fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEy)[index] << endl;
          //fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEz)[index] << endl;
        }
      }
    }
    fs5Ex.close();
    fs5Ey.close();
    fs5Ez.close();
    */







  //end = omp_get_wtime();
  //cout << "Communitation Time " << end-start << " seconds." << endl;
  ////////////////////////////////////////////////////////////////////////////////
  // End : Communication of E
  ////////////////////////////////////////////////////////////////////////////////


  updateE();

  /*
  if(timestep == 2) {
    stringstream fnEx;
    stringstream fnEy;
    stringstream fnEz;
    fnEx << "datatest_Ex_" << m_iRank << ".txt";
    fnEy << "datatest_Ey_" << m_iRank << ".txt";
    fnEz << "datatest_Ez_" << m_iRank << ".txt";
    ofstream fsEx(fnEx.str().c_str());
    ofstream fsEy(fnEy.str().c_str());
    ofstream fsEz(fnEz.str().c_str());

    fsEx << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsEx << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsEx << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsEx << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsEx << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsEx << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsEy << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsEy << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsEy << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsEy << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsEy << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsEy << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsEz << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsEz << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsEz << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsEz << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsEz << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsEz << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;



    double difference_ex;
    double difference_ey;
    double difference_ez;
    double index;
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY ; j++) {
        for(i=0 ; i<m_iGridSizeX ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          difference_ex = fabs( (*m_vCurrentEx)[index] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) );
          difference_ey = fabs( (*m_vCurrentEy)[index] - getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time) );
          difference_ez = fabs( (*m_vCurrentEz)[index] - getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );
          fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_ex << endl;
          fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_ey << endl;
          fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_ez << endl;

          //fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEx)[index] << "\t" << getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) << endl;
          //fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEy)[index] << "\t" << getExactEy(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time) << endl;
          //fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEz)[index] << "\t" << getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) << endl;
        }
      }
    }
    fsEx.close();
    fsEy.close();
    fsEz.close();
  }
  //*/




  // Solving H.
  //start = omp_get_wtime();
  #pragma omp parallel for private(i,j,k,index)
  for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
    for(j=0 ; j<m_iGridSizeY-1 ; j++) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++) {
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = (*m_vCurrentHx)[index] + (m_dt/m_fMu) * ( ((*m_vCurrentEy)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEy)[index])/m_dz - ((*m_vCurrentEz)[index+m_iGridSizeX] - (*m_vCurrentEz)[index])/m_dy );
        (*m_vNewHy)[index] = (*m_vCurrentHy)[index] + (m_dt/m_fMu) * ( ((*m_vCurrentEz)[index+1] - (*m_vCurrentEz)[index])/m_dx - ((*m_vCurrentEx)[index+m_iGridSizeX*m_iGridSizeY] - (*m_vCurrentEx)[index])/m_dz );
        (*m_vNewHz)[index] = (*m_vCurrentHz)[index] + (m_dt/m_fMu) * ( ((*m_vCurrentEx)[index+m_iGridSizeX] - (*m_vCurrentEx)[index])/m_dy - ((*m_vCurrentEy)[index+1] - (*m_vCurrentEy)[index])/m_dx );
      }
    }
  }



    /*
    if(timestep == 2) {
      stringstream fnNewHy;
      stringstream fnCurrentHy;
      stringstream fnCurrentEx;
      stringstream fnCurrentEz;
      stringstream fnCurrentEx3;
      stringstream fnCurrentEz1;
      
      fnNewHy << "special_datatest_newHy_" << m_iRank << ".txt";
      fnCurrentHy << "special_datatest_currentHy_" << m_iRank << ".txt";
      fnCurrentEx << "special_datatest_currentEx_" << m_iRank << ".txt";
      fnCurrentEz << "special_datatest_currentEz_" << m_iRank << ".txt";
      fnCurrentEx3 << "special_datatest_currentEx3_" << m_iRank << ".txt";
      fnCurrentEz1 << "special_datatest_currentEz1_" << m_iRank << ".txt";

      ofstream fsNewHy(fnNewHy.str().c_str());
      ofstream fsCurrentHy(fnCurrentHy.str().c_str());
      ofstream fsCurrentEx(fnCurrentEx.str().c_str());
      ofstream fsCurrentEz(fnCurrentEz.str().c_str());
      ofstream fsCurrentEx3(fnCurrentEx3.str().c_str());
      ofstream fsCurrentEz1(fnCurrentEz1.str().c_str());

      fsNewHy << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
      fsNewHy << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
      fsNewHy << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
      fsNewHy << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
      fsNewHy << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
      fsNewHy << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

      fsCurrentHy << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
      fsCurrentHy << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
      fsCurrentHy << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
      fsCurrentHy << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
      fsCurrentHy << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
      fsCurrentHy << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

      fsCurrentEx << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
      fsCurrentEx << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
      fsCurrentEx << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
      fsCurrentEx << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
      fsCurrentEx << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
      fsCurrentEx << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

      fsCurrentEz << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
      fsCurrentEz << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
      fsCurrentEz << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
      fsCurrentEz << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
      fsCurrentEz << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
      fsCurrentEz << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

      fsCurrentEx3 << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
      fsCurrentEx3 << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
      fsCurrentEx3 << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
      fsCurrentEx3 << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
      fsCurrentEx3 << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
      fsCurrentEx3 << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

      fsCurrentEz1 << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
      fsCurrentEz1 << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
      fsCurrentEz1 << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
      fsCurrentEz1 << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
      fsCurrentEz1 << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
      fsCurrentEz1 << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;


      double difference_new_hy;
      double difference_current_hy;
      double difference_current_ex;
      double difference_current_ez;
      double difference_current_ex3;
      double difference_current_ez1;

      double index1,index2,index3;
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
            index1 = index + 1;
            index3 = index + m_iGridSizeX*m_iGridSizeY;

            difference_new_hy = fabs( (*m_vNewHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
            difference_current_hy = fabs( (*m_vCurrentHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time - 0.5*m_dt) );
            double difference_current_ex = fabs( (*m_vCurrentEx)[index] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz, time) );
            double difference_current_ez = fabs( (*m_vCurrentEz)[index] - getExactEz(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );
            double difference_current_ex3 = fabs( (*m_vCurrentEx)[index3] - getExactEx(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + (k+1)*m_dz, time) );
            double difference_current_ez1 = fabs( (*m_vCurrentEz)[index1] - getExactEz(m_fLowerX + (i+1)*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time) );

            fsNewHy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_new_hy << endl;
            fsCurrentHy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_current_hy << endl;
            fsCurrentEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_current_ex << endl;
            fsCurrentEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_current_ez << endl;
            fsCurrentEx3 << "["<<index3<<"]\t["<<i<<"]["<<j<<"]["<<(k+1)<<"]\t" << difference_current_ex3 << endl;
            fsCurrentEz1 << "["<<index1<<"]\t["<<(i+1)<<"]["<<j<<"]["<<k<<"]\t" << difference_current_ez1 << endl;
          }
        }
      }
      fsNewHy.close();
      fsCurrentHy.close();
      fsCurrentEx.close();
      fsCurrentEz.close();
      fsCurrentEx3.close();
      fsCurrentEz1.close();
    }
    //*/








        
        

  MPI_Barrier(MPI_COMM_WORLD);

  //*
  ////////////////////////////////////////////////////////////////////////////////
  // Start : Communication of H
  ////////////////////////////////////////////////////////////////////////////////
  // Start : Sending Hx, Hy, Hz to Lower(Left) Cell.
  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    buf_upper_x_Hx = new double[buf_size_x_H];
    buf_upper_x_Hy = new double[buf_size_x_H];
    buf_upper_x_Hz = new double[buf_size_x_H];
    MPI_Irecv(buf_upper_x_Hx, buf_size_x_H, MPI_DOUBLE, m_iUpperXNeighbor, 0, MPI_COMM_WORLD, &requestXHx);
    MPI_Irecv(buf_upper_x_Hy, buf_size_x_H, MPI_DOUBLE, m_iUpperXNeighbor, 1, MPI_COMM_WORLD, &requestXHy);
    MPI_Irecv(buf_upper_x_Hz, buf_size_x_H, MPI_DOUBLE, m_iUpperXNeighbor, 2, MPI_COMM_WORLD, &requestXHz);
  }


  if(m_iLowerXNeighbor >= 0) { // Not Lower(Left) Boundary along x-coordinates.
    buf_lower_x_Hx = new double[buf_size_x_H];
    buf_lower_x_Hy = new double[buf_size_x_H];
    buf_lower_x_Hz = new double[buf_size_x_H];

    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index = j + (m_iGridSizeY)*k;
        index2 = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_x_Hx[index] = (*m_vNewHx)[index2];
        buf_lower_x_Hy[index] = (*m_vNewHy)[index2];
        buf_lower_x_Hz[index] = (*m_vNewHz)[index2];
      }
    }

    MPI_Send(buf_lower_x_Hx, buf_size_x_H, MPI_DOUBLE, m_iLowerXNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_x_Hy, buf_size_x_H, MPI_DOUBLE, m_iLowerXNeighbor, 1, MPI_COMM_WORLD);
    MPI_Send(buf_lower_x_Hz, buf_size_x_H, MPI_DOUBLE, m_iLowerXNeighbor, 2, MPI_COMM_WORLD);
  }


  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    MPI_Wait(&requestXHx, &statusXHx);
    MPI_Wait(&requestXHy, &statusXHy);
    MPI_Wait(&requestXHz, &statusXHz);
    
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index2 = j + (m_iGridSizeY)*k;
        index = (m_iGridSizeX - 1) + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = buf_upper_x_Hx[index2];
        (*m_vNewHy)[index] = buf_upper_x_Hy[index2];
        (*m_vNewHz)[index] = buf_upper_x_Hz[index2];
      }
    }
  }



  MPI_Barrier(MPI_COMM_WORLD);



  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    buf_upper_y_Hx = new double[buf_size_y_H];
    buf_upper_y_Hy = new double[buf_size_y_H];
    buf_upper_y_Hz = new double[buf_size_y_H];
    MPI_Irecv(buf_upper_y_Hx, buf_size_y_H, MPI_DOUBLE, m_iUpperYNeighbor, 0, MPI_COMM_WORLD, &requestYHx);
    MPI_Irecv(buf_upper_y_Hy, buf_size_y_H, MPI_DOUBLE, m_iUpperYNeighbor, 1, MPI_COMM_WORLD, &requestYHy);
    MPI_Irecv(buf_upper_y_Hz, buf_size_y_H, MPI_DOUBLE, m_iUpperYNeighbor, 2, MPI_COMM_WORLD, &requestYHz);
  }


  if(m_iLowerYNeighbor >= 0) { // Not Lower(Left) Boundary along y-coordinates.
    buf_lower_y_Hx = new double[buf_size_y_H];
    buf_lower_y_Hy = new double[buf_size_y_H];
    buf_lower_y_Hz = new double[buf_size_y_H];

    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + (m_iGridSizeX)*k;
        index2 = i + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_y_Hx[index] = (*m_vNewHx)[index2];
        buf_lower_y_Hy[index] = (*m_vNewHy)[index2];
        buf_lower_y_Hz[index] = (*m_vNewHz)[index2];
      }
    }

    MPI_Send(buf_lower_y_Hx, buf_size_y_H, MPI_DOUBLE, m_iLowerYNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_y_Hy, buf_size_y_H, MPI_DOUBLE, m_iLowerYNeighbor, 1, MPI_COMM_WORLD);
    MPI_Send(buf_lower_y_Hz, buf_size_y_H, MPI_DOUBLE, m_iLowerYNeighbor, 2, MPI_COMM_WORLD);
  }


  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    MPI_Wait(&requestYHx, &statusYHx);
    MPI_Wait(&requestYHy, &statusYHy);
    MPI_Wait(&requestYHz, &statusYHz);
    
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + (m_iGridSizeX)*k;
        index = i + m_iGridSizeX*(m_iGridSizeY-1) + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vNewHx)[index] = buf_upper_y_Hx[index2];
        (*m_vNewHy)[index] = buf_upper_y_Hy[index2];
        (*m_vNewHz)[index] = buf_upper_y_Hz[index2];
      }
    }
  }



  MPI_Barrier(MPI_COMM_WORLD);



  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    buf_upper_z_Hx = new double[buf_size_z_H];
    buf_upper_z_Hy = new double[buf_size_z_H];
    buf_upper_z_Hz = new double[buf_size_z_H];
    MPI_Irecv(buf_upper_z_Hx, buf_size_z_H, MPI_DOUBLE, m_iUpperZNeighbor, 0, MPI_COMM_WORLD, &requestZHx);
    MPI_Irecv(buf_upper_z_Hy, buf_size_z_H, MPI_DOUBLE, m_iUpperZNeighbor, 1, MPI_COMM_WORLD, &requestZHy);
    MPI_Irecv(buf_upper_z_Hz, buf_size_z_H, MPI_DOUBLE, m_iUpperZNeighbor, 2, MPI_COMM_WORLD, &requestZHz);
  }

  
  if(m_iLowerZNeighbor >= 0) { // Not Lower(Left) Boundary along z-coordinates.
    buf_lower_z_Hx = new double[buf_size_z_H];
    buf_lower_z_Hy = new double[buf_size_z_H];
    buf_lower_z_Hz = new double[buf_size_z_H];

    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + (m_iGridSizeX)*j;
        index2 = i + m_iGridSizeX*j;
        buf_lower_z_Hx[index] = (*m_vNewHx)[index2];
        buf_lower_z_Hy[index] = (*m_vNewHy)[index2];
        buf_lower_z_Hz[index] = (*m_vNewHz)[index2];
      }
    }

    MPI_Send(buf_lower_z_Hx, buf_size_z_H, MPI_DOUBLE, m_iLowerZNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_z_Hy, buf_size_z_H, MPI_DOUBLE, m_iLowerZNeighbor, 1, MPI_COMM_WORLD);
    MPI_Send(buf_lower_z_Hz, buf_size_z_H, MPI_DOUBLE, m_iLowerZNeighbor, 2, MPI_COMM_WORLD);
  }

  
  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    MPI_Wait(&requestZHx, &statusZHx);
    MPI_Wait(&requestZHy, &statusZHy);
    MPI_Wait(&requestZHz, &statusZHz);
    
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + (m_iGridSizeX)*j;
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        (*m_vNewHx)[index] = buf_upper_z_Hx[index2];
        (*m_vNewHy)[index] = buf_upper_z_Hy[index2];
        (*m_vNewHz)[index] = buf_upper_z_Hz[index2];
      }
    }
  }
  
  // End : Sending Hx, Hy, Hz to Lower(Left) Cell.




  if(buf_lower_x_Hx != 0) delete[] buf_lower_x_Hx;
  if(buf_lower_x_Hy != 0) delete[] buf_lower_x_Hy;
  if(buf_lower_x_Hz != 0) delete[] buf_lower_x_Hz;

  if(buf_upper_x_Hx != 0) delete[] buf_upper_x_Hx;
  if(buf_upper_x_Hy != 0) delete[] buf_upper_x_Hy;
  if(buf_upper_x_Hz != 0) delete[] buf_upper_x_Hz;

  if(buf_lower_y_Hx != 0) delete[] buf_lower_y_Hx;
  if(buf_lower_y_Hy != 0) delete[] buf_lower_y_Hy;
  if(buf_lower_y_Hz != 0) delete[] buf_lower_y_Hz;

  if(buf_upper_y_Hx != 0) delete[] buf_upper_y_Hx;
  if(buf_upper_y_Hy != 0) delete[] buf_upper_y_Hy;
  if(buf_upper_y_Hz != 0) delete[] buf_upper_y_Hz;

  if(buf_lower_z_Hx != 0) delete[] buf_lower_z_Hx;
  if(buf_lower_z_Hy != 0) delete[] buf_lower_z_Hy;
  if(buf_lower_z_Hz != 0) delete[] buf_lower_z_Hz;

  if(buf_upper_z_Hx != 0) delete[] buf_upper_z_Hx;
  if(buf_upper_z_Hy != 0) delete[] buf_upper_z_Hy;
  if(buf_upper_z_Hz != 0) delete[] buf_upper_z_Hz;
  ////////////////////////////////////////////////////////////////////////////////
  // End : Communication of H
  ////////////////////////////////////////////////////////////////////////////////
  //*/


  updateH();


    /*
    stringstream fn4Ex;
    stringstream fn4Ey;
    stringstream fn4Ez;
    fn4Ex << "after_m_vCurrent_Ex_" << m_iRank << ".txt";
    fn4Ey << "after_m_vCurrent_Ey_" << m_iRank << ".txt";
    fn4Ez << "after_m_vCurrent_Ez_" << m_iRank << ".txt";
    ofstream fs4Ex(fn4Ex.str().c_str());
    ofstream fs4Ey(fn4Ey.str().c_str());
    ofstream fs4Ez(fn4Ez.str().c_str());


    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
          fs4Ex << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEx->at(index) << endl;
          fs4Ey << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEy->at(index) << endl;
          fs4Ez << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vCurrentEz->at(index) << endl;

          //fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEx)[index] << endl;
          //fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEy)[index] << endl;
          //fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vNewEz)[index] << endl;
        }
      }
    }
    fs4Ex.close();
    fs4Ey.close();
    fs4Ez.close();
    */
  
  /*
  if(timestep == 1) {
    stringstream fnHx;
    stringstream fnHy;
    stringstream fnHz;
    fnHx << "datatest_Hx_" << m_iRank << ".txt";
    fnHy << "datatest_Hy_" << m_iRank << ".txt";
    fnHz << "datatest_Hz_" << m_iRank << ".txt";
    ofstream fsHx(fnHx.str().c_str());
    ofstream fsHy(fnHy.str().c_str());
    ofstream fsHz(fnHz.str().c_str());

    fsHx << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsHx << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsHx << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsHx << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsHx << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsHx << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsHy << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsHy << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsHy << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsHy << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsHy << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsHy << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;

    fsHz << "m_iLowerXNeighbor : " << m_iLowerXNeighbor << endl;
    fsHz << "m_iUpperXNeighbor : " << m_iUpperXNeighbor << endl;
    fsHz << "m_iLowerYNeighbor : " << m_iLowerYNeighbor << endl;
    fsHz << "m_iUpperYNeighbor : " << m_iUpperYNeighbor << endl;
    fsHz << "m_iLowerZNeighbor : " << m_iLowerZNeighbor << endl;
    fsHz << "m_iUpperZNeighbor : " << m_iUpperZNeighbor << endl;



    double difference_hx;
    double difference_hy;
    double difference_hz;
    double index;
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY ; j++) {
        for(i=0 ; i<m_iGridSizeX ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          difference_hx = fabs( (*m_vCurrentHx)[index] - getExactHx(m_fLowerX + i*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
          difference_hy = fabs( (*m_vCurrentHy)[index] - getExactHy(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy, m_fLowerZ + k*m_dz + 0.5*m_dz, time + 0.5*m_dt) );
          difference_hz = fabs( (*m_vCurrentHz)[index] - getExactHz(m_fLowerX + i*m_dx + 0.5*m_dx, m_fLowerY + j*m_dy + 0.5*m_dy, m_fLowerZ + k*m_dz, time + 0.5*m_dt) );
          fsHx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_hx << endl;
          fsHy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_hy << endl;
          fsHz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << difference_hz << endl;

          //fsHx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentHx)[index] << endl;
          //fsHy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentHy)[index] << endl;
          //fsHz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentHz)[index] << endl;
        }
      }
    }
    fsHx.close();
    fsHy.close();
    fsHz.close();
  }
  //*/






  //MPI_Barrier(MPI_COMM_WORLD);



  return 0;
}



void MPIFieldSolver::updateE() {
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




void MPIFieldSolver::updateH() {
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




int MPIFieldSolver::updateCurrent(double time, int cellX, int cellY, int cellZ, double jx_y0z0, double jx_y1z0, double jx_y0z1, double jx_y1z1, double jy_x0z0, double jy_x1z0, double jy_x0z1, double jy_x1z1, double jz_x0y0, double jz_x1y0, double jz_x0y1, double jz_x1y1) {
  int index = cellX + m_iGlobalGridSizeX*cellY + m_iGlobalGridSizeX*m_iGlobalGridSizeY*cellZ;

  m_vGlobalEx->at(index) -= m_dt*jx_y0z0/getEpsilon(cellX*m_dx + 0.5*m_dx, cellY*m_dy, cellZ*m_dz, time);
  m_vGlobalEx->at(index + m_iGlobalGridSizeX) -= m_dt*jx_y1z0/getEpsilon(cellX*m_dx + 0.5*m_dx, (cellY+1)*m_dy, cellZ*m_dz, time);
  m_vGlobalEx->at(index + m_iGlobalGridSizeX*m_iGlobalGridSizeY) -= m_dt*jx_y0z1/getEpsilon(cellX*m_dx + 0.5*m_dx, cellY*m_dy, (cellZ+1)*m_dz, time);
  m_vGlobalEx->at(index + m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY) -= m_dt*jx_y1z1/getEpsilon(cellX*m_dx + 0.5*m_dx, (cellY+1)*m_dy, (cellZ+1)*m_dz, time);

  m_vGlobalEy->at(index) -= m_dt*jy_x0z0/getEpsilon(cellX*m_dx, cellY*m_dy + 0.5*m_dy, cellZ*m_dz, time);
  m_vGlobalEy->at(index + 1) -= m_dt*jy_x1z0/getEpsilon((cellX+1)*m_dx, cellY*m_dy + 0.5*m_dy, cellZ*m_dz, time);
  m_vGlobalEy->at(index + m_iGlobalGridSizeX*m_iGlobalGridSizeY) -= m_dt*jy_x0z1/getEpsilon(cellX*m_dx, cellY*m_dy + 0.5*m_dy, (cellZ+1)*m_dz, time);
  m_vGlobalEy->at(index + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY) -= m_dt*jy_x1z1/getEpsilon((cellX+1)*m_dx, cellY*m_dy + 0.5*m_dy, (cellZ+1)*m_dz, time);

  m_vGlobalEz->at(index) -= m_dt*jz_x0y0/getEpsilon(cellX*m_dx, cellY*m_dy, cellZ*m_dz + 0.5*m_dz, time);
  m_vGlobalEz->at(index + 1) -= m_dt*jz_x1y0/getEpsilon((cellX+1)*m_dx, cellY*m_dy, cellZ*m_dz + 0.5*m_dz, time);
  m_vGlobalEz->at(index + m_iGlobalGridSizeX) -= m_dt*jz_x0y1/getEpsilon(cellX*m_dx, (cellY+1)*m_dy, cellZ*m_dz + 0.5*m_dz, time);
  m_vGlobalEz->at(index + 1 + m_iGlobalGridSizeX) -= m_dt*jz_x1y1/getEpsilon((cellX+1)*m_dx, (cellY+1)*m_dy, cellZ*m_dz + 0.5*m_dz, time);


  /*
  cout << "index : " << index << endl;
  cout << "Ez at <index> : " << m_vCurrentEz->at(index) << endl;
  cout << "Ez at <index + 1> : " << m_vCurrentEz->at(index + 1) << endl;
  cout << "Ez at <index + m_iGridSizeX> : " << m_vCurrentEz->at(index + m_iGridSizeX) << endl;
  cout << "Ez at <index + 1 + m_iGridSizeX> : " << m_vCurrentEz->at(index + 1 + m_iGridSizeX) <<"\n"<<endl;
  */

  return EM_SUCCESS;
} 




double MPIFieldSolver::getMu(double x, double y, double z, double time) {
  return 0.0;
}

double MPIFieldSolver::getEpsilon(double x, double y, double z, double time) {
  return 0.0;
}

double MPIFieldSolver::getSigma(double x, double y, double z, double time) {
  return 0.0;
}



void MPIFieldSolver::initGlobalFieldDomain() {

  m_vGolableFieldData.clear();

  size_t localMemSizeOfFieldComponent = (m_iGridSizeX - 1)*(m_iGridSizeY - 1)*(m_iGridSizeZ - 1);
  size_t blockMemSize = localMemSizeOfFieldComponent*6;

  BlockDomainField blockField;

  m_vGolableFieldData.resize(m_iTotalProcNum, blockField);


  if(m_vSendBuffer != 0) {
    delete[] m_vRecvBuffer;
  }
  if(m_vRecvBuffer != 0) {
    delete[] m_vRecvBuffer;
  }

  m_vSendBuffer = new double[blockMemSize];
  //cout << "Rank : " << m_iRank << " : m_vSendBuffer address : " << m_vSendBuffer << endl;
  m_vRecvBuffer = new double[m_iTotalProcNum*blockMemSize];
  //cout << "Rank : " << m_iRank << " : m_vRecvBuffer address : " << m_vRecvBuffer << endl;
  

  for(int i=0 ; i<m_iTotalProcNum ; ++i) {
    BlockDomainField& field = m_vGolableFieldData[i];
    field.ex = m_vRecvBuffer + i*blockMemSize;
    field.ey = m_vRecvBuffer + i*blockMemSize + localMemSizeOfFieldComponent;
    field.ez = m_vRecvBuffer + i*blockMemSize + localMemSizeOfFieldComponent*2;

    field.bx = m_vRecvBuffer + i*blockMemSize + localMemSizeOfFieldComponent*3;
    field.by = m_vRecvBuffer + i*blockMemSize + localMemSizeOfFieldComponent*4;
    field.bz = m_vRecvBuffer + i*blockMemSize + localMemSizeOfFieldComponent*5;
  }
  
  /*
  for(int i=0 ; i<m_iTotalProcNum ; ++i) {
    BlockDomainField& field = m_vGolableFieldData[i];
    cout << "Rank : " << m_iRank << ", m_vGolableFieldData["<<i<<"].ex address : " << m_vGolableFieldData[i].ex << endl;
    cout << "Rank : " << m_iRank << ", field.ex address : " << field.ex << endl;
  }
  */

}





void MPIFieldSolver::getGlobalFieldData(size_t timestep) {
  /*
  if(m_iRank == 0) {
    cout << "timestep=" << timestep << ", " << "Rank : " << m_iRank << ", " << "m_vGolableFieldData[0].ex address : " << (m_vGolableFieldData[0].ex) << endl;
  }
  */
  
  
  /*
  if(timestep == 2) {
    
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
    fnEx << "local_current_Ex_" << m_iRank << ".txt";
    fnEy << "local_current_Ey_" << m_iRank << ".txt";
    fnEz << "local_current_Ez_" << m_iRank << ".txt";
    fnBx << "local_current_Hx_" << m_iRank << ".txt";
    fnBy << "local_current_Hy_" << m_iRank << ".txt";
    fnBz << "local_current_Hz_" << m_iRank << ".txt";
    fn2Ex << "local_new_Ex_" << m_iRank << ".txt";
    fn2Ey << "local_new_Ey_" << m_iRank << ".txt";
    fn2Ez << "local_new_Ez_" << m_iRank << ".txt";
    fn2Bx << "local_new_Hx_" << m_iRank << ".txt";
    fn2By << "local_new_Hy_" << m_iRank << ".txt";
    fn2Bz << "local_new_Hz_" << m_iRank << ".txt";
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

    int i,j,k,index;
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
  
  
  

  size_t index, index2;
  

  int arrayLen = (m_iGridSizeX - 1)*(m_iGridSizeY - 1)*(m_iGridSizeZ - 1);
  int sendBufferSize = arrayLen*6;
  //size_t recvBufferSize = m_iTotalNumProc*sendBufferSize;

  //cout << "sendBufferSize : " << sendBufferSize << endl;

  for(int i=0 ; i<m_iGridSizeX-1 ; ++i) {
    for(int j=0 ; j<m_iGridSizeY-1 ; ++j) {
      for(int k=0 ; k<m_iGridSizeZ-1 ; ++k) {
        index = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
        index2 = i + j*(m_iGridSizeX-1) + k*(m_iGridSizeX-1)*(m_iGridSizeY-1);
        
        //cout << "i=" << i << ", j=" << j << ", k=" << k <<", index : " << index << ", index2 : " << index2 << endl;

        m_vSendBuffer[index2] = m_vCurrentEx->at(index);
        m_vSendBuffer[index2 + arrayLen] = m_vCurrentEy->at(index);
        m_vSendBuffer[index2 + arrayLen*2] = m_vCurrentEz->at(index);
        m_vSendBuffer[index2 + arrayLen*3] = m_fMu*0.5*(m_vCurrentHx->at(index) + m_vNewHx->at(index));
        m_vSendBuffer[index2 + arrayLen*4] = m_fMu*0.5*(m_vCurrentHy->at(index) + m_vNewHy->at(index));
        m_vSendBuffer[index2 + arrayLen*5] = m_fMu*0.5*(m_vCurrentHz->at(index) + m_vNewHz->at(index));
      }
    }
  }




  //cout << "here? - 1" << endl;
  MPI_Allgather(m_vSendBuffer, sendBufferSize, MPI_DOUBLE, m_vRecvBuffer, sendBufferSize, MPI_DOUBLE, MPI_COMM_WORLD);
  //cout << "here? - 2" << endl;
  
  int blockGridNumX = m_iGridSizeX - 1;
  int blockGridNumY = m_iGridSizeY - 1;
  int blockGridNumZ = m_iGridSizeZ - 1;
  int startX, startY, startZ;
  //int endX, endY, endZ;
  int block;
  int srcIndex, targetIndex;

  for(int i=0 ; i<m_iProcNumX ; i++) {
    for(int j=0 ; j<m_iProcNumY ; j++) {
      for(int k=0 ; k<m_iProcNumZ ; k++) {
        block = i + j*m_iProcNumX + k*m_iProcNumX*m_iProcNumY;
        startX = i*blockGridNumX;
        startY = j*blockGridNumY;
        startZ = k*blockGridNumZ;

        BlockDomainField& field = m_vGolableFieldData[block];
        //cout << "Rank : " << m_iRank << ", block=" << block << ", field.ex address : " << field.ex << endl;
        //cout << "Rank : " << m_iRank << ", " << "m_vGolableFieldData["<<block<<"].ex address : " << (m_vGolableFieldData[block].ex) << endl;
        for(int l=0 ; l<blockGridNumX ; ++l) {
          for(int m=0 ; m<blockGridNumY ; ++m) {
            for(int n=0 ; n<blockGridNumZ ; ++n) {
              srcIndex = l + m*blockGridNumX + n*blockGridNumX*blockGridNumY;
              //cout << "srcIndex : " << srcIndex << endl;
              targetIndex = (startX + l) + (startY + m)*m_iGlobalGridSizeX + (startZ + n)*m_iGlobalGridSizeX*m_iGlobalGridSizeY;
              //cout << "targetIndex : " << targetIndex << endl;
              m_vGlobalEx->at(targetIndex) = field.ex[srcIndex];
              m_vGlobalEy->at(targetIndex) = field.ey[srcIndex];
              m_vGlobalEz->at(targetIndex) = field.ez[srcIndex];
              m_vGlobalBx->at(targetIndex) = field.bx[srcIndex];
              m_vGlobalBy->at(targetIndex) = field.by[srcIndex];
              m_vGlobalBz->at(targetIndex) = field.bz[srcIndex];
            }
          }
        }

      }
    }
  }


    /*
    //int i,j,k;
    if(timestep == 2) {
      stringstream fnBx;
      stringstream fnBy;
      stringstream fnBz;
      fnBx << "global_Bx_" << m_iRank << ".txt";
      fnBy << "global_By_" << m_iRank << ".txt";
      fnBz << "global_Bz_" << m_iRank << ".txt";
      ofstream fsBx(fnBx.str().c_str());
      ofstream fsBy(fnBy.str().c_str());
      ofstream fsBz(fnBz.str().c_str());
  
      
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
            fsBx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBx->at(index) << endl;
            fsBy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBy->at(index) << endl;
            fsBz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBz->at(index) << endl;
          }
        }
      }
      fsBx.close();
      fsBy.close();
      fsBz.close();
      
      //MPI_Finalize();
    }
    //*/
    
    
    
    
    
    
    
    
    
    /*
    if(timestep == 3) {
      if(m_iRank == 0) {
        index = 14 + 12*m_iGridSizeX + 7*m_iGridSizeX*m_iGridSizeY;
        printf("m_vCurrentEy[%d] : %.20f\n", index+m_iGridSizeX*m_iGridSizeY, (*m_vCurrentEy)[index+m_iGridSizeX*m_iGridSizeY]);
        printf("m_vCurrentEy[%d] : %.20f\n", index, (*m_vCurrentEy)[index]);
        printf("m_vCurrentEz[%d] : %.20f\n", index+m_iGridSizeX, (*m_vCurrentEz)[index+m_iGridSizeX]);
        printf("m_vCurrentEz[%d] : %.20f\n", index, (*m_vCurrentEz)[index]);
        
        //cout << "m_vCurrentEy["<<index+m_iGridSizeX*m_iGridSizeY<<"] : " << (*m_vCurrentEy)[index+m_iGridSizeX*m_iGridSizeY] << endl;
        //cout << "m_vCurrentEy["<<index<<"] : " << (*m_vCurrentEy)[index] << endl;
        //cout << "m_vCurrentEz["<<index+m_iGridSizeX<<"] : " << (*m_vCurrentEz)[index+m_iGridSizeX] << endl;
        //cout << "m_vCurrentEz["<<index<<"] : " << (*m_vCurrentEz)[index] << endl;
        
        //printf("m_vNewEy[9122] : %.10f\n", (*m_vNewEy)[9122]);
        //printf("m_vNewEy[8033] : %.10f\n", (*m_vNewEy)[8033]);
        //printf("m_vNewEz[8066] : %.10f\n", (*m_vNewEz)[8066]);
        //printf("m_vNewEz[8033] : %.10f\n", (*m_vNewEz)[8033]);
      }
      
      
      stringstream fnEx;
      stringstream fnEy;
      stringstream fnEz;
      fnEx << "global_Ex_" << m_iRank << ".txt";
      fnEy << "global_Ey_" << m_iRank << ".txt";
      fnEz << "global_Ez_" << m_iRank << ".txt";
      ofstream fsEx(fnEx.str().c_str());
      ofstream fsEy(fnEy.str().c_str());
      ofstream fsEz(fnEz.str().c_str());
  
      
      for(k=0 ; k<m_iGlobalGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGlobalGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGlobalGridSizeX-1 ; i++) {
            index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
            fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalEx->at(index) << endl;
            fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalEy->at(index) << endl;
            fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalEz->at(index) << endl;
          }
        }
      }
      fsEx.close();
      fsEy.close();
      fsEz.close();
      
      
      
      
      
      
      stringstream fnBx;
      stringstream fnBy;
      stringstream fnBz;
      fnBx << "global_Bx_" << m_iRank << ".txt";
      fnBy << "global_By_" << m_iRank << ".txt";
      fnBz << "global_Bz_" << m_iRank << ".txt";
      ofstream fsBx(fnBx.str().c_str());
      ofstream fsBy(fnBy.str().c_str());
      ofstream fsBz(fnBz.str().c_str());
  
  
  
      for(k=0 ; k<m_iGlobalGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGlobalGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGlobalGridSizeX-1 ; i++) {
            index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
            fsBx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBx->at(index) << endl;
            fsBy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBy->at(index) << endl;
            fsBz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBz->at(index) << endl;
          }
        }
      }
      fsBx.close();
      fsBy.close();
      fsBz.close();
      
      
      MPI_Finalize();
    }
    //*/
    
    
    
    
    
    
    
    
    
    
    
    /*
    if(timestep == 5) {
      
      stringstream FNsendBufferEx;
      stringstream FNsendBufferEy;
      stringstream FNsendBufferEz;
      stringstream FNsendBufferBx;
      stringstream FNsendBufferBy;
      stringstream FNsendBufferBz;
      
      FNsendBufferEx << "sendBuffer_Ex_" << m_iRank << ".txt";
      FNsendBufferEy << "sendBuffer_Ey_" << m_iRank << ".txt";
      FNsendBufferEz << "sendBuffer_Ez_" << m_iRank << ".txt";
      FNsendBufferBx << "sendBuffer_Bx_" << m_iRank << ".txt";
      FNsendBufferBy << "sendBuffer_By_" << m_iRank << ".txt";
      FNsendBufferBz << "sendBuffer_Bz_" << m_iRank << ".txt";
      
      ofstream sendBufferEx(FNsendBufferEx.str().c_str());
      ofstream sendBufferEy(FNsendBufferEy.str().c_str());
      ofstream sendBufferEz(FNsendBufferEz.str().c_str());
      ofstream sendBufferBx(FNsendBufferBx.str().c_str());
      ofstream sendBufferBy(FNsendBufferBy.str().c_str());
      ofstream sendBufferBz(FNsendBufferBz.str().c_str());
  
  
  
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            //index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
            index2 = i + j*(m_iGridSizeX-1) + k*(m_iGridSizeX-1)*(m_iGridSizeY-1);
            
            sendBufferEx << "["<<index2<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vSendBuffer[index2] << endl;
            sendBufferEy << "["<<index2<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vSendBuffer[index2 + arrayLen] << endl;
            sendBufferEz << "["<<index2<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vSendBuffer[index2 + arrayLen*2] << endl;
            sendBufferBx << "["<<index2<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vSendBuffer[index2 + arrayLen*3] << endl;
            sendBufferBy << "["<<index2<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vSendBuffer[index2 + arrayLen*4] << endl;
            sendBufferBz << "["<<index2<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vSendBuffer[index2 + arrayLen*5] << endl;
          }
        }
      }
      
      sendBufferEx.close();
      sendBufferEy.close();
      sendBufferEz.close();
      sendBufferBx.close();
      sendBufferBy.close();
      sendBufferBz.close();
      
      
      
      
      
      
      
      
      stringstream fn2Bx;

      fn2Bx << "m_vGolableFieldData_Bx_" << m_iRank << ".txt";

      ofstream fs2Bx(fn2Bx.str().c_str());
  
      BlockDomainField& field = m_vGolableFieldData[m_iRank];
      fs2Bx << "Rank : " << m_iRank << ", " << "field.ex address : " << (field.ex) << endl;
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + (m_iGridSizeX-1)*j + (m_iGridSizeX-1)*(m_iGridSizeY-1)*k;
            fs2Bx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << field.bx[index] << endl;
          }
        }
      }
      fs2Bx.close();
      
      
      
      
      
      
      
      stringstream fn4Bx;

      fn4Bx << "m_vGolableFieldData_direct_Bx_" << m_iRank << ".txt";

      ofstream fs4Bx(fn4Bx.str().c_str());
  
      size_t localMemSizeOfFieldComponent = (m_iGridSizeX - 1)*(m_iGridSizeY - 1)*(m_iGridSizeZ - 1);
      size_t blockMemSize = localMemSizeOfFieldComponent*6;
      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + (m_iGridSizeX-1)*j + (m_iGridSizeX-1)*(m_iGridSizeY-1)*k;
            fs4Bx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vRecvBuffer[index + m_iRank*blockMemSize + 3*localMemSizeOfFieldComponent] << endl;
          }
        }
      }
      fs4Bx.close();
      
      
      
      
      
      
      
      
      stringstream fn3Bx;

      fn3Bx << "original_Bx_" << m_iRank << ".txt";

      ofstream fs3Bx(fn3Bx.str().c_str());

      for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGridSizeX-1 ; i++) {
            index = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
            fs3Bx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_fMu*0.5*(m_vCurrentHx->at(index) + m_vNewHx->at(index)) << endl;
          }
        }
      }
      fs3Bx.close();
      
      
      
      
      
      
      
      
      
      
      
      MPI_Finalize();
    }
    //*/




    /*
    if(timestep == 1) {
      stringstream fn2Ex;
      stringstream fn2Ey;
      stringstream fn2Ez;
      fn2Ex << "field_Ex_" << m_iRank << ".txt";
      fn2Ey << "field_Ey_" << m_iRank << ".txt";
      fn2Ez << "field_Ez_" << m_iRank << ".txt";
      ofstream fs2Ex(fn2Ex.str().c_str());
      ofstream fs2Ey(fn2Ey.str().c_str());
      ofstream fs2Ez(fn2Ez.str().c_str());
  
  
      for(k=0 ; k<m_iGlobalGridSizeZ-1 ; k++) {
        for(j=0 ; j<m_iGlobalGridSizeY-1 ; j++) {
          for(i=0 ; i<m_iGlobalGridSizeX-1 ; i++) {
            index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
            fs2Ex << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGolableFieldData[block].ex[index] << endl;
            fs2Ey << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGolableFieldData[block].ey[index] << endl;
            fs2Ez << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGolableFieldData[block].ez[index] << endl;
          }
        }
      }
      fs2Ex.close();
      fs2Ey.close();
      fs2Ez.close();
    }
    //*/



    //
    //exit(0);

}




void MPIFieldSolver::updateGlobalElectricIntensityByCurrent(size_t timestep) {
  
  /*
  if(timestep == 2) {
    size_t i,j,k,index;
    //stringstream fnEx;
    //stringstream fnEy;
    stringstream fnEz;
    //fnEx << "before_updateCurrent_global_Ex_" << m_iRank << ".txt";
    //fnEy << "before_updateCurrent_global_Ey_" << m_iRank << ".txt";
    fnEz << "before_updateCurrent_global_Ez_" << m_iRank << ".txt";
    //FILE* fsEx = fopen(fnEx.str().c_str(), "w");
    //FILE* fsEy = fopen(fnEy.str().c_str(), "w");
    FILE* fsEz = fopen(fnEz.str().c_str(), "w");
    
    fprintf(fsEz, "timestep : %d\n", timestep);
    
    for(k=0 ; k<m_iGlobalGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGlobalGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGlobalGridSizeX-1 ; i++) {
          index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
          //fprintf(fsEx, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, m_vGlobalEx->at(index));
          //fprintf(fsEy, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, m_vGlobalEy->at(index));
          fprintf(fsEz, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, m_vGlobalEz->at(index));
        }
      }
    }
    //fclose(fsEx);
    //fclose(fsEy);
    fclose(fsEz);
    
    
    //MPI_Finalize();
  }
  //*/
    
  int index, index2, block, targetIndex, srcIndex;
  int startX, startY, startZ;
  
  
  //MPI_Request sendRequest;
  //MPI_Request recvRequest;
  //MPI_Status recvStatus;
    
    
  int blockGridNumX = m_iGridSizeX - 1;
  int blockGridNumY = m_iGridSizeY - 1;
  int blockGridNumZ = m_iGridSizeZ - 1;
    
  int arrayLen = blockGridNumX*blockGridNumY*blockGridNumZ;
  
  
    
  // Data copy.
  for(int i=0 ; i<m_iProcNumX ; i++) {
    for(int j=0 ; j<m_iProcNumY ; j++) {
      for(int k=0 ; k<m_iProcNumZ ; k++) {
        block = i + j*m_iProcNumX + k*m_iProcNumX*m_iProcNumY;

        BlockDomainField& field = m_vGolableFieldData[block];// reuse of m_vGolableFieldData
        startX = i*blockGridNumX;
        startY = j*blockGridNumY;
        startZ = k*blockGridNumZ;

        for(int l=0 ; l<blockGridNumX ; ++l) {
          for(int m=0 ; m<blockGridNumY ; ++m) {
            for(int n=0 ; n<blockGridNumZ ; ++n) {
              targetIndex = l + m*blockGridNumX + n*blockGridNumX*blockGridNumY;
              srcIndex = (startX + l) + (startY + m)*m_iGlobalGridSizeX + (startZ + n)*m_iGlobalGridSizeX*m_iGlobalGridSizeY;
              field.ex[targetIndex] = m_vGlobalEx->at(srcIndex);
              field.ey[targetIndex] = m_vGlobalEy->at(srcIndex);
              field.ez[targetIndex] = m_vGlobalEz->at(srcIndex);
            }
          }
        }


      }
    }
  }
  







  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Reduction
  int sendBufferSize = arrayLen*3;
  double* sendBuffer; 
  double* recvBuffer;

  for(int i=0 ; i<m_iTotalProcNum ; ++i) {
    sendBuffer = m_vGolableFieldData[block].ex;
    recvBuffer = sendBuffer;

    if(m_iRank == 0) {
      MPI_Reduce(MPI_IN_PLACE, recvBuffer, sendBufferSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(sendBuffer, 0, sendBufferSize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Scattering
  int* displs = new int[m_iTotalProcNum];
  for(int i=0 ; i<m_iTotalProcNum ; ++i) {
    displs[i] = i*(sendBufferSize*2);
  }
  int* sendcounts = new int[m_iTotalProcNum];
  for(int i=0 ; i<m_iTotalProcNum ; ++i) {
    sendcounts[i] = sendBufferSize;
  }

  sendBuffer = m_vGolableFieldData[0].ex;
  recvBuffer = m_vGolableFieldData[m_iRank].ex;
  MPI_Scatterv(sendBuffer, sendcounts, displs, MPI_DOUBLE, recvBuffer, sendBufferSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  delete[] displs;
  delete[] sendcounts;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data copy
  BlockDomainField& field = m_vGolableFieldData[m_iRank];
  for(int i=0 ; i<blockGridNumX ; ++i) {
    for(int j=0 ; j<blockGridNumY ; ++j) {
      for(int k=0 ; k<blockGridNumZ ; ++k) {
        index = i + j*blockGridNumX + k*blockGridNumX*blockGridNumY;
        index2 = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
        (*m_vCurrentEx)[index2] += field.ex[index];
        (*m_vCurrentEy)[index2] += field.ey[index];
        (*m_vCurrentEz)[index2] += field.ez[index];
      }
    }
  }


















  
  /*
  // Self copy.
  BlockDomainField& field = m_vGolableFieldData[m_iRank];
  for(int i=0 ; i<blockGridNumX ; ++i) {
    for(int j=0 ; j<blockGridNumY ; ++j) {
      for(int k=0 ; k<blockGridNumZ ; ++k) {
        index = i + j*blockGridNumX + k*blockGridNumX*blockGridNumY;
        index2 = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
        (*m_vCurrentEx)[index2] += field.ex[index];
        (*m_vCurrentEy)[index2] += field.ey[index];
        (*m_vCurrentEz)[index2] += field.ez[index];
      }
    }
  }

  size_t sendBufferSize = arrayLen*3;

  // Receive
  for(int i=0 ; i<m_iTotalProcNum-1 ; ++i) {
    MPI_Irecv(recvBuffer, sendBufferSize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvRequest);
  }

  // Send
  for(int i=0 ; i<m_iProcNumX ; i++) {
    for(int j=0 ; j<m_iProcNumY ; j++) {
      for(int k=0 ; k<m_iProcNumZ ; k++) {
        block = i + j*m_iProcNumX + k*m_iProcNumX*m_iProcNumY;
        if(m_iRank != block) {
          sendBuffer = m_vGolableFieldData[block].ex;// reuse of m_vGolableFieldData

          MPI_Isend(sendBuffer, sendBufferSize, MPI_DOUBLE, block, 0, MPI_COMM_WORLD, &sendRequest);

        } // if

      }
    }
  }
  

  // Wait of Receive
  //field = m_vGolableFieldData[block];
  field = m_vGolableFieldData[m_iRank];
  //double* recvBuffer = m_vGolableFieldData[m_iRank].ex;
  for(int h=0 ; h<m_iTotalProcNum-1 ; ++h) {
    MPI_Wait(&recvRequest, &recvStatus);
    
    for(int i=0 ; i<blockGridNumX ; ++i) {
      for(int j=0 ; j<blockGridNumY ; ++j) {
        for(int k=0 ; k<blockGridNumZ ; ++k) {
          index = i + j*blockGridNumX + k*blockGridNumX*blockGridNumY;
          index2 = i + j*m_iGridSizeX + k*m_iGridSizeX*m_iGridSizeY;
          (*m_vCurrentEx)[index2] += field.ex[index];
          (*m_vCurrentEy)[index2] += field.ey[index];
          (*m_vCurrentEz)[index2] += field.ez[index];
        }
      }
    }

  }
  */
  
  
  /*
  if(timestep == 2) {
    size_t i,j,k,index;
    //stringstream fnEx;
    //stringstream fnEy;
    stringstream fnEz;
    //fnEx << "middle1_updateCurrent_global_Ex_" << m_iRank << ".txt";
    //fnEy << "middle1_updateCurrent_global_Ey_" << m_iRank << ".txt";
    fnEz << "middle1_updateCurrent_m_vCurrentEz_" << m_iRank << ".txt";
    //FILE* fsEx = fopen(fnEx.str().c_str(), "w");
    //FILE* fsEy = fopen(fnEy.str().c_str(), "w");
    FILE* fsEz = fopen(fnEz.str().c_str(), "w");
    
    fprintf(fsEz, "timestep : %d\n", timestep);
    
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fprintf(fsEz, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, (*m_vCurrentEz)[index]);
        }
      }
    }
    //fclose(fsEx);
    //fclose(fsEy);
    fclose(fsEz);
    
    
    //MPI_Finalize();
  }
  //*/


  MPI_Barrier(MPI_COMM_WORLD);

  

  ////////////////////////////////////////////////////////////////////////////////
  // Start : Communication of E
  ////////////////////////////////////////////////////////////////////////////////
  //start = omp_get_wtime();

  double *buf_lower_x_Ey = 0;
  double *buf_lower_x_Ez = 0;
  double *buf_upper_x_Ey = 0;
  double *buf_upper_x_Ez = 0;
  double *buf_lower_y_Ex = 0;
  double *buf_lower_y_Ez = 0;
  double *buf_upper_y_Ex = 0;
  double *buf_upper_y_Ez = 0;
  double *buf_lower_z_Ex = 0;
  double *buf_lower_z_Ey = 0;
  double *buf_upper_z_Ex = 0;
  double *buf_upper_z_Ey = 0;

  int buf_size_x_Ey = (m_iGridSizeY-1)*m_iGridSizeZ;
  int buf_size_x_Ez = m_iGridSizeY*(m_iGridSizeZ-1);
  int buf_size_y_Ex = (m_iGridSizeX-1)*m_iGridSizeZ;
  int buf_size_y_Ez = m_iGridSizeX*(m_iGridSizeZ-1);
  int buf_size_z_Ex = (m_iGridSizeX-1)*m_iGridSizeY;
  int buf_size_z_Ey = m_iGridSizeX*(m_iGridSizeY-1);

  MPI_Request requestXEy;
  MPI_Request requestXEz;
  MPI_Request requestYEx;
  MPI_Request requestYEz;
  MPI_Request requestZEx;
  MPI_Request requestZEy;
  MPI_Status statusXEy;
  MPI_Status statusXEz;
  MPI_Status statusYEx;
  MPI_Status statusYEz;
  MPI_Status statusZEx;
  MPI_Status statusZEy;

  int i,j,k;
  
  // Start : Sending Ex, Ey, Ez to Lower(Left) Cell.
  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    buf_upper_x_Ey = new double[buf_size_x_Ey];
    buf_upper_x_Ez = new double[buf_size_x_Ez];
    MPI_Irecv(buf_upper_x_Ey, buf_size_x_Ey, MPI_DOUBLE, m_iUpperXNeighbor, 0, MPI_COMM_WORLD, &requestXEy);
    MPI_Irecv(buf_upper_x_Ez, buf_size_x_Ez, MPI_DOUBLE, m_iUpperXNeighbor, 1, MPI_COMM_WORLD, &requestXEz);
  }

  if(m_iLowerXNeighbor >= 0) { // Not Lower(Left) Boundary along x-coordinates.
    buf_lower_x_Ey = new double[buf_size_x_Ey];
    buf_lower_x_Ez = new double[buf_size_x_Ez];

    // X : Ey
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index = j + (m_iGridSizeY-1)*k;
        index2 = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_x_Ey[index] = (*m_vCurrentEy)[index2];
      }
    }
    // X  : Ez
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index = j + m_iGridSizeY*k;
        index2 = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_x_Ez[index] = (*m_vCurrentEz)[index2];
      }
    }

    MPI_Send(buf_lower_x_Ey, buf_size_x_Ey, MPI_DOUBLE, m_iLowerXNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_x_Ez, buf_size_x_Ez, MPI_DOUBLE, m_iLowerXNeighbor, 1, MPI_COMM_WORLD);
  }

  
  if(m_iUpperXNeighbor >= 0) { // Not Upper(Right) Boundary along x-coordinates.
    MPI_Wait(&requestXEy, &statusXEy);
    MPI_Wait(&requestXEz, &statusXEz);
    
    
    // X : Ey
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
        index2 = j + (m_iGridSizeY-1)*k;
        index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        //cout << "index : " << index << endl;
        //cout << "index2 : " << index2 << endl;
        //EM_LOG("here? 7-2-2-1");
        //double temp;
        (*m_vCurrentEy)[index] = buf_upper_x_Ey[index2];
        //temp = buf_upper_x_Ey[index2];
        //cout << "temp : " << temp << endl;
        //EM_LOG("here? 7-2-2-2");
        //(*m_vCurrentEy)[index] = temp;
        //temp = (*m_vCurrentEy)[index];
        //EM_LOG("here? 7-2-2-3");
      }
    }
    
    // X  : Ez
    #pragma omp parallel for private(j,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(j=0 ; j<m_iGridSizeY ; j++ ) {
        index2 = j + m_iGridSizeY*k;
        index = m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vCurrentEz)[index] = buf_upper_x_Ez[index2];
      }
    }

  }




  /*
  if(timestep == 2) {
    size_t i,j,k,index;
    stringstream fnEz;
    fnEz << "middle2_updateCurrent_m_vCurrentEz_" << m_iRank << ".txt";
    FILE* fsEz = fopen(fnEz.str().c_str(), "w");
    fprintf(fsEz, "timestep : %d\n", timestep);
    
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fprintf(fsEz, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, (*m_vCurrentEz)[index]);
        }
      }
    }
    fclose(fsEz);
  }
  //*/
  
  

  MPI_Barrier(MPI_COMM_WORLD);




  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    buf_upper_y_Ex = new double[buf_size_y_Ex];
    buf_upper_y_Ez = new double[buf_size_y_Ez];
    MPI_Irecv(buf_upper_y_Ex, buf_size_y_Ex, MPI_DOUBLE, m_iUpperYNeighbor, 0, MPI_COMM_WORLD, &requestYEx);
    MPI_Irecv(buf_upper_y_Ez, buf_size_y_Ez, MPI_DOUBLE, m_iUpperYNeighbor, 1, MPI_COMM_WORLD, &requestYEz);
  }

  
  if(m_iLowerYNeighbor >= 0) { // Not Lower(Left) Boundary along y-coordinates.
    buf_lower_y_Ex = new double[buf_size_y_Ex];
    buf_lower_y_Ez = new double[buf_size_y_Ez];

    // Y : Ex
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + (m_iGridSizeX-1)*k;
        index2 = i + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_y_Ex[index] = (*m_vCurrentEx)[index2];
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*k;
        index2 = i + m_iGridSizeX*m_iGridSizeY*k;
        buf_lower_y_Ez[index] = (*m_vCurrentEz)[index2];
      }
    }

    MPI_Send(buf_lower_y_Ex, buf_size_y_Ex, MPI_DOUBLE, m_iLowerYNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_y_Ez, buf_size_y_Ez, MPI_DOUBLE, m_iLowerYNeighbor, 1, MPI_COMM_WORLD);
  }


  
  if(m_iUpperYNeighbor >= 0) { // Not Upper(Right) Boundary along y-coordinates.
    MPI_Wait(&requestYEx, &statusYEx);
    MPI_Wait(&requestYEz, &statusYEz);

    // Y : Ex
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ ; k++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index2 = i + (m_iGridSizeX-1)*k;
        index = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vCurrentEx)[index] = buf_upper_y_Ex[index2];
      }
    }
    // Y : Ez
    #pragma omp parallel for private(i,k,index,index2)
    for(k=0 ; k<m_iGridSizeZ-1 ; k++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + m_iGridSizeX*k;
        index = i + m_iGridSizeX*(m_iGridSizeY - 1) + m_iGridSizeX*m_iGridSizeY*k;
        (*m_vCurrentEz)[index] = buf_upper_y_Ez[index2];
      }
    }
  }





  /*
  if(timestep == 2) {
    size_t i,j,k,index;
    stringstream fnEz;
    fnEz << "middle3_updateCurrent_m_vCurrentEz_" << m_iRank << ".txt";
    FILE* fsEz = fopen(fnEz.str().c_str(), "w");
    fprintf(fsEz, "timestep : %d\n", timestep);
    
    for(k=0 ; k<m_iGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGridSizeX-1 ; i++) {
          index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
          fprintf(fsEz, "[%d]\t[%d][%d][%d]\t%.60f\n", index, i, j, k, (*m_vCurrentEz)[index]);
        }
      }
    }
    fclose(fsEz);
  }
  //*/

  MPI_Barrier(MPI_COMM_WORLD);




  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    buf_upper_z_Ex = new double[buf_size_z_Ex];
    buf_upper_z_Ey = new double[buf_size_z_Ey];
    MPI_Irecv(buf_upper_z_Ex, buf_size_z_Ex, MPI_DOUBLE, m_iUpperZNeighbor, 0, MPI_COMM_WORLD, &requestZEx);
    MPI_Irecv(buf_upper_z_Ey, buf_size_z_Ey, MPI_DOUBLE, m_iUpperZNeighbor, 1, MPI_COMM_WORLD, &requestZEy);
  }


  if(m_iLowerZNeighbor >= 0) { // Not Lower(Left) Boundary along z-coordinates.
    buf_lower_z_Ex = new double[buf_size_z_Ex];
    buf_lower_z_Ey = new double[buf_size_z_Ey];

    // Z : Ex
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index = i + (m_iGridSizeX-1)*j;
        index2 = i + m_iGridSizeX*j;
        buf_lower_z_Ex[index] = (*m_vCurrentEx)[index2];
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index = i + m_iGridSizeX*j;
        index2 = i + m_iGridSizeX*j;
        buf_lower_z_Ey[index] = (*m_vCurrentEy)[index2];
      }
    }

    MPI_Send(buf_lower_z_Ex, buf_size_z_Ex, MPI_DOUBLE, m_iLowerZNeighbor, 0, MPI_COMM_WORLD);
    MPI_Send(buf_lower_z_Ey, buf_size_z_Ey, MPI_DOUBLE, m_iLowerZNeighbor, 1, MPI_COMM_WORLD);
  }


  if(m_iUpperZNeighbor >= 0) { // Not Upper(Right) Boundary along z-coordinates.
    MPI_Wait(&requestZEx, &statusZEx);
    MPI_Wait(&requestZEy, &statusZEy);
    
    // Z : Ex
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY ; j++ ) {
      for(i=0 ; i<m_iGridSizeX-1 ; i++ ) {
        index2 = i + (m_iGridSizeX-1)*j;
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        (*m_vCurrentEx)[index] = buf_upper_z_Ex[index2];
      }
    }
    // Z : Ey
    #pragma omp parallel for private(i,j,index,index2)
    for(j=0 ; j<m_iGridSizeY-1 ; j++ ) {
      for(i=0 ; i<m_iGridSizeX ; i++ ) {
        index2 = i + m_iGridSizeX*j;
        index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*(m_iGridSizeZ-1);
        (*m_vCurrentEy)[index] = buf_upper_z_Ey[index2];
      }
    }
  }
  // End : Sending Ex, Ey, Ez to Upper(Right) Cell.



  //MPI_Barrier(MPI_COMM_WORLD);


  if(buf_lower_x_Ey != 0) delete[] buf_lower_x_Ey;
  if(buf_lower_x_Ez != 0) delete[] buf_lower_x_Ez;
  if(buf_upper_x_Ey != 0) delete[] buf_upper_x_Ey;
  if(buf_upper_x_Ez != 0) delete[] buf_upper_x_Ez;

  if(buf_lower_y_Ex != 0) delete[] buf_lower_y_Ex;
  if(buf_lower_y_Ez != 0) delete[] buf_lower_y_Ez;
  if(buf_upper_y_Ex != 0) delete[] buf_upper_y_Ex;
  if(buf_upper_y_Ez != 0) delete[] buf_upper_y_Ez;

  if(buf_lower_z_Ex != 0) delete[] buf_lower_z_Ex;
  if(buf_lower_z_Ey != 0) delete[] buf_lower_z_Ey;
  if(buf_upper_z_Ex != 0) delete[] buf_upper_z_Ex;
  if(buf_upper_z_Ey != 0) delete[] buf_upper_z_Ey;




  //end = omp_get_wtime();
  //cout << "Communitation Time " << end-start << " seconds." << endl;
  ////////////////////////////////////////////////////////////////////////////////
  // End : Communication of E
  ////////////////////////////////////////////////////////////////////////////////



    /*
    if(timestep == 2) {
      //stringstream fn2Ex;
      //stringstream fn2Ey;
      stringstream fn2Ez;
      //fn2Ex << "after_update_current_m_vCurrentEx_" << m_iRank << ".txt";
      //fn2Ey << "after_update_current_m_vCurrentEy_" << m_iRank << ".txt";
      fn2Ez << "after_update_current_m_vCurrentEz_" << m_iRank << ".txt";
      //ofstream fs2Ex(fn2Ex.str().c_str());
      //ofstream fs2Ey(fn2Ey.str().c_str());
      ofstream fs2Ez(fn2Ez.str().c_str());
      //fs2Ex.precision(60);
      //fs2Ey.precision(60);
      fs2Ez.precision(60);
      
      fs2Ez << "timestep : " << timestep << endl;
  
      for(k=0 ; k<m_iGridSizeZ ; k++) {
        for(j=0 ; j<m_iGridSizeY ; j++) {
          for(i=0 ; i<m_iGridSizeX ; i++) {
            index = i + m_iGridSizeX*j + m_iGridSizeX*m_iGridSizeY*k;
            //fs2Ex << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEx)[index] << endl;
            //fs2Ey << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEy)[index] << endl;
            fs2Ez << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << (*m_vCurrentEz)[index] << endl;
          }
        }
      }
      //fs2Ex.close();
      //fs2Ey.close();
      fs2Ez.close();
      
      //MPI_Finalize();
      //_exit(0);
    }
    //*/



    

}



void MPIFieldSolver::resetGlobalElectricIntensity() {
  int index;
  for(int k=0 ; k<m_iGlobalGridSizeZ ; k++) {
    for(int j=0 ; j<m_iGlobalGridSizeY ; j++) {
      for(int i=0 ; i<m_iGlobalGridSizeX ; i++) {
        index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
        m_vGlobalEx->at(index) = 0.0;
        m_vGlobalEy->at(index) = 0.0;
        m_vGlobalEz->at(index) = 0.0;
      }
    }
  }
}


int MPIFieldSolver::getElectricIntensityAt(double time, size_t timestep, double x, double y, double z, double& ex, double& ey, double& ez) {

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

  /*
  int globalGridSizeX = m_vGlobalEx->size();
  int globalGridSizeY = m_vGlobalEy->size();
  int globalGridSizeZ = m_vGlobalEz->size();
  */

  i = cellX + cellY*m_iGlobalGridSizeX + cellZ*m_iGlobalGridSizeX*m_iGlobalGridSizeY;
  int x1 = i + 1;
  int y1 = i + m_iGlobalGridSizeX;
  int z1 = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
  int x1y1 = x1 + m_iGlobalGridSizeX;
  int x1z1 = z1 + 1;
  int y1z1 = z1 + m_iGlobalGridSizeX;

  if(cellX == 0) {
    lxlylzEx = getBoundaryXLowerEx(m_fGlobalLowerX, m_fGlobalLowerY + m_dy*cellY, m_fGlobalLowerZ + m_dz*cellZ, time);
    lxuylzEx = getBoundaryXLowerEx(m_fGlobalLowerX, m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalLowerZ + m_dz*cellZ, time);
    lxlyuzEx = getBoundaryXLowerEx(m_fGlobalLowerX, m_fGlobalLowerY + m_dy*cellY, m_fGlobalLowerZ + m_dz*(cellZ+1), time);
    lxuyuzEx = getBoundaryXLowerEx(m_fGlobalLowerX, m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalLowerZ + m_dz*(cellZ+1), time);
  } else {
    lxlylzEx = 0.5*(m_vGlobalEx->at(i) + m_vGlobalEx->at(i-1));
    lxuylzEx = 0.5*(m_vGlobalEx->at(y1) + m_vGlobalEx->at(y1-1));
    lxlyuzEx = 0.5*(m_vGlobalEx->at(z1) + m_vGlobalEx->at(z1-1));
    lxuyuzEx = 0.5*(m_vGlobalEx->at(y1z1) + m_vGlobalEx->at(y1z1-1));
  }

  if(cellX == m_iGridSizeX - 1) {
    uxlylzEx = getBoundaryXUpperEx(m_fGlobalUpperX, m_fGlobalLowerY + m_dy*cellY, m_fGlobalLowerZ + m_dz*cellZ, time);
    uxuylzEx = getBoundaryXUpperEx(m_fGlobalUpperX, m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalLowerZ + m_dz*cellZ, time);
    uxlyuzEx = getBoundaryXUpperEx(m_fGlobalUpperX, m_fGlobalLowerY + m_dy*cellY, m_fGlobalLowerZ + m_dz*(cellZ+1), time);
    uxuyuzEx = getBoundaryXUpperEx(m_fGlobalUpperX, m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalLowerZ + m_dz*(cellZ+1), time);
  } else {
    uxlylzEx = 0.5*(m_vGlobalEx->at(i) + m_vGlobalEx->at(i+1));
    uxuylzEx = 0.5*(m_vGlobalEx->at(y1) + m_vGlobalEx->at(y1+1));
    uxlyuzEx = 0.5*(m_vGlobalEx->at(z1) + m_vGlobalEx->at(z1+1));
    uxuyuzEx = 0.5*(m_vGlobalEx->at(y1z1) + m_vGlobalEx->at(y1z1+1));
  }
  


  if(cellY == 0) {
    lxlylzEy = getBoundaryYLowerEy(m_fGlobalLowerX + m_dx*cellX, m_fGlobalLowerY, m_fGlobalLowerZ + m_dz*cellZ, time);
    uxlylzEy = getBoundaryYLowerEy(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalLowerY, m_fGlobalLowerZ + m_dz*cellZ, time);
    lxlyuzEy = getBoundaryYLowerEy(m_fGlobalLowerX + m_dx*cellX, m_fGlobalLowerY, m_fGlobalLowerZ + m_dz*(cellZ+1), time);
    uxlyuzEy = getBoundaryYLowerEy(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalLowerY, m_fGlobalLowerZ + m_dz*(cellZ+1), time);
  } else {
    lxlylzEy = 0.5*(m_vGlobalEy->at(i) + m_vGlobalEy->at(i-m_iGlobalGridSizeX));
    uxlylzEy = 0.5*(m_vGlobalEy->at(x1) + m_vGlobalEy->at(x1-m_iGlobalGridSizeX));
    lxlyuzEy = 0.5*(m_vGlobalEy->at(z1) + m_vGlobalEy->at(z1-m_iGlobalGridSizeX));
    uxlyuzEy = 0.5*(m_vGlobalEy->at(x1z1) + m_vGlobalEy->at(x1z1-m_iGlobalGridSizeX));
  }

  if(cellY == m_iGridSizeY - 1) {
    lxuylzEy = getBoundaryYUpperEy(m_fGlobalLowerX + m_dx*cellX, m_fGlobalUpperY, m_fGlobalLowerZ + m_dz*cellZ, time);
    uxuylzEy = getBoundaryYUpperEy(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalUpperY, m_fGlobalLowerZ + m_dz*cellZ, time);
    lxuyuzEy = getBoundaryYUpperEy(m_fGlobalLowerX + m_dx*cellX, m_fGlobalUpperY, m_fGlobalLowerZ + m_dz*(cellZ+1), time);
    uxuyuzEy = getBoundaryYUpperEy(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalUpperY, m_fGlobalLowerZ + m_dz*(cellZ+1), time);
  } else {
    lxuylzEy = 0.5*(m_vGlobalEy->at(i) + m_vGlobalEy->at(i+m_iGlobalGridSizeX));
    uxuylzEy = 0.5*(m_vGlobalEy->at(x1) + m_vGlobalEy->at(x1+m_iGlobalGridSizeX));
    lxuyuzEy = 0.5*(m_vGlobalEy->at(z1) + m_vGlobalEy->at(z1+m_iGlobalGridSizeX));
    uxuyuzEy = 0.5*(m_vGlobalEy->at(x1z1) + m_vGlobalEy->at(x1z1+m_iGlobalGridSizeX));
  }



  if(cellZ == 0) {
    lxlylzEz = getBoundaryZLowerEz(m_fGlobalLowerX + m_dx*cellX, m_fGlobalLowerY + m_dy*cellY, m_fGlobalLowerZ, time);
    uxlylzEz = getBoundaryZLowerEz(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalLowerY + m_dy*cellY, m_fGlobalLowerZ, time);
    lxuylzEz = getBoundaryZLowerEz(m_fGlobalLowerX + m_dx*cellX, m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalLowerZ, time);
    uxuylzEz = getBoundaryZLowerEz(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalLowerZ, time);
  } else {
    lxlylzEz = 0.5*(m_vGlobalEz->at(i) + m_vGlobalEz->at(i-m_iGlobalGridSizeX*m_iGlobalGridSizeY));
    uxlylzEz = 0.5*(m_vGlobalEz->at(x1) + m_vGlobalEz->at(x1-m_iGlobalGridSizeX*m_iGlobalGridSizeY));
    lxuylzEz = 0.5*(m_vGlobalEz->at(y1) + m_vGlobalEz->at(y1-m_iGlobalGridSizeX*m_iGlobalGridSizeY));
    uxuylzEz = 0.5*(m_vGlobalEz->at(x1y1) + m_vGlobalEz->at(x1y1-m_iGlobalGridSizeX*m_iGlobalGridSizeY));
  }

  if(cellZ == m_iGridSizeZ - 1) {
    lxlyuzEz = getBoundaryZUpperEz(m_fGlobalLowerX + m_dx*cellX, m_fGlobalLowerY + m_dy*cellY, m_fGlobalUpperZ, time);
    uxlyuzEz = getBoundaryZUpperEz(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalLowerY + m_dy*cellY, m_fGlobalUpperZ, time);
    lxuyuzEz = getBoundaryZUpperEz(m_fGlobalLowerX + m_dx*cellX, m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalUpperZ, time);
    uxuyuzEz = getBoundaryZUpperEz(m_fGlobalLowerX + m_dx*(cellX+1), m_fGlobalLowerY + m_dy*(cellY+1), m_fGlobalUpperZ, time);
  } else {
    lxlyuzEz = 0.5*(m_vGlobalEz->at(i) + m_vGlobalEz->at(i+m_iGlobalGridSizeX*m_iGlobalGridSizeY));
    uxlyuzEz = 0.5*(m_vGlobalEz->at(x1) + m_vGlobalEz->at(x1+m_iGlobalGridSizeX*m_iGlobalGridSizeY));
    lxuyuzEz = 0.5*(m_vGlobalEz->at(y1) + m_vGlobalEz->at(y1+m_iGlobalGridSizeX*m_iGlobalGridSizeY));
    uxuyuzEz = 0.5*(m_vGlobalEz->at(x1y1) + m_vGlobalEz->at(x1y1+m_iGlobalGridSizeX*m_iGlobalGridSizeY));
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

int MPIFieldSolver::getMagneticDensityAt(double time, size_t timestep, double x, double y, double z, double& bx, double& by, double& bz) {
  
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

  /*
  int globalGridSizeX = m_vGlobalEx->size();
  int globalGridSizeY = m_vGlobalEy->size();
  int globalGridSizeZ = m_vGlobalEz->size();
  */

  i = cellX + cellY*m_iGlobalGridSizeX + cellZ*m_iGlobalGridSizeX*m_iGlobalGridSizeX;

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
  //double mucoeff = getMu(x,y,z,time)*0.5;

  if( (x >= m_fGlobalLowerX + 0.5*m_dx) && (x <= m_fGlobalUpperX - 0.5*m_dx) && (y >= m_fGlobalLowerY + 0.5*m_dy) && (y <= m_fGlobalUpperY - 0.5*m_dy) && (z >= m_fGlobalLowerZ + 0.5*m_dz) && (z <= m_fGlobalUpperZ - 0.5*m_dz) ) {
    // Hx
    if( (dy < 0.5) && (dz < 0.5) ) {
      lxlylz = i - m_iGlobalGridSizeX - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxlyuz = i - m_iGlobalGridSizeX;
      lxuyuz = i;
      uxlylz = i - m_iGlobalGridSizeX - m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      uxuylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      uxlyuz = i - m_iGlobalGridSizeX + 1;
      uxuyuz = i + 1;
      rdx = dx;
      rdy = dy + 0.5;
      rdz = dz + 0.5;
    } else if ( (dy < 0.5) && (dz >= 0.5) ) {
      lxlylz = i - m_iGlobalGridSizeX;
      lxuylz = i;
      lxlyuz = i - m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlylz = i - m_iGlobalGridSizeX + 1;
      uxuylz = i + 1;
      uxlyuz = i - m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      uxuyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      rdx = dx;
      rdy = dy + 0.5;
      rdz = dz - 0.5;
    } else if ( (dy >= 0.5) && (dz < 0.5) ) {
      lxlylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuylz = i + m_iGlobalGridSizeX - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxlyuz = i;
      lxuyuz = i + m_iGlobalGridSizeX;
      uxlylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      uxuylz = i + m_iGlobalGridSizeX - m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      uxlyuz = i + 1;
      uxuyuz = i + m_iGlobalGridSizeX + 1;
      rdx = dx;
      rdy = dy - 0.5;
      rdz = dz + 0.5;
    } else {
      lxlylz = i;
      lxuylz = i + m_iGlobalGridSizeX;
      lxlyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuyuz = i + m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlylz = i + 1;
      uxuylz = i + m_iGlobalGridSizeX + 1;
      uxlyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      uxuyuz = i + m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY + 1;
      rdx = dx;
      rdy = dy - 0.5;
      rdz = dz - 0.5;
    }
    
    lxlylzBx = m_vGlobalBx->at(lxlylz);
    uxlylzBx = m_vGlobalBx->at(uxlylz);
    lxuylzBx = m_vGlobalBx->at(lxuylz);
    uxuylzBx = m_vGlobalBx->at(uxuylz);
    lxlyuzBx = m_vGlobalBx->at(lxlyuz);
    uxlyuzBx = m_vGlobalBx->at(uxlyuz);
    lxuyuzBx = m_vGlobalBx->at(lxuyuz);
    uxuyuzBx = m_vGlobalBx->at(uxuyuz);
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
      lxlylz = i - 1 - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxlyuz = i - 1;
      uxlyuz = i;
      lxuylz = i - 1 - m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      uxuylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      lxuyuz = i - 1 + m_iGlobalGridSizeX;
      uxuyuz = i + m_iGlobalGridSizeX;
      rdx = dx + 0.5;
      rdy = dy;
      rdz = dz + 0.5;
    } else if ( (dx < 0.5) && (dz >= 0.5) ) {
      lxlylz = i - 1;
      uxlylz = i;
      lxlyuz = i - 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuylz = i - 1 + m_iGlobalGridSizeX;
      uxuylz = i + m_iGlobalGridSizeX;
      lxuyuz = i - 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      uxuyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      rdx = dx + 0.5;
      rdy = dy;
      rdz = dz - 0.5;
    } else if ( (dx >= 0.5) && (dz < 0.5) ) {
      lxlylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlylz = i + 1 - m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxlyuz = i;
      uxlyuz = i + 1;
      lxuylz = i - m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      uxuylz = i + 1 - m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      lxuyuz = i + m_iGlobalGridSizeX;
      uxuyuz = i + 1 + m_iGlobalGridSizeX;
      rdx = dx - 0.5;
      rdy = dy;
      rdz = dz + 0.5;
    } else {
      lxlylz = i;
      uxlylz = i + 1;
      lxlyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlyuz = i + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuylz = i + m_iGlobalGridSizeX;
      uxuylz = i + 1 + m_iGlobalGridSizeX;
      lxuyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      uxuyuz = i + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY + m_iGlobalGridSizeX;
      rdx = dx - 0.5;
      rdy = dy;
      rdz = dz - 0.5;
    }
    
    lxlylzBy = m_vGlobalBy->at(lxlylz);
    uxlylzBy = m_vGlobalBy->at(uxlylz);
    lxuylzBy = m_vGlobalBy->at(lxuylz);
    uxuylzBy = m_vGlobalBy->at(uxuylz);
    lxlyuzBy = m_vGlobalBy->at(lxlyuz);
    uxlyuzBy = m_vGlobalBy->at(uxlyuz);
    lxuyuzBy = m_vGlobalBy->at(lxuyuz);
    uxuyuzBy = m_vGlobalBy->at(uxuyuz);
    by = lxlylzBy*(1-rdx)*(1-rdy)*(1-rdz) + uxlylzBy*rdx*(1-rdy)*(1-rdz) + lxuylzBy*(1-rdx)*rdy*(1-rdz) + uxuylzBy*rdx*rdy*(1-rdz) + lxlyuzBy*(1-rdx)*(1-rdy)*rdz + uxlyuzBy*rdx*(1-rdy)*rdz + lxuyuzBy*(1-rdx)*rdy*rdz + uxuyuzBy*rdx*rdy*rdz;

    // Hz
    if( (dx < 0.5) && (dy < 0.5) ) {
      lxlylz = i - m_iGlobalGridSizeX - 1;
      uxlylz = i - m_iGlobalGridSizeX;
      lxuylz = i - 1;
      uxuylz = i;
      lxlyuz = i - m_iGlobalGridSizeX - 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlyuz = i - m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuyuz = i - 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxuyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      rdx = dx + 0.5;
      rdy = dy + 0.5;
      rdz = dz;
    } else if ( (dx < 0.5) && (dy >= 0.5) ) {
      lxlylz = i - 1;
      uxlylz = i;
      lxuylz = i - 1 + m_iGlobalGridSizeX;
      uxuylz = i + m_iGlobalGridSizeX;
      lxlyuz = i - 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuyuz = i - 1 + m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxuyuz = i + m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      rdx = dx + 0.5;
      rdy = dy - 0.5;
      rdz = dz;
    } else if ( (dx >= 0.5) && (dy < 0.5) ) {
      lxlylz = i - m_iGlobalGridSizeX;
      uxlylz = i - m_iGlobalGridSizeX + 1;
      lxuylz = i;
      uxuylz = i + 1;
      lxlyuz = i - m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlyuz = i - m_iGlobalGridSizeX + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxuyuz = i + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      rdx = dx - 0.5;
      rdy = dy + 0.5;
      rdz = dz;
    } else {
      lxlylz = i;
      uxlylz = i + 1;
      lxuylz = i + m_iGlobalGridSizeX;
      uxuylz = i + m_iGlobalGridSizeX + 1;
      lxlyuz = i + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxlyuz = i + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      lxuyuz = i + m_iGlobalGridSizeX + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      uxuyuz = i + m_iGlobalGridSizeX + 1 + m_iGlobalGridSizeX*m_iGlobalGridSizeY;
      rdx = dx - 0.5;
      rdy = dy - 0.5;
      rdz = dz;
    }
    
    lxlylzBz = m_vGlobalBz->at(lxlylz);
    uxlylzBz = m_vGlobalBz->at(uxlylz);
    lxuylzBz = m_vGlobalBz->at(lxuylz);
    uxuylzBz = m_vGlobalBz->at(uxuylz);
    lxlyuzBz = m_vGlobalBz->at(lxlyuz);
    uxlyuzBz = m_vGlobalBz->at(uxlyuz);
    lxuyuzBz = m_vGlobalBz->at(lxuyuz);
    uxuyuzBz = m_vGlobalBz->at(uxuyuz);
    bz = lxlylzBz*(1-rdx)*(1-rdy)*(1-rdz) + uxlylzBz*rdx*(1-rdy)*(1-rdz) + lxuylzBz*(1-rdx)*rdy*(1-rdz) + uxuylzBz*rdx*rdy*(1-rdz) + lxlyuzBz*(1-rdx)*(1-rdy)*rdz + uxlyuzBz*rdx*(1-rdy)*rdz + lxuyuzBz*(1-rdx)*rdy*rdz + uxuyuzBz*rdx*rdy*rdz;
  } else { // in the case that poisiton is located at boundary side.
    // TODO : implementation.
    return -1;
  }



  /*
  if(timestep == 2) {
    
    stringstream fnBx;
    stringstream fnBy;
    stringstream fnBz;
    stringstream fnEx;
    stringstream fnEy;
    stringstream fnEz;
    fnBx << "global_Bx_" << m_iRank << ".txt";
    fnBy << "global_By_" << m_iRank << ".txt";
    fnBz << "global_Bz_" << m_iRank << ".txt";
    fnEx << "global_Ex_" << m_iRank << ".txt";
    fnEy << "global_Ey_" << m_iRank << ".txt";
    fnEz << "global_Ez_" << m_iRank << ".txt";
    ofstream fsBx(fnBx.str().c_str());
    ofstream fsBy(fnBy.str().c_str());
    ofstream fsBz(fnBz.str().c_str());
    ofstream fsEx(fnEx.str().c_str());
    ofstream fsEy(fnEy.str().c_str());
    ofstream fsEz(fnEz.str().c_str());
    
    fsBx.precision(60);
    fsBy.precision(60);
    fsBz.precision(60);
    fsEx.precision(60);
    fsEy.precision(60);
    fsEz.precision(60);


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
    fsBx << "timestep : " << timestep << ", mu=" << m_fMu << endl;
    fsBx << "timestep : " << timestep << ", m_vGlobalBx->at(8065)=" << m_vGlobalBx->at(8065) << endl;
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

    int i,j,k,index;
    for(k=0 ; k<m_iGlobalGridSizeZ-1 ; k++) {
      for(j=0 ; j<m_iGlobalGridSizeY-1 ; j++) {
        for(i=0 ; i<m_iGlobalGridSizeX-1 ; i++) {
          index = i + m_iGlobalGridSizeX*j + m_iGlobalGridSizeX*m_iGlobalGridSizeY*k;
          fsBx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBx->at(index) << endl;
          fsBy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBy->at(index) << endl;
          fsBz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalBz->at(index) << endl;
          fsEx << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalEx->at(index) << endl;
          fsEy << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalEy->at(index) << endl;
          fsEz << "["<<index<<"]\t["<<i<<"]["<<j<<"]["<<k<<"]\t" << m_vGlobalEz->at(index) << endl;
        }
      }
    }
    fsBx.close();
    fsBy.close();
    fsBz.close();
    fsEx.close();
    fsEy.close();
    fsEz.close();
    //_exit(0);
  }
  //*/

  return EM_SUCCESS;
}



int MPIFieldSolver::findCellHavingPosition(double x, double y, double z, int& cellX, int& cellY, int& cellZ, double& dx, double& dy, double& dz) {
  if( (x < m_fGlobalLowerX) || (x > m_fGlobalUpperX) || (y < m_fGlobalLowerY) || (y > m_fGlobalUpperY) || (z < m_fGlobalLowerZ) || (z > m_fGlobalUpperZ) ) {
    return EM_ERR_COM_DOM_SCOPE;
  }

  double distanceX = x - m_fGlobalLowerX;
  double distanceY = y - m_fGlobalLowerY;
  double distanceZ = z - m_fGlobalLowerZ;

  cellX = static_cast<int>(floor(distanceX/m_dx));
  cellY = static_cast<int>(floor(distanceY/m_dy));
  cellZ = static_cast<int>(floor(distanceZ/m_dz));

  dx = distanceX/m_dx - cellX; // (distanceX - cellX*m_dx)/m_dx
  dy = distanceY/m_dy - cellY; // (distanceY - cellX*m_dy)/m_dy
  dz = distanceZ/m_dz - cellZ; // (distanceZ - cellX*m_dz)/m_dz

  if(cellX == m_vGlobalEx->size()) {
    cellX--;
    dx = 1.0;
  }
  if(cellY == m_vGlobalEy->size()) {
    cellY--;
    dy = 1.0;
  }
  if(cellZ == m_vGlobalEz->size()) {
    cellZ--;
    dz = 1.0;
  }

  return EM_SUCCESS;
}


////////////////////////////////////////////////////////////////////////////////
// End : MPIFieldSolver
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// Start : MPITestSolver
////////////////////////////////////////////////////////////////////////////////
//*
const double MPITestSolver::mu = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
const double MPITestSolver::epsilon = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
const double MPITestSolver::frequency = 2.45e9;
const double MPITestSolver::omega = 2*PI*MPITestSolver::frequency;
//*/

MPITestSolver::MPITestSolver(int rank, int procNumX, int procNumY, int procNumZ, double dt, double leftX, double rightX, int gridNumX, double leftY, double rightY, int gridNumY, double leftZ, double rightZ, int gridNumZ)
: MPIFieldSolver(rank, procNumX, procNumY, procNumZ, dt, leftX, rightX, gridNumX, leftY, rightY, gridNumY, leftZ, rightZ, gridNumZ)
{
}



void MPITestSolver::initializeSolver() {

  double a = m_fGlobalUpperX - m_fGlobalLowerX;
  double b = m_fGlobalUpperY - m_fGlobalLowerY;

  m_fPIa = PI/a;
  m_fPIb = PI/b;

  m_fBeta = sqrt(omega*omega*mu*epsilon - m_fPIa*m_fPIa - m_fPIb*m_fPIb);
  m_fHsquare = m_fPIa*m_fPIa + m_fPIb*m_fPIb;
  m_fH = sqrt(m_fHsquare);

  MPIFieldSolver::initializeSolver();
}





double MPITestSolver::getBoundaryXLowerEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*m_fLowerX) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryXLowerEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*m_fLowerX) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryXLowerEz(double x, double y, double z, double t) {
  return sin(m_fPIa*m_fLowerX) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryXUpperEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*m_fUpperX) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryXUpperEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*m_fUpperX) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryXUpperEz(double x, double y, double z, double t) {
  return sin(m_fPIa*m_fUpperX) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryYLowerEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*m_fLowerY) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryYLowerEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*m_fLowerY) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryYLowerEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*m_fLowerY) * cos(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryYUpperEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*m_fUpperY) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryYUpperEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*m_fUpperY) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryYUpperEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*m_fUpperY) * cos(omega*t - m_fBeta*z);
}

double MPITestSolver::getBoundaryZLowerEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*m_fLowerZ);
}

double MPITestSolver::getBoundaryZLowerEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*m_fLowerZ);
}

double MPITestSolver::getBoundaryZLowerEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*m_fLowerZ);
}

double MPITestSolver::getBoundaryZUpperEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*m_fUpperZ);
}

double MPITestSolver::getBoundaryZUpperEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*m_fUpperZ);
}

double MPITestSolver::getBoundaryZUpperEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*m_fUpperZ);
}






double MPITestSolver::getInitialEx(double x, double y, double z) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(-m_fBeta*z);
}

double MPITestSolver::getInitialEy(double x, double y, double z) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(-m_fBeta*z);
}

double MPITestSolver::getInitialEz(double x, double y, double z) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(-m_fBeta*z);
}

double MPITestSolver::getInitialHx(double x, double y, double z) {
  //return -(omega*epsilon/m_fHsquare)*m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*0.5*m_dt - m_fBeta*z);
  return -(omega*epsilon/m_fHsquare)*m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(- m_fBeta*z);
}

double MPITestSolver::getInitialHy(double x, double y, double z) {
  //return (omega*epsilon/m_fHsquare)*m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*0.5*m_dt - m_fBeta*z);
  return (omega*epsilon/m_fHsquare)*m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(- m_fBeta*z);
}

double MPITestSolver::getInitialHz(double x, double y, double z) {return 0.0;}





double MPITestSolver::getExactEx(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getExactEy(double x, double y, double z, double t) {
  return m_fBeta/m_fHsquare * m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getExactEz(double x, double y, double z, double t) {
  return sin(m_fPIa*x) * sin(m_fPIb*y) * cos(omega*t - m_fBeta*z);
}

double MPITestSolver::getExactHx(double x, double y, double z, double t) {
  return -(omega*epsilon/m_fHsquare)*m_fPIb * sin(m_fPIa*x) * cos(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getExactHy(double x, double y, double z, double t) {
  return (omega*epsilon/m_fHsquare)*m_fPIa * cos(m_fPIa*x) * sin(m_fPIb*y) * sin(omega*t - m_fBeta*z);
}

double MPITestSolver::getExactHz(double x, double y, double z, double t) {return 0.0;}



void MPITestSolver::evaluateError(double time, double& Ex, double& Ey, double& Ez, double& Hx, double& Hy, double& Hz, double threshold) {


  

  size_t index;

#ifndef _MSC_VER
  size_t i,j,k;
#else
  int i,j,k;
#endif 

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
// End : MPITestSolver
////////////////////////////////////////////////////////////////////////////////





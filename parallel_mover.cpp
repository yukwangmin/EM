/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Tue Aug. 07 2012
 * 
 *      
 */





#include <list>
#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <mpi.h>


#include "parallel_mover.h"
#include "solver.h"
#include "em_error.h"

using namespace std;


int MPIParticleMover::addParticle(Particle& p, int group) {
  map< int, list<Particle> >::iterator it = m_mParticlesGroup.find(group);
  if(it == m_mParticlesGroup.end()) {
    EM_ERROR(STR_ERR_PARTICLE_WRONG_PARTICLEGROUP);
    return EM_ERR_PARTICLE_WRONG_PARTICLEGROUP;
  } else {
    it->second.push_back(p);
  }
  return EM_SUCCESS;
}

size_t MPIParticleMover::getNumberOfParticles(int group) {
  map< int, list<Particle> >::iterator it = m_mParticlesGroup.find(group);
  if(it == m_mParticlesGroup.end()) {
    return 0;
  } else {
    return it->second.size();
  }
}

size_t MPIParticleMover::getTotalNumberOfParticles() {
  size_t totalNum = 0;
  for(map< int, list<Particle> >::iterator it2 = m_mParticlesGroup.begin() ; it2 != m_mParticlesGroup.end() ; ++it2) {
    totalNum += it2->second.size();
  }
  return totalNum;
}

list<Particle>* MPIParticleMover::getParticles(int group) {
  map< int, list<Particle> >::iterator it = m_mParticlesGroup.find(group);
  if(it == m_mParticlesGroup.end()) {
    return 0;
  } else {
    return &(it->second);
  }
}

int MPIParticleMover::moveParticles(double time, size_t timestep) {
  
  int err;

  m_pFieldSolver->resetGlobalElectricIntensity();

  for(map< int, list<Particle> >::iterator it2 = m_mParticlesGroup.begin() ; it2 != m_mParticlesGroup.end() ; ++it2) {
    list<Particle>& particles = it2->second;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      double x,y,z;
      //double dx,dy,dz; // relative distance ration in a cell. 0 <= value < 1
      //int cellX, cellY, cellZ;
  
      double new_x, new_y, new_z;
      //double new_dx, new_dy, new_dz; // relative distance ration in a cell. 0 <= value < 1
      //int new_cellX, new_cellY, new_cellZ;
  
      //double ex,ey,ez,bx,by,bz;
  
      
  
      (*iter).getPosition(x,y,z);
      (*iter).updatePosition(m_dt);
      (*iter).getPosition(new_x,new_y,new_z);
  
      //*
      err = updateCurrentDensity(time, (*iter).getCharge(), x, y, z, new_x, new_y, new_z);
      if(err) {
        EM_ERROR("Updating Current Density Error");
        return err;
      }
      //*/
    }
  }
  
  m_pFieldSolver->updateGlobalElectricIntensityByCurrent(timestep);

 
  m_pFieldSolver->getGlobalFieldData(timestep);

  for(map< int, list<Particle> >::iterator it2 = m_mParticlesGroup.begin() ; it2 != m_mParticlesGroup.end() ; ++it2) {
    list<Particle>& particles = it2->second;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      double ex,ey,ez,bx,by,bz;
      double x,y,z;
  
      (*iter).getPosition(x,y,z);
      
      
      err = m_pFieldSolver->getElectricIntensityAt(time,timestep,x,y,z,ex,ey,ez);
      if(err) {
        EM_ERROR("Getting Electric Intensity at a Particle Position Error");
        return err;
      }
      err = m_pFieldSolver->getMagneticDensityAt(time,timestep,x,y,z,bx,by,bz);
      if(err) {
        EM_ERROR("Getting Magnetic Density at a Particle Position Error");
        return err;
      }
      /*
      {
        printf("timestep : %d, x : %.50f\n", timestep, x);
        printf("timestep : %d, y : %.50f\n", timestep, y);
        printf("timestep : %d, z : %.50f\n", timestep, z);
        //printf("timestep : %d, new_x : %.50f\n", timestep, new_x);
        //printf("timestep : %d, new_y : %.50f\n", timestep, new_y);
        //printf("timestep : %d, new_z : %.50f\n", timestep, new_z);
        printf("timestep : %d, ex : %.50f\n", timestep, ex);
        printf("timestep : %d, ey : %.50f\n", timestep, ey);
        printf("timestep : %d, ez : %.50f\n", timestep, ez);
        printf("timestep : %d, bx : %.50f\n", timestep, bx);
        printf("timestep : %d, by : %.50f\n", timestep, by);
        printf("timestep : %d, bz : %.50f\n", timestep, bz);
        if(timestep == 5) _exit(0);
      }
      //*/
  
  
      (*iter).updateElectricIntensity(ex,ey,ez);
      (*iter).updateMagneticDensity(bx,by,bz);
  
      (*iter).updateVelocity(m_dt);
      
      
    }
  }

  


  return 0;
}


void MPIParticleMover::calculateInitialVelocity() {
  double x,y,z;
  double ex, ey, ez;
  double hx, hy, hz;
  double bx, by, bz;
  double mu = m_pFieldSolver->getMu(0,0,0,0);
  
  for(map< int, list<Particle> >::iterator it2 = m_mParticlesGroup.begin() ; it2 != m_mParticlesGroup.end() ; ++it2) {
    list<Particle>& particles = it2->second;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      Particle& p = (*iter);
  
      p.getPosition(x,y,z);
      ex = m_pFieldSolver->getInitialEx(x,y,z);
      ey = m_pFieldSolver->getInitialEy(x,y,z);
      ez = m_pFieldSolver->getInitialEz(x,y,z);
      hx = m_pFieldSolver->getInitialHx(x,y,z);
      hy = m_pFieldSolver->getInitialHy(x,y,z);
      hz = m_pFieldSolver->getInitialHz(x,y,z);
      bx = mu*hx;
      by = mu*hy;
      bz = mu*hz;
      
      p.updateElectricIntensity(ex,ey,ez);
      p.updateMagneticDensity(bx,by,bz);
      p.updateInitialVelocity(m_dt);
    }
  }
}



int MPIParticleMover::updateCurrentDensity(double time, double charge, double x, double y, double z, double new_x, double new_y, double new_z) {
  double dx,dy,dz; // relative distance ration in a cell. 0 <= value < 1
  int cellX, cellY, cellZ;

  double ndx, ndy, ndz; // relative distance ration in a cell. 0 <= value < 1
  int ncellX, ncellY, ncellZ;

  // directional vectors
  double dirX = new_x - x;
  double dirY = new_y - y;
  double dirZ = new_z - z;

  // dx in a grid cell in FieldSolver discretization.
  double gDX = m_pFieldSolver->getdx();
  double gDY = m_pFieldSolver->getdy();
  double gDZ = m_pFieldSolver->getdz();

  double lowerX = m_pFieldSolver->getGlobalLowerX();
  double lowerY = m_pFieldSolver->getGlobalLowerY();
  double lowerZ = m_pFieldSolver->getGlobalLowerZ();
  
  // Since x,y,z were new_x,new_y,new_z in the last time step, we don't need to chech their domain.
  m_pFieldSolver->findCellHavingPosition(x,y,z,cellX,cellY,cellZ,dx,dy,dz);
  int error = m_pFieldSolver->findCellHavingPosition(new_x,new_y,new_z,ncellX,ncellY,ncellZ,ndx,ndy,ndz);
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

  int diffX = ncellX-cellX;
  int diffY = ncellY-cellY;
  int diffZ = ncellZ-cellZ;

  if((abs(diffX) > 1) || (abs(diffY) > 1) || (abs(diffZ) > 1)) {
    EM_ERROR(STR_ERR_MAT_CFL_GEN);
    return EM_ERR_MAT_CFL_GEN;
  }

  if((diffX == -1) && (dx == 0.0)) {
    cellX--;
    diffX = 0;
  }
  if((diffX == 1) && (ndx == 0.0)) {
    ncellX--;
    diffX = 0;
  }
  if((diffY == -1) && (dy == 0.0)) {
    cellY--;
    diffY = 0;
  }
  if((diffY == 1) && (ndy == 0.0)) {
    ncellY--;
    diffY = 0;
  }
  if((diffZ == -1) && (dz == 0.0)) {
    cellZ--;
    diffZ = 0;
  }
  if((diffZ == 1) && (ndz == 0.0)) {
    ncellZ--;
    diffZ = 0;
  }

  int passCnt = 0;
  if(diffX != 0) passCnt++;
  if(diffY != 0) passCnt++;
  if(diffZ != 0) passCnt++;


  //double tmp_x, tmp_y, tmp_z;

  if( passCnt == 0 ) { // 1 cell
    updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, new_x, new_y, new_z, ndx, ndy, ndz);



  } else if ( passCnt == 1 ) { // 2 cell
    double tmpT; // temporary time.
    double tmpX, tmpY, tmpZ;
    double tmpdx, tmpdy, tmpdz;

    if(diffX != 0) { // 2 cell : x-direction
      if(diffX < 0) {
        tmpX = cellX*gDX + lowerX;
        tmpdx = 0.0;
      } else {
        tmpX = (cellX+1)*gDX + lowerX;
        tmpdx = 1.0;
      }
      tmpT = (tmpX - x)/dirX;
      tmpY = y + dirY*tmpT;
      tmpZ = z + dirZ*tmpT;
      tmpdy = (tmpY - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;
      tmpdz = (tmpZ - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;
      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX, tmpY, tmpZ, tmpdx, tmpdy, tmpdz);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX, tmpY, tmpZ, 1.0 - tmpdx, tmpdy, tmpdz, new_x, new_y, new_z, ndx, ndy, ndz);
    } else if(diffY != 0) { // 2 cell : y-direction
      if(diffY < 0) {
        tmpY = cellY*gDY + lowerY;
        tmpdy = 0.0;
      } else {
        tmpY = (cellY+1)*gDY + lowerY;
        tmpdy = 1.0;
      }
      tmpT = (tmpY - y)/dirY;
      tmpX = x + dirX*tmpT;
      tmpZ = z + dirZ*tmpT;
      tmpdx = (tmpX - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
      tmpdz = (tmpZ - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;
      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX, tmpY, tmpZ, tmpdx, tmpdy, tmpdz);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX, tmpY, tmpZ, tmpdx, 1.0 - tmpdy, tmpdz, new_x, new_y, new_z, ndx, ndy, ndz);
    } else { // 2 cell : z-direction
      if(diffZ < 0) {
        tmpZ = cellZ*gDZ + lowerZ;
        tmpdz = 0.0;
      } else {
        tmpZ = (cellZ+1)*gDZ + lowerZ;
        tmpdz = 1.0;
      }
      tmpT = (tmpZ - z)/dirZ;
      tmpX = x + dirX*tmpT;
      tmpY = y + dirY*tmpT;
      tmpdx = (tmpX - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
      tmpdy = (tmpY - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;
      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX, tmpY, tmpZ, tmpdx, tmpdy, tmpdz);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX, tmpY, tmpZ, tmpdx, tmpdy, 1.0 - tmpdz, new_x, new_y, new_z, ndx, ndy, ndz);
    }



  } else if ( passCnt == 2 ) {
    int tmpCellX = cellX;
    int tmpCellY = cellY;
    int tmpCellZ = cellZ;

    double tmpT1, tmpT2;
    double tmpX1, tmpY1, tmpZ1, tmpX2, tmpY2, tmpZ2;
    double tmpdx1, tmpdy1, tmpdz1, tmpdx2, tmpdy2, tmpdz2;

    if(diffX == 0) { // 3 cell : y,z-direction
      tmpY1 = (diffY < 0) ? (cellY*gDY + lowerY) : ((cellY+1)*gDY + lowerY);
      tmpZ1 = (diffZ < 0) ? (cellZ*gDZ + lowerZ) : ((cellZ+1)*gDZ + lowerZ);

      tmpT1 = (tmpY1 - y)/dirY;
      tmpT2 = (tmpZ1 - z)/dirZ;

      if(tmpT1 >= tmpT2) {
        double tmp = tmpT1;
        tmpT1 = tmpT2;
        tmpT2 = tmp;
        if(diffZ < 0) {
          tmpdz1 = 0.0;
          tmpCellZ--;
        } else {
          tmpdz1 = 1.0;
          tmpCellZ++;
        }
        tmpX1 = x + dirX*tmpT1;
        tmpY1 = y + dirY*tmpT1;
        tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
        tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;

        updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);

        if(diffY < 0) {
          tmpY2 = cellY*gDY + lowerY;
          tmpdy2 = 0.0;
        } else {
          tmpY2 = (cellY+1)*gDY + lowerY;
          tmpdy2 = 1.0;
        }

        tmpX2 = x + dirX*tmpT2;
        tmpZ2 = z + dirZ*tmpT2;
        tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX; // (tmpX - (tmpCellX*gDX + lowerX))/gDX;
        tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ; // (tmpZ - (tmpCellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, tmpCellX, tmpCellY, tmpCellZ, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, 1.0 - tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);
        updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX2, tmpY2, tmpZ2, tmpdx2, 1.0 - tmpdy2, tmpdz2, new_x, new_y, new_z, ndx, ndy, ndz);

      } else { // tmpT1 < tmpT2
        if(diffY < 0) {
          tmpdy1 = 0.0;
          tmpCellY--;
        } else {
          tmpdy1 = 1.0;
          tmpCellY++;
        }
        tmpX1 = x + dirX*tmpT1;
        tmpZ1 = z + dirZ*tmpT1;
        tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
        tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);

        if(diffZ < 0) {
          tmpZ2 = cellZ*gDZ + lowerZ;
          tmpdz2 = 0.0;
        } else {
          tmpZ2 = (cellZ+1)*gDZ + lowerZ;
          tmpdz2 = 1.0;
        }

        tmpX2 = x + dirX*tmpT2;
        tmpY2 = y + dirY*tmpT2;
        tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX; // (tmpX - (tmpCellX*gDX + lowerX))/gDX;
        tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY; // (tmpY - (tmpCellY*gDY + lowerY))/gDY;

        updateCurrentDensityInOneCell(time, charge, tmpCellX, tmpCellY, tmpCellZ, tmpX1, tmpY1, tmpZ1, tmpdx1, 1.0 - tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);
        updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, 1.0 - tmpdz2, new_x, new_y, new_z, ndx, ndy, ndz);
      }

    } else if(diffY == 0) { // 3 cell : x,z-direction
      tmpX1 = (diffX < 0) ? (cellX*gDX + lowerX) : ((cellX+1)*gDX + lowerX);
      tmpZ1 = (diffZ < 0) ? (cellZ*gDZ + lowerZ) : ((cellZ+1)*gDZ + lowerZ);

      tmpT1 = (tmpX1 - x)/dirX;
      tmpT2 = (tmpZ1 - z)/dirZ;

      if(tmpT1 >= tmpT2) {
        double tmp = tmpT1;
        tmpT1 = tmpT2;
        tmpT2 = tmp;
        if(diffZ < 0) {
          tmpdz1 = 0.0;
          tmpCellZ--;
        } else {
          tmpdz1 = 1.0;
          tmpCellZ++;
        }
        tmpX1 = x + dirX*tmpT1;
        tmpY1 = y + dirY*tmpT1;
        tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
        tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;

        updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);

        if(diffX < 0) {
          tmpX2 = cellX*gDX + lowerX;
          tmpdx2 = 0.0;
        } else {
          tmpX2 = (cellX+1)*gDX + lowerX;
          tmpdx2 = 1.0;
        }

        tmpY2 = y + dirY*tmpT2;
        tmpZ2 = z + dirZ*tmpT2;
        tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY; // (tmpY - (tmpCellY*gDY + lowerY))/gDY;
        tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ; // (tmpZ - (tmpCellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, tmpCellX, tmpCellY, tmpCellZ, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, 1.0 - tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);
        updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX2, tmpY2, tmpZ2, 1.0 - tmpdx2, tmpdy2, tmpdz2, new_x, new_y, new_z, ndx, ndy, ndz);

      } else { // tmpT1 < tmpT2
        if(diffX < 0) {
          tmpdx1 = 0.0;
          tmpCellX--;
        } else {
          tmpdx1 = 1.0;
          tmpCellX++;
        }
        tmpY1 = y + dirY*tmpT1;
        tmpZ1 = z + dirZ*tmpT1;
        tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;
        tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);

        if(diffZ < 0) {
          tmpZ2 = cellZ*gDZ + lowerZ;
          tmpdz2 = 0.0;
        } else {
          tmpZ2 = (cellZ+1)*gDZ + lowerZ;
          tmpdz2 = 1.0;
        }

        tmpX2 = x + dirX*tmpT2;
        tmpY2 = y + dirY*tmpT2;
        tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX; // (tmpX - (tmpCellX*gDX + lowerX))/gDX;
        tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY; // (tmpY - (tmpCellY*gDY + lowerY))/gDY;

        updateCurrentDensityInOneCell(time, charge, tmpCellX, tmpCellY, tmpCellZ, tmpX1, tmpY1, tmpZ1, 1.0 - tmpdx1, tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);
        updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, 1.0 - tmpdz2, new_x, new_y, new_z, ndx, ndy, ndz);
      }

    } else { // diffZ == 0 // 3 cell : x,y-direction 
      tmpX1 = (diffX < 0) ? (cellX*gDX + lowerX) : ((cellX+1)*gDX + lowerX);
      tmpY1 = (diffY < 0) ? (cellY*gDY + lowerY) : ((cellY+1)*gDY + lowerY);

      tmpT1 = (tmpX1 - x)/dirX;
      tmpT2 = (tmpY1 - y)/dirY;

      if(tmpT1 >= tmpT2) {
        double tmp = tmpT1;
        tmpT1 = tmpT2;
        tmpT2 = tmp;
        if(diffY < 0) {
          tmpdy1 = 0.0;
          tmpCellY--;
        } else {
          tmpdy1 = 1.0;
          tmpCellY++;
        }
        tmpX1 = x + dirX*tmpT1;
        tmpZ1 = z + dirZ*tmpT1;
        tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
        tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);

        if(diffX < 0) {
          tmpX2 = cellX*gDX + lowerX;
          tmpdx2 = 0.0;
        } else {
          tmpX2 = (cellX+1)*gDX + lowerX;
          tmpdx2 = 1.0;
        }

        tmpY2 = y + dirY*tmpT2;
        tmpZ2 = z + dirZ*tmpT2;
        tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY; // (tmpY - (tmpCellY*gDY + lowerY))/gDY;
        tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ; // (tmpZ - (tmpCellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, tmpCellX, tmpCellY, tmpCellZ, tmpX1, tmpY1, tmpZ1, tmpdx1, 1.0 - tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);
        updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX2, tmpY2, tmpZ2, 1.0 - tmpdx2, tmpdy2, tmpdz2, new_x, new_y, new_z, ndx, ndy, ndz);

      } else { // tmpT1 < tmpT2
        if(diffX < 0) {
          tmpdx1 = 0.0;
          tmpCellX--;
        } else {
          tmpdx1 = 1.0;
          tmpCellX++;
        }
        tmpY1 = y + dirY*tmpT1;
        tmpZ1 = z + dirZ*tmpT1;
        tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;
        tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);

        if(diffY < 0) {
          tmpY2 = cellY*gDY + lowerY;
          tmpdy2 = 0.0;
        } else {
          tmpY2 = (cellY+1)*gDY + lowerY;
          tmpdy2 = 1.0;
        }

        tmpX2 = x + dirX*tmpT2;
        tmpZ2 = z + dirZ*tmpT2;
        tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX; // (tmpX - (tmpCellX*gDX + lowerX))/gDX;
        tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ; // (tmpZ - (tmpCellZ*gDZ + lowerZ))/gDZ;

        updateCurrentDensityInOneCell(time, charge, tmpCellX, tmpCellY, tmpCellZ, tmpX1, tmpY1, tmpZ1, 1.0 - tmpdx1, tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);
        updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX2, tmpY2, tmpZ2, tmpdx2, 1.0 - tmpdy2, tmpdz2, new_x, new_y, new_z, ndx, ndy, ndz);
      }
      
    }



  } else { // if ( (diffX != 0) && (diffY != 0) && (diffZ != 0) ) // 4 cell : all direction
    int tmpCellX1 = cellX;
    int tmpCellY1 = cellY;
    int tmpCellZ1 = cellZ;

    double tmpT1, tmpT2, tmpT3;
    double tmpX1, tmpY1, tmpZ1, tmpX2, tmpY2, tmpZ2, tmpX3, tmpY3, tmpZ3;
    double tmpdx1, tmpdy1, tmpdz1, tmpdx2, tmpdy2, tmpdz2, tmpdx3, tmpdy3, tmpdz3;

    tmpX1 = (diffX < 0) ? (cellX*gDX + lowerX) : ((cellX+1)*gDX + lowerX);
    tmpY1 = (diffY < 0) ? (cellY*gDY + lowerY) : ((cellY+1)*gDY + lowerY);
    tmpZ1 = (diffZ < 0) ? (cellZ*gDZ + lowerZ) : ((cellZ+1)*gDZ + lowerZ);

    tmpT1 = (tmpX1 - x)/dirX;
    tmpT2 = (tmpY1 - y)/dirY;
    tmpT3 = (tmpZ1 - z)/dirZ;

    if((tmpT1 <= tmpT2) && (tmpT2 <= tmpT3)) {
      if(diffX < 0) {
        tmpdx1 = 0.0;
        tmpCellX1--;
      } else {
        tmpdx1 = 1.0;
        tmpCellX1++;
      }
      tmpY1 = y + dirY*tmpT1;
      tmpZ1 = z + dirZ*tmpT1;
      tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;
      tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);
      
      int tmpCellX2 = tmpCellX1;
      int tmpCellY2 = tmpCellY1;
      int tmpCellZ2 = tmpCellZ1;

      if(diffY < 0) {
        tmpY2 = cellY*gDY + lowerY;
        tmpdy2 = 0.0;
        tmpCellY2--;
      } else {
        tmpY2 = (cellY+1)*gDY + lowerY;
        tmpdy2 = 1.0;
        tmpCellY2++;
      }

      tmpX2 = x + dirX*tmpT2;
      tmpZ2 = z + dirZ*tmpT2;
      tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX1; // (tmpX - (tmpCellX1*gDX + lowerX))/gDX;
      tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ1; // (tmpZ - (tmpCellZ1*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, tmpCellX1, tmpCellY1, tmpCellZ1, tmpX1, tmpY1, tmpZ1, 1.0 - tmpdx1, tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);

      if(diffZ < 0) {
        tmpZ3 = cellZ*gDZ + lowerZ;
        tmpdz3 = 0.0;
      } else {
        tmpZ3 = (cellZ+1)*gDZ + lowerZ;
        tmpdz3 = 1.0;
      }

      tmpX3 = x + dirX*tmpT3;
      tmpY3 = y + dirY*tmpT3;
      tmpdx3 = (tmpX3 - lowerX)/gDX - tmpCellX2; // (tmpX - (tmpCellX2*gDX + lowerX))/gDX;
      tmpdy3 = (tmpY3 - lowerY)/gDY - tmpCellY2; // (tmpY - (tmpCellY2*gDY + lowerY))/gDY;

      updateCurrentDensityInOneCell(time, charge, tmpCellX2, tmpCellY2, tmpCellZ2, tmpX2, tmpY2, tmpZ2, tmpdx2, 1.0 - tmpdy2, tmpdz2, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, tmpdz3);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, 1.0 - tmpdz3, new_x, new_y, new_z, ndx, ndy, ndz);

    } else if((tmpT1 <= tmpT3) && (tmpT3 <= tmpT2)) {
      if(diffX < 0) {
        tmpdx1 = 0.0;
        tmpCellX1--;
      } else {
        tmpdx1 = 1.0;
        tmpCellX1++;
      }
      tmpY1 = y + dirY*tmpT1;
      tmpZ1 = z + dirZ*tmpT1;
      tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;
      tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);
      
      int tmpCellX2 = tmpCellX1;
      int tmpCellY2 = tmpCellY1;
      int tmpCellZ2 = tmpCellZ1;

      if(diffZ < 0) {
        tmpZ2 = cellZ*gDZ + lowerZ;
        tmpdz2 = 0.0;
        tmpCellZ2--;
      } else {
        tmpZ2 = (cellZ+1)*gDZ + lowerZ;
        tmpdz2 = 1.0;
        tmpCellZ2++;
      }

      tmpX2 = x + dirX*tmpT3;
      tmpY2 = y + dirY*tmpT3;
      tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX1; // (tmpX - (tmpCellX1*gDX + lowerX))/gDX;
      tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY1; // (tmpY - (tmpCellY1*gDY + lowerY))/gDY;

      updateCurrentDensityInOneCell(time, charge, tmpCellX1, tmpCellY1, tmpCellZ1, tmpX1, tmpY1, tmpZ1, 1.0 - tmpdx1, tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);

      if(diffY < 0) {
        tmpY3 = cellY*gDY + lowerY;
        tmpdy3 = 0.0;
      } else {
        tmpY3 = (cellY+1)*gDY + lowerY;
        tmpdy3 = 1.0;
      }

      tmpX3 = x + dirX*tmpT2;
      tmpZ3 = z + dirZ*tmpT2;
      tmpdx3 = (tmpX3 - lowerX)/gDX - tmpCellX2; // (tmpX - (tmpCellX2*gDX + lowerX))/gDX;
      tmpdz3 = (tmpZ3 - lowerZ)/gDZ - tmpCellZ2; // (tmpZ - (tmpCellZ2*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, tmpCellX2, tmpCellY2, tmpCellZ2, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, 1.0 - tmpdz2, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, tmpdz3);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX3, tmpY3, tmpZ3, tmpdx3, 1.0 - tmpdy3, tmpdz3, new_x, new_y, new_z, ndx, ndy, ndz);

    } else if((tmpT2 <= tmpT1) && (tmpT1 <= tmpT3)) {
      if(diffY < 0) {
        tmpdy1 = 0.0;
        tmpCellY1--;
      } else {
        tmpdy1 = 1.0;
        tmpCellY1++;
      }
      tmpX1 = x + dirX*tmpT2;
      tmpZ1 = z + dirZ*tmpT2;
      tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
      tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);
      
      int tmpCellX2 = tmpCellX1;
      int tmpCellY2 = tmpCellY1;
      int tmpCellZ2 = tmpCellZ1;

      if(diffX < 0) {
        tmpX2 = cellX*gDX + lowerX;
        tmpdx2 = 0.0;
        tmpCellX2--;
      } else {
        tmpX2 = (cellX+1)*gDX + lowerX;
        tmpdx2 = 1.0;
        tmpCellX2++;
      }

      tmpY2 = y + dirY*tmpT1;
      tmpZ2 = z + dirZ*tmpT1;
      tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY1; // (tmpY - (tmpCellY1*gDY + lowerY))/gDY;
      tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ1; // (tmpZ - (tmpCellZ1*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, tmpCellX1, tmpCellY1, tmpCellZ1, tmpX1, tmpY1, tmpZ1, tmpdx1, 1.0 - tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);

      if(diffZ < 0) {
        tmpZ3 = cellZ*gDZ + lowerZ;
        tmpdz3 = 0.0;
      } else {
        tmpZ3 = (cellZ+1)*gDZ + lowerZ;
        tmpdz3 = 1.0;
      }

      tmpX3 = x + dirX*tmpT3;
      tmpY3 = y + dirY*tmpT3;
      tmpdx3 = (tmpX3 - lowerX)/gDX - tmpCellX2; // (tmpX - (tmpCellX2*gDX + lowerX))/gDX;
      tmpdy3 = (tmpY3 - lowerY)/gDY - tmpCellY2; // (tmpY - (tmpCellY2*gDY + lowerY))/gDY;

      updateCurrentDensityInOneCell(time, charge, tmpCellX2, tmpCellY2, tmpCellZ2, tmpX2, tmpY2, tmpZ2, 1.0 - tmpdx2, tmpdy2, tmpdz2, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, tmpdz3);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, 1.0 - tmpdz3, new_x, new_y, new_z, ndx, ndy, ndz);

    } else if((tmpT2 <= tmpT3) && (tmpT3 <= tmpT1)) {
      if(diffY < 0) {
        tmpdy1 = 0.0;
        tmpCellY1--;
      } else {
        tmpdy1 = 1.0;
        tmpCellY1++;
      }
      tmpX1 = x + dirX*tmpT2;
      tmpZ1 = z + dirZ*tmpT2;
      tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
      tmpdz1 = (tmpZ1 - lowerZ)/gDZ - cellZ; // (tmpZ - (cellZ*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);
      
      int tmpCellX2 = tmpCellX1;
      int tmpCellY2 = tmpCellY1;
      int tmpCellZ2 = tmpCellZ1;

      if(diffZ < 0) {
        tmpZ2 = cellZ*gDZ + lowerZ;
        tmpdz2 = 0.0;
        tmpCellZ2--;
      } else {
        tmpZ2 = (cellZ+1)*gDZ + lowerZ;
        tmpdz2 = 1.0;
        tmpCellZ2++;
      }

      tmpX2 = x + dirX*tmpT3;
      tmpY2 = y + dirY*tmpT3;
      tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX1; // (tmpX - (tmpCellX1*gDX + lowerX))/gDX;
      tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY1; // (tmpY - (tmpCellY1*gDY + lowerY))/gDY;
      
      updateCurrentDensityInOneCell(time, charge, tmpCellX1, tmpCellY1, tmpCellZ1, tmpX1, tmpY1, tmpZ1, tmpdx1, 1.0 - tmpdy1, tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);

      if(diffX < 0) {
        tmpX3 = cellX*gDX + lowerX;
        tmpdx3 = 0.0;
      } else {
        tmpX3 = (cellX+1)*gDX + lowerX;
        tmpdx3 = 1.0;
      }
      
      tmpY3 = y + dirY*tmpT1;
      tmpZ3 = z + dirZ*tmpT1;
      tmpdy3 = (tmpY3 - lowerY)/gDY - tmpCellY2; // (tmpY - (tmpCellY2*gDY + lowerY))/gDY;
      tmpdz3 = (tmpZ3 - lowerZ)/gDZ - tmpCellZ2; // (tmpZ - (tmpCellZ2*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, tmpCellX2, tmpCellY2, tmpCellZ2, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, 1.0 - tmpdz2, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, tmpdz3);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX3, tmpY3, tmpZ3, 1.0 - tmpdx3, tmpdy3, tmpdz3, new_x, new_y, new_z, ndx, ndy, ndz);

    } else if((tmpT3 <= tmpT1) && (tmpT1 <= tmpT2)) {
      if(diffZ < 0) {
        tmpdz1 = 0.0;
        tmpCellZ1--;
      } else {
        tmpdz1 = 1.0;
        tmpCellZ1++;
      }
      tmpX1 = x + dirX*tmpT3;
      tmpY1 = y + dirY*tmpT3;
      tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
      tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;

      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);
      
      int tmpCellX2 = tmpCellX1;
      int tmpCellY2 = tmpCellY1;
      int tmpCellZ2 = tmpCellZ1;

      if(diffX < 0) {
        tmpX2 = cellX*gDX + lowerX;
        tmpdx2 = 0.0;
        tmpCellX2--;
      } else {
        tmpX2 = (cellX+1)*gDX + lowerX;
        tmpdx2 = 1.0;
        tmpCellX2++;
      }
      
      tmpY2 = y + dirY*tmpT1;
      tmpZ2 = z + dirZ*tmpT1;
      tmpdy2 = (tmpY2 - lowerY)/gDY - tmpCellY1; // (tmpY - (tmpCellY1*gDY + lowerY))/gDY;
      tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ1; // (tmpZ - (tmpCellZ1*gDZ + lowerZ))/gDZ;
      
      updateCurrentDensityInOneCell(time, charge, tmpCellX1, tmpCellY1, tmpCellZ1, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, 1.0 - tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);

      if(diffY < 0) {
        tmpY3 = cellY*gDY + lowerY;
        tmpdy3 = 0.0;
      } else {
        tmpY3 = (cellY+1)*gDY + lowerY;
        tmpdy3 = 1.0;
      }
      
      tmpX3 = x + dirX*tmpT2;
      tmpZ3 = z + dirZ*tmpT2;
      tmpdx3 = (tmpX3 - lowerX)/gDX - tmpCellX2; // (tmpX - (tmpCellX2*gDX + lowerX))/gDX;
      tmpdz3 = (tmpZ3 - lowerZ)/gDZ - tmpCellZ2; // (tmpZ - (tmpCellZ2*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, tmpCellX2, tmpCellY2, tmpCellZ2, tmpX2, tmpY2, tmpZ2, 1.0 - tmpdx2, tmpdy2, tmpdz2, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, tmpdz3);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX3, tmpY3, tmpZ3, tmpdx3, 1.0 - tmpdy3, tmpdz3, new_x, new_y, new_z, ndx, ndy, ndz);

    } else { // (tmpT3 <= tmpT2) && (tmpT2 <= tmpT1)
      if(diffZ < 0) {
        tmpdz1 = 0.0;
        tmpCellZ1--;
      } else {
        tmpdz1 = 1.0;
        tmpCellZ1++;
      }
      tmpX1 = x + dirX*tmpT3;
      tmpY1 = y + dirY*tmpT3;
      tmpdx1 = (tmpX1 - lowerX)/gDX - cellX; // (tmpX - (cellX*gDX + lowerX))/gDX;
      tmpdy1 = (tmpY1 - lowerY)/gDY - cellY; // (tmpY - (cellY*gDY + lowerY))/gDY;

      updateCurrentDensityInOneCell(time, charge, cellX, cellY, cellZ, x, y, z, dx, dy, dz, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, tmpdz1);
      
      int tmpCellX2 = tmpCellX1;
      int tmpCellY2 = tmpCellY1;
      int tmpCellZ2 = tmpCellZ1;

      if(diffY < 0) {
        tmpY2 = cellY*gDY + lowerY;
        tmpdy2 = 0.0;
        tmpCellY2--;
      } else {
        tmpY2 = (cellY+1)*gDY + lowerY;
        tmpdy2 = 1.0;
        tmpCellY2++;
      }
      
      tmpX2 = x + dirX*tmpT2;
      tmpZ2 = z + dirZ*tmpT2;
      tmpdx2 = (tmpX2 - lowerX)/gDX - tmpCellX1; // (tmpX - (tmpCellX1*gDX + lowerX))/gDX;
      tmpdz2 = (tmpZ2 - lowerZ)/gDZ - tmpCellZ1; // (tmpZ - (tmpCellZ1*gDZ + lowerZ))/gDZ;
      
      updateCurrentDensityInOneCell(time, charge, tmpCellX1, tmpCellY1, tmpCellZ1, tmpX1, tmpY1, tmpZ1, tmpdx1, tmpdy1, 1.0 - tmpdz1, tmpX2, tmpY2, tmpZ2, tmpdx2, tmpdy2, tmpdz2);

      if(diffX < 0) {
        tmpX3 = cellX*gDX + lowerX;
        tmpdx3 = 0.0;
      } else {
        tmpX3 = (cellX+1)*gDX + lowerX;
        tmpdx3 = 1.0;
      }
      
      tmpY3 = y + dirY*tmpT1;
      tmpZ3 = z + dirZ*tmpT1;
      tmpdy3 = (tmpY3 - lowerY)/gDY - tmpCellY2; // (tmpY - (tmpCellY2*gDY + lowerY))/gDY;
      tmpdz3 = (tmpZ3 - lowerZ)/gDZ - tmpCellZ2; // (tmpZ - (tmpCellZ2*gDZ + lowerZ))/gDZ;

      updateCurrentDensityInOneCell(time, charge, tmpCellX2, tmpCellY2, tmpCellZ2, tmpX2, tmpY2, tmpZ2, tmpdx2, 1.0 - tmpdy2, tmpdz2, tmpX3, tmpY3, tmpZ3, tmpdx3, tmpdy3, tmpdz3);
      updateCurrentDensityInOneCell(time, charge, ncellX, ncellY, ncellZ, tmpX3, tmpY3, tmpZ3, 1.0 - tmpdx3, tmpdy3, tmpdz3, new_x, new_y, new_z, ndx, ndy, ndz);

    }


  }


  return 0;

}


void MPIParticleMover::updateCurrentDensityInOneCell(double time, double charge, int cellX, int cellY, int cellZ, double x, double y, double z, double dx, double dy, double dz, double new_x, double new_y, double new_z, double ndx, double ndy, double ndz) {
  if( !((dx >= 0.0) && (dx <= 1.0) && (dy >= 0.0) && (dy <= 1.0) && (dz >= 0.0) && (dz <= 1.0) && (ndx >= 0.0) && (ndx <= 1.0) && (ndy >= 0.0) && (ndy <= 1.0) && (ndz >= 0.0) && (ndz <= 1.0)) ) {
    EM_DEBUG("Wrong dx,dy,dz,ndx,ndy or ndz. At least one of them is not in the ragne between 0.0 and 1.0!");
  }
  double delX = ndx - dx;
  double delY = ndy - dy;
  double delZ = ndz - dz;

  double mdx = 0.5*(ndx + dx);
  double mdy = 0.5*(ndy + dy);
  double mdz = 0.5*(ndz + dz);

  double appX = delY*delZ/12;
  double appY = delX*delZ/12;
  double appZ = delX*delY/12;
  
  double volume = m_pFieldSolver->getdx() * m_pFieldSolver->getdy() * m_pFieldSolver->getdz();
  double chargedensity = charge/volume;

  // J
  double jx_y0z0 = chargedensity*(new_x-x)*((1 - mdy)*(1 - mdz) + appX);
  double jx_y1z0 = chargedensity*(new_x-x)*(mdy*(1 - mdz) - appX);
  double jx_y0z1 = chargedensity*(new_x-x)*((1 - mdy)*mdz - appX);
  double jx_y1z1 = chargedensity*(new_x-x)*(mdz*mdz + appX);

  double jy_x0z0 = chargedensity*(new_y-y)*((1 - mdx)*(1 - mdz) + appY);
  double jy_x1z0 = chargedensity*(new_y-y)*(mdx*(1 - mdz) - appY);
  double jy_x0z1 = chargedensity*(new_y-y)*((1 - mdx)*mdz - appY);
  double jy_x1z1 = chargedensity*(new_y-y)*(mdx*mdz + appY);

  double jz_x0y0 = chargedensity*(new_z-z)*((1 - mdx)*(1 - mdy) + appZ);
  double jz_x1y0 = chargedensity*(new_z-z)*(mdx*(1 - mdy) - appZ);
  double jz_x0y1 = chargedensity*(new_z-z)*((1 - mdx)*mdy - appZ);
  double jz_x1y1 = chargedensity*(new_z-z)*(mdx*mdy + appZ);

  /*
  {
    printf("Current from (%f, %f, %.20f) to (%f, %f, %.20f)\n", x, y, z, new_x, new_y, new_z);

    printf("jx_y0z0 = %.50f\n", jx_y0z0);
    printf("jx_y1z0 = %.50f\n", jx_y1z0);
    printf("jx_y0z0 = %.50f\n", jx_y0z0);
    printf("jx_y0z1 = %.50f\n", jx_y0z1);

    printf("jy_x0z0 = %.50f\n", jy_x0z0);
    printf("jy_x1z0 = %.50f\n", jy_x1z0);
    printf("jy_x0z1 = %.50f\n", jy_x0z1);
    printf("jy_x1z1 = %.50f\n", jy_x1z1);

    printf("jz_x0y0 = %.50f\n", jz_x0y0);
    printf("jz_x1y0 = %.50f\n", jz_x1y0);
    printf("jz_x0y1 = %.50f\n", jz_x0y1);
    printf("jz_x1y1 = %.50f\n", jz_x1y1);
  }
  //*/

  m_pFieldSolver->updateCurrent(time, cellX, cellY, cellZ, jx_y0z0, jx_y1z0, jx_y0z1, jx_y1z1, jy_x0z0, jy_x1z0, jy_x0z1, jy_x1z1, jz_x0y0, jz_x1y0, jz_x0y1, jz_x1y1);
}







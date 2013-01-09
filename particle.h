/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Tue Aug. 07 2012
 * 
 *      
 */


#ifndef __PARTICLE_H__
#define __PARTICLE_H__




class Particle {
private:
  //double m_dt;

  double m_fPositionX;
  double m_fPositionY;
  double m_fPositionZ;

  double m_fVelocityU;
  double m_fVelocityV;
  double m_fVelocityW;

  //double m_fCurrentDensityX;
  //double m_fCurrentDensityY;
  //double m_fCurrentDensityZ;

  double m_fCharge;
  double m_fMass;

  double m_fEx;
  double m_fEy;
  double m_fEz;
  double m_fBx;
  double m_fBy;
  double m_fBz;


public:
  Particle(double x, double y, double z, double u, double v, double w, double charge, double mass);

public:
  double getMass() {return m_fMass;}
  double getCharge() {return m_fCharge;}
  int getPosition(double& x, double& y, double& z) {x = m_fPositionX; y = m_fPositionY; z = m_fPositionZ; return 0;}
  int getVelocity(double& u, double& v, double& w) {u = m_fVelocityU; v = m_fVelocityV; w = m_fVelocityW; return 0;}
  int updatePosition(double dt);
  int updateVelocity(double dt);
  int updateInitialVelocity(double dt);
  int updateElectricIntensity(double ex, double ey, double ez) {m_fEx = ex; m_fEy = ey; m_fEz = ez; return 0;}
  int updateMagneticDensity(double bx, double by, double bz) {m_fBx = bx; m_fBy = by; m_fBz = bz; return 0;}
  //int updateCurrentDensityByMoving();



};





#endif //__PARTICLE_H__


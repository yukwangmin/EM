/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Tue Aug. 07 2012
 * 
 *      
 */


#include <iostream>


#include "particle.h"
#include "em_error.h"


Particle::Particle(double x, double y, double z, double u, double v, double w, double charge, double mass/*, double dt*/) 
: m_fPositionX(x), m_fPositionY(y), m_fPositionZ(z), m_fVelocityU(u), m_fVelocityV(v), m_fVelocityW(w), m_fCharge(charge), m_fMass(mass)//, m_dt(dt)
{}


int Particle::updatePosition(double dt) {

  //cout << "velocity : ("<<m_fVelocityU<<","<<m_fVelocityV<<","<<m_fVelocityW<<")" << endl;
  m_fPositionX += dt*m_fVelocityU;
  m_fPositionY += dt*m_fVelocityV;
  m_fPositionZ += dt*m_fVelocityW;
  
  return EM_SUCCESS;
}


int Particle::updateVelocity(double dt) {

  double alpha = m_fCharge*dt/m_fMass;

  // u^-
  double half_u = m_fVelocityU + 0.5*alpha*m_fEx;
  double half_v = m_fVelocityV + 0.5*alpha*m_fEy;
  double half_w = m_fVelocityW + 0.5*alpha*m_fEz;

  double t_x = alpha*m_fBx;
  double t_y = alpha*m_fBy;
  double t_z = alpha*m_fBz;

  // t^2
  double t2 = t_x*t_x + t_y*t_y + t_y*t_y;

  //
  double s = 2/(1 + t2);

  //u^- + u^- X t
  double ut_x = half_u + half_v*t_z - half_w*t_y;
  double ut_y = half_v + half_w*t_x - half_u*t_z;
  double ut_z = half_w + half_u*t_y - half_v*t_x;

  // s * (u^- + u^- X t) X t
  double sutt_x = s * (ut_y*t_z - ut_z*t_y);
  double sutt_y = s * (ut_z*t_x - ut_x*t_z);
  double sutt_z = s * (ut_x*t_y - ut_y*t_x);

  m_fVelocityU += (alpha*m_fEx + sutt_x);
  m_fVelocityV += (alpha*m_fEy + sutt_y);
  m_fVelocityW += (alpha*m_fEz + sutt_z);

  //std::cout << "velocity : " << m_fVelocityW <<"\n"<< std::endl;

  return EM_SUCCESS;

}


int Particle::updateInitialVelocity(double dt) {

  double alpha = (0.5*m_fCharge*dt)/m_fMass;
  double u = m_fVelocityU;
  double v = m_fVelocityV;
  double w = m_fVelocityW;

  m_fVelocityU = u + alpha*(m_fEx + (v*m_fBz - w*m_fBy));
  m_fVelocityV = v + alpha*(m_fEy + (w*m_fBx - u*m_fBz));
  m_fVelocityW = w + alpha*(m_fEz + (u*m_fBy - v*m_fBx));


  return EM_SUCCESS;
}

/*
int Particle::updateCurrentDensityByMoving() {

  return 0;
}
*/





/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 13 2012
 * 
 *      
 */


#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "resultview.h"


using namespace std;


#define FILE_DIGIT_SIZE 10
#define PROC_DIGIT_SIZE 4
#define GROUP_DIGIT_SIZE 2



std::string getFilename(const string& name, int digit, size_t number) {
  stringstream ss;
  ss << name;

  ss << "_";

  int num_digit = 1;
  if(number > 0) {
    num_digit = (int)(floor(log10((double)number))) + 1;
  }

  for(int i=0 ; i<digit-num_digit ; i++) {
    ss << "0";
  }

  ss << number;

  return ss.str();
}



std::string getFilenameWithProc(const string& name, int digit, size_t number, int proc_digit_size, size_t procNo) {
  stringstream ss;
  ss << name;

  ss << "_proc";

  int proc_digit = 1;
  if(procNo > 0) {
    proc_digit = (int)(floor(log10((double)procNo))) + 1;
  }

  for(int i=0 ; i<proc_digit_size-proc_digit ; i++) {
    ss << "0";
  }

  ss << procNo << "_";

  int num_digit = 1;
  if(number > 0) {
    num_digit = (int)(floor(log10((double)number))) + 1;
  }

  for(int i=0 ; i<digit-num_digit ; i++) {
    ss << "0";
  }

  ss << number;

  return ss.str();
}


std::string getFilenameWithGroup(const string& name, int digit, size_t number, int group_digit_size, size_t groupNo) {
  stringstream ss;
  ss << name;

  ss << "_grp";

  int grp_digit = 1;
  if(groupNo > 0) {
    grp_digit = (int)(floor(log10((double)groupNo))) + 1;
  }

  for(int i=0 ; i<group_digit_size-grp_digit ; i++) {
    ss << "0";
  }

  ss << groupNo << "_";

  int num_digit = 1;
  if(number > 0) {
    num_digit = (int)(floor(log10((double)number))) + 1;
  }

  for(int i=0 ; i<digit-num_digit ; i++) {
    ss << "0";
  }

  ss << number;

  return ss.str();
}



std::string getFilenameWithProcGroup(const string& name, int digit, size_t number, int proc_digit_size, size_t procNo, int group_digit_size, size_t groupNo) {
  stringstream ss;
  ss << name;

  ss << "_grp";

  int grp_digit = 1;
  if(groupNo > 0) {
    grp_digit = (int)(floor(log10((double)groupNo))) + 1;
  }

  for(int i=0 ; i<group_digit_size-grp_digit ; i++) {
    ss << "0";
  }

  ss << groupNo << "_proc";

  int proc_digit = 1;
  if(procNo > 0) {
    proc_digit = (int)(floor(log10((double)procNo))) + 1;
  }

  for(int i=0 ; i<proc_digit_size-proc_digit ; i++) {
    ss << "0";
  }

  ss << procNo << "_";

  int num_digit = 1;
  if(number > 0) {
    num_digit = (int)(floor(log10((double)number))) + 1;
  }

  for(int i=0 ; i<digit-num_digit ; i++) {
    ss << "0";
  }

  ss << number;

  return ss.str();
}




////////////////////////////////////////////////////////////////////////////////
// Start : VTKFieldViewer
////////////////////////////////////////////////////////////////////////////////

VTKFieldViewer::VTKFieldViewer(FieldSolver* solver, string outputname, int precision, int mode)
: m_pSolver(solver), m_iMode(mode), m_sOutputName(outputname), m_iPrecision(precision)
{
  m_fLowerX = m_pSolver->getLowerX();
  m_fUpperX = m_pSolver->getUpperX();
  m_fLowerY = m_pSolver->getLowerY();
  m_fUpperY = m_pSolver->getUpperY();
  m_fLowerZ = m_pSolver->getLowerZ();
  m_fUpperZ = m_pSolver->getUpperZ();
}

VTKFieldViewer::VTKFieldViewer(FieldSolver* solver, string outputname, double lowerX, double upperX, double lowerY, double upperY, double lowerZ, double upperZ, int precision, int mode)
: m_pSolver(solver), m_iMode(mode), m_sOutputName(outputname), m_iPrecision(precision)
{
  m_fLowerX = (m_pSolver->getLowerX() >= lowerX) ? m_pSolver->getLowerX() : lowerX;
  m_fUpperX = (m_pSolver->getUpperX() <= upperX) ? m_pSolver->getUpperX() : upperX;
  m_fLowerY = (m_pSolver->getLowerY() >= lowerY) ? m_pSolver->getLowerY() : lowerY;
  m_fUpperY = (m_pSolver->getUpperY() <= upperY) ? m_pSolver->getUpperY() : upperY;
  m_fLowerZ = (m_pSolver->getLowerZ() >= lowerZ) ? m_pSolver->getLowerZ() : lowerZ;
  m_fUpperZ = (m_pSolver->getUpperZ() <= upperZ) ? m_pSolver->getUpperZ() : upperZ;
}



int VTKFieldViewer::writeResult(double time, int timestep, int proc) {

  if( (m_fLowerX >= m_fUpperX) || (m_fLowerY >= m_fUpperY) || (m_fLowerZ >= m_fUpperZ) ) {
    return 1;
  }

  string outputFileName;
  if(proc < 0) outputFileName = getFilename(m_sOutputName, FILE_DIGIT_SIZE, timestep);
  else outputFileName = getFilenameWithProc(m_sOutputName, FILE_DIGIT_SIZE, timestep, PROC_DIGIT_SIZE, proc);
  
  //double dt = m_pSolver->getdt();
  //double time = dt*timestep;

  double dx = m_pSolver->getdx();
  double dy = m_pSolver->getdy();
  double dz = m_pSolver->getdz();

  double lowerX = m_pSolver->getLowerX();
  double lowerY = m_pSolver->getLowerY();
  double lowerZ = m_pSolver->getLowerZ();

  int startXindex = static_cast<int>(ceil((m_fLowerX - lowerX) / dx));
  int endXindex = static_cast<int>(floor((m_fUpperX - lowerX) / dx));
  int startYindex = static_cast<int>(ceil((m_fLowerY - lowerY) / dy));
  int endYindex = static_cast<int>(floor((m_fUpperY - lowerY) / dy));
  int startZindex = static_cast<int>(ceil((m_fLowerZ - lowerZ) / dz));
  int endZindex = static_cast<int>(floor((m_fUpperZ - lowerZ) / dz));

  int gridSizeX = endXindex - startXindex + 1;
  int gridSizeY = endYindex - startYindex + 1;
  int gridSizeZ = endZindex - startZindex + 1;

  if( (gridSizeX <= 1) || (gridSizeX <= 1) || (gridSizeX <= 1) ) {
    return 1;
  }

  int totalGridSizeX = m_pSolver->getGridSizeX();
  int totalGridSizeY = m_pSolver->getGridSizeY();
  //int totalGridSizeZ = m_pSolver->getGridSizeZ();

  vector<double>* ex;
  vector<double>* ey;
  vector<double>* ez;
  vector<double>* hx;
  vector<double>* hy;
  vector<double>* hz;

  m_pSolver->getCurrentEx(&ex);
  m_pSolver->getCurrentEy(&ey);
  m_pSolver->getCurrentEz(&ez);
  m_pSolver->getCurrentHx(&hx);
  m_pSolver->getCurrentHy(&hy);
  m_pSolver->getCurrentHz(&hz);
  
  int point_size = (gridSizeX-1)*(gridSizeY-1)*(gridSizeZ-1);


  switch(m_iMode) {

    case SCALAR:
      {
        string fnameEx = outputFileName + "_Ex_scalar.vtk";
        ofstream fsEx(fnameEx.c_str());
        fsEx.precision(m_iPrecision);
        string fnameEy = outputFileName + "_Ey_scalar.vtk";
        ofstream fsEy(fnameEy.c_str());
        fsEy.precision(m_iPrecision);
        string fnameEz = outputFileName + "_Ez_scalar.vtk";
        ofstream fsEz(fnameEz.c_str());
        fsEz.precision(m_iPrecision);
        string fnameHx = outputFileName + "_Hx_scalar.vtk";
        ofstream fsHx(fnameHx.c_str());
        fsHx.precision(m_iPrecision);
        string fnameHy = outputFileName + "_Hy_scalar.vtk";
        ofstream fsHy(fnameHy.c_str());
        fsHy.precision(m_iPrecision);
        string fnameHz = outputFileName + "_Hz_scalar.vtk";
        ofstream fsHz(fnameHz.c_str());
        fsHz.precision(m_iPrecision);
        
        stringstream header;
        header << "# vtk DataFile Version 3.0\n" << time << "\n" << "ASCII\n" << "DATASET STRUCTURED_POINTS\n" << "DIMENSIONS " << gridSizeX-1 << " " << gridSizeY-1 << " " << gridSizeZ-1 << "\n";
        header << "ORIGIN " << (lowerX + startXindex*dx + 0.5*dx) << " " << (lowerY + startYindex*dy + 0.5*dx) << " " << (lowerZ + startZindex*dz + 0.5*dx) << "\nSPACING " << dx << " " << dy << " " << dz << "\n";

        
        fsEx << scientific << header.str() << endl;
        fsEy << scientific << header.str() << endl;
        fsEz << scientific << header.str() << endl;
        fsHx << scientific << header.str() << endl;
        fsHy << scientific << header.str() << endl;
        fsHz << scientific << header.str() << endl;



        stringstream dataheaderEx;
        dataheaderEx << "POINT_DATA " << point_size << "\nSCALARS electric_field_intensity_x double\nLOOKUP_TABLE default";

        stringstream dataheaderEy;
        dataheaderEy << "POINT_DATA " << point_size << "\nSCALARS electric_field_intensity_y double\nLOOKUP_TABLE default";

        stringstream dataheaderEz;
        dataheaderEz << "POINT_DATA " << point_size << "\nSCALARS electric_field_intensity_z double\nLOOKUP_TABLE default";

        stringstream dataheaderHx;
        dataheaderHx << "POINT_DATA " << point_size << "\nSCALARS magnetic_field_intensity_x double\nLOOKUP_TABLE default";

        stringstream dataheaderHy;
        dataheaderHy << "POINT_DATA " << point_size << "\nSCALARS magnetic_field_intensity_y double\nLOOKUP_TABLE default";

        stringstream dataheaderHz;
        dataheaderHz << "POINT_DATA " << point_size << "\nSCALARS magnetic_field_intensity_z double\nLOOKUP_TABLE default";


        fsEx << dataheaderEx.str() << endl;
        fsEy << dataheaderEy.str() << endl;
        fsEz << dataheaderEz.str() << endl;
        fsHx << dataheaderHx.str() << endl;
        fsHy << dataheaderHy.str() << endl;
        fsHz << dataheaderHz.str() << endl;

        int index_x1;
        int index_y1;
        int index_z1;
        int index_x1y1;
        int index_y1z1;
        int index_x1z1;
        int i,j,k;
        int index;
        for(k=startZindex ; k<endZindex ; k++) {
          for(j=startYindex ; j<endYindex ; j++) {
            for(i=startXindex ; i<endXindex ; i++) {
              index = i + totalGridSizeX*j + totalGridSizeX*totalGridSizeY*k;

              index_x1 = index + 1;
              index_y1 = index + totalGridSizeX;
              index_z1 = index + totalGridSizeX*totalGridSizeY;
              index_x1y1 = index_x1 + totalGridSizeX;
              index_y1z1 = index_y1 + totalGridSizeX*totalGridSizeY;
              index_x1z1 = index_x1 + totalGridSizeX*totalGridSizeY;

              fsEx << (( ex->at(index) + ex->at(index_y1) + ex->at(index_z1) + ex->at(index_y1z1) )/4) << endl;
              fsEy << (( ey->at(index) + ey->at(index_x1) + ey->at(index_z1) + ey->at(index_x1z1) )/4) << endl;
              fsEz << (( ez->at(index) + ez->at(index_x1) + ez->at(index_y1) + ez->at(index_x1y1) )/4) << endl;

              fsHx << (( hx->at(index) + hx->at(index_x1) )/2) << endl;
              fsHy << (( hy->at(index) + hy->at(index_y1) )/2) << endl;
              fsHz << (( hz->at(index) + hz->at(index_z1) )/2) << endl;

              //(*m_pOStream) << ex->at(index) << " " << 0 << " " << 0 << endl; // Ex value.
              //(*m_pOStream) << 0 << " " << ey->at(index) << " " << 0 << endl; // Ey value.
              //(*m_pOStream) << 0 << " " << 0 << " " << ez->at(index) << endl; // Ez value.
              //(*m_pOStream) << hx->at(index) << " " << 0 << " " << 0 << endl; // Hx value.
              //(*m_pOStream) << 0 << " " << hy->at(index) << " " << 0 << endl; // Hy value.
              //(*m_pOStream) << 0 << " " << 0 << " " << hz->at(index) << endl; // Hz value.
            }
          }
        }



        fsEx.close();
        fsEy.close();
        fsEz.close();
        fsHx.close();
        fsHy.close();
        fsHz.close();

      }
    break;

    case VECTOR_SEPARATED:
      {
        string fnameE = outputFileName + "_E_vector.vtk";
        ofstream fsE(fnameE.c_str());
        fsE.precision(m_iPrecision);
        string fnameH = outputFileName + "_H_vector.vtk";
        ofstream fsH(fnameH.c_str());
        fsH.precision(m_iPrecision);
        
        stringstream header;
        header << "# vtk DataFile Version 3.0\n" << time << "\n" << "ASCII\n" << "DATASET STRUCTURED_POINTS\n" << "DIMENSIONS " << gridSizeX-1 << " " << gridSizeY-1 << " " << gridSizeZ-1 << "\n";
        header << "ORIGIN " << (lowerX + startXindex*dx + 0.5*dx) << " " << (lowerY + startYindex*dy + 0.5*dx) << " " << (lowerZ + startZindex*dz + 0.5*dx) << "\nSPACING " << dx << " " << dy << " " << dz << "\n";

        
        fsE << scientific << header.str() << endl;
        fsH << scientific << header.str() << endl;
        


        stringstream dataheaderE;
        dataheaderE << "POINT_DATA " << point_size << "\nVECTORS electric_field_intensity double";

        stringstream dataheaderH;
        dataheaderH << "POINT_DATA " << point_size << "\nVECTORS magnetic_field_intensity double";

        fsE << dataheaderE.str() << endl;
        fsH << dataheaderH.str() << endl;
        
        int index_x1;
        int index_y1;
        int index_z1;
        int index_x1y1;
        int index_y1z1;
        int index_x1z1;
        int i,j,k;
        int index;
        for(k=startZindex ; k<endZindex ; k++) {
          for(j=startYindex ; j<endYindex ; j++) {
            for(i=startXindex ; i<endXindex ; i++) {
              index = i + totalGridSizeX*j + totalGridSizeX*totalGridSizeY*k;

              index_x1 = index + 1;
              index_y1 = index + totalGridSizeX;
              index_z1 = index + totalGridSizeX*totalGridSizeY;
              index_x1y1 = index_x1 + totalGridSizeX;
              index_y1z1 = index_y1 + totalGridSizeX*totalGridSizeY;
              index_x1z1 = index_x1 + totalGridSizeX*totalGridSizeY;

              fsE << (( ex->at(index) + ex->at(index_y1) + ex->at(index_z1) + ex->at(index_y1z1) )/4) << "\t" << (( ey->at(index) + ey->at(index_x1) + ey->at(index_z1) + ey->at(index_x1z1) )/4) << "\t" << (( ez->at(index) + ez->at(index_x1) + ez->at(index_y1) + ez->at(index_x1y1) )/4) << endl;
              fsH << (( hx->at(index) + hx->at(index_x1) )/2) << "\t" << (( hy->at(index) + hy->at(index_y1) )/2) << "\t" << (( hz->at(index) + hz->at(index_z1) )/2) << endl;
            }
          }
        }


        fsE.close();
        fsH.close();
      }

    break;



    case VECTOR_MERGED:
      {
        string fname = outputFileName + "_merged_vector.vtk";
        ofstream fs(fname.c_str());
        fs.precision(m_iPrecision);

        
        stringstream header;
        header << "# vtk DataFile Version 3.0\n" << time << "\n" << "ASCII\n" << "DATASET STRUCTURED_POINTS\n" << "DIMENSIONS " << gridSizeX-1 << " " << gridSizeY-1 << " " << gridSizeZ-1 << "\n";
        header << "ORIGIN " << (lowerX + startXindex*dx + 0.5*dx) << " " << (lowerY + startYindex*dy + 0.5*dx) << " " << (lowerZ + startZindex*dz + 0.5*dx) << "\nSPACING " << dx << " " << dy << " " << dz << "\n";

        
        fs << scientific << header.str() << endl;
        


        stringstream dataheaderE;
        dataheaderE << "POINT_DATA " << point_size << "\nVECTORS electric_field_intensity double";



        fs << dataheaderE.str() << endl;
        
        
        int index_x1;
        int index_y1;
        int index_z1;
        int index_x1y1;
        int index_y1z1;
        int index_x1z1;
        int i,j,k;
        int index;
        for(k=startZindex ; k<endZindex ; k++) {
          for(j=startYindex ; j<endYindex ; j++) {
            for(i=startXindex ; i<endXindex ; i++) {
              index = i + totalGridSizeX*j + totalGridSizeX*totalGridSizeY*k;

              index_x1 = index + 1;
              index_y1 = index + totalGridSizeX;
              index_z1 = index + totalGridSizeX*totalGridSizeY;
              index_x1y1 = index_x1 + totalGridSizeX;
              index_y1z1 = index_y1 + totalGridSizeX*totalGridSizeY;
              index_x1z1 = index_x1 + totalGridSizeX*totalGridSizeY;

              fs << (( ex->at(index) + ex->at(index_y1) + ex->at(index_z1) + ex->at(index_y1z1) )/4) << "\t" << (( ey->at(index) + ey->at(index_x1) + ey->at(index_z1) + ey->at(index_x1z1) )/4) << "\t" << (( ez->at(index) + ez->at(index_x1) + ez->at(index_y1) + ez->at(index_x1y1) )/4) << endl;
              
            }
          }
        }



        stringstream dataheaderH;
        dataheaderH << "VECTORS magnetic_field_intensity double";

        fs << dataheaderH.str() << endl;


        for(k=startZindex ; k<endZindex ; k++) {
          for(j=startYindex ; j<endYindex ; j++) {
            for(i=startXindex ; i<endXindex ; i++) {
              index = i + totalGridSizeX*j + totalGridSizeX*totalGridSizeY*k;

              index_x1 = index + 1;
              index_y1 = index + totalGridSizeX;
              index_z1 = index + totalGridSizeX*totalGridSizeY;
              index_x1y1 = index_x1 + totalGridSizeX;
              index_y1z1 = index_y1 + totalGridSizeX*totalGridSizeY;
              index_x1z1 = index_x1 + totalGridSizeX*totalGridSizeY;

              
              fs << (( hx->at(index) + hx->at(index_x1) )/2) << "\t" << (( hy->at(index) + hy->at(index_y1) )/2) << "\t" << (( hz->at(index) + hz->at(index_z1) )/2) << endl;
            }
          }
        }


        fs.close();

      }

    break;


  }


  return 0;
}



////////////////////////////////////////////////////////////////////////////////
// End : VTKFieldViewer
////////////////////////////////////////////////////////////////////////////////









////////////////////////////////////////////////////////////////////////////////
// Start : VTKParticleViewer
////////////////////////////////////////////////////////////////////////////////

/*
VTKParticleViewer::VTKParticleViewer(ParticleMover* mover, int mode)
: m_pMover(mover), m_iMode(mode)
{
  m_fLowerX = m_pMover->getLowerX();
  m_fUpperX = m_pMover->getUpperX();
  m_fLowerY = m_pMover->getLowerY();
  m_fUpperY = m_pMover->getUpperY();
  m_fLowerZ = m_pMover->getLowerZ();
  m_fUpperZ = m_pMover->getUpperZ();
}
*/


VTKParticleViewer::VTKParticleViewer(ParticleMover* mover, string outputname, double lowerX, double upperX, double lowerY, double upperY, double lowerZ, double upperZ, int precision)
: m_pMover(mover), m_sOutputName(outputname), m_iPrecision(precision), m_fLowerX(lowerX), m_fUpperX(upperX), m_fLowerY(lowerY), m_fUpperY(upperY), m_fLowerZ(lowerZ), m_fUpperZ(upperZ)
{
}


int VTKParticleViewer::writeResult(double time, int timestep, int proc) {

  /*
  if( (m_fLowerX >= m_fUpperX) || (m_fLowerY >= m_fUpperY) || (m_fLowerZ >= m_fUpperZ) ) {
    return 1;
  }
  */

  map< int, list<Particle> >& particlesGroup = m_pMover->getParticlesGroup();
    

  for(map< int, list<Particle> >::iterator it2 = particlesGroup.begin() ; it2 != particlesGroup.end() ; ++it2) {
    list<Particle> particles = it2->second;
    int group = it2->first;

    string outputFileName;
    if(proc < 0) outputFileName = getFilenameWithGroup(m_sOutputName, FILE_DIGIT_SIZE, timestep, GROUP_DIGIT_SIZE, group);
    else outputFileName = getFilenameWithProcGroup(m_sOutputName, FILE_DIGIT_SIZE, timestep, PROC_DIGIT_SIZE, proc, GROUP_DIGIT_SIZE, group);

    string fname = outputFileName + ".vtk";
    ofstream fs(fname.c_str());
    fs.precision(m_iPrecision);

    size_t num_particles = particles.size();

    fs << scientific << "# vtk DataFile Version 3.0\n" << time << "\n" << "ASCII\n" << "DATASET POLYDATA\n" << "POINTS " << num_particles << " double" << endl;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      double x,y,z;
      (*iter).getPosition(x,y,z);
      fs << x << " " << y << " " << z << endl;
    }

    fs << "POINT_DATA " << num_particles << endl;


    fs << "SCALARS charge double\nLOOKUP_TABLE default" << endl;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      fs << (*iter).getCharge() << endl;
    }

    fs << "SCALARS mass double\nLOOKUP_TABLE default" << endl;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      fs << (*iter).getMass() << endl;
    }

    fs << "VECTORS velocity double" << endl;
    for(list<Particle>::iterator iter = particles.begin() ; iter!=particles.end() ; iter++) {
      double u,v,w;
      (*iter).getVelocity(u,v,w);
      fs << u << " " << v << " " << w << endl;
    }

    fs.close();

  }


  return 0;
}



////////////////////////////////////////////////////////////////////////////////
// End : VTKParticleViewer
////////////////////////////////////////////////////////////////////////////////


/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \date Mon Jul. 13 2012
 * 
 *      
 */


#ifndef __EM_ERROR_H__
#define __EM_ERROR_H__





#ifndef _MSC_VER
#include <unistd.h>
#else
#include <process.h>
#endif

#include <string>
#include <cstdio>


using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////
// Start : Error Codes
////////////////////////////////////////////////////////////////////////////////////////////


#define EM_SUCCESS                                       0x00000000

#define EM_ERR_UNKNOWN                                   -1

#define EM_ERR_SYS_GEN                                   0x00000001

#define EM_ERR_SYS_INITIAL_SETTING_GEN                   0x00000100
#define EM_ERR_SYS_INITIAL_SETTING_FILE                  0x00000110
#define EM_ERR_SYS_INITIAL_SETTING_PARTICLE              0x00000120

#define EM_ERR_SYS_OBJECT_CREATION                       0x00000200
#define EM_ERR_SYS_OBJECT_CREATION_FIELDSOLVER           0x00000201
#define EM_ERR_SYS_OBJECT_CREATION_PARTICLEMOVER         0x00000202
#define EM_ERR_SYS_OBJECT_CREATION_FIELDVIEWER           0x00000203
#define EM_ERR_SYS_OBJECT_CREATION_PARTICLEVIEWER        0x00000204
#define EM_ERR_SYS_OBJECT_CREATION_USERSPECPARTICLE      0x00000205



#define EM_ERR_MAT_GEN                                   0x00001000
#define EM_ERR_MAT_CFL_GEN                               0x00001100

#define EM_ERR_PHY_GEN                                   0x00002000

#define EM_ERR_COM_DOM_SCOPE                             0x00003000
#define EM_ERR_COM_DOM_SCOPE_X                           0x00003001
#define EM_ERR_COM_DOM_SCOPE_Y                           0x00003002
#define EM_ERR_COM_DOM_SCOPE_Z                           0x00003003

#define EM_ERR_FIELD                                     0x00004000 // reserved

#define EM_ERR_PARTICLE                                  0x00005000
#define EM_ERR_PARTICLE_WRONG_PARTICLEGROUP              0x00005001

#define EM_ERR_EXTERNAL_LIB_INIT                         0x00021000

#define EM_ERR_EXTERNAL_LIB_RUNNING                      0x00021100
#define EM_ERR_EXTERNAL_LIB_RUNNING_FILEOPEN             0x00021101

////////////////////////////////////////////////////////////////////////////////////////////


#define STR_ERR_UNKNOWN                                  "Unknown Error"

#define STR_ERR_SYS_GEN                                  "General System Error"

#define STR_ERR_SYS_INITIAL_SETTING_GEN                  "Initial Setting : General Error"
#define STR_ERR_SYS_INITIAL_SETTING_FILE                 "Initial Setting : File Error"
#define STR_ERR_SYS_INITIAL_SETTING_PARTICLE             "Initial Setting : Particle Error"

#define STR_ERR_SYS_OBJECT_CREATION                      "System Object Creation Error"
#define STR_ERR_SYS_OBJECT_CREATION_FIELDSOLVER          "System Object Creation Error : FieldSolver"
#define STR_ERR_SYS_OBJECT_CREATION_PARTICLEMOVER        "System Object Creation Error : ParticleMover"
#define STR_ERR_SYS_OBJECT_CREATION_FIELDVIEWER          "System Object Creation Error : FieldViewer"
#define STR_ERR_SYS_OBJECT_CREATION_PARTICLEVIEWER       "System Object Creation Error : ParticleViewer"
#define STR_ERR_SYS_OBJECT_CREATION_USERSPECPARTICLE     "System Object Creation Error : UserSpecificInitialParticleSetting"


#define STR_ERR_MAT_GEN                                  "General Mathematical Error"
#define STR_ERR_MAT_CFL_GEN                              "Violation of CFL Condition Error"


#define STR_ERR_PHY_GEN                                  "General Physical Error"


#define STR_ERR_COM_DOM_SCOPE                            "Out of Computational Domain Error"
#define STR_ERR_COM_DOM_SCOPE_X                          "Out of Computational Domain Error Especially along x-axis"
#define STR_ERR_COM_DOM_SCOPE_Y                          "Out of Computational Domain Error Especially along y-axis"
#define STR_ERR_COM_DOM_SCOPE_Z                          "Out of Computational Domain Error Especially along z-axis"


#define STR_ERR_PARTICLE                                 "Particle Error"
#define STR_ERR_PARTICLE_WRONG_PARTICLEGROUP             "Wrong Particle Group Error"


#define STR_ERR_EXTERNAL_LIB_INIT                        "External Library Initialization Error"


#define STR_ERR_EXTERNAL_LIB_RUNNING                     "External Library Running Error"
#define STR_ERR_EXTERNAL_LIB_RUNNING_FILEOPEN            "External Library Running Error : File Open"

////////////////////////////////////////////////////////////////////////////////////////////
// End : Error Codes
////////////////////////////////////////////////////////////////////////////////////////////












#ifndef _MSC_VER
  #define EM_LOG(CONTENTS) EM_logger(std::clog, __FILE__, __func__, __LINE__, getpid(), (CONTENTS))
#else
  #define EM_LOG(CONTENTS) EM_logger(std::clog, __FILE__, __FUNCTION__, __LINE__, _getpid(), (CONTENTS))
#endif


#ifndef _MSC_VER
  #define EM_ERROR(CONTENTS) EM_logger(std::cerr, __FILE__, __func__, __LINE__, getpid(), (CONTENTS))
#else
  #define EM_ERROR(CONTENTS) EM_logger(std::cerr, __FILE__, __FUNCTION__, __LINE__, _getpid(), (CONTENTS))
#endif




#ifndef RELEASE

  #ifndef _MSC_VER
    #define EM_DEBUG(CONTENTS) EM_logger(std::cout, __FILE__, __func__, __LINE__, getpid(), (CONTENTS))
  #else
    #define EM_DEBUG(CONTENTS) EM_logger(std::cout, __FILE__, __FUNCTION__, __LINE__, _getpid(), (CONTENTS))
  #endif


#else
  
  #define SPH_DEBUG(CONTENTS)

#endif// RELEASE


void EM_logger(const ostream &stream, const char filename[], const char functionname[], size_t linenumber, int processID, const string& contents);

//void EM_logger(const ostream &stream, const char filename[], const char functionname[], size_t linenumber, int processID, char contents[]);









#endif

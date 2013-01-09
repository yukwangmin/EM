

#include <iostream>
#include <cmath>


#include "initial_setting.h"
#include "em_error.h"


#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>



using namespace xercesc;
using namespace std;


string& trim(string& s) {
  s.erase(s.find_last_not_of(" \n\r\t")+1); 
  s.erase(0, s.find_first_not_of(" \n\r\t")); 
  /*
  see #include <boost/algorithm/string.hpp> 
  trim(fname);
  */
  return s;
}





void parseDomain(const DOMElement* element, Domain& domain) {

  XMLCh* ATTR_lower = XMLString::transcode("lower");
  XMLCh* ATTR_upper = XMLString::transcode("upper");
  XMLCh* ATTR_gridnumber = XMLString::transcode("gridnumber");

  DOMNodeList *domainChildrenNodeList = element->getChildNodes();
  const XMLSize_t domainNodeCount = domainChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<domainNodeCount ; ++l) {
    DOMNode* domainCurrentNode = domainChildrenNodeList->item(l);
    if(domainCurrentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* domainCurrentElement = dynamic_cast<DOMElement*>( domainCurrentNode );
      const char* tagName2 = XMLString::transcode(domainCurrentElement->getTagName());
      if( XMLString::equals(tagName2, "x") ) {
        const char* xmlch_lower = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_lower));
        const char* xmlch_upper = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_upper));
        const char* xmlch_gridnumber = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_gridnumber));
        domain.xLower = atof(xmlch_lower);
        domain.xUpper = atof(xmlch_upper);
        domain.xGridNumber = atoi(xmlch_gridnumber);
      } else if( XMLString::equals(tagName2, "y") ) {
        const char* xmlch_lower = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_lower));
        const char* xmlch_upper = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_upper));
        const char* xmlch_gridnumber = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_gridnumber));
        domain.yLower = atof(xmlch_lower);
        domain.yUpper = atof(xmlch_upper);
        domain.yGridNumber = atoi(xmlch_gridnumber);
      } else if( XMLString::equals(tagName2, "z") ) {
        const char* xmlch_lower = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_lower));
        const char* xmlch_upper = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_upper));
        const char* xmlch_gridnumber = XMLString::transcode(domainCurrentElement->getAttribute(ATTR_gridnumber));
        domain.zLower = atof(xmlch_lower);
        domain.zUpper = atof(xmlch_upper);
        domain.zGridNumber = atoi(xmlch_gridnumber);
      }
    }
  }


  XMLString::release(&ATTR_lower);
  XMLString::release(&ATTR_upper);
  XMLString::release(&ATTR_gridnumber);

}




void parseVisualization(const DOMElement* element, Visualization& visualization) {

  const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("use")));
  if(XMLString::equals(xmlch_value, "y")) {
    visualization.use = true;
  } else {
    visualization.use = false;
  }

  visualization.timeinterval = atof(XMLString::transcode(element->getAttribute(XMLString::transcode("timeinterval"))));
  visualization.classname = XMLString::transcode(element->getAttribute(XMLString::transcode("classname")));
  visualization.filename = XMLString::transcode(element->getAttribute(XMLString::transcode("filename")));
  visualization.precision = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("precision"))));
  if(visualization.precision < 1) visualization.precision = 6;


  DOMNodeList *visualizationChildrenNodeList = element->getChildNodes();
  const XMLSize_t visualizationNodeCount = visualizationChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<visualizationNodeCount ; ++l) {
    DOMNode* visualizationChildNode = visualizationChildrenNodeList->item(l);
    if(visualizationChildNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* visualizationChildElement = dynamic_cast<DOMElement*>( visualizationChildNode );
      const char* tagName2 = XMLString::transcode(visualizationChildElement->getTagName());
      if( XMLString::equals(tagName2, "domain") ) {
        parseDomain(visualizationChildElement, visualization.domain);
      }
    }
  }


}




void parseField(const DOMElement* element, FieldSetting& field) {


  DOMNodeList *fieldChildrenNodeList = element->getChildNodes();
  const XMLSize_t fieldNodeCount = fieldChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<fieldNodeCount ; ++l) {
    DOMNode* fieldChildNode = fieldChildrenNodeList->item(l);
    if(fieldChildNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* fieldChildElement = dynamic_cast<DOMElement*>( fieldChildNode );
      const char* tagName2 = XMLString::transcode(fieldChildElement->getTagName());
      if( XMLString::equals(tagName2, "solver") ) {
        field.solver = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("classname")));
      } else if( XMLString::equals(tagName2, "epsilon") ) {
        field.epsilon = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("value")));
      } else if( XMLString::equals(tagName2, "mu") ) {
        field.mu = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("value")));
      } else if( XMLString::equals(tagName2, "initial_distribution") ) {
        field.initial_distribution_ex = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("ex")));
        field.initial_distribution_ey = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("ey")));
        field.initial_distribution_ez = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("ez")));
        field.initial_distribution_hx = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("hx")));
        field.initial_distribution_hy = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("hy")));
        field.initial_distribution_hz = XMLString::transcode(fieldChildElement->getAttribute(XMLString::transcode("hz")));
      } else if( XMLString::equals(tagName2, "visualization") ) {
        parseVisualization(fieldChildElement, field.visualization);
      }
    }

  }


}




void parsePhysicalData(DOMElement* element, PhysicalData& physical_data) {

  XMLCh* ATTR_value = XMLString::transcode("value");


  DOMNodeList *physicaDataChildrenNodeList = element->getChildNodes();
  const XMLSize_t physicalDataNodeCount = physicaDataChildrenNodeList->getLength();
  
  for(XMLSize_t m=0 ; m<physicalDataNodeCount ; ++m) {
    DOMNode* physicalDataCurrentNode = physicaDataChildrenNodeList->item(m);
    if(physicalDataCurrentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* physicalDataCurrentElement = dynamic_cast<DOMElement*>( physicalDataCurrentNode );
      const char* tagName3 = XMLString::transcode(physicalDataCurrentElement->getTagName());

      
      if( XMLString::equals(tagName3, "charge") ) {
        string value = string(XMLString::transcode(physicalDataCurrentElement->getAttribute(ATTR_value)));
        physical_data.charge = trim(value);
      } else if( XMLString::equals(tagName3, "mass") ) {
        string value = string(XMLString::transcode(physicalDataCurrentElement->getAttribute(ATTR_value)));
        physical_data.mass = trim(value);
      } else if( XMLString::equals(tagName3, "velocity") ) {
        string value = string(XMLString::transcode(physicalDataCurrentElement->getAttribute(XMLString::transcode("u"))));
        physical_data.velocity_u = trim(value);
        value = string(XMLString::transcode(physicalDataCurrentElement->getAttribute(XMLString::transcode("v"))));
        physical_data.velocity_v = trim(value);
        value = string(XMLString::transcode(physicalDataCurrentElement->getAttribute(XMLString::transcode("w"))));
        physical_data.velocity_w = trim(value);
      }

    }
  }


  XMLString::release(&ATTR_value);

}




void parseUserSpecificDistribution(const DOMElement* element, UserSpecificSetting& userSetting) {
  const char* xmlch_value;

  xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("use")));
  if(XMLString::equals(xmlch_value, "y")) {
    userSetting.use = true;
  } else {
    userSetting.use = false;
  }

  userSetting.group = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("group"))));

  xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("mode")));
  userSetting.mode = xmlch_value;

  xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("classname")));
  userSetting.classname = xmlch_value;

  xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("dx")));
  userSetting.dx = atof(xmlch_value);

  DOMNodeList *userSpecificChildrenNodeList = element->getChildNodes();
  const XMLSize_t userSpecificNodeCount = userSpecificChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<userSpecificNodeCount ; ++l) {
    DOMNode* userSpecificChildNode = userSpecificChildrenNodeList->item(l);
    if(userSpecificChildNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* userSpecificChildElement = dynamic_cast<DOMElement*>( userSpecificChildNode );
      const char* tagName2 = XMLString::transcode(userSpecificChildElement->getTagName());
      if( XMLString::equals(tagName2, "physical_data") ) {
        parsePhysicalData(userSpecificChildElement, userSetting.physical_data);
      }
    }
  }

}


void parseDirectSetting(const DOMElement* element, list<InitialParticle>& particleList) {
  const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("use")));
  if(!XMLString::equals(xmlch_value, "y")) {
    return;
  }

  DOMNodeList *directSettingChildrenNodeList = element->getChildNodes();
  const XMLSize_t directSettingNodeCount = directSettingChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<directSettingNodeCount ; ++l) {
    DOMNode* directSettingChildNode = directSettingChildrenNodeList->item(l);
    if(directSettingChildNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* directSettingChildElement = dynamic_cast<DOMElement*>( directSettingChildNode );
      const char* tagName2 = XMLString::transcode(directSettingChildElement->getTagName());
      if( XMLString::equals(tagName2, "particle") ) {
        InitialParticle ip;
        ip.x = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("x")));
        ip.y = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("y")));
        ip.z = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("z")));
        ip.physical_data.velocity_u = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("u")));
        ip.physical_data.velocity_v = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("v")));
        ip.physical_data.velocity_w = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("w")));
        ip.physical_data.charge = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("charge")));
        ip.physical_data.mass = XMLString::transcode(directSettingChildElement->getAttribute(XMLString::transcode("mass")));
        particleList.push_back(ip);
      }
    }
  }
}





void parseBunch(DOMElement* element, Bunch& bunch) {

  const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("use")));
  if(XMLString::equals(xmlch_value, "y")) {
    bunch.use = true;
  } else {
    bunch.use = false;
  }

  bunch.mode = string(XMLString::transcode(element->getAttribute(XMLString::transcode("mode"))));

  bunch.group = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("group"))));

  DOMNodeList *regionChildrenNodeList = element->getChildNodes();
  const XMLSize_t regionNodeCount = regionChildrenNodeList->getLength();
  for(XMLSize_t m=0 ; m<regionNodeCount ; ++m) {
    DOMNode* regionCurrentNode = regionChildrenNodeList->item(m);
    if(regionCurrentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* regionCurrentElement = dynamic_cast<DOMElement*>( regionCurrentNode );
      const char* tagName3 = XMLString::transcode(regionCurrentElement->getTagName());

      if( XMLString::equals(tagName3, "domain") ) {
        parseDomain(regionCurrentElement, bunch.domain);
      } else if( XMLString::equals(tagName3, "level_function") ) {
        const char* xmlch_range = XMLString::transcode(regionCurrentElement->getAttribute(XMLString::transcode("range")));  
        const char* xmlch_equality = XMLString::transcode(regionCurrentElement->getAttribute(XMLString::transcode("equality")));  
        //cout << "xmlch_range = " << xmlch_range << endl;
        LevelFunction levelfunction;

        if(XMLString::equals(xmlch_range, "more")) {
          levelfunction.range = '>';
        } else {
          levelfunction.range = '<';
        }

        if(XMLString::equals(xmlch_equality, "y")) {
          levelfunction.equality = true;
        } else {
          levelfunction.equality = false;
        }

        DOMNodeList *level_functionChildrenNodeList = regionCurrentElement->getChildNodes();
        const XMLSize_t level_functionNodeCount = level_functionChildrenNodeList->getLength();
        for(XMLSize_t n=0 ; n<level_functionNodeCount ; ++n) {
          DOMNode* level_functionCurrentNode = level_functionChildrenNodeList->item(n);
          if(level_functionCurrentNode->getNodeType() == DOMNode::TEXT_NODE) {
            DOMText* level_functionCurrentText = dynamic_cast<DOMText*>( level_functionCurrentNode );
            const char* function = XMLString::transcode(level_functionCurrentText->getData());
            string fname = string(function);

            trim(fname);

            //printf("level function : %s\n", fname.c_str());
            levelfunction.function = fname;
          }
        }
        bunch.functions.push_back(levelfunction);
      } else if( XMLString::equals(tagName3, "physical_data") ) {
        parsePhysicalData(regionCurrentElement, bunch.physical_data);

      }
    }
  }
  //printf("The number of level functions : %d\n", region.functions.size());
                            

}




void parseNormalRandomBunch(DOMElement* element, NormalRandomBunch& normal_bunch) {

  const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("use")));
  if(XMLString::equals(xmlch_value, "y")) {
    normal_bunch.use = true;
  } else {
    normal_bunch.use = false;
  }

  normal_bunch.num_particles = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("num_particles"))));

  normal_bunch.group = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("group"))));

  DOMNodeList *bunchChildrenNodeList = element->getChildNodes();
  const XMLSize_t bunchNodeCount = bunchChildrenNodeList->getLength();
  for(XMLSize_t m=0 ; m<bunchNodeCount ; ++m) {
    DOMNode* bunchCurrentNode = bunchChildrenNodeList->item(m);
    if(bunchCurrentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* bunchCurrentElement = dynamic_cast<DOMElement*>( bunchCurrentNode );
      const char* tagName3 = XMLString::transcode(bunchCurrentElement->getTagName());

      if( XMLString::equals(tagName3, "mu") ) {
        normal_bunch.mu_x = atof(XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("x"))));
        normal_bunch.mu_y = atof(XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("y"))));
        normal_bunch.mu_z = atof(XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("z"))));
      } else if( XMLString::equals(tagName3, "sigma") ) {
        normal_bunch.sigma_x = atof(XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("x"))));
        normal_bunch.sigma_y = atof(XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("y"))));
        normal_bunch.sigma_z = atof(XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("z"))));
      } else if( XMLString::equals(tagName3, "domain") ) {
        parseDomain(bunchCurrentElement, normal_bunch.domain);
      } else if( XMLString::equals(tagName3, "physical_data") ) {
        parsePhysicalData(bunchCurrentElement, normal_bunch.physical_data);
      }
    }
  }
  //printf("The number of level functions : %d\n", region.functions.size());
                            

}





void parseUniformRandomBunch(DOMElement* element, UniformRandomBunch& uniform_bunch) {

  const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("use")));
  if(XMLString::equals(xmlch_value, "y")) {
    uniform_bunch.use = true;
  } else {
    uniform_bunch.use = false;
  }

  uniform_bunch.num_particles = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("num_particles"))));

  uniform_bunch.group = atoi(XMLString::transcode(element->getAttribute(XMLString::transcode("group"))));

  DOMNodeList *bunchChildrenNodeList = element->getChildNodes();
  const XMLSize_t bunchNodeCount = bunchChildrenNodeList->getLength();
  for(XMLSize_t m=0 ; m<bunchNodeCount ; ++m) {
    DOMNode* bunchCurrentNode = bunchChildrenNodeList->item(m);
    if(bunchCurrentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* bunchCurrentElement = dynamic_cast<DOMElement*>( bunchCurrentNode );
      const char* tagName3 = XMLString::transcode(bunchCurrentElement->getTagName());

      if( XMLString::equals(tagName3, "domain") ) {
        parseDomain(bunchCurrentElement, uniform_bunch.domain);
      } else if( XMLString::equals(tagName3, "level_function") ) {
        const char* xmlch_range = XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("range")));  
        const char* xmlch_equality = XMLString::transcode(bunchCurrentElement->getAttribute(XMLString::transcode("equality")));  
        //cout << "xmlch_range = " << xmlch_range << endl;
        LevelFunction levelfunction;

        if(XMLString::equals(xmlch_range, "more")) {
          levelfunction.range = '>';
        } else {
          levelfunction.range = '<';
        }

        if(XMLString::equals(xmlch_equality, "y")) {
          levelfunction.equality = true;
        } else {
          levelfunction.equality = false;
        }

        DOMNodeList *level_functionChildrenNodeList = bunchCurrentElement->getChildNodes();
        const XMLSize_t level_functionNodeCount = level_functionChildrenNodeList->getLength();
        for(XMLSize_t n=0 ; n<level_functionNodeCount ; ++n) {
          DOMNode* level_functionCurrentNode = level_functionChildrenNodeList->item(n);
          if(level_functionCurrentNode->getNodeType() == DOMNode::TEXT_NODE) {
            DOMText* level_functionCurrentText = dynamic_cast<DOMText*>( level_functionCurrentNode );
            const char* function = XMLString::transcode(level_functionCurrentText->getData());
            string fname = string(function);

            trim(fname);

            //printf("level function : %s\n", fname.c_str());
            levelfunction.function = fname;
          }
        }
        uniform_bunch.functions.push_back(levelfunction);
      } else if( XMLString::equals(tagName3, "physical_data") ) {
        parsePhysicalData(bunchCurrentElement, uniform_bunch.physical_data);
      }
    }
  }
  //printf("The number of level functions : %d\n", region.functions.size());
                            

}






void parseParticleInitialDistribution(const DOMElement* element, ParticleInitialDistribution& distribution) {
  distribution.userSetting.use = false;
  distribution.directSetting_group = 0;

  const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("buffer_size")));
  distribution.buffer_size = atoi(xmlch_value);

  DOMNodeList *initialDistributionChildrenNodeList = element->getChildNodes();
  const XMLSize_t initialDistributionNodeCount = initialDistributionChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<initialDistributionNodeCount ; ++l) {
    DOMNode* initialDistributionChildNode = initialDistributionChildrenNodeList->item(l);
    if(initialDistributionChildNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* initialDistributionChildElement = dynamic_cast<DOMElement*>( initialDistributionChildNode );
      const char* tagName2 = XMLString::transcode(initialDistributionChildElement->getTagName());
      if( XMLString::equals(tagName2, "direct_setting") ) {
        distribution.directSetting_group = atoi(XMLString::transcode(initialDistributionChildElement->getAttribute(XMLString::transcode("group"))));
        parseDirectSetting(initialDistributionChildElement, distribution.directSetting);
      } else if( XMLString::equals(tagName2, "bunch") ) {
        Bunch bunch;
        parseBunch(initialDistributionChildElement, bunch);
        distribution.bunchs.push_back(bunch);
      } else if( XMLString::equals(tagName2, "normal_random_bunch") ) {
        NormalRandomBunch normal_bunch;
        parseNormalRandomBunch(initialDistributionChildElement, normal_bunch);
        distribution.normal_bunchs.push_back(normal_bunch);
      } else if( XMLString::equals(tagName2, "uniform_random_bunch") ) {
        UniformRandomBunch uniform_bunch;
        parseUniformRandomBunch(initialDistributionChildElement, uniform_bunch);
        distribution.uniform_bunchs.push_back(uniform_bunch);
      } else if( XMLString::equals(tagName2, "user_specific_distribution") ) {
        parseUserSpecificDistribution(initialDistributionChildElement, distribution.userSetting);
      }
    }
  }

}



void parseParticle(const DOMElement* element, ParticleSetting& particle) {
  //const char* xmlch_value = XMLString::transcode(element->getAttribute(XMLString::transcode("gamma")));
  //particle.gamma = atof(xmlch_value);

  DOMNodeList *particleChildrenNodeList = element->getChildNodes();
  const XMLSize_t particleNodeCount = particleChildrenNodeList->getLength();
  for(XMLSize_t l=0 ; l<particleNodeCount ; ++l) {
    DOMNode* particleChildNode = particleChildrenNodeList->item(l);
    if(particleChildNode->getNodeType() == DOMNode::ELEMENT_NODE) {
      DOMElement* particleChildElement = dynamic_cast<DOMElement*>( particleChildNode );
      const char* tagName2 = XMLString::transcode(particleChildElement->getTagName());
      if( XMLString::equals(tagName2, "initial_distribution") ) {
        parseParticleInitialDistribution(particleChildElement, particle.distribution);
      } else if( XMLString::equals(tagName2, "visualization") ) {
        parseVisualization(particleChildElement, particle.visualization);
      }
    }

  }

}



int loadInitialSetting(const char* filename, InitialSetting* is) {

  //cout << "loadInitialSetting() : filename : " << filename << endl;

  try {
    XMLPlatformUtils::Initialize();
  //} catch (const XMLException& toCatch) {
  } catch (...) {
    // Do your failure processing here
    EM_ERROR("XMLPlatformUtils::Initialize() error");
    return EM_ERR_EXTERNAL_LIB_INIT;
  }


  XercesDOMParser *parser = new XercesDOMParser;

  if(!parser) {
    EM_ERROR("DOM Parser creation failed");
    return EM_ERR_EXTERNAL_LIB_RUNNING;
  } 


  parser->setDoSchema(false);
  parser->setDoNamespaces(false);
  parser->setValidationScheme(XercesDOMParser::Val_Auto);
  parser->setCreateEntityReferenceNodes(false);


  try {
    parser->parse(filename);
  } catch(...) {
    string str = "Opening initial setting file : ";
    str += filename;
    EM_ERROR(str);
    return EM_ERR_EXTERNAL_LIB_RUNNING_FILEOPEN;
  }

  
  try {
    DOMDocument *xmlDoc = parser->getDocument();

    DOMElement* elementRoot = xmlDoc->getDocumentElement();  
    if( !elementRoot ) {
      EM_ERROR("Empty XML docmuent");
      return EM_ERR_EXTERNAL_LIB_RUNNING;
    }

    const char* xmlch_value = 0;
    
    
    xmlch_value = XMLString::transcode(elementRoot->getAttribute(XMLString::transcode("dim")));  
    is->dim = atoi(xmlch_value);

    xmlch_value = XMLString::transcode(elementRoot->getAttribute(XMLString::transcode("omp_num_threads")));  
    is->omp_num_threads = atoi(xmlch_value);

    DOMNodeList *systemChildrenNodeList = elementRoot->getChildNodes();
    const XMLSize_t systemNodeCount = systemChildrenNodeList->getLength();
    //printf("node count : %d\n", systemNodeCount);

    for( XMLSize_t i = 0; i < systemNodeCount; ++i ) {  

      DOMNode* systemCurrentNode = systemChildrenNodeList->item(i);  

      if( systemCurrentNode->getNodeType() &&  // true is not NULL  
          systemCurrentNode->getNodeType() == DOMNode::ELEMENT_NODE ) // is element   
      //if( systemCurrentNode->getNodeType() ) 
      {  
        
        // Found node which is an Element. Re-cast node as element  
        DOMElement* systemCurrentElement = dynamic_cast<DOMElement*>( systemCurrentNode );  
        //*
        if(!systemCurrentElement->getTagName()) {
          EM_DEBUG(XMLString::transcode(systemCurrentElement->getTagName()));
        }
        //*/
        const char* tagName = XMLString::transcode(systemCurrentElement->getTagName());
        //cout << "Tag name : " << tagName << endl;
        
        if( XMLString::equals(tagName, "time") ) {
          xmlch_value = XMLString::transcode(systemCurrentElement->getAttribute(XMLString::transcode("start")));  
          is->time.start = atof(xmlch_value);
          xmlch_value = XMLString::transcode(systemCurrentElement->getAttribute(XMLString::transcode("end")));  
          is->time.end = atof(xmlch_value);
          xmlch_value = XMLString::transcode(systemCurrentElement->getAttribute(XMLString::transcode("dt")));  
          is->time.dt = atof(xmlch_value);
        } else if ( XMLString::equals(tagName, "domain") ) {
          parseDomain(systemCurrentElement, is->domain);
        } else if ( XMLString::equals(tagName, "controller") ) {
          xmlch_value = XMLString::transcode(systemCurrentElement->getAttribute(XMLString::transcode("classname")));
          is->controller = xmlch_value;
        } else if ( XMLString::equals(tagName, "field") ) {
          parseField(systemCurrentElement, is->field);
        } else if ( XMLString::equals(tagName, "particle") ) {
          parseParticle(systemCurrentElement, is->particle);
        } else if ( XMLString::equals(tagName, "looger") ) {
          // Not Implemented Yet
        } else {
          // TODO : Error
        }


      } else {
        /*
        string error = "Node parsing error : ";
        error += "TAG name = ";
        error += XMLString::transcode(systemCurrentNode->getNodeName());
        error += ", TYPE = ";
        //error += XMLString::transcode(systemCurrentNode->getNodeType());
        EM_ERROR(error);
        return 1;
        */
      }

    }  




  } catch(...) {
    //cout << "parsing failed" << endl;
    EM_ERROR("Opening initial setting file error or parsing the file!");
    return EM_ERR_EXTERNAL_LIB_RUNNING;
  }

  

  delete parser;



  XMLPlatformUtils::Terminate();



  return 0;





}












// return : The size of vectors.
size_t getInitialHexagonalPacking3D(double localXXmin, double localYYmin, double localZZmin, double localXXmax, double localYYmax, double localZZmax, double dz, std::vector<double>& coord_x, std::vector<double>& coord_y, std::vector<double>& coord_z) {
      
  double h_r = 0.5*dz;

  int l0,l1;
  int m0_odd, m1_odd, m0_even, m1_even, n0_odd, n1_odd, n0_even, n1_even;
  int nn0_odd, nn1_odd, nn0_even, nn1_even; 


  double range_x_min = localXXmin;
  double range_x_max = localXXmax;
  double range_y_min = localYYmin;
  double range_y_max = localYYmax;
  double range_z_min = localZZmin;
  double range_z_max = localZZmax;


			
  //layers
  if (localZZmin-range_z_min-h_r < 0.5*(2.0*sqrt(6.0)/3.0*h_r)) l0 = 0;
  else l0 = (int) ((localZZmin-range_z_min-h_r-0.4999999*(2.0*sqrt(6.0)/3.0*h_r))/(2.0*sqrt(6.0)/3.0*h_r))+1;
  
  l1 = (int) ((localZZmax-range_z_min-h_r-0.4999999*(2.0*sqrt(6.0)/3.0*h_r))/(2.0*sqrt(6.0)/3.0*h_r))+1; 

  //rows (for odd-numbered layers)
  if (localYYmin-range_y_min-h_r < 0.5*(sqrt(3.0)*h_r)) m0_odd = 0;
  else m0_odd = (int) ((localYYmin-range_y_min-h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1;

  m1_odd = (int) ((localYYmax-range_y_min-h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1; 
            
  //columns (for odd-numbered layers & odd-numbered rows)
  if (localXXmin-range_x_min-2.0*h_r < 0.5*(2.0*h_r)) n0_odd = 0;
  else n0_odd = (int) ((localXXmin-range_x_min-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;
  
  n1_odd = (int) ((localXXmax-range_x_min-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
            
  //columns (for odd-numbered layers & even-numbered rows)       
  if (localXXmin-range_x_min-h_r < 0.5*(2.0*h_r)) n0_even = 0;
  else n0_even = (int) ((localXXmin-range_x_min-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;

  n1_even = (int) ((localXXmax-range_x_min-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
                 
  //rows (for even-numbered layers)
  if (localYYmin-range_y_min-sqrt(3.0)/3.0*h_r < 0.5*(sqrt(3.0)*h_r)) m0_even = 0;
  else m0_even = (int) ((localYYmin-range_y_min-sqrt(3.0)/3.0*h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1;

  m1_even = (int) ((localYYmax-range_y_min-sqrt(3.0)/3.0*h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1; 
            
  //columns (for even-numbered layers & odd-numbered rows)       
  if (localXXmin-range_x_min-h_r < 0.5*(2.0*h_r)) nn0_odd = 0;
  else nn0_odd = (int) ((localXXmin-range_x_min-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;
  
  nn1_odd = (int) ((localXXmax-range_x_min-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 
                 
  //columns (for even-numbered layers & even-numbered rows)
  if (localXXmin-range_x_min-2.0*h_r < 0.5*(2.0*h_r)) nn0_even = 0;
  else nn0_even = (int) ((localXXmin-range_x_min-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;
  
  nn1_even = (int) ((localXXmax-range_x_min-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1; 




  //compute the location of particles 
  for (int i=l0; i<l1; i++) { 
    if ((i+1)%2 != 0) { //odd-numbered layers

      for (int j=m0_odd; j<m1_odd; j++) { 
        if ((j+1)%2 != 0) { //odd-numbered rows 
          for (int k=n0_odd; k<n1_odd; k++) { 
            coord_x.push_back(range_x_min+2.0*h_r+k*2.0*h_r);
            coord_y.push_back(range_y_min+h_r+j*sqrt(3.0)*h_r);
            coord_z.push_back(range_z_min+h_r+i*2.0*sqrt(6.0)/3.0*h_r);
          }
        } else { //even-numbered rows
          for (int k=n0_even; k<n1_even; k++) {
            coord_x.push_back(range_x_min+h_r+k*2.0*h_r);
            coord_y.push_back(range_y_min+h_r+j*sqrt(3.0)*h_r);
            coord_z.push_back(range_z_min+h_r+i*2.0*sqrt(6.0)/3.0*h_r);
          }
        }
	    }
			
    } else { //even-numbered layers
      for (int j=m0_even; j<m1_even; j++) { 
        if ((j+1)%2 != 0) { //odd-numbered rows
          for (int k=nn0_odd; k<nn1_odd; k++) { 
            coord_x.push_back(range_x_min+h_r+k*2.0*h_r);
            coord_y.push_back(range_y_min+sqrt(3.0)/3.0*h_r+j*sqrt(3.0)*h_r);
            coord_z.push_back(range_z_min+h_r+i*2.0*sqrt(6.0)/3.0*h_r);
          }
        } else { //even-numbered rows
          for (int k=nn0_even; k<nn1_even; k++) {
            coord_x.push_back(range_x_min+2.0*h_r+k*2.0*h_r);
            coord_y.push_back(range_y_min+sqrt(3.0)/3.0*h_r+j*sqrt(3.0)*h_r);
            coord_z.push_back(range_z_min+h_r+i*2.0*sqrt(6.0)/3.0*h_r);
          }
        }
	    }
    }    
  }

  return coord_x.size();

}


size_t getInitialRectangularPacking3D(double localXXmin, double localYYmin, double localZZmin, double localXXmax, double localYYmax, double localZZmax, double dz, std::vector<double>& coord_x, std::vector<double>& coord_y, std::vector<double>& coord_z) {

  double range_x_min = localXXmin;
  double range_x_max = localXXmax;
  double range_y_min = localYYmin;
  double range_y_max = localYYmax;
  double range_z_min = localZZmin;
  double range_z_max = localZZmax;
  
  double dx = dz;
  double dy = dz;


  size_t max_x_step = (size_t)((range_x_max - range_x_min)/dx + 1);
  size_t max_y_step = (size_t)((range_y_max - range_y_min)/dy + 1);
  size_t max_z_step = (size_t)((range_z_max - range_z_min)/dz + 1);

  for(size_t i=0 ; i<=max_x_step ; ++i) {
    for(size_t j=0 ; j<=max_y_step ; ++j) {
      for(size_t k=0 ; k<=max_z_step ; ++k) {
        coord_x.push_back(range_x_min + i*dx);
        coord_y.push_back(range_y_min + j*dy);
        coord_z.push_back(range_z_min + k*dz);
      } // for k
    } // for j
  } // for i

  return coord_x.size();
}



// return : The size of vectors.
size_t getInitialHexagonalPacking2D(double localXXmin, double localZZmin, double localXXmax, double localZZmax, double dx, std::vector<double>& coord_x, std::vector<double>& coord_z) {
  
  double h_r = 0.5*dx;

  int l0,l1;
	int n0_odd, n1_odd, n0_even, n1_even;


  double range_x_min = localXXmin;
  double range_x_max = localXXmax;
  double range_z_min = localZZmin;
  double range_z_max = localZZmax;


	if (localZZmin-range_z_min-h_r < 0.5*(sqrt(3.0)*h_r))
		l0 = 0;
	else
		l0 = (int) ((localZZmin-range_z_min-h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1;
	l1 = (int) ((localZZmax-range_z_min-h_r-0.4999999*(sqrt(3.0)*h_r))/(sqrt(3.0)*h_r))+1;

	if (localXXmin-range_x_min-2.0*h_r < 0.5*(2.0*h_r))
		n0_odd = 0;
	else
		n0_odd = (int) ((localXXmin-range_x_min-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;
	n1_odd = (int) ((localXXmax-range_x_min-2.0*h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;

	if (localXXmin-range_x_min-h_r < 0.5*(2.0*h_r))
		n0_even = 0;
	else
		n0_even = (int) ((localXXmin-range_x_min-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;
	n1_even = (int) ((localXXmax-range_x_min-h_r-0.4999999*(2.0*h_r))/(2.0*h_r))+1;


	//compute the location of particles
	for (int i=l0; i<l1; i++) {
    if ((i+1)%2 != 0) { //odd-numbered rows 
			for (int j=n0_odd; j<n1_odd; j++) {
				coord_x.push_back(range_x_min+2.0*h_r+j*2.0*h_r);
        coord_z.push_back(range_z_min+h_r+i*sqrt(3.0)*h_r);
			}
    } else { //even-numbered rows
			for (int j=n0_even; j<n1_even; j++) {
				coord_x.push_back(range_x_min+h_r+j*2.0*h_r);
        coord_z.push_back(range_z_min+h_r+i*sqrt(3.0)*h_r);
			}
		}
	}


  return coord_x.size();
}



size_t getInitialRectangularPacking2D(double localXXmin, double localZZmin, double localXXmax, double localZZmax, double dx, std::vector<double>& coord_x, std::vector<double>& coord_z) {
  double range_x_min = localXXmin;
  double range_x_max = localXXmax;
  double range_z_min = localZZmin;
  double range_z_max = localZZmax;
  
  double dz = dx;


  size_t max_x_step = (size_t)((range_x_max - range_x_min)/dx + 1);
  size_t max_z_step = (size_t)((range_z_max - range_z_min)/dz + 1);

  for(size_t i=0 ; i<=max_x_step ; ++i) {
    for(size_t k=0 ; k<=max_z_step ; ++k) {
      coord_x.push_back(range_x_min + i*dx);
      coord_z.push_back(range_z_min + k*dz);
    } // for k
  } // for i

  return coord_x.size();

}
















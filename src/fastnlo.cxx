/*

  proxy for fastnlov1.cxx or fastnlov2.cxx

*/

#include "appl_grid/fastnlo.h"

#include <iostream>
#include <fstream>
#include <sys/stat.h>


void fastnlo::readgrid( const std::string& filename ) {
  
  // can we open the fastnlo grid file?
  struct stat stfileinfo;
  if ( stat(filename.c_str(),&stfileinfo) )   {
    std::cerr << "fastnlo::fastnlo() cannot open fastnlo grid file " << filename << std::endl;
    return;
  }

  // read the fastnlo grid
  std::ifstream faststream( filename.c_str() );

  // skip the first number (magic number)
  std::string dummy; faststream >> dummy; 

  // the version number is the second
  long long tabversion;  faststream >> tabversion;
  
  if (tabversion < 20000) {
    std::cout << "fastnlo::fastnlo() reading fastnlo grid file " << filename
	      << ", version " << tabversion << std::endl;
    constructv1(filename);
  } else if (20000 <= tabversion && tabversion <= 23600) {
    std::cout << "fastnlo::fastnlo() reading fastnlo grid file " << filename
	      << ", version " << tabversion << std::endl;
    constructv2(filename);
  } else {
    std::cerr << "fastnlo::fastnlo() too recent version of fastnlo grid file " << filename << std::endl;
  }
  
  faststream.close();
  
}







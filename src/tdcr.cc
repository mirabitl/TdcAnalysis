#include "tdcrb.hh"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/dir.h>  
#include <sys/param.h>  
#include <stdio.h>  
#include <stdlib.h>  
#include <unistd.h>  
#include <string.h>
#include <stdint.h>
#include <fcntl.h>
#include <iostream>
#include <sstream>
#include <map>
#include <bitset>
#include "TApplication.h"
#include "TCanvas.h"
#include <dirent.h>
#include <fnmatch.h>

int main(int argc, char** argv )
{
  static tdcrb bs("/tmp");
  TApplication theApp("tapp", &argc, argv);
  bs.geometry("gifpp_geom.json");

  std::stringstream spat;
  int runask=atol(argv[1]);
  spat<<"SMM*"<<atol(argv[1])<<"*.dat";
  //spat<<"SMM*"<<argv[1]<<"*.dat";
 #define UNTEST
 #ifdef UNTEST
  struct dirent **namelist;
  int n;
  std::cout<<"Pattern "<<spat.str()<<std::endl;
  //std::string dirp="/data/srv02/RAID6/TDCRET";
  std::string dirp="/data/local/GIFPP";
  //dirp=".";
  n = scandir(dirp.c_str(), &namelist, NULL, alphasort);
  if (n < 0)
    perror("scandir");
  else {
    while (n--) {

      if (fnmatch(spat.str().c_str(), namelist[n]->d_name, 0)==0)
	{
	  printf("%s %d \n", namelist[n]->d_name,fnmatch(spat.str().c_str(), namelist[n]->d_name, 0));
	  printf("found\n");
	  std::stringstream sf;
	  sf<<dirp<<"/"<< namelist[n]->d_name;
	  bs.addRun(runask,sf.str());
	}
      free(namelist[n]);
    }
    free(namelist);
  }
#endif
  
  DCHistogramHandler rootHandler;
  //DHCalEventReader  dher;


  #ifdef STREAMOUT
  std::stringstream filename("");    
  char dateStr [64];
            
  time_t tm= time(NULL);
  strftime(dateStr,20,"SMO_%d%m%y_%H%M%S",localtime(&tm));
  filename<<"/tmp/"<<dateStr<<"_"<<runask<<".dat";
  int32_t _fdOut= ::open(filename.str().c_str(),O_CREAT| O_RDWR | O_NONBLOCK,S_IRWXU);
  if (_fdOut<0)
    {
      perror("No way to store to file :");
      //std::cout<<" No way to store to file"<<std::endl;
      exit(0);
    }
  else
    bs.setOutFileId(_fdOut);
#endif
  bs.setOutFileId(-1);
  #ifndef UNTEST
  
  bs.addRun(735986,"/tmp/SMO_220517_143427_735986.dat");
  #endif
  bs.Read();
  bs.setRun(runask);
  bs.end();
}

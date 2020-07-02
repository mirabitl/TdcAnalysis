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
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc, char **argv )
{
  static tdcrb bs("/tmp");
  DCHistogramHandler* rh=DCHistogramHandler::instance();
  lmana::Analyzer* a= new lmana::TdcAnalyzer(rh);
  bs.setAnalyzer(a);
  TApplication theApp("tapp", &argc, argv);
  std::string geom_file;
  std::string dirp="/home/acqilc/backup_aug2018";
  int32_t runask=0;
  char c;
   while ( (c = getopt(argc, argv, "g:r:d:hv")) != -1 ) {
     fprintf(stderr,"%c read\n",c);
     
        switch ( c ) {
            case 'g':
                /* Chaque fois que l'option -v est utilisée,
                 * on augmente le degré de verbosité. */
	      
	      geom_file.assign(optarg);
	      fprintf(stderr, "Geometry %s \n",optarg);
                break;
	case 'd':
                /* Chaque fois que l'option -v est utilisée,
                 * on augmente le degré de verbosité. */
	      
	      dirp.assign(optarg);
	      fprintf(stderr, "Directory %s \n",optarg);
                break;


            case 'r':
                /* Ici, on convertit la valeur de optarg
                 * en entier. Il conviendrait de gérer les
                 * cas où la conversion est impossible (si
                 * l'utilisateur a fait n'importe quoi...). */
                runask = atoi(optarg);
		fprintf(stderr, "Run %s \n",optarg);

                break;

            case 'h':
                /* On affiche l'aide et on termine. */
                fprintf(stderr, "Usage: tdcr [options] \n"
                                "Options:\n"
                                "  -v         Augmente la verbosité.\n"
                                "  -V n       Définit le niveau de verbosité à n.\n"
                                "  -h         Affiche ce message d'aide.\n"
                                "  -f FICHIER Écrit dans FICHIER plutôt que sur la sortie standard.\n"
                       );
                return 0;

            case 'v':
                /* On stocke le pointeur vers l'argument
                 * (car sa valeur pourrait être écrasée
                 * au prochain tour de boucle). */
	      fprintf(stderr, "On serait verbose Usage: tdcr [options] \n");
                break;

            case '?':
                /* getopt renvoie (par défaut) '?' en cas
                 * d'erreur, si une option non acceptée
                 * est utilisée ou si une option attendant
                 * un argument n'est pas suivi de son argument. */
                fprintf(stderr, "Ligne de commande incorrecte.\n");
                return 1;
        }
    }







  
   //bs.geometry("gifpp_geom.json");
   std::cout<<geom_file<<" " <<runask<<std::endl;
   bs.geometry(geom_file);
   //getchar();
  std::stringstream spat;
  //int runask=atol(argv[1]);
  spat<<"SMM*R"<<runask<<".dat";
  //spat<<"SMM*"<<argv[1]<<"*.dat";
 #define UNTEST
 #ifdef UNTEST
  struct dirent **namelist;
  int n;
  std::cout<<"Pattern "<<spat.str()<<std::endl;
  //std::string dirp="/data/srv02/RAID6/Dome0718";

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

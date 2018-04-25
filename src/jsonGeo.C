#include "jsonGeo.hh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
std::string itoa(int k)
{
  char buf[24];
  memset(buf,0,24);
  sprintf(buf,"%d",k);
  return std::string(buf);
}
jsonGeo::jsonGeo(std::string config)
{
    std::ifstream ifs (config.c_str(), std::ifstream::in);
  
 // Let's parse it  

    Json::Reader reader;
    
    bool parsedSuccess = reader.parse(ifs, 
				      _jroot, 
				      false);
    
    if(not parsedSuccess)
      {
	// Report failures and their locations 
	// in the document.
	std::cout<<"Failed to parse JSON file ["<<config<<"]"<<std::endl 
		 <<reader.getFormattedErrorMessages()
		 <<std::endl;

      }
    for (int32_t i=0;i<255;i++)
      {
	if (_jroot["difs"].isMember(itoa(i)))
	  {
	    _jDif[i].chamber=_jroot["difs"][itoa(i)]["chamber"].asUInt();
	    _jDif[i].di=_jroot["difs"][itoa(i)]["di"].asDouble();
	    _jDif[i].dj=_jroot["difs"][itoa(i)]["dj"].asDouble();
	    _jDif[i].poli=_jroot["difs"][itoa(i)]["poli"].asDouble();
	    _jDif[i].polj=_jroot["difs"][itoa(i)]["polj"].asDouble();
	    _jDif[i].type=_jroot["difs"][itoa(i)]["type"].asString();
	  }
	if (_jroot["chambers"].isMember(itoa(i)))
	  {
	    _jChamber[i].plan=_jroot["chambers"][itoa(i)]["plan"].asUInt();
	    _jChamber[i].x0=_jroot["chambers"][itoa(i)]["x0"].asDouble();
	    _jChamber[i].y0=_jroot["chambers"][itoa(i)]["y0"].asDouble();
	    _jChamber[i].z0=_jroot["chambers"][itoa(i)]["z0"].asDouble();
	    _jChamber[i].x1=_jroot["chambers"][itoa(i)]["x1"].asDouble();
	    _jChamber[i].y1=_jroot["chambers"][itoa(i)]["y1"].asDouble();
	    _jChamber[i].z1=_jroot["chambers"][itoa(i)]["z1"].asDouble();
	    
	  }
      }

  }
void jsonGeo::convert(uint32_t difid,uint32_t asicid,uint32_t ipad,ROOT::Math::XYZPoint* p)
{
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Upper view (ASICs SIDE)
    //Asics mapping (channel 0 to 63) on small chambers  with HARDROC I
    //                            0		    5		  10		15	      20        25        30        35        40        45        50        55        60
    const unsigned short MapJSmallHR1[64]={2,2,2,7,5,1,3,1,4,1,6,0,3,0,3,0,7,5,4,0,6,7,2,5,3,1,4,4,6,7,5,6,7,5,6,4,3,7,1,0,5,2,6,7,0,5,4,6,0,6,0,4,4,3,3,3,5,7,1,1,2,1,2,2};
    const unsigned short MapISmallHR1[64]={7,6,5,7,7,7,7,5,7,6,7,7,5,5,6,6,6,6,6,4,6,5,4,5,4,4,5,4,5,4,4,4,3,3,3,3,3,2,3,3,2,3,2,1,0,1,2,1,2,0,1,1,0,1,2,0,0,0,2,1,2,0,1,0};
    //(J Axis)
    //
    // 7	|11|05|00|06|08|04|10|03|  	
    // 6	|15|09|01|14|18|17|20|16|  	     
    // 5	|13|07|02|12|26|23|28|21|  	     
    // 4	|19|25|22|24|27|30|31|29|  	     
    // 3	|39|38|41|36|35|33|34|32|  	     TOP VIEW (ASICs SIDE)
    // 2	|48|58|60|54|46|40|42|37|  	
    // 1	|50|59|62|53|51|45|47|43| 
    // 0	|44|61|63|55|52|56|49|57|   
    //		  0  1  2  3  4  5  6  7    (I Axis)  ----->




    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Upper view (ASICs SIDE)
    //ASIC 90 deg rotation in the trigonometric way compare to the MapISmallHR1/MapJSmallHR1.
    //Asics mapping on large chambers with HARDROC I
    //                            0		    5		  10		15	      20        25        30        35        40        45        50        55        60
    const unsigned short MapJLargeHR1[64]={0,1,2,0,0,0,0,2,0,1,0,0,2,2,1,1,1,1,1,3,1,2,3,2,3,3,2,3,2,3,3,3,4,4,4,4,4,5,4,4,5,4,5,6,7,6,5,6,5,7,6,6,7,6,5,7,7,7,5,6,5,7,6,7};
    const unsigned short MapILargeHR1[64]={2,2,2,7,5,1,3,1,4,1,6,0,3,0,3,0,7,5,4,0,6,7,2,5,3,1,4,4,6,7,5,6,7,5,6,4,3,7,1,0,5,2,6,7,0,5,4,6,0,6,0,4,4,3,3,3,5,7,1,1,2,1,2,2};
    //(J Axis)
    //
    // 7	|03|16|21|29|32|37|43|57|  	
    // 6	|10|20|28|31|34|42|47|49|  	     
    // 5	|04|17|23|30|33|40|45|56|  	     
    // 4	|08|18|26|27|35|46|51|52|  	     
    // 3	|06|14|12|24|36|54|53|55|  	     TOP VIEW (ASICs SIDE)
    // 2	|00|01|02|22|41|60|62|63|  	
    // 1	|05|09|07|25|38|58|59|61| 
    // 0	|11|15|13|19|39|48|50|44|   
    //		  0  1  2  3  4  5  6  7    (I Axis)  ----->
    //				|	 |
    //				|DIFF|
    //				|____|	


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Upper view (ASICs SIDE)
    //90 deg rotation in the trigonometric way compare to the MapILargeHR1/MapJLargeHR1.
    //180 deg rotation in the trigonometric way compare to the MapISmallHR1/MapJSmallHR1.
    //Asics mapping on large chambers with HARDROC II & IIB
    //								 0		   5		 10		   15	     20        25        30        35        40        45        50        55        60
    const unsigned short MapJLargeHR2[64]={1,1,2,2,3,3,4,4,5,5,6,6,7,7,4,3,2,0,0,1,0,5,6,7,7,6,5,4,3,2,1,0,0,1,2,3,4,5,6,7,7,6,0,0,1,1,0,2,3,4,5,7,7,6,6,5,5,4,4,3,3,2,2,1};
    const unsigned short MapILargeHR2[64]={1,0,1,0,1,0,1,0,0,1,0,1,0,1,2,2,2,0,1,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,6,5,6,7,5,5,5,5,6,7,6,7,6,7,6,7,7,6,7,6,7};
    //(J Axis)
    //
    // 7	|46|63|61|59|58|56|54|52|  	
    // 6	|43|45|62|60|57|55|53|51|  	     
    // 5	|42|44|47|48|49|50|41|40|  	     
    // 4	|32|33|34|35|36|37|38|39|  	     
    // 3	|31|30|29|28|27|26|25|24|  	     TOP VIEW (ASICs SIDE)
    // 2	|20|19|16|15|14|21|22|23|  	
    // 1	|18|00|02|04|06|09|11|13| 
    // 0	|17|01|03|05|07|08|10|12|   
    //		  0  1  2  3  4  5  6  7    (I Axis)  ----->
    //				|	 |
    //				|DIFF|
    //				|____|	

    unsigned short AsicShiftI[49]={	0,	
					0,	0,	0,	0,	
					8,	8,	8,	8,	
					16,	16,	16,	16,	
					24,	24,	24,	24,	
					32,	32,	32,	32,	
					40,	40,	40,	40,	
					48,	48,	48,	48,	
					56,	56,	56,	56,	
					64,	64,	64,	64,	
					72,	72,	72,	72,	
					80,	80,	80,	80,	
					88,	88,	88,	88};
    unsigned short AsicShiftJ[49]={	0,		
					0,	8,	16,	24,
					24,	16,	8,	0,	
					0,	8,	16,	24,	
					24,	16,	8,	0,	
					0,	8,	16,	24,	
					24,	16,	8,	0,	
					0,	8,	16,	24,	
					24,	16,	8,	0,	
					0,	8,	16,	24,	
					24,	16,	8,	0,	
					0,	8,	16,	24,	
					24,	16,	8,	0

    };

    jsonDifInfo dif=_jDif[difid];
    jsonChamberInfo ch=_jChamber[dif.chamber];
    uint32_t I=0,J=0;
    double pad=1.04125;
    //Small chambers  HR
    if(dif.type.compare("HR0PAD")==0){
      I = MapISmallHR1[ipad]+AsicShiftI[asicid]+1;
      J = MapJSmallHR1[ipad]+AsicShiftJ[asicid]+1;
    }
 
    //First square meter HR1
    if(dif.type.compare("HR1PAD")==0){
      I = MapILargeHR1[ipad]+AsicShiftI[asicid]+1;
      J = MapJLargeHR1[ipad]+AsicShiftJ[asicid]+1;
    }


    //Second Square meter HR2
    if(dif.type.compare("HR2PAD")==0){
      I = MapILargeHR2[ipad]+AsicShiftI[asicid]+1;
      J = MapJLargeHR2[ipad]+AsicShiftJ[asicid]+1;
      J=33-J;
    }
    if (dif.type.compare("MR2PAD")==0){
      int jligne = 3-(asicid-1)/12;
      int icol = (asicid-1)%12;
      int ispad= ((ipad-1)/8)+1;
      int jspad=(ipad-1)%8+1;
      I = icol*8+ispad;
      J =  jligne*8+jspad;
      //j =  (jligne+1)*8-jspad;
    }

    // from DIF (I,J) to chamber
    double fI,fJ;
    if (dif.poli>0)
      fI=I+dif.di;
    else
      fI=(96-I+1)+dif.di;
    if (dif.polj>0)
      fJ=J+dif.dj;
    else
      fJ=(32-J+1)+dif.dj;

    // Now go global
    if (p!=NULL)
      p->SetXYZ(fI*pad+ch.x0,fJ*pad+ch.y0,ch.z0);
    
}

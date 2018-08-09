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
#include "TCanvas.h"
#include "TdcMapping.hh"
#include <fstream>


using namespace zdaq;
tdcrb::tdcrb(std::string dire) : _directory(dire),_run(0),_started(false),_fdIn(-1),_totalSize(0),_event(0),_geo(NULL),_t0(2E50),_t(0),_tspill(0)
			       ,_readoutTotalTime(0),_numberOfMuon(0),_numberOfShower(0),_runType(0),_dacSet(0),_fdOut(-1),_bxId0(0)
{_rh=DCHistogramHandler::instance();
_analyzer= new lydaq::TdcAnalyzer(_rh);
 _mezMap.clear();
 for (uint32_t i=1;i<255;i++)
    { 
      _mezMap.insert(std::pair<uint32_t,std::vector<lydaq::TdcChannel> >(i,std::vector<lydaq::TdcChannel>()));}
}

void tdcrb::geometry(std::string name)
{
  _geo=new jsonGeo(name);
  _analyzer->setGeometry(_geo);
}
void tdcrb::open(std::string filename)
{
  if (_geo==NULL)
    {
      std::cout<<"Please speicfy a geometry"<<std::endl;
      exit(0);
    }
  _fdIn= ::open(filename.c_str(), O_RDONLY | O_NONBLOCK,S_IRWXU);
  if (_fdIn<0)
    {
      perror("Ici No way to store to file :");
      //std::cout<<" No way to store to file"<<std::endl;
      return;
    }  
  _event=0;
  _started=true;
}
void tdcrb::close()
{
  _started=false;
  ::sleep(1);
  if (_fdIn>0)
    {
      ::close(_fdIn);
      _fdIn=-1;
    }


}
uint32_t tdcrb::totalSize(){return _totalSize;}
uint32_t tdcrb::eventNumber(){return _event;}
uint32_t tdcrb::runNumber(){return _run;}
void tdcrb::Read()
{
  int nfile=0;
  _nread=0;
  for (std::vector<std::pair<uint32_t,std::string> >::iterator it=_files.begin();it!=_files.end();it++)
    {
      std::cout<<"NEW File "<<it->first<<" "<<it->second<<std::endl;
      _run=it->first;
      std::stringstream sff;
      sff<<"sudo chmod o+r "<<it->second;
      system(sff.str().c_str());
      this->open(it->second);
      if (_fdOut<=0)
	this->read();
      //else
      //this->streamout(4);
      
      this->close();
      // std::stringstream sroot;
      // sroot<<"/tmp/histo"<<_run<<"_"<<nfile<<".root";
      //  _rh->writeHistograms(sroot.str());

      //  std::cout<<"Writing histos "<<sroot.str()<<std::endl;
      //getchar();
       nfile++;
    }
}
#define FEBCMS
#ifndef FEBCMS
void tdcrb::read()
{
  
  zdaq::buffer b(0x100000);
  int last=-1;
  uint64_t _eventChannel[4096*8];
  std::vector<lydaq::TdcChannel> _vAll;
  uint32_t _eventChannels;
  _geo->fillFebs(_run);
  while (_started)
    {
      if (!_started) return;
      uint32_t theNumberOfDIF=0;
      // To be implemented
      if (_fdIn>0)
	{
	  _idx=0;

	  
	  int ier=::read(_fdIn,&_event,sizeof(uint32_t));
	  _nread++;
	  if (ier<0 || last==_event)
	    {
	      printf("Cannot read Event anymore %d %d %d \n ",ier,last,_event);return;
	    }
	  //else
	  last=_event;
	  if (_event%100==0)
	    INFO_PRINTF("Event read %d \n",_event);
      
	  ier=::read(_fdIn,&theNumberOfDIF,sizeof(uint32_t));
	  if (ier<0)
	    {
	      printf("Cannot read anymore number of DIF %d \n ",ier);return;
	    }
	  else
	    if (_event%100==0)
	      DEBUG_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
	  DEBUG_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
	  uint32_t difFound[256];
	  memset(difFound,0,256*sizeof(uint32_t));
	  uint32_t trigFound[256];
	  memset(trigFound,0,256*sizeof(uint32_t));
	  _analyzer->clear();
	  _vAll.clear();
	  memset(_eventChannel,0,4096*8*sizeof(uint64_t));
	  _eventChannels=0;
	  for (uint32_t idif=0;idif<theNumberOfDIF;idif++) 
	    {
	      uint32_t tbcid=0;
	      //DEBUG_PRINTF("\t writing %d bytes",idata[SHM_BUFFER_SIZE]);
	      //(*iv)->compress();
	      uint32_t bsize=0;
	      // _totalSize+=bsize;
	      ier=::read(_fdIn,&bsize,sizeof(uint32_t));
	      if (ier<0)
		{
		  printf("Cannot read anymore  DIF Size %d \n ",ier);return;
		}
	      else
		if (_event%100==0)
		  DEBUG_PRINTF("\t DIF size %d \n",bsize);
	  
	      ier=::read(_fdIn,b.ptr(),bsize);
	      if (ier<0)
		{
		  printf("Cannot read anymore Read data %d \n ",ier);return;
		}
	      b.setPayloadSize(bsize-(3*sizeof(uint32_t)+sizeof(uint64_t)));
	      b.uncompress();
	      memcpy(&_buf[_idx], b.payload(),b.payloadSize());
	      b.setDetectorId(b.detectorId()&0xFF);
	      DEBUG_PRINTF("\t \t %d %d %d %x %d %d %d\n",b.detectorId()&0XFF,b.dataSourceId(),b.eventId(),b.bxId(),b.payloadSize(),bsize,_idx);
	      _bxId=b.bxId();
	      if (_bxId0==0) _bxId0=_bxId;
	      uint32_t _detId=b.detectorId()&0xFF;
	      //DEBUG_PRINTF("DETID %d \n",_detId);
	      // getchar();
	      if (_detId==255)
		{
		  uint32_t* buf=(uint32_t*) b.payload();
		  printf("NEW RUN %d \n",_event);
		  _run=_event;


		  for (int i=0;i<b.payloadSize()/4;i++)
		    {
		      printf("%d ",buf[i]);
		    }
		  _difId=b.dataSourceId();
		  _runType=buf[0];
		  if (_runType==1)
		    _dacSet=buf[1];
		  if (_runType==2)
		    _vthSet=buf[1];
		  printf("\n Run type %d DAC set %d VTH set %d \n",_runType,_dacSet,_vthSet);
		  // getchar();

		}
	      if (_detId==120)
		{
		  uint32_t* ibuf=(uint32_t*) b.payload();
	       
		   // for (int i=0;i<7;i++)
		   //   {
		   //     printf("%d ",ibuf[i]);
		   //  }
		  uint32_t nch=ibuf[6];
		  //printf("\n channels -> %d \n",nch);
		  _mezzanine=ibuf[4];
		  _difId=(ibuf[5]>>24)&0xFF;
		  _gtc=ibuf[1];

		  DEBUG_PRINTF("\t \t \t %d %d %d %d \n",_mezzanine,_difId,_gtc,nch);
		  std::vector<lydaq::TdcChannel> vch;
		  vch.clear();
		  _analyzer->setInfo(_difId,_run,_event,_gtc,_bxId,TDC_TRIGGER_CHANNEL,_vthSet,_dacSet);
		  if (ibuf[6]>=0)
		    {
		      for (uint32_t i=1;i<255;i++)
			if (_mezMap[i].size()>0) _mezMap[i].clear();
		      uint8_t* cbuf=( uint8_t*)&ibuf[7];
		      bool tfound=false;
		      for (int i=0;i<nch;i++)
			{
			  
			  //for (int j=0;j<8;j++)
			  //  DEBUG_PRINTF("\t %.2x ",cbuf[i*8+j]);
			  //  DEBUG_PRINTF("\n");
			  memcpy(&_eventChannel[_eventChannels],&cbuf[8*i],8*sizeof(uint8_t));
			  lydaq::TdcChannel c(&cbuf[8*i],_difId&0xFF);
			  lydaq::TdcChannel ca((uint8_t*) &_eventChannel[_eventChannels],_difId&0xFF);
			  _eventChannels++;
			  if (c.channel()==24)
			    {
			      tfound=true;
			      trigFound[_difId]++;
			      tbcid=c.bcid();
			    }
			  //_mezMap[_difId].push_back(c);
			  vch.push_back(c);
			  _vAll.push_back(ca);
			  
			}
		      //if (!tfound && _runType==0 && _event%10000!=0 ) continue;
		      if (_runType==1) _analyzer->pedestalAnalysis(_difId,vch);
		      if (_runType==2) _analyzer->scurveAnalysis(_difId,vch);
		      //if (_runType==0) _analyzer->normalAnalysis(_difId,vch);
		    }
		  difFound[ _difId]+=vch.size();
		  if (_event%100==0)
		    DEBUG_PRINTF("Time Acquisition => \t %lu %lu %10.3f \n",_bxId0,_bxId,(_bxId-_bxId0)*2E-7);
		  
		
		  if (_gtc%100==0)
		    DEBUG_PRINTF("\t \t \t ========>Oops Type %d Mez %d DIF %d %x  channels %d Event %d %d\n",_runType,_mezzanine,_difId,ibuf[5],difFound[_difId],_gtc,tbcid);


		  // for (uint32_t i=1;i<255;i++)
		  //   if (_mezMap[i].size()>0)
		  //     {
		  // 	if (_runType==1) _analyzer->pedestalAnalysis(i,_mezMap[i]);
		  // 	if (_runType==2) _analyzer->scurveAnalysis(i,_mezMap[i]);
		  // 	if (_runType==0) _analyzer->normalAnalysis(i,_mezMap[i]);
		  //     }
		
		}

     
	    }
	  //_analyzer->fullAnalysis(_vAll);
	  _analyzer->multiChambers(_vAll);
	  if (_analyzer->trigger())
	    INFO_PRINTF("EVENT SUMMARY \t \t ========>Oops %d total %d, %d, triggers %d %f \n",_event,_eventChannels,_vAll.size(),_analyzer->triggers(),_analyzer->acquisitionTime());
	  //getchar();

	  // bool stop=false;
	  // uint32_t ndf=0;
	  // for (int i=0;i<256;i++)
	  //   if (difFound[i]!=0)
	  //     {DEBUG_PRINTF("(%d:%d) ",i,difFound[i]);
	  // 	ndf++;
	  // 	if (difFound[i]!=1) stop=true;
	  //     }
	  // DEBUG_PRINTF("\n");
	  // if (ndf!=theNumberOfDIF || stop) getchar();
	}

    }
}
#else
void tdcrb::read()
{
  
  zdaq::buffer b(0x100000);
  int last=-1;
  uint64_t _eventChannel[4096*8];
  std::vector<lydaq::TdcChannel> _vAll;
  uint32_t _eventChannels;
  _geo->fillFebs(_run);
  while (_started)
    {
      if (!_started) return;
      uint32_t theNumberOfDIF=0;
      // To be implemented
      if (_fdIn>0)
	{
	  _idx=0;

	  
	  int ier=::read(_fdIn,&_event,sizeof(uint32_t));
	  _nread++;
if (ier<0 || ((last==_event)&_nread>200))
	    {
	      printf("Cannot read Event anymore %d  %d %d read %d \n ",ier,last,_event,_nread);return;
	    }
	  //else
	  last=_event;
	  if (_event%100==0)
	    DEBUG_PRINTF("Event read %d \n",_event);
      
	  ier=::read(_fdIn,&theNumberOfDIF,sizeof(uint32_t));
	  if (ier<0)
	    {
	      printf("Cannot read anymore number of DIF %d \n ",ier);return;
	    }
	  else
	    if (_event%100==0)
	      DEBUG_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
	  INFO_PRINTF("================> %d %d Event %d Number of DIF found %d \n",ier,_fdIn,_event,theNumberOfDIF);
	  uint32_t difFound[256];
	  memset(difFound,0,256*sizeof(uint32_t));
	  uint32_t trigFound[256];
	  memset(trigFound,0,256*sizeof(uint32_t));
	  _analyzer->clear();
	  _vAll.clear();
	  memset(_eventChannel,0,4096*8*sizeof(uint64_t));
	  _eventChannels=0;
	  for (uint32_t idif=0;idif<theNumberOfDIF;idif++) 
	    {
	      uint32_t tbcid=0;
	      //DEBUG_PRINTF("\t writing %d bytes",idata[SHM_BUFFER_SIZE]);
	      //(*iv)->compress();
	      uint32_t bsize=0;
	      // _totalSize+=bsize;
	      ier=::read(_fdIn,&bsize,sizeof(uint32_t));
	      if (ier<0)
		{
		  printf("Cannot read anymore  DIF Size %d \n ",ier);return;
		}
	      else
		if (_event%100<100)
		  INFO_PRINTF("\t DIF size %d \n",bsize);
	  
	      ier=::read(_fdIn,b.ptr(),bsize);
	      if (ier<0)
		{
		  printf("Cannot read anymore Read data %d \n ",ier);return;
		}
	      b.setPayloadSize(bsize-(3*sizeof(uint32_t)+sizeof(uint64_t)));
	      b.uncompress();
	      memcpy(&_buf[_idx], b.payload(),b.payloadSize());
	      b.setDetectorId(b.detectorId()&0xFF);
	      INFO_PRINTF("\t \t det %d source %d event %d bx %x payload %d size  %d %d\n",b.detectorId()&0XFF,b.dataSourceId(),b.eventId(),b.bxId(),b.payloadSize(),bsize,_idx);
	      _bxId=b.bxId();
	      if (_bxId0==0) _bxId0=_bxId;
	      uint32_t _detId=b.detectorId()&0xFF;
	      DEBUG_PRINTF("DUMP DETID %d \n",_detId);
	      uint8_t* cc=(uint8_t*) b.payload();
	      
	      // for (int ib=0;ib<b.payloadSize();ib++)
	      // 	{
	      // 	  printf(" %.2x ",cc[ib]);
	      // 	  if (ib%16==15 &&ib>0) printf("\n");
	      // 	}


	      
	      // getchar();
	      if (_detId==255)
		{
		  uint32_t* buf=(uint32_t*) b.payload();
		  printf("NEW RUN %d \n",_event);
		  _run=_event;


		  for (int i=0;i<b.payloadSize()/4;i++)
		    {
		      printf("%d ",buf[i]);
		    }
		  _difId=b.dataSourceId();
		  _runType=buf[0];
		  if (_runType==1)
		    _dacSet=buf[1];
		  if (_runType==2)
		    _vthSet=buf[1];
		  printf("\n Run type %d DAC set %d VTH set %d \n",_runType,_dacSet,_vthSet);
		  // getchar();

		}
	      if (_detId==130)
		{
		  uint32_t* ibuf=(uint32_t*) b.payload();
	       
		   // for (int i=0;i<7;i++)
		   //   {
		   //     printf("%d ",ibuf[i]);
		   //  }
		  uint32_t nch=ibuf[6];
		  printf("\n channels -> %d \n",nch);
		  _mezzanine=ibuf[4];
		  _difId=(ibuf[5]>>24)&0xFF;
		  _gtc=ibuf[1];

		  INFO_PRINTF("\t \t \t %d %d GTC %d NCH %d \n",_mezzanine,_difId,_gtc,nch);
		  std::vector<lydaq::TdcChannel> vch;
		  vch.clear();
		  _analyzer->setInfo(_difId,_run,_event,_gtc,_bxId,TDC_TRIGGER_CHANNEL,_vthSet,_dacSet);
		  if (ibuf[6]>=0)
		    {
		      for (uint32_t i=1;i<255;i++)
			if (_mezMap[i].size()>0) _mezMap[i].clear();
		      uint8_t* cbuf=( uint8_t*)&ibuf[7];
		      bool tfound=false;
		      for (int i=0;i<nch;i++)
			{
			  
			  // for (int j=0;j<6;j++)
			  //    INFO_PRINTF("\t %.2x ",cbuf[i*6+j]);
			  // INFO_PRINTF("\n");
			  memcpy(&_eventChannel[_eventChannels],&cbuf[6*i],6*sizeof(uint8_t));
			  lydaq::TdcChannel c(&cbuf[6*i],_difId&0xFF);
			  lydaq::TdcChannel ca((uint8_t*) &_eventChannel[_eventChannels],_difId&0xFF);
			  //c.dump();
			 
			  _eventChannels++;
			  if (c.channel()==0)
			    {
			      tfound=true;
			      trigFound[_difId]++;
			      tbcid=c.bcid();
			    }
			  //_mezMap[_difId].push_back(c);
			  vch.push_back(c);
			  _vAll.push_back(ca);
			  
			}
		      // if (nch>0)
		      // getchar();
		      //if (!tfound && _runType==0 && _event%10000!=0 ) continue;
		      if (_runType==1) _analyzer->pedestalAnalysis(_difId,vch);
		      if (_runType==2) _analyzer->scurveAnalysis(_difId,vch);
		      if (_runType==0) _analyzer->normalAnalysis(_difId,vch);
		    }
		  difFound[ _difId]+=vch.size();
		  if (_event%100==0)
		    DEBUG_PRINTF("Time Acquisition => \t %lu %lu %10.3f \n",_bxId0,_bxId,(_bxId-_bxId0)*2E-7);
		  
		
		  if (_gtc%100==0)
		    DEBUG_PRINTF("\t \t \t ========>Oops Type %d Mez %d DIF %d %x  channels %d Event %d %d\n",_runType,_mezzanine,_difId,ibuf[5],difFound[_difId],_gtc,tbcid);

		
		}

     
	    }
	  //_analyzer->fullAnalysis(_vAll);
	  if (_runType==0) _analyzer->multiChambers(_vAll);
	  if (_analyzer->trigger())
	    INFO_PRINTF("EVENT SUMMARY \t \t ========>Oops %d total %d, %d, triggers %d %f \n",_event,_eventChannels,_vAll.size(),_analyzer->triggers(),_analyzer->acquisitionTime());
	}

    }
}
#endif
void tdcrb::end()
{

  _analyzer->end();
  std::stringstream sr;
  sr<<"./histo"<<_run<<"_0.root";
  
  _rh->writeHistograms(sr.str());


  
}
#ifdef TESTMAINEXAMPLE
int main()
{
  levbdim::tdcrb bs("/tmp");
  bs.geometry("/home/acqilc/SDHCAL/SDHCAL_EventReader/pluggins/m3_avril2015.json");
  //bs.open("/data/NAS/June2016/SMM_160616_163121_732786.dat");
  // bs.open("/data/NAS/June2016/SMM_160616_110612_732783.dat");
  //bs.open("/data/NAS/June2016/SMM_170616_052256_732795.dat");
  //bs.open("/data/NAS/June2016/SMM_170616_092331_732799.dat");
  // bs.open("/data/NAS/June2016/SMM_160616_110612_732783.dat");
  /*
    /data/NAS/Oct2016/SMM_071016_123856_733633.dat
    /data/NAS/Oct2016/SMM_071016_124907_733636.dat
    /data/NAS/Oct2016/SMM_071016_125306_733637.dat
    /data/NAS/Oct2016/SMM_071016_153619_733641.dat
    /data/NAS/Oct2016/SMM_071016_154539_733642.dat
    /data/NAS/Oct2016/SMM_071016_155358_733643.dat
    /data/NAS/Oct2016/SMM_071016_155937_733644.dat
    /data/NAS/Oct2016/SMM_071016_165755_733645.dat
    /data/NAS/Oct2016/SMM_071016_170520_733646.dat
    /data/NAS/Oct2016/SMM_071016_173728_733647.dat
    /data/NAS/Oct2016/SMM_071016_193047_733650.dat
    /data/NAS/Oct2016/SMM_071016_205435_733653.dat
    /data/NAS/Oct2016/SMM_071016_210657_733654.dat
    /data/NAS/Oct2016/SMM_071016_232430_733655.dat
    /data/NAS/Oct2016/SMM_071016_233612_733655.dat
    /data/NAS/Oct2016/SMM_081016_012606_733656.dat
    /data/NAS/Oct2016/SMM_081016_015222_733656.dat
    /data/NAS/Oct2016/SMM_081016_033908_733658.dat
    /data/NAS/Oct2016/SMM_081016_035422_733659.dat
    /data/NAS/Oct2016/SMM_081016_035811_733660.dat
    /data/NAS/Oct2016/SMM_081016_054542_733660.dat
    /data/NAS/Oct2016/SMM_081016_082718_733665.dat
    /data/NAS/Oct2016/SMM_081016_092637_733665.dat
    /data/NAS/Oct2016/SMM_081016_110948_733665.dat
    /data/NAS/Oct2016/SMM_081016_160612_733675.dat
    /data/NAS/Oct2016/SMM_081016_171614_733678.dat
    /data/NAS/Oct2016/SMM_081016_173543_733679.dat
    /data/NAS/Oct2016/SMM_091016_010348_733680.dat
    /data/NAS/Oct2016/SMM_091016_041603_733680.dat
    /data/NAS/Oct2016/SMM_091016_065059_733680.dat
    /data/NAS/Oct2016/SMM_091016_072800_733683.dat
    /data/NAS/Oct2016/SMM_091016_104731_733686.dat
    /data/NAS/Oct2016/SMM_091016_122843_733688.dat
    /data/NAS/Oct2016/SMM_091016_124430_733688.dat
    /data/NAS/Oct2016/SMM_091016_140032_733689.dat
    /data/NAS/Oct2016/SMM_091016_154423_733689.dat
    /data/NAS/Oct2016/SMM_091016_163928_733692.dat
    /data/NAS/Oct2016/SMM_091016_164335_733693.dat
    /data/NAS/Oct2016/SMM_091016_182544_733693.dat
    /data/NAS/Oct2016/SMM_091016_184359_733696.dat
    /data/NAS/Oct2016/SMM_091016_202828_733696.dat
    /data/NAS/Oct2016/SMM_091016_211028_733698.dat
    /data/NAS/Oct2016/SMM_091016_223807_733698.dat
    /data/NAS/Oct2016/SMM_091016_232004_733699.dat
    /data/NAS/Oct2016/SMM_101016_001144_733700.dat
    /data/NAS/Oct2016/SMM_101016_001759_733701.dat
    /data/NAS/Oct2016/SMM_101016_004241_733702.dat
    /data/NAS/Oct2016/SMM_101016_012948_733705.dat
    /data/NAS/Oct2016/SMM_101016_013700_733707.dat
    /data/NAS/Oct2016/SMM_101016_014138_733708.dat
    /data/NAS/Oct2016/SMM_101016_015433_733710.dat
    /data/NAS/Oct2016/SMM_101016_020437_733711.dat
    /data/NAS/Oct2016/SMM_101016_033112_733711.dat
    /data/NAS/Oct2016/SMM_101016_045829_733711.dat
    /data/NAS/Oct2016/SMM_101016_054606_733718.dat
    /data/NAS/Oct2016/SMM_101016_063034_733718.dat
    /data/NAS/Oct2016/SMM_101016_075846_733718.dat
    /data/NAS/Oct2016/SMM_101016_092544_733719.dat
    /data/NAS/Oct2016/SMM_101016_103528_733720.dat
    /data/NAS/Oct2016/SMM_101016_111905_733722.dat
    /data/NAS/Oct2016/SMM_101016_113444_733723.dat
    /data/NAS/Oct2016/SMM_101016_120745_733724.dat
    /data/NAS/Oct2016/SMM_101016_132143_733724.dat
    /data/NAS/Oct2016/SMM_101016_144929_733725.dat
    /data/NAS/Oct2016/SMM_101016_151024_733728.dat
    /data/NAS/Oct2016/SMM_101016_175002_733738.dat
    /data/NAS/Oct2016/SMM_101016_175744_733740.dat
    /data/NAS/Oct2016/SMM_101016_181441_733741.dat
    /data/NAS/Oct2016/SMM_101016_181943_733742.dat
    /data/NAS/Oct2016/SMM_101016_224057_733743.dat
    /data/NAS/Oct2016/SMM_111016_004226_733743.dat
    /data/NAS/Oct2016/SMM_111016_005214_733748.dat
    /data/NAS/Oct2016/SMM_111016_022420_733750.dat
    /data/NAS/Oct2016/SMM_111016_024904_733750.dat
    /data/NAS/Oct2016/SMM_111016_044205_733750.dat
    /data/NAS/Oct2016/SMM_111016_063534_733750.dat
    /data/NAS/Oct2016/SMM_111016_072721_733754.dat
    /data/NAS/Oct2016/SMM_111016_084258_733754.dat
    /data/NAS/Oct2016/SMM_111016_105343_733756.dat
    /data/NAS/Oct2016/SMM_111016_113540_733756.dat
    /data/NAS/Oct2016/SMM_111016_135741_733757.dat
    /data/NAS/Oct2016/SMM_111016_140953_733758.dat
    /data/NAS/Oct2016/SMM_111016_151355_733759.dat
    acqilc@lyosdhcal9:~$ 

  */
  bs.open("/data/srv02/RAID6/Oct2016/SMM_101016_224057_733743.dat");
  bs.read();
}
#endif

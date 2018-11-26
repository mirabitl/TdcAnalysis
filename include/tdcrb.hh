#ifndef _zdaq_tdcrb_h
#define _zdaq_tdcrb_h

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include "jsonGeo.hh"
#include "DCHistogramHandler.hh"
#include "TdcChannel.hh"
#include "TdcAnalyzer.hh"
#include "zmBuffer.hh"
#ifndef FEBCMS
#define TDC_TRIGGER_CHANNEL 24
#else
#define TDC_TRIGGER_CHANNEL 0
#endif
class tdcrb
{
public:
  tdcrb(std::string dire="/tmp");
  void Read();
  void end();
  void open(std::string name);
  void close();
  void read();
   

  void geometry(std::string name);
    
  uint32_t totalSize();
  uint32_t eventNumber();
  uint32_t runNumber();
  void addRun(uint32_t r,std::string name) { _files.push_back(std::pair<uint32_t,std::string>(r,name));}
  void setRun(int r){_run=r;}
  void setOutFileId(int32_t fid){_fdOut=fid;}
  void setAnalyzer(lmana::Analyzer* a){_analyzer=a;}
  void clearDataSet(){_files.clear();}
  void findDataSet(std::string dir,uint32_t run);
  void stop(){_started=false;}
  void processRawEvent(uint64_t idx);
  void clearShm();
  void monitor();
  void pull(std::string name,zdaq::buffer* buf,std::string sourcedir);
  uint32_t numberOfDataSource();

private:
  std::vector<std::pair<uint32_t,std::string> > _files;
  uint64_t _bxId,_bxId0;
  uint32_t _gtc;
  double _t,_t0,_tspill;
  std::string _directory;
  uint32_t _run,_event,_totalSize,_nread;
  uint32_t _nevt,_ntrigger,_nfound,_nbside;
  int32_t _fdIn,_fdOut;
  bool _started;
  unsigned char _buf[32*1024*1024];
  uint32_t _idx;
  jsonGeo* _geo;
  std::map<uint32_t,std::bitset<64> > _tcount;
  lmana::Analyzer* _analyzer;
  double _readoutTime,_readoutTotalTime;
  uint32_t _numberOfShower,_numberOfMuon;
  DCHistogramHandler* _rh;
  uint32_t _runType,_dacSet,_vthSet,_mezzanine,_difId;
  
  std::map<uint32_t,std::vector<lydaq::TdcChannel> > _mezMap;
  std::map<uint64_t,std::vector<zdaq::buffer*> > _eventMap;

};
#endif

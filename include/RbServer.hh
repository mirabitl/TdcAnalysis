#ifndef _RbServer_h

#define _RbServer_h
#include <iostream>

#include <string.h>
#include<stdio.h>
#include "baseApplication.hh"
#include "ReadoutLogger.hh"

#include "tdcrb.hh"

using namespace std;
#include <sstream>



  class RbServer : public zdaq::baseApplication
  {
  public:
    RbServer(std::string name);
    void configure(zdaq::fsmmessage* m);
    void destroy(zdaq::fsmmessage* m);
    void process(zdaq::fsmmessage* m);
    void startMonitor(zdaq::fsmmessage* m);
    void processRun();
    void monitorRun();
    void stop(zdaq::fsmmessage* m);
    void c_status(Mongoose::Request &request, Mongoose::JsonResponse &response);
    void c_histolist(Mongoose::Request &request, Mongoose::JsonResponse &response);
    void c_histo(Mongoose::Request &request, Mongoose::JsonResponse &response);
    
  private:

    zdaq::fsmweb* _fsm;
    boost::thread_group _g_group;
    tdcrb* _bs;
    DCHistogramHandler* _rh;
    std::string _directory;
    std::string _geofile;
    std::string _analyzer;
    uint32_t _run;
  };
#endif


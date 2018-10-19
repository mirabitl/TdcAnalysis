
#include "RbServer.hh"
using namespace zdaq;
using namespace lmana;

void RbServer::configure(zdaq::fsmmessage* m)
{
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" CMD: "<<m->content());

  if (m->content().isMember("directory"))
    { 
      _directory=m->content()["directory"].asString();
      this->parameters()["directory"]=m->content()["directory"];
    }
  else
    _directory=this->parameters()["directory"].asString();



  if (m->content().isMember("geometry"))
    { 
      _geofile=m->content()["geometry"].asString();
      this->parameters()["geometry"]=m->content()["geometry"];
    }
  else
    _geofile=this->parameters()["geometry"].asString();

  if (m->content().isMember("analyzer"))
    { 
      _analyzer=m->content()["analyzer"].asString();
      this->parameters()["analyzer"]=m->content()["analyzer"];
    }
  else
    _analyzer=this->parameters()["analyzer"].asString();

  _bs= new tdcrb("/tmp");


  if (_analyzer.compare("TdcAnalyzer")==0)
    {
      lmana::Analyzer* a = new TdcAnalyzer(_rh);
      _bs->setAnalyzer(a);
    }
  if (_analyzer.compare("RecoAnalyzer")==0)
    {
      lmana::Analyzer* a = new RecoAnalyzer(_rh);
      _bs->setAnalyzer(a);
    }
    _bs->setOutFileId(-1);
    _bs->geometry(_geofile);
}
void RbServer::destroy(zdaq::fsmmessage* m)
{
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" CMD: "<<m->command());
  if (_bs==NULL)
    {
       LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<"Please open Configure first");
       return;
    }
  //_bs->end();
  delete _bs;
}
void RbServer::processRun()
{
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" Calling read "<<_run);
  _bs->Read();
}

void RbServer::stop(zdaq::fsmmessage* m)
{
  _bs->stop();
  _g_group.join_all();
  _bs->end();
}

void RbServer::process(zdaq::fsmmessage* m)
{
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" CMD: "<<m->command());
  if (_bs==NULL)
    {
       LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<"Please open Configure first");
       return;
    }
  if (m->content().isMember("run"))
    { 
      _run=m->content()["run"].asUInt();

    }
  else
    {
    LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<"Please provide a run number");
    return;
    }
  LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<"RUn"<<_run);
  _bs->clearDataSet();
  _bs->findDataSet(_directory,_run);
  _bs->setRun(_run);
  LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<" launching thread"<<_run);
  // No launch the process
  _g_group.create_thread(boost::bind(&RbServer::processRun, this));
  //_g_group.join_all();
  return;

}
void RbServer::c_status(Mongoose::Request &request, Mongoose::JsonResponse &response)
{

  response["STATUS"]="DONE";
  if (_bs!=NULL)
    {
      response["run"]=_bs->runNumber();
      response["event"]=_bs->eventNumber();
    }
  else
    {
      response["run"]=-1;
      response["event"]=-1;
    }
    
}


void RbServer::c_histolist(Mongoose::Request &request, Mongoose::JsonResponse &response)
{
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" Histo List called ");

  if (_rh==NULL)
    {
       LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<"No Hisogram handler");
       response["STATUS"]="FAILED";
       return;
    }

  response["STATUS"]="DONE";
  response["list"]=_rh->getXMLHistoList();
  
}
void RbServer::c_histo(Mongoose::Request &request, Mongoose::JsonResponse &response)
{
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" Histo called ");

  if (_rh==NULL)
    {
       LOG4CXX_ERROR(_logMDCC,__PRETTY_FUNCTION__<<"No Hisogram handler");
       response["STATUS"]="FAILED";
       return;
    }
  std::string sh=request.get("histo","none");
  LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" Histo called for "<<sh);
  response["STATUS"]="DONE";
  if (sh.compare("none")!=0)
    response["histo"]=_rh->getJSONHisto(sh);
  else
    response["STATUS"]="MISSED";
  
}

RbServer::RbServer(std::string name) : zdaq::baseApplication(name),_bs(NULL),_rh(NULL)
{

  
 

  //_fsm=new zdaq::fsm(name);
  _fsm=this->fsm();

  
// Register state
  _fsm->addState("CONFIGURED");
  _fsm->addState("RUNNING");


  _fsm->addTransition("CONFIGURE","CREATED","CONFIGURED",boost::bind(&RbServer::configure, this,_1));
  _fsm->addTransition("PROCESS","CONFIGURED","RUNNING",boost::bind(&RbServer::process, this,_1));
  _fsm->addTransition("STOP","RUNNING","CONFIGURED",boost::bind(&RbServer::stop, this,_1));
  _fsm->addTransition("DESTROY","CONFIGURED","CREATED",boost::bind(&RbServer::destroy, this,_1));
  
 _fsm->addCommand("STATUS",boost::bind(&RbServer::c_status,this,_1,_2));

 _fsm->addCommand("HISTOLIST",boost::bind(&RbServer::c_histolist,this,_1,_2));
 _fsm->addCommand("HISTO",boost::bind(&RbServer::c_histo,this,_1,_2));

 _rh=DCHistogramHandler::instance();

  char* wp=getenv("WEBPORT");
  if (wp!=NULL)
    {
      LOG4CXX_INFO(_logMDCC,__PRETTY_FUNCTION__<<" Service "<<name<<" is starting on "<<atoi(wp));
      //      std::cout<<"Service "<<name<<" started on port "<<atoi(wp)<<std::endl;
    _fsm->start(atoi(wp));
    }


}



#include "tdcreadbinary.hh"
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
#include "RecoHit.hh"
#include "recoPoint.hh"
#include "recoTrack.hh"
#include "rCluster.hh"
#include "TCanvas.h"

class framePoint : public recoPoint
{
public:
  framePoint(uint32_t p) : _plan(p){;}
  virtual double dX(){return 0.3;}
  virtual double dY(){return 0.3;}
  virtual uint32_t plan(){ return _plan;}
private:
  uint32_t _plan;
};
using namespace levbdim;
tdcreadbinary::tdcreadbinary(std::string dire) : _directory(dire),_run(0),_started(false),_fdIn(-1),_totalSize(0),_event(0),_geo(NULL),_t0(2E50),_t(0),_tspill(0)
					       ,_readoutTotalTime(0),_numberOfMuon(0),_numberOfShower(0),_runType(0),_dacSet(0)
{_rh=DCHistogramHandler::instance();}

void tdcreadbinary::geometry(std::string name)
{
  _geo=new jsonGeo(name);
}
void tdcreadbinary::open(std::string filename)
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
void tdcreadbinary::close()
{
  _started=false;
  ::sleep(1);
  if (_fdIn>0)
    {
      ::close(_fdIn);
      _fdIn=-1;
    }


}
uint32_t tdcreadbinary::totalSize(){return _totalSize;}
uint32_t tdcreadbinary::eventNumber(){return _event;}
uint32_t tdcreadbinary::runNumber(){return _run;}
void tdcreadbinary::Read()
{
  for (std::vector<std::pair<uint32_t,std::string> >::iterator it=_files.begin();it!=_files.end();it++)
    {
      std::cout<<"NEW File "<<it->first<<" "<<it->second<<std::endl;
      _run=it->first;
      std::stringstream sff;
      sff<<"sudo chmod o+r "<<it->second;
      system(sff.str().c_str());
      this->open(it->second);
      this->read();
      this->close();
    }
}
void tdcreadbinary::read()
{
  
  levbdim::buffer b(0x100000);
  while (_started)
    {
      if (!_started) return;
      uint32_t theNumberOfDIF=0;
      // To be implemented
      if (_fdIn>0)
	{
	  _idx=0;
	  // Clear previous event
	  for (uint32_t i=0;i<theDIFPtrList_.size();i++) 
	    {   
	      //theDIFPtrList_[i]->dumpDIFInfo();
	      delete theDIFPtrList_[i];
	    }
	  theDIFPtrList_.clear();
 
	  int ier=::read(_fdIn,&_event,sizeof(uint32_t));
	  if (ier<=0)
	    {
	      printf("Cannot read anymore %d \n ",ier);return;
	    }
	  //else
	  //	printf("Event read %d \n",_event);
      
	  ier=::read(_fdIn,&theNumberOfDIF,sizeof(uint32_t));
	  if (ier<=0)
	    {
	      printf("Cannot read anymore number of DIF %d \n ",ier);return;
	    }
	  // else
	  //   printf("Number of DIF found %d \n",theNumberOfDIF);

	  for (uint32_t idif=0;idif<theNumberOfDIF;idif++) 
	    {
	      //printf("\t writing %d bytes",idata[SHM_BUFFER_SIZE]);
	      //(*iv)->compress();
	      uint32_t bsize=0;
	      // _totalSize+=bsize;
	      ier=::read(_fdIn,&bsize,sizeof(uint32_t));
	      if (ier<=0)
		{
		  printf("Cannot read anymore  DIF Size %d \n ",ier);return;
		}
	      //else
	      //  printf("\t DIF size %d \n",bsize);
	  
	      ier=::read(_fdIn,b.ptr(),bsize);
	      if (ier<=0)
		{
		  printf("Cannot read anymore Read data %d \n ",ier);return;
		}
	      b.setPayloadSize(bsize-(3*sizeof(uint32_t)+sizeof(uint64_t)));
	      b.uncompress();
	      memcpy(&_buf[_idx], b.payload(),b.payloadSize());
	      //printf("\t \t %d %d %d %x %d %d %d\n",b.detectorId(),b.dataSourceId(),b.eventId(),b.bxId(),b.payloadSize(),bsize,_idx);
	      if (b.detectorId()==255)
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

		}
	      if (b.detectorId()==110)
		{
		  uint32_t* ibuf=(uint32_t*) b.payload();
	       
		  // for (int i=0;i<7;i++)
		  //   {
		  //     printf("%d ",ibuf[i]);
		  //  }
		  uint32_t nch=ibuf[6];
		  //printf("\n channels -> %d \n",nch);
		  _channels.clear();	  
		  if (ibuf[6]>0)
		    {

		      uint8_t* cbuf=( uint8_t*)&ibuf[7];
		      for (int i=0;i<nch;i++)
			{
			  /*
			    for (int j=0;j<8;j++)
			    printf("\t %.2x ",cbuf[i*8+j]);
			    printf("\n");
			  */
			  TdcChannel c(&cbuf[8*i]);
			  _channels.push_back(c);
			}
		    }
		      _mezzanine=ibuf[4];
		      _difId=(ibuf[5]>>24)&0xFF;
		      _gtc=ibuf[1];

		      //printf("Type %d Mez %d DIF %d %x  channels %d Event %d \n",_runType,_mezzanine,_difId,ibuf[5],_channels.size(),_gtc);
		      if (_runType==1) this->pedestalAnalysis();
		      if (_runType==2) this->scurveAnalysis();
		      if (_runType==0) this->normalAnalysis();
		
		}
	      if (b.detectorId()!=100) continue;
	      if (idif==0)
		{
		  _bxId=b.bxId();
		  _gtc=b.eventId();
	      
		  if (_t0>1E50)
		    _t0=b.bxId()*2E-7;
		  double ct=b.bxId()*2E-7-_t0;
		  if ((ct-_t)>5.)
		    {
		      _tspill=ct;
		      std::cout<<" New Spill====>"<<_tspill<<std::endl;
		    }
		  _t=ct;
		}
	      DIFPtr* d= new DIFPtr(&_buf[_idx],b.payloadSize());
	      _idx+=b.payloadSize();
	      uint8_t* buf=(uint8_t*) b.payload();
	      theDIFPtrList_.push_back(d);
#ifdef BIFDUMP
	      if (b.dataSourceId()==3)
		{
		  for (int i=0;i<b.payloadSize();i++)
		    {
		      printf("%.2x ",buf[i]);
		    }
		  printf("\n %d \n",b.payloadSize());

		}
#endif	  
	    }
	  // Minimal histos
	  continue;
	  TH1* hacqtim= _rh->GetTH1("/BR/AcquistionTime");	
	  if (hacqtim==NULL)
	    {
	      hacqtim=_rh->BookTH1("/BR/AcquistionTime",1000.,0.,5.);
	    }
	  // Filling frame and selecting events
	  _tframe.clear();
	  _tcount.clear();
	  for (std::vector<recoTrack*>::iterator it=_vtk.begin();it!=_vtk.end();it++) delete (*it);

	  _vtk.clear();
	  std::map<uint32_t,std::vector<std::pair<DIFPtr*,uint32_t> > >::iterator ifm=_tframe.end();
	  std::map<uint32_t,std::bitset<64> >::iterator im=_tcount.end();
	  int32_t window=_geo->cuts()["timeWindow"].asInt();
	  for (std::vector<DIFPtr*>::iterator it = theDIFPtrList_.begin();it!=theDIFPtrList_.end();it++) // Loop on DIF
	    {
	      DIFPtr* d = (*it);
	      //uint32_t chid= getChamber(d->getID());
	      // LMTest      uint32_t bc = rint(f->getBunchCrossingTime()/DCBufferReader::getDAQ_BC_Period());
	      uint32_t chid = _geo->difGeo(d->getID())["chamber"].asUInt();
    
	      // Loop on Frames
	      for (uint32_t ifra=0;ifra<d->getNumberOfFrames();ifra++)
		{
		  uint32_t bc=d->getFrameTimeToTrigger(ifra);
		  bool found=false;
		  for (int dt=-window;dt<=window;dt++)
		    {
		      im=_tcount.find(bc+dt);
	  

		      if (im!=_tcount.end())
			{
			  im->second.set(chid);
			  ifm=_tframe.find(bc+dt);
			  ifm->second.push_back(std::pair<DIFPtr*,uint32_t>(d,ifra));
			  found=true;
			  break;
			}
		    }
		  if (found) continue;
		  std::bitset<64> v(0);
		  v.set(chid);
		  std::pair<uint32_t,std::bitset<64> > p(bc,v);
		  _tcount.insert(p);
		  std::vector<std::pair<DIFPtr*,uint32_t> > vf;
		  vf.push_back(std::pair<DIFPtr*,uint32_t>(d,ifra));
		  std::pair<uint32_t,std::vector<std::pair<DIFPtr*,uint32_t> > > pf(bc,vf);
		  _tframe.insert(pf);
       
		}
	    }
    
	  std::vector<uint32_t> seeds;seeds.clear();
	  uint32_t npmin=_geo->cuts()["minPlans"].asUInt();
	  uint32_t itmin=0xFFFFFFFF;
	  uint32_t itmax=0;
 
	  for (std::map<uint32_t,std::bitset<64> >::iterator im=_tcount.begin();im!=_tcount.end();im++)
	    {
	      //
	      //if (im->second[57] || im->second.count()>12)
	      //      std::cout<<im->first<<" "<<im->second.count()<<" "<<im->second<<std::endl;

	      if (im->first<itmin) itmin=im->first;
	      if (im->first>itmax) itmax=im->first;
	      if (im->second.count()>=npmin)
		{
		  this->processEvent(im->first);
		  seeds.push_back(im->first);
		}
	    }
	  //getchar();
	  _readoutTime=(itmax-itmin)*2E-7;
	  hacqtim->Fill(_readoutTime);
	  if (_readoutTime>0 && _readoutTime<5.)
	    _readoutTotalTime+=_readoutTime;
	  std::cout<<"seeds : "<<seeds.size()<<" readout time :"<<_readoutTime<<" absolute "<<_t<<" in spill "<<_t-_tspill<<" Sh (n:rate) "<<_numberOfShower<<":"<<_numberOfShower/_readoutTotalTime<<" Mu(n:rate) "<<_numberOfMuon<<":"<<_numberOfMuon/_readoutTotalTime<<" total time "<<_readoutTotalTime<<std::endl;
	  if (_event%100==0 &&_event>0)
	    _rh->writeSQL();
     
	}
    }

}

void tdcreadbinary::pedestalAnalysis()
{

  //if (_gtc[_mezzanine-1]
  std::cout<<"Mezzanine "<<_mezzanine<<"Event "<<_event<<" GTC"<<_gtc<<" hits"<<_channels.size()<<std::endl;

  // Analyze
  std::stringstream sr;
  sr<<"/run"<<_run<<"/TDC"<<_mezzanine<<"/";

  uint32_t dac =_dacSet;
  for (int ich=0;ich<28;ich++)
    {
 
      std::stringstream src;
      src<<sr.str()<<"dac"<<ich;
      TH1* hdac=_rh->GetTH1(src.str());
      if (hdac==NULL)
	{
	 
	  hdac=_rh->BookTH1(src.str(),64,0.,64.);
	}
      bool found=false;
      double lastf=0;
      for (auto x :_channels)
	{
	  if (x.channel()==ich) {
	    //printf("%d %d %f \n",x.channel(),x.bcid(),x.tdcTime());
	    double dt=x.tdcTime()-lastf;
	    lastf=x.tdcTime();
	    if (dt>25 || dt<0)
	      hdac->Fill(dac*1.);
	  }//break;}
	}
    }

}
void tdcreadbinary::scurveAnalysis()
{

  //if (_gtc[_mezzanine-1]
  std::cout<<"Mezzanine "<<_difId<<"Event "<<_event<<" GTC"<<_gtc<<" hits"<<_channels.size()<<std::endl;

  // Analyze
  std::stringstream sr;
  sr<<"/run"<<_run<<"/TDC"<<_difId<<"/";

  uint32_t vth =_vthSet;
  for (int ich=0;ich<28;ich++)
    {
 
      std::stringstream src;
      src<<sr.str()<<"vth"<<ich;
      TH1* hvth=_rh->GetTH1(src.str());
      if (hvth==NULL)
	{
	 
	  hvth=_rh->BookTH1(src.str(),900,0.,900.);
	}
      bool found=false;
      double lastf=0;
      for (auto x :_channels)
	{
	  if (x.channel()==ich) {
	    //printf("%d %d %f \n",x.channel(),x.bcid(),x.tdcTime());
	    double dt=x.tdcTime()-lastf;
	    lastf=x.tdcTime();
	    //  if (dt>25 || dt<0)
	    hvth->Fill(vth*1.);
	    break;
	  }
	}
    }

}
void tdcreadbinary::normalAnalysis()
{
  this->LmAnalysis();
}
void tdcreadbinary::processEvent(uint32_t iseed)
{
  std::vector<recoPoint*> point;
  //double bestX1=0.049985, bestY1 =0.049985,bestZ1 =0.049985,bestX2 =2.63409e-05, bestY2=3.00341e-05,bestZ2=7.91877e-05,bestX3=-1.75883e-08, bestY3=-9.99385e-09,bestZ3 =-9.99385e-09;
  double bestX1=0.033; double bestY1 =0.0905167;double bestZ1 =0.136332;double bestX2 =2.69635e-05; double bestY2=9.13179e-06;double bestZ2=3.056545e-05;double bestX3=-1.7897e-08; double bestY3= 1.49317e-09;double bestZ3 =2.599392e-08;
  double par[9];
  par[0]= 2.859E-02;// 1.889E-02 
  par[1]= 4.376E-05;// 4.794E-05 
  par[2]= -2.248E-08;// 2.660E-08 
  par[3]= 1.533E-01;// 6.956E-02 
  par[4]= -3.566E-04;// 1.819E-04 
  par[5]= 2.174E-07;// 1.020E-07 
  par[6]= -2.665E-01;// 8.000E-02 
  par[7]= 1.456E-03;// 2.061E-04 
  par[8]= -8.541E-07;// 1.186E-07 
  par[0]= 4.609E-02;// 5.274E-03 
  par[1]= -4.682E-06;// 7.414E-06 
  par[2]= -1.366E-09;// 2.253E-09


  par[3]= 3.461E-02;// 1.744E-02 
  par[4]= 4.376E-05;// 2.096E-05 
  par[5]= -1.419E-08;// 2.931E-08

 
  par[6]= 7.173E-02;// 1.256E-02 
  par[7]= 3.515E-04;// 1.073E-05 
  par[8]= -1.723E-07;// 1.379E-08 

  uint32_t chbif=_geo->cuts()["bifChamber"].asUInt();
  //std::cout<<chbif<<std::endl;
  TH1* hnhit= _rh->GetTH1("/BR/nhit");
  TH1* hhitparasic= _rh->GetTH1("/BR/hitparasic");
  TH1* hnhitm= _rh->GetTH1("/BR/nhitmuon");
  TH1* hnhits= _rh->GetTH1("/BR/nhitshower");
  TH1* hmsi= _rh->GetTH1("/BR/sizemuon");
  TH1* hssi= _rh->GetTH1("/BR/sizeshower");

  TH2* hhitvsratio=_rh->GetTH2("/BR/hitvsratio");
  TH1* hen= _rh->GetTH1("/BR/energy");
  TH1* henb= _rh->GetTH1("/BR/energyb");

  TProfile* hr0 = (TProfile*) _rh->GetTH1("/BR/r0");
  TProfile* hr1 = (TProfile*) _rh->GetTH1("/BR/r1");
  TProfile* hr2 = (TProfile*) _rh->GetTH1("/BR/r2");
  TH2* hshe=_rh->GetTH2("/BR/showerentry");
    
  
  // if (hhtx==NULL)
  //   {
  //     hhtx = _rh->BookProfile("HoughTransformX",100,0.,100,-50.,50.);
  if (hnhit==NULL)
    {
      hnhit=_rh->BookTH1("/BR/nhit",1000.,0.,2000.);
      hnhitm=_rh->BookTH1("/BR/nhitmuon",1000.,0.,2000.);
      hnhits=_rh->BookTH1("/BR/nhitshower",1000.,0.,2000.);
    
      hmsi=_rh->BookTH1("/BR/sizemuon",128.,0.1,128.1);
      hssi=_rh->BookTH1("/BR/sizeshower",128.,0.1,128.1);

      hhitparasic=_rh->BookTH1("/BR/hitparasic",512.,0.,64.);
      hhitvsratio=_rh->BookTH2("/BR/hitvsratio",512,0.,10.,512,0.,0.5);
      hen=_rh->BookTH1("/BR/energy",500.,0.,150.);
      henb=_rh->BookTH1("/BR/energyb",500.,0.,150.);
      hr0 = _rh->BookProfile("/BR/r0",30,0.,5.,0.,5000.);
      hr1 = _rh->BookProfile("/BR/r1",30,0.,5.,0.,5000.);
      hr2 = _rh->BookProfile("/BR/r2",30,0.,5.,0.,5000.);
      hshe=_rh->BookTH2("/BR/showerentry",100,0.,100.,100,0.,100.);

    }
  std::map<uint32_t,std::vector<std::pair<DIFPtr*,uint32_t> > >::iterator ifm=_tframe.find(iseed);
  if (ifm==_tframe.end())
    return;
  std::map<uint32_t,std::bitset<64> >::iterator im=_tcount.find(iseed);
  if (im==_tcount.end())
    return;
  uint32_t nshower=0,nmuon=0;	
  uint32_t nhit=0,nh0=0,nh1=0,nh2=0;
  uint32_t np[60];memset(np,0,60*4);
  for (std::vector<std::pair<DIFPtr*,uint32_t> >::iterator itf=ifm->second.begin();itf!=ifm->second.end();itf++)
    {
      DIFPtr* d=itf->first;
      uint32_t ifra=itf->second;
      //std::cout<<hex<<(uint64_t) d<<" "<<dec<<ifra<<std::endl;
      uint32_t chid = _geo->difGeo(d->getID())["chamber"].asUInt();
      np[chid]++;
      for (uint32_t j=0;j<64;j++)
	{
	  if (!(d->getFrameLevel(ifra,j,0) || d->getFrameLevel(ifra,j,1))) continue;
	  nhit++;
	  if (d->getFrameLevel(ifra,j,0) && d->getFrameLevel(ifra,j,1)) 
	    nh2++;
	  else
	    if (d->getFrameLevel(ifra,j,0)) 
	      nh1++;
	    else
	      if (d->getFrameLevel(ifra,j,1)) 
		nh0++;
      
	  framePoint* p= new framePoint(chid);
      
	  _geo->convert(d->getID(),d->getFrameAsicHeader(ifra),j,p);
	  std::stringstream s;
	  s<<"/SH/CH"<<chid;
	  TH2* hpos=_rh->GetTH2(s.str()+"/pos");
	  if (hpos==NULL)
	    {
	      hpos=_rh->BookTH2(s.str()+"/pos",120.,-10.,110.,120,-10.,110.);
	    }
	  hpos->Fill(p->X(),p->Y());
	  point.push_back(p);
	}
      //std::cout<<ifra<<" "<<nhit<<std::endl;
    }
  if (point.size()<4) return;
  pcaComponents cp=RecoHit::calculateComponents(point);
  //std::cout<<cp[3]<<" "<<cp[4]<<" "<<cp[5]<<std::endl;
  double bx=cp[0]-cp[6]*cp[2]/cp[8],by=cp[1]-cp[7]*cp[2]/cp[8];
  hshe->Fill(bx,by);
  hnhit->Fill(nhit*1.);
  double hitsparasic=(nh1+nh2+nh0)*1./ifm->second.size();
  hitsparasic=0;
  int nplane=0;
  for (int i=0;i<60;i++)
    {
      if (np[i]==0) continue;
      hitsparasic+=np[i];
      nplane++;
    }
  if (nplane>0) hitsparasic/=nplane;
  hhitparasic->Fill(hitsparasic);
  hhitvsratio->Fill(hitsparasic,(cp[3])/cp[5]);
  
  // Build clusters
  std::vector<rCluster<recoPoint>*> clusters;
  for (std::vector<recoPoint*>::iterator it=point.begin();it!=point.end();it++)
    {
      bool found=false;
      for (std::vector<rCluster<recoPoint>*>::iterator iv=clusters.begin();iv!=clusters.end();iv++)
	if ((*iv)->Append((*it),3.)){found=true;break;}
      if (found) continue;
      clusters.push_back(new rCluster<recoPoint>((*it)));
    }
  uint32_t nbad=0,ngood=0;
  for (std::vector<rCluster<recoPoint>*>::iterator iv=clusters.begin();iv!=clusters.end();iv++)
    {
      if ((*iv)->plan()==1)
	{
	  if ((*iv)->X()>50 && (*iv)->X()<70  &&(*iv)->Y()>45 && (*iv)->X()<60)
	    ngood++;
	  else
	    nbad++;
	}
    }
  // std::cout<<"Clusters :"<<clusters.size()<<" Ratio "<<point.size()*1./clusters.size()<<std::endl;
  
  /** Track */
   

  std::vector<recoPoint*> vrp;
  float zmin=1E20,zmax=-1E20;
   
  for (std::vector<rCluster<recoPoint>*>::iterator iv=clusters.begin();iv!=clusters.end();iv++) 
    {
      uint32_t nv=0;
      for (std::vector<rCluster<recoPoint>*>::iterator jv=clusters.begin();jv!=clusters.end();jv++) 
	{
	  if ((*iv)==(*jv)) continue;
	  ROOT::Math::XYZVector d1=(*(*iv))-(*(*jv));
	  if (d1.Mag2()>100) continue;
	  nv++;				    
	}
       
      //std::cout<<"Size "<<(*iv)->size()<<" vois "<<nv<<std::endl; 
      if ((*iv)->size()<3 || ((*iv)->size()<5 && nv<5) )
	vrp.push_back(*iv);
      //hmsi->Fill((*iv)->size()*1.);
      else
	{
	  if ((*iv)->Z()<zmin) zmin=(*iv)->Z();
	  if ((*iv)->Z()>zmax) zmax=(*iv)->Z();
	}
    }
  recoTrack::combinePoint(vrp,_geo,_vtk);
  //std::cout<<_vtk.size()<<" tracks found "<<std::endl;
  uint32_t ntk=_vtk.size();
    
    
  bool bif=false;
  if (chbif>0)
    {
      std::map<uint32_t,std::bitset<64> >::iterator ick=_tcount.find(iseed-5);
      if (ick!=_tcount.end())
	if (ick->second[57]!=0) bif=true;
	  
      ick=_tcount.find(iseed-6);
      if (ick!=_tcount.end()) 
	if (ick->second[57]!=0) bif=true;
      ick=_tcount.find(iseed-7); 
      if (ick!=_tcount.end())
	if (ick->second[57]!=0) bif=true;
    }
  bif=bif || (chbif==0);
    
  // if (bif)
  //	this->draw(point);
  // Select shower and muons
  if (hitsparasic>1.5 && (cp[3])/cp[5]>1E-2 &&  point.size()*1./clusters.size()>3. && nbad<5 && bx>50 && bx< 70 && by>45 && by<60 && ntk>0 && zmin>2 && zmax<130.)
    {nshower++;_numberOfShower++;

      printf("%d,%d,%d,%lld,%d,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%d\n",_run,_event,_gtc,_bxId,iseed,_t,_tspill,im->second.count(),ifm->second.size(),nhit,nh0,nh1,nh2,clusters.size(),ntk,zmin,zmax,bif);
    
      double n=nhit;
      double n2=n*n;
      double fe =(nh0*(par[0]+par[1]*n+par[2]*n2)+
		  nh1*(par[3]+par[4]*n+par[5]*n2)+
		  nh2*(par[6]+par[7]*n+par[8]*n2));

    
      //printf("%d %d %f %f %f \n",ngood,nbad,zmin,zmax,f);
    
      /** 733754 80 GeV */
      //  nh0=nh0+int((_t-_tspill-iseed*2E-7)*5.5);
      //nh1=nh1+int((_t-_tspill-iseed*2E-7)*3.1);
      //nh2=nh2+int((_t-_tspill-iseed*2E-7)*3.);

      //std::cout<<nhit<<" "<<nh0+nh1+nh2<<std::endl;

      /** 733724  40  GeV
	  nh0=nh0+int((_t-_tspill-iseed*2E-7)*1.2);
	  nh1=nh1+int((_t-_tspill-iseed*2E-7)*1.28);
	  nh2=nh2+int((_t-_tspill-iseed*2E-7)*1.26);
      */
      hr0->Fill(_t-_tspill-iseed*2E-7,nh0);
      hr1->Fill(_t-_tspill-iseed*2E-7,nh1);
      hr2->Fill(_t-_tspill-iseed*2E-7,nh2);
      double en=(bestX1+bestX2*nhit+bestX3*nhit*nhit)*nh0;
 
      en=en+(bestY1+bestY2*nhit+bestY3*nhit*nhit)*nh1;
      en=en+(bestZ1+bestZ2*nhit+bestZ3*nhit*nhit)*nh2;
      en=en*1.;
      for (std::vector<rCluster<recoPoint>*>::iterator iv=clusters.begin();iv!=clusters.end();iv++) 
	{
	  hssi->Fill((*iv)->size()*1.);
	}
      //double en=0.057*(nh0*0.2+(1+nhit/15000.)*nh1+5*(1-nhit/15000)*nh2);
      //  std::cout<<" CANDIDATE seed :"<<iseed<<" -> chambers: "<<im->second.count()<<" -> frames: "<<ifm->second.size()<<" -> Hits: "<<nhit<<" -> ratio:"<<hitsparasic<<" "<<nh0<<" "<<nh1<<" "<<nh2<<" "<<en<<" " <<ngood<<" "<<nbad<<" tk "<<ntk<<std::endl;
      hnhits->Fill(nh0+nh1+nh2);
      hen->Fill(fe);
      if (bif) henb->Fill(fe);
      //this->draw(point);

 
    }
  else
    {nmuon++;_numberOfMuon++;
      hnhitm->Fill(nhit);
      for (std::vector<rCluster<recoPoint>*>::iterator iv=clusters.begin();iv!=clusters.end();iv++) 
	{
  
	  hmsi->Fill((*iv)->size()*1.);
	}
      //if (ntk==0 && point.size()>30)
      //  this->draw(point);
    }

  for (std::vector<recoPoint*>::iterator it=point.begin();it!=point.end();it++) delete (*it);
  for (std::vector<rCluster<recoPoint>*>::iterator iv=clusters.begin();iv!=clusters.end();iv++) delete((*iv));
  //std::cout<<" CANDIDATE seed :"<<im->first<<" -> chambers: "<<im->second.count()<<" -> frames: "<<ifm->second.size()<<" -> Hits: "<<nhit<<" -> ratio:"<<hitsparasic<<std::endl;
}

static TCanvas* TCHits=NULL;
void tdcreadbinary::draw(std::vector<recoPoint*> vp)
{
  TH2* hzx=_rh->GetTH2("/BR/hzx");
  TH2* hzy=_rh->GetTH2("/BR/hzy");
  if (hzx==NULL)
    {
      hzx=_rh->BookTH2("/BR/zx",150,0.,150.,100,0.,100.);
      hzy=_rh->BookTH2("/BR/zy",150,0.,150.,100,0.,100.);
    }
  else
    {
      hzx->Reset();
      hzy->Reset();
    }
  for (std::vector<recoPoint*>::iterator it=vp.begin();it!=vp.end();it++)
    {
      hzx->Fill((*it)->Z(),(*it)->X());
      hzy->Fill((*it)->Z(),(*it)->Y());
    }
  if (TCHits==NULL)
    {
      TCHits=new TCanvas("TCHits","tChits1",900,900);
      TCHits->Modified();
      TCHits->Draw();
      TCHits->Divide(1,2);
    }
  TCHits->cd(1);
  hzx->Draw("COLZ");
#define drawtk
#ifdef drawtk
  std::vector<TLine*> vl;
  for (std::vector<recoTrack*>::iterator itk=_vtk.begin();itk!=_vtk.end();itk++)
    {
      ROOT::Math::XYZPoint pmin=(*itk)->extrapolate((*itk)->zmin());
      ROOT::Math::XYZPoint pmax=(*itk)->extrapolate((*itk)->zmax());

      TLine* l = new TLine(pmin.Z(),pmin.X(),pmax.Z(),pmax.X());
      l->SetLineColor(2);
      l->Draw("SAME");
      vl.push_back(l);
      std::cout<<pmin.X()<<" "<<pmax.X()<<std::endl;
    }
#endif
  TCHits->Modified();
  TCHits->cd(2);
  hzy->Draw("COLZ");
  for (std::vector<recoTrack*>::iterator itk=_vtk.begin();itk!=_vtk.end();itk++)
    {
      ROOT::Math::XYZPoint pmin=(*itk)->extrapolate((*itk)->zmin());
      ROOT::Math::XYZPoint pmax=(*itk)->extrapolate((*itk)->zmax());

      TLine* l = new TLine(pmin.Z(),pmin.Y(),pmax.Z(),pmax.Y());
      l->SetLineColor(2);
      l->Draw("SAME");
      vl.push_back(l);
      std::cout<<pmin.X()<<" "<<pmax.X()<<std::endl;
    }
  TCHits->Modified();
  TCHits->Draw();
  TCHits->Update();
  getchar();
  for (std::vector<TLine*>::iterator il=vl.begin();il!=vl.end();il++) delete (*il);
  std::cout<<"fini "<<std::endl;
}

void tdcreadbinary::end()
{

  if (_runType==1)
    {
      for (int mez=1;mez<=2;mez++)
	{
	  std::stringstream sr;
	  
	  sr<<"/run"<<_run<<"/TDC"<<mez<<"/";

	  int ipr=0;
	  for (int ich=0;ich<32;ich++)
	    {

	      if (ich%2==0)
		ipr=ich/2;
	      else
		ipr=31-ich/2;
	      std::stringstream src;
	      src<<sr.str()<<"dac"<<ich;
	      TH1* hdac=_rh->GetTH1(src.str());
	      int ped=31;
	      if (hdac!=NULL)
		{
		  printf("Mezzanine %d Channel %d Mean %f RMS %f \n",mez,ich,hdac->GetMean(),hdac->GetRMS());
		  ped=int(hdac->GetMean());
		  if (hdac->GetRMS()>6)
		    {
		      printf("\t \t ======================>ILL %d \n",ipr);
		      ped-=int(hdac->GetRMS());
		    }
		  // if (mez==1)
		  //   _s1.set6bDac(ipr,ped&0xFF);
		  // else
		  //   _s2.set6bDac(ipr,ped&0xFF);
		}
	      // else
	      // 	if (mez==1)
	      // 	  _s1.set6bDac(ipr,31);
	      // 	else
	      // 	  _s2.set6bDac(ipr,31);

	      if (ped==0)
		{printf("\t \t ======================>DEAD %d \n",ipr);
		  ped=31;
		}
	      printf("\t %d %d \n",ipr,ped);
	    }
	}
      // _s1.toJson();
      // _s1.dumpJson();
      // _s2.toJson();
      // _s2.dumpJson();

    }
  



  std::stringstream sr;
  sr<<"/tmp/tdcb"<<_run<<".root";
  
  _rh->writeHistograms(sr.str());


  
}
#define FW29

void tdcreadbinary::timeAnalysis()
{

  if (_channels.size()>512) return;
  uint32_t run=_run;
  // Analyze
  uint16_t chtrg;
  if (_mezzanine==1)
    chtrg=0x10;
  else
    chtrg=0x1c;
  chtrg=0x18; //24
  chtrg=0x1c; //28
  std::stringstream sr;
  sr<<"/run"<<run<<"/TDC"<<_mezzanine<<"/TimeAnalysis/";
  TH2* hpos=_rh->GetTH2(sr.str()+"Position");
  TH1* hdt=_rh->GetTH1(sr.str()+"DeltaT");
  TH1* hdtau=_rh->GetTH1(sr.str()+"DeltaTU");
  TH1* hdc=_rh->GetTH1(sr.str()+"DeltaCluster");

  TH2* hdtr=_rh->GetTH2(sr.str()+"DeltaTr");
  TH1* hns=_rh->GetTH1(sr.str()+"NChannel");
  TH1* hfin=_rh->GetTH1(sr.str()+"Fine");
  TH1* heff=_rh->GetTH1(sr.str()+"Efficiency");
  TH1* hstrip=_rh->GetTH1(sr.str()+"Strips");
  TH1* hstript=_rh->GetTH1(sr.str()+"Stript");
  TH1* hnstrip=_rh->GetTH1(sr.str()+"NStrips");
  TH1* hxp=_rh->GetTH1(sr.str()+"XP");
  TH1* hti=_rh->GetTH1(sr.str()+"time");
  TH1* hra=_rh->GetTH1(sr.str()+"rate");
  if (hpos==NULL)
    {
      hpos=_rh->BookTH2(sr.str()+"Position",100,0.,20.,300,-300.,300.);
      hdt=_rh->BookTH1(sr.str()+"DeltaT",500,-10.,10.);
      hdtau=_rh->BookTH1(sr.str()+"DeltaTU",500,-10.,10.);
      hdc=_rh->BookTH1(sr.str()+"DeltaCluster",500,-10.,10.);
      hdtr=_rh->BookTH2(sr.str()+"DeltaTr",32,0,32.,3000,-300.,300.);
      hns=_rh->BookTH1(sr.str()+"NChannel",1024,0.,1024.);
      hfin=_rh->BookTH1(sr.str()+"Fine",257,0.,257.);
      heff=_rh->BookTH1(sr.str()+"Efficiency",32,0.,32.);
      hstrip=_rh->BookTH1(sr.str()+"Strips",32,0.,32.);
      hstript=_rh->BookTH1(sr.str()+"Stript",32,0.,32.);
      hnstrip=_rh->BookTH1(sr.str()+"NStrips",32,0.,32.);
      hxp=_rh->BookTH1(sr.str()+"XP",400,0.,10.);
      hti=_rh->BookTH1(sr.str()+"time",400,0.,0.2);
      hra=_rh->BookTH1(sr.str()+"rate",750,0.,2000.);

    }
  hns->Fill(_channels.size()*1.);
  // Firmware 100 ps
#ifdef OLDFIRMWARE
  double shift1[32]={-0.61,2.61,-0.742,-0.665,-0.311,0.71,0.198,0.,24*0};
#else
  // firmware 30 ps 1/2 16 pistes
  //double shift1[32]={0.,-1.152,-0.637,-2.124,-0.880,-0.731,-1.072,0.,24*0};
  // double shift1[32]={32*0};
#ifdef FW24
  double shift1[32]={-3.07,-3.66,-5.441,-4.61,-5.008,-5.491,-5.666,-5.873,
		     -5.558,-6.298,-5.339,-4.911,20*0};
#endif
#ifdef FW29OLD
  double shift1[32]={-3.07,-3.66,-4.719,-4.80,-5.036,-5.299,-5.249,-5.925,
		     -6.491,-6.248,-4.971,-5.319,-5.592,-5.725,18*0};
  double shift2[32]={0,0,6.413,5.292,5.614,7.932,6.453,6.407,3.822,9.33,8.923,3.406,5.729,19*0};
#endif
#ifdef FW29
  double shift1[32]={32*0};
		    
  double shift2[32]={32*0};

#endif
 
#endif
  _nevt++;
  heff->Fill(0.1);
  uint32_t ndec=0;
  //printf("Event %d %d %d \n",_gtc,_mezzanine,_channels.size());
  float ti=0,tmax=0;
  // if (_channels.size()>0)
  //   ti=_channels.begin()->bcid()*2.E-7;
  uint32_t lbcid=0,bcidshift=0;
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++){

    if (it->bcid()<lbcid) bcidshift+=65535;
    lbcid=it->bcid();
    float t=((int) it->bcid()*2E-7)+(bcidshift*2.E-7)-ti;
    if (t>tmax) tmax=t;
    if (it->channel()==chtrg) ndec++;
    //printf("%d %d %x %x %x \n",_gtc,it->channel(),it->coarse(),it->fine(),it->bcid());
  }
  if (tmax==0 && _channels.size()>0)
    {
      for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
	std::cout<<(int) it->channel()<<" "<<it->coarse()*2.5E-9<<" "<<it->bcid()*2E-7<<endl;
      //getchar();
    }
  else
    std::cout<<"T MAX "<<tmax<<std::endl;
  hti->Fill(tmax);
  if (tmax>0) hra->Fill(_channels.size()/tmax);
  if (ndec>1) return;
  _trigger=_channels.end();

  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
    if (it->channel()==chtrg) {_trigger=it;it->setUsed(true);break;}
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
    {
      hstrip->Fill(it->channel()*1.);
      for (int j=0;j<8;j++)
	{
	  if (_mezzanine==2 )
	    shift1[j]=0;
	  uint8_t fi=it->fine();
	  //printf("%d %d %d %d \n",fi,j,(fi>>j),((fi>>j)&0x1)); 
	  if (((fi>>j)&0x1) !=0)
	    {
	      //  printf("--->ok\n");
	      hfin->Fill(j+0.1);
	    }
	}

      // Mask channel 0/1 fro FW29 and mezzanine 1
#ifdef FW29
      //if (_mezzanine==1 && it->channel()/2==0) it->setUsed(true);
#endif
    }
  if (_trigger==_channels.end()) return;
  printf("TDC %d  GTC %d   TRIGGER %x %x \n",_mezzanine,_gtc,_trigger->bcid(),_trigger->coarse()/80);
  //_trigger=_channels.begin();
  _ntrigger++;
  heff->Fill(1.1);
  bool found=false;
  uint32_t ngood=0;
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
    {

      if (it->used()) continue;
      int32_t dbcid=it->bcid()-_trigger->bcid();
      if ((dbcid>-4 && dbcid<4))
	{found=true; ngood++; hstript->Fill(it->channel()*1.);}
      else
	it->setUsed(true);
    }
  if (!found) return;
  _nfound++;
  heff->Fill(2.1);
  // if (ngood>6) return;
  bool bside=false;
  bool sused[32]={32*0};
#ifdef FW29
  //if (_mezzanine==1) sused[0]=true;
#endif
  typedef std::pair<uint8_t,double> tdchit;
  typedef std::vector<tdchit> tdcl;
  std::vector<tdchit> hits;
  hits.clear();
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
    {
      printf("%d %d %f \n",it->channel(),it->bcid(),it->tdcTime());
      if (it->used()) continue;
      hdtr->Fill(it->channel(),it->tdcTime()-_trigger->tdcTime());
      if (it->channel()%2!=0) continue;
     
      if (it->tdcTime()-_trigger->tdcTime()>-135.) continue;
      int str=it->channel()/2;
      if (sused[str]) continue;
      for (std::vector<TdcChannel>::iterator jt=_channels.begin();jt!=_channels.end();jt++)
	{

	  if (it==jt) continue;
	  if (jt->used()) continue;
	  if (jt->channel()!=(it->channel()+1)) continue;
	  double t0=it->tdcTime();
	  double t1=jt->tdcTime();
	  ///if (it->coarse()!=jt->coarse()-1) continue;
	  if (abs(t1-t0)>100) continue;
	  it->setUsed(true);
	  jt->setUsed(true);
	  //if (str==0) std::cout<<"oops"<<(int) it->channel()<<" "<<(int) jt->channel()<<std::endl;

	  sused[str]=true;
	  double tsh=0;
	  if (_mezzanine==1)
	    tsh=shift1[str];
	  else
	    tsh=shift2[str];
	  if (_mezzanine==1 && str>=4&&str<=12)
	    hdt->Fill((t0-t1)-tsh);
	  if (_mezzanine==2 && str>=7&&str<=12)
	    hdt->Fill((t0-t1)-tsh);
	  hits.push_back(tdchit((uint8_t) str&0xFF,(t0-t1)-tsh));
	  std::stringstream s;
	  s<<"hdts"<<str;
	  TH1* hdts=_rh->GetTH1(sr.str()+s.str());
	  if (hdts==NULL)
	    {
	      hdts=_rh->BookTH1(sr.str()+s.str(),500,-20.,20.);
	    }
	  hdts->Fill(t0-t1-tsh);
          
	  double x=str*0.4+0.2;
	  double y=((t0-t1)-tsh)*100./7.;
	  //std::cout<<x<<" "<<y<<std::endl;
	  bside=true;
	  hpos->Fill(x,y);
	}     
    }
  // Loop on hits
  hnstrip->Fill(hits.size()*1.);
  if (hits.size()<=8)
    {
   
      std::vector<tdchit> hitg;
      for (int i=0;i<32;i++)
	{
	  int np=0;
	  for (auto x:hits)
	    if (x.first==i) np++;
	  if (np!=1) continue;
     
	  printf(" strip %d np %d \n",i,np);
	  for (auto x:hits)
	    if (x.first==i)
	      {
		hitg.push_back(x);
		std::stringstream s;
		s<<"hdtu"<<i;
		TH1* hdtu=_rh->GetTH1(sr.str()+s.str());
		if (hdtu==NULL)
		  {
		    hdtu=_rh->BookTH1(sr.str()+s.str(),500,-20.,20.);
		  }
		hdtu->Fill(x.second);
		if (x.first>=4&&x.first<=12)
		  hdtau->Fill(x.second);

	      }
    
	}
      tdcl v;
      for (auto x:hitg)
	{
	  bool found=false;
	  for (auto c:v)
	    {
	      if (abs(c.first-x.first)<2)
		{
		  found=true;
		  v.push_back(x);
		  break;
		}
	    }
	  if (!found && v.size()==0) v.push_back(x);
	}
      float dt=0,xp=0;
      float nt=0;
      for (auto x:v)
	{
	  nt++;
	  dt+=x.second;
	  xp+=x.first;
	  //hxp->Fill(x.first*0.4+0.2);
	}
      if (nt>0&&nt<7) {
	hdc->Fill(dt/nt);
	xp=xp/nt;
       
	hxp->Fill(xp*0.4+0.2);
      }
    }
  if (bside) {_nbside++;heff->Fill(3.1);}
  printf("%d-%d %d  #evt %d #trig %d #found %d  #time %d \n",_run,_event,_gtc,_nevt,_ntrigger,_nfound,_nbside); 
}
void tdcreadbinary::LmAnalysis()
{

  if (_channels.size()>1024) return;
  uint32_t run=_run;
  // Analyze
  uint16_t chtrg;
  if (_mezzanine==1)
    chtrg=0x10;
  else
    chtrg=0x1c;
  chtrg=0x18; //24
  chtrg=16; //28
  std::stringstream sr;
  sr<<"/run"<<run<<"/TDC"<<_difId<<"/LmAnalysis/";
  TH2* hpos=_rh->GetTH2(sr.str()+"Position");
  TH1* hdt=_rh->GetTH1(sr.str()+"DeltaT");
  TH1* hdtau=_rh->GetTH1(sr.str()+"DeltaTU");
  TH1* hdc=_rh->GetTH1(sr.str()+"DeltaCluster");

  TH2* hdtr=_rh->GetTH2(sr.str()+"DeltaTr");
  TH1* hns=_rh->GetTH1(sr.str()+"NChannel");
  TH1* hnst=_rh->GetTH1(sr.str()+"NChannelTrigger");
  TH1* hfin=_rh->GetTH1(sr.str()+"Fine");
  TH1* heff=_rh->GetTH1(sr.str()+"Efficiency");
  TH1* hstrip=_rh->GetTH1(sr.str()+"Strips");
  TH1* hstript=_rh->GetTH1(sr.str()+"Stript");
  TH1* hnstrip=_rh->GetTH1(sr.str()+"NStrips");
  TH1* hxp=_rh->GetTH1(sr.str()+"XP");
  TH1* hti=_rh->GetTH1(sr.str()+"time");
  TH1* hra=_rh->GetTH1(sr.str()+"rate");
  if (hpos==NULL)
    {
      hpos=_rh->BookTH2(sr.str()+"Position",100,0.,20.,300,-300.,300.);
      hdt=_rh->BookTH1(sr.str()+"DeltaT",500,-10.,10.);
      hdtau=_rh->BookTH1(sr.str()+"DeltaTU",500,-10.,10.);
      hdc=_rh->BookTH1(sr.str()+"DeltaCluster",500,-10.,10.);
      hdtr=_rh->BookTH2(sr.str()+"DeltaTr",32,0,32.,400,-1000.,-800.);
      hns=_rh->BookTH1(sr.str()+"NChannel",1024,0.,1024.);
      hnst=_rh->BookTH1(sr.str()+"NChannelTrigger",1024,0.,1024.);
      hfin=_rh->BookTH1(sr.str()+"Fine",257,0.,257.);
      heff=_rh->BookTH1(sr.str()+"Efficiency",32,0.,32.);
      hstrip=_rh->BookTH1(sr.str()+"Strips",32,0.,32.);
      hstript=_rh->BookTH1(sr.str()+"Stript",32,0.,32.);
      hnstrip=_rh->BookTH1(sr.str()+"NStrips",32,0.,32.);
      hxp=_rh->BookTH1(sr.str()+"XP",400,0.,10.);
      hti=_rh->BookTH1(sr.str()+"time",400,0.,0.2);
      hra=_rh->BookTH1(sr.str()+"rate",750,0.,200000.);

    }
  hns->Fill(_channels.size()*1.);
  // Firmware 100 ps
#ifdef OLDFIRMWARE
  double shift1[32]={-0.61,2.61,-0.742,-0.665,-0.311,0.71,0.198,0.,24*0};
#else
  // firmware 30 ps 1/2 16 pistes
  //double shift1[32]={0.,-1.152,-0.637,-2.124,-0.880,-0.731,-1.072,0.,24*0};
  // double shift1[32]={32*0};
#ifdef FW24
  double shift1[32]={-3.07,-3.66,-5.441,-4.61,-5.008,-5.491,-5.666,-5.873,
		     -5.558,-6.298,-5.339,-4.911,20*0};
#endif
#ifdef FW29OLD
  double shift1[32]={-3.07,-3.66,-4.719,-4.80,-5.036,-5.299,-5.249,-5.925,
		     -6.491,-6.248,-4.971,-5.319,-5.592,-5.725,18*0};
  double shift2[32]={0,0,6.413,5.292,5.614,7.932,6.453,6.407,3.822,9.33,8.923,3.406,5.729,19*0};
#endif
#ifdef FW29
  double shift1[32]={32*0};
		    
  double shift2[32]={32*0};

#endif
 
#endif
  _nevt++;
  heff->Fill(0.1);
  uint32_t ndec=0;
  //printf("Event %d %d %d \n",_gtc,_mezzanine,_channels.size());
  float ti=0,tmax=0;
  // if (_channels.size()>0)
  //   ti=_channels.begin()->bcid()*2.E-7;
  uint32_t lbcid=0,bcidshift=0,bcidmax=0;
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++){
    if (it->bcid()>bcidmax) bcidmax=it->bcid();
    // if (it->bcid()<lbcid) {
    //   printf("lbcid %d bcid %d  Shift %d \n",lbcid,it->bcid(),b
    //   bcidshift+=65535;
    // }
    lbcid=it->bcid();
    float t=((int) it->bcid()*2E-7)+(bcidshift*2.E-7)-ti;
    if (t>tmax) tmax=t;
    if (it->channel()==chtrg) ndec++;
    //printf("%d %d %x %x %x \n",_gtc,it->channel(),it->coarse(),it->fine(),it->bcid());
  }
  //printf("BCID max %d Bcidshift %d Tmax %f \n",bcidmax,bcidshift,tmax);
  if (tmax==0 && _channels.size()>0)
    {
      for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
	std::cout<<(int) it->channel()<<" "<<it->coarse()*2.5E-9<<" "<<it->bcid()*2E-7<<endl;
      //getchar();
    }
  
  hti->Fill(tmax);
  if (tmax>0) hra->Fill(_channels.size()/tmax);
  // Accept events with only one trigger
  if (ndec!=1) return;
  // Find the trigger


  std::map<uint32_t,std::vector<TdcChannel> > _trmap;
  _trmap.clear();
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
    {
      if (it->channel()==chtrg)
	{
	it->setUsed(true);
	std::vector<TdcChannel> vc;
	vc.push_back(*it);

	
	for (auto x:_channels)
	  {
	    if (x.used()) continue;
	    if (x.channel() == chtrg) continue;
	    if (x.bcid()>(it->bcid()-3) || x.bcid()<(it->bcid()-6)) continue;
	    //if ((x.bcid()-it->bcid())!=-4) continue;
	    vc.push_back(x);
	    x.setUsed(true);
	  }
	std::pair<uint32_t,std::vector<TdcChannel> > p(it->bcid(),vc);
	_trmap.insert(p);
	}

      }//break;}


  // Now loop on all channel of the event
  for (std::vector<TdcChannel>::iterator it=_channels.begin();it!=_channels.end();it++)
    {
      hstrip->Fill(it->channel()*1.);
    }

  printf("TDC %d  GTC %d   Number %d \n",_difId,_gtc,_trmap.size());
  bool found=false;
   bool bside=false;
 
  for (auto t:_trmap)
    {
      //printf("Trigger %d \n",t.first);
      double trigtime=0,trigbcid=0;
      bool chused[32]={32*0};
      bool sused[32]={32*0};
      for (auto x:t.second)
	if (x.channel()==chtrg) {chused[chtrg]=1;trigtime=x.tdcTime();trigbcid=x.bcid()*200;}
      hnst->Fill(t.second.size()*1.);
      if (t.second.size()>20) continue;
      for (auto x:t.second)
	{
	  printf("Chan %d bcid %d Time %f %f \n",x.channel(),x.bcid(),x.tdcTime(),trigtime);
	  if (x.channel()!=chtrg)
	    found=true;
	  else
	    continue;
	  if (chused[x.channel()]) continue;
	  if (x.channel()%2!=0) continue;
	  if (x.tdcTime()-trigtime>-861) continue;
	  if (x.tdcTime()-trigtime<-900) continue;
	  hdtr->Fill(x.channel()+1.,x.tdcTime()-trigtime);
	  //hdtr->Fill(x.channel(),x.bcid()*200-trigbcid);

	   int str=x.channel()/2;
	   if (sused[str]) continue;
	  for (auto y:t.second)
	    {
	      if (chused[y.channel()]) continue;
	  if (y.channel()!=(x.channel()+1)) continue;
	  double t0=x.tdcTime();
	  double t1=y.tdcTime();
	  ///if (it->coarse()!=jt->coarse()-1) continue;
	  //if (abs(t1-t0)>100) continue;
	  chused[x.channel()]=1;
	  chused[y.channel()]=1;
	  //if (str==0) std::cout<<"oops"<<(int) it->channel()<<" "<<(int) jt->channel()<<std::endl;

	  sused[str]=true;
	  double tsh=0;
	  //if (_mezzanine==1 && str>=4&&str<=12)
	  //  hdt->Fill((t0-t1)-tsh);
	  //if (_mezzanine==2 && str>=7&&str<=12)
	    hdt->Fill((t0-t1)-tsh);
	    //hits.push_back(tdchit((uint8_t) str&0xFF,(t0-t1)-tsh));
	  std::stringstream s;
	  s<<"hdts"<<str;
	  TH1* hdts=_rh->GetTH1(sr.str()+s.str());
	  if (hdts==NULL)
	    {
	      hdts=_rh->BookTH1(sr.str()+s.str(),500,-20.,20.);
	    }
	  std::stringstream sx;
	  sx<<"hdtr"<<(int) x.channel();
	  TH1* hdtrx=_rh->GetTH1(sr.str()+sx.str());
	  if (hdtrx==NULL)
	    {
	      hdtrx=_rh->BookTH1(sr.str()+sx.str(),400,-900.,-860.);
	    }
	  hdts->Fill(t0-t1-tsh);
	  hdtrx->Fill(t0-trigtime);
	  std::stringstream sry;
	  sry<<"hdtr"<<(int) y.channel();
	  TH1* hdtry=_rh->GetTH1(sr.str()+sry.str());
	  if (hdtry==NULL)
	    {
	      hdtry=_rh->BookTH1(sr.str()+sry.str(),400,-900.,-860.);
	    }
	  hdtry->Fill(t1-trigtime);
	  double x=str*0.4+0.2;
	  double yp=((t0-t1)-tsh)*100./7.;
	  //std::cout<<x<<" "<<y<<std::endl;
	  bside=true;
	  hpos->Fill(x,yp);
	    }     
	}
    }
  //if (_difId==15 ) getchar();
  //_trigger=_channels.begin();
  _ntrigger++;
  heff->Fill(1.1);
  if (!found) return;
  _nfound++;
  // update efficency
  heff->Fill(2.1);
  if (bside) {_nbside++;heff->Fill(3.1);}
  printf("%d-%d %d  #evt %d #trig %d #found %d  #time %d \n",_run,_event,_gtc,_nevt,_ntrigger,_nfound,_nbside); 
}


#ifdef TESTMAINEXAMPLE
int main()
{
  levbdim::tdcreadbinary bs("/tmp");
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

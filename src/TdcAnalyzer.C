#include "TdcAnalyzer.hh"
#include "jsonGeo.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <bitset>
#include <TCanvas.h>
static TCanvas* TCHits=NULL;

using namespace lydaq;
using namespace lmana;
void lmana::TdcAnalyzer::drawHits(int ch)
{
  
 
  TH2* hpx = rh()->GetTH2("realx");
  TH2* hpy = rh()->GetTH2("realy");
 
  
  if (hpx==NULL)
    {
 
      hpx =rh()->BookTH2("Return",96,0.1,48.1,256,-15.,15.);
      hpy =rh()->BookTH2("Coaxial",96,0.1,48.1,256,-15.,15.);
      hpx->SetMarkerStyle(25);
      hpx->SetMarkerColor(kRed);
      hpy->SetMarkerStyle(25);
      hpy->SetMarkerColor(kBlue);

    }
  else
    {
      if (ch==1)
	hpx->Reset();
      else
	hpy->Reset();

    }

  if (TCHits==NULL)
    {
      TCHits=new TCanvas("TCHits","tChits1",900,900);
      TCHits->Modified();
      TCHits->Draw();
      TCHits->Divide(1,2);
    }
  TCHits->cd(3-ch);
  for (auto x:_strips)
    {
      if (ch==1)
	hpx->Fill(x.xpos()-70,x.ypos());
      else
	{
	  float dx=48-(x.xpos()-70);
	  std::cout<<dx<<std::endl;
	  hpy->Fill(dx,x.ypos());
	}
    }
  if (ch==1)
    hpx->Draw("P");
  else
    {
      hpy->Draw("P");
      TCHits->Modified();
      TCHits->Draw();
      TCHits->Update();
    }
  
}
/**
 * \file evt.C
 * \brief Main analysis class
 * \author L.Mirabito
 * \version 0.1
 * \date 15 septembre 2017
 *
 * Simple example reading SDHCAL H2 Spetember 2017 data
 *
 */
using namespace std;
lmana::TdcAnalyzer::TdcAnalyzer(DCHistogramHandler*r ) : Analyzer(r),_pedestalProcessed(false),_nevt(0),_ntrigger(0),_nfound(0),_nbside(0),_triggerFound(false),_display(false),_noise(false)
{
}
void lmana::TdcAnalyzer::clear()
{_strips.clear();}
void lmana::TdcAnalyzer::processFEB(uint32_t feb,std::vector<lydaq::TdcChannel>& vChannel)
{
  std::cout<<"Process FEB "<<feb<<" "<<jEvent()["runtype"].asUInt()<<std::endl;
  if (jEvent()["runtype"].asUInt()==1) pedestalAnalysis(feb,vChannel);
  if (jEvent()["runtype"].asUInt()==2) scurveAnalysis(feb,vChannel);

}
void lmana::TdcAnalyzer::processChannels(std::vector<lydaq::TdcChannel>& vChannel)
{
  for (int i=0;i<255;i++)
    {
      //printf("FEB %d %d \n",i,geo()->feb(i).id);
      if (geo()->feb(i).isEmpty()) continue;
      processFEB(i,vChannel);
    }
  if (jEvent()["runtype"].asUInt()==0)
    this->multiChambers(vChannel);
}
void lmana::TdcAnalyzer::setInfo(uint32_t dif,uint32_t run,uint32_t ev,uint32_t gt,uint64_t ab,uint16_t trgchan,uint32_t vth,uint32_t dac)
{_dif=dif;
  _run=run;
  _event=ev;
  _gtc=gt;
  _abcid=ab;
  _triggerChannel=trgchan;
  _vthSet=vth;
  _dacSet=dac;
  if (_abcid0==0 || _abcid<_abcid0) _abcid0=_abcid;
  _display=geo()->general()["display"].asUInt()==1;
  _noise=geo()->general()["noise"].asUInt()==1;
  _jEvent["run"]=_run;
  _jEvent["event"]=_event;
  _jEvent["gtc"]=_gtc;
  _jEvent["abcid"]=Json::Value((Json::Value::UInt64)_abcid);
  _jEvent["triggerChannel"]=trgchan;
  _jEvent["vthset"]=vth;
  _jEvent["dacset"]=dac;
  
  //std::cout<<_jEvent<<std::endl;
  //getchar();
}
bool lmana::TdcAnalyzer::noiseStudy(std::vector<lydaq::TdcChannel>& vChannel,std::string subdir)
{
  float ch1_dt[128];
  float ch2_dt[128];

  memset(ch2_dt,0,128*sizeof(float));
  memset(ch1_dt,0,128*sizeof(float));

#define RUN743134
#ifdef RUN743065
  float alg1[45]={0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -5.951, -5.836, -5.586, 0.000, -5.932, -5.518, -5.538, -5.583, 0.000, -6.518, -6.195, -5.831, -5.572, -6.000, -6.076, -5.975, -6.086, -5.966, -5.547, -7.015, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};

  for (int i=0;i<45;i++) {
    ch1_dt[72+i]=alg1[i];}
  float alg2[49]={0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.970, 1.697, 0.485, 0.729, -0.419, -1.690, -2.601, -1.783, -1.922, -2.199, -2.046, -2.089, -2.064, -2.621, -4.695, -4.766, -1.886, -2.313, -2.396, -1.830, -2.271, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
  for (int i=0;i<49;i++) {
    ch2_dt[72+i]=alg2[i];}
#endif
#ifdef RUN743079
  float alg1[45]={0.000, 0.000, -6.572, -5.284, -4.974, -5.085, 0.000, -4.329, -4.858, -5.385, -4.898, -6.042, -6.749, -6.021, -6.035, -5.759, 0.000, -5.863, -5.516, -5.449, -5.467, 0.000, -6.745, -6.107, -5.797, -5.555, -5.856, -5.925, -5.829, -6.056, -5.897, -5.669, -7.332, -4.162, -3.824, -4.705, -4.881, 0.000, -5.893, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
  float alg2[49]={0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 7.167, 1.371, 1.293, 1.741, 0.559, 0.756, -0.408, -1.673, -2.624, -1.695, -1.853, -2.149, -1.987, -1.961, -1.817, -2.609, -4.628, -5.021, -2.207, -2.201, -2.281, -1.775, -2.401, -2.072, 0.000, -4.127, -4.165, -3.657, -4.816, -4.676, -5.599, -6.948, -7.514, -7.008, -6.912, -6.900, 0.000};
  for (int i=0;i<45;i++) {   ch1_dt[72+i]=alg1[i];}
  for (int i=0;i<49;i++) {   ch2_dt[72+i]=alg2[i];}
#endif
#ifdef RUN743134
  float alg1[45]={0.000, 0.000, -6.161, -5.713, -5.275, -4.911, 0.000, -4.309, -4.803, -5.357, -4.827, -5.955, -6.649, -5.904, -5.986, -5.782, 0.000, -5.835, -5.494, -5.430, -5.425, -4.701, -6.648, -6.014, -5.770, -5.537, -5.851, -5.887, -5.830, -6.019, -5.924, -5.680, -7.326, -3.991, -3.796, -4.641, -4.825, 0.000, -5.961, -5.593, -6.831, -6.694, -5.650, 0.000, 0.000};

  float alg2[49]={0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 7.251, 1.260, 1.334, 1.705, 0.474, 0.727, -0.513, -1.777, -2.623, -1.778, -1.943, -2.188, -2.119, -1.867, -2.075, -2.640, -4.647, -5.065, -2.271, -2.221, -2.303, -1.887, -2.338, -2.136, -29.730, -4.138, -4.045, -3.788, -4.752, -4.816, -5.459, -6.912, -7.614, -7.012, -7.248, -7.764, 0.000};
  for (int i=0;i<45;i++) {   ch1_dt[72+i]=alg1[i];}
  for (int i=0;i<49;i++) {   ch2_dt[72+i]=alg2[i];}
#endif
  uint32_t triggerChannel=0;
  float dtmin=-615,dtmax=-585;
  bool noisy=false;
  //dtmin+=100;dtmax+=100;
  for (uint32_t chamber=1;chamber<=2;chamber++)
    {
      _strips.clear();
      std::vector<TdcChannel*> c_strip[128];
      for (int i=0;i<128;i++) c_strip[i].clear();
      float maxtime=0,mttime=0;uint32_t nch=0,ntrg=0;
      float ttime[24];
      memset(ttime,0,24*sizeof(float));
      for (auto x=vChannel.begin();x!=vChannel.end();x++)
	{
	  if (geo()->feb(x->feb()).chamber!=chamber) continue;
	  if (x->channel()!=triggerChannel) continue;
	  ttime[x->feb()]=x->tdcTime();
	  mttime+=x->tdcTime();
	  ntrg++;
	}
      // for (int i=1;i<24;i++)
      // 	printf("%d -> %f \n",i,ttime[i]);
      //if (chamber==1 && ntrg!=4) return true;
      mttime=mttime/ntrg;
      if (mttime<abs(dtmin)+50) return true;
      std::bitset<24> stfeb(0);
      for (auto x=vChannel.begin();x!=vChannel.end();x++)
	{
	  if (geo()->feb(x->feb()).chamber!=chamber) continue;
	  if (x->channel()==triggerChannel) continue;
	  // if (chamber==1 && x->side(geo()->feb(x->feb()))==0)
	  //   {dtmin=-614;dtmax=-594;}
	  // if (chamber==1 && x->side(geo()->feb(x->feb()))==1)
	  //   {dtmin=-606;dtmax=-586;}
	  // if (chamber==2 )
	  //   {
	  //     dtmin=-615.; dtmax=-585;
	  //   }
	  //dtmin=dtm[x->feb()][ x->side(geo()->feb(x->feb()))]-10.;
	  //dtmax=dtm[x->feb()][ x->side(geo()->feb(x->feb()))]+10.;
	  dtmin=geo()->feb(x->feb()).dt[x->side(geo()->feb(x->feb()))]-10.;
	  dtmax=geo()->feb(x->feb()).dt[x->side(geo()->feb(x->feb()))]+10.;

	  if (x->tdcTime()-ttime[x->feb()]>dtmin-200 && x->tdcTime()-ttime[x->feb()]<dtmax-200)
	    {
	      stfeb.set(x->feb(),1);
	      //std::cout<<"================================> bit set \n";
	      std::stringstream sraw;
	      sraw<<"/run"<<_run<<"/"<<subdir<<"/Chamber"<<chamber<<"/Raw/";
	      TH1* hchan=rh()->GetTH1(sraw.str()+"Channels");
	      TH1* hstrips=rh()->GetTH1(sraw.str()+"Strips");
	      if (hchan==NULL)
		{
		  hchan=rh()->BookTH1(sraw.str()+"Channels",24*16,0.,24.*16);
		  hstrips=rh()->BookTH1(sraw.str()+"Strips",96,0,96);
		}
	      hchan->Fill(x->feb()*24+x->channel());
	      hstrips->Fill( x->side(geo()->feb(x->feb()))*48+x->detectorStrip(geo()->feb(x->feb()))-72);
	    }
	  // printf("%f %f \n",dtmin,dtmax);
	  // getchar();
	  // Book and fill time to trigger
	  std::stringstream src;
	  src<<"/run"<<_run<<"/"<<subdir<<"/Chamber"<<chamber<<"/FEB/"<<x->feb()<<"/Side"<<(int) x->side(geo()->feb(x->feb()))<<"/channel"<<(int) x->channel();
	  std::stringstream srcp;
	  srcp<<"/run"<<_run<<"/"<<subdir<<"/Chamber"<<chamber<<"/FEB/"<<x->feb()<<"/Side"<<(int) x->side(geo()->feb(x->feb()))<<"/";
		  
	  TH1* hdt=rh()->GetTH1(src.str());
	  TH2* hdtr=rh()->GetTH2(srcp.str()+"DeltaTrigger");
	  TH1* hdtrp=rh()->GetTH1(srcp.str()+"DTall");
	  //TH1* hfi=rh()->GetTH1(srcp.str()+"Fine");
	  
	  if (hdt==NULL)
	    {
	      
	      hdt=rh()->BookTH1(src.str(),50,-25,25);
	      hdtr=rh()->BookTH2(srcp.str()+"DeltaTrigger",4000,-2000.,1500.,32,0.,32.);
	      hdtrp=rh()->BookTH1(srcp.str()+"DTall",2000,-1000.,0.);
	      //hfi=rh()->BookTH1(srcp.str()+"Fine",100,-2.5,2.5);

	    }
	  // hfi->Fill(x->fine()/256.0*TDC_COARSE_TIME);
	  // printf("%d %d %d %d %f \n",x->feb(),x->channel(),x->coarse(),x->fine(),x->tdcTime());
	  // getchar();
	  float dt=geo()->feb(x->feb()).dtc[x->channel()];
	  hdtr->Fill(x->tdcTime()-ttime[x->feb()]-dt,x->channel());
	  hdtrp->Fill(x->tdcTime()-ttime[x->feb()]-dt);
	  //printf("%f %f %f \n",x->tdcTime(),ttime[x->feb()],x->tdcTime()-ttime[x->feb()]);
	  if (_noise)
	    {
	      dtmin-=200;
	      dtmax-=200;
	    }

	  if (x->tdcTime()>maxtime) maxtime=x->tdcTime();
	  if (x->tdcTime()-ttime[x->feb()]<dtmin) continue;
	  if (x->tdcTime()-ttime[x->feb()]>dtmax) continue;
	  hdt->Fill(x->tdcTime()-ttime[x->feb()]-
		    geo()->feb(x->feb()).dt[x->side(geo()->feb(x->feb()))]);	  
	  //dtm[x->feb()][ 
	  c_strip[x->detectorStrip(geo()->feb(x->feb()))].push_back(&(*x));

	  	      std::stringstream sraw;
	      sraw<<"/run"<<_run<<"/"<<subdir<<"/Chamber"<<chamber<<"/Raw/";
	      TH1* hchan=rh()->GetTH1(sraw.str()+"SelChannels");
	      TH1* hstrips=rh()->GetTH1(sraw.str()+"SelStrips");
	      if (hchan==NULL)
		{
		  hchan=rh()->BookTH1(sraw.str()+"SelChannels",24*16,0.,24.*16);
		  hstrips=rh()->BookTH1(sraw.str()+"SelStrips",96,0,96);
		}
	      hchan->Fill(x->feb()*24+x->channel());
	      hstrips->Fill( x->side(geo()->feb(x->feb()))*48+x->detectorStrip(geo()->feb(x->feb()))-72);



	  
	  nch++;
	}

      std::stringstream srcc;
      srcc<<"/run"<<_run<<"/"<<subdir<<"/Chamber"<<chamber<<"/";
      
      TH1* hfr=rh()->GetTH1(srcc.str()+"FebCount");
      TH1* hfrs=rh()->GetTH1(srcc.str()+"FebCountSel");

      
      if (hfr==NULL)
	{
	      
	  hfr=rh()->BookTH1(srcc.str()+"FebCount",25,0.,25.);
	  hfrs=rh()->BookTH1(srcc.str()+"FebCountSel",25,0.,25.);

	}
      hfr->Fill(24.);
      hfrs->Fill(24.);
      for (int i=0;i<24;i++)
	if (stfeb[i]>0) hfr->Fill(i*1.);

      
      maxtime=maxtime*1E-9;

      fprintf(stderr," Maxtime %d %f %d %f \n",chamber,maxtime,nch,nch/maxtime/6500);
      //getchar();
      bool dostop=false;int nstrip=0;
      uint16_t febc[24];
      memset(febc,0,48);
      std::bitset<49> stb(0);
      for (int i=0;i<128;i++)
	{
	  if (c_strip[i].size()>0)
	    {
	      //fprintf(stderr,"Chamber %d Strip %d # %d \n",chamber,i,c_strip[i].size());
	      nstrip++;

	      stb.set(i-70,1);
	    }
	  if (c_strip[i].size()>2) dostop=true;

	  if (c_strip[i].size()==2)
	    {
	      double t0=-1,t1=-1;
	      for (auto x:c_strip[i])
		{

		  //fprintf(stderr,"\t %d %d %f %f \n",x->channel(), x->side(geo()->feb(x->feb())),x->tdcTime(),x->tdcTime()-ttime);
		  double dt=geo()->feb(x->feb()).dtc[x->channel()];
		  if (t0<0 &&  x->side(geo()->feb(x->feb()))==0)
		    {
		      t0=x->tdcTime()-dt;

		      //printf("T0 %d %d %d %d %f %f dt=%f \n",x->feb(),x->channel(),x->coarse(),x->fine(),x->tdcTime(),t0,dt);
		    }
		  if (t1<0 &&  x->side(geo()->feb(x->feb()))==1)
		    {
		      t1=x->tdcTime()-dt;
		      //printf("T1 %d %d %d %d %f %f \n",x->feb(),x->channel(),x->coarse(),x->fine(),x->tdcTime(),t1);
		    }
		  if(t0>0 && t1>0 )
		    {
		      febc[x->feb()]++;
		      //std::cout<<x->feb()<<" FEBC "<< febc[x->feb()]<<std::endl;
		      if (geo()->feb(x->feb()).polarity==-1)
			{
			  double tt=t1;
			  t1=t0;
			  t0=tt;
			}
		    
		      //lmana::TdcStrip ts(geo()->feb(x->feb()).chamber,x->feb(),x->detectorStrip(geo()->feb(x->feb())),t0,t1,geo()->feb(x->feb()).timePedestal[x->detectorStrip( geo()->feb(x->feb()))-70]);
		      if (chamber==1)
			{
			  lmana::TdcStrip ts(geo()->feb(x->feb()).chamber,x->feb(),x->detectorStrip(geo()->feb(x->feb())),t0,t1,ch1_dt[x->detectorStrip( geo()->feb(x->feb()))+1]);
			  _strips.push_back(ts);
			}
		      else
			{
			  lmana::TdcStrip ts(geo()->feb(x->feb()).chamber,x->feb(),x->detectorStrip(geo()->feb(x->feb())),t0,t1,ch2_dt[x->detectorStrip( geo()->feb(x->feb()))+1]);
			  _strips.push_back(ts);
			}

		    }
		}



	    
	    }
	}
      for (int i=0;i<24;i++)
	if (febc[i]>0) hfrs->Fill(i*1.);

      if (dostop) return true;
      if (stb.count()>24) return true;
      noisy=(stb.count()>24);
      //for (int i=0;i<24;i++)
      // if (febc[i]>=10) return true;
      //std::cout<<stb<<std::endl;
      std::vector<lmana::TdcCluster> vclus;
      vclus.clear();
      float step=2.;
      if (chamber==1) step=4.;
      for (auto it=_strips.begin();it!=_strips.end();it++)
	{
	  //fprintf(stderr,"%d %d %f %f %f %f \n",it->chamber(),it->strip(),it->xpos(),it->ypos(),it->t0(),it->t1());
	  if (it->chamber()!=chamber) continue;
	  //	      if (it->ypos()<-10 || it->ypos()>-0.2) continue;
	  bool found=false;
	  for (auto ic=vclus.begin();ic!=vclus.end();ic++)
	    {
	      if (ic->isAdjacent((*it),step))
		{
		  ic->addStrip((*it));
		  found=true;
		  break;
		}
	    }
	  if(!found)
	    {
	      lmana::TdcCluster c;
	      c.addStrip((*it));
	      vclus.push_back(c);
	    }
	}

      // Merge adjacent cluster
      bool merged=false;
      //printf("vclus size %d \n",vclus.size());
      for (auto it=vclus.begin();it!=vclus.end();it++)
	{
	  for (auto jt=it+1;jt!=vclus.end();)
	    {
	      bool adj=false;
	      for (int i=0;i<jt->size();i++)
		{
		  if (it->isAdjacent(jt->strip(i),step))
		    {adj=true;break;}
		}
	      if (adj)
		{
		  merged=true;
		  printf("Merigng cluster \n");
		  for (int i=0;i<jt->size();i++)
		    {
		      it->addStrip(jt->strip(i));
		    }
		  vclus.erase(jt);
		}
	      else
		++jt;
	    }
	}
      //printf("vclus size after %d \n",vclus.size());
      //if (merged)
      //	getchar();
      std::stringstream src;
      src<<"/run"<<_run<<"/"<<subdir<<"/Chamber"<<chamber<<"/ClusterNew/";
		  
      TH2* hposs=rh()->GetTH2(src.str()+"XYStrip");
      TH2* hposc=rh()->GetTH2(src.str()+"XY");
      TH2* hposc1=rh()->GetTH2(src.str()+"XY1");
      TH2* hposcm=rh()->GetTH2(src.str()+"XYMore");
      TH2* hposcma=rh()->GetTH2(src.str()+"XYMax");
      TH2* hposx=rh()->GetTH2(src.str()+"XYX");
      TH1* hncl=rh()->GetTH1(src.str()+"Clusters");
      TH1* hmulc=rh()->GetTH1(src.str()+"ClusterSize");
      TH1* hmulc1=rh()->GetTH1(src.str()+"ClusterSize1");
      TH1* hns=rh()->GetTH1(src.str()+"nstrip");
      TH1* hns2=rh()->GetTH1(src.str()+"nstrip2");
      TH1* htoa=rh()->GetTH1(src.str()+"TOA");

      if (hposc==NULL)
	{
	      
	  hposc=rh()->BookTH2(src.str()+"XY",128,0.,128.,256,-15.,15.);
	  hposs=rh()->BookTH2(src.str()+"XYStrip",128,0.,128.,256,-15.,15.);
	  hposc1=rh()->BookTH2(src.str()+"XY1",128,0.,128.,256,-15.,15.);
	  hposcm=rh()->BookTH2(src.str()+"XYMore",128,0.,128.,256,-15.,15.);
	  hposcma=rh()->BookTH2(src.str()+"XYMax",128,0.,128.,3000,-15.,15.);
	  hposx=rh()->BookTH2(src.str()+"XYX",128,0.,128.,600,-160.,160.);
	  hncl=rh()->BookTH1(src.str()+"Clusters",32,0.,32.);
	  hmulc=rh()->BookTH1(src.str()+"ClusterSize",32,0.,32.);
	  hmulc1=rh()->BookTH1(src.str()+"ClusterSize1",32,0.,32.);
	  hns=rh()->BookTH1(src.str()+"nstrip",48,0.,48.);
	  hns2=rh()->BookTH1(src.str()+"nstrip2",48,0.,48.);
	  htoa=rh()->BookTH1(src.str()+"TOA",100,-570.,-600.);

	}
      hns->Fill(nstrip);
      hns2->Fill(_strips.size());
      printf(" ===> %d  %d strips , Number of clusters %d \n",chamber,_strips.size(),vclus.size());
      //if (vclus.size()>0)
      hncl->Fill(vclus.size()*1.);
      uint32_t maxs=0;
      for (auto x:vclus)
	{
	  // if (vclus.size()>1)
	  //   {
	  if (_display)
	    {
	      fprintf(stderr,"\t %f %f %d \n",x.X(),x.Y(),x.size());
	      for (int i=0;i<x.size();i++)
		{
		  fprintf(stderr,"\t \t  %5.1f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n",x.strip(i).xpos(),x.strip(i).ypos(),x.strip(i).shift(),x.strip(i).t0(),x.strip(i).t1(),(x.strip(i).t0()+x.strip(i).t1())/2.-ttime[x.strip(i).dif()],ttime[x.strip(i).dif()]);
		}
	    }
	  //   }
	  if (vclus.size()==1)
	    for (int i=0;i<x.size();i++)
	      {
		htoa->Fill(
			   (x.strip(i).t0()+x.strip(i).t1())/2.-ttime[x.strip(i).dif()]);
	      }

	  if (x.size()>30) continue;
	  hposc->Fill(x.X(),x.Y());
	  if (vclus.size()==1)
	    {hposc1->Fill(x.X(),x.Y());
	      hmulc1->Fill(x.size()*1.);
	    }
	  
	  
	  hmulc->Fill(x.size()*1.);
	  //printf("%d %f %f \n",x.size(),x.X(),x.Y());
	  if (x.size()>maxs) maxs=x.size();
	}
      //if (vclus.size()>1) getchar();
      int nc=0;
      for (auto x:vclus)
	{
	  nc++;
	  if (x.size()>16) continue;
	  if (x.size()==maxs)
	    {
	      hposcma->Fill(x.X(),x.Y());
	      float L=160.;
	      float v=160./8.7;
	      float xl=(L-x.Y()*v)/2.-L/2.;
	      hposx->Fill(x.X(),xl);

	      for (int i=0;i<x.size();i++)
		{
		  hposs->Fill(x.strip(i).xpos(),x.strip(i).ypos());
		  std::stringstream srcs;
		  srcs<<src.str()<<"align/strip"<<int(x.strip(i).xpos());
		  //std:cout<<srcs.str()<<std::endl;
		  TH1* hdts=rh()->GetTH1(srcs.str());
		  if (hdts==NULL)
		    hdts=rh()->BookTH1(srcs.str(),100,-20,20.);
		  hdts->Fill(x.strip(i).ypos());
		}
	      break;
	    }
	}
      int ncp=0;
      for (auto x:vclus)
	{
	  ncp++;
	  if (x.size()>16) continue;

	  if (ncp==nc) continue;
	  hposcm->Fill(x.X(),x.Y());  
	}
      if (_display)
	{
	  this->drawHits(chamber);
	  if (chamber==2 ) getchar();
	}

    }
  return false;
}

void lmana::TdcAnalyzer::multiChambers(std::vector<lydaq::TdcChannel>& vChannel)
{
  _noise=true;
  this->noiseStudy(vChannel,"OffTime");
  _noise=false;
  if (this->noiseStudy(vChannel,"InTime")) return;
  //if (_noise) return;
  return;
  std::stringstream sr;
  

  uint32_t ndifread=7;
  uint32_t triggerChannel=0;
  
  std::bitset<16> btrg(0);uint32_t ntrg=0,atbcid=0;
  for (auto x:vChannel)
    {
      if (x.channel()==triggerChannel) {
	//printf("Trigger found %d  %d %f\n",x.feb(),x.bcid(),x.tdcTime());
	btrg.set(x.feb(),1);
	atbcid=x.bcid();
	ntrg++;
      }
    }
  //getchar();
  _triggerFound=(btrg.count()%ndifread==0 && ntrg>0);
  // if (btrg.count()==ndifread) _ntrigger++;
  // if (ntrg!=ndifread) return;
  if (!_triggerFound)
    {            return;}
  if (atbcid<2) return;

  for (uint32_t chamber=1;chamber<=2;chamber++)
    {
      _strips.clear();
      double dtmin,dtmax,dtmean;
      std::stringstream sr;
      sr.clear();
      sr.str("");

      sr<<"/run"<<_run<<"/Chamber"<<chamber<<"/";
      TH2* hdtr=rh()->GetTH2(sr.str()+"DeltaTrigger");
      TH2* hdtr0=rh()->GetTH2(sr.str()+"DeltaTrigger0");
      TH2* hdtr1=rh()->GetTH2(sr.str()+"DeltaTrigger1");
      TH2* hdtrt=rh()->GetTH2(sr.str()+"DeltaTriggerSel");
      
      TH2* hdtrt0=rh()->GetTH2(sr.str()+"DeltaTriggerSel0");
      TH2* hdtrt1=rh()->GetTH2(sr.str()+"DeltaTriggerSel1");
      TH1* hnstrip=rh()->GetTH1(sr.str()+"NSTRIP");
      TH1* heff=rh()->GetTH1(sr.str()+"Efficiency");
      TH1* hbp2=rh()->GetTH1(sr.str()+"BeamProfile");
      TH1* hti=rh()->GetTH1(sr.str()+"MaxTime");
      TH1* htti=rh()->GetTH1(sr.str()+"TriggerTime");
      TH1* hnch=rh()->GetTH1(sr.str()+"Channels");
      TH1* hrate=rh()->GetTH1(sr.str()+"Rate");
	    
      TH2* hpos=rh()->GetTH2(sr.str()+"DeltaTvsStrip");
      if (hdtr==NULL)
	{
	  hdtr=rh()->BookTH2(sr.str()+"DeltaTrigger",4000,-2000.,1500.,48,71.,119.);
	  hdtr0=rh()->BookTH2(sr.str()+"DeltaTrigger0",4000,-2000.,1500.,48,71.,119.);
	  hdtr1=rh()->BookTH2(sr.str()+"DeltaTrigger1",4000,-2000.,1500.,48,71.,119.);
	  hdtrt=rh()->BookTH2(sr.str()+"DeltaTriggerSel",4000,dtmin,dtmax,48,71.,119.);
	  hdtrt0=rh()->BookTH2(sr.str()+"DeltaTriggerSel0",4000,dtmin,dtmax,48,71.,119.);
	  hdtrt1=rh()->BookTH2(sr.str()+"DeltaTriggerSel1",4000,dtmin,dtmax,48,71.,119.);
	  hpos=rh()->BookTH2(sr.str()+"DeltaTvsStrip",3000,-30.,30.,48,71.,119.);
	  hnstrip=rh()->BookTH1(sr.str()+"NSTRIP",24,-0.1,23.9);
	  heff=rh()->BookTH1(sr.str()+"Efficiency",10,-0.1,9.9);
	  hbp2=rh()->BookTH1(sr.str()+"BeamProfile",128,0.1,128.1);
	  hti=rh()->BookTH1(sr.str()+"MaxTime",200000,0.,0.25);
	  htti=rh()->BookTH1(sr.str()+"TriggerTime",20000,0.,0.25);
	  hnch=rh()->BookTH1(sr.str()+"Channels",4096,-0.1,4095.9);
	  hrate=rh()->BookTH1(sr.str()+"Rate",20000,0.,10000.);

	}

      std::bitset<128> side[2];
      side[0].reset();
      side[1].reset();

      double febbcid[128];
      memset(febbcid,0,128*sizeof(double));
      heff->Fill(1.1);

      uint32_t maxbcid=0,nch=0;
      for (auto x:vChannel)
	{
	  if (geo()->feb(x.feb()).chamber!=chamber) continue;
	  if (x.channel()==triggerChannel)
	    {
	      htti->Fill(x.tdcTime()*1E-9);
	      continue;
	  
	    }
	  

	  if (x.bcid()>maxbcid) maxbcid=x.bcid();
	  nch++;
      
	}
      if (nch>0)
	{hti->Fill(maxbcid*2e-7);hrate->Fill(nch*1./(maxbcid*2E-7)/7250/2.);
	  std::cout<<nch*1./(maxbcid*2E-7)<<std::endl;}
      hnch->Fill(nch*1.);
      //if (maxbcid*2E-7<1E-5) continue;
      heff->Fill(2.1);

      for (int idif=1;idif<=24;idif++)
	{
	  //      std::cout<<"idif "<<idif<<" "<<geo()->feb(idif).id<<" "<<geo()->feb(idif).chamber<<std::endl;
	  //getchar();
	  if (geo()->feb(idif).id!=idif) continue;
	  if (geo()->feb(idif).chamber!=chamber) continue;

	  std::stringstream srd;
	  srd.clear();
	  srd.str("");
	  srd<<sr.str()<<"FPGATDC"<<idif<<"/";
	  TH1* hbp2d1=rh()->GetTH1(srd.str()+"BeamProfile1");
	  TH1* hbp2d0=rh()->GetTH1(srd.str()+"BeamProfile0");
	  if (hbp2d0==NULL)
	    {
	      hbp2d0=rh()->BookTH1(srd.str()+"BeamProfile0",128,0.1,128.1);
	      hbp2d1=rh()->BookTH1(srd.str()+"BeamProfile1",128,0.1,128.1);
	      std::cout<<"Booking "<<srd.str()+"BeamProfile0"<<std::endl;
	      std::cout<<"Booking "<<srd.str()+"BeamProfile1"<<std::endl;
	      //getchar();
	    }




	  
	  double tbcid=0;
	  dtmin=geo()->feb(idif).triggerMin,dtmax=geo()->feb(idif).triggerMax,dtmean=geo()->feb(idif).triggerMean;uint32_t nchfeb=0;
	  for (auto x:vChannel)
	    {
	      if (x.feb()!=idif) continue;
	      nchfeb++;
	      if (x.channel()!=triggerChannel) continue;
	  
	      tbcid=x.tdcTime();
	      //printf("TRIGGER FOUND %f \n",tbcid);
	      //getchar();
	      break;
	    }
	  //printf("FEB %d %d \n",idif,nchfeb);
	  febbcid[idif]=tbcid;
	  for (auto x:vChannel)
	    {
	      if (x.feb()!=idif) continue;
	      if (x.channel()==triggerChannel) continue;
	      //printf("strip %f %f %d \n",x.tdcTime(),tbcid,x.strip());
	      hdtr->Fill(x.tdcTime()-tbcid,x.detectorStrip(geo()->feb(idif)));
	      if (x.side(geo()->feb(idif)))
		hdtr1->Fill(x.tdcTime()-tbcid,x.detectorStrip(geo()->feb(idif)));
	      else
		hdtr0->Fill(x.tdcTime()-tbcid,x.detectorStrip(geo()->feb(idif)));
	    }
	  // Now fill those in good range
	  for (auto x:vChannel)
	    {
	      if (x.channel()==triggerChannel) continue;
	      if (x.feb()!=idif) continue;
	      if ((x.tdcTime()-tbcid)<dtmin) continue;
	      if ((x.tdcTime()-tbcid)>dtmax) continue;

	      
	      if (chamber==112)
		{
		  printf("strip time  %f  bcid %f strip %d detStrip %d channel %d feb %d \n",x.tdcTime(),tbcid,x.strip(geo()->feb(idif)),x.detectorStrip(geo()->feb(idif)),x.channel(),x.feb());
		  geo()->feb(idif).dump();
		  getchar();
		}
	      hdtrt->Fill(x.tdcTime()-tbcid,x.detectorStrip(geo()->feb(idif)));
	      side[x.side(geo()->feb(idif))].set(x.detectorStrip(geo()->feb(idif)),1);
	      // printf("FEB %d bcid %d \n",x.feb(),x.bcid()); 
	      if (x.side(geo()->feb(idif)))
		{
		  hdtrt1->Fill(x.tdcTime()-tbcid,x.detectorStrip(geo()->feb(idif)));
		  hbp2d1->Fill(x.detectorStrip(geo()->feb(idif)));

		}
	      else
		{
		  hdtrt0->Fill(x.tdcTime()-tbcid,x.detectorStrip(geo()->feb(idif)));
		  hbp2d0->Fill(x.detectorStrip(geo()->feb(idif)));
		}
	    }
  
	}

      //std::cout<<side[0]<<std::endl;
      // std::cout<<side[1]<<std::endl;
      if (side[0].count()!=0 || side[1].count()!=0)
	heff->Fill(3.1);
      //else
      //  getchar();
      uint32_t nstrip=0;
      for (int i=0;i<128;i++)
	if (side[0][i]&&side[1][i])
	  {
	    nstrip++;
	    hbp2->Fill(i*1.);

	  }
      nstrip=0;
      for (int i=0;i<128;i++)
	if (side[0][i]&&side[1][i])
	  {
	    double t0=-1,t1=-1;
	    for (auto x:vChannel)
	      {
		double tbcid=febbcid[x.feb()];
		if (tbcid==0) continue;
		if ((x.tdcTime()-tbcid)<dtmin) continue;
		if ((x.tdcTime()-tbcid)>dtmax) continue;
		if (x.detectorStrip(geo()->feb(x.feb()))!=i) continue;
		if (x.side(geo()->feb(x.feb()))==0 && t0==-1) 	    t0=x.tdcTime();
		if (x.side(geo()->feb(x.feb()))==1 && t1==-1) 	    t1=x.tdcTime();
		//if (x.side(geo()->feb(x.feb()))==0 ) 	    t0=x.tdcTime();
		//if (x.side(geo()->feb(x.feb()))==1 ) 	    t1=x.tdcTime();
	      
	
		if(t0>0 && t1>0 && abs(t0-t1)<11.5)
		  {
		    nstrip++;
		    //printf("%f %f %f\n",t0,t1,t1-t0);
		    hpos->Fill(t1-t0,x.detectorStrip(geo()->feb(x.feb())));
		    std::stringstream s;
		    s<<"Timing/All/hdtpos"<<(int) x.detectorStrip(geo()->feb(x.feb()));
		    TH1* hdts=rh()->GetTH1(sr.str()+s.str());
		    if (hdts==NULL)
		      {
			hdts=rh()->BookTH1(sr.str()+s.str(),300,-25.,25.);
		      }
		    hdts->Fill(t1-t0);
		    s.str("");
		    s.clear();
		    s<<"Timing/Side0/hdtr_"<<(int) x.detectorStrip(geo()->feb(x.feb()));
		    TH1* hdts0=rh()->GetTH1(sr.str()+s.str());
		    if (hdts0==NULL)
		      {
			hdts0=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
		      }
		    hdts0->Fill(t0-tbcid-dtmean);
		    s.str("");
		    s.clear();
		    s<<"Timing/Side1/hdtr_"<<(int) x.detectorStrip(geo()->feb(x.feb()));
		    TH1* hdts1=rh()->GetTH1(sr.str()+s.str());
		    if (hdts1==NULL)
		      {
			hdts1=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
		      }
		    hdts1->Fill(t1-tbcid-dtmean);
		    if (geo()->feb(x.feb()).polarity==-1)
		      {
			double tt=t1;
			t1=t0;
			t0=tt;
		      }
		    
		    lmana::TdcStrip ts(geo()->feb(x.feb()).chamber,x.feb(),x.detectorStrip(geo()->feb(x.feb())),t0,t1,geo()->feb(x.feb()).timePedestal[x.detectorStrip( geo()->feb(x.feb()))-70]);

		    //if (chamber==2)
		    // {
		    //  printf("New strip %d %d %d %f %f %f \n",geo()->feb(x.feb()).chamber,x.feb(),x.detectorStrip(x.feb()),t0,t1,geo()->feb(x.feb()).timePedestal[x.detectorStrip( geo()->feb(x.feb()))-70]);
		    //}
		    _strips.push_back(ts);

		    if (nstrip==1)
		      {
			std::stringstream s;
			s<<"Timing/OneStrip/hdtpos"<<(int) x.detectorStrip(geo()->feb(x.feb()));
			TH1* hdts=rh()->GetTH1(sr.str()+s.str());
			if (hdts==NULL)
			  {
			    hdts=rh()->BookTH1(sr.str()+s.str(),300,-25.,25.);
			  }
			hdts->Fill(t1-t0);
		      }
		    if (nstrip<=3)
		      {
			std::stringstream s;
			s.str("");
			s.clear();
			s<<"Timing/OneStrip/Side0/hdtr_"<<(int) x.detectorStrip(geo()->feb(x.feb()));
			TH1* hdts0=rh()->GetTH1(sr.str()+s.str());
			if (hdts0==NULL)
			  {
			    hdts0=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
			  }
			hdts0->Fill(t0-tbcid-dtmean);
			s.str("");
			s.clear();
			s<<"Timing/OneStrip/Side1/hdtr_"<<(int) x.detectorStrip(geo()->feb(x.feb()));
			TH1* hdts1=rh()->GetTH1(sr.str()+s.str());
			if (hdts1==NULL)
			  {
			    hdts1=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
			  }
			hdts1->Fill(t1-tbcid-dtmean);
		      }
		    break;
		  }
	      }
	  }
      //std::cout<<"NSTRIP "<<nstrip<<std::endl;
      //getchar();

      if (nstrip>=1) {heff->Fill(4.1);  hnstrip->Fill(nstrip*1.);}

      uint32_t nst[8];
      memset(nst,0,8*4);

      if (_strips.size()>0)
	{
	  DEBUG_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
	  DEBUG_PRINTF(" ======================================> Strips \n");

	  std::vector<lmana::TdcCluster> vclus;
	  vclus.clear();
	  float step=2.;
	  if (chamber==1) step=3.;
	  for (auto it=_strips.begin();it!=_strips.end();it++)
	    {
	      if (it->chamber()!=chamber) continue;
	      //	      if (it->ypos()<-10 || it->ypos()>-0.2) continue;
	      bool found=false;
	      for (auto ic=vclus.begin();ic!=vclus.end();ic++)
		{
		  if (ic->isAdjacent((*it),step))
		    {
		      ic->addStrip((*it));
		      found=true;
		      break;
		    }
		}
	      if(!found)
		{
		  lmana::TdcCluster c;
		  c.addStrip((*it));
		  vclus.push_back(c);
		}
	    }
	  std::stringstream src;
	  src<<"/run"<<_run<<"/Chamber"<<chamber<<"/Cluster/";
		  
	  TH2* hposc=rh()->GetTH2(src.str()+"XY");
	  TH2* hposc1=rh()->GetTH2(src.str()+"XY1");
	  TH2* hposcm=rh()->GetTH2(src.str()+"XYMore");
	  TH2* hposcma=rh()->GetTH2(src.str()+"XYMax");
	  TH1* hncl=rh()->GetTH1(src.str()+"Clusters");
	  TH1* hmulc=rh()->GetTH1(src.str()+"ClusterSize");

	  if (hposc==NULL)
	    {
	      
	      hposc=rh()->BookTH2(src.str()+"XY",128,0.,128.,256,-15.,15.);
	      hposc1=rh()->BookTH2(src.str()+"XY1",128,0.,128.,256,-15.,15.);
	      hposcm=rh()->BookTH2(src.str()+"XYMore",128,0.,128.,256,-15.,15.);
	      hposcma=rh()->BookTH2(src.str()+"XYMax",128,0.,128.,256,-15.,15.);
	      hncl=rh()->BookTH1(src.str()+"Clusters",16,0.,16.);
	      hmulc=rh()->BookTH1(src.str()+"ClusterSize",16,0.,16.);

	    }
	  printf(" %d  %d strips , Number of clusters %d \n",chamber,_strips.size(),vclus.size());
	  hncl->Fill(vclus.size()*1.);
	  uint32_t maxs=0;
	  for (auto x:vclus)
	    {
	      if (x.size()>8) continue;
	      hposc->Fill(x.X(),x.Y());
	      if (vclus.size()==1)
		hposc1->Fill(x.X(),x.Y());
	      else
		hposcm->Fill(x.X(),x.Y());
	      hmulc->Fill(x.size()*1.);
	      //printf("%d %f %f \n",x.size(),x.X(),x.Y());
	      if (x.size()>maxs) maxs=x.size();
	    }
	  for (auto x:vclus)
	    {
	      if (x.size()>8) continue;
	      if (x.size()==maxs)
		hposcma->Fill(x.X(),x.Y());
	    }

	  
	  //getchar();
	    
	  for (auto it=_strips.begin();it!=_strips.end();it++)
	    {
	      lmana::TdcStrip& x=(*it);
	      if (chamber==2)
		DEBUG_PRINTF("\t STRIP %d %d %f %f pos %f %f \n",x.dif(),x.strip(),x.t0(),x.t1(),x.xpos(),x.ypos());
	      nst[x.dif()/2]++;
	      std::stringstream sr;
	      uint32_t ich=x.chamber();
	      if (ich!=chamber) continue;
	      sr<<"/run"<<_run<<"/Chamber"<<chamber<<"/";
	  
	      TH2* hpos=rh()->GetTH2(sr.str()+"XY");
	      if (hpos==NULL)
		{
	      
		  hpos=rh()->BookTH2(sr.str()+"XY",128,0.,128.,128,-10.,10.);
	      
	      
		}
	      hpos->Fill(x.xpos(),x.ypos());
	    }
	  // getchar();

	  std::stringstream sr;

	  sr<<"/run"<<_run<<"/Chamber"<<chamber<<"/ChamberAll/";
		  
	  TH2* hpos=rh()->GetTH2(sr.str()+"XY");
	  if (hpos==NULL)
	    {
		      
	      hpos=rh()->BookTH2(sr.str()+"XY",128,0.,128.,256,-15.,15.);


	    }
	  if (_strips.size()==1)
	    hpos->Fill(_strips[0].xpos(),_strips[0].ypos());
	  if (_strips.size()==2)
	    hpos->Fill((_strips[0].xpos()+_strips[1].xpos())/2.,(_strips[0].ypos()+_strips[1].ypos())/2.);
	  if (_strips.size()==3)
	    hpos->Fill(_strips[1].xpos(),_strips[1].ypos());
	  if (_strips.size()==4)
	    hpos->Fill((_strips[2].xpos()+_strips[1].xpos())/2.,(_strips[2].ypos()+_strips[1].ypos())/2.);

	  if (_strips.size()>=5 && _strips.size()<=8)
	    {
	      double x=0,y=0;
	      for (int i=2;i<_strips.size()-2;i++)
		{
		  x+=_strips[i].xpos();
		  y+=_strips[i].ypos();
		}
	      x/=(_strips.size()-4);
	      y/=(_strips.size()-4);
	      hpos->Fill(x,y);
	    }
		

	  for (int i=71;i<71+48;i++)
	    {
	      float ti=-2000.,tj=-2000.;
	      for (auto it=_strips.begin();it!=_strips.end();it++)
		{
		  lmana::TdcStrip& x=(*it);
		  //if (x.dif()!=8) continue;
		  if (x.strip()==i) {ti=x.ypos();}
		  if (x.strip()==i+1) {tj=x.ypos();}
		}

	      if (ti>-100 && tj>-100 && _strips.size()<12)
		{
		  std::stringstream sd;
		      
		  sd<<"/run"<<_run<<"/ChamberDif/dif"<<i<<"_"<<i+1;

		  TH1* hpos=rh()->GetTH1(sd.str());
		  if (hpos==NULL)
		    {
		  
		      hpos=rh()->BookTH1(sd.str(),128,-10.,10.);


		    }
		  hpos->Fill(tj-ti);

		}
	    }

	      
	}
    }
#ifdef AFAIRE
	 
#endif
}

void lmana::TdcAnalyzer::fullAnalysis(std::vector<lydaq::TdcChannel>& vChannel)
{

  double fe1_2tr[128];
  memset(fe1_2tr,0,128*sizeof(double));
  fe1_2tr[71]=0.538;
  fe1_2tr[72]=2.684;
  fe1_2tr[73]=4.380;
  fe1_2tr[74]=6.970;
  fe1_2tr[75]=4.854;
  fe1_2tr[76]=3.512;
  fe1_2tr[77]=3.841;
  fe1_2tr[78]=2.302;
  fe1_2tr[79]=0.000;
  fe1_2tr[80]=3.909;
  fe1_2tr[81]=0.658;
  fe1_2tr[82]=0.000;
  fe1_2tr[83]=-0.868;
  fe1_2tr[84]=2.402;
  fe1_2tr[85]=4.748;
  fe1_2tr[86]=5.251;
  fe1_2tr[87]=3.620;
  fe1_2tr[88]=0.697;
  fe1_2tr[89]=1.292;
  fe1_2tr[90]=0.465;
  fe1_2tr[91]=3.092;
  fe1_2tr[92]=2.013;
  fe1_2tr[93]=-1.252;
  fe1_2tr[94]=0.000;
  fe1_2tr[95]=-3.725;
  fe1_2tr[96]=-0.851;
  fe1_2tr[97]=0.862;
  fe1_2tr[98]=2.577;
  fe1_2tr[99]=0.960;
  fe1_2tr[100]=-1.688;
  fe1_2tr[101]=-0.950;
  fe1_2tr[102]=-1.264;
  fe1_2tr[103]=0.236;
  fe1_2tr[104]=-0.445;
  fe1_2tr[105]=-3.270;
  fe1_2tr[106]=0.000;
  fe1_2tr[107]=-8.650;
  fe1_2tr[108]=-3.841;
  fe1_2tr[109]=-1.670;
  fe1_2tr[110]=-0.204;
  fe1_2tr[111]=-1.868;
  fe1_2tr[112]=-3.775;
  fe1_2tr[113]=-3.519;
  fe1_2tr[114]=-3.850;
  fe1_2tr[115]=-1.360;
  fe1_2tr[116]=-2.633;
  fe1_2tr[117]=-5.976;
  fe1_2tr[118]=0.000;

  //740473
  fe1_2tr[71]=-0.729;
  fe1_2tr[72]=5.578;
  fe1_2tr[73]=5.346;
  fe1_2tr[74]=6.979;
  fe1_2tr[75]=5.484;
  fe1_2tr[76]=2.787;
  fe1_2tr[77]=3.493;
  fe1_2tr[78]=5.086;
  fe1_2tr[79]=0.000;
  fe1_2tr[80]=4.209;
  fe1_2tr[81]=0.596;
  fe1_2tr[82]=0.000;
  fe1_2tr[83]=6.131;
  fe1_2tr[84]=4.099;
  fe1_2tr[85]=3.636;
  fe1_2tr[86]=5.273;
  fe1_2tr[87]=1.221;
  fe1_2tr[88]=0.899;
  fe1_2tr[89]=3.930;
  fe1_2tr[90]=0.601;
  fe1_2tr[91]=3.246;
  fe1_2tr[92]=4.602;
  fe1_2tr[93]=-1.179;
  fe1_2tr[94]=0.000;
  fe1_2tr[95]=-2.724;
  fe1_2tr[96]=-0.747;
  fe1_2tr[97]=1.005;
  fe1_2tr[98]=2.766;
  fe1_2tr[99]=1.104;
  fe1_2tr[100]=-1.629;
  fe1_2tr[101]=-1.153;
  fe1_2tr[102]=-1.796;
  fe1_2tr[103]=0.899;
  fe1_2tr[104]=-0.291;
  fe1_2tr[105]=-3.507;
  fe1_2tr[106]=0.000;
  fe1_2tr[107]=-7.014;
  fe1_2tr[108]=-3.733;
  fe1_2tr[109]=-4.378;
  fe1_2tr[110]=-2.629;
  fe1_2tr[111]=-1.741;
  fe1_2tr[112]=-4.752;
  fe1_2tr[113]=-3.914;
  fe1_2tr[114]=-4.568;
  fe1_2tr[115]=-1.864;
  fe1_2tr[116]=-3.137;
  fe1_2tr[117]=-6.686;
  fe1_2tr[118]=0.000;

  // 740478
  fe1_2tr[71]=0.000;
  fe1_2tr[72]=0.000;
  fe1_2tr[73]=0.000;
  fe1_2tr[74]=0.000;
  fe1_2tr[75]=0.000;
  fe1_2tr[76]=0.000;
  fe1_2tr[77]=0.000;
  fe1_2tr[78]=0.000;
  fe1_2tr[79]=0.000;
  fe1_2tr[80]=8.305;
  fe1_2tr[81]=0.000;
  fe1_2tr[82]=0.000;
  fe1_2tr[83]=2.105;
  fe1_2tr[84]=0.000;
  fe1_2tr[85]=0.000;
  fe1_2tr[86]=8.945;
  fe1_2tr[87]=3.062;
  fe1_2tr[88]=5.103;
  fe1_2tr[89]=8.276;
  fe1_2tr[90]=3.793;
  fe1_2tr[91]=6.835;
  fe1_2tr[92]=8.472;
  fe1_2tr[93]=2.664;
  fe1_2tr[94]=0.000;
  fe1_2tr[95]=0.865;
  fe1_2tr[96]=3.147;
  fe1_2tr[97]=4.950;
  fe1_2tr[98]=7.080;
  fe1_2tr[99]=5.983;
  fe1_2tr[100]=2.446;
  fe1_2tr[101]=2.696;
  fe1_2tr[102]=2.037;
  fe1_2tr[103]=4.759;
  fe1_2tr[104]=3.647;
  fe1_2tr[105]=0.141;
  fe1_2tr[106]=0.000;
  fe1_2tr[107]=0.000;
  fe1_2tr[108]=-0.161;
  fe1_2tr[109]=-1.075;
  fe1_2tr[110]=0.997;
  fe1_2tr[111]=1.908;
  fe1_2tr[112]=-1.301;
  fe1_2tr[113]=-0.374;
  fe1_2tr[114]=-1.463;
  fe1_2tr[115]=0.753;
  fe1_2tr[116]=-0.498;
  fe1_2tr[117]=-3.602;
  fe1_2tr[118]=0.000;

  // 740507
  fe1_2tr[71]=0.000;
  fe1_2tr[72]=0.000;
  fe1_2tr[73]=0.000;
  fe1_2tr[74]=0.000;
  fe1_2tr[75]=0.000;
  fe1_2tr[76]=0.000;
  fe1_2tr[77]=0.000;
  fe1_2tr[78]=0.000;
  fe1_2tr[79]=0.000;
  fe1_2tr[80]=5.018;
  fe1_2tr[81]=2.981;
  fe1_2tr[82]=0.000;
  fe1_2tr[83]=6.550;
  fe1_2tr[84]=2.139;
  fe1_2tr[85]=1.201;
  fe1_2tr[86]=2.231;
  fe1_2tr[87]=1.144;
  fe1_2tr[88]=0.898;
  fe1_2tr[89]=1.089;
  fe1_2tr[90]=0.557;
  fe1_2tr[91]=0.463;
  fe1_2tr[92]=2.280;
  fe1_2tr[93]=1.013;
  fe1_2tr[94]=0.000;
  fe1_2tr[95]=-1.166;
  fe1_2tr[96]=-0.423;
  fe1_2tr[97]=-1.491;
  fe1_2tr[98]=-0.280;
  fe1_2tr[99]=-1.530;
  fe1_2tr[100]=-1.665;
  fe1_2tr[101]=-1.579;
  fe1_2tr[102]=-1.863;
  fe1_2tr[103]=-2.033;
  fe1_2tr[104]=-0.032;
  fe1_2tr[105]=-1.390;
  fe1_2tr[106]=0.000;
  fe1_2tr[107]=0.000;
  fe1_2tr[108]=-3.692;
  fe1_2tr[109]=-4.508;
  fe1_2tr[110]=-3.352;
  fe1_2tr[111]=-4.442;
  fe1_2tr[112]=-4.755;
  fe1_2tr[113]=-4.512;
  fe1_2tr[114]=-4.733;
  fe1_2tr[115]=-4.545;
  fe1_2tr[116]=-2.749;
  fe1_2tr[117]=-3.738;
  fe1_2tr[118]=0.000;
  std::stringstream sr;
  
  sr<<"/run"<<_run<<"/";
  TH2* hdtr=rh()->GetTH2(sr.str()+"DeltaTrigger");
  TH2* hdtrt=rh()->GetTH2(sr.str()+"DeltaTriggerSel");

  TH2* hdtrt0=rh()->GetTH2(sr.str()+"DeltaTriggerSel0");
  TH2* hdtrt1=rh()->GetTH2(sr.str()+"DeltaTriggerSel1");
  TH1* hnstrip=rh()->GetTH1(sr.str()+"NSTRIP");
  TH1* heff=rh()->GetTH1(sr.str()+"Efficiency");
  TH1* hbp2=rh()->GetTH1(sr.str()+"BeamProfile");
  TH2* hpos=rh()->GetTH2(sr.str()+"DeltaTvsStrip");
  if (hdtr==NULL)
    {
      hdtr=rh()->BookTH2(sr.str()+"DeltaTrigger",4000,-1000.,0.,48,71.,119.);
      hdtrt=rh()->BookTH2(sr.str()+"DeltaTriggerSel",4000,-600.,-560,48,71.,119.);
      hdtrt0=rh()->BookTH2(sr.str()+"DeltaTriggerSel0",4000,-600.,-560,48,71.,119.);
      hdtrt1=rh()->BookTH2(sr.str()+"DeltaTriggerSel1",4000,-600.,-560,48,71.,119.);
      hpos=rh()->BookTH2(sr.str()+"DeltaTvsStrip",3000,-30.,30.,48,71.,119.);
      hnstrip=rh()->BookTH1(sr.str()+"NSTRIP",24,-0.1,23.9);
      heff=rh()->BookTH1(sr.str()+"Efficiency",10,-0.1,9.9);
      hbp2=rh()->BookTH1(sr.str()+"BeamProfile",128,0.1,128.1);

    }



  
  std::bitset<16> btrg(0);
  for (auto x:vChannel)
    {
      if (x.channel()==24) {
	//printf("Trigger found %d %d %f\n",x.feb(),x.channel(),x.tdcTime());
	btrg.set(x.feb(),1);
      }
    }
  _triggerFound=(btrg.count()==8);
  if (btrg.count()==8) _ntrigger++;
  heff->Fill(1.1);
  if (!_triggerFound) return;
  heff->Fill(2.1);

  std::bitset<128> side[2];
  side[0].reset();
  side[1].reset();
  double dtmin=-640.,dtmax=-590.,dtmean=-640.;
  double febbcid[128];
  memset(febbcid,0,128*sizeof(double));
  for (int idif=0;idif<24;idif++)
    {
      if (FEB2STRIP[idif]==255) continue;
      double tbcid=0;

      for (auto x:vChannel)
	{
	  if (x.feb()!=idif) continue;
	  if (x.channel()!=24) continue;
	  tbcid=x.tdcTime();
	  //printf("TRIGGER FOUND %f \n",tbcid);
	  //getchar();
	  break;
	}
      febbcid[idif]=tbcid;
      for (auto x:vChannel)
	{
	  if (x.feb()!=idif) continue;
	  //printf("strip %f %f %d \n",x.tdcTime(),tbcid,x.strip());
	  hdtr->Fill(x.tdcTime()-tbcid,x.detectorStrip(x.feb()));
	}

      // Now fill those in good range
      for (auto x:vChannel)
	{
	  if (x.feb()!=idif) continue;
	  if ((x.tdcTime()-tbcid)<dtmin) continue;
	  if ((x.tdcTime()-tbcid)>dtmax) continue;
	  
	  //printf("strip %f %f %d \n",x.tdcTime(),tbcid,x.strip());
	  hdtrt->Fill(x.tdcTime()-tbcid,x.detectorStrip(x.feb()));
	  side[x.side()].set(x.detectorStrip(x.feb()),1);
	  if (x.side())
	    {
	      hdtrt1->Fill(x.tdcTime()-tbcid,x.detectorStrip(x.feb()));

	    }
	  else
	    hdtrt0->Fill(x.tdcTime()-tbcid,x.detectorStrip(x.feb()));
	}
  
    }
  //std::cout<<side[0]<<std::endl;
  // std::cout<<side[1]<<std::endl;
  if (side[0].count()!=0 || side[1].count()!=0) heff->Fill(3.1);
  uint32_t nstrip=0;
  for (int i=0;i<128;i++)
    if (side[0][i]&&side[1][i])
      {
	nstrip++;
	hbp2->Fill(i*1.);
      }
  for (int i=0;i<128;i++)
    if (side[0][i]&&side[1][i])
      {
	double t0=-1,t1=-1;
	for (auto x:vChannel)
	  {
	    double tbcid=febbcid[x.feb()];
	    if (tbcid==0) continue;
	    if ((x.tdcTime()-tbcid)<dtmin) continue;
	    if ((x.tdcTime()-tbcid)>dtmax) continue;
	    if (x.detectorStrip(x.feb())!=i) continue;
	    if (x.side()==0 && t0==-1) 	    t0=x.tdcTime();
	    if (x.side()==1 && t1==-1) 	    t1=x.tdcTime();
	      
	
	    if(t0>0 && t1>0)
	      {
		//printf("%f %f %f\n",t0,t1,t1-t0);
		hpos->Fill(t1-t0,x.detectorStrip(x.feb()));
		std::stringstream s;
		s<<"Timing/All/hdtpos"<<(int) x.detectorStrip(x.feb());
		TH1* hdts=rh()->GetTH1(sr.str()+s.str());
		if (hdts==NULL)
		  {
		    hdts=rh()->BookTH1(sr.str()+s.str(),300,-25.,25.);
		  }
		hdts->Fill(t1-t0);
		s.str("");
		s.clear();
		s<<"Timing/Side0/hdtr_"<<(int) x.detectorStrip(x.feb());
		TH1* hdts0=rh()->GetTH1(sr.str()+s.str());
		if (hdts0==NULL)
		  {
		    hdts0=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
		  }
		hdts0->Fill(t0-tbcid-dtmean);
		s.str("");
		s.clear();
		s<<"Timing/Side1/hdtr_"<<(int) x.detectorStrip(x.feb());
		TH1* hdts1=rh()->GetTH1(sr.str()+s.str());
		if (hdts1==NULL)
		  {
		    hdts1=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
		  }
		hdts1->Fill(t1-tbcid-dtmean);

		lmana::TdcStrip ts(x.feb(),x.detectorStrip(x.feb()),t0,t1,fe1_2tr[x.detectorStrip(x.feb())]);
		_strips.push_back(ts);

		if (nstrip==1)
		  {
		    std::stringstream s;
		    s<<"Timing/OneStrip/hdtpos"<<(int) x.detectorStrip(x.feb());
		    TH1* hdts=rh()->GetTH1(sr.str()+s.str());
		    if (hdts==NULL)
		      {
			hdts=rh()->BookTH1(sr.str()+s.str(),300,-25.,25.);
		      }
		    hdts->Fill(t1-t0);
		  }
		if (nstrip<=3)
		  {
		    std::stringstream s;
		    s.str("");
		    s.clear();
		    s<<"Timing/OneStrip/Side0/hdtr_"<<(int) x.detectorStrip(x.feb());
		    TH1* hdts0=rh()->GetTH1(sr.str()+s.str());
		    if (hdts0==NULL)
		      {
			hdts0=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
		      }
		    hdts0->Fill(t0-tbcid-dtmean);
		    s.str("");
		    s.clear();
		    s<<"Timing/OneStrip/Side1/hdtr_"<<(int) x.detectorStrip(x.feb());
		    TH1* hdts1=rh()->GetTH1(sr.str()+s.str());
		    if (hdts1==NULL)
		      {
			hdts1=rh()->BookTH1(sr.str()+s.str(),100,dtmin-dtmean,dtmax-dtmean);
		      }
		    hdts1->Fill(t1-tbcid-dtmean);
		  }
		break;
	      }
	  }
      }
  //std::cout<<"NSTRIP "<<nstrip<<std::endl;
  //getchar();

  if (nstrip>=1) {heff->Fill(4.1);  hnstrip->Fill(nstrip*1.);}



  uint32_t nst[8];
  memset(nst,0,8*4);

  if (_strips.size()>0)
    {
      DEBUG_PRINTF("================> Event %d Number of DIF found %d \n",_event,theNumberOfDIF);
      DEBUG_PRINTF(" ======================================> Strips \n");
      
      
      for (auto it=_strips.begin();it!=_strips.end();it++)
	{
	  lmana::TdcStrip& x=(*it);
	  DEBUG_PRINTF("\t STRIP %d %d %f %f pos %f %f \n",x.dif(),x.strip(),x.t0(),x.t1(),x.xpos(),x.ypos());
	  nst[x.dif()/2]++;
	  std::stringstream sr;
	  uint32_t ich=x.dif();
	  sr<<"/run"<<_run<<"/Chamber"<<ich<<"/";
	  
	  TH2* hpos=rh()->GetTH2(sr.str()+"XY");
	  if (hpos==NULL)
	    {
	      
	      hpos=rh()->BookTH2(sr.str()+"XY",128,0.,128.,128,-10.,10.);
	      
	      
	    }
	  hpos->Fill(x.xpos(),x.ypos());
	}
      // getchar();

      std::stringstream sr;

      sr<<"/run"<<_run<<"/ChamberAll/";
		  
      TH2* hpos=rh()->GetTH2(sr.str()+"XY");
      if (hpos==NULL)
	{
		      
	  hpos=rh()->BookTH2(sr.str()+"XY",128,0.,128.,256,-15.,15.);


	}
      if (_strips.size()==1)
	hpos->Fill(_strips[0].xpos(),_strips[0].ypos());
      if (_strips.size()==2)
	hpos->Fill((_strips[0].xpos()+_strips[1].xpos())/2.,(_strips[0].ypos()+_strips[1].ypos())/2.);
      if (_strips.size()==3)
	hpos->Fill(_strips[1].xpos(),_strips[1].ypos());
      if (_strips.size()==4)
	hpos->Fill((_strips[2].xpos()+_strips[1].xpos())/2.,(_strips[2].ypos()+_strips[1].ypos())/2.);

      if (_strips.size()>=5 && _strips.size()<=12)
	{
	  double x=0,y=0;
	  for (int i=2;i<_strips.size()-2;i++)
	    {
	      x+=_strips[i].xpos();
	      y+=_strips[i].ypos();
	    }
	  x/=(_strips.size()-4);
	  y/=(_strips.size()-4);
	  hpos->Fill(x,y);
	}
		

      for (int i=71;i<71+48;i++)
	{
	  float ti=-2000.,tj=-2000.;
	  for (auto it=_strips.begin();it!=_strips.end();it++)
	    {
	      lmana::TdcStrip& x=(*it);
	      //if (x.dif()!=8) continue;
	      if (x.strip()==i) {ti=x.ypos();}
	      if (x.strip()==i+1) {tj=x.ypos();}
	    }

	  if (ti>-100 && tj>-100 && _strips.size()<12)
	    {
	      std::stringstream sd;
		      
	      sd<<"/run"<<_run<<"/ChamberDif/dif"<<i<<"_"<<i+1;

	      TH1* hpos=rh()->GetTH1(sd.str());
	      if (hpos==NULL)
		{
		  
		  hpos=rh()->BookTH1(sd.str(),128,-10.,10.);


		}
	      hpos->Fill(tj-ti);

	    }
	}

	      
    }
	 

}
void lmana::TdcAnalyzer::pedestalAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel)
{
  _pedestalProcessed=true;

  std::cout<<"Mezzanine "<<mezId<<"Event "<<_event<<" GTC"<<_gtc<<" hits"<<vChannel.size()<<" Trigger channel "<<_triggerChannel<<std::endl;

  // Analyze
  std::stringstream sr;
  sr<<"/run"<<_run<<"/TDC"<<mezId<<"/";

  uint32_t dac =_dacSet;
  for (int ich=0;ich<_triggerChannel+1;ich++)
    {
 
      std::stringstream src;
      src<<sr.str()<<"dac"<<ich;
      TH1* hdac=rh()->GetTH1(src.str());
      if (hdac==NULL)
	{
	 
	  hdac=rh()->BookTH1(src.str(),64,0.,64.);
	}
      bool found=false;
      double lastf=0;
      for (std::vector<lydaq::TdcChannel>::iterator x=vChannel.begin();x!=vChannel.end();x++)
	{
	  if (x->channel()==ich) {
	    //printf("%d %d %f \n",x.channel(),x.bcid(),x.tdcTime());
	    double dt=x->tdcTime()-lastf;
	    lastf=x->tdcTime();
	    if (dt>25 || dt<0)
	      hdac->Fill(dac*1.);
	  }
	}
    }

}
void lmana::TdcAnalyzer::scurveAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel)
{

  //if (gtc[mezId-1]
  uint32_t maxChannels=64;
  std::cout<<"Mezzanine "<<mezId<<"Event "<<_event<<" GTC"<<_gtc<<" hits"<<vChannel.size()<<" Vth set "<<_vthSet<<" Trigger channel "<<maxChannels<<std::endl;

  // Analyze
  std::stringstream sr;
  sr<<"/run"<<_run<<"/TDC"<<mezId<<"/";
  
  uint32_t vth =_vthSet;
  for (int ich=0;ich<maxChannels+1;ich++)
    {
 
      bool found=false;
      double lastf=0;
      std::stringstream src;
      src<<sr.str()<<"vth"<<ich;

      for (std::vector<lydaq::TdcChannel>::iterator x=vChannel.begin();x!=vChannel.end();x++)
	{
	  if (x->feb()!=mezId) continue;
	  if (x->channel()==ich) {

	    TH1* hvth=rh()->GetTH1(src.str());
	    if (hvth==NULL)
	      {
		
		hvth=rh()->BookTH1(src.str(),900,0.,900.);
		printf("Booking %s \n",src.str().c_str());
	      }

	    printf("Found %d  %d \n",x->channel(),vth);
	    double dt=x->tdcTime()-lastf;
	    lastf=x->tdcTime();
	    hvth->Fill(vth*1.);
	    break;
	  }
	}
    }
  std::stringstream srt;
  srt<<sr.str()<<"maxtime";
  TH1* hti=rh()->GetTH1(srt.str());
  if (hti==NULL)
    {
	 
      hti=rh()->BookTH1(srt.str(),90000,0.,90000.);
      printf("Booking %s \n",srt.str().c_str());
    }
  double maxt=0;
  for (int ich=0;ich<maxChannels+1;ich++)
    {
 
      std::stringstream src;
      src<<sr.str()<<"vthc"<<ich;
      TH1* hvthc=rh()->GetTH1(src.str());
      std::stringstream srd;
      srd<<sr.str()<<"vthd"<<ich;
      TH1* hvthd=rh()->GetTH1(srd.str());

      bool found=false;
      double lastf=0;
      for (std::vector<lydaq::TdcChannel>::iterator x=vChannel.begin();x!=vChannel.end();x++)
	{
	  // x->dump();
	  if (x->feb()!=mezId) continue;
	  if (x->tdcTime()>maxt) maxt=x->tdcTime();
	  if (x->channel()==ich) {



	    //printf("%d \n",x->channel());
	    double dt=x->tdcTime()-lastf;
	    lastf=x->tdcTime();
	    if (x->falling())
	      {
		if (hvthd==NULL)
		  {
		    
		    hvthd=rh()->BookTH1(srd.str(),900,0.,900.);

		  }

	      hvthd->Fill(vth*1.);
	      }
	    else
	      {
		if (hvthc==NULL)
		  {
		
		    hvthc=rh()->BookTH1(src.str(),900,0.,900.);
		
		    printf("Booking %s \n",src.str().c_str());
		  }

	      hvthc->Fill(vth*1.);
	      }
	  }
	}
    }
  hti->Fill(maxt);
  printf("MAXTIME %f %f \n",maxt,maxt*2.5E-9);
  //  if (maxt>0)
  //  getchar();
}
void lmana::TdcAnalyzer::normalAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel)
{
  this->LmAnalysis(mezId,vChannel);
}

void lmana::TdcAnalyzer::end()
{

  if (_pedestalProcessed)
    {
      for (int mez=1;mez<=255;mez++)
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
	      TH1* hdac=rh()->GetTH1(src.str());
	      if (hdac==NULL) continue;
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
	       
		  if (ped==0)
		    {printf("\t \t ======================>DEAD %d \n",ipr);
		      ped=31;
		    }
		  printf("\t %d %d \n",ipr,ped);
		}
	    }


	}
    }




  // std::stringstream sr;
  // sr<<"/tmp/toto"<<_run<<".root";
  
  // rh()->writeHistograms(sr.str());


}

void lmana::TdcAnalyzer::LmAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel)
{
  //if (vChannel.size()==254) return;
  //printf("%d %d %d \n",_event,mezId,vChannel.size());
  double fe2_shift[128];
  double fe1_shift[128];
  fe1_shift[71]=-1.36;
  fe1_shift[72]=-0.138;
  fe1_shift[73]=-0.496;
  fe1_shift[74]=-0.415;
  fe1_shift[75]=-0.076;
  fe1_shift[76]=-0.297;
  fe1_shift[77]=-0.032;
  fe1_shift[78]=-0.133;
  fe1_shift[79]=-0.121;
  fe1_shift[80]=-0.073;
  fe1_shift[81]=0.410;
 
  double fe1_2tr[32];
  memset(fe1_2tr,0,32*sizeof(double));
  memset(fe1_shift,0,128*sizeof(double));

  double fe8_739331[12]={5.000000000000002, 6.805555555555557, 9.333333333333334, 10.634057971014494, 9.129901960784316, 6.8024344569288395, 7.720125786163522, 6.672979797979799, 9.41666666666667, 7.990740740740743, 4.895833333333334, 0};

  for (int i=0;i<12;i++){fe1_shift[71+i]=fe8_739331[i];}


  // From TDC8 (FE#5) run 739303
  fe1_2tr[0]=24.279;
  fe1_2tr[1]=13.385;
  fe1_2tr[2]=23.355;
  fe1_2tr[3]=12.594;
  fe1_2tr[4]=23.170;
  fe1_2tr[5]=16.227;
  fe1_2tr[6]=22.260;
  fe1_2tr[7]=13.219;
  fe1_2tr[8]=23.474;
  fe1_2tr[9]=13.109;
  fe1_2tr[10]=20.214;
  fe1_2tr[11]=12.218;
  fe1_2tr[12]=23.346;
  fe1_2tr[13]=12.962;
  fe1_2tr[14]=14.287;
  fe1_2tr[15]=12.735;
  fe1_2tr[16]=19.923;
  fe1_2tr[17]=14.467;
  fe1_2tr[18]=22.203;
  fe1_2tr[19]=15.252;
  fe1_2tr[20]=22.474;
  fe1_2tr[21]=14.479;
  fe1_2tr[22]=19.608;
  fe1_2tr[23]=11.500;

  fe1_shift[71+0]=-0.000;
  fe1_shift[71+1]=1.750;
  fe1_shift[71+2]=1.083;
  fe1_shift[71+3]=1.392;
  fe1_shift[71+4]=1.188;
  fe1_shift[71+5]=0.797;
  fe1_shift[71+6]=0.515;
  fe1_shift[71+7]=0.656;
  fe1_shift[71+8]=0.611;
  fe1_shift[71+9]=1.031;
  fe1_shift[71+10]=0.431;
  fe1_shift[71+11]=0.000;

  for (int i=0;i<12;i++){fe1_shift[71+i]+=fe8_739331[i];}

  /* Run jusqu'au 739543
     fe1_shift[75]+=0.2148;
     fe1_shift[76]+=0.2148-0.1172;
     fe1_shift[77]+=0.2148-0.1172-0.3095;
     fe1_shift[78]+=0.2148-0.1172-0.3095+0.08156;
     fe1_shift[79]+=0.2148-0.1172-0.3095+0.08156+0.1106;
     fe1_shift[80]+=0.2148-0.1172-0.3095+0.08156+0.1106-0.168;
     fe1_shift[81]+=0.2148-0.1172-0.3095+0.08156+0.1106-0.168+0.4425;
  */
  // fe1_shift[72]+=-2.5;
  // fe1_shift[73]+=-2.5-2.5;
  // fe1_shift[74]+=-2.5-2.5+0.1225;
  // fe1_shift[75]+=-2.5-2.5+0.1225+0.066;
  // fe1_shift[76]+=-2.5-2.5+0.1225+0.066+2.525;
  // fe1_shift[77]+=-2.5-2.5+0.1225+0.066+2.525-0.4326;
  // fe1_shift[78]+=-2.5-2.5+0.1225+0.066+2.525-0.4326-0.0352;
  // fe1_shift[79]+=-2.5-2.5+0.1225+0.066+2.525-0.4326-0.0352-2.579;
  // fe1_shift[80]+=-2.5-2.5+0.1225+0.066+2.525-0.4326-0.0352-2.579+2.587;
  // fe1_shift[81]+=-2.5-2.5+0.1225+0.066+2.525-0.4326-0.0352-2.579+2.93;
  memset(fe1_shift,0,128*sizeof(double));

  /* 739560 center to -1.4 
     fe1_shift[72]=5.785;
     fe1_shift[73]=5.785-0.681;
     fe1_shift[74]=5.785-0.681+0.463;
     fe1_shift[75]=5.785-0.681+0.463-4.13;
     fe1_shift[76]=5.785-0.681+0.463-4.13-0.27;
     fe1_shift[77]=5.785-0.681+0.463-4.13-0.27+2.647;
     fe1_shift[78]=5.785-0.681+0.463-4.13-0.27+2.467-3.237;
     fe1_shift[79]=5.785-0.681+0.463-4.13-0.27+2.467-3.237+2.73;
     fe1_shift[80]=5.785-0.681+0.463-4.13-0.27+2.467-3.237+2.73-1.035;
     fe1_shift[81]=5.785-0.681+0.463-4.13-0.27+2.467-3.237+2.73-1.035-3.352;
  */
  double fe1_diff[128];
  fe1_diff[72]=5.474;
  fe1_diff[73]=-0.5078;
  fe1_diff[74]=1.567;
  fe1_diff[75]=-4.277;
  fe1_diff[76]=-0.2433;
  fe1_diff[77]=2.71;
  fe1_diff[78]=-3305;
  fe1_diff[79]=2.710;
  fe1_diff[80]=-1.048;
  fe1_diff[81]=-3.354;

  for (int i=72;i<81;i++)
    {fe1_shift[i]=0; for (int j=72;j<=i;j++) fe1_shift[i]+=fe1_diff[j];}

  memset(fe1_shift,0,128*sizeof(double));

 
 
 
 

 
  memset(fe1_2tr,0,32*sizeof(double));
  if (vChannel.size()>4096) return; // Skip noise event
  double tmi=1E33,tma=-1E33;
  for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();it++)
    {
      if (it->tdcTime()<tmi) tmi=it->tdcTime();
      if (it->tdcTime()>tma) tma=it->tdcTime();
    }
  // Make a cut on maximal time
  //if (abs(tma-tmi)>0.5) return;

  // Analyze

  std::stringstream sr;
  std::stringstream difname;
  std::stringstream runname;
  difname<<mezId;
  runname<<_run;
  sr<<"/run"<<_run<<"/TDC"<<mezId<<"/LmAnalysis/";
  TH2* hpos=rh()->GetTH2(sr.str()+"Position");
  TH2* hpost=rh()->GetTH2(sr.str()+"Post");
  TH1* hdt=rh()->GetTH1(sr.str()+"DeltaT");
  TH1* hnti=rh()->GetTH1(sr.str()+"nvstime");
  TH1* hmti=rh()->GetTH1(sr.str()+"muvstime");
  TH1* hdtr0=rh()->GetTH1(sr.str()+"DeltaTr0");
  
  TH2* hdtr=rh()->GetTH2(sr.str()+"DeltaTr");
  TH2* hcor=rh()->GetTH2(sr.str()+"corr");
  TH2* hdtra=rh()->GetTH2(sr.str()+"DeltaTrAll");
  TH1* hns=rh()->GetTH1(sr.str()+"NChannel");
  TH1* hnst=rh()->GetTH1(sr.str()+"NChannelTrigger");
  TH1* hfin=rh()->GetTH1(sr.str()+"Fine");
  TH1* heff=rh()->GetTH1(sr.str()+"Efficiency");
  TH1* hstrip=rh()->GetTH1(sr.str()+"Strips");
  TH1* hstripo=rh()->GetTH1(sr.str()+"Stripo");
  TH1* hstript=rh()->GetTH1(sr.str()+"Stript");
  TH1* hstriptb=rh()->GetTH1(sr.str()+"Striptb");
  TH1* hnstrip=rh()->GetTH1(sr.str()+"NStrips");
  TH1* hstript2=rh()->GetTH1(sr.str()+"Stript2");
  TH1* hstript1=rh()->GetTH1(sr.str()+"Stript1");
  TH1* hnstrip2=rh()->GetTH1(sr.str()+"NStrips2");

  TH1* hxp=rh()->GetTH1(sr.str()+"XP");
  TH1* hti=rh()->GetTH1(sr.str()+"time");
  TH1* hra=rh()->GetTH1(sr.str()+"rate");
  if (hpos==NULL)
    {
      hpos=rh()->BookTH2(sr.str()+"Position",40,0.,50.,200,-100.,100.);
      hpost=rh()->BookTH2(sr.str()+"Post",128,0.,128.,1000,-25.,+25.);
      hdt=rh()->BookTH1(sr.str()+"DeltaT",1500,-30.,30.);

      hcor=rh()->BookTH2(sr.str()+"corr",32,0.1,32.1,32,0.1,32.1);
      hdtr=rh()->BookTH2(sr.str()+"DeltaTr",32,0,32.,4000,-400.,400.);
      hnti=rh()->BookProfile(sr.str()+"nvstime",50000,0,14400.,0.,20000.);
      hmti=rh()->BookProfile(sr.str()+"muvstime",50000,0,14400.,0.,20000.);
      
      hdtr0=rh()->BookTH1(sr.str()+"DeltaTr0",10000,-1000.,2500.);
      hdtra=rh()->BookTH2(sr.str()+"DeltaTrAll",32,0,32.,500,-500.,500.);
      hns=rh()->BookTH1(sr.str()+"NChannel",1024,0.,1024.);
      hnst=rh()->BookTH1(sr.str()+"NChannelTrigger",1024,0.,1024.);
      hfin=rh()->BookTH1(sr.str()+"Fine",257,0.,257.);
      heff=rh()->BookTH1(sr.str()+"Efficiency",32,0.,32.);
      hstrip=rh()->BookTH1(sr.str()+"Strips",32,0.,32.);
      hstript=rh()->BookTH1(sr.str()+"Stript",32,0.,32.);
      hstriptb=rh()->BookTH1(sr.str()+"Striptb",64,0.,64.);
      hstripo=rh()->BookTH1(sr.str()+"Stripo",64,0.,64.);
      hnstrip=rh()->BookTH1(sr.str()+"NStrips",32,0.,32.);
      hstript2=rh()->BookTH1(sr.str()+"Stript2",32,0.,32.);
      hstript1=rh()->BookTH1(sr.str()+"Stript1",32,0.,32.);
      hnstrip2=rh()->BookTH1(sr.str()+"NStrips2",32,0.,32.);
      hxp=rh()->BookTH1(sr.str()+"XP",400,0.,10.);
      hti=rh()->BookTH1(sr.str()+"time",90000,0.,0.5);
      hra=rh()->BookTH1(sr.str()+"rate",750,0.,200000.);

    }

  // Number of channel raw
  hns->Fill((vChannel.size()+1)*1.);

  _nevt++;
  heff->Fill(0.1); // First bin number of window
  
  // Profile histo
  hnti->Fill(acquisitionTime(),vChannel.size());
  uint32_t ndeclenchement=0;
  
  double ti=0,tmax=0;
  uint32_t lbcid=0,bcidshift=0,bcidmax=0;
  uint32_t tbcid=0;
  double ttime=0;
  _triggerChannel=0;
  if (_event%1000==0)
    printf("Event %d DIF %d GTC %d ABCID %lu Size %d %10.3f \n",_event,mezId,_gtc,_abcid,vChannel.size(),acquisitionTime());
  for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();it++)
    {
      hstrip->Fill(it->channel()*1.+0.5);
      if (it->channel()!=_triggerChannel) hstripo->Fill(32*it->side()+it->strip()-70);
      if (it->bcid()>bcidmax) bcidmax=it->bcid();

	  
      // if (it->bcid()<lbcid) {
      //   printf("lbcid %d bcid %d  Shift %d \n",lbcid,it->bcid(),b
      //   bcidshift+=65535;
      // }
      lbcid=it->bcid();
      //float t=((int) it->bcid()*2E-7)+(bcidshift*2.E-7)-ti;
      double t =it->tdcTime();
      if (t>tmax) tmax=t;
      // Find trigger bcid ( last one)
      if (it->channel()==_triggerChannel) {
	ndeclenchement++; tbcid=it->bcid();ttime=it->tdcTime();
	INFO_PRINTF("Event %d DIF %d GTC %d ABCID %lu BCID trigger %d # %d %f \n",_event,mezId,_gtc,_abcid,tbcid,ndeclenchement,ttime);
      }
      //printf("%d %d %d %d %d  %f \n",_gtc,it->channel(),it->coarse(),it->fine(),it->bcid(),it->tdcTime());
    }
  // //getchar();
  // fprintf(stderr,"BCID max %d Bcidshift %d Tmax %f \n",bcidmax,bcidshift,tmax);
  // for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();it++)
  //   {
  //     if (it->bcid()==bcidmax)
  // 	{
  // 	  fprintf(stderr,"c %d f %d t %f bc %d \n",it->coarse(),it->fine(),it->tdcTime(),it->bcid()); 
  // 	}
  //     if (it->tdcTime()==tmax)
  // 	{
  // 	  fprintf(stderr,"c %d f %d t %f bc %d \n",it->coarse(),it->fine(),it->tdcTime(),it->bcid()); 
  // 	}
  //   }

  // if (tmax==0 && vChannel.size()>0)
  //   {
  //     for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();it++)
  // 	std::cout<<(int) it->channel()<<" "<<it->coarse()*2.5E-9<<" "<<it->bcid()*2E-7<<endl;
  //     //getchar();
  //   }
  // Strore the maximum time of acquisition
  if (tmax==0) tmax=0.5E9;
  if (ndeclenchement==0) hti->Fill(tmax*1.E-9);
  // Calculate channel occupancy

  float ncx=vChannel.size();
  if (ncx==0) ncx=1;
  if (ndeclenchement==0) hra->Fill(ncx/tmax/1E-9);
  // Accept events with only one trigger
  if (ndeclenchement!=1) return;
  // Find the trigger
  hmti->Fill(acquisitionTime(),1.);
  heff->Fill(1.1);


  // Map of tdc channels per declenchement
  std::map<uint32_t,std::vector<lydaq::TdcChannel> > _trmap;
  
  _trmap.clear();
  // Fill bcid distance of hits wrt to trigger
  for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();it++)
    {
      it->setUsed(false);
      //printf(" channel %d bcid %d Time %f \n",it->channel(),it->bcid(),it->tdcTime());
      if (it->channel()!=_triggerChannel) {

	hdtra->Fill(it->channel(),(it->bcid()-tbcid)*1.);
	hdtr->Fill(it->channel(),(it->tdcTime()-ttime)*1.);
      }
    }
  //getchar();
  for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();++it)
    {
      // Find trigger channel
      if (it->channel()==_triggerChannel)
	{
	  it->setUsed(true);
	  std::vector<lydaq::TdcChannel> vc;
	  vc.push_back(*it);
	  
	  // Loop on hits and find nearby channels hits
	  
	  for (std::vector<lydaq::TdcChannel>::iterator x=vChannel.begin();x!=vChannel.end();++x)
	    {
	      if (x->used()) continue;
	      if (x->channel() == _triggerChannel) continue;
	      //if (x->bcid()>(it->bcid()-2) || x->bcid()<(it->bcid()-3)) continue;
	      if ((x->tdcTime()-it->tdcTime())<-220) continue;
	      if ((x->tdcTime()-it->tdcTime())>-160) continue;

	      
#ifdef TIMECORRCERN
	      if ((x->tdcTime()-it->tdcTime())<-220) continue;
	      if ((x->tdcTime()-it->tdcTime())>-160) continue;
#endif
	      //if ((x->tdcTime()-it->tdcTime())<-165) continue;
	      //if ((x->tdcTime()-it->tdcTime())>-135) continue;

	      //printf("\t TDC %d LEMO %d STRIP %d  SIDE %d   time %f \n",x->channel(),x->lemo(),x->strip(),x->side(),  (x->tdcTime()-it->tdcTime()));
	      vc.push_back((*x));
	      x->setUsed(true);
	    }
	  // Insert bcid, vector of hits in the trigger map
	  std::pair<uint32_t,std::vector<lydaq::TdcChannel> > p(it->bcid(),vc);
	  _trmap.insert(p);
	  //getchar();

	}
      
    }

  // Now loop on all channel of the event
  for (std::vector<lydaq::TdcChannel>::iterator it=vChannel.begin();it!=vChannel.end();it++)
    {
      hstrip->Fill(it->channel()*1.);
    }

  if (_trmap.size()>0) DEBUG_PRINTF("TDC %d  GTC %d   Number %d \n",mezId,_gtc,_trmap.size());
  bool found=false;
  bool bside=false;
  int nt=0;
  for ( std::map<uint32_t,std::vector<lydaq::TdcChannel> >::iterator t=_trmap.begin();t!=_trmap.end();t++)
    {
      //if (nt) break;
      nt++;
      DEBUG_PRINTF("Trigger %d => channels %d  \n",t->first,t->second.size());
      double trigtime=0,trigbcid=0;
      bool chused[32]={32*0};
      bool sused[96]={96*0};
      for (int i=0;i<32;i++) {chused[i]=false;}
      for (int i=0;i<96;i++){sused[i]=false;}
      std::sort(t->second.begin(),t->second.end());
      double tev0=0;
      std::vector<lydaq::TdcChannel>::iterator ittrig=t->second.end();
      std::vector<lydaq::TdcChannel>::iterator itchan=t->second.end();
      //hnstrip->Fill(t->second.size()-1.);

      for (std::vector<lydaq::TdcChannel>::iterator x=t->second.begin();x!=t->second.end();++x)
	{

	  

	  
	  x->setUsed(false);
	  if (x->channel()==_triggerChannel)
	    {chused[_triggerChannel]=1;trigtime=x->tdcTime();trigbcid=x->bcid()*200;x->setUsed(true);ittrig=x;
	    }
	  else
	    {
	      if (tev0==0 && x->side()==0 ) {tev0=x->tdcTime();itchan=x;}
	      //hstript->Fill(x->channel()*1.+0.5);
	    }
	}
      std::bitset<32> spat;
      spat.reset();
      std::bitset<64> spatb;
      spatb.reset();
      std::bitset<32> spat2;
      spat2.reset();
      for (std::vector<lydaq::TdcChannel>::iterator x=t->second.begin();x!=t->second.end();++x)
	{
	  INFO_PRINTF("\t TDC %d LEMO %d STRIP %d  SIDE %d   time %f \n",x->channel(),x->lemo(),x->strip(),x->side(),  (x->tdcTime()-trigtime));
	  if (x->channel()!=_triggerChannel) spat.set(x->strip()-70,1);
	  if (x->channel()!=_triggerChannel) spatb.set(32*x->side()+x->strip()-70,1);
	}
      std::cout<<_event<<" "<<_dif<<" "<<spat<<" "<<spat.count()<<std::endl;
      
      hnstrip->Fill(spat.count()*1.);
      for (int i=0;i<64;i++)
	if (spatb[i])
	  {
	    hstriptb->Fill(i+0.5);
	  }
      for (int i=0;i<32;i++)
	if (spat[i])
	  {
	    hstript->Fill(i+0.5);
	    uint8_t str=i+70;
	    double t0=-1000,t1=-1000,sh0=0,sh1=0;
	    
	    for (std::vector<lydaq::TdcChannel>::iterator x=t->second.begin();x!=t->second.end();++x)
	      {
		if (x->strip()==str && x->side()==0 && t0<0)
		  {
		    t0=x->tdcTime()+fe1_2tr[x->channel()];
		    //sh0=fe1_2tr[x->channel()];
		    std::stringstream s;
		    s<<"Timing/hdt2tr"<<(int) x->channel();
		    TH1* hdts=rh()->GetTH1(sr.str()+s.str());
		    if (hdts==NULL)
		      {
			hdts=rh()->BookTH1(sr.str()+s.str(),120,-655.,-445.);
		      }
		    hdts->Fill(t0-trigtime);
		  }
		if (x->strip()==str && x->side()==1 && t1<0 )
		  {
		    t1=x->tdcTime()+fe1_2tr[x->channel()];
		    //sh1=fe1_2tr[x->channel()];
		    std::stringstream s;
		    s<<"Timing/hdt2tr"<<(int) x->channel();
		    TH1* hdts=rh()->GetTH1(sr.str()+s.str());
		    if (hdts==NULL)
		      {
			hdts=rh()->BookTH1(sr.str()+s.str(),120,-655.,-445.);
		      }
		    hdts->Fill(t1-trigtime);

		  }
	      }
	    if (t0>0 && t1>0)
	      {
		lmana::TdcStrip ts(_dif,str+FEB2STRIP[_dif],t0,t1,fe1_shift[str]);
		_strips.push_back(ts);
		std::cout<<"Adding strip "<<str+FEB2STRIP[_dif]<<std::endl;
		spat2.set(i,1);
		//fe1_shift[str]=sh0-sh1;
		hpost->Fill(str,t0-t1-fe1_shift[str]);
		if (spat.count()<8)
		  hdt->Fill(t0-t1-fe1_shift[str]);
		std::stringstream s;
		s<<"Timing/hdtpos"<<(int) str;
		TH1* hdts=rh()->GetTH1(sr.str()+s.str());
		if (hdts==NULL)
		  {
		    hdts=rh()->BookTH1(sr.str()+s.str(),300,-25.,25.);
		  }
		//		 hdts->Fill(t0-t1-fe2_shift[str]);
		if (spat.count()<12)
		  hdts->Fill(t0-t1-fe1_shift[str]);
		 

	      }
	  }

      std::cout<<_dif<<" "<<spat2<<" "<<spat2.count()<<std::endl;
      hnstrip2->Fill(spat2.count()*1.);
      for (int i=0;i<32;i++)
	if (spat2[i])
	  {
	    hstript2->Fill(i+0.5);
	  }
	else
	  if (spat[i])
	    hstript1->Fill(i+0.5);


      if (spat2.count()==1)
	for (int i=0;i<32;i++)
	  if (spat2[i])
	    {

	      uint8_t str=i+70;
	      double t0=-1000,t1=-1000,sh0=0,sh1=0;
	    
	      for (std::vector<lydaq::TdcChannel>::iterator x=t->second.begin();x!=t->second.end();++x)
		{
		  if (x->strip()==str && x->side()==0 && t0<0)
		    {
		      t0=x->tdcTime()+fe1_2tr[x->channel()];
		    }
		  if (x->strip()==str && x->side()==1 && t1<0 )
		    {
		      t1=x->tdcTime()+fe1_2tr[x->channel()];

		    }
		}
	      if (t0>0 && t1>0)
		{
		  std::stringstream s1;
		  s1<<"Timing/hdtone"<<(int) str;
		  TH1* hdts1=rh()->GetTH1(sr.str()+s1.str());
		  if (hdts1==NULL)
		    {
		      hdts1=rh()->BookTH1(sr.str()+s1.str(),300,-25.,25.);
		    }
		  //		 hdts->Fill(t0-t1-fe2_shift[str]);

		  hdts1->Fill(t0-t1-fe1_shift[str]);

		}
	    }

      //getchar();
      uint32_t nch=0;
      for (std::vector<lydaq::TdcChannel>::iterator it=t->second.begin();it!=t->second.end();++it)
	
	{

	  if ((it->tdcTime()-trigtime<-600) &&
	      (it->tdcTime()-trigtime>-500) )
	    {
	      for (std::vector<lydaq::TdcChannel>::iterator jt=t->second.begin();jt!=t->second.end();++jt)
		if ((jt->tdcTime()-trigtime<-600) &&
		    (jt->tdcTime()-trigtime>-500) )
		  {
		    if (it->channel()==jt->channel()) continue;
		    if (it->side()==jt->side()) continue;
		    hcor->Fill(it->lemo()+1.,jt->lemo()+1.);
		  }
	    }
	  
	  // stringstream s;
	  // s<<"hdco"<<(int) it->channel();

	  
	  // TH1* hdco=rh()->GetTH1(sr.str()+s.str());
	  // if (hdco==NULL)
	  //   {
	  //     hdco=rh()->BookTH1(sr.str()+s.str(),65536,-32767,32767);
	  //   }
	  // if (it->channel()!=itchan->channel())
	  //   {
	  //     hdco->Fill(it->coarse()-itchan->coarse());
	  //   }

	  
	  if (abs(it->tdcTime()-tev0)>5000) it->setUsed(true);
	  if (it->used()) continue;
	  //DEBUG_PRINTF("%d %u %u %u %f \n",x->channel(),x->coarse(),x->fine(),x->bcid(),x->tdcTime()-tev0);
	  nch++;
	}
      //      getchar();
      if (itchan==t->second.end()) continue;
      if (ittrig==t->second.end()) continue;
      DEBUG_PRINTF("Trigtime %f %u %d tev0 %f %f %u %d \n",trigtime,ittrig->coarse(),ittrig->fine(),tev0,trigtime-tev0,itchan->coarse(),itchan->fine());
      //getchar();
      hdtr0->Fill(tev0-trigtime);
      hnst->Fill(nch*1.);
      if (nch>=1)  heff->Fill(2.1);
      if (nch>=1)  heff->Fill(4.1);
      //DEBUG_PRINTF(" Effective TDC %d  GTC %d   Number %d \n",mezId,_gtc,t->second.size());
      if (t->second.size()>2000) continue; // Use trigger with less than  20 strip

    }
  //if (mezId==15 ) getchar();
  //getchar();
  _ntrigger++;
  //heff->Fill(1.1);
  if (!found) return;
  _nfound++;
  // update efficency
  //  heff->Fill(2.1);
  if (bside) {_nbside++;heff->Fill(3.1);}
  DEBUG_PRINTF("%d-%d %d  #evt %d #dif %d #trig %d #found %d  #time %d \n",_run,_event,_gtc,_nevt,mezId,_ntrigger,_nfound,_nbside); 
}


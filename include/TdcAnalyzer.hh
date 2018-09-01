#ifndef TDC_ANALYZER_HH
#define TDC_ANALYZER_HH
#include "TdcChannel.hh"
#include "DCHistogramHandler.hh"
#include "TdcMapping.hh"
#include "jsonGeo.hh"

namespace lydaq
{
  class TdcStrip
  {
  public:
    TdcStrip() :_ch(0),_dif(0), _str(0),_t0(0),_t1(0),_shift(0){;}
    TdcStrip(uint16_t dif,uint16_t st,double t0,double t1,double shift=0) :_dif(dif), _str(st),_t0(t0),_t1(t1),_shift(shift),_ch(1) {;}
    TdcStrip(uint16_t ch,uint16_t dif,uint16_t st,double t0,double t1,double shift=0) :_ch(ch),_dif(dif), _str(st),_t0(t0),_t1(t1),_shift(shift) {;}
    inline uint16_t strip() const {return _str;}
    inline uint16_t chamber() const {return _ch;}
    inline uint16_t dif() const {return _dif;}
    inline double t0() const {return _t0;}
    inline double t1() const {return _t1;}
    inline double shift() const {return _shift;}
    #ifdef SMALLPCB
    inline double ypos() const {return (_t0-_t1-_shift)/0.125;}
    inline double xpos() const {
      if (_dif%2==1) return -1*(_str*0.4+0.2);
      else return +1*(_str*0.4+0.2);
    }
    #else
    inline double ypos() const {return (_t0-_t1-_shift)/1.;}
    inline double xpos() const {
      return _str*1.0;
    }

    #endif
  private:
    uint16_t _dif,_str,_ch;
    double _t0,_t1,_shift;

  };
  class TdcCluster
  {
  public:
    TdcCluster(): _x(0),_y(0) { _strips.clear();}
    bool isAdjacent(TdcStrip& s,float step=2)
    {

      for (auto x:_strips)
	{
	//if (abs(x.xpos()-s.xpos())<step && abs(x.ypos()-s.ypos())<2)
	float dta=2.5;
	if (x.dif()!=s.dif()) dta=3*2.5;
	if (abs(x.xpos()-s.xpos())<step && abs((x.t0()+x.t1())/2-(s.t0()+s.t1())/2)<dta && abs(x.ypos()-s.ypos())<1.5)
	  {
	    return true;
	  }
	}
      return false;
    }
    void addStrip(TdcStrip& s){
	_strips.push_back(s);
	this->calcpos();}
    void calcpos()
    {
      if (_strips.size()==1)	{_x=_strips[0].xpos(); _y=_strips[0].ypos();}
      if (_strips.size()==2){_x=(_strips[0].xpos()+_strips[1].xpos())/2.;_y=(_strips[0].ypos()+_strips[1].ypos())/2.;}
      if (_strips.size()==3)       if (_strips.size()==1)	{_x=_strips[1].xpos(); _y=_strips[1].ypos();}
      if (_strips.size()==4){_x=(_strips[2].xpos()+_strips[1].xpos())/2.;_y=(_strips[2].ypos()+_strips[1].ypos())/2.;}
      if (_strips.size()>=5 && _strips.size()<=12)
	{
	  _x=0;_y=0;
	  for (int i=2;i<_strips.size()-2;i++)
	    {
	      _x+=_strips[i].xpos();
	      _y+=_strips[i].ypos();
	    }
	  _x/=(_strips.size()-4);
	  _y/=(_strips.size()-4);

	}
    }
    double X(){return _x;}
    double Y(){return _y;}
    uint32_t size(){return _strips.size();}
    lydaq::TdcStrip& strip(int n) {return _strips[n];}
  private:
    double _x,_y;
    std::vector<lydaq::TdcStrip> _strips;

  };
  class TdcAnalyzer {
  public:
    TdcAnalyzer(DCHistogramHandler* r);
    void pedestalAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void scurveAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void normalAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void LmAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void fullAnalysis(std::vector<lydaq::TdcChannel>& vChannel);
    void multiChambers(std::vector<lydaq::TdcChannel>& vChannel);
    bool noiseStudy(std::vector<lydaq::TdcChannel>& vChannel,std::string stubdir="InTime");
    void drawHits(int nch);
    void end();
    void setInfo(uint32_t dif,uint32_t run,uint32_t ev,uint32_t gt,uint64_t ab,uint16_t trgchan,uint32_t vth,uint32_t dac);
    double acquisitionTime(){ return (_abcid-_abcid0)*2E-7;}
    void clear(){_strips.clear();}
    std::vector<lydaq::TdcStrip>& strips(){return _strips;}
    uint32_t triggers(){return _ntrigger;}
    bool trigger(){return _triggerFound;}
    uint32_t gtc(){return _gtc;}
    uint64_t abcid(){return _abcid;}
    jsonGeo* geometry(){return _geo;}
    void setGeometry(jsonGeo* g){_geo=g;}
  private:
    DCHistogramHandler* _rh;
    std::vector<lydaq::TdcChannel>::iterator _trigger;
    std::vector<lydaq::TdcStrip> _strips;
    uint32_t _dif,_run,_event,_gtc,_vthSet,_dacSet,_nevt,_ntrigger,_nfound,_nbside;
    uint64_t _abcid,_abcid0;
    uint16_t _triggerChannel;
    bool _pedestalProcessed,_triggerFound,_display,_noise;
    jsonGeo* _geo;
  };
};
#endif

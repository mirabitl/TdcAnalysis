#ifndef TDC_ANALYZER_HH
#define TDC_ANALYZER_HH
#include "TdcChannel.hh"
#include "DCHistogramHandler.hh"
#include "TdcMapping.hh"
namespace lydaq
{
  class TdcStrip
  {
  public:
    TdcStrip() :_dif(0), _str(0),_t0(0),_t1(0),_shift(0){;}
    TdcStrip(uint16_t dif,uint16_t st,double t0,double t1,double shift=0) :_dif(dif), _str(st),_t0(t0),_t1(t1),_shift(shift){;}
    inline uint8_t strip() const {return _str;}
    inline uint8_t dif() const {return _dif;}
    inline double t0() const {return _t0;}
    inline double t1() const {return _t1;}
    #ifdef SMALLPCB
    inline double ypos() const {return (_t0-_t1-_shift)/0.125;}
    inline double xpos() const {
      if (_dif%2==1) return -1*(_str*0.4+0.2);
      else return +1*(_str*0.4+0.2);
    }
    #else
    inline double ypos() const {return (_t0-_t1-_shift)/1.;}
    inline double xpos() const {
      return _str;
    }

    #endif
  private:
    uint16_t _dif,_str;
    double _t0,_t1,_shift;

  };
  class TdcAnalyzer {
  public:
    TdcAnalyzer(DCHistogramHandler* r);
    void pedestalAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void scurveAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void normalAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void LmAnalysis(uint32_t mezId,std::vector<lydaq::TdcChannel>& vChannel);
    void fullAnalysis(std::vector<lydaq::TdcChannel>& vChannel);

    void end();
    void setInfo(uint32_t dif,uint32_t run,uint32_t ev,uint32_t gt,uint64_t ab,uint16_t trgchan,uint32_t vth,uint32_t dac);
    double acquisitionTime(){ return (_abcid-_abcid0)*2E-7;}
    void clear(){_strips.clear();}
    std::vector<lydaq::TdcStrip>& strips(){return _strips;}
    uint32_t triggers(){return _ntrigger;}
    bool trigger(){return _triggerFound;}
    uint32_t gtc(){return _gtc;}
    uint64_t abcid(){return _abcid;}

  private:
    DCHistogramHandler* _rh;
    std::vector<lydaq::TdcChannel>::iterator _trigger;
    std::vector<lydaq::TdcStrip> _strips;
    uint32_t _dif,_run,_event,_gtc,_vthSet,_dacSet,_nevt,_ntrigger,_nfound,_nbside;
    uint64_t _abcid,_abcid0;
    uint16_t _triggerChannel;
    bool _pedestalProcessed,_triggerFound;
  };
};
#endif

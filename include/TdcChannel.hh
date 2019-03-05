#ifndef _TDCCHANNEL_HH
#define _TDCCHANNEL_HH
#define TDC_COARSE_TIME 2.5
#define TDC_FINE_TIME 0.009765625
#include <stdint.h>
#include <stdio.h>
#include "TdcMapping.hh"
#include "jsonGeo.hh"
using namespace std;
namespace lydaq {
#define FEBCMS
#undef SUCHAIN
#ifdef SUCHAIN
  class TdcChannel
{
public:
  TdcChannel() :_fr(NULL),_used(false) {;}
  TdcChannel(uint8_t*  b,uint8_t feb=0) :_fr(b),_used(false),_feb(feb) {;}
  inline uint8_t channel() {return  (_fr[0]&0XFF);}
  #ifdef OLDMAPPING
  inline uint16_t pr() {return  TDC2PR[channel()];}
  inline uint16_t lemo() {return  PR2LEMO[pr()];}
  inline uint16_t side() {return  lemo()/16;}
  inline uint16_t strip() {return  LEMO2STRIP[lemo()];}
  inline uint16_t feb(){return _feb;}
  #else
  inline uint16_t pr() {return  TDC2PR[channel()];}
  inline uint16_t lemo() {return  strip()+12*side();}
  inline uint16_t side() {return  SIDE[channel()];}
  inline uint16_t strip() {return  70+STRIP[channel()];}
  inline uint16_t lemo(jsonFebInfo& f) {return  strip(f)+12*side(f);}
  inline uint16_t side(jsonFebInfo& f) {return  f.tdc2side[channel()];}
  inline uint16_t strip(jsonFebInfo& f) {return  70+f.tdc2strip[channel()];}
  inline uint16_t feb(){return _feb;}

  #endif
  inline uint16_t detectorStrip(uint32_t feb) {return  strip()+FEB2STRIP[feb];}
  inline uint16_t detectorStrip(jsonFebInfo& f) {return  strip(f)+f.stripShift;}
  inline uint8_t length(){return 8;}
  inline uint64_t coarse() const {return ((uint64_t)_fr[6])|((uint64_t)_fr[5]<<8)|((uint64_t)_fr[4]<<16)|((uint64_t)_fr[3]<<24);}
  inline uint8_t fine() const {return _fr[7];}
  #ifdef BCIDFROMCOARSE
  inline uint16_t bcid(){return (uint16_t) (coarse()*TDC_COARSE_TIME/200);}
  #else
  inline  uint16_t bcid() const {return (_fr[2]|(_fr[1]<<8));}
  #endif
  //inline  double tdcTime() const { return (coarse()*TDC_COARSE_TIME+fine()*TDC_FINE_TIME);}
  inline  double tdcTime() const { return (coarse()+fine()/256.0)*TDC_COARSE_TIME;}
  inline uint8_t* frame(){ return _fr;}
  inline bool used(){return _used;}
  inline void setUsed(bool t){_used=t;}
  bool operator<(const TdcChannel &ipaddr){
    if( coarse() < ipaddr.coarse())
    return true;
  else
    return false;
    
}
  void dump()
  {
    for (int i=0;i<8;i++)
      printf("%.2x ",_fr[i]);
    printf("%d %d %d %f \n",channel(),coarse(),fine(),tdcTime());
    printf("\n");
  }
private:
  uint8_t* _fr;
  bool _used;
  uint8_t _feb;
};

  #endif
#define LASTCHAN 44
  #ifdef FEBCMS
  class TdcChannel
{
public:
  TdcChannel() :_fr(NULL),_used(false) {;}
  TdcChannel(uint8_t*  b,uint8_t feb=0) :_fr(b),_used(false),_feb(feb) {;}
  inline uint8_t channel() {
    if (falling())
      return (_fr[0]-0x80)-LASTCHAN;
    else
      return _fr[0];}
  inline bool falling() {return  ((_fr[0]&0X80)>>7)==1;}
  inline uint16_t pr() {return  TDC2PR[channel()];}
  inline uint16_t lemo() {return  strip()+12*side();}
  inline uint16_t side() {return  SIDE[channel()];}
  inline uint16_t strip() {return  70+STRIP[channel()];}
  inline uint16_t lemo(jsonFebInfo& f) {return  strip(f)+12*side(f);}
  inline uint16_t side(jsonFebInfo& f) {return  f.tdc2side[channel()];}
  inline uint16_t strip(jsonFebInfo& f) {return  70+f.tdc2strip[channel()];}
  inline uint16_t feb(){return _feb;}


  inline uint16_t detectorStrip(uint32_t feb) {return  strip()+FEB2STRIP[feb];}
  inline uint16_t detectorStrip(jsonFebInfo& f) {return  strip(f)+f.stripShift;}
  inline uint8_t length(){return 6;}
  inline uint64_t coarse() const {return ((uint64_t)_fr[4])|((uint64_t)_fr[3]<<8)|((uint64_t)_fr[2]<<16)|((uint64_t)_fr[1]<<24);}
  inline uint8_t fine() const {return _fr[5];}

  inline uint32_t bcid(){return (uint32_t) (coarse()*TDC_COARSE_TIME/200);}
  inline  double tdcTime() const { return (coarse()+fine()/256.0)*TDC_COARSE_TIME;}
  inline uint8_t* frame(){ return _fr;}
  inline bool used(){return _used;}
  inline void setUsed(bool t){_used=t;}
  bool operator<(const TdcChannel &ipaddr){
    if( coarse() < ipaddr.coarse())
    return true;
  else
    return false;
    
}
  void dump()
  {
    for (int i=0;i<6;i++)
      printf("%.2x ",_fr[i]);
    printf("%d %d %d %f \n",channel(),coarse(),fine(),tdcTime());
    printf("\n");
  }
private:
  uint8_t* _fr;
  bool _used;
  uint8_t _feb;
};

  #endif
};
#endif

#ifndef _RECOFRAME_HH
#define _RECOFRAME_HH
#include "jsonGeo.hh"
#include <stdint.h>
#include <Math/Point3Dfwd.h>
#include <bitset>
#include <Math/Point3Dfwd.h>
#include <TPrincipal.h>
#include <cfloat>
using namespace ROOT::Math;
class pcaComponents
{
public:
  pcaComponents() {memset(_components,0,21*sizeof(double));}
  double& operator[](uint32_t i) { return _components[i];}
private:
  double _components[21];
};

class RecoFrame : public ROOT::Math::XYZPoint
  {
  public:
    enum Type {THR0=0,THR1=1,THR2=2,CORE=3,EDGE=4,ISOLATED=5,HOUGH=6,MIP=7};
    RecoFrame(){_frame=0;_geo=0;}
    RecoFrame(uint64_t f){_frame=f;_geo=0;}
    RecoFrame(uint64_t f,jsonGeo* g){_frame=f;_geo=g;}
    RecoFrame(uint8_t dif,uint8_t asic,uint8_t channel,uint8_t threshold,uint32_t bc){ this->set(dif,asic,channel,threshold,bc);}
    void initialise(jsonGeo* g) {_geo=g;g->convert(this->dif(),this->asic(),this->channel(),this);}
    void initialise(uint64_t f,jsonGeo* g) {_frame=f;_geo=g;g->convert(this->dif(),this->asic(),this->channel(),this);}
      
    inline uint16_t I(){return int(X()/1.04125);}
    inline uint16_t  J(){return int(Y()/1.04125);}
    inline uint16_t chamber(){return (_geo==NULL)?0:_geo->difInfo(this->dif()).chamber;}
    inline uint32_t dif(){return (_frame&0xFF);}
    inline uint32_t asic( ){return ((_frame>>8)&0x3F);}
    inline uint32_t channel( ){return ((_frame>>14)&0x3F);}
    inline uint32_t threshold(){return ((_frame>>20)&0x3);}
    inline uint32_t bc(){ return ((_frame>>32)&0xFFFFFFFF);}
    inline void set(uint8_t dif,uint8_t asic,uint8_t channel,uint8_t threshold,uint32_t bc)
    {
      uint64_t d=dif,a=asic,c=channel,t=threshold,b=bc;
      _frame=((b&0xFFFFFFFF)<<32)|((t&0x3)<<20)|((c&0x3f)<<14)|((a&0x3F)<<8)|
        (d&0xFF);}
    inline uint64_t value() const {return _frame;}
    inline void set( uint64_t f){_frame=f;}


    inline uint16_t Tag(){return _tag;}
    inline void setTag(RecoFrame::Type t,bool v){_tag=(v)?_tag|(1<<t):_tag&~(1<<t);}
    inline bool isTagged(RecoFrame::Type t){return (_tag&(1<<t))!=0;}
    inline void clearTag(){_tag=0;}
    inline jsonGeo* geo() const {return _geo;}
    inline void setWeight(float t)  { _weight=t;}
    inline float weight() const {return _weight;}
    inline void setUse(bool t)  { _used=t;}
    inline bool isUsed() const {return _used;}
      template <class T>
  static pcaComponents calculateComponents(std::vector<T*> vnear_)
{
  pcaComponents c;

  if (vnear_.size()<3) return c;
  TPrincipal tp(3,"");
  double xp[3];
  uint32_t nh=0;
  double xb=0,yb=0,zb=0;
  double wt=0.;
  
  double fp=DBL_MAX;
  double lp=-DBL_MAX;
  double fx=DBL_MAX;
  double lx=-DBL_MAX;
  double fy=DBL_MAX;
  double ly=-DBL_MAX;


  //INFO_PRINT("%d vector size\n",v.size());
  for (typename std::vector<T*>::iterator it=vnear_.begin();it!=vnear_.end();it++)
    {
      T* iht=(*it);
 
      double w=1.;
      xb+=iht->X()*w;
      yb+=iht->Y()*w;
      zb+=iht->Z()*w;
      wt+=w;
      xp[0]=iht->X();
      xp[1]=iht->Y();
      xp[2]=iht->Z();
      tp.AddRow(xp);
      if (iht->Z()<fp) fp=iht->Z();
      if (iht->Z()>lp) lp=iht->Z();
      if (iht->X()<fx) fx=iht->X();
      if (iht->X()>lx) lx=iht->X();
      if (iht->Y()<fy) fy=iht->Y();
      if (iht->Y()>ly) ly=iht->Y();
      nh++;
    }
  
  if (nh<3) return c;
  tp.MakePrincipals();
  c[0]=(*tp.GetMeanValues())[0];
  c[1]=(*tp.GetMeanValues())[1];
  c[2]=(*tp.GetMeanValues())[2];



  c[3]=(*tp.GetEigenValues())[2];
  c[4]=(*tp.GetEigenValues())[1];
  c[5]=(*tp.GetEigenValues())[0];


  c[6]=(*tp.GetEigenVectors())(0,0);
  c[7]=(*tp.GetEigenVectors())(1,0);
  c[8]=(*tp.GetEigenVectors())(2,0);	
  c[9]=(*tp.GetEigenVectors())(0,1);
  c[10]=(*tp.GetEigenVectors())(1,1);
  c[11]=(*tp.GetEigenVectors())(2,1);
  c[12]=(*tp.GetEigenVectors())(0,2);
  c[13]=(*tp.GetEigenVectors())(1,2);
  c[14]=(*tp.GetEigenVectors())(2,2);

  // Store First and last Z
  c[15]=fp;
  c[16]=lp;
  c[17]=fx;
  c[18]=lx;
  c[19]=fy;
  c[20]=ly;
  
  return c;
}
  private:
    uint64_t _frame;
    jsonGeo* _geo;
    uint16_t _tag;
    float _weight;
    bool _used;
  };
#endif

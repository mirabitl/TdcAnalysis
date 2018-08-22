import sqlite3
from ROOT import *
import math
import time
from array import array

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self


defped=[31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31]
#defped=[30, 30, 34, 29, 32, 32, 37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 33, 27, 31, 29, 25, 35, 33, 9]
#defped=[39,38,42,37,41,40,46,0,43,48,37,49,0,0,0,0,0,0,0,0,39,39,36,35,41,35,40,34,32,44,43,34]

defped=[38,38,41,39,41,40,45,0,42,47,37,47,0,0,0,0,0,0,0,0,38,38,35,36,42,36,40,36,33,43,42,35]
def calP(p0,alt):
  return p0*(1-0.0065*alt/288.15)**5.255

def calV(V,P,T):
  print 1-(0.2+0.8*P/990.*293./T)
  return  V/(0.2+0.8*P/990.*293./T)

def calValice(V,P,T):
  print 1-(1*P/990.*293./T)
  return  V/(1.0*P/990.*293./T)

def calApp(V,P,T):
  print 1-(0.2+0.8*P/990.*293./T)
  return  V*(0.2+0.8*P/990.*293./T)

def getdt(run,chamber,feb):
  f82=TFile("./histo%d_0.root" % run);
  r=[]
  for i in range(0,25):
    r.append(0)
  for i in range(1,25):
    #print "/run%d/Chamber%d/FEB/%d/Side0/channel%d" % (run,chamber,feb,i)
    f82.cd("/run%d/Chamber%d/FEB/%d/Side0/" % (run,chamber,feb))
    hch=f82.Get("/run%d/Chamber%d/FEB/%d/Side0/channel%d" % (run,chamber,feb,i))
    #print hch
    if (hch!=None):
      #print i,hch.GetEntries(),hch.GetMean();
      if (hch.GetEntries()>5):
        r[i]=hch.GetMean()
      continue
    f82.cd("/run%d/Chamber%d/FEB/%d/Side1/" % (run,chamber,feb))
    hch1=f82.Get("/run%d/Chamber%d/FEB/%d/Side1/channel%d" % (run,chamber,feb,i))
    #print hch1
    if (hch1!=None):
      #print i,hch1.GetEntries(),hch1.GetMean();
      if (hch1.GetEntries()>5):
        r[i]=hch1.GetMean()
      continue

  #print r
  r = map(prettyfloat, r)
  print '"dtc":',r
def fitProfile(run,sel=92):
  f82=TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/ChamberAll" % (run));
  c1=TCanvas();
  gStyle.SetOptFit();
  pos0=[]
  pmean0=[]
  pos1=[]
  pmean1=[]
  hall=f82.Get("/run%d/ChamberAll/XY" % (run));
  for i in range(sel,sel+1):
    pos0.append(0)
    pmean0.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=hall.ProjectionY("strip%d" % (i),(i),i+1)
    if (hstrip.GetEntries()<100):
      continue
    scfit=TF1("scfit","gaus",hstrip.GetMean()-2.*hstrip.GetRMS(),hstrip.GetMean()+2.*hstrip.GetRMS())
    hstrip.Fit("scfit","Q","",hstrip.GetMean()-2.*hstrip.GetRMS(),hstrip.GetMean()+2.*hstrip.GetRMS())
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    print run,i,dtmean,dtres,dtmean*80./8.3,dtres*80/8.3
    c1.Update()
    #val = raw_input()
    time.sleep(2)


def fitdif(run):
  f82=TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/ChamberDif" % (run));
  c1=TCanvas();
  gStyle.SetOptFit();
  pos0=[]
  pmean0=[]
  pos1=[]
  pmean1=[]
  for i in range(48):
    pos0.append(0)
    pmean0.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=f82.Get("/run%d/Timing/OneStrip/Side0/hdtr_%d" % (run,i+71));
    if (hstrip == None):
      continue
    if (hstrip.GetEntries()<20):
      continue
    hstrip.Rebin(2)
    
    hstrip.Draw()
    c1.Update()

    print "Enter min max"
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=hstrip.GetMean()-5.*hstrip.GetRMS()
    hmax=hstrip.GetMean()+5.*hstrip.GetRMS()
    print hmin,hmax
    #scfit=TF1("scfit","gaus",hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    #hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    
    scfit=TF1("scfit","gaus",hmin,hmax)
    hstrip.GetXaxis().SetRangeUser(hmin,hmax)
    hstrip.Fit("scfit","Q");
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    hstrip.Draw()
    c1.Update()
    #c1.SaveAs("Run%d_Strip_pos.png" % (run));

    val = raw_input()
    pmean0[i]=hstrip.GetMean()
    if (dtres<hstrip.GetRMS()):
      pos0[i]=dtmean
    else:
      pos0[i]=hstrip.GetMean()
      
  print pos0
  print pmean0
  for i in range(48):
    pos1.append(0)
    pmean1.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=f82.Get("/run%d/Timing/OneStrip/Side1/hdtr_%d" % (run,i+71));
    if (hstrip == None):
      continue
    if (hstrip.GetEntries()<20):
      continue
    hstrip.Rebin(2)
    
    hstrip.Draw()
    c1.Update()

    print "Enter min max"
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=hstrip.GetMean()-5.*hstrip.GetRMS()
    hmax=hstrip.GetMean()+5.*hstrip.GetRMS()
    print hmin,hmax
    #scfit=TF1("scfit","gaus",hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    #hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    
    scfit=TF1("scfit","gaus",hmin,hmax)
    hstrip.GetXaxis().SetRangeUser(hmin,hmax)
    hstrip.Fit("scfit","Q");
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    hstrip.Draw()
    c1.Update()
    #c1.SaveAs("Run%d_Strip_pos.png" % (run));

    val = raw_input()
    pmean1[i]=hstrip.GetMean()
    if (dtres<hstrip.GetRMS()):
      pos1[i]=dtmean
    else:
      pos1[i]=hstrip.GetMean()
      
  print pos0
  print pmean0





  
  dt=0
  for i in range(48):
    print "fe1_2tr[%d]=%5.3f;" % (i+71,pos0[i]-pos1[i])

  #print pos

def fitpos(run,tdc):
  f82=TFile("./histo%d_0.root" % run);
  #f82.cd("/run%d/TDC%d/LmAnalysis/Timing" % (run,tdc));
  f82.cd("/run%d/Chamber%d/Timing" % (run,tdc));
  c1=TCanvas();
  gStyle.SetOptFit();
  pos=[]
  pmean=[]
  for i in range(48):
    pos.append(0)
    pmean.append(0)
    #hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/Timing/hdtpos%d" % (run,tdc,i+71));
    hstrip=f82.Get("/run%d/Chamber%d/Timing/All/hdtpos%d" % (run,tdc,i+71));
    if (hstrip == None):
      continue
    hstrip.Rebin(2)
    
    hstrip.Draw()
    c1.Update()

    print "Enter min max"
    #hmin = float(raw_input())
    #hmax = float(raw_input())
    hmin=-20.
    hmax=20.
    print hmin,hmax
    #scfit=TF1("scfit","gaus",hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    #hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-3.*hstrip.GetRMS(),hstrip.GetMean()+3.*hstrip.GetRMS())
    
    scfit=TF1("scfit","gaus",hmin,hmax)
    hstrip.GetXaxis().SetRangeUser(hmin,hmax)
    hstrip.Fit("scfit","Q");
    dtmean=scfit.GetParameter(1)
    dtres=scfit.GetParameter(2)
    hstrip.Draw()
    c1.Update()
    #c1.SaveAs("Run%d_Strip_pos.png" % (run));

    val = raw_input()
    pmean[i]=hstrip.GetMean()
    if (dtres>1. and  dtres<4):
      pos[i]=dtmean
    else:
      pos[i]=hstrip.GetMean()
  print pos
  print pmean
  
  for i in range(12):
    print "fe1_2tr[%d]=%5.3f;" % (i,pos[i])
  #print pos

def calcefn(run,chamber,hv=0,dirp="."):
  f82=TFile("%s/histo%d_0.root" % (dirp,run));
  f82.cd("/run%d/Chamber%d/" % (run,chamber));
  c1=TCanvas();
  gStyle.SetOptFit();

  hns=f82.Get("/run%d/Chamber%d/Efficiency" % (run,chamber));
  #hns2=f82.Get("/run%d/TDC%d/LmAnalysis/NStrips2" % (run,tdc));
  hstrip=f82.Get("/run%d/Chamber%d/NSTRIP" % (run,chamber));
  hrate=f82.Get("/run%d/Chamber%d/Rate" % (run,chamber));

  hclusters=f82.Get("/run%d/Chamber%d/ClusterNew/Clusters" % (run,chamber));
  hclusterm=f82.Get("/run%d/Chamber%d/ClusterNew/ClusterSize1" % (run,chamber));
  #hstrip.Rebin(2)
  csize=0.1
  effc=0.0
  deffc=0.0
  ncev=0
  ncl=0
  if (hclusterm!=None):
    csize=hclusterm.GetMean()
    ncev=hclusters.GetEntries()
    nc=ncev-hclusters.GetBinContent(1)
    ncl=0
    nevcl=0
    for i in range(2,16):
        x=i-1.;
        y=hclusters.GetBinContent(i)
        nevcl=nevcl+y
        ncl=ncl+x*y
    if (nevcl>0):
        ncl=ncl*1./nevcl
    effc=nc*1./ncev
    deffc=math.sqrt(effc*(1-effc)/ncev)
  ntrg=hns.GetBinContent(3)
  nall=hns.GetBinContent(4)
  nxy=hns.GetBinContent(5)

  hstrip.GetXaxis().SetRangeUser(0.5,9.5)
  hrate.GetXaxis().SetRangeUser(0.1,60000.)
  mul=hstrip.GetMean()

  
  #print ntrg,hns.GetBinContent(1),hns2.GetBinContent(1)
  eff=nall*1./ntrg
  deff=math.sqrt(eff*(1-eff)/ntrg)
  effp=nxy*1./ntrg
  deffp=math.sqrt(effp*(1-effp)/ntrg)
  
  print "|%d|%d|%7.1f|%d|%d|%d|%5.2f|%5.2f|%5.2f|%5.2f|%5.1f|%7.1f|%d|%5.2f|%5.2f|%5.2f|%5.2f|" % (run,chamber,hv,int(ntrg),int(nall),int(nxy),eff*100,deff*100,effp*100,deffp*100,mul,hrate.GetMean(),ncev,effc*100,deffc*100,csize,ncl)
  #hstrip.Draw()
  #c1.Update()
  #c1.SaveAs("Run%d_Strip_pos.png" % (run));

  #val = raw_input()
  r=(run,chamber,hv,int(ntrg),int(nall),int(nxy),eff*100,deff*100,effp*100,deffp*100,mul,hrate.GetMean(),ncev,effc*100,deffc*100,csize,ncl)
  #r = map(prettyfloat, r)
  return r

def calceff(run,tdc,strip=71):
  f82=TFile("./histo%d_0.root" % run);
  f82.cd("/run%d/TDC%d/LmAnalysis" % (run,tdc));
  c1=TCanvas();
  gStyle.SetOptFit();

  heff=f82.Get("/run%d/TDC%d/LmAnalysis/Efficiency" % (run,tdc));
  hstrip=f82.Get("/run%d/TDC%d/LmAnalysis/hdts%d" % (run,tdc,strip));
  hstrip.Rebin(2)
  scfit=TF1("scfit","gaus",hstrip.GetMean()-6.*hstrip.GetRMS(),hstrip.GetMean()+6.*hstrip.GetRMS())
  hstrip.GetXaxis().SetRangeUser(hstrip.GetMean()-6.*hstrip.GetRMS(),hstrip.GetMean()+6.*hstrip.GetRMS())
  hstrip.Fit("scfit","Q");
  dtmean=scfit.GetParameter(1)
  dtres=scfit.GetParameter(2)

  ntrg=heff.GetBinContent(2)
  n1=heff.GetBinContent(3)
  n2=heff.GetBinContent(4)
  eff=n1/ntrg
  deff=math.sqrt(eff*(1-eff)/ntrg)
  effp=n2/ntrg
  deffp=math.sqrt(effp*(1-effp)/ntrg)
  
  print "|%d|150|%d|%d|%d|%5.2f|%5.2f|%5.2f|%5.2f|%5.2f|%5.3f|" % (run,int(ntrg),int(n1),int(n2),eff*100,deff*100,effp*100,deffp*100,dtmean,dtres)
  hstrip.Draw()
  c1.Update()
  c1.SaveAs("Run%d_Strip%d_pos.png" % (run,strip));

  val = raw_input()

def fitped(run,tdc,vthmin,vthmax,old=defped):
  ped=[]
  for i in range(32):
    ped.append(0)
  f82=TFile("./histo%d_0.root" % run);
  f82.cd("/run%d/TDC%d" % (run,tdc));
  c1=TCanvas();
  #c2=TCanvas("c2","Test",1400,900);
  #c2.cd()
  #c2.Divide(6,4)
  #c2.Draw()
  #c2.Update()
  #val = raw_input()
  #c2.Draw()
  fout=open("summary_pedestal_%d_tdc%d.txt" % (run,tdc),"w");
  fout.write("+--+-----+-----+-----+ \n");
  gStyle.SetOptFit();
  hmean=TH1F("hmean","Summary %d %d " %(run,tdc),vthmax-vthmin+1,vthmin,vthmax)
  hnoise=TH1F("hnoise","Summary noise %d %d " %(run,tdc),100,0.,30.)
  hpmean=TH1F("hpmean","Summary %d %d " %(run,tdc),32,0.,32.);
  hpnoise=TH1F("hpnoise","Summary noise %d %d " %(run,tdc),32,0.,32.);
  scfit=TF1("scfit","[0]*TMath::Erfc((x-[1])/[2])",vthmin+1,vthmax);
  
  for ip in range(1,25):
      #c2.cd()
      hs=f82.Get("/run%d/TDC%d/vthc%d" % (run,tdc,ip));
      if (hs.GetEntries()==0):
        continue
      nmax=0
      for i in range(1,hs.GetNbinsX()):
        if (hs.GetBinContent(i)==0):
              if (hs.GetBinContent(i-1)!=0 and hs.GetBinContent(i+1)!=0):
                  hs.SetBinContent(i,(hs.GetBinContent(i-1)+hs.GetBinContent(i+1))/2.)
        else:
          if (hs.GetBinContent(i)>nmax):
              nmax=hs.GetBinContent(i)

      hs.GetXaxis().SetRangeUser(vthmin-1,vthmax);
      icolor= ip%4 +1
      istyle= ip/4+1
      hs.SetLineColor(icolor)
      hs.SetLineStyle(istyle)
      hs.SetLineWidth(2)
      c1.cd()
      c1.Draw()
      



      if (ip==0):
        hs.Draw()
      else:
        hs.Draw("SAME")
  c1.Update()
  c1.SaveAs("Run%d_AllStrip%d.root" % (run,tdc));
  c1.SaveAs("Run%d_AllStrip%d.png" % (run,tdc));

  val = raw_input()
      
  for ip in range(0,25):
      #c2.cd()
      hs=f82.Get("/run%d/TDC%d/vthc%d" % (run,tdc,ip));
      if (hs.GetEntries()==0):
        continue
      hder=TH1F("hder%d" % ip,"derivative",900,0.,900.)	
      #hs.Rebin(2)
      nmax=0
      for i in range(1,hs.GetNbinsX()):
          if (hs.GetBinContent(i)==0):
              if (hs.GetBinContent(i-1)!=0 and hs.GetBinContent(i+1)!=0):
                  hs.SetBinContent(i,(hs.GetBinContent(i-1)+hs.GetBinContent(i+1))/2.)
          else:
            if (hs.GetBinContent(i)>nmax):
              nmax=hs.GetBinContent(i)
      for i in range(1,hs.GetNbinsX()):
        if (hs.GetBinContent(i)-hs.GetBinContent(i+1)>-10):
          hder.SetBinContent(i,hs.GetBinContent(i)-hs.GetBinContent(i+1))
      hder.Rebin(4)
      hder.GetXaxis().SetRangeUser(vthmin-1,vthmax);
      scfit.SetParameter(0,nmax/2.);
      scfit.SetParameter(1,hder.GetMean());
      scfit.SetParameter(2,hder.GetRMS());



      hs.GetXaxis().SetRangeUser(vthmin-1,vthmax);
      hs.Fit("scfit","Q","",vthmin+2,vthmax);
      #hs.GetXaxis().SetRangeUser(vthmin-1,scfit.GetParameter(1)+60);
      #gPad.SetLogy();
      rped=scfit.GetParameter(1)
      c1.cd()
      c1.Draw()
      

      hder.Draw()
      c1.Update()
      #val = raw_input()

      print "heho ",rped,hder.GetMean()
      rped=hder.GetMean()
      hs.Draw()
      

      c1.cd()
      c1.Draw()
      c1.Update()

      fout.write("|%2d|%5.1f|%5.1f|%5.2f| \n" % (ip,scfit.GetParameter(0),rped,scfit.GetParameter(2)));
      ipr=0
      if (ip%2==1):
        ipr=ip/2
      else:
        ipr=31-ip/2
      firmwaret=[31,29,27,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6]
      #firmwaco=[31,29,27,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4]
      firmware2=[24,5,3,1,0,2,4,6,7,8,9,10,26,28,30,31,29,27,25,23,22,21,20,19]
      firmware=firmware2
      if (ip>0):
        ipr=firmware[ip-1]
      else:
        ipr=0
      ped[ipr]=rped
      print ip,ipr,ped[ipr]
      hmean.Fill(rped)
      hnoise.Fill(scfit.GetParameter(2))
      hpmean.SetBinContent(ip+1,rped);
      hpnoise.SetBinContent(ip+1,scfit.GetParameter(2))
      #c1.SaveAs("Run%d_Strip%d.root" % (run,ip));
      val = raw_input()

      #hder.Draw()
      
      #c1.Update()
      #val = raw_input()
  c1.cd()
  hmean.Draw()
  hpmean.GetYaxis().SetRangeUser(vthmin,vthmax)
  hpmean.Draw()
  c1.Update()
  c1.SaveAs("Summary_%d_TDC%d.png" % (run,tdc));
  val = raw_input()
  hnoise.Draw()
  hpnoise.Draw()
  c1.Update()
  c1.SaveAs("Summary_Noise_%d_TDC%d.png" % (run,tdc));

  fout.write("+--+-----+-----+-----+ \n");
  fout.close()
  print ped
  val = raw_input()
  med=5550.0
  for i in range(32):
    if (ped[i]==0):
      continue;
    if (ped[i] < med):
      #print med,ped[i]
      med=ped[i]

  med=med+5
  med=480
  print "Alignment to :",med
  dac=ped
  for i in range(32):
    if (ped[i]==0):
      continue;
    old[i]=0
    dac[i]=int(round(old[i]+(med-ped[i])*1./2.97))
  print "cor%d=" % tdc,dac
  return dac

def calcped(oldpr,ped,median):
  print "Alignment to :",median
  dac=[]
  for i in range(32):
    dac.append(0)

  for i in range(32):
    if (ped[i]<=0):
      continue;
    dac[i]=int(round(oldpr[i]+(median-ped[i])*1./2.97))
  print dac
  return dac
import os
def process(runs,proc=True):
  if (proc):  
    for run in runs:
      os.system("./bin/tdcr -ggifpp_geom.json -r%d " % run)

  v=6800
  for run in runs:
    calcefn(run,1,v)
    v=v+200

  v=6800
  for run in runs:
    calcefn(run,2,v)
    v=v+200

def proclist(first,last,proc=True,vf=6700,step=100):
  if (proc):  
    for run in range(first,last+1):
      os.system("./bin/tdcr -ggifpp_geom.json -r%d " % run)

  v=vf
  for run in range(first,last+1):
    calcefn(run,1,v)
    v=v+step

  v=vf
  for run in range(first,last+1):
    calcefn(run,2,v)
    v=v+step

def processDCS(fdb,webdcs,proc=True,diro="."):
    conn = sqlite3.connect(fdb)
    conn.text_factory = str
    curs = conn.cursor()
    curs.execute("SELECT RUN,HV,ATT,DEAD,TRET,POSITION FROM runs WHERE DCS=%d" % webdcs)
    v=curs.fetchall()
    #print v[0]
    #print v[0][0]
    #return
    if (proc):  
        for x in v:
            os.system("./bin/tdcr -ggifpp_geom.json -r%d " % x[0])
            
    res=[]        
    for x in v:
        res.append(calcefn(x[0],1,x[1],dirp=diro))
    for x in v:
        res.append(calcefn(x[0],2,x[1],dirp=diro))
    return res

def storeResults(fdbi,fdbo,webdcs,histod="."):
    conn = sqlite3.connect(fdbo)
    conn.text_factory = str
    curs = conn.cursor()
    res=processDCS(fdbi,webdcs,False,diro=histod)
    #run,chamber,hv,int(ntrg),int(nall),int(nxy),eff*100,deff*100,effp*100,deffp*100,mul,hrate.GetMean(),ncev,effc*100,deffc*100,csize,ncl]
    for r in res:
        print len(r)
        sql_ins='''INSERT INTO INTIME(run,chamber,hv,ntrg,nall,nxy,eff,deff,effp,deffp,mul,rate,ncev,effc,deffc,csize,ncl,dcs) VALUES(%d,%d,%d,%d,%d,%d,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%7.1f,%d,%5.2f,%5.2f,%5.2f,%5.2f,%d)''' % (r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10],r[11],r[12],r[13],r[14],r[15],r[16],webdcs)
        #print sql_ins
        curs.execute(sql_ins)
    conn.commit()
    
def buildTGraph(title,vx,dvx,vy,dvy,tx,ty):
    r_ti, r_p = array( 'd' ), array( 'd' )
    er_ti, er_p = array( 'd' ), array( 'd' )
    for x in vx:
      r_ti.append(x)
    for y in vy:
      r_p.append(y)
    for x in dvx:
      er_ti.append(x)
    for y in dvy:
      er_p.append(y)

    gr = TGraphErrors( len(vx), r_ti, r_p,er_ti,er_p )
    gr.SetLineColor( 1 )
    gr.SetLineWidth( 1 )
    gr.SetMarkerColor( 2 )
    gr.SetMarkerStyle( 21 )
    gr.SetMarkerSize( 0.4 )
    gr.SetTitle(title)
    gr.GetXaxis().SetTitle(tx )
    gr.GetYaxis().SetTitle( ty )
    return gr

def  drawDCS(fdbi,fdbo,webdcs,chamber,c=None):
    if (c==None):
        c=TCanvas()
    gStyle.SetOptFit(1)

    conn = sqlite3.connect(fdbi)
    conn.text_factory = str
    curs = conn.cursor()
    curs.execute("select DISTINCT(DCS),ATT,TRET,DEAD from runs WHERE DCS=%d" % webdcs);
    v=curs.fetchall()
    if (len(v)==0):
        return
    att=v[0][1]
    dthr=v[0][2]-500
    dead=v[0][3]
    print att,dthr,dead
    connr = sqlite3.connect(fdbo)
    connr.text_factory = str
    cursr = connr.cursor()
    cursr.execute("select hv,effc,deffc from INTIME where DCS=%d AND CHAMBER=%d" % (webdcs,chamber));
    vo=cursr.fetchall()
    if (len(vo)<6):
        return;
    hv=[]
    eff=[]
    dhv=[]
    deff=[]
    for x in vo:
        hv.append(x[0])
        dhv.append(10.)
        eff.append(x[1])
        deff.append(x[2])
    title="DCS%d_ATT%3.1f_THR%d_DT%d_CH%d" % (webdcs,att,dthr,dead,chamber)
    gr=buildTGraph("effi",hv,dhv,eff,deff,"HV effective (V)","efficiency (%)")

    func = TF1("func", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", 6500,8200)
    func.SetParameters(100, 9.E-3, 7000)
    gr.Fit(func,"","",6700,8200)
    hv95=func.GetX(func.GetParameter(0)*0.95)
    print "HV95",hv95,hv95+150
    title=title+"_HV95_%4.0f_WP_%4.0f" % (hv95,hv95+150)
    gr.SetTitle(title)
    gStyle.SetStatX(0.85)
    gStyle.SetStatY(0.7)
    c.cd()
    gr.Draw("AP")
    c.Update()
    c.SaveAs("%s.png" % title)
    val=raw_input()

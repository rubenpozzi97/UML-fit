import argparse
parser = argparse.ArgumentParser(description="")
parser.add_argument("bin", help = "choose q2 bin range", default = -1, type = float)
parser.add_argument("year", help = "choose year [format:2016, 20162017]", default = '2016')
args = parser.parse_args()

import ROOT
from ROOT import RooFit
from math import sqrt
import math
from collections import OrderedDict
ROOT.gROOT.SetBatch(True)


class fit_pars(object):
    '''
    '''
    def __init__(self):
        self.Reset()

    def Reset(self):
        self.Fl          = -10
        self.P1          = -10
        self.P2          = -10
        self.P3          = -10
        self.P4p         = -10
        self.P5p         = -10
        self.P6p         = -10
        self.P8p         = -10

    def __str__(self):
        toWrite = ' Fl  = .3%f \n \
                    P1  = .3%f \n \
                    P2  = .3%f \n \
                    P3  = .3%f \n \
                    P4p = .3%f \n \
                    P5p = .3%f \n \
                    P6p = .3%f \n \
                    P8p = .3%f \n '.format(self.Fl, self.P1, self.P2,self.P3, self.P4p, self.P5p, self.P6p, self.P8p)
        return toWrite 

    ### Physical Region: should be < 0
    def ctL4phi1(self):
      return (self.P4p*self.P4p + \
              self.P5p*self.P5p + \
              self.P6p*self.P6p + \
              self.P8p*self.P8p - \
              2 + 2*abs( 2*self.P2 - self.P4p*self.P5p +self.P6p*self.P8p ) < 0 )
             
    def ctK2(self):
        return (self.P1*self.P1 + \
                4*self.P2*self.P2 + \
                4*self.P3*self.P3 -1  < 0 )

    def ctL2phi1(self):
        return (self.P5p*self.P5p*(1-self.P1) + \
                self.P6p*self.P6p*(1+self.P1) - \
                4*self.P3*self.P5p*self.P6p - \
                1 + self.P1*self.P1 + 4*self.P3*self.P3 < 0)
    def ctL2phi2(self):      
        return (self.P6p*self.P6p - 1 + self.P1 < 0)
    def ctL2phi3(self):      
        return (self.P5p*self.P5p - 1 - self.P1 < 0)
    
    def numerical(self):
    
        a0 = 1 - self.P1*self.P1 - self.P6p*self.P6p*(1+self.P1) - self.P8p*self.P8p*(1-self.P1) - 4*self.P2*self.P2 - 4*self.P2*self.P6p*self.P8p; 
        a4 = 1 - self.P1*self.P1 - self.P4p*self.P4p*(1+self.P1) - self.P5p*self.P5p*(1-self.P1) - 4*self.P2*self.P2 + 4*self.P2*self.P4p*self.P5p; 
        a1 = 4*self.P3*self.P8p*self.P8p - 4*self.P3*self.P6p*self.P6p - 8*self.P1*self.P3 + 2*self.P5p*self.P6p*(1+self.P1) - 2*self.P4p*self.P8p*(1-self.P1) - 4*self.P2*self.P4p*self.P6p + 4*self.P2*self.P5p*self.P8p;
        a3 = 4*self.P3*self.P4p*self.P4p - 4*self.P3*self.P5p*self.P5p + 8*self.P1*self.P3 + 2*self.P5p*self.P6p*(1-self.P1) - 2*self.P4p*self.P8p*(1+self.P1) - 4*self.P2*self.P4p*self.P6p + 4*self.P2*self.P5p*self.P8p;
        a2 = 2 + 2*self.P1*self.P1 - 8*self.P2*self.P2 - 16*self.P3*self.P3 - (self.P4p*self.P4p+self.P6p*self.P6p)*(1-self.P1) - (self.P5p*self.P5p+self.P8p*self.P8p)*(1+self.P1) + 4*self.P2*self.P4p*self.P5p - 4*self.P2*self.P6p*self.P8p + 8*self.P3*self.P4p*self.P8p + 8*self.P3*self.P5p*self.P6p;
        b0 = self.P8p*self.P8p - 1 - self.P1 + 2*self.P2 + self.P6p*self.P8p; 
        b2 = self.P4p*self.P4p - 1 + self.P1 + 2*self.P2 - self.P4p*self.P5p; 
        b1 = self.P4p*self.P8p - 2*self.P3 + 0.5 * ( self.P4p*self.P6p - self.P5p*self.P8p );
        c0 = self.P8p*self.P8p - 1 - self.P1 - 2*self.P2 - self.P6p*self.P8p; 
        c2 = self.P4p*self.P4p - 1 + self.P1 - 2*self.P2 + self.P4p*self.P5p; 
        c1 = self.P4p*self.P8p - 2*self.P3 - 0.5 * ( self.P4p*self.P6p - self.P5p*self.P8p );
        
        nSteps = 100
        y = 0
        for step in range(nSteps+1):
            phi = 3.14159 * step / nSteps;
            sin2   = math.sin(phi)*math.sin(phi);
            sincos = math.sin(phi)*math.cos(phi);
            cos2   = math.cos(phi)*math.cos(phi);
            
            ctL5p = b0*sin2 + b1*sincos + b2*cos2;
            if ( ctL5p >= 0 ): continue
            ctL5m = c0*sin2 + c1*sincos + c2*cos2;
            if ( ctL5m >= 0 ): continue
            ctL1 = a0*sin2*sin2 + a1*sin2*sincos + a2*sin2*cos2 + a3*sincos*cos2 + a4*cos2*cos2;
            if ( ctL1 >= 0 ): continue
         
            tmp_y = min(ctL5p, ctL5m, ctL1)
            if (tmp_y < y):  y = tmp_y
        return (y >=0)



par_list = {
  "Fl" : [0, 1] ,
  "P1" : [-1, 1],
  "P2" : [-0.5, 0.5],
  "P3" : [-0.5, 0.5],
  "P4p": [-1*sqrt(2),sqrt(2)],
  "P5p": [-1*sqrt(2),sqrt(2)],
  "P6p": [-1*sqrt(2),sqrt(2)],
  "P8p": [-1*sqrt(2),sqrt(2)],
}
par_list = OrderedDict(sorted(par_list.items(), key=lambda t: t[0]))

bin_range = [0,1,2,3,5,7]
if args.bin != -1:
    bin_range = [int(args.bin)]  

for ibin in bin_range:
  fname = 'simFitResultsExtConstr/simFitResult_recoMC_fullAngular%s_dataStat_b%s.root'%(args.year,ibin)
  try:
    f = ROOT.TFile(fname,'open')
  except:
    print ('file %s not found'%fname)
    exit(0) 

  keys = f.GetListOfKeys()

  hist_list = []
  hist_list_phys_reg = []
  for i in par_list.keys():
    hist_list         .append( ROOT.TH1F('h_%s'%i     , '%s; %s; fitted value'%(i,i)     , 40, par_list[i][0], par_list[i][1] ))
    hist_list_phys_reg.append( ROOT.TH1F('h_%s_good'%i, '%s_good; %s; fitted value'%(i,i), 40, par_list[i][0], par_list[i][1] ))
  
  cond1 = ROOT.TH1F('cond1', 'ctL4phi1' , 2, 0, 2)
  cond2 = ROOT.TH1F('cond2', 'ctK2'     , 2, 0, 2)
  cond3 = ROOT.TH1F('cond3', 'ctL2phi1' , 2, 0, 2)
  cond4 = ROOT.TH1F('cond4', 'ctL2phi2' , 2, 0, 2)
  cond5 = ROOT.TH1F('cond5', 'ctL2phi3' , 2, 0, 2)
  cond6 = ROOT.TH1F('cond6', 'phi cond.', 2, 0, 2)
  
  lls = {}
  pars = fit_pars()
  for i in keys:
      if 'simFit' in i.GetName(): 
          fitres = f.Get(i.GetName()) 
          if fitres.status() == 0 and fitres.covQual()==3:
              lls[i.GetName()] = fitres.minNll()
  
              for ipar in par_list.keys():
                  setattr(pars, ipar, fitres.floatParsFinal().find( ipar).getVal() )
  
              for ih in hist_list:
                ih.Fill(  getattr(pars,  ih.GetName().strip("h_") ) )  
              for ih in hist_list_phys_reg:
                if (pars.ctL4phi1()*pars.ctK2()*pars.ctL2phi1()*pars.ctL2phi2()*pars.ctL2phi3()*pars.numerical()>0):
                  ih.Fill(  getattr(pars,  ih.GetName().strip("h_").strip("_good") ) )  
  
          cond1.Fill(pars.ctL4phi1()  )
          cond2.Fill(pars.ctK2()  )
          cond3.Fill(pars.ctL2phi1() )
          cond4.Fill(pars.ctL2phi2() )
          cond5.Fill(pars.ctL2phi3() )
          cond6.Fill(pars.numerical() )
          
  canv = ROOT.TCanvas('canv', 'canv')
  canv.Divide(3,3)
  for i, ih in enumerate(hist_list):
    canv.cd(i+1)
    ih.Draw()
    hist_list_phys_reg[i].SetLineColor(ROOT.kGreen)
    hist_list_phys_reg[i].Draw('same')
    
    leg = ROOT.TLegend(0.6,0.7,0.8,0.9, "", "brNDC")
    leg.AddEntry(ih                   , 'all converging fits', 'l')
    leg.AddEntry(hist_list_phys_reg[i], 'inside boundary', 'l')
    leg.Draw()
    
  canv.SaveAs('compare_datalike_results_b%s.pdf'%ibin)
  
  canv2 = ROOT.TCanvas('canv2', 'canv2', 700,400)
  canv2.Divide(3,2)
  canv2.cd(1)
  cond1.Draw()
  canv2.cd(2)
  cond2.Draw()
  canv2.cd(3)
  cond3.Draw()
  canv2.cd(4)
  cond4.Draw()
  canv2.cd(5)
  cond5.Draw()
  canv2.cd(6)
  cond6.Draw()
  
  canv2.SaveAs('compare_datalike_results_conditions_b%s.pdf'%ibin)

# minll = min(lls.values())
# for k,v in lls.items():
#     if minll==v:
#         f.Get(k).Print()
#         
        
        
 

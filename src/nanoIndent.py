##
# @file
# @brief Classes to evaluate indentation data and indenter tip
#
# Indentation data: indentation experiment
# - Methods: iso (can have multiple unloading segments), csm
# - Vendor: Agilent, Hysitron
#
# Indenter tip: shape of indenter tip and gantry stiffness (that what you calibrate)
#
# UNITS: one should use mSB units in this code, since Agilent area function is unit-dependent
# mSB: [mN], [um], [GPa] (force, length, stress)
#
##### Variables: differentiate different length ########
# array of full length: force, time, depth, validMask, ...  [used for plotting]
# array of valid length: E,H,Ac,hc, ... [only has the length where these values are valid]
#    force[validMask] = pMax
#    all these are vectors: OliverPharr et al methods are only vector functions
#
# Coding rules:
# - Change all variables: do not keep original-depth as can be reread and makes code less readable

# 2. clean frame stiffness, additional compilance  : differentiate between both
#    - fitting unloading curve: assume as intial guess m=1.5
import math, io, lmfit, re, os, traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from enum import Enum
from scipy.optimize import differential_evolution, fmin_tnc, fmin_l_bfgs_b, curve_fit, OptimizeResult, newton
from scipy import ndimage, interpolate
from scipy.signal import savgol_filter
from zipfile import ZipFile


# enum classes: make code more readable
class Method(Enum):
  ISO = 1    #one unloading
  MULTI = 2  #multiple ISO unloadings
  CSM = 3    #CSM method: number of unloadings ~ all data-points

class Vendor(Enum):
  Agilent = 1
  HysiHLD = 2
  HysiTXT = 3
  Micromaterials = 4
  FischerScope = 5

class FileType(Enum):
  Single = 1  #single test in file
  Multi  = 2  #multiple tests in file

class Indentation:
  def __init__(self, fileName, nuMat= 0.3, tip=None, verbose=2):
    """
    Initialize indentation experiment data

    Args:
       fileName: fileName to open (.xls, .hld)
       nuMat: material's Poisson ratio.
       tip:  tip class to use; None=perfect
       verbose: the higher, the more information printed: 2=default, 1=minimal, 0=print nothing
    """
    self.nuMat = nuMat
    self.nuIndent = 0.07
    self.EIndent  = 1140                                    #GPa from Oliver,Pharr Method paper
    self.beta = 0.75                                        #beta: contact depth coefficient
    self.verbose = verbose
    self.method    = Method.ISO                             #iso default: csm uses different methods
    self.onlyLoadingSegment = False                         #use all data by default
    if tip is None: tip = Tip()
    self.tip = tip
    self.iLHU   = [ [-1,-1,-1,-1] ]                         #indicies of Load-Hold-Unload cycles (StartLoad-StartHold-StartUnload-EndLoad)
    self.iDrift = [-1,-1]                                   #start and end indicies of drift segment
    self.meta = {}                                          #some results come from input file, others are added by analysis
    self.slope= []

    #initialize and load first data set
    #set default parameters
    success = False
    if not os.path.exists(fileName) and fileName!='':
      print("*ERROR* __init__: file does not exist",fileName)
      return None
    if fileName.endswith(".xls") or fileName.endswith(".xlsx"):
      # KLA, Agilent, Keysight, MTS
      self.vendor = Vendor.Agilent
      self.fileType = FileType.Multi
      self.unloadPMax = 0.999
      self.unloadPMin = 0.5
      success = self.loadAgilent(fileName)
    if (fileName.endswith(".hld") or fileName.endswith(".txt")) and not success:
      # Hysitron
      if fileName.endswith(".hld"): self.vendor = Vendor.HysiHLD
      else                        : self.vendor = Vendor.HysiTXT
      self.fileType = FileType.Single
      self.unloadPMax = 0.95
      self.unloadPMin = 0.4
      success = self.loadHysitron(fileName)
    if fileName.endswith(".txt") and not success:
      # Micromaterials
      self.vendor = Vendor.Micromaterials
      self.fileType = FileType.Single
      self.unloadPMax = 0.99
      self.unloadPMin = 0.5
      success = self.loadMicromaterials(fileName)
    if fileName.endswith(".zip") and not success:
      # Micromaterials
      self.vendor = Vendor.Micromaterials
      self.fileType = FileType.Multi
      self.unloadPMax = 0.99
      self.unloadPMin = 0.5
      success = self.loadMicromaterialsZip(fileName)
    if fileName.endswith(".txt") and not success:
      # Fischer Scope
      self.vendor = Vendor.FischerScope
      self.fileType = FileType.Multi
      self.unloadPMax = 0.95
      self.unloadPMin = 0.21
      success = self.loadFischerScope(fileName)
    return

  ##
  # @name CONVENTIONAL NANOINDENTATION FUNCTIONS: area, E,.
  # only area prefactors and stiffness are class variables
  #@{
  def YoungsModulus(self, redE, nuThis=-1):
    """
    Calculate the Youngs modulus from the reduced Youngs modulus

    Args:
       redE: reduced Youngs modulus [GPa]
       nuThis: use a non-standard Young's modulus
    """
    nu = self.nuMat
    if nuThis>0:
      nu = nuThis
    E = (1.0-nu*nu) / ( 1.0/redE - (1.0-self.nuIndent*self.nuIndent)/self.EIndent )
    return E


  def ReducedModulus(self, E, nuThis=-1):
    """
    Calculate the reduced modulus from the Youngs modulus

    Args:
       E: Youngs modulus [GPa]
       nuThis: use a non-standard Young's modulus
    """
    nu = self.nuMat
    if nuThis>0:
      nu = nuThis
    Ered =  1.0/(  (1.0-nu*nu)/E + (1.0-self.nuIndent*self.nuIndent)/self.EIndent )
    return Ered


  def OliverPharrMethod(self, S, P, h):
    """
    Conventional Oliver-Pharr indentation method to calculate reduced Modulus E*

    The following equations are used in that order:<br>
      h_c = h-beta P/S<br>
      A = h_c(prefactors)<br>
      S = 2/sqrt(pi) sqrt(A) E*<br>
      A the contact area, h_c the contact depth<br>

    Args:
       S: stiffness = slope dP/dh
       P: maximal force
       h: total penetration depth

    Returns:
       as list: reducedModulus E*, A, h_c
    """
    threshA = 1.e-12  #units in um: threshold = 1pm^2
    h_c = h - self.beta*P/S
    A   = self.tip.areaFunction(h_c)
    A[A< threshA] = threshA  # prevent zero or negative area that might lock sqrt
    E   = S / (2.0*np.sqrt(A)/np.sqrt(np.pi))
    return [E, A, h_c]


  def inverseOliverPharrMethod(self, S, P, E):
    """
    Inverse Oliver-Pharr indentation method to calculate contact area A

    equations and variable definitions given above; order in reverse order

    Args:
       S: slope dP/dh
       P: maximal force
       E: reducedModulus E*

    Returns:
       h penetration depth
    """
    A = math.pow( S / (2.0*E/math.sqrt(math.pi))  ,2)
    h_cGuess = math.sqrt(A / 24.494) # first guess: perfect Berkovich
    # print("  DEBUG A,self.beta,P,S,h_c0", A, self.beta, P, S, h_cGuess)
    h_c = self.tip.areaFunctionInverse(A, h_c0=h_cGuess)
    h = h_c + self.beta*P/S
    return h.flatten()


  def UnloadingPowerFunc(self,h_,B,hf,m):
    """
    internal function describing the unloading regime
    """
    import warnings
    warnings.filterwarnings('error') #catch runtime warning and raise exception which is catched in the caller-function
    value = B*np.power(h_-hf,m)
    return value


  def stiffnessFromUnloading(self, P, h, plot=False):
    """
    Calculate single unloading stiffness from Unloading
    see G200 manual, p7-6

    Args:
       P: vector of forces
       h: vector of depth
       plot: plot results

    Results:
       stiffness, validMask [values of P,h where stiffness is determined], mask, optimalVariables, powerlawFit-success
    """
    debug = False
    if self.method== Method.CSM:
      print("*ERROR* Should not land here: CSM method")
      return None,None,None,None
    if self.verbose>1: print("Number of unloading segments:"+str(len(self.iLHU))+"  Method:"+str(self.method))
    S, mask, opt, powerlawFit = [], None, None, []
    maxDeltaP = -0.01
    t = self.t - np.min(self.t)  #undo resetting zero during init
    validMask = np.zeros_like(P, dtype=bool)
    if plot:
      plt.plot(h,P,'-k')
    for cycleNum, cycle in enumerate(self.iLHU):
      loadStart, loadEnd, unloadStart, unloadEnd = cycle
      if loadStart>loadEnd or loadEnd>unloadStart or unloadStart>unloadEnd:
        print('*ERROR* stiffnessFromUnloading: indicies not in order:',cycle)
      maskSegment = np.zeros_like(h, dtype=bool)
      maskSegment[unloadStart:unloadEnd+1] = True
      maskForce   = np.logical_and(P<P[loadEnd]*self.unloadPMax, P>P[loadEnd]*self.unloadPMin)
      mask        = np.logical_and(maskSegment,maskForce)
      if plot:
        plt.plot(h[mask],P[mask],'ob')
      B0  = (P[mask][-1]-P[mask][0])/(h[mask][-1]-h[mask][0])
      hf0 = h[mask][0] - P[mask][0]/B0
      try:
        opt, _ = curve_fit(self.UnloadingPowerFunc, h[mask],P[mask], p0=[B0,hf0,1.5], ftol=1e-4, maxfev=1000 )#set ftol to 1e-4 if accept more and fail less
        B,hf,m = opt
        if np.isnan(B):
          raise ValueError("NAN after fitting")
        powerlawFit.append(True)
      except:
        if self.verbose>0:
          print("stiffnessFrommasking: #",cycleNum," Fitting failed. use linear")
        B  = (P[mask][-1]-P[mask][0])/(h[mask][-1]-h[mask][0])
        hf = h[mask][0] -P[mask][0]/B
        m  = 1.
        opt= (B,hf,m)
        powerlawFit.append(False)
      S_ = B*m*math.pow( (h[unloadStart]-hf), m-1)
      S.append(S_)
      validMask[unloadStart]=True
      if plot:
        plt.plot(h[mask],   self.UnloadingPowerFunc(h[mask],B,hf,m),'m-')
        Sn= P[unloadStart]-S_*h[unloadStart]
        plt.plot(h[mask],   S_*h[mask]+Sn, 'r--', lw=3)
    if plot:
      plt.xlim(left=0)
      plt.ylim(bottom=0)
      plt.title('magenta: power function, red: linear slope')
      plt.xlabel(r'depth [$\mu m$]')
      plt.ylabel(r'force [$mN$]')
      plt.show()
    return S,validMask, mask, opt, powerlawFit


  def popIn(self, correctH=True, plot=True, removeInitialNM=2.):
    """
    Search for pop-in by jump in depth rate

    Certainty:
    - deltaSlope: higher is better (difference in elastic - plastic slope). Great indicator
    - prefactor: higher is better (prefactor of elastic curve). Great indicator
    - secondRate: lower is better (height of second largest jump). Nice indicator 0.25*deltaRate
    - covElast: lower is better. bad indicator
    - deltaH: higher is better (delta depth in jump). bad indicator
    - deltaRate: higher is better (depth rate during jump). bad indicator

    Future: iterate over largest, to identify best

    Args:
       correctH: correct depth such that curves aligned
       plot: plot pop-in curve
       removeInitialNM: remove initial nm from data as they have large scatter

    Returns:
       pop-in force, dictionary of certainty
    """
    maxPlasticFit = 150
    minElasticFit = 0.01

    mask = (self.h[self.valid]-np.min(self.h[self.valid]))  >removeInitialNM/1.e3
    h = self.h[self.valid][mask]
    p = self.p[self.valid][mask]

    depthRate = (h[1:]-h[:-1])
    x_        = np.arange(len(depthRate))
    fits      = np.polyfit(x_,depthRate,2)  #substract 2nd order fit b/c depthRate increases over time
    depthRate-= np.polyval(fits,x_)
    iJump     = np.argmax(depthRate)
    iMax      = min(np.argmax(p), iJump+maxPlasticFit)      #max for fit: 150 data-points or max. of curve
    iMin      = np.min(np.where(p>minElasticFit))
    fitPlast  = np.polyfit(h[iJump+1:iMax],p[iJump+1:iMax],2) #does not have to be parabola, just close fit
    slopePlast= np.polyder(np.poly1d(fitPlast))(h[iJump+1] )
    def funct(depth, prefactor, h0):
      diff           = depth-h0
      if type(diff)==np.float64:
        diff = max(diff,0.0)
      else:
        diff[diff<0.0] = 0.0
      return prefactor* (diff)**(3./2.)
    fitElast, pcov = curve_fit(funct, h[iMin:iJump], p[iMin:iJump], p0=[100.,0.])
    slopeElast= (funct(h[iJump],*fitElast) - funct(h[iJump]*0.9,*fitElast)) / (h[iJump]*0.1)
    fPopIn    = p[iJump]
    certainty = {"deltaRate":depthRate[iJump], "prefactor":fitElast[0], "h0":fitElast[1], \
                 "deltaSlope": slopeElast-slopePlast, 'deltaH':h[iJump+1]-h[iJump],\
                 "covElast":pcov[0,0] }
    listDepthRate = depthRate.tolist()
    iJump2 = np.argmax(listDepthRate)
    while ((iJump2-iJump)<3):
      del listDepthRate[iJump2]
      iJump2 = np.argmax(listDepthRate)
    certainty["secondRate"] = np.max(listDepthRate)
    if plot:
      _, ax1 = plt.subplots()
      ax2 = ax1.twinx()
      ax1.plot(self.h,self.p)
      h_ = np.linspace(self.h[iJump+1],self.h[iMax])
      ax1.plot(h_, np.polyval(fitPlast,h_))
      ax1.plot(self.h[iMin:iJump], funct(self.h[iMin:iJump],*fitElast))
      ax2.plot(h[:-1],depthRate,'r')
      ax1.axvline(h[iJump], color='k', linestyle='dashed')
      ax1.axhline(fPopIn, color='k', linestyle='dashed')
      ax1.set_xlim(right=4.*self.h[iJump])
      ax1.set_ylim(top=4.*self.p[iJump], bottom=0)
      plt.show()
    if correctH:
      self.h -= certainty["h0"]
    return fPopIn, certainty


  #@}
  ##
  # @name Calculate YoungsModulus, Hardess and deterimine area function
  # Access to class variables
  #@{
  def calcYoungsModulus(self, minDepth=-1, plot=False):
    """
    Calculate and plot Young's modulus as a function of the depth
    -  use corrected h, s (do not recalculate)

    Args:
       minDepth: minimum depth for fitting horizontal; if negative: no line is fitted
       plot: plot comparison this calculation to data read from file

    Returns:
       average Youngs, minDepth>0
    """
    self.modulusRed, self.A_c, self.h_c = self.OliverPharrMethod(self.slope, self.p[self.valid], self.h[self.valid])
    E = self.YoungsModulus(self.modulusRed)
    if minDepth>0:
      #eAve = np.average(       self.modulusRed[ self.h>minDepth ] )
      eAve = np.average( E[  np.bitwise_and(E>0, self.h[self.valid]>minDepth) ] )
      eStd = np.std(     E[  np.bitwise_and(E>0, self.h[self.valid]>minDepth) ] )
      print("Average and StandardDeviation of Young's Modulus",round(eAve,1) ,round(eStd,1) ,' [GPa]')
    else:
      eAve, eStd = -1, 0
    if plot:
      if not self.modulus is None:
        plt.plot(self.h[self.h>minDepth], self.modulus[self.h>minDepth], '-r', lw=3, label='read')
      plt.plot(self.h[self.h>minDepth], E[self.h>minDepth], '-b', label='calc')
      if minDepth>0:
        plt.axhline(eAve, color='k')
        plt.axhline(eAve+eStd, color='k', linestyle='dashed')
        plt.axhline(eAve-eStd, color='k', linestyle='dashed')
        plt.ylim([eAve-4*eStd,eAve+4*eStd])
      plt.xlabel(r'depth [$\mu m$]')
      plt.ylim(ymin=0)
      plt.ylabel('Youngs modulus [GPa]')
      plt.legend(loc=0)
      plt.show()
    self.modulus = E
    return eAve


  def calcHardness(self, minDepth=-1, plot=False):
    """
    Calculate and plot Hardness as a function of the depth

    Args:
       minDepth: minimum depth for fitting horizontal; if negative: no line is fitted
       plot: plot comparison this calculation to data read from file
    """
    H = self.p[self.valid] / self.OliverPharrMethod(self.slope, self.p[self.valid], self.h[self.valid])[1] #use area function
    if plot:
      plt.plot(self.h, H , '-b', label='calc')
      if not self.hardness is None:
        plt.plot(self.h, self.hardness, '-r', label='readFromFile')
      if minDepth>0:
        HAve = np.average( H[  np.bitwise_and(H>0, self.h>minDepth) ] )
        HStd = np.std(     H[  np.bitwise_and(H>0, self.h>minDepth) ] )
        print("Average and StandardDeviation of Hardness",round(HAve,1) ,round(HStd,1) ,' [GPa]')
        plt.axhline(HAve, color='b')
        plt.axhline(HAve+HStd, color='b', linestyle='dashed')
        plt.axhline(HAve-HStd, color='b', linestyle='dashed')
      plt.xlabel(r'depth [$\mu m$]')
      plt.ylabel(r'hardness [$GPa$]')
      plt.legend(loc=0)
      plt.show()
    self.hardness = H
    return


  def calcStiffness2Force(self, minDepth=0.01, plot=True, calibrate=False):
    """
    Calculate and plot stiffness squared over force as a function of the depth

    Args:
       minDepth: minimum depth for fitting line
       plot: plot curve and slope
       calibrate: calibrate additional stiffness and save value
    """
    compliance0 = self.compliance
    prefactors = None
    def errorFunction(compliance):
      s   = 1./(1./self.sRaw-compliance)
      s2f = np.divide(np.multiply(s,s),self.p)
      h   = self.hRaw-compliance*self.p
      h_ = h[ h>minDepth ]
      s2f_  = s2f[ h>minDepth ]
      if len(h_)>4:
        prefactors = np.polyfit(h_,s2f_,1)
        print(compliance,"Fit f(x)=",prefactors[0],"*x+",prefactors[1])
        return np.abs(prefactors[0])
      else:
        print("*WARNING*: too short vector",len(h_))
        return 9999999.
    if calibrate:
      result = fmin_l_bfgs_b(errorFunction, compliance0, bounds=[(-0.1,0.1)], approx_grad=True, epsilon=0.000001, factr=1e11)
      print("  Best values   ",result[0], "\tOptimum residual:",np.round(result[1],3))
      print('  Number of function evaluations~size of globalData',result[2]['funcalls'])
      self.compliance = result[0]
      compliance0 = self.compliance
      #self.correct_H_S()
    if plot:
      s = 1./(1./self.sRaw-self.compliance)
      s2f = np.divide(np.multiply(s,s),self.p)
      h   = self.hRaw-compliance0*self.p
      h_ = h[ h>minDepth ]
      s2f_  = s2f[ h>minDepth ]
      prefactors = np.polyfit(h_,s2f_,1)
      plt.plot(h,s2f, 'b-')
      s2fFit = np.polyval(prefactors,h)
      plt.plot(h, s2fFit, 'r-', lw=3)
      plt.xlabel(r'depth [$\mu m$]')
      plt.ylabel(r'stiffness2/force [$GPa$]')
      plt.show()
    return prefactors


  def tareDepthForce(self, slopeThreshold=100, compareRead=False, plot=False):
    """
    Calculate surface contact (by slope being larger than threshold)
    and offset depth,force,time by the surface

    Future improvements:
    - surface identification in future
    - handle more cases

    Args:
       slopeThreshold: threshold slope in P-h curve for contact: 200,300
       compareRead: compare new results to the ones from the file
       plot: plot comparison new data and data from file
    """
    if self.vendor!=Vendor.Agilent or self.method==Method.CSM:
      print("tareDepthForce only valid for ISO method of Agilent at this moment")
      return
    iSurface = np.min(np.where(self.pVsHSlope>slopeThreshold))#determine point of contact
    h = self.hRaw   - self.hRaw[iSurface]                   #tare to point of contact
    p = self.pRaw   - self.pRaw[iSurface]
    t = self.tTotal - self.tTotal[iSurface]
    h-= p/self.frameStiffness                               #compensate depth for instrument deflection
    maskDrift = np.zeros_like(h, dtype=bool)
    maskDrift[self.iDrift[0]:self.iDrift[1]]   =  True
    tMiddle = (t[self.iDrift[1]]+t[self.iDrift[0]])/2
    maskDrift = np.logical_and(maskDrift, t>=tMiddle)
    iDriftS, iDriftE = np.where(maskDrift)[0][0],np.where(maskDrift)[0][-1]
    driftRate        = (h[iDriftE]-h[iDriftS])/(t[iDriftE]-t[iDriftS])  #calc. as rate between last and first point
                                #according to plot shown in J.Hay Univerisity part 3; fitting line would be different
    print("Drift rate: %.3f nm/s"%(driftRate*1e3))
    h-= driftRate*t                                          #compensate thermal drift
    p-= self.slopeSupport*(self.hRaw-self.hRaw[iSurface])    #compensate supporting mechanism (use original data since h changed)
    if compareRead:
      mask = self.h>0.010                                    #10nm
      error = (h[mask]-self.h[mask])/self.h[mask]
      print("Error in h: {0:.2f}%".format(np.linalg.norm(error)/len(error)*100.) )
      error = (p[mask]-self.p[mask])/self.p[mask]
      print("Error in p: {0:.2f}%".format(np.linalg.norm(error)/len(error)*100.) )
      error = (t[mask]-self.t[mask])/self.t[mask]
      print("Error in t: {0:.2f}%".format(np.linalg.norm(error)/len(error)*100.) )
    if plot:
      fig, ax1 = plt.subplots()
      ax2 = ax1.twinx()
      ax1.plot(t,p, label='new')
      ax1.plot(self.t,self.p, label='read')
      ax1.axhline(0, linestyle='dashed', c='k')
      ax1.axvline(0, linestyle='dashed', c='k')
      ax1.legend(loc=2)
      ax2.plot(self.t, self.pVsHSlope, "C2--", label='pVsHSlope')
      ax1.set_xlabel(r"depth [$\mu m$]")
      ax1.set_ylabel(r"force [$mN$]")
      ax2.set_ylabel(r"depth [$mN/\mu m$]", color='C2')
      plt.show()
    #set newly obtained data
    self.h, self.p, self.t = h, p, t
    return


  def analyse(self):
    """
    update slopes/stiffness, Young's modulus and hardness after displacement correction by:
    - drift change
    - compliance change

    ONLY DO ONCE AFTER LOADING FILE: if this causes issues introduce flag analysed which is toggled during loading and analysing
    """
    self.h -= self.tip.compliance*self.p
    if self.method == Method.CSM:
      self.slope = 1./(1./self.slope-self.tip.compliance)
    else:
      self.slope, self.valid, _, _ , _= self.stiffnessFromUnloading(self.p, self.h)
      self.slope = np.array(self.slope)
    self.k2p = self.slope*self.slope/self.p[self.valid]
    self.calcYoungsModulus()
    self.calcHardness()
    return


  def analysePreDrift(self, plot=True):
    """
    Analyse drift segment by Hysitron before the test

    Args:
       plot: plot drift data

    Results:
       drift in um/s
    """
    if not self.vendor == Vendor.HysiHLD: return None
    t = self.dataDrift[:,0]
    h = self.dataDrift[:,1]
    rate = np.zeros_like(t)
    rate[:] = np.nan
    for idxEnd,ti in enumerate(t):
      if ti<20: continue
      idxStart = np.argmin(np.abs( t[:]-(ti-20.) )  )
      rate[idxEnd] = (h[idxEnd]-h[idxStart])/(t[idxEnd]-t[idxStart])
    idxEnd = np.argmin(np.abs( t[:]-(40.) )  )
    drift = rate[idxEnd]
    print("Drift:", round(drift*1.e3,3),"nm/s")
    self.meta['drift'] = drift
    if plot:
      _, ax1 = plt.subplots()
      ax2 = ax1.twinx()
      ax1.plot(t,rate*1.e3,'r')
      ax1.axhline(0.05,c='r',linestyle='dashed')
      ax2.plot(t,h*1.e3,'b')
      ax1.set_xlabel('time [$s$]')
      ax1.set_ylabel('drift rate [$nm/s$]', color='r')
      ax2.set_ylabel('depth [$nm$]', color='b')
      plt.show()
    return drift

  def identifyLoadHoldUnload(self,plot=False):
    """
    internal method: identify ALL load - hold - unload segments in data

    Args:
       plot: verify by plotting
    """
    if self.method==Method.CSM:
      self.identifyLoadHoldUnloadCSM()
      return
    rate = self.p[1:]-self.p[:-1]
    #using histogram, define masks for loading and unloading
    hist,bins= np.histogram(rate , bins=1000) #TODO Better algorithm: 1000 is good for FischerScope, but leads to other failures
    if self.vendor==Vendor.HysiTXT:
      hist = ndimage.filters.gaussian_filter1d(hist,2)
    binCenter = (bins[1:]+bins[:-1])/2
    peaks = np.where(hist>10)[0]                  #peaks with more than 10 items
    try:
      zeroID = np.argmin(np.abs(binCenter[peaks]))  #id which is closest to zero
      zeroValue = binCenter[peaks][zeroID]
    except:
      self.iLHU = []
      return False
    ## Better algorithm: look for closest zero historgram-peak to zeroValue; take that to calculate delta
    zeroPeaks = np.logical_and(hist<0.3, binCenter<zeroValue)
    zeroPeaks = np.where(zeroPeaks)[0]
    firstZero = np.argmin(np.abs(zeroValue-binCenter[zeroPeaks]))
    zeroDelta = abs(binCenter[zeroPeaks][firstZero]-zeroValue)
    ## Old algorithm: find next binCenter and calculate its distance; depends on histogram step size; not good
    # zeroDelta = max( binCenter[zeroID+1]-binCenter[zeroID],\
    #                  binCenter[zeroID]-binCenter[zeroID-1])/2
    # zeroDelta = max( zeroDelta, 0.002)    #0.002mN/s seems to be accuracy limit
    loadMask  = rate>(zeroValue+zeroDelta)
    unloadMask= rate<(zeroValue-zeroDelta)
    if plot:     # verify visually
      plt.plot(binCenter,hist,'o')#, width=0.001)
      plt.axvline(zeroValue, c='k')
      plt.axvline(binCenter[zeroPeaks][firstZero], c='r')
      plt.axvline(zeroValue+zeroDelta, c='k', linestyle='dashed')
      plt.axvline(zeroValue-zeroDelta, c='k', linestyle='dashed')
      plt.ylabel('count []')
      plt.xlabel('rate [$mN/sec$]')
      plt.show()
      plt.plot(self.t[1:],rate)
      plt.axhline(zeroValue, c='k')
      plt.axhline(zeroValue+zeroDelta, c='k', linestyle='dashed')
      plt.axhline(zeroValue-zeroDelta, c='k', linestyle='dashed')
      plt.xlabel('time incr. []')
      plt.ylabel('rate [$mN/sec$]')
      plt.show()
    #clean small fluctuations
    size = 7
    loadMask = ndimage.binary_closing(loadMask, structure=np.ones((size,)) )
    unloadMask = ndimage.binary_closing(unloadMask, structure=np.ones((size,)))
    loadMask = ndimage.binary_opening(loadMask, structure=np.ones((size,)))
    unloadMask = ndimage.binary_opening(unloadMask, structure=np.ones((size,)))
    #find index where masks are changing from true-false
    loadMask  = np.r_[False,loadMask,False] #pad with false on both sides
    unloadMask= np.r_[False,unloadMask,False]
    loadIdx   = np.flatnonzero(loadMask[1:]   != loadMask[:-1])
    unloadIdx = np.flatnonzero(unloadMask[1:] != unloadMask[:-1])
    if plot:     # verify visually
      plt.plot(self.p)
      plt.plot(loadIdx[::2],  self.p[loadIdx[::2]],  'o',label='load',markersize=12)
      plt.plot(loadIdx[1::2], self.p[loadIdx[1::2]], 'o',label='hold',markersize=10)
      plt.plot(unloadIdx[::2],self.p[unloadIdx[::2]],'o',label='unload',markersize=8)
      plt.plot(unloadIdx[1::2],self.p[unloadIdx[1::2]],'o',label='unload-end',markersize=6)
      plt.legend(loc=0)
      plt.title("BEFORE Cleaning")
      plt.xlabel('time incr. []')
      plt.ylabel('force [$mN$]')
      plt.show()
    while len(loadIdx)<len(unloadIdx) and unloadIdx[0]<loadIdx[0]:
      print("*WARNING* identifyLoadHoldUnload: cut two from front of unloadIdx: UNDESIRED")
      unloadIdx = unloadIdx[2:]
    while len(loadIdx)<len(unloadIdx) and unloadIdx[-3]>loadIdx[-1]:
      print("*WARNING* identifyLoadHoldUnload: cut two from end of unloadIdx: UNDESIRED")
      unloadIdx = unloadIdx[:-2]
    while len(loadIdx)>len(unloadIdx) and loadIdx[3]<unloadIdx[1]:
      print("*WARNING* identifyLoadHoldUnload: cut two from front of loadIdx: UNDESIRED")
      loadIdx = loadIdx[2:]
    if plot:     # verify visually
      plt.plot(self.p)
      plt.plot(loadIdx[::2],  self.p[loadIdx[::2]],  'o',label='load',markersize=12)
      plt.plot(loadIdx[1::2], self.p[loadIdx[1::2]], 'o',label='hold',markersize=10)
      plt.plot(unloadIdx[::2],self.p[unloadIdx[::2]],'o',label='unload',markersize=8)
      plt.plot(unloadIdx[1::2],self.p[unloadIdx[1::2]],'o',label='unload-end',markersize=6)
      plt.legend(loc=0)
      plt.title("AFTER Cleaning")
      plt.xlabel('time incr. []')
      plt.ylabel('force [$mN$]')
      plt.show()
    if len(loadIdx)!=len(unloadIdx):
      print("*WARNING* identifyLoadHoldUnload: Repair required")
      #TODO Repair possibly that this is not required
      self.identifyLoadHoldUnloadCSM()
      return
    self.iLHU = []
    for i,_ in enumerate(loadIdx[::2]):
      self.iLHU.append([loadIdx[::2][i],loadIdx[1::2][i],unloadIdx[::2][i],unloadIdx[1::2][i]])
    if len(self.iLHU)>1:
      self.method=Method.MULTI
    if self.verbose>1:
      print("Number Unloading segments",len(self.iLHU))
    #drift segments
    iDriftS = unloadIdx[1::2][i]+1
    iDriftE = len(self.p)-1
    if iDriftS+1>iDriftE:
      iDriftS=iDriftE-1
    self.iDrift = [iDriftS,iDriftE]
    return True


  def identifyLoadHoldUnloadCSM(self):
    """
    internal method: identify load - hold - unload segment in CSM data

    Backup: if identifyLoadHoldUnload fails
    """
    iSurface = np.min(np.where( self.h>=0                     ))
    iLoad    = np.min(np.where( self.p-np.max(self.p)*0.999>0 ))
    if iLoad<len(self.p)-1:
      iHold    = np.max(np.where( self.p-np.max(self.p)*0.999>0 ))
      hist,bins= np.histogram( self.p[iHold:] , bins=1000)
      pDrift   = bins[np.argmax(hist)+1]
      pCloseToDrift = np.logical_and(self.p>pDrift*0.999,self.p<pDrift/0.999)
      pCloseToDrift[:iHold] = False
      if len(pCloseToDrift[pCloseToDrift])>3:
        iDriftS  = np.min(np.where( pCloseToDrift ))
        iDriftE  = np.max(np.where( pCloseToDrift ))
      else:
        iDriftS   = len(self.p)-2
        iDriftE   = len(self.p)-1
      if not (iSurface<iLoad and iLoad<iHold and iHold<iDriftS and iDriftS<iDriftE and iDriftE<len(self.h)):
        print("*ERROR* identifyLoadHoldUnloadCSM in identify load-hold-unloading cycles")
        print(iSurface,iLoad,iHold,iDriftS,iDriftE, len(self.h))
    else:  #This part is required
      if self.method != Method.CSM:
        print("*WARNING*: no hold or unloading segments in data")
      iHold     = len(self.p)-3
      iDriftS   = len(self.p)-2
      iDriftE   = len(self.p)-1
    self.iLHU   = [[iSurface,iLoad,iHold,iDriftS]]
    self.iDrift = [iDriftS,iDriftE]
    return



  #@}
  ##
  # @name FILE HANDLING; PLOTTING
  # Access to class variables
  #@{
  def nextTest(self):
    """
    Wrapper for all next test for all vendors
    """
    if self.vendor == Vendor.Agilent:
      success = self.nextAgilentTest()
    elif self.vendor == Vendor.Micromaterials:
      success = self.nextMicromaterialsTest()
    elif self.vendor == Vendor.FischerScope:
      success = self.nextFischerScopeTest()
    else:
      print("No multiple tests in file")
      success = False
    return success


  def loadAgilent(self, fileName):
    """
    Initialize G200 excel file for processing

    Args:
       fileName: file name
    """
    self.testList = []
    self.fileName = fileName    #one file can have multiple tests
    self.indicies = {}
    #TODO STEFFEN
    # wb = pd.read_excel(fileName,sheet_name='Required Inputs')
    # self.meta.update( dict(wb.iloc[-1]) )
    # if self.meta['Poissons Ratio']!=self.nuMat:
    #   print("*WARNING*: your Poisson Ratio is different than in file.",self.nuMat,self.meta['Poissons Ratio'])
    ## TODO check if CSM can be deduced form other sheets
    self.datafile = pd.read_excel(fileName, sheet_name=None)
    tagged = []
    code = {"Load On Sample":"p", "Force On Surface":"p", "LOAD":"p", "_Load":"pRaw", "Raw Load":"pRaw"\
                ,"Displacement Into Surface":"h", "DEPTH":"h", "_Displacement":"hRaw", "Raw Displacement": "hRaw"\
                ,"Time On Sample":"t", "Time in Contact":"t", "TIME":"t", "Time":"tTotal"\
                ,"Load vs Disp Slope":"pVsHSlope", "_Column": "Column", "_Frame": "Frame", "Stiffness":"slope" \
                ,"Contact Area":"A_c", "Contact Depth":"h_c"\
                ,"Harmonic Displacement":"hHarmonic", "Harmonic Load":"pHarmonic","Phase Angle":"phaseAngle"\
                ,"Harmonic Stiffness":"slopeInvalid", "Harmonic Contact Stiffness":"slope", "STIFFNESS":"slope", "Stiffness Squared Over Load":"k2p"\
                , "Support Spring Stiffness":"slopeSupport", "Frame Stiffness": "frameStiffness"\
                ,"Hardness":"hardness", "H_IT Channel":"hardness", "Modulus": "modulus", "E_IT Channel": "modulus","Reduced Modulus":"modulusRed"\
                ,"Scratch Distance": "s", "XNanoPosition": "x", "YNanoPosition": "y", "X Position": "xCoarse", "Y Position": "yCoarse"\
                ,"TotalLateralForce": "L", "X Force": "pX", "_XForce": "pX", "Y Force": "pY", "_YForce": "pY"\
                ,"_XDeflection": "Ux", "_YDeflection": "Uy" }
    self.fullData = ['h','p','t','pVsHSlope','hRaw','pRaw','tTotal','slopeSupport']
    if self.verbose>1: print("=============  "+fileName+"  ============")
    for dfName in self.datafile.keys():
      df    = self.datafile.get(dfName)
      if "Test " in dfName and not "Tagged" in dfName and not "Test Inputs" in dfName:
        self.testList.append(dfName)
        #print "  I should process sheet |",sheet.name,"|"
        if len(self.indicies)==0:               #find index of colums for load, etc
          for cell in df.columns:
            if cell in code:
              self.indicies[code[cell]] = cell
              if self.verbose>1: print("     %-30s : %-20s "%(cell,code[cell]) )
            else:
              if self.verbose>1: print(" *** %-30s NOT USED"%cell)
            if "Harmonic" in cell:
              self.method = Method.CSM
          #reset to ensure default values are set
          if not "p" in self.indicies: self.indicies['p']=self.indicies['pRaw']
          if not "h" in self.indicies: self.indicies['h']=self.indicies['hRaw']
          if not "t" in self.indicies: self.indicies['t']=self.indicies['tTotal']
          #if self.verbose: print("   Found column names: ",sorted(self.indicies))
      if "Tagged" in dfName: tagged.append(dfName)
    if len(tagged)>0 and self.verbose>1: print("Tagged ",tagged)
    if not ("t" in self.indicies) or not ("p" in self.indicies) or \
       not ("h" in self.indicies) or not ("slope" in self.indicies)  :
          print("*WARNING*: INDENTATION: Some index is missing (t,p,h,slope) should be there")
    self.meta['measurementType'] = 'MTS, Agilent Indentation XLS'
    self.nextAgilentTest()
    return True


  def nextAgilentTest(self, newTest=True):
    """
    Go to next sheet in worksheet and prepare indentation data

    Data: _Raw: without frame stiffness correction, _Frame:  with frame stiffness correction (remove postscript finally)
    - only affects/applies directly depth (h) and stiffness (s)
    - modulus, hardness and k2p always only use the one with frame correction

    Args:
      newTest: take next sheet (default)
    """
    if self.vendor!=Vendor.Agilent: return False #cannot be used
    if len(self.testList)==0: return False   #no sheet left
    if newTest:
      self.testName = self.testList.pop(0)

    #read data and identify valid data points
    df     = self.datafile.get(self.testName)
    h       = np.array(df[self.indicies['h'    ]][1:-1], dtype=np.float)
    self.validFull = np.isfinite(h)
    if 'slope' in self.indicies:
      slope   = np.array(df[self.indicies['slope']][1:-1], dtype=np.float)
      self.valid =  np.isfinite(slope)
      self.valid[self.valid] = slope[self.valid] > 0.0  #only valid points if stiffness is positiv
    else:
      self.valid = self.validFull
    for index in self.indicies:
      data = np.array(df[self.indicies[index]][1:-1], dtype=np.float)
      mask = np.isfinite(data)
      mask[mask] = data[mask]<1e99
      self.valid = np.logical_and(self.valid, mask)                       #adopt/reduce mask continously

    #Run through all items again and crop to only valid data
    for index in self.indicies:
      data = np.array(df[self.indicies[index]][1:-1], dtype=np.float)
      if not index in self.fullData:
        data = data[self.valid]
      else:
        data = data[self.validFull]
      setattr(self, index, data)
      # print(index, len(data))
    self.valid = self.valid[self.validFull]
    #  now all fields (incl. p) are full and defined

    success = self.identifyLoadHoldUnload()
    if self.onlyLoadingSegment and self.method==Method.CSM:
      # print("Length test",len(self.valid), len(self.h[self.valid]), len(self.p[self.valid])  )
      iMin, iMax = 2, self.iLHU[0][1]
      self.valid[iMax:] = False
      self.valid[:iMin] = False
      self.slope = self.slope[iMin:np.sum(self.valid)+iMin]

    #correct data and evaluate missing
    self.h /= 1.e3 #from nm in um
    if "A_c" in self.indicies         : self.A_c /= 1.e6  #from nm in um
    if "slope" in self.indicies       : self.slope /= 1.e3 #from N/m in mN/um
    if "slopeSupport" in self.indicies: self.slopeSupport /= 1.e3 #from N/m in mN/um
    if 'h_c' in self.indicies         : self.h_c /= 1.e3  #from nm in um
    if 'hRaw' in self.indicies        : self.hRaw /= 1.e3  #from nm in um
    if not "k2p" in self.indicies and 'slope' in self.indicies:
      self.k2p = self.slope * self.slope / self.p[self.valid]
    return success


  def loadHysitron(self, fileName, plotContact=False):
    """
    Load Hysitron hld or txt file for processing, only contains one test

    Args:
       fileName: file name
       plotContact: plot intial contact identification (use this method for access)
    """
    from io import StringIO
    self.fileName = fileName
    inFile = open(self.fileName, 'r',encoding='iso-8859-1')
    #### HLD FILE ###
    if self.fileName.endswith('.hld'):
      line = inFile.readline()
      if not "File Version: Hysitron" in line:
        #not a Hysitron file
        inFile.close()
        return False
      if self.verbose>1: print("Open Hysitron file: "+self.fileName)

      #read meta-data
      prefact = [0]*6
      segmentTime = []
      segmentDeltaP = []
      segmentPoints = []
      while True:
        line = inFile.readline()
        label = line.split(":")[0]
        try:
          data = line.split(":")[1].split(" ")
          value = float(data[1])
          #if len(data)>2: unit  = data[2]
          #else:           unit  = ""
        except:
          value = line.split(":")[1].rstrip()
          #unit  = ""
        if label == "Sample Approach Data Points": break
        if label == "Machine Comp": self.compliance = value #assume nm/uN = um/mN
        if label == "Tip C0":       prefact[0] = value #nm^2/nm^2
        if label == "Tip C1":       prefact[1] = value #nm^2/nm
        if label == "Tip C2":       prefact[2] = value #nm^2/nm^0.5
        if label == "Tip C3":       prefact[3] = value #nm^2/nm^0.25
        if label == "Tip C4":       prefact[4] = value #nm^2/nm^0.125
        if label == "Tip C5":       prefact[5] = value #nm^2/nm^0.0625
        if label == "Contact Threshold": forceTreshold = value/1.e3 #uN
        if label == "Drift Rate":   self.meta['drift_rate'] = value/1.e3 #um/s
        if label == "Number of Segments"  : numSegments  = value
        if label == "Segment Begin Time"  : segmentTime.append(value)
        if label == "Segment Begin Demand": pStart     = value
        if label == "Segment End Demand"  : segmentDeltaP.append( (value-pStart)/1.e3 ) #to mN
        if label == "Segment Points"      : segmentPoints.append(int(value))
        if label == "Time Stamp"          : self.timeStamp = ":".join(line.rstrip().split(":")[1:])
      self.tip.prefactors = prefact
      self.tip.prefactors.append('iso')
      if (numSegments!=len(segmentTime)) or (numSegments!=len(segmentDeltaP)):
        print("*ERROR*", numSegments,len(segmentTime),len(segmentDeltaP ) )
      segmentDeltaP = np.array(segmentDeltaP)
      segmentPoints = np.array(segmentPoints)
      segmentTime   = np.array(segmentTime)

      #read approach data
      line = inFile.readline() #Time_s  MotorDisp_mm    Piezo Extension_nm"
      data = ""
      for idx in range(int(value)):
        data +=inFile.readline()
      if (len(data)>1):
        dataApproach = np.loadtxt( StringIO(str(data))  )

      #read drift data
      value = inFile.readline().split(":")[1]
      line = inFile.readline()  #Time_s	Disp_nm",value
      data = ""
      for idx in range(int(value)):
        data +=inFile.readline()
      if (len(data)>1):
        self.dataDrift = np.loadtxt( StringIO(str(data))  )
        self.dataDrift[:,1] /= 1.e3  #into um

      #read test data
      value = inFile.readline().split(":")[1]
      line = inFile.readline() #Time_s	Disp_nm	Force_uN	LoadCell_nm	PiezoDisp_nm	Disp_V	Force_V	Piezo_LowV
      data = ""
      for idx in range(int(value)):
        data +=inFile.readline()
      dataTest = np.loadtxt( StringIO(str(data))  )
      #store data
      self.t = dataTest[:,0]
      self.h = dataTest[:,1]/1.e3
      self.p = dataTest[:,2]/1.e3

      # create loading-holding-unloading cycles
      listLoading = np.where(segmentDeltaP>0.1 )[0]
      listUnload  = np.where(segmentDeltaP<-0.1)[0]
      segmentPoints  -= 1                     #since the first / last point of each segment are double in both segments
      segmentPoints[0]+=1
      segPnts   = np.cumsum(segmentPoints)
      self.iLHU = []
      for idx in range(len(listLoading)):
        iSurface = segPnts[listLoading[idx]-1]+1
        iLoad    = segPnts[listLoading[idx]]
        iHold    = segPnts[listUnload[idx]-1]+1
        iUnload  = segPnts[listUnload[idx]]
        self.iLHU.append( [iSurface,iLoad,iHold,iUnload] )

    #### TXT FILE ###
    if self.fileName.endswith('.txt'):
      line0 = inFile.readline()
      line1 = inFile.readline()
      line2 = inFile.readline()
      line3 = inFile.readline()
      self.meta = {'measurementType': 'Hysitron Indentation TXT', 'dateMeasurement':line0.strip()}
      if line1 != "\n" or "Number of Points" not in line2 or not "Depth (nm)" in line3:
        inFile.close()
        return False #not a Hysitron file
      if self.verbose>1: print("Open Hysitron file: "+self.fileName)
      dataTest = np.loadtxt(inFile)
      #store data
      self.t = dataTest[:,2]
      self.h = dataTest[:,0]/1.e3
      self.p = dataTest[:,1]/1.e3
      #set unknown values
      forceTreshold = 0.25 #250uN
      self.identifyLoadHoldUnload()

    #correct data
    #flatten intial section of retraction
    idxMinH  = np.argmin(self.h)
    self.p[:idxMinH] = self.p[idxMinH]
    self.h[:idxMinH] = self.h[idxMinH]
    idxMask = int( np.where(self.p>forceTreshold)[0][0])
    fractionMinH = 0.5
    hFraction    = (1.-fractionMinH)*self.h[idxMinH]+fractionMinH*self.h[idxMask]
    idxMask = np.argmin(np.abs(self.h-hFraction))
    if idxMask>2:
      mask     = np.zeros_like(self.h, dtype=bool)
      mask[:idxMask] = True
      fit = np.polyfit(self.h[mask],self.p[mask],1)
      self.p -= np.polyval(fit,self.h)

      #use force signal and its threshold to identify surface
      #Option: use lowpass-filter and then evaluate slope: accurate surface identifaction possible, however complicated
      #Option: use lowpass-filter and then use force-threshold: accurate surface identifaction, however complicated
      #Best: use medfilter or wiener on force signal and then use force-threshold: accurate and easy
      #see also Bernado_Hysitron/FirstTests/PhillipRWTH/testSignal.py
      pZero    = np.average(self.p[mask])
      pNoise   = max(pZero-np.min(self.p[mask]), np.max(self.p[mask])-pZero )
      #from initial loading: back-extrapolate to zero force
      maskInitLoad = np.logical_and(self.p>pZero+pNoise*2. , self.p<forceTreshold)
      maskInitLoad[np.argmax(self.p):] = False
      fitInitLoad  = np.polyfit(self.p[maskInitLoad],self.h[maskInitLoad],2)  #inverse h-p plot -> next line easier
      hZero        = np.polyval(fitInitLoad, pZero)
      ## idx = np.where(  self.p>(pZero+pNoise)  )[0][0] OLD SYSTEM NOT AS ACCURATE, better fitInitLoad
      if plotContact:
        """
        attempt to use dp/dt for contact identification
        from scipy import signal
        h = signal.wiener(self.h,5)
        p = signal.wiener(self.p,5)
        dpdt = (p[1:]-p[:-1])/(h[1:]-h[:-1])
        _, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax2.plot(h[:-1],dpdt)
        ax1.plot(h,p,'r')
        ax2.set_ylim([0,1000])
        plt.show()
        """
        plt.axhline(pZero,c='g', label='pZero')
        plt.axhline(pZero+pNoise,c='g',linestyle='dashed', label='pNoise')
        plt.axvline(self.h[idxMinH],c='k')
        plt.axvline(self.h[idxMask],c='k')
        plt.plot(self.h,self.p)
        plt.plot(self.h[mask],self.p[mask], label='used for pNoise')
        plt.plot(self.h[maskInitLoad],self.p[maskInitLoad], label='used for backfit')
        plt.axvline(hZero,c='r',linestyle='dashed',label='Start')
        plt.plot(hZero,pZero, "ro", label="Start")
        plt.legend(loc=0)
        plt.ylim(np.min(self.p), self.p[idx]+forceTreshold )
        plt.xlim(-0.1,self.h[idx]+0.05)
        plt.show()
        print ("Debug pZero and pNoise:",pZero,pNoise)
    else:
      print("Error", forceTreshold,np.where(self.p>forceTreshold)[0][:10])
      pZero = 0
      idx   = 0
    if True:
      ## self.t -= self.t[idx] #do not offset time since segment times are given
      self.h -= hZero
      self.p -= pZero
    inFile.close()
    return True



  def loadMicromaterials(self, fileName, plotContact=False):
    """
    Load Micromaterials txt file for processing, contains only one test

    Args:
       fileName: file name or file-content
       plotContact: plot intial contact identification (use this method for access)
    """
    try:            #file-content given
      dataTest = np.loadtxt(fileName)  #exception caught
      if not isinstance(fileName, io.TextIOWrapper):
        self.fileName = fileName
        if self.verbose>1: print("Open Micromaterials file: "+self.fileName)
        self.meta = {'measurementType': 'Micromaterials Indentation TXT'}
    except:
      if self.verbose>1:
        print("Is not a Micromaterials file")
      return False

    #store data
    self.t = dataTest[:,0]
    self.h = dataTest[:,1]/1.e3
    self.p = dataTest[:,2]
    self.valid = np.ones_like(self.t, dtype=bool)
    #set unknown values
    forceTreshold = 0.25 #250uN
    self.identifyLoadHoldUnload()
    return True


  def loadMicromaterialsZip(self, fileName):
    """
    Initialize zip-file file for processing

    Args:
       fileName: file name
    """
    if self.verbose>1: print("Open Micromaterials zip-file: "+fileName)
    self.datafile = ZipFile(fileName)
    self.testList = self.datafile.namelist()
    self.fileName = fileName
    self.nextMicromaterialsTest()
    self.meta = {'measurementType': 'Micromaterials Indentation ZIP'}
    return True


  def nextMicromaterialsTest(self):
    """
    Go to next file in zip-file
    """
    if self.vendor!=Vendor.Micromaterials: return False #cannot be used
    if len(self.testList)==0: return False   #no sheet left
    self.testName = self.testList.pop(0)
    with self.datafile.open(self.testName) as myFile:
      txt = io.TextIOWrapper(myFile, encoding="utf-8")
      success = self.loadMicromaterials(txt)
    return success


  def loadFischerScope(self,fileName):
    """
    Initialize txt-file from Fischer-Scope for processing

    Args:
      fileName: file name
    """
    self.meta = {'date':[], 'shape correction':[], 'coordinate x':[], 'coordinate y':[],
            'work elastic':[], 'work nonelastic':[], 'EIT/(1-vs^2) [GPa]':[], 'HIT [N/mm]':[],
            'HUpl [N/mm]': [], 'hr [um]':[], 'hmax [um]':[], 'Compliance [um/N]':[],
            'epsilon':[], 'fit range': []}
    self.workbook = []
    self.testList = []
    self.fileName = fileName
    block = None
    with open(fileName,'r',encoding='iso-8859-1') as fIn:
      # read initial lines and initialialize
      line = fIn.readline()
      if ".hap	Name of the application" not in line:
        print("Not a Fischer Scope")
        return False
      identifier = line.split()[0]
      temp = fIn.readline()
      self.meta['Indent_Type'] = fIn.readline().split()[0]
      self.meta['Indent_F'] = ' '.join( fIn.readline().split()[2:] )
      self.meta['Indent_C'] = ' '.join( fIn.readline().split()[2:] )
      self.meta['Indent_R'] = ' '.join( fIn.readline().split()[2:] )
      #read all lines after initial lines
      for line in fIn:
        pattern = identifier+r"   \d\d\.\d\d\.\d\d\d\d  \d\d:\d\d:\d\d"
        if re.match(pattern, line) is not None:
          ## finish old individual measurement
          if block is not None:
            df = pd.DataFrame(np.array(block), columns=['F','h','t'] )
            self.workbook.append(df)
          ## start new  individual measurement
          block = []
          self.meta['date'] += [' '.join(line.split()[-2:])]
          self.testList.append('_'.join(line.split()[-2:]))
        elif line.startswith('Indenter shape correction:'):
          self.meta['shape correction'] += [line.split()[-1]]
        elif 'x=  ' in line and 'y=  ' in line:
          self.meta['coordinate x'] += [float(line.split()[1])]
          self.meta['coordinate y'] += [float(line.split()[3])]
        elif line.startswith('We	['):
          self.meta['work elastic'] += [line.split()[-1]]
        elif line.startswith('Wr	['):
          self.meta['work nonelastic'] += [line.split()[-1]]
        elif line.startswith('EIT/(1-vs^2)	[GPa]'):
          self.meta['EIT/(1-vs^2) [GPa]'] += [float(line.split()[-1])]
        elif line.startswith('HIT	[N/mm'):
          self.meta['HIT [N/mm]'] += [float(line.split()[-1])]
        elif line.startswith('HUpl	[N/mm'):
          self.meta['HUpl [N/mm]'] += [float(line.split()[-1])]
        elif line.startswith('hr	['):
          self.meta['hr [um]'] += [float(line.split()[-1])]
        elif line.startswith('hmax	['):
          self.meta['hmax [um]'] += [float(line.split()[-1])]
        elif line.startswith('Compliance	['):
          self.meta['Compliance [um/N]'] += [float(line.split()[-1])]
        elif 'Epsilon =' in line:
          self.meta['epsilon'] += [float(line.split()[-1])]
          self.meta['fit range'] += [' '.join(line.split()[:-3])]
        elif re.match(r'^\d+,\d+\s\d+,\d+\s\d+,\d+$', line):
          line = line.replace(',','.').strip()
          block.append( [float(i) for i in line.split()] )
      ## add last dataframe
      df = pd.DataFrame(np.array(block), columns=['F','h','t'] )
      self.workbook.append(df)
    if self.verbose>1:
      print("Meta information:",self.meta)
      print("Number of measurements read:",len(self.workbook))
    self.meta['measurementType'] = 'Fischer-Scope Indentation TXT'
    if self.meta['Indent_F'].startswith('ESP'):
      self.method = Method.MULTI
    else:
      self.method = Method.ISO
    self.nextFischerScopeTest()
    return True


  def nextFischerScopeTest(self):
    """
    Go to next dataframe
    """
    df = self.workbook.pop(0)
    self.testName = self.testList.pop(0)
    self.t = np.array(df['t'])
    self.h = np.array(df['h'])
    self.p = np.array(df['F'])
    self.valid = np.ones_like(self.t, dtype=bool)
    self.identifyLoadHoldUnload()
    return


  #@}
  ##
  # @name Output and plotting
  #@{
  def getDataframe(self):
    """
    previous name: getResult
    save all results to dataframe variable for easy proccessing
    """
    dfAll = pd.DataFrame()
    if self.method == Method.CSM:
        i = -1 # only last value is saved
        results = {"S_mN/um":self.slope[i], "hMax_um":self.h[self.valid][i], "pMax_mN":self.p[self.valid][i],\
                  "redE_GPa":self.modulusRed[i], "A_um2":self.A_c[i], "hc_um":self.h_c[i], "E_GPa":self.modulus[i],\
                  "H_GPa":self.hardness[i],"segment":str(i+1)}
        df = pd.DataFrame(results, index=[self.fileName])
        dfAll = dfAll.append(df)
    else:
      for i in range(len(self.slope)):
        results = {"S_mN/um":self.slope[i], "hMax_um":self.h[self.valid][i], "pMax_mN":self.p[self.valid][i],\
                  "redE_GPa":self.modulusRed[i], "A_um2":self.A_c[i], "hc_um":self.h_c[i], "E_GPa":self.modulus[i],\
                  "H_GPa":self.hardness[i],"segment":str(i+1)}
        df = pd.DataFrame(results, index=[self.fileName+'-'+self.testName])
        dfAll = dfAll.append(df)
    dfAll['method'] = "Python"
    if "timeStamp" in dfAll.columns:
      dfAll['timeStamp'] = self.timeStamp
    dfAll = dfAll.rename(index=str, columns={'drift':'drift_um/s'})
    return dfAll


  def getDictionary(self):
    """
    save last result to dictionary variable for easy proccessing
    """
    if len(self.slope)>1:
      print("Error in getDictionary")
      return None
    i = -1
    results = {"S_mN/um":self.slope[i], "hMax_um":self.h[self.valid][i], "pMax_mN":self.p[self.valid][i],\
              "redE_GPa":self.modulusRed[i], "A_um2":self.A_c[i], "hc_um":self.h_c[i], "E_GPa":self.modulus[i],\
              "H_GPa":self.hardness[i],"segment":str(i+1)}
    results.update(self.meta)
    return results


  def plotTestingMethod(self):
    """
    plot testing method
    """
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(self.t, self.p,'C0')
    ax2.plot(self.t, self.h,'C1')
    for mask in self.iLHU:
      ax1.plot(self.t[mask][0], self.p[mask][0], 'C0s')
      ax1.plot(self.t[mask][1], self.p[mask][1], 'C0x')
      ax1.plot(self.t[mask][2], self.p[mask][2], 'C0+')
      ax1.plot(self.t[mask][3], self.p[mask][3], 'C0o')
    ax1.plot(self.t[self.iDrift], self.p[self.iDrift], 'k.')
    ax1.axhline(0,color='C0', linestyle='dashed')
    ax2.axhline(0,color='C1', linestyle='dashed')
    ax1.set_xlabel(r"time [$s$]")
    ax2.set_ylabel(r"depth [$\mu m$]", color='C1', fontsize=14)
    ax1.set_ylabel(r"force [$mN$]", color='C0', fontsize=14)
    plt.grid()
    plt.show()
    return


  def plot(self, saveFig=False, show=True):
    """
    Plot force-depth curve with all data

    Args:
      saveFig: save plot to file [use known filename plus extension png]
      show: show figure, else do not show
    """
    if len(self.slope)==1 and self.verbose>1:
      print("Stiffness:"+str(round(self.slope[0],1))     +"mN/um   hMax:"+str(round(self.h[self.valid][0],4))+"um    pMax:"+str(round(self.p[self.valid][0],2))+"mN")
      print("E*:       "+str(round(self.modulusRed[0],1))+"GPa     A:   "+str(round(self.A_c[0],4))+          "um2    hc: "+str(round(self.h_c[0],4))+"um")
      print("E:        "+str(round(self.modulus[0],1))   +"GPa     H:   "+str(round(self.hardness[0],1))+     "GPa")
    f, ax = plt.subplots()
    ax.axhline(0,ls="dashed",c='k')
    ax.axvline(0,ls="dashed",c='k')
    ax.plot(self.h,self.p)
    if self.method != Method.CSM:
      _, _, maskUnload, optPar, _ = self.stiffnessFromUnloading(self.p, self.h)
      h_, p_ = self.h[maskUnload], self.p[maskUnload]
      ax.plot(self.h[maskUnload], self.UnloadingPowerFunc(self.h[maskUnload],*optPar), 'C1', label='fit powerlaw' )
      ax.plot(self.h[self.valid],self.p[self.valid],"or",label="max", markersize=10)
      ax.plot(self.h_c, np.zeros_like(self.h_c),"ob", label="h_c", markersize=10)
      if len(self.h_c)<2:
        ax.plot(h_[0],p_[0],'og',)
        ax.plot(h_[-1],p_[-1],'og', label="fit domain")
        Sn= self.p[self.valid]-self.slope*self.h[self.valid]
        h_ = np.linspace(self.h_c,self.h[self.valid],10)
        ax.plot(h_,   self.slope*h_+Sn, 'r--', lw=2, label='stiffness')
      ax.legend(loc=0, numpoints=1)
    else:
      ax.plot(self.h[self.iLHU[0]],self.p[self.iLHU[0]],"or",label="specific", markersize=10)
    ax.set_xlim(left=-0.03)
    ax.set_xlabel(r"depth [$\mu m$]")
    ax.set_ylabel(r"force [$mN$]")
    if saveFig:
      plt.savefig(self.fileName.split('.')[0]+".png", dpi=150, bbox_inches='tight')
    if show:
      plt.show()
    return ax


  def plotAsDepth(self, property, saveFig=False, hvline=None):
    """
    Plot as function of depth either Young's modulus, hardness,
    stiffnessSquaredForce, ContactDepth, Contact Area, reducedModulus

    Args:
      property: what to plot on y-axis [E,H,K,K2P,h_c,A_c,Ered]
      saveFig: save plot to file [use known filename plus extension png]
    """
    if not isinstance(property, str):
      print("**ERROR plotAsDepth: property=[E,H,K,K2P,h_c,A_c,Ered]")
      return
    if hvline is not None:
      plt.axhline(hvline, c='k')
    if   property == "E":
      plt.plot(self.h[self.valid], self.modulus, "o")
      plt.ylabel("Young's modulus [GPa]")
    elif property == "Ered":
      plt.plot(self.h[self.valid], self.modulusRed, "o")
      plt.ylabel("reduced Young's modulus [GPa]")
    elif property == "H":
      plt.plot(self.h[self.valid], self.hardness, "o")
      plt.ylabel("Hardness [GPa]")
    elif property == "K":
      plt.plot(self.h[self.valid], self.slope, "o")
      plt.ylabel("Stiffness [kN/m]")
    elif property == "K2P":
      if not hasattr(self, 'k2p'):
        self.k2p = np.array(self.slope)*np.array(self.slope)/np.array(self.p[self.valid])
      plt.plot(self.h[self.valid], self.k2p, "C0o")
      mask = self.h[self.valid]>0.1
      fit = np.polyfit(self.h[self.valid][mask], self.k2p[mask],1)
      plt.plot(self.h[self.valid], np.polyval(fit,self.h[self.valid]), 'C1-')
      plt.axvline(0.1, linestyle='dashed',color='C1')
      plt.ylabel(r"Stiffness Squared Over Load [GPa]")
    elif property == "h_c":
      plt.plot(self.h[self.valid], self.h_c, "o")
      plt.ylabel(r"Contact depth [$\mu m$]")
    elif property == "A_c":
      plt.plot(self.h[self.valid], self.A_c, "o")
      plt.ylabel(r"Contact area [$\mu m^2$]")
    else:
      print("Unknown property")
      return
    plt.xlabel(r"depth "+r'$[\mu m]$')
    plt.show()
    return


  #@}
  ##
  # @name CALIBRATION METHOD
  #@{
  def calibration(self,eTarget=72.0,numPolynomial=3,critDepth=1.0,critForce=0.0,plotStiffness=False,plotTip=False, **kwargs):
    """
    Calibrate by first frame-stiffness and then area-function calibration

    Args:
       eTarget: target Young's modulus (not reduced), nu is known
       numPolynomial: number of area function polynomial; if None: return interpolation function
       critDepth: frame stiffness: what is the minimum depth of data used; area function: what is the maximum depth of data used
                  (if deep data is used for area function, this data can scew the area function)
       critForce: frame stiffness: what is the minimum force used for fitting
       plotStiffness: plot stiffness graph with compliance
       pltTip: plot tip shape after fitting
    """
    constantTerm = kwargs.get('constantTerm', False)
    frameCompliance, res = self.calibrateStiffness(critDepth=critDepth,critForce=critForce,plotStiffness=plotStiffness)

    print("\nStart tip-area function fitting by reading the data and using compliance")
    ## re-create data-frame of all files
    temp = {'method': self.method, 'onlyLoadingSegment': self.onlyLoadingSegment}
    self.__init__(self.fileName, nuMat=self.nuMat, verbose=self.verbose)
    self.tip.compliance = frameCompliance
    for item in temp:
      setattr(self, item, temp[item])
    if self.method==Method.CSM:
      self.nextAgilentTest(newTest=False)  #rerun to ensure that onlyLoadingSegment used
      slope = None
      while True:
        self.analyse()
        if slope is None:
          slope = self.slope
          h     = self.h[self.valid]
          p     = self.p[self.valid]
        else:
          slope = np.hstack((slope, self.slope))
          h     = np.hstack((h,     self.h[self.valid]))
          p     = np.hstack((p,     self.p[self.valid]))
        if not self.testList: break
        self.nextTest()
    else:
      dfNew = pd.DataFrame()
      while True:
        self.analyse()
        dfNew = dfNew.append(self.getDataframe())
        if len(self.testList)==0: break
        self.nextTest()
      slope = np.array(dfNew['S_mN/um'])
      ## verify (only makes sense for non-CSM because CSM calculates stiffness like this)
      totalCompliance = 1./dfAll['S_mN/um']
      contactStiffness = 1./(totalCompliance-frameCompliance) #mN/um
      print("Info: difference direct-indirect stiffness %",round(np.linalg.norm((slope-contactStiffness)/slope)*100,2),"%. Should be small")
      slope = np.array(dfNew['S_mN/um'])
      h     = dfNew['hMax_um']
      p     = dfNew['pMax_mN']

    ## fit shape function
    #reverse OliverPharrMethod to determine area function
    EredGoal = self.ReducedModulus(eTarget, self.nuMat)
    A = np.array( np.power( slope  / (2.0*EredGoal/np.sqrt(np.pi))  ,2))
    h_c = np.array( h - self.beta*p/slope )
    #calculate shape function as interpolation of 30 points (log-spacing)
    #  first calculate the  savgol-average using a adaptive window-size
    if numPolynomial is None:
      # use interpolation function using random points
      data = np.vstack((h_c,A))
      data = data[:, data[0].argsort()]
      windowSize = int(len(A)/20) if int(len(A)/20)%2==1 else int(len(A)/20)-1
      output = savgol_filter(data,windowSize,3)
      interpolationFunct = interpolate.interp1d(output[0,:],output[1,:])
      h_c_ = np.logspace(np.log(0.0001),np.log(np.max(output[0,:])),num=30,base=np.exp(1))
      A_c_ = interpolationFunct(h_c_)
      interpolationFunct = interpolate.interp1d(h_c_, A_c_)
      del output, data
    else:
      #It is possible to crop only interesting contact depth: h_c>1nm
      # A = A[h_c>0.001]
      # h_c = h_c[h_c>0.001]
      if constantTerm:
        appendix = 'isoPlusConstant'
      else:
        appendix = 'iso'
      def fitFunct(params):     #error function
        self.tip.prefactors = [params[x].value for x in params]+[appendix]
        A_temp = self.tip.areaFunction(h_c)                #use all datapoints as critDepth is for compliance plot
        residual     = np.abs(A-A_temp)/len(A)             #normalize by number of points
        return residual
      # Parameters, 'value' = initial condition, 'min' and 'max' = boundaries
      params = lmfit.Parameters()
      params.add('m0', value= 24.3, min=10.0, max=60.0)
      for idx in range(1,numPolynomial):
        startVal = np.power(100,idx)
        params.add('m'+str(idx), value= startVal/1000, min=-startVal*100, max=startVal*100)
      if constantTerm:
        params.add('c',  value= 20, min=0.5, max=300.0) ##all prefactors are in nm, this has to be too
      # do fit, here with leastsq model; args=(h_c, A)
      result = lmfit.minimize(fitFunct, params, max_nfev=10000)
      self.tip.prefactors = [result.params[x].value for x in result.params]+[appendix]
      print("\nTip shape:")
      print("iterated prefactors",[round(i,1) for i in self.tip.prefactors[:-1]])
      stderr = [result.params[x].stderr for x in result.params]
      print("  standard error",['NaN' if x is None else round(x,2) for x in stderr])
      res['tip prefact factors and std.error'] = [self.tip.prefactors[:-1],stderr]

    if plotTip:
      if numPolynomial is None:
        self.tip.setInterpolationFunction(interpolationFunct)
      rNonPerfect = np.sqrt(A/np.pi)
      plt.plot(rNonPerfect, h_c,'.')
      self.tip.plotIndenterShape(maxDepth=1.5)
      #Error plot
      plt.plot(h_c,(A-self.tip.areaFunction(h_c))/A,'o',markersize=2)
      plt.axhline(0,color='k',linewidth=2)
      plt.xlabel("Depth [$\mu m$]")
      plt.ylabel("Relative area error")
      plt.ylim([-0.1,0.1])
      plt.xlim(left=0)
      plt.yticks([-0.1,-0.05,0,0.05,0.1])
      plt.show()

    #rerun everything with calibrated area function to see
    self.__init__(self.fileName, nuMat=self.nuMat, verbose=self.verbose, tip=self.tip)
    for item in temp:
      setattr(self, item, temp[item])
    if self.method==Method.CSM:
      self.nextAgilentTest(newTest=False)  #rerun to ensure that onlyLoadingSegment used
    else:
      ## create data-frame of all files
      dfAll = pd.DataFrame()
      while True:
        self.analyse()
        dfAll = dfAll.append(self.getDataframe())
        if not self.testList: break
        self.nextTest()
      ## output representative values
      maskPrint = dfAll['pMax_mN'] > 0.95*np.max(dfAll['pMax_mN'])
      res['End Max depth: ave,stderr']=[dfAll['hMax_um'][maskPrint].mean(),dfAll['hMax_um'][maskPrint].std()/dfAll['hMax_um'][maskPrint].count()]
      res['End MeasStiff: ave,stderr']=[dfAll['S_mN/um'][maskPrint].mean(),dfAll['S_mN/um'][maskPrint].std()/dfAll['S_mN/um'][maskPrint].count()]
      res['End A_c: ave,stderr']=[      dfAll['A_um2'][maskPrint].mean(),  dfAll['A_um2'][maskPrint].std()/  dfAll['A_um2'][maskPrint].count()]
      res['End h_c: ave,stderr']=[      dfAll['hc_um'][maskPrint].mean(),  dfAll['hc_um'][maskPrint].std()/  dfAll['hc_um'][maskPrint].count()]
      res['End E: ave,stderr']=[  dfAll['E_GPa'].mean(),   dfAll['E_GPa'].std()/   dfAll['E_GPa'].count()]
      res['End E_r: ave,stderr']=[dfAll['redE_GPa'].mean(),dfAll['redE_GPa'].std()/dfAll['redE_GPa'].count()]
      res['End H: ave,stderr']=[  dfAll['H_GPa'].mean(),   dfAll['H_GPa'].std()/   dfAll['H_GPa'].count()]
    if numPolynomial is None:
      return res, interpolationFunct
    return res


  def calibrateStiffness(self,critDepth=1.0,critForce=0.0,plotStiffness=True, returnAxis=False):
    """
    Calibrate by first frame-stiffness from K^2/P of individual measurement

    Args:
       critDepth: frame stiffness: what is the minimum depth of data used
       critForce: frame stiffness: what is the minimum force used for fitting
       plotStiffness: plot stiffness graph with compliance
       returnAxis: return axis of plot
    """
    print("Start compliance fitting")
    ## output representative values
    res = {}
    if self.method==Method.CSM:
      x, y, mask = None, None, None
      while True:
        self.analyse()
        if x is None:
          x    = 1./np.sqrt(self.p[self.valid])
          y    = 1./self.slope
          mask = self.h[self.valid]
        else:
          x =    np.hstack((x,    1./np.sqrt(self.p[self.valid]) ))
          y =    np.hstack((y,    1./self.slope))
          mask = np.hstack((mask, self.h[self.valid]))
        if not self.testList: break
        self.nextTest()
      mask = np.logical_and(mask>critDepth, x<1./np.sqrt(critForce))
      maskPrint = []
    else:
      ## create data-frame of all files
      dfAll = pd.DataFrame()
      while True:
        self.analyse()
        dfAll = dfAll.append(self.getDataframe())
        if not self.testList: break
        self.nextTest()
      maskPrint = dfAll['pMax_mN'] > 0.95*np.max(dfAll['pMax_mN'])
      res['Input Max force: ave,stderr'] = [dfAll['pMax_mN'][maskPrint].mean(),dfAll['pMax_mN'][maskPrint].std()/dfAll['pMax_mN'][maskPrint].count()]
      res['Input Max depth: ave,stderr'] = [dfAll['hMax_um'][maskPrint].mean(),dfAll['hMax_um'][maskPrint].std()/dfAll['hMax_um'][maskPrint].count()]
      res['Input MeasStiff: ave,stderr'] = [dfAll['S_mN/um'][maskPrint].mean(),dfAll['S_mN/um'][maskPrint].std()/dfAll['S_mN/um'][maskPrint].count()]
      ## determine compliance by intersection of 1/sqrt(p) -- compliance curve
      x = 1./np.sqrt(dfAll['pMax_mN'])
      y = 1./dfAll['S_mN/um']
      mask = dfAll['hMax_um'] > critDepth
      mask = np.logical_and(mask, dfAll['pMax_mN'] > critForce)
    if(len(mask[mask])==0):
      print("ERROR too much filtering, no data left. Decrease critForce and critDepth")
      return None

    param, covM = np.polyfit(x[mask],y[mask],1, cov=True)
    print("fit f(x)=",round(param[0],5),"*x+",round(param[1],5))
    frameStiff = 1./param[1]
    frameCompliance = param[1]
    print("  frame compliance: %8.4e um/mN = %8.4e m/N"%(frameCompliance,frameCompliance/1000.))
    stderrPercent = np.abs( np.sqrt(np.diag(covM)[1]) / param[1] * 100. )
    print("  compliance and stiffness standard error in %:",round(stderrPercent,2) )
    res['Stiffness and error in %']=[frameStiff,stderrPercent]
    print("  frame stiffness: %6.0f mN/um = %6.2e N/m"%(frameStiff,1000.*frameStiff))
    self.tip.compliance = frameCompliance

    if plotStiffness:
      f, ax = plt.subplots()
      ax.plot(x[~mask], y[~mask], 'o', color='#165480', fillstyle='none', markersize=1, label='excluded')
      ax.plot(x[mask], y[mask],   'C0o', markersize=5, label='for fit')
      x_ = np.linspace(0, np.max(x)*1.1, 50)
      y_ = np.polyval(param, x_)
      ax.plot(x_,y_,'w-')
      ax.plot(x_,y_,'C0--')
      ax.plot([0,np.min(x)/2],[frameCompliance,frameCompliance],'k')
      ax.text(np.min(x)/2,frameCompliance,'frame compliance')
      ax.set_xlabel(r"1/sqrt(p) [$mN^{-1/2}$]")
      ax.set_ylabel(r"meas. compliance [$\mu m/mN$]")
      ax.legend(loc=4)
      ax.set_ylim([0,np.max(y[mask])*1.5])
      ax.set_xlim([0,np.max(x[mask])*1.5])
      if returnAxis:
        return ax
      plt.show()
    return [frameCompliance, res]


  #@}
  ##
  # @name VERIFY METHODS
  #@{
  def verifyOneData(self):
    """
    Test one data set to ensure everything still working: OliverPharrMethod and area functions (normal and inverse)
    """
    self.tip.prefactors = [32.9049, -6418.303798, 288484.8518, -989287.0625, 103588.5588, 675977.3345, "iso"]
    print("Test CSM method, and area functions (normal and inverse)")
    #values from time=78.47sec of Cu_500muN_Creep_Vergleich
    harmStiff   = 159111.704268288/1000.   #  159111.704268288N/m
    load        = 0.491297865144331        #  0.491297865144331mN
    totalDepth  = 111.172457420282/1000.   #  111.901346020458nm
    print("   Set Poisson's ratio 0.35")
    self.nuMat = 0.35
    print("   From Agilent software")
    print("      harmStiff   = 159111.704268288 N/m")
    print("      load        = 0.49129786514433 mN")
    print("      totalDepth  = 111.901346020458 nm")
    print("      H           = 0.82150309678705 GPa")
    print("      E           = 190.257729329881 GPa")
    print("      redE        = 182.338858733495 GPa")
    print("      Stiffness Squared Over Load=51529.9093101531 GPa")
    print("      ContactArea = 598047.490101769 nm^2")
    [redE, A_c, _]  = self.OliverPharrMethod(np.array([harmStiff]), np.array([load]), np.array([totalDepth]))
    print("   Evaluated by this python method")
    print("      reducedModulus [GPa] =",round(redE[0],4),"  with error=", round((redE[0]-182.338858733495)*100/182.338858733495,4),'%')
    print("      ContactArea    [um2] =",round(A_c[0],4),"  with error=", round((A_c[0]-598047.490101769/1.e6)*100/598047.490101769/1.e6,4),'%')
    E = self.YoungsModulus(redE)
    print("      Youngs Modulus [GPa] =",round(E[0],4),"  with error=", round((E[0]-190.257729329881)*100/190.257729329881,4),'%')
    totalDepth2 = self.inverseOliverPharrMethod(np.array([harmStiff]), np.array([load]), redE)
    print("      By using inverse methods: total depth h=",totalDepth2[0], "[um]  with error=", round((totalDepth2[0]-totalDepth)*100/totalDepth,4),'%')
    print("End Test")
    return

  def verifyOneData1(self):
    """
    Test one data set to ensure everything still working: OliverPharrMethod and area functions (normal and inverse)
    """
    self.tip.prefactors = [24.8204,402.507,-3070.91,3699.87,"iso" ]
    self.tip.compliance = 1000.0/9.2358e6
    print("Test CSM method, and area functions (normal and inverse)")
    #values from time=78.47sec of FS.xls in 5Materials
    harmStiff   = 25731.8375827836/1000.   #  25731.8375827836 N/m
    load        = 0.987624311132669        #  0.987624311132669 mN
    totalDepth  = 88.2388854303261/1000.   #  88.2388854303261 nm
    print("   Set Poisson's ratio 0.18")
    self.nuMat = 0.18
    print("   From Agilent software")
    print("      harmContactStiff = 25731.8375827836 N/m")
    print("      load             = 0.987624311132669 mN")
    print("      totalDepth       = 88.2388854303261 nm")
    print("      H           = 10.0514655820034 GPa")
    print("      E           = 75.1620054287519 GPa")
    print("      Stiffness Squared Over Load=670.424429535749 GPa")
    [redE, _, _]  = self.OliverPharrMethod(np.array([harmStiff]), np.array([load]), np.array([totalDepth]))
    E = self.YoungsModulus(redE)
    print("      Youngs Modulus [GPa] =",E[0],"  with error=", round((E[0]-75.1620054287519)*100/75.1620054287519,4),'%'  )
    totalDepth2 = self.inverseOliverPharrMethod(np.array([harmStiff]), np.array([load]), redE)
    print("      By using inverse methods: total depth h=",totalDepth2[0], "[um]  with error=", round((totalDepth2[0]-totalDepth)*100/totalDepth,4), '%')
    print("End Test")
    return

  def verifyReadCalc(self, plot=True):
    redE,A_c,h_c = self.OliverPharrMethod(self.slope, self.p[self.valid], self.h[self.valid])
    E = self.YoungsModulus(redE)
    H = self.p[self.valid] / A_c
    if self.method==Method.CSM:
      if plot:
        plt.plot(self.t[self.valid],self.h_c,'o',label='read')
        plt.plot(self.t[self.valid],h_c,label='calc')
        plt.legend(loc=0)
        plt.xlim(left=0)
        plt.ylim([0,np.max(self.h_c)])
        plt.title("Error in hc: {0:.2e}".format(np.linalg.norm(h_c-self.h_c)) )
        plt.show()
      else:
        print("Error in hc: {0:.2e}".format(np.linalg.norm(h_c-self.h_c)) )
    else:
      print("Error in hc: %.3e %% between %.3e and %.3e" %(abs(h_c-self.h_c)*100./h_c, h_c, self.h_c) )
    if self.method==Method.CSM:
      if plot:
        plt.plot(self.t[self.valid],self.A_c,'o',label='read')
        plt.plot(self.t[self.valid],A_c,label='calc')
        plt.legend(loc=0)
        plt.xlim(left=0)
        plt.ylim([0,np.max(self.A_c)])
        plt.title("Error in Ac: {0:.2e}".format(np.linalg.norm((A_c-self.A_c))) )
        plt.show()
      else:
        print("Error in Ac: {0:.2e}".format(np.linalg.norm((A_c-self.A_c))))
    else:
      print("Error in Ac: %.3e %% between %.3e and %.3e" %(abs(A_c-self.A_c)*100./A_c,A_c,self.A_c) )
    if self.method==Method.CSM:
      if plot:
        plt.plot(self.t[self.valid],self.modulusRed,'o',label='read')
        plt.plot(self.t[self.valid],redE,label='calc')
        plt.legend(loc=0)
        plt.xlim(left=0)
        plt.ylim([0,np.max(self.modulusRed)])
        plt.title("Error in E*: {0:.2e}".format(np.linalg.norm((redE-self.modulusRed))) )
        plt.show()
      else:
        print("Error in E*: {0:.2e}".format(np.linalg.norm((redE-self.modulusRed))))
    else:
      print("Error in E*: %.3e %% between %.3e and %.3e" %(abs(redE-self.modulusRed)*100./redE,redE,self.modulusRed) )
    if self.method==Method.CSM:
      if plot:
        plt.plot(self.t[self.valid],self.modulus,'o',label='read')
        plt.plot(self.t[self.valid],E,label='calc')
        plt.legend(loc=0)
        plt.xlim(left=0)
        plt.ylim([0,np.max(self.modulus)])
        plt.title("Error in E: {0:.2e}".format(np.linalg.norm((E-self.modulus))) )
        plt.show()
      else:
        print("Error in E: {0:.2e}".format(np.linalg.norm((E-self.modulus))))
    else:
      print("Error in E:  %.3e %% between %.3e and %.3e" %(abs(E-self.modulus)*100./E, E,self.modulus) )
    if self.method==Method.CSM:
      if plot:
        plt.plot(self.t[self.valid],self.hardness,'o',label='read')
        plt.plot(self.t[self.valid],H,label='calc')
        plt.legend(loc=0)
        plt.xlim(left=0)
        plt.ylim([0,np.max(self.hardness)])
        plt.title("Error in H: {0:.2e}".format(np.linalg.norm((H-self.hardness))) )
        plt.show()
      else:
        print("Error in H: {0:.2e}".format(np.linalg.norm((H-self.hardness))))
    else:
      print("Error in H:  %.3e %% between %.3e and %.3e" %(abs(H-self.hardness)*100./H, H,self.hardness) )
    return
  #@}




class Tip:
  def __init__(self, shape="perfect", interpFunction=None, compliance=0.0, plot=False, verbose=0):
    """
    Initialize indenter shape

    Args:
       shape: list of prefactors, default="perfect"
       interpFunction: interpolation function A_c(h_c): if given superseeds other information
       compliance: additional compliance in um/mN. sensible values: 0.0001..0.01
       plot: plot indenter shape
       verbose: how much output
    """
    #define indenter shape: could be overwritten
    if callable(interpFunction):
      self.prefactors = None
      self.interpFunction = interpFunction
    elif shape[-1]=="sphere" or shape[-1]=="iso":
      self.prefactors = shape
    elif type(shape)==list:  #assume iso
      self.prefactors = shape
      self.prefactors.append("iso")
    else:
      self.prefactors = ["perfect"]
    self.compliance = compliance

    #verify and set default values
    if self.compliance > 0.01 or self.compliance < 0.0000001:
      if compliance == 0:
        if verbose>1:
          print("*WARNING*: stiffness outside domain 1e5...1e10 N/m: infinite")
      else:
        if verbose>1:
          print("*WARNING*: stiffness outside domain 1e5...1e10 N/m:",round(1000./self.compliance) )
    if plot:
      self.plotIndenterShape()
    return

  def __repr__(self):
    outString = 'compliance: '+str(self.compliance)+';   '
    if self.prefactors is None:
      outString+= 'with interpolation function with '+str(len(self.interpFunction.x))+' points'
    else:
      outString+= 'prefactors: '+str(self.prefactors)
    return outString


  def setInterpolationFunction(self,interpFunction):
    self.interpFunction = interpFunction
    self.prefactors = None
    return

  def areaFunction(self, h):
    """
    AREA FUNCTION: from contact depth h_c calculate area

    all functions inside are using [nm]; the outside of this function uses [um]; hence at the start and end there is conversion

    prefactors:
    - "iso" type area function A=ax^2+bx^1+cx^0.5..., requires !!nm!! units
    - "perfect" type area function of a perfect Berkovich A=3*sqrt(3)*tan(65.27)^2 h_c^2 = 24.494 h_c^2
    - "sphere" type: A=pi(2Rh-h^2) h=depth, R indenter radius; for small h-> h^2=0<br>
               prefactors [-pi, 2piR] R in nm<br>
               does not account for cone at top; hence here other approach<br>

   Args:
       h [array]: contact depth in um

    Returns:
       area projected contact area in [um2]
    """
    h = h* 1000.   #starting here: all is in nm
    threshH = 1.e-3 #1pm
    h[h< threshH] = threshH
    area = np.zeros_like(h)
    if self.prefactors is None:
      self.interpFunction.bounds_error=False
      self.interpFunction.fill_value='extrapolate'
      return self.interpFunction(h/1000.)
    elif self.prefactors[-1]=='iso':
      for i in range(0, len(self.prefactors)-1):
        exponent = 2./math.pow(2,i)
        area += self.prefactors[i]*np.power(h,exponent)
        #print(i, self.prefactors[i], h,exponent, area)
    elif self.prefactors[-1]=='isoPlusConstant':
      h += self.prefactors[-2]
      for i in range(0, len(self.prefactors)-2):
        exponent = 2./math.pow(2,i)
        area += self.prefactors[i]*np.power(h,exponent)
    elif self.prefactors[-1]=='perfect':
      area = 24.494*np.power(h,2)
    elif self.prefactors[-1]=='sphere':
      radius = self.prefactors[0]*1000.
      openingAngle = self.prefactors[1]
      cos      = math.cos(openingAngle/180.0*math.pi)
      sin      = math.sin(openingAngle/180.0*math.pi)
      tan      = math.tan(openingAngle/180.0*math.pi)
      mask     = radius-h > radius*sin
      rA       = np.zeros_like(h)
      rA[mask] = np.sqrt(radius**2 - (radius-h[mask])**2 )  #spherical section
      deltaY = radius / cos			 #tapered section
      deltaX = radius-h[~mask]
      rA[~mask] = deltaY - tan*deltaX
      area = math.pi * rA * rA
    else:
      print("*ERROR*: prefactors last value does not contain type")
    area[area<0] = 0.0
    return area/1.e6


  def areaFunctionInverse(self, area, h_c0=70):
    """
    INVERSE AREA FUNCTION: from area calculate contact depth h_c

    using Newton iteration with initial guess contact depth h_c0<br>

    prefactors:
    -  "iso" type area function A=ax^2+bx^1+cx^0.5..., requires !!nm!! units
    -  "perfect" type area function of a perfect Berkovich A=3*sqrt(3)*tan(65.27)^2 h_c^2 = 24.494 h_c^2

    Args:
       area: projected contact area
       h_c0: initial guess contact depth

    Returns:
       h depth
    """
    ## define function in form f(x)-y=0
    def function(height):
      return self.areaFunction(np.array([height]))-area
    ## solve
    if self.prefactors[-1]=="iso":
      h = newton(function, h_c0)
    elif self.prefactors[-1]=="perfect":
      h = math.sqrt(area / 24.494)
    else:
      print("*ERROR*: prefactors last value does not contain type")
    return h


  def plotIndenterShape(self, maxDepth=1, steps=2000, show=True, tipLabel=None, fileName=None):
    """
    check indenter shape: plot shape function against perfect Berkovich

    analytical: perfect shape is 2.792254*x

    Args:
       maxDepth: maximum depth [um] to plot; default=10um
       steps: number of steps for plotting
       show: show figure
       tipLabel: label for this tip
       fileName: if given, save to file
    """
    zoom = 0.5
    h_c = np.linspace(0, maxDepth, steps)
    rNonPerfect = np.sqrt( self.areaFunction(h_c)/math.pi)
    rPerfect  = 2.792254*h_c
    if tipLabel is None:  tipLabel = 'this tip'
    plt.plot(rNonPerfect, h_c, '-', label=tipLabel)
    plt.plot(rPerfect,h_c, '-k', label='Berkovich')
    plt.plot(np.tan(np.radians(60.0))*h_c,h_c, '--k', label='$60^o$')
    plt.legend(loc="best")
    plt.ylabel(r'contact depth [$\mu m$]')
    plt.xlabel(r'contact radius [$\mu m$]')
    plt.xlim([0,maxDepth*4./3./zoom])
    plt.ylim([0,maxDepth/zoom])
    if show:
      plt.grid()
      if fileName is not None:
        plt.savefig(fileName, dpi=150, bbox_inches='tight')
      plt.show()
    return

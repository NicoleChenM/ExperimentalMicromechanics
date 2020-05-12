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
##### Variables: differentiate ########
# array of full length: force, time, depth, validMask, ...  [used for plotting]
# array of valid length: E,H,Ac,hc, ... [only has the length where these values are valid]
#    force[validMask] = pMax
#    all these are vectors: OliverPharr et al methods are only vector functions
#
# Coding rules:
# - Change all variables: do not keep original-depth as can be reread and makes code less readable
# TODO
# 1. see TODOs in this file
# 2. clean frame stiffness, additional compilance  : differentiate between both
#    - fitting unloading curve: assume as intial guess m=1.5
# 3. see separate file nanoIndentCalibrationMethods
#    - fit: besser gleich A_c = f(h_c) direct zu fitten weil beide bekannt
#    - better calibration methods
import numpy as np
import pandas as pd
import math, logging
import matplotlib.pyplot as plt
from enum import Enum
from scipy.optimize import differential_evolution, fmin_tnc, fmin_l_bfgs_b, curve_fit, OptimizeResult, newton

class Method(Enum):
  ISO = 1
  MULTI = 2  #multiple ISO unloadings
  CSM = 3

class Vendor(Enum):
  Agilent = 1
  HysiHLD = 2
  HysiTXT = 3
  Micromaterials = 4

class Indentation:
  def __init__(self, fileName, tip=None, verbose=1):
    """
    Initialize indentation experiment data

    Args:
       fileName: fileName to open (.xls, .hld)
       tip:  tip class to use; None=perfect
       verbose: the higher, the more information printed: 1=default, 0=print nothing
    """
    self.nuMat = 0.3
    self.nuIndent = 0.07
    self.EIndent  = 1140                                    #GPa from Oliver,Pharr Method paper
    self.beta = 0.75                                        #beta: shape coeffcient
    self.verbose = verbose
    self.method    = Method.ISO                             #iso default: csm uses different methods
    if tip is None: tip = Tip()
    self.tip = tip
    self.iLHU   = [ [-1,-1,-1,-1] ]                         #indicies of Load-Hold-Unload cycles (StartLoad-StartHold-StartUnload-EndLoad)
    self.iDrift = [-1,-1]                                   #start and end indicies of drift segment
    self.meta = {}                                          #some results come from input file, others are added by analysis
    self.slope= []

    #initialize and load first data set
    #set default parameters
    success = False
    if fileName.endswith(".xls"):
      # KLA, Agilent, Keysight, MTS
      self.vendor = Vendor.Agilent
      self.unloadPMax = 0.999
      self.unloadPMin = 0.5
      success = self.initAgilentXLS(fileName)
    if (fileName.endswith(".hld") or fileName.endswith(".txt")) and not success:
      # Hysitron
      if fileName.endswith(".hld"): self.vendor = Vendor.HysiHLD
      else                        : self.vendor = Vendor.HysiTXT
      self.unloadPMax = 0.95
      self.unloadPMin = 0.4
      success = self.initHysitron(fileName)
    if fileName.endswith(".txt") and not success:
      # Micromaterials
      self.unloadPMax = 0.95
      self.unloadPMin = 0.4
      self.vendor = Vendor.Micromaterials
      success = self.initMicromaterials(fileName)
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
    return B*np.power( (h_-hf),m)


  def stiffnessFromUnloading(self, P, h, plot=False):
    """
    Calculate single unloading stiffness from Unloading
    see G200 manual, p7-6

    Args:
       P: vector of forces
       h: vector of depth
       plot: plot results

    Results:
       stiffness, validMask [values of P,h where stiffness is determined], mask, optimalVariables
    """
    debug = False
    if self.method== Method.CSM:
      print("***WARNING Cannot calculate unloading data from CSM method")
      return None,None,None,None
    logging.info("Number of unloading segments:"+str(len(self.iLHU)))
    S         = []
    maxDeltaP = -0.01
    t = self.t - np.min(self.t)  #undo resetting zero during init
    validMask = np.zeros_like(P, dtype=np.bool)
    if plot:
      plt.plot(h,P,'-k')
    for cycle in self.iLHU:
      loadStart, loadEnd, unloadStart, unloadEnd = cycle
      maskSegment = np.zeros_like(h, dtype=np.bool)
      maskSegment[unloadStart:unloadEnd+1] = True
      maskForce   = np.logical_and(P<P[loadEnd]*self.unloadPMax, P>P[loadEnd]*self.unloadPMin)
      mask        = np.logical_and(maskSegment,maskForce)
      if plot:
        plt.plot(h[mask],P[mask],'ob')
      B0  = (P[mask][-1]-P[mask][0])/(h[mask][-1]-h[mask][0])
      hf0 = h[mask][0] - P[mask][0]/B0
      try:
        opt, _ = curve_fit(self.UnloadingPowerFunc, h[mask],P[mask], p0=[B0,hf0,1.] )
        B,hf,m = opt
      except: #if np.isnan(B):
        print("stiffnessFrommasking - Multi Unload - Fitting failed. use linear")
        B  = (P[mask][-1]-P[mask][0])/(h[mask][-1]-h[mask][0])
        hf = h[mask][0] -P[mask][0]/B
        m  = 1.
        opt= (B,hf,m)
      S_ = B*m*math.pow( (h[unloadStart]-hf), m-1)
      S.append(S_)
      validMask[unloadStart]=True
      if plot:
        Sn= P[unloadStart]-S_*h[unloadStart]
        plt.plot(h[mask],   S_*h[mask]+Sn, 'm--', lw=4)
    if plot:
      plt.show()
    return S,validMask, mask, opt


  def popIn(self, correctH=True, plot=True, removeInitialNM=2.):
    """
    Search for pop-in

    Certainty:
    - deltaSlope: higher is better (difference in elastic - plastic slope). Great indicator
    - prefactor: higher is better (prefactor of elastic curve). Great indicator
    - secondRate: lower is better (height of second largest jump). Nice indicator 0.25*deltaRate
    - covElast: lower is better. bad indicator
    - deltaH: higher is better (delta depth in jump). bad indicator
    - deltaRate: higher is better (depth rate during jump). bad indicator

    Args:
       correctH: correct depth such that curves aligend
       plot: plot pop-in curve
       removeInitialNM: remove initial nm from data as they have large scatter

    Returns:
       pop-in force, dictonary of certainty
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
      plt.xlabel('depth [um]')
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
      plt.xlabel('depth [um]')
      plt.ylabel('hardness [GPa]')
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
        print("*** WARNING: too short vector",len(h_))
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
      plt.xlabel('depth [um]')
      plt.ylabel('stiffness2/force [GPa]')
      plt.show()
    return prefactors


  def tareDepthForce(self, slopeThreshold=100, compareRead=False, plot=False):
    """
    Calculate depth, force, time from raw data

    TODO Improve surface identifation in future

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
    maskDrift = np.zeros_like(h, dtype=np.bool)
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
      ax1.set_xlabel(r"depth $[\mu m]$")
      ax1.set_ylabel(r"force $[mN]$")
      ax2.set_ylabel(r"depth $[mN/\mu m]$", color='C2')
      plt.show()
    #set newly obtained data
    self.h, self.p, self.t = h, p, t
    return


  def updateSlopes(self):
    """
    update slopes/stiffness, Young's modulus and hardness after displacement correction by:
    - drift change
    - compliance change
    """
    self.h -= self.tip.compliance*self.p
    self.slope, self.valid, _, _ = self.stiffnessFromUnloading(self.p, self.h)
    self.slope = np.array(self.slope)
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
    self.results['drift'] = drift
    if plot:
      _, ax1 = plt.subplots()
      ax2 = ax1.twinx()
      ax1.plot(t,rate*1.e3,'r')
      ax1.axhline(0.05,c='r',linestyle='dashed')
      ax2.plot(t,h*1.e3,'b')
      ax1.set_xlabel('time [s]')
      ax1.set_ylabel('drift rate [nm/s]', color='r')
      ax2.set_ylabel('depth [nm]', color='b')
      plt.show()
    return drift

  def identifyLoadHoldUnload(self):
    #identify load-hold-unloading cycles
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
        print("**ERROR in identify load-hold-unloading cycles")
        print(iSurface,iLoad,iHold,iDriftS,iDriftE, len(self.h))
    else:
      print("No hold or unloading segments")
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
  def initAgilentXLS(self, fileName):
    """
    Initialize G200 excel file for processing

    TODO Read Poisson ratio from file

    Args:
       fileName: file name
       plot: plot curve (not implemented)
    """
    self.testList = []
    self.fileName = fileName    #one file can have multiple tests
    self.indicies = {}
    self.workbook = pd.read_excel(fileName, sheet_name=None)
    tagged = []
    code = {"Load On Sample":"p", "Force On Surface":"p", "_Load":"pRaw", "Raw Load":"pRaw"\
                ,"Displacement Into Surface":"h", "_Displacement":"hRaw", "Raw Displacement": "hRaw"\
                ,"Time On Sample":"t", "Time in Contact":"t", "Time":"tTotal"\
                ,"Load vs Disp Slope":"pVsHSlope", "_Column": "Column", "_Frame": "Frame", "Stiffness":"slope" \
                ,"Contact Area":"A_c", "Contact Depth":"h_c"\
                ,"Harmonic Displacement":"hHarmonic", "Harmonic Load":"pHarmonic","Phase Angle":"phaseAngle"\
                ,"Harmonic Stiffness":"slopeInvalid", "Harmonic Contact Stiffness":"slope", "Stiffness Squared Over Load":"k2p"\
                , "Support Spring Stiffness":"slopeSupport", "Frame Stiffness": "frameStiffness"\
                ,"Hardness":"hardness", "H_IT Channel":"hardness", "Modulus": "modulus", "E_IT Channel": "modulus","Reduced Modulus":"modulusRed"\
                ,"Scratch Distance": "s", "XNanoPosition": "x", "YNanoPosition": "y", "X Position": "xCoarse", "Y Position": "yCoarse"\
                ,"TotalLateralForce": "L", "X Force": "pX", "_XForce": "pX", "Y Force": "pY", "_YForce": "pY"\
                ,"_XDeflection": "Ux", "_YDeflection": "Uy" }
    self.fullData = ['h','p','t','pVsHSlope','hRaw','pRaw','tTotal','slopeSupport']
    if self.verbose: print("=============  "+fileName+"  ============")
    for dfName in self.workbook.keys():
      df    = self.workbook.get(dfName)
      if "Test " in dfName and not "Tagged" in dfName:
        self.testList.append(dfName)
        #print "  I should process sheet |",sheet.name,"|"
        if len(self.indicies)==0:               #find index of colums for load, etc
          for cell in df.columns:
            if cell in code:
              self.indicies[code[cell]] = cell
              if self.verbose: print("     %-30s : %-20s "%(cell,code[cell]) )
            else:
              if self.verbose: print(" *** %-30s NOT USED"%cell)
          for item in self.indicies:
            if "Harmonic" in item: self.method = Method.CSM
          #reset to ensure default values are set
          if not "p" in self.indicies: self.indicies['p']=self.indicies['pRaw']
          if not "h" in self.indicies: self.indicies['h']=self.indicies['hRaw']
          if not "t" in self.indicies: self.indicies['t']=self.indicies['tTotal']
          #if self.verbose: print("   Found column names: ",sorted(self.indicies))
      if "Tagged" in dfName: tagged.append(dfName)
    if len(tagged)>0 and self.verbose: print("Tagged ",tagged)
    self.nextTest()
    return True


  def nextTest(self):
    """
    Go to next sheet in worksheet and prepare indentation data

    Data: _Raw: without frame stiffness correction, _Frame:  with frame stiffness correction (remove postscript finally)
    - only affects/applies directly depth (h) and stiffness (s)
    - modulus, hardness and k2p always only use the one with frame correction
    """
    if not self.vendor == Vendor.Agilent: return False #cannot be used
    if len(self.testList)==0: return False   #no sheet left
    if not ("t" in self.indicies) or \
        not ("p" in self.indicies) or \
        not ("h" in self.indicies) or \
        not ("slope" in self.indicies)  :
          print("WARNING: INDENTATION: Some index is missing?", self.indicies)
    self.testName = self.testList.pop(0)

    #read data and identify valid data points
    df     = self.workbook.get(self.testName)
    slope   = np.array(df[self.indicies['slope']][1:-1], dtype=np.float)
    h       = np.array(df[self.indicies['h'    ]][1:-1], dtype=np.float)
    self.validFull = np.isfinite(h)
    self.valid =  np.isfinite(slope)
    self.valid[self.valid] = slope[self.valid] > 0.0  #only valid points if stiffness is positiv
    for index in self.indicies:
      data = np.array(df[self.indicies[index]][1:-1], dtype=np.float)
      mask = np.isfinite(data)
      mask[mask] = data[mask]<1e99
      self.valid = np.logical_and(self.valid, mask)                       #adopt/reduce mask continously

    #crop to only valid data
    for index in self.indicies:
      data = np.array(df[self.indicies[index]][1:-1], dtype=np.float)
      if not index in self.fullData:
        data = data[self.valid]
      else:
        data = data[self.validFull]
      setattr(self, index, data)
    self.valid = self.valid[self.validFull]
      #print(index, len(data))

    #correct data and evaluate missing
    self.h /= 1.e3 #from nm in um
    if "A_c" in self.indicies         : self.A_c /= 1.e6  #from nm in um
    if "slope" in self.indicies       : self.slope /= 1.e3 #from N/m in mN/um
    if "slopeSupport" in self.indicies: self.slopeSupport /= 1.e3 #from N/m in mN/um
    if 'h_c' in self.indicies         : self.h_c /= 1.e3  #from nm in um
    if 'hRaw' in self.indicies        : self.hRaw /= 1.e3  #from nm in um
    if not "k2p" in self.indicies     : self.k2p = self.slope * self.slope / self.p[self.valid]
    self.identifyLoadHoldUnload()
    return True


  def initHysitron(self, fileName, plotContact=False):
    """
    Initialize Hysitron hld or txt file for processing

    Args:
       fileName: file name
       plotContact: plot intial contact identification (use this method for access)
    """
    from io import StringIO
    self.fileName = fileName
    logging.info("Open file: "+self.fileName)
    inFile = open(self.fileName, 'r',encoding='iso-8859-1')
    #### HLD FILE ###
    if self.fileName.endswith('.hld'):
      line = inFile.readline()
      if not "File Version: Hysitron" in line:
         #not a Hysitron file
        return False

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
        if label == "Drift Rate":   self.results['drift'] = value/1.e3 #um/s
        if label == "Number of Segments"  : numSegments  = value
        if label == "Segment Begin Time"  : segmentTime.append(value)
        if label == "Segment Begin Demand": pStart     = value
        if label == "Segment End Demand"  : segmentDeltaP.append( (value-pStart)/1.e3 ) #to mN
        if label == "Segment Points"      : segmentPoints.append(int(value))
        if label == "Time Stamp"          : self.timeStamp = ":".join(line.rstrip().split(":")[1:])
      self.tip.prefactors = prefact
      self.tip.prefactors.append("iso")
      if (numSegments!=len(segmentTime)) or (numSegments!=len(segmentDeltaP)):
        print("ERROR", numSegments,len(segmentTime),len(segmentDeltaP ) )
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
        return False #not a Hysitron file
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
      mask     = np.zeros_like(self.h, dtype=np.bool)
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
    return True



  def initMicromaterials(self, fileName, plotContact=False):
    """
    Initialize Micromaterials txt file for processing

    Args:
       fileName: file name
       plotContact: plot intial contact identification (use this method for access)
    """
    self.fileName = fileName
    logging.info("Open file: "+self.fileName)
    try:
      inFile = open(self.fileName, 'r',encoding='iso-8859-1')
      dataTest = np.loadtxt(inFile)
    except:
      print("Is not a Micromaterials file")
      return False
    #store data
    self.t = dataTest[:,0]
    self.h = dataTest[:,1]/1.e3
    self.p = dataTest[:,2]
    #set unknown values
    forceTreshold = 0.25 #250uN
    self.identifyLoadHoldUnload()
    self.meta = {'measurementType': 'Micromaterials Indentation TXT'}
    return True



  #@}
  ##
  # @name Output and plotting
  #@{
  def getDataframe(self):
    """
    save all results to dataframe variable for easy proccessing
    """
    dfAll = pd.DataFrame()
    if self.method == Method.CSM:
        i = -1 # TODO improve, use last segment to fit
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
        results.update(self.results)
        df = pd.DataFrame(results, index=[self.fileName])
        dfAll = dfAll.append(df)
    dfAll['method'] = "Python"
    if "timeStamp" in dfAll.columns:
      dfAll['timeStamp'] = self.timeStamp
    dfAll = dfAll.rename(index=str, columns={'drift':'drift_um/s'})
    return dfAll


  def getDictionary(self):
    """
    save all results to dictionary variable for easy proccessing
    """
    if len(self.slope)>1:
      print("Error in getDictionary")
      return None
    i = 0
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
    ax1.set_xlabel("time [sec]")
    ax2.set_ylabel("depth [um]", color='C1', fontsize=14)
    ax1.set_ylabel("force [mN]", color='C0', fontsize=14)
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
    if len(self.slope)==1:
      logging.info("Stiffness:"+str(round(self.slope[0],1))     +"mN/um   hMax:"+str(round(self.h[self.valid][0],4))+"um    pMax:"+str(round(self.p[self.valid][0],2))+"mN")
      logging.info("E*:       "+str(round(self.modulusRed[0],1))+"GPa     A:   "+str(round(self.A_c[0],4))+          "um2    hc: "+str(round(self.h_c[0],4))+"um")
      logging.info("E:        "+str(round(self.modulus[0],1))   +"GPa     H:   "+str(round(self.hardness[0],1))+     "GPa")
    f, ax = plt.subplots()
    ax.axhline(0,ls="dashed",c='k')
    ax.axvline(0,ls="dashed",c='k')
    ax.plot(self.h,self.p)
    if self.method == Method.ISO:
      _, _, maskUnload, optPar = self.stiffnessFromUnloading(self.p, self.h)
      h_, p_ = self.h[maskUnload], self.p[maskUnload]
      ax.plot(self.h[maskUnload], self.UnloadingPowerFunc(self.h[maskUnload],*optPar), 'C1', label='fit powerlaw' )
      ax.plot(self.h[self.valid],self.p[self.valid],"or",label="max", markersize=10)
      ax.plot(self.h_c, np.zeros_like(self.h_c),"ob", label="h_c")
      if len(self.h_c)<2:
        ax.plot(h_[0],p_[0],'og',)
        ax.plot(h_[-1],p_[-1],'og', label="fit domain")
        Sn= self.p[self.valid]-self.slope*self.h[self.valid]
        h_ = np.linspace(self.h_c,self.h[self.valid],10)
        ax.plot(h_,   self.slope*h_+Sn, 'r--', lw=2, label='stiffness')
      ax.legend(loc=0, numpoints=1)
    ax.set_xlim(left=-0.03)
    ax.set_xlabel("depth [um]")
    ax.set_ylabel("force [mN]")
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
      fit = np.polyfit(self.h[self.valid], self.k2p,1)
      plt.plot(self.h[self.valid], np.polyval(fit,self.h[self.valid]), 'C1-')
      plt.plot(self.h[self.valid], self.k2p, "C0o")
      plt.ylabel("Stiffness Squared Over Load [GPa]")
    elif property == "h_c":
      plt.plot(self.h[self.valid], self.h_c, "o")
      plt.ylabel("Contact depth [um]")
    elif property == "A_c":
      plt.plot(self.h[self.valid], self.A_c, "o")
      plt.ylabel("Contact area [um^2]")
    else:
      print("Unknown property")
      return
    plt.xlabel("depth "+r'$[\mu m]$')
    plt.show()
    return

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
    print("      reducedModulus [GPa] =",redE[0],"  with error=", (redE-182.338858733495)/182.338858733495)
    print("      ContactArea    [um2] =",A_c[0],"  with error=", (A_c-598047.490101769/1.e6)/598047.490101769/1.e6)
    E = self.YoungsModulus(redE)
    print("      Youngs Modulus [GPa] =",E[0],"  with error=", (E-190.257729329881)/190.257729329881)
    totalDepth2 = self.inverseOliverPharrMethod(np.array([harmStiff]), np.array([load]), redE)
    print("      By using inverse methods: total depth h=",totalDepth2[0], "[um]  with error=", (totalDepth2[0]-totalDepth)/totalDepth)
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
    print("      Youngs Modulus [GPa] =",E,"  with error=", (E-75.1620054287519)/75.1620054287519)
    totalDepth2 = self.inverseOliverPharrMethod(np.array([harmStiff]), np.array([load]), redE)
    print("      By using inverse methods: total depth h=",totalDepth2[0], "[um]  with error=", (totalDepth2[0]-totalDepth)/totalDepth)
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
  def __init__(self, shape="perfect", compliance=0.0, plot=False, verbose=0):
    """
    Initialize indenter shape

    Args:
       shape: list of prefactors, default="perfect"
       compliance: additional compliance in um/mN. sensible values: 0.0001..0.01
       plot: plot indenter shape
       verbose: how much output
    """
    #define indenter shape: could be overwritten
    if (shape[-1]=="sphere" or shape[-1]=="iso" ):
      self.prefactors = shape
    elif (type(shape)==list):  #assume iso
      self.prefactors = shape
      self.prefactors.append("iso")
    else:
      self.prefactors = ["perfect"]
    self.compliance = compliance

    #verify and set default values
    if self.compliance > 0.01 or self.compliance < 0.0000001:
      if compliance == 0:
        if verbose>1:
          print("***Warning: stiffness outside domain 1e5...1e10 N/m: infinite")
      else:
        if verbose>1:
          print("***Warning: stiffness outside domain 1e5...1e10 N/m:",round(1000./self.compliance) )
    if plot:
      self.plotIndenterShape()
    return


  def areaFunction(self, hc):
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
       hc: depth in um

    Returns:
       area projected contact area in [um2]
    """
    h = hc* 1000.   #starting here: all is in nm
    threshH = 1.e-3 #1pm
    area = np.zeros_like(h)
    h[h< threshH] = threshH
    if self.prefactors[-1]=="iso":
      for i in range(0, len(self.prefactors)-1):
        exponent = 2./math.pow(2,i)
        area += self.prefactors[i]*np.power(h,exponent)
        #print i, self.prefactors[i], h,exponent, area
    elif self.prefactors[-1]=="perfect":
      area = 24.494*np.power(h,2)
    elif self.prefactors[-1]=="sphere":
      radius = self.prefactors[0]  #BUG: *1000.
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
      print("*** ERROR: prefactors last value does not contain type")
    if type(h)==np.ndarray:  area[area<0] = 0.0
    else:                    area = max(area,0.0)
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
      print("*** ERROR: prefactors last value does not contain type")
    return h


  def plotIndenterShape(self, maxDepth=1, steps=200, show=True, tipLabel=None, fileName=None):
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
    plt.legend(loc="best")
    plt.ylabel('depth '+r'$[\mu m]$')
    plt.xlabel('contact radius '+r'$[\mu m]$')
    plt.xlim([0,maxDepth*4./3./zoom])
    plt.ylim([0,maxDepth/zoom])
    if show:
      plt.grid()
      if fileName is not None:
        plt.savefig(fileName, dpi=150, bbox_inches='tight')
      plt.show()
    return

#! /usr/bin/python -W ignore
#####################################################################
# GLSpyf.py                                                         #
# Periodogram analysis using LS, GLS or GLS including a linear term.#
# Fortran library with LS, GLS and CGLS routines is required        #
#                                                                   #
# Author: J.C. Morales                                              #
# Last update: 06-Mar-2017                                          #
#####################################################################
#
#Loading Libraries
print ('*** PERIODOGRAM ANALYSIS ***')
print (' Loading minimal libraries...\n')

import numpy as np
import matplotlib.pyplot as plt
import sys
#plt.style.use('dark_background')
#
#Definitions
#Reading the table with the radial velocity
def read_table(filename):
   '''
   Function reading a file using numpy.loadtxt
   Parameters:	filename: {string}
                   Name of the file to read, numerical values only
   Returns:     x: {float, array}
                   Data
   '''
   x=np.loadtxt(filename,unpack=True,skiprows=1,delimiter=',')
   return x
#
#
#  Generalised Lomb-Scargle periodogram
def periodGLS(x,y,dy,oversampling=10,frequency=None,createplot=True,errorbarplot=True,verbose=True,saveoutput=True):
   '''
   Compute the periodogram and fits the following function to data
      y=c0+c1*cos(x)+c2*sin(x)
   Amplitude of the radial velocity = c1**2+c2**2
   Phase of the radial velocity: x=2*pi/Period*(t-tref)=2*pi*frequency*(t-tref)
   Window is computed following Lomb-Scargle definition
   Power definition (see Zechmeister and Kurster 2009)
   Parameters: x: {float, array}
                  Array of independent variable
               y: {float, array}
                  Array of dependent variable
               dy: {float, array}
                  Array of dependent variable uncertainties
               frequency: {float, array}
                  Array of frequencies to test
               tref: {float}
                  Reference time to avoid large numeric values on times
               mfactor: {float}
                  Factor for analytic FAP
   Returns: power: {float, array}
               Array with the power of each frequency
            window: (float, array}
               Array with values of the window function
            solution: {float, array}
               Array with the solution for each frequency. The solution has 4
               elements: chi2, c1, c2, c0
            dsolution: {float, array}
               Array with the solution formal uncertainties for each frequency. The solution has 4
               elements: chi_0, c1, c2, c0
   '''
   #If frequency array is not provided, construct one based on observations
   if(frequency==None):
      #####################################################################                                                                   #
      # Basically this is the inizialization of the array of frequencies  #
      # to test. As it is done here, it depends on the input data but it  #
      # can be defined in a different way. Formally Nyquist frequency is  #
      # is the largest period that can be sample and it depends on the    #
      # dataset as 1/2*dt where dt is the time between observations       #
      #                                                                   #
      fmin=1.0/(2.0*(np.max(x)-np.min(x)))       #Minimum frequency       #
      dx=x[1:]-x[:-1]                            #f_Nyquist: Max freq is  #
      # x array must be sorted                   #equal to the inverse of #
      #                                          #the median value of     #
      #                                          #the time between        #
      fmax=1.0/np.median(dx)                     #observations            #
      #                                                                   #
      n_oversampling=oversampling   #Number of frequencies per interval   #
      #                                                                   #
      #                                                                   #
      #
      fmin=0.001
      fmax=0.95
      frequency=np.arange(fmin,fmax,fmin/n_oversampling)                  #
      #####################################################################
   # File name to save plots and results                               #
   filename='GLS'                                                      #
   # Reference time (to avoid use large HJD numeric numbers)           #
   tref=np.min(x)                                                      #
   # M factor to compute analytic FAP (Zechmeister and Kurster 2009).  #
   # This factor depends on the time span of data and the array of     #
   # frequencies tested (frequency_range/frequency_step)               #
   m_factor=(np.max(frequency)-np.min(frequency))*(np.max(x)-np.min(x))
   #Print input data parameters
   if(verbose):
      print ('\nInput data')
      print ('Time domain:')
      print ('\tNumber of observations: ',len(t))
      print ('\tTime span of observations: ',np.max(t)-np.min(t))
      print ('\tTime between observations:',1.0/fmax,'(median)')
      print ('\tReference time subtracted: ',tref)
      print ('Frequency domain:')
      print ('\tf_min:',fmin,'(P:',1.0/fmin,')')
      print ('\tf_max:',fmax,'(P:',1.0/fmax,')')
      print ('\tFrequency step: ',fmin/n_oversampling,'(oversampling=',str(n_oversampling)+')')
      print ('\tNumber of frequencies to test: ',len(frequency))
      print ('\tFAP M factor: ',m_factor)
   #Compute the periodogram
   #Compute the weights of each measurement
   w=1.0/(dy*dy)
   #Initialize variables with results
   power=[]
   window=[]
   solution=[]
   dsolution=[]
   #For each frequency in the frequency array:
   for i in range(len(frequency)):
      #Compute the angular frequency array corresponding to the array of frequencies
      wfreq=2.0*np.pi*frequency[i]
      #Compute the cosx and sinx terms
      cosx=np.cos(wfreq*(x-tref))
      sinx=np.sin(wfreq*(x-tref))
      #Compute values of the matrix of equations
      sum_y=np.sum(y*w)
      sum_y2=np.sum(y*y*w)
      sum_cos=np.sum(cosx*w)
      sum_sin=np.sum(sinx*w)
      sum_ycos=np.sum(y*cosx*w)
      sum_ysin=np.sum(y*sinx*w)
      sum_cos2=np.sum(cosx*cosx*w)
      sum_sin2=np.sum(sinx*sinx*w)
      sum_cossin=np.sum(cosx*sinx*w)
      weight=np.sum(w)
      #Compute fitting values
      chi_0=weight*(sum_y2/weight-(sum_y/weight)*(sum_y/weight))
      yy=sum_y2/weight-(sum_y/weight)*(sum_y/weight)
      yc=sum_ycos/weight-(sum_y/weight)*(sum_cos/weight)
      ys=sum_ysin/weight-(sum_y/weight)*(sum_sin/weight)
      cc=sum_cos2/weight-(sum_cos/weight)*(sum_cos/weight)
      ss=sum_sin2/weight-(sum_sin/weight)*(sum_sin/weight)
      cs=sum_cossin/weight-(sum_cos/weight)*(sum_sin/weight)
      dd=cc*ss-cs*cs
      #Solve the system of equations
      c1=(yc*ss-ys*cs)/dd
      c2=(ys*cc-yc*cs)/dd
      c0=sum_y/weight-c1*sum_cos/weight-c2*sum_sin/weight
      chi2=weight*(yy-c1*yc-c2*ys)
      solution.append([chi2,c1,c2,c0])
      #Compute the uncertainties inverting the matrix of normal equations
      a11=weight
      a12=sum_cos
      a13=sum_sin
      a21=sum_cos
      a22=sum_cos2
      a23=sum_cossin
      a31=sum_sin
      a32=sum_cossin
      a33=sum_sin2
      det_a=a11*a22*a33+a12*a23*a31+a21*a32*a13-a13*a22*a31-a32*a23*a11-a12*a21*a33
      a_inv11=(a22*a33-a32*a23)/det_a
      a_inv12=-(a12*a33-a32*a13)/det_a
      a_inv13=(a12*a23-a22*a13)/det_a
      a_inv21=-(a21*a33-a31*a23)/det_a
      a_inv22=(a11*a33-a31*a13)/det_a
      a_inv23=-(a11*a23-a21*a13)/det_a
      a_inv31=(a21*a32-a31*a22)/det_a
      a_inv32=-(a11*a32-a31*a12)/det_a
      a_inv33=(a11*a22-a21*a12)/det_a
      dsolution.append([chi_0,np.sqrt(a_inv22),np.sqrt(a_inv33),np.sqrt(a_inv11)])
      #Compute the power of the periodogram and the LS window
      power.append((ss*yc*yc+cc*ys*ys-2.0*cs*yc*ys)/(yy*dd))
      window.append(np.sqrt(sum_cos*sum_cos+sum_sin*sum_sin)/weight)
   #Convert lists to arrays
   power=np.array(power)
   window=np.array(window)
   solution=np.array(solution)
   dsolution=np.array(dsolution)
   #Get the index of the highest periodogram peak
   i_best=np.argmax(power)
   #Compute residuals
   residuals=compute_residuals(x,y,solution[i_best],frequency[i_best],tref)
   #
   #Parameters used for plots and results
   peak_width=1.0/(np.max(x)-np.min(x))
   #Compute the FAP value for the best frequency
   probability_power=(1.0-power[i_best])**((float(len(x))-3.0)/2.0)
   FAP=1.0-(1.0-probability_power)**m_factor
   #Compute the power for which FAP=0.01 (1%) and FAP=0.001 (0.1%)
   FAPlim1=0.01
   prob_FAP1=1.0-np.exp(np.log(1.0-FAPlim1)/m_factor)
   power_FAP1=1.0-np.exp(np.log(prob_FAP1)/((float(len(x))-3.0)/2.0))
   FAPlim2=0.001
   prob_FAP2=1.0-np.exp(np.log(1.0-FAPlim2)/m_factor)
   power_FAP2=1.0-np.exp(np.log(prob_FAP2)/((float(len(x))-3.0)/2.0))
   power_lim=[power_FAP1,power_FAP2]
   #Show the results
   if(verbose):
      results_periodogram(x,frequency[i_best],power[i_best],solution[i_best],dsolution[i_best],residuals,m_factor)
   #
   #Plot figures the periodogram
   if(createplot):
      plot_periodogram(frequency,power,i_best,peak_width,frequency,window,power_lim,filename)
      if(errorbarplot):
         plot_phase_errorbar(x,y,dy,residuals,frequency[i_best],tref,solution[i_best],filename)
      else:
         plot_phase(x,y,residuals,frequency[i_best],tref,solution[i_best],filename)
   #
   #Save output datasets the periodogram
   if(saveoutput):
      #Save periodogram
      save_periodogram(x,frequency,power,solution,m_factor,filename)
      #Save window function
      save_window(frequency,window,filename)
      #Save the residuals and plot the phase curve
      save_residuals(x,residuals,dy,filename)
   return
#
#Function to show the results
def results_periodogram(x,freq,power,solution,dsolution,residuals,mfactor):
   '''
   Show the results of the best fitted frequency.
   Parameters: x: {array}
                  Array of independent variable
               freq: {float}
                  Frequency of the best fit
               power: {float}
                  Power of the best fit
               solution: {array}
                  Array with the best fit parameters
               dsolution: {array}
                  Array with the best fit parameters uncertainties
               residuals: {array}
                  Array of residuals of the best fit
               mfactor: {float}
                  Factor for analytic FAP
   '''
   #Compute the FAP value for the best frequency
   probability_power=(1.0-power)**((float(len(x))-3.0)/2.0)
   FAP=1.0-(1.0-probability_power)**mfactor
   #Compute the power for which FAP=0.01 (1%) and FAP=0.001 (0.1%)
   FAPlim1=0.01
   prob_FAP1=1.0-np.exp(np.log(1.0-FAPlim1)/mfactor)
   power_FAP1=1.0-np.exp(np.log(prob_FAP1)/((float(len(x))-3.0)/2.0))
   FAPlim2=0.001
   prob_FAP2=1.0-np.exp(np.log(1.0-FAPlim2)/mfactor)
   power_FAP2=1.0-np.exp(np.log(prob_FAP2)/((float(len(x))-3.0)/2.0))
   power_lim=[power_FAP1,power_FAP2]
   #Compute the theoretical width of the best frequency peak
   peak_width=1.0/(np.max(x)-np.min(x))
   #Compute the residuals
   residuals2=[]
   for i in range(len(residuals)):
      residuals2.append(residuals[i]*residuals[i])
   #Compute the amplitude
   amplitude=np.sqrt(solution[1]*solution[1]+solution[2]*solution[2])
   damplitude=np.sqrt(solution[1]*solution[1]*dsolution[1]*dsolution[1]+solution[2]*solution[2]*dsolution[2]*dsolution[2])/amplitude
   #Compute the phase angle:
   costheta=solution[1]
   sintheta=-solution[2]
   theta=np.arctan(sintheta/costheta)
   dtheta=np.sqrt(solution[2]*solution[2]*dsolution[1]*dsolution[1]+solution[1]*solution[1]*dsolution[2]*dsolution[2])/(solution[1]*solution[1]+solution[2]*solution[2])
   if(costheta<0):
      theta=theta+np.pi
   #Writing results
   print (' ')
   print ('Best fit of the periodogram:')
   print ('\tBest frequency:',freq,'width:',peak_width)
   print ('\tBest period:',1.0/freq,'width:',((1.0/freq)**2.0)*peak_width)
   print ('\tFAP:',FAP*100.0,' %')
   print ('\tFit parameters:')
   print ('\t\tRV amplitude fitted:',amplitude,'+/-',damplitude)
   print ('\t\tPhase:',theta,'+/-',dtheta,'rad','('+str(round(theta*180.0/np.pi,3))+'+/-'+str(round(dtheta*180.0/np.pi,3))+' deg)')
   if(len(solution)>=4):
      print( '\t\tRV zero value:',solution[3],'+/-',dsolution[3])
   if(len(solution)>=5):
      print ('\t\tRV linear parameter:',solution[4],'+/-',dsolution[4])
   print ('\tFit properties:')
   print ('\t\tchi^2:',solution[0],'(chi^2(mean)='+str(round(dsolution[0],3))+')')
   print ('\t\trms residuals:',np.sqrt(np.average(residuals2)))
   print ('\t\tmean residuals:',np.average(residuals))
   print ('\t\tsigma residuals:',np.std(residuals))
   print (' ')
   return

#Function to compute the residuals of the best fit
def compute_residuals(x,y,best_fit,freq_fit,tref):
   '''
   Computation of the residuals of a fit:
   Parameters: x: {float,array}
                  Array of independent variable
               y: {float,array}
                  Array of dependent variable
               best_fit: {float,array}
                  Array with the parameters of the fit
               freq_fit: {float}
                  Frequency of the fit
               tref: {float}
                  Reference time of the fit
   Output: residuals: {float,array}
                  Residuals of the fit
   '''
   residuals=[]
   for i in range(len(x)):
      x_fit=2.0*np.pi*(x[i]-tref)*freq_fit
      sinx=np.sin(x_fit)
      cosx=np.cos(x_fit)
      if(len(best_fit)==3):
#        Residuals of the LS fit
         y_fit=best_fit[1]*cosx+best_fit[2]*sinx
      elif(len(best_fit)>=4):
#        Residuals of the sinusoidal component of the fits
         y_fit=best_fit[3]+best_fit[1]*cosx+best_fit[2]*sinx
      residuals.append(y[i]-y_fit)
   return residuals
#
#Function to plot the periodogram
def plot_periodogram(x,y,i_fit,peak_width,xwindow,window,ylim,plotname):
   '''
   Plot the computed periodogram
   Parameters: x: {float, array}
                  Array of frequencies tested
               y: {float, array}
                  Array of power for each frequency
               i_fit: {int}
                  Index of the best fit
               x_window: {float, array}
                  Array with frequencies of the window function
               window: {float, array}
                  Array with the value of the window function
               ylim: {float, array}
                  Array with FAP limits
               plotname: {float, array}
                  Name of the plot (saved as pdf)
   '''
   xmin=np.min([np.min(x),np.min(xwindow)])
   xmax=np.max([np.max(x),np.max(xwindow)])
   ymin=0
   ymax=np.max(y)*1.10
#  Figure with the periodogram and window function
   fig=plt.figure(figsize=(8, 6))
   t=fig.text(0.5,0.9,f'Best period: {1.0/x[i_fit]:.2g} days',horizontalalignment='right',fontsize=11)
#  Periodogram
   panel1=fig.add_axes([0.10,0.10,0.85,0.50])
   panel1.minorticks_on()
   f1=panel1.plot(x, y, color='cyan', linestyle='-')
   fl1=panel1.plot([xmin,xmax],[ylim[0],ylim[0]], color='cyan', linestyle='-')
   fl2=panel1.plot([xmin,xmax],[ylim[1],ylim[1]], color='cyan', linestyle='-')
   panel1.axvline(x=x[i_fit],linestyle='-.',color='g',lw=2)
   panel1.set_xlabel('Frequency', fontsize=12)
   panel1.set_ylabel('Power', fontsize=12)
   panel1.set_xlim(xmin,xmax)
   
   panel1.set_ylim(ymin,ymax)
   plt.tight_layout()
   #plt.savefig(plotname+'s53_periodogram.png',orientation='landscape',dpi=400)
   
#  Window function
   winmin=0
   winmax=1.10*np.max(window)
   panel2=fig.add_axes([0.10,0.60,0.85,0.35])
   panel2.minorticks_on()
   f2=panel2.plot(xwindow,window, color='cyan', linestyle='-')
   panel2.axvline(x=x[i_fit],linestyle='-.',color='r',lw=2)
   panel2.set_ylabel('Window function', fontsize=12)
   panel2.set_xlim(xmin,xmax)
   panel2.set_xticklabels([])
   panel2.set_ylim(winmin,winmax)
   #panel2.legend(loc='center')
   plt.savefig(plotname+'s53_periodogram.png',orientation='landscape',dpi=400)
   plt.show()
#
#  Figure with a zoom around the best period found
   fig=plt.figure()
   t=fig.text(0.5,0.96,'Best period: '+str(1.0/x[i_fit]),horizontalalignment='center')
#  Periodogram
   panel1=fig.add_axes([0.10,0.10,0.85,0.50])
   panel1.minorticks_on()
   f1=panel1.plot(x,y,'k-')
   fl1=panel1.plot([xmin,xmax],[ylim[0],ylim[0]],'k--')
   fl2=panel1.plot([xmin,xmax],[ylim[1],ylim[1]],'k-.')
   panel1.axvline(x=x[i_fit],linestyle='-.',color='r')
   panel1.set_xlabel('Frequency')
   panel1.set_ylabel('Power')
   panel1.set_xlim(x[i_fit]-10.0*peak_width,x[i_fit]+10.0*peak_width)
   panel1.set_ylim(ymin,ymax)
#  Window function
   panel2=fig.add_axes([0.10,0.60,0.85,0.35])
   panel2.minorticks_on()
   f2=panel2.plot(xwindow,window,'k-')
   f2p=panel2.plot(xwindow+x[i_fit],window,'k-.')
   f2m=panel2.plot(-xwindow+x[i_fit],window,'k-.')
   panel2.axvline(x=x[i_fit],linestyle='-.',color='r')
   panel2.set_ylabel('Window func.')
   panel2.set_xlim(x[i_fit]-10.0*peak_width,x[i_fit]+10.0*peak_width)
   panel2.set_xticklabels([])
   panel2.set_ylim(winmin,winmax)
   plt.savefig(plotname+'_peak.png',orientation='landscape')
   print ('Periodogram plot saved as',plotname+'_periodogram.png and',plotname+'_peak.png')
   plt.show()
#
#Function to plot the phase curve
def plot_phase(x,yreal,residu,freq,tref,best_fit,filename):   
   '''
   Plot the phase curve of the best fit
   Parameters: x: {float, array}
                  Array of independent variables
               yreal: {float, array}
                  Array of dependent variables
               freq: {float}
                  Frequency of the best fit
               best_fit: {float,array}
                  Array with the parameters of the best fit
               filename: {str}
                  Root name of the file with the plot saved as pdf
   '''
   #Compute phase of the curve
   phase=[]
   y1=[]
   y2=[]
   for k in range(len(x)):
      ph=(x[k]-tref)*freq-int((x[k]-tref)*freq)
      if(ph>=1.0):
         ph=ph-1.0
      if(ph<0.0):
         ph=ph+1.0
      phase.append(ph)
      y1.append(yreal[k])
      y2.append(residu[k])
      if(ph>0.8):
         phase.append(ph-1.0)
         y1.append(yreal[k])
         y2.append(residu[k])
      if(ph<0.2):
         phase.append(ph+1.0)
         y1.append(yreal[k])
         y2.append(residu[k])
   y1=np.array(y1)
   y2=np.array(y2)
   #Binning
   if(len(x)>1000):
      dy=np.array([1.0]*len(phase))
      x1bin,y1bin,dy1bin=phase_bin(phase,y1,dy,-0.2,1.2,0.01)
      x2bin,y2bin,dy2bin=phase_bin(phase,y2,dy,-0.2,1.2,0.01)
   #Compute the theoretical curve and the best fit curve
   ph=np.arange(-0.2,1.201,0.01)
   fit=[]
   for k in range(len(ph)):
      sinx=np.sin(2.0*np.pi*ph[k])
      cosx=np.cos(2.0*np.pi*ph[k])
      if(len(best_fit)==3):
#        Residuals of the LS fit
         y_fit=best_fit[1]*cosx+best_fit[2]*sinx
      elif(len(best_fit)>=4):
#        Residuals of the sinusoidal component of the fits
         y_fit=best_fit[3]+best_fit[1]*cosx+best_fit[2]*sinx
      fit.append(y_fit)
   #Plot the figure
   rv_amp_fit=np.sqrt(best_fit[1]*best_fit[1]+best_fit[2]*best_fit[2])
   period_fit=1.0/freq
   label_fit='P='+str(round(period_fit,2))+' K='+str(round(rv_amp_fit,4))
   #Phase range
   xmin=-0.2 #Minimum phase
   xmax=1.2  #Maximum phase
   ticks=np.arange(xmin,xmax,0.1)
   mticks=np.arange(xmin,xmax,0.02)
   #Initialize figure
   fig=plt.figure()
   t=fig.text(0.5,0.96,'Best period: '+str(1.0/freq),horizontalalignment='center')
   #Panel 1: radial velocity curve
   #Y-scale
   ymin=np.min([np.min(y1),np.min(fit)])
   ymax=np.max([np.max(y1),np.max(fit)])
   dy=(ymax-ymin)/10.0
   ymin=ymin-dy
   ymax=ymax+dy
   #Plot of data
   panel1=fig.add_axes([0.10,0.10,0.85,0.65])
   panel1.minorticks_on()
   f11=panel1.scatter(phase,y1,c='k')
   #Plot binned curve
   if(len(x)>1000):
      f1b=panel1.scatter(x1bin,y1bin,c='c',zorder=2)
   f12=panel1.plot(ph,fit,'r-',label=label_fit)
   panel1.set_xlim(xmin,xmax)           #X-axis limits
   panel1.set_xticks(ticks)             #X-axis major ticks (only with linear labels)
   panel1.set_xticks(mticks,minor=True) #X-axis minor ticks (only with linear labels)
   panel1.set_xlabel('Phase')
   panel1.set_ylim(ymin,ymax)           #Y-axis limits reversed -> panel1.set_ylim(ymax,ymin)
   panel1.set_ylabel('Variable')
   panel1.legend(loc=0,frameon=False,shadow=None)
   #Panel 2: O-C values
   #Y scale: definned to be symmetric
   ymax=np.max([np.abs(np.min(y2)),np.max(y2)])
   ymin=-ymax
   dy=(ymax-ymin)/10.0
   ymin=ymin-dy
   ymax=ymax+dy
   #Plotting data
   panel2=fig.add_axes([0.10,0.75,0.85,0.2])
   panel2.locator_params(axis = 'y', nbins = 4)
   panel2.minorticks_on()
   f21=panel2.scatter(phase,y2,c='k')
   if(len(x)>1000):
      f21b=panel2.scatter(x2bin,y2bin,c='c',zorder=1)
   f22=panel2.plot([xmin,xmax],[0,0],'r-')
   panel2.set_xlim(xmin,xmax)           #X-axis limits
   panel2.set_xticks(ticks)             #X-axis major ticks (only with linear labels)
   panel2.set_xticks(mticks,minor=True) #X-axis minor ticks (only with linear labels)
   panel2.set_xticklabels([])           #Delete labels
   panel2.set_ylim(ymin,ymax)
   panel2.set_ylabel('O-C')
   #Save figure to an .eps file
   print ('Phase plot saved as',filename+'_phase.ng')
   plt.savefig(filename+'_phase.png',orientation='landscape')
   plt.show(block=True)
#  End plotting function
#
def phase_bin(x,y,dy,start,end,size):
   xbins=np.arange(start,end+size,size)
   xbin=[]
   ybin=[]
   dybin=[]
   for i in range(len(xbins)-1):
      value_x=[]
      value_y=[]
      value_dy=[]
      for j in range(len(x)):
         if(xbins[i]<=x[j]<xbins[i+1]):
            value_x.append(x[j])
            value_y.append(y[j])
            value_dy.append(dy[j])
      value_x=np.array(value_x)
      value_y=np.array(value_y)
      value_dy=np.array(value_dy)
      if(len(value_x)>0):
         xbin.append(np.average(value_x,weights=1.0/(value_dy)))
         ybin.append(np.average(value_y,weights=1.0/(value_dy)))
         dybin.append(np.std(value_y))
   return np.array(xbin),np.array(ybin),np.array(dybin)
#phase_bin(phase,flux,eflux,phase_start,phase_end,phase_bin_size)
#Function to plot the phase curve with errorbars
def plot_phase_errorbar(x,yreal,dyerror,residu,freq,tref,best_fit,filename):   
   '''
   Plot the phase curve of the best fit
   Parameters: x: {float, array}
                  Array of independent variables
               yreal: {float, array}
                  Array of dependent variables
               dy: {float, array}
                  Array of dependent variables uncertainties
               freq: {float}
                  Frequency of the best fit
               best_fit: {float,array}
                  Array with the parameters of the best fit
               filename: {str}
                  Root name of the file with the plot saved as pdf
   '''
   #Compute phase of the curve
   phase=[]
   y1=[]
   y2=[]
   dy1=[]
   for k in range(len(x)):
      ph=(x[k]-tref)*freq-int((x[k]-tref)*freq)
      if(ph>=1.0):
         ph=ph-1.0
      if(ph<0.0):
         ph=ph+1.0
      phase.append(ph)
      y1.append(yreal[k])
      y2.append(residu[k])
      dy1.append(dyerror[k])
      if(ph>0.8):
         phase.append(ph-1.0)
         y1.append(yreal[k])
         y2.append(residu[k])
         dy1.append(dyerror[k])
      if(ph<0.2):
         phase.append(ph+1.0)
         y1.append(yreal[k])
         y2.append(residu[k])
         dy1.append(dyerror[k])
   y1=np.array(y1)
   y2=np.array(y2)
   dy1=np.array(dy1)
   #Binning
   if(len(x)>1000):
      x1bin,y1bin,dy1bin=phase_bin(phase,y1,dy1,-0.2,1.2,0.01)
      x2bin,y2bin,dy2bin=phase_bin(phase,y2,dy1,-0.2,1.2,0.01)
   #Compute the theoretical curve and the best fit curve
   ph=np.arange(-0.2,1.201,0.01)
   fit=[]
   for k in range(len(ph)):
      sinx=np.sin(2.0*np.pi*ph[k])
      cosx=np.cos(2.0*np.pi*ph[k])
      if(len(best_fit)==3):
#        Residuals of the LS fit
         y_fit=best_fit[1]*cosx+best_fit[2]*sinx
      elif(len(best_fit)>=4):
#        Residuals of the sinusoidal component of the fits
         y_fit=best_fit[3]+best_fit[1]*cosx+best_fit[2]*sinx
      fit.append(y_fit)
   #Plot the figure
   rv_amp_fit=np.sqrt(best_fit[1]*best_fit[1]+best_fit[2]*best_fit[2])
   period_fit=1.0/freq
   label_fit='P='+str(round(period_fit,2))+' K='+str(round(rv_amp_fit,4))
   #Phase range
   xmin=-0.2 #Minimum phase
   xmax=1.2  #Maximum phase
   ticks=np.arange(xmin,xmax,0.1)
   mticks=np.arange(xmin,xmax,0.02)
   #Initialize figure
   fig=plt.figure()
   t=fig.text(0.5,0.96,'Best period: '+str(1.0/freq),horizontalalignment='center')
   #Panel 1: radial velocity curve
   #Y-scale
   ymin=np.min([np.min(y1),np.min(fit)])
   ymax=np.max([np.max(y1),np.max(fit)])
   dy=(ymax-ymin)/10.0
   ymin=ymin-dy
   ymax=ymax+dy
   #Plot of data
   panel1=fig.add_axes([0.10,0.10,0.85,0.65])
   panel1.minorticks_on()
   f11=panel1.errorbar(phase,y1,yerr=dy1,fmt='o',c='k',capthick=2)
   #Plot binned curve
   if(len(x)>1000):
      f11b=panel1.errorbar(x1bin,y1bin,yerr=dy1bin,fmt='o',c='c',ecolor='c',markeredgecolor='none',capsize=2,zorder=1)
   f12=panel1.plot(ph,fit,'r-',label=label_fit)
   panel1.set_xlim(xmin,xmax)           #X-axis limits
   panel1.set_xticks(ticks)             #X-axis major ticks (only with linear labels)
   panel1.set_xticks(mticks,minor=True) #X-axis minor ticks (only with linear labels)
   panel1.set_xlabel('Phase')
   panel1.set_ylim(ymin,ymax)           #Y-axis limits reversed -> panel1.set_ylim(ymax,ymin)
   panel1.set_ylabel('Variable')
   panel1.legend(loc=0,frameon=False,shadow=None)
   #Panel 2: O-C values
   #Y scale: definned to be symmetric
   ymax=np.max([np.abs(np.min(y2)),np.max(y2)])
   ymin=-ymax
   dy=(ymax-ymin)/10.0
   ymin=ymin-dy
   ymax=ymax+dy
   #Plotting data
   panel2=fig.add_axes([0.10,0.75,0.85,0.2])
   panel2.locator_params(axis = 'y', nbins = 4)
   panel2.minorticks_on()
   f21=panel2.scatter(phase,y2,c='k')
   #Plot binned residuals
   if(len(x)>1000):
      f21b=panel2.errorbar(x2bin,y2bin,yerr=dy2bin,fmt='o',c='c',ecolor='c',markeredgecolor='none',capsize=2,zorder=1)
   f22=panel2.plot([xmin,xmax],[0,0],'r-')
   panel2.set_xlim(xmin,xmax)           #X-axis limits
   panel2.set_xticks(ticks)             #X-axis major ticks (only with linear labels)
   panel2.set_xticks(mticks,minor=True) #X-axis minor ticks (only with linear labels)
   panel2.set_xticklabels([])           #Delete labels
   panel2.set_ylim(ymin,ymax)
   panel2.set_ylabel('O-C')
   #Save figure to an .eps file
   print ('Phase plot saved as',filename+'_phase.png')
   plt.savefig(filename+'_phase.png',orientation='landscape',dpi=400)
   plt.show(block=True)
#  End plotting function
#
# Function to save the periodogram results
def save_periodogram(x,freq,power,solution,mfactor,filename):
   '''
   Save the periodogram results
   Parameters: x: {float, array}
                  Array of independent variables
               freq: {float,array}
                  Array of frequencies
               power: {float, array}
                  Array of periodogram power
               solution: {float,array}
                  Array with the parameters of the best fit for each frequency
               filename: {str}
                  Root name of the file with results stored as plain text
   '''
   #Save periodogram results
   file_out=open(filename+'_periodogram.dat','w')
   if(len(solution[0])==3):
      #LS data only
      file_out.write('#Frequency   Power   FAP(%)   chi^2   RV_ampl  ampl_cos   ampl_sin\n')
   else:
      #GLS additional parameters are also stored
      file_out.write('#Frequency   Power   FAP(%)   chi^2   RV_ampl  ampl_cos   ampl_sin   GLS_terms\n')
   #Save data for each frequency
   for i in range(len(freq)):
      probability_power=(1.0-power[i])**((float(len(x))-3.0)/2.0)
      FAP=1.0-(1.0-probability_power)**mfactor
      rv_fit=np.sqrt(solution[i][1]*solution[i][1]+solution[i][2]*solution[i][2])
      textline=str(freq[i])+' '+str(power[i])+' '+str(FAP*100.0)+' '+str(solution[i][0])+' '+str(rv_fit)+' '+str(solution[i][1])+' '+str(solution[i][2])
      #Add parameters of the generalised LS periodogram
      if(len(solution[i])>=4):
         for param in solution[i][3:]:
            textline=textline+' '+str(param)
      file_out.write(textline+'\n')
   file_out.close()
   print ('Periodogram data saved as',filename+'_periodogram.dat')
   return
#
# Function to save the window function
def save_window(freq,window,filename):
   '''
   Save the LS window function
   Parameters: freq: {float,array}
                  Array of frequencies
               window: {float, array}
                  Array of values of the window functions
               filename: {str}
                  Root name of the file with window function stored as plain text
   '''
   #Save periodogram results
   file_out=open(filename+'_window.dat','w')
   file_out.write('#Frequency   window  \n')
   for i in range(len(freq)):
      file_out.write(str(freq[i])+' '+str(window[i])+'\n')
   file_out.close()
   print ('Window function data saved as',filename+'_window.dat')
   return
#
# Function to compute the residuals
def save_residuals(x,residuals,dy,filename):
   '''
   Save the residuals of the best fit so they can be analyzed again
   Parameters: x: {float, array}
                  Array of independent variables
               residuals: {float,array}
                  Array of residuals of the best fit
               dy: {float, array}
                  Array of uncertainties
               filename: {str}
                  Root name of the file with residuals stored as plain text
   '''
   #Save periodogram results
   file_out=open(filename+'_residuals.dat','w')
   file_out.write('#Time  residual  error\n')
   for i in range(len(x)):
      file_out.write(str(x[i])+' '+str(residuals[i])+' '+str(dy[i])+'\n')
   file_out.close()
   print ('Residuals saved as',filename+'_residuals.dat')
   return
#
#Main code
#Load input data
filename=sys.argv[1]
t,y,dy=read_table(filename)
#Call periodogram function
periodGLS(t,y,dy,verbose=True)


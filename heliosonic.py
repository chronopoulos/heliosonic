"""
heliosonic.py
Tools for the sonification of helioseismic data.
Chris Chronopoulos, 01-22-2013
chronopoulos.chris@gmail.com
"""

import os, pylab, pyfits
import numpy as np
import scipy.io.wavfile as wav

def MeshGrid(x,y):
   xx = np.outer(x,np.ones(len(y)))
   yy = np.outer(np.ones(len(x)),y)
   return xx,yy

def SpatialApod(ny,nx):
   mask = np.zeros((ny,nx),dtype='float32')
   if nx <= 16:
      rInner = 0.8750
   else:
      rInner = 0.9375
   rOuter = 1.0
   xx,yy = MeshGrid(np.linspace(-1.,1.,num=nx), np.linspace(-1.,1.,num=ny))
   rr = np.sqrt(xx**2+yy**2)
   indInner = np.where( rr <= rInner )
   indBetween = np.where( (rr > rInner) * (rr < rOuter) )
   indOuter = np.where( rr >= rOuter )
   mask[indInner] = 1.
   rs = (rr[indBetween]-rInner)/(rOuter-rInner)
   mask[indBetween] = (1.0 - rs**2)**2
   mask[indOuter] = 0.
   return mask

def FullApod(nt,ny,nx):
   """
   Thanks to Ben Greer for this function.
   """
   t = np.linspace(-1.,1.,num=nt)
   t = t.astype('float32') # that's annoying
   tOuter = 1.
   tInner = 0.96875
   indInner = np.where( np.abs(t) <= tInner )
   indBetween = np.where( (np.abs(t) > tInner) * (np.abs(t)<tOuter) )
   indOuter = np.where( np.abs(t) >= tOuter )
   tMask = np.zeros(nt, dtype='float32')
   tMask[indInner] = 1.0
   rs = (abs(t[indBetween])-tInner)/(tOuter-tInner)
   tMask[indBetween] = (1.0 - rs**2)**2
   tMask[indOuter] = 0.
   rMask=SpatialApod(ny,nx)
   mask=np.zeros((nt,ny,nx),dtype='float32')
   for i in range(nt):
      mask[i,:,:] = tMask[i]*rMask
   return mask

###

class DataCube():
   """
   Encapsulates the raw data cube (images over time) and its manipulations.
   Takes an SDO/HMI FITS file as an argument.
   Example Instantiation:
   >>> fft = DataCube('full_fft.npy')
   """

   def __init__(self,filename):
      print 'Reading in file '+filename
      fitsFile = pyfits.open(filename)
      self.data = fitsFile[0].data
      self.nt, self.ny, self.nx = pylab.shape(self.data)
      # Remove the DC offset due to the spacecraft's motion
      for i in range(self.nt):
         self.data[i,:,:] -= np.mean(self.data[i,:,:])

   def SpatialFFT(self,*arg):
      """
      Compute the spatial FFT of each timestep of the data set.
      Uses a tapered-tophat 2D apodization mask.
      If called without an argument, return the result.
      If called with a string argument ending in .npy, then save to that file in NumPy format.
      """
      spatialMask=SpatialApod(self.ny,self.nx)
      data_fft = np.empty( (self.nt,self.ny,self.nx), dtype='complex64' )
      for i in range(self.nt):
         print 'Timestep '+str(i)+' of '+str(self.nt)
         data_fft[i,:,:] = np.fft.fft2(spatialMask*self.data[i,:,:])
      if len(arg)==0:
         return data_fft[:,0:self.ny//2+1,0:self.nx//2+1]
      elif os.path.splitext(arg[0])[1]=='.npy':
         print 'Writing to file '+arg[0]
         np.save(arg[0], data_fft[:,0:self.ny//2+1,0:self.nx//2+1])
         print 'Done writing.'
      else:
         raise Exception('Invalid filename argument')

   def FullFFT(self,*arg):
      """
      Compute the full FFT of the data set - in lat, lon, and time.
      Uses a tapered-tophat apodization mask (see function fullApod).
      If called without an argument, return the result.
      If called with a string argument ending in .npy, then save to that file in NumPy format.
      """
      mask = FullApod(self.nt, self.ny, self.nx)
      data_fft = np.empty((self.nt,self.ny,self.nx),dtype='complex64')
      # confused: rfftn only treats the last axis as real??
      data_fft[:] = np.fft.fftn(mask*self.data)
      if len(arg)==0:
         return data_fft[0:self.nt//2+1, 0:self.ny//2+1, 0:self.nx//2+1]
      elif os.path.splitext(arg[0])[1]=='.npy':
         print 'Writing to file '+arg[0]
         np.save(arg[0], data_fft[0:self.nt//2+1, 0:self.ny//2+1, 0:self.nx//2+1])
         print 'Done writing.'
      else:
         raise Exception('Invalid filename argument')

   def MakeMovie(self,imgDir=''):
      """
      Write out 2D colormaps of each timestep, to be used as frames in a movie.
      Optional argument imgDir specifies the directory in which to place the files.
      """
      for i in range(self.nt):
         imgPlot = pylab.imshow(self.data[i,:,:])
         imgPlot.set_cmap('hot')
         numstr = str(i)
         imgfilename = imgDir + numstr.zfill( int( np.log10(self.nt) ) ) + '.png'
         print 'Saving ' + imgfilename
         pylab.savefig(imgfilename)

class FullFFT():
   """
   Encapsulates the full, 3D Fourier Transform of the data, and its manipulations.
   Takes a NumPy .npy file containing the full FFT data as the argument.
   Example Instantiation:
   >>> fft = FullFFT('fullFFT.npy')
   """

   def __init__(self,filename):
      print 'Reading in file '+filename
      self.data = np.load(filename)
      self.nomega, self.nky, self.nkx = np.shape(self.data)
      self.omega = np.arange(self.nomega)
      self.kx, self.ky = MeshGrid(np.arange(self.nkx),np.arange(self.nky))
      self.kr = np.sqrt(self.kx**2 + self.ky**2)
      #self.ktheta = np.arctan(self.ky/self.kx)

   def CutKy(self, iky):
         spectrum = np.absolute(self.data[:,iky,:])**2
         imgPlot = pylab.imshow(np.log10(spectrum))
         imgPlot.set_cmap('cool_r')
         pylab.show()


class SpatialFFT():
   """
   Encapsulates the spatial Fourier Transform of the data, and its manipulations.
   Takes a NumPy .npy file containing the spatial FFT data as the argument.
   Example Instantiation:
   >>> fft = SpatialFFT('spatialFFT.npy')
   """

   def __init__(self,filename):
      print 'Reading in file '+filename
      self.data = np.load(filename)
      self.nt, self.nky, self.nkx = np.shape(self.data)
      self.kx, self.ky = MeshGrid(np.arange(self.nkx),np.arange(self.nky))
      self.kr = np.sqrt(self.kx**2 + self.ky**2)
      self.ktheta = np.arctan2(self.ky,self.kx)

   def GetPointSignal(self,ikx,iky):
      signal = self.data[:,iky,ikx]
      return signal

   def GetRingSignal(self,ikr,dikr=2):
      """
      Display the signal and spectrum of all modes in a ring of radius kr.
      Arguments: index ikr, ring width dikr (optional)
      """
      indRing = np.where( (self.kr >= ikr-dikr)*(self.kr <= ikr+dikr) )
      signal=np.zeros(self.nt)
      for i in range(np.shape(indRing)[1]):
         signal += self.GetPointSignal(indRing[0][i], indRing[1][i])
      return signal

   def GetRingAmpSignal(self,ikr,dikr=2):
      """
      Display the signal and spectrum of all modes in a ring of radius kr.
      Arguments: index ikr, ring width dikr (optional)
      """
      indRing = np.where( (self.kr >= ikr-dikr)*(self.kr <= ikr+dikr) )
      signal=np.zeros(self.nt)
      for i in range(np.shape(indRing)[1]):
         signal += np.absolute( self.GetPointSignal(indRing[0][i], indRing[1][i]) )
      return signal

   def GetNormalizedSignal(self,signal):
      """
      Given a signal, remove its DC offset,
      and normalize to a peak value of 1.
      """
      signal = np.real(signal)
      signal = signal - np.mean(signal)
      signal = signal / max(abs(signal))
      return signal

   def PlotSignalSpectrum(self,signal):
      """
      Display a plot of the signal and its spectrum.
      Argument: 1D signal array
      """
      fft = np.fft.fft(signal)
      spectrum = np.absolute(fft[0:self.nt/2-1])**2
      pylab.figure(figsize=(14,6))
      pylab.subplot(121, axisbg='k')
      pylab.plot(signal,'y-')
      pylab.xlim([0,2047])
      pylab.title('Signal')
      pylab.subplot(122, axisbg='k')
      #pylab.plot(spectrum,'g-')
      pylab.plot(spectrum,'g-')
      pylab.xlim([100,500])
      pylab.title('Spectrum')
      pylab.show()

   def DisplayPoint(self,ikx,iky):
      """
      Display the signal and spectrum of a single mode (kx,ky).
      Arguments: indeces ikx, iky
      """
      signal=self.GetPointSignal(ikx,iky)
      self.PlotSignalSpectrum(signal)

   def DisplayRing(self, ikr, dikr=2):
      signal = self.GetNormalizedSignal( self.GetRingSignal(ikr,dikr) )
      self.PlotSignalSpectrum(signal)

   def ExploreWAV(self):
      satisfied=False
      while not satisfied:
         ikr = input('Enter ikr: ')
         signal = self.GetNormalizedSignal(self.GetRingSignal(ikr))
         self.PlotSignalSpectrum(signal)
         response = raw_input('Satisfied? y/n: ')
         satisfied = (response=='y')
      filename = 'ringsample_'+str(ikr)+'.wav'
      response = raw_input('Write to '+filename+'? y/n: ')
      if response=='y':
         signal_int16 = np.int16(signal*2**15*0.99)
         wav.write(filename,44100,signal_int16)


"""
heliosonic.py
Tools for the sonficiation of helioseismic data.
Chris Chronopoulos, 01-22-2013
chronopoulos.chris@gmail.com
"""

import os, pylab, pyfits
import numpy as np
import scikits.audiolab as audiolab

def meshGrid(x,y):
   xx = np.outer(x,np.ones(len(y)))
   yy = np.outer(np.ones(len(x)),y)
   return xx,yy

def apodMask(nx,ny):
   mask = np.zeros((nx,ny),dtype='float32')
   if nx <= 16:
      rInner = 0.8750
   else:
      rInner = 0.9375
   rOuter = 1.0
   xx,yy = meshGrid(np.linspace(-1.,1.,num=nx), np.linspace(-1.,1.,num=ny))
   rr = np.sqrt(xx**2+yy**2)
   indInner = np.where( rr <= rInner )
   indBetween = np.where( (rr > rInner) * (rr < rOuter) )
   indOuter = np.where( rr >= rOuter )
   mask[indInner] = 1.
   rs = (rr[indBetween]-rInner)/(rOuter-rInner)
   mask[indBetween] = (1.0 - rs**2)**2
   mask[indOuter] = 0.
   return mask

###

class imgObject():
   """
   Encapsulates the raw data (image, as opposed to its transform) and its manipulations.
   Takes an SDO/HMI FITS file as an argument.
   """

   def __init__(self,filename):
      print 'Reading in file '+filename
      fitsFile = pyfits.open(filename)
      self.data = fitsFile[0].data
      self.nt, self.ny, self.nx = pylab.shape(self.data)

   def computeFFT(self,*arg):
      """
      Compute the spatial FFT of each timestep of the data set.
      Uses a tapered-tophat 2D apodization mask.
      If called without an argument, return the result.
      If called with a string argument ending in .npy, then save to that file in NumPy format.
      """
      spatialMask=apodMask(self.ny,self.nx)
      data_fft = pylab.zeros( (self.nt,self.ny,self.nx), dtype='complex64' )
      for i in range(self.nt):
         print 'Timestep '+str(i)+' of '+str(self.nt)
         data_fft[i,:,:] = np.fft.fft2(spatialMask*self.data[i,:,:])
      if len(arg)==0:
         return data_fft
      elif os.path.splitext(arg)[1]=='.py':
         np.save(arg, data_fft[:,0:nx/2-1,0:nx/2-1])
      else:
         raise Exception('Invalid filename argument')

   def makeMovie(self,imgDir=''):
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
      

class fftObject():
   """
   Encapsulates the spatial Fourier Transform of the data, and its manipulations.
   Takes a NumPy .npy file containing the FFT data as the argument.
   """

   def __init__(self,filename):
      print 'Reading in file '+filename
      self.data = np.load(filename)
      self.nt, self.nky, self.nkx = np.shape(self.data)
      self.kx, self.ky = meshGrid(np.arange(self.nkx),np.arange(self.nky))
      self.kr = np.sqrt(self.kx**2 + self.ky**2)
      self.ktheta = np.arctan(self.ky/self.kx)

   def plotSignalSpectrum(self,signal):
      """
      Display a plot of the signal and its spectrum.
      Argument: 1D signal array
      """
      fft = np.fft.fft(signal)
      spectrum = np.absolute(fft[0:self.nt/2-1])**2
      pylab.figure(figsize=(14,6))
      pylab.subplot(121, axisbg='k')
      pylab.plot(signal,'y-')
      pylab.title('Signal')
      pylab.subplot(122, axisbg='k')
      pylab.semilogy(spectrum,'g-')
      pylab.title('Spectrum')
      pylab.show()

   def displayMode(self,ikx,iky):
      """
      Display the signal and spectrum of a single mode (kx,ky).
      Arguments: indeces ikx, iky
      """
      signal=self.data[:,iky,ikx]
      self.plotSignalSpectrum(signal)

   def displayRing(self,ikr,dikr=2):
      """
      Display the signal and spectrum of all modes in a ring of radius kr.
      Arguments: index ikr, ring width dikr (optional)
      """
      indRing = np.where( (self.kr >= ikr-dikr)*(self.kr <= ikr+dikr) )
      signal=np.zeros(self.nt)
      for i in range(np.shape(indRing)[1]):
         signal += self.data[ :, indRing[1][i], indRing[0][i] ]
      self.plotSignalSpectrum(signal)



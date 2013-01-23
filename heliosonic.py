import pylab
import pyfits
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


class imgObject():

   def __init__(self,filename):
      print 'Reading in file '+filename
      fitsFile = pyfits.open(filename)
      self.data = fitsFile[0].data
      self.nt, self.ny, self.nx = pylab.shape(self.data)

   def computeFFT(self,filename):
      spatialMask=apodMask(self.ny,self.nx)
      data_fft = pylab.zeros( (self.nt,self.ny,self.nx), dtype='complex64' )
      for i in range(self.nt):
         print 'Timestep '+str(i)+' of '+str(self.nt)
         data_fft[i,:,:] = np.fft.fft2(spatialMask*self.data[i,:,:])
      np.save(filename, data_fft[:,0:nx/2-1,0:nx/2-1])

   def makeMovie(self,filename,imgDir='/home/chrono/data/blork/helio/img/'):
      for i in range(self.nt):
         imgPlot = pylab.imshow(self.data[i,:,:])
         imgPlot.set_cmap('hot')
         numstr = str(i)
         imgfilename = imgDir + numstr.zfill( int( np.log10(self.nt) ) ) + '.png'
         print 'Saving ' + imgfilename
         pylab.savefig(imgfilename)
      

class fftObject():

   def __init__(self,filename):
      print 'Reading in file '+filename
      self.data = np.load(filename)
      self.nt, self.nky, self.nkx = np.shape(self.data)
      self.kx, self.ky = meshGrid(np.arange(self.nkx),np.arange(self.nky))
      self.kr = np.sqrt(self.kx**2 + self.ky**2)
      self.ktheta = np.arctan(self.ky/self.kx)

   def plotSignalSpectrum(self,signal,spectrum):
      pylab.figure(figsize=(14,6))
      pylab.subplot(121, axisbg='k')
      pylab.plot(signal,'y-')
      pylab.title('Signal')
      pylab.subplot(122, axisbg='k')
      pylab.semilogy(spectrum,'g-')
      pylab.title('Spectrum')
      pylab.show()

   def displaySingleMode(self,ikx,iky):
      signal=self.data[:,iky,ikx]
      fft = np.fft.fft(signal)
      spectrum = np.absolute(fft[0:self.nt/2-1])**2
      self.plotSignalSpectrum(signal,spectrum)

   def displaySingleKr(self,ikr,dikr=2):
      indRing = np.where( (self.kr >= ikr-dikr)*(self.kr <= ikr+dikr) )
      signal=np.zeros(self.nt)
      for i in range(np.shape(indRing)[1]):
         signal += self.data[ :, indRing[1][i], indRing[0][i] ]
      fft = np.fft.fft(signal)
      spectrum = np.absolute(fft[0:self.nt/2-1])**2
      self.plotSignalSpectrum(signal,spectrum)


# Main Program
if __name__=='__main__':
   tile16_img=imgObject('/home/chrono/data/blork/helio/tile16_+285.0000_+0.0000.fits')
   (nt,ny,nx) = pylab.shape(data_img)
   spatialMask=apodMask(ny,nx)
   data_fft = pylab.zeros( (nt,ny,nx), dtype='complex64' )
   for i in range(nt):
      print 'Timestep '+str(i)+' of '+str(nt)
      data_fft[i,:,:] = pylab.fft2(spatialMask*data_img[i,:,:])
   np.save('/home/chrono/data/blork/helio/data_fft.npy', data_fft[:,0:nx/2-1,0:nx/2-1])

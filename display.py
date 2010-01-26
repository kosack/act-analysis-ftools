import pyfits
from kapteyn import maputils
from matplotlib import pyplot as plt
import sys

fname = sys.argv[1]
print "*** ",fname

im = pyfits.open( fname )[0]

f = maputils.FITSimage( externalheader=im.header, externaldata=im.data)

fig = plt.figure()
frame = fig.add_subplot(1,1,1)

img = f.Annotatedimage(frame)
colorbar = img.Colorbar()

img.Image()
img.Graticule()

# galactic overlay:
#gr2 = img.Graticule(skyout="galactic", visible=True)
#gr2.setp_plotaxis(("top","right"), mode=maputils.native, color='g', visible=False)
#gr2.setp_tick(wcsaxis=(0,1), color='g')
#gr2.setp_gratline(wcsaxis=(0,1), color='g')

img.plot()
img.interact_imagecolors()
img.interact_toolbarinfo()
img.interact_writepos()

plt.title( fname )
plt.show()



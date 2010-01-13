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
img.plot()
img.interact_imagecolors()
img.interact_toolbarinfo()



plt.title( fname )
plt.show()



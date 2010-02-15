import pyfits
from kapteyn import maputils
from matplotlib import pyplot as plt
import sys
import math

nfigs = int(len(sys.argv)-1)
fig = plt.figure()
ifig = 1

nx = 1
while (nfigs > nx**2):
    nx+=1

print "nx=",nx,"nfigs=",nfigs

while (1):

    if len(sys.argv) == 1:
        break

    fname = sys.argv.pop(1)

    print "Loading:",fname
    ff = pyfits.open( fname )
    foundit=0
    for ii in range(100):
        im = ff[ii]
        if im.header["NAXIS"]>1:
            print "Found data at extension:",im.name
            foundit=1
            break
    if foundit==0:
        print "Couldn't find image in: ", fname
        ifig += 1
        continue

    f = maputils.FITSimage( externalheader=im.header, externaldata=im.data)


    frame = fig.add_subplot(nx,nx,ifig)
    ifig += 1

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


print "SHOWING..."
plt.show()



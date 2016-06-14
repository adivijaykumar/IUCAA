#################--IMPORTS--#####################

import numpy as np
from scipy import optimize
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.mplot3d import Axes3D
import os
import astropy
from astropy.wcs import WCS
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import SkyCoord
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('multipage.pdf')


##################--FUNCTIONS-RELATED-TO-PLOTTING-GAUSSIAN--#################
# These codes were taken and edited from http://scipy.github.io/old-wiki/pages/Cookbook/FittingData#Fitting_a_2D_gaussian #

matplotlib.rcParams.update({'font.size': 6}) #Updating font size for plotting

def gaussian(height, center_x, center_y, width_x, width_y,constant):
    """Returns a gaussian function with the given parameters"""
    
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2) + constant


def normal_gaussian(center_x,center_y,width_x,width_y):
    """Returns a normalized gaussian function with the given parameters"""
    
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)


def moments(data):
    """Returns (height, x, y, width_x, width_y,constant) - the gaussian parameters of a 2D distribution by calculating its moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    constant = np.mean(data)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    constant = np.median(data)
    return [height, x, y, width_x, width_y, constant]


def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y,constant)
    the gaussian parameters of a 2D distribution found by a least squares fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

########################--MAIN--############################

image_concat_plain = []
image_concat_subtracted = []
f = open('filenames.list') #filenames.list is a file containing the names of all files we need to process
f1 = open('info.list') #info.list contains RA, DEC details corresponding to RUN, RERUN, CAMCOL

for name in f:

    if name[0] == '#':
        continue
        # Do not include file in filenames.list whose first character is '#'. Done for flexibility of choosing sources


    #Obtain RUN, RERUN, CAMCOL and other details from the file name.

    l = name.split('-')
    x = l[2]
    g = l[3].split('.')
    from_file = [int(l[1]),int(x[1]),int(g[0])]
    trigger = 0

    # Obtain RA, DEC corresponding to RUN, RERUN, CAMCOL and other details. This and the last section might be buggy
    # and would have to changed depending on data sets and naming conventions

    f1 = open('info.list')
    for source in f1:
        a = source.split(',')
        from_info = [int(a[2]),int(a[4]),int(a[5])]
        if  from_file == from_info:
            ra = float(a[0])
            dec = float(a[1])
            trigger = 1 
            break
    
    f1.close()

    if trigger == 1: #ie if a match is detected
        
        ## Get Pixel Coordinates
        header = fits .getheader('./fits/'+name[:-1])
        a = WCS(header)
        c = FK5(float(ra),float(dec), unit='deg')
        output = astropy.wcs.utils.skycoord_to_pixel(c,a) #

        pix1 = np.floor(output[0])
        pix2 = np.floor(output[1])
        
        # Check if the peak coincides with centre and shift if it does not.

        temp = fits.getdata('./fits/'+name[:-1])
        temp = temp[pix2-5:pix2+6,pix1-5:pix1+6]
        new_pix = np.argmax(temp)
        new_pix = np.unravel_index(new_pix,(11,11))
        pix2 = pix2 + (new_pix[0] - 5)
        pix1 = pix1 + (new_pix[1] - 5)
        
        # Now with given peak value, form an image cutout
        image_from_file = fits.getdata('./fits/'+name[:-1])
        image_from_file = image_from_file[pix2-25:pix2+26,pix1-25:pix1+26]
        image_concat_plain.append(image_from_file) # Append to stack for final Stacking 
        #plt.imshow(image_from_file, cmap='Greys_r')

        # Execute program to obtain PSF. It is obtained as a 51x51 FIT image
        print '/home/aditya/Downloads/readAtlasImages-v5_4_11/read_PSF -v ./PS/psField-%s-%s-%s 3 %s %s foo.fit'%(l[1],l[2],l[3][:-1],pix1,pix2)
        #print output, name[:-1], ra, dec
        os.system('/home/aditya/Downloads/readAtlasImages-v5_4_11/read_PSF -v ./PS/psField-%s-%s-%s 3 %s %s foo.fit'%(l[1],l[2][1],l[3][:-1]   ,pix1,pix2))
        image = fits.getdata('foo.fit')
        dim = image.shape
        
        # Fit a Gaussian onto the PSF obtained, and get the parameters and image of the Gaussian        
        params = fitgaussian(image)
        #params = moments(image)
        #print 'The Parameters are : \nHeight = %s\nmu_x = %s\nmu_y = %s\nsigma_x = %s\nsigma_y = %s\nconstant = %s\n'%(params[0],params[1],params[2],params[3],params[4],params[5])
        normal_params = [params[1],params[2],params[3],params[4]]
        fit = normal_gaussian(*normal_params)
        Z = fit(*np.indices(image_from_file.shape))
        
        # Scale the Gaussian to coincide with the central peak of the image
        scale = np.amax(image_from_file)
        Z = (scale-np.median(image_from_file))*Z
        Z = Z + np.median(image_from_file)
        
        image1 = image_from_file - Z # Subtract reconstructed Gaussian from original image
        image_concat_subtracted.append(image1) # Append to stack for final Stacking 


#######--PERFORM-STACKING--########

final_image_plain = np.median(image_concat_plain,axis=0) 
final_image_subtracted = np.median(image_concat_subtracted,axis=0)

###########--PLOT-AS-REQUIRED--############

plt.imshow(final_image_plain,cmap='Greys_r')
plt.show()

params = fitgaussian(final_image_plain)
normal_params = [params[1],params[2],params[3],params[4]]
fit = normal_gaussian(*normal_params)
Z = fit(*np.indices(final_image_plain.shape))


scale = np.amax(final_image_plain)
Z = (scale-np.median(final_image_plain))*Z
Z = Z + np.median(final_image_plain)

final_image_plain = final_image_plain - Z

plt.imshow(final_image_plain,cmap='Greys_r')
plt.show()

plt.imshow(final_image_subtracted,cmap='Greys_r')
plt.show()


########--SAVE-AS-FITS-FILES-AS-REQUIRED--##########

print 'Making Outfile'
outfile = 'fits_plain.fits'
hdu = fits.PrimaryHDU(final_image_plain)
hdu.writeto(outfile, clobber=True)
print 'done!'

print 'Making Outfile'
outfile = 'fits_sub.fits'
hdu = fits.PrimaryHDU(final_image_subtracted)
hdu.writeto(outfile, clobber=True)
print 'done!'
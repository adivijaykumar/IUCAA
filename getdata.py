import astropy.io.fits as fits


f = open('list.txt','w')


hdulist = fits.open('dr12q_noabs.fits')
data = hdulist[1].data
RUN = data.field('RUN_NUMBER')
COL = data.field('COL_NUMBER')
FIELD = data.field('FIELD_NUMBER')
RERUN = data.field('RERUN_NUMBER')
RA = data.field('RA')
DEC = data.field('dec')


for i in range(len(RUN)):

	f.write('%s %s %s %s %s %s\n'%(RUN[i],COL[i],FIELD[i],RERUN[i],str(RA[i]),str(DEC[i])))
import astropy.io.fits as fits

for name in ['cop1.fits','cop2.fits','cop3.fits']:#['op3.fits','op4.fits']:#,'op5.fits']:#,'eop1.fits','eop2.fits']:
	
	print name

	f = open('input'+name[0]+name[-6]+'.txt','w')

	hdulist = fits.open(name,ignore_missing_end=True)
	for element in hdulist[1].data:
		f.write(str(element[0])+' '+str(element[1])+' '+str(element[5])+'\n')


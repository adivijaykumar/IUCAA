import numpy as np
import math
import matplotlib.pyplot as plt
import astropy.io.fits as fits

def ellipse(a,e):

	step = 0.2
	#sma = np.full((a/step)*2 +1,a)
	sma = a
	#print sma
	#ecc = np.full((a/step)*2 +1,e)
	ecc = e 
	#print ecc
	xdata = np.linspace((-1)*(1.0*a//step)*step,(1.0*a//step)*step,num = (a/step)*2 +1)
	#print xdata.shape
	xdata = np.insert(xdata,0,(-1)*sma)
	xdata = np.insert(xdata,xdata.shape[0],sma)
	#print xdata.shape
	#print xdata
	ydata1 = np.sqrt((1.0-ecc**2)*(sma**2 - xdata**2))
	ydata2 = (-1.0)*np.sqrt((1.0-ecc**2)*(sma**2 - xdata**2))
	return [xdata,ydata1,ydata2]


for name in ['inputc1.txt','inputc2.txt','inputc3.txt']:#,'input4.txt','input5.txt']:
	#data = fits.getdata('fits_sub_op1_1.fits')
	plt.figure()
	#plt.imshow(data,cmap='Greys_r')
	count = 0
	
	g=open(name[5]+name[6]+'.con','w')
	f= open(name,'r')
	pix_list = []
	for line in f:
		
		count +=1 
		#print count

		l = line.split()
		a = float(l[0])
		e = float(l[2])
		#print a,e
		if float(l[1]) > 0.001:

			pix_list.append(float(l[0]))
			#print l
			x, y1, y2 = ellipse(a,e)
			x = x + 25
			y1 = y1 + 25
			y2 = y2 + 25
		
			x_in = x+1
			y1_in = y1 +1
			y2_in = y2 + 1
			
			
			
			for i in range(x.shape[0]):
				g.write(' '+str(x_in[i])+' '+str(y1_in[i])+'\n')
		
			g.write('\n')
			for i in range(x.shape[0]):
				g.write(' '+str(x_in[i])+' '+str(y2_in[i])+'\n')
			g.write('\n')

	print name[5]+name[6], max(pix_list)
	g.close()

	#plt.show()

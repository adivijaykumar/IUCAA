import numpy as np

for name in ['inputc1.txt','inputc2.txt','inputc3.txt']:
	print name
	f = open(name,'r')
	g = open(name+'.r','w')

	for line in f:
		l = line.split()
		intensity = float(l[1])
		if  intensity > 0:

			temp = 23.9 - 2.5*np.log10(intensity*3.631)
			g.write('%s %s\n'%(l[0],str(temp)))

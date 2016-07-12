f = open('list.txt','r')
g = open('to-download-PSF','w')

no =0

for line in f:
	no +=1
	print no
	l = line.split()

	no_1 = '0'*(6-len(l[0]))+l[0]
	no_2 = '0'*(4-len(l[2]))+l[2]

#	g.write('%s,%s,%s\n'%(l[2],l[4],l[5]))
	#g.write('psField-%s-%s-%s.fit\n'%(no_1,l[4],no_
	g.write('http://data.sdss3.org/sas/dr12/boss/photo/redux/301/%s/objcs/%s/psField-%s-%s-%s.fit\n'%(l[0],l[1],no_1,l[1],no_2))

#	for color in ['r']:#,'g','r','i','z']:
#		g.write('http://dr12.sdss3.org/sas/dr12/boss/photoObj/frames/%s/%s/%s/frame-%s-%s-%s-%s.fits.bz2\n'%(l[3],l[0],l[1],color,no_1,l[1],no_2))
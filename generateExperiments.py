import random

for i in xrange(500):
	mutswappb = random.uniform(0.0001,0.025)
	cxswappb = random.uniform(0.1,1.0)
	cxop = random.choice(['upmx','pbx'])
	popsize = random.randint(5,100)
	tournSel = random.choice([True,False])
	total = random.choice([True,False])
	smallstart = False if total else random.choice([True,False])
	uniformsize = False if total or smallstart else random.choice([True,False])
	print "./optnetalign --net1 ~/Dropbox/optnetalign/tests/dmela.net --net2 ~/Dropbox/optnetalign/tests/hsapi.net --bitscores ~/Dropbox/optnetalign/tests/dm-hs.sim --ics --mutswappb " + str(mutswappb) + " --cxswappb " + str(cxswappb) + " --popsize " + str(popsize) + " " + ("--tournsel " if tournsel else " ") + ("--smallstart " if smallstart else " ") + (" --uniformsize " if uniformsize else " ") + ">> experimentalData.csv"
	print "rm *.aln"
	print "rm *.info"
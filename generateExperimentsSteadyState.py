import random

for i in xrange(500):
	mutswappb = random.uniform(0.0001,0.001)
	mutrate = random.uniform(0.01, 0.5)
	cxswappb = random.uniform(0.1,1.0)
	cxrate = random.uniform(0.1,0.9)
	hillclimbiters = int(random.uniform(100, 10000))
	dynparams = random.sample([True,False],1)
	print "./steadystateMOGA --net1 ~/Dropbox/optnetalign/optnetalign/tests/dmela.net --net2 ~/Dropbox/optnetalign/optnetalign/tests/hsapi.net --annotations1 ~/Dropbox/optnetalign/optnetalign/tests/dmela.annos --annotations2 ~/Dropbox/optnetalign/optnetalign/tests/hsapi.annos --s3 " + " --hillclimbiters " + str(hillclimbiters) + " --mutswappb " + str(mutswappb) + " --mutrate " + str(mutrate) + " --cxswappb " + str(cxswappb) + " --cxrate " + str(cxrate) + " --popsize 500 " + (" --dynparams " if dynparams else " ")  + "--nthreads 16 --nooutput --generations 1000000 --timelimit 60 --finalstats >> experimentalDataSteadyState.csv"
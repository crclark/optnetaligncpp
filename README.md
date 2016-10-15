optnetaligncpp
==============

Optimizing network aligner

### Quick start

Install the C++ Boost libraries on your system, then do the following:

```
make optnetalignubuntu
./optnetalign --net1 scere.net --net2 dmela.net --total --s3 --bitscores sc-dm.sim --blastsum --annotations1 scere.annos --annotations2 dmela.annos --cxrate 0.05 --cxswappb 0.75 --mutrate 0.05 --mutswappb 0.0001 --oneobjrate 0.75 --dynparams --popsize 200 --generations 1000000000 --hillclimbiters 10000  --timelimit 720 --outprefix sc-dm1 --finalstats >> sc-dm1.finalstats
```

For more information, see manual.tex

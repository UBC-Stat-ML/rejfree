Summary 
-------

Implementation of the [Bouncy Particle Sampler (BPS) algorithm](http://www.stat.ubc.ca/~bouchard/bps/) (focussing on the priority queue-based local BPS algorithm). 

Please include the following citation:

```
@article{bps,
	Author = {Bouchard-C\^{o}t\'{e}, A. and Vollmer, S.J. and Doucet, A.},
	Title = {The Bouncy particle sampler: a non-reversible rejection-free {M}arkov chain {M}onte {C}arlo method},
	Note = {Technical report arxiv:1510.02451},
	Year = {2015}
}
```


Installation
------------

There are several options available to install the package:

### Integrate to a gradle script

Simply add the following lines (replacing 1.0.0 by the current version (see git tags)):

```groovy
repositories {
 mavenCentral()
 jcenter()
 maven {
    url "http://www.stat.ubc.ca/~bouchard/maven/"
  }
}

dependencies {
  compile group: 'ca.ubc.stat', name: 'rejfree', version: '1.0.0'
}
```

### Compile using the provided gradle script

- Check out the source ``git clone git@github.com:alexandrebouchard/rejfree.git``
- Compile using ``gradle installApp``
- Add the jars in ``build/install/rejfree/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone git@github.com:alexandrebouchard/rejfree.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates
  

Usage
-----

We document here two usage scenarios: reproducing results from our paper, and creation of new models.


### Reproducing results

Complete reproducibility information is stored in separate git repos (containing precise git commits, plotting scripts, etc), but for most purpose the following information should be sufficient (those not listed here are implemented in Julia, and documented elsewhere; replace ``@@{x,y,z}`` by the cross product of all argument assignments or use [westrun](https://github.com/alexandrebouchard/westrun):

- Figure 4: 
    - main class: ``rejfree.models.normal.CompareStanRFOnNormalModel``
    - arguments (note: ignore the ``restrictVelocityNorm`` option as 
      the organization of this command line argument has been changed 
      since the time this experiment was ran):
    
```
-useLocal @@{true,false} \
-offDiag @@{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9} \
-nPairs 1000 \
-timeMilli 60000 \
-collectRate 0 \
-refreshRate @@{0.01,0.1,1,10} \
-restrictVelocityNorm @@{true,false}
```

- Figure 5:
    - main class: ``rejfree.models.normal.CompareStanRFOnNormalModel``
    - arguments:

```
-useLocal true \
-nPairs @@{100,1000} \
-timeMilli 30000 \
-collectRate 0 \
-refreshRate @@{0.01,0.1,1,10} \
-refreshmentMethod @@{GLOBAL,LOCAL,RESTRICTED,PARTIAL}
```

- Figure 6:
    - main class: ``rejfree.models.normal.CompareStanRFOnRN``
    - arguments:

```
++plans/singleRunOnRNOptions/options-@@{1,2,3,4}.map \
-stanHome /home/bouchard/bin/cmdstan-2.7.0 \
-restrictVelocityNorm @@{true,false} \
-usePartialRefreshment false
```

- Figure 7 and 8:
    - main class: ``rejfree.models.normal.CompareStanRFOnNormalModel``
    - arguments:

```
-recordPartialSums true \
-stanHome /Users/bouchard/Documents/workspace-cpp/cmdstan-2.7.0 \
-nPairs @@{10,100,1000} \
-nRepeats 40 \
-variableMonitorInterval 100
```

- Figure 9:
    - main class: ``rejfree.models.expfam.MRFMain``
    - arguments:

```
-stanHome /home/bouchard/bin/cmdstan-2.7.0 \
-nStanIters @@{2^[4--12]} \
-randomForBothAlgorithms @@{1--10} \
-collectRate 0.0 \
-useLocalRefreshment true \
-restrictVelocityNorm false
```

### Creating a new model

**Note:** we are in the process of making this more user-friendly and better documented.

- Create the factor graph
    - See ``rejfree.models.expfam.MRF`` for an example, in particular the inner class ``rejfree.models.expfam.MRF$ModelSpec``
    - The main step is to create an object where the factors (modelled by classes implementing ``rejfree.local.CollisionFactor``) are accessible in fields with the annotation ``@DefineFactor``. The structure of the factor graph is built using reflection by looking at which factors have shared accessibility of some objects. 
- To create and run the local BPS from the factor graph:
    - See ``rejfree.models.expfam.MRFMain`` for an example
    - The key snippet there is:
    
```java
ModelSpec modelSpec = mrf.newModelSpecFromGenerativeProcess(new Random(generateRandom));
LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);
rfRunner.init(modelSpec);
rfRunner.addMomentRayProcessor();
rfRunner.run();
```

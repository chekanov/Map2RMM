# Map2RMM

This is a library for mapping collision data to the RMM matrix format.
It uses input Monte Carlo data in the form of ProMC files from the [HepSim Monte Carlo repository](http://atlaswww.hep.anl.gov/hepsim/). The example builds anti-KT jets, and fill ROOT histograms.


 1. Install ProMC (http://atlaswww.hep.anl.gov/asc/promc/), ROOT and FastJet 
 2. Check the installation. The variables: 

```
   echo $PROMC
   echo $ROOTSYS
   echo $FASTJET
```
  they should return the installation paths. 

 3. Go to "map2rmm/" and compile the library as "make"
    This create 2 libraries :

```
      lib/libmap2rmm.so
      lib/libmap2rmm_static.a
```

 4. Go to the upper level and compile the example.cc  "make". The compilation will link the above library

 5. Download ProMC files from HepSim and put them to the "data" directory. Use hs-tools as: 
  
``` 
   hs-get tev100_higgs_ttbar_mg5 data
```
   See the HepSim documentation. 

 6. Process all files inside the directory "data" using the command "./example".

 7. Loook at the output root file "output.root" with histograms.
    The RMM data are stored as a tree "inputNN" with "proj" branch.

S.Chekanov (ANL) 


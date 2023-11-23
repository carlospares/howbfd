How to run:
* use python2
* python2 main.py -c config/config_file.py -N=npts
* atm only upwind flux is completed with the gf part, only for burgers
* added time to a bunch of functions types (funH, euler, RKk, BCs,nm_upwing, odiintegrator etc) 
* added MMS case to burgers !be careful: you have to change S in eq_burgers otherwise 


* Grid convergence  with   MMS 
* Grid convergence on steady state:
   ** initialise with steady exact
   **  run finite time and compute error
   ** do this with WENOk-AMp -- error == p+1=?
* Perturbations for fun


Then we pass to shallow water

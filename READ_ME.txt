How to run:
* use python2
* python2 main.py -c config/config_file.py -N=npts
* atm only upwind flux is completed with the gf part, only for burgers
* added time to a bunch of functions types (funH, euler, RKk, BCs,nm_upwing, odiintegrator etc) 
* added MMS case to burgers !be careful: you have to change S in eq_burgers otherwise 


DONE
* Grid convergence  with   MMS 
* Grid convergence on steady state:
   ** initialise with steady exact
   **  run finite time and compute error
   ** do this with WENOk-AMp -- error == p+1=?
* Perturbations for fun

To do 
* Latex the above


* add option to chose AB/AM outside and the steps outside
* add option to chose AB/AM or ABSW/AMSW outside
* With numerical  B'(x) show machine error Lake at rest 
* With numerical B'(x) grid convergence for super/subcritical (check if what we have is with analytical or numerical)
* perturbations of above and also super.sub critical (compare analytical and numerical b(x) )
* MS for shallow water

* discontinuous H: first order ok, missing AMK and AB
* discontinuous H: first order ok, missing extrapolation
* discontinuous H: Carlo's stuff
* discontinuous H: shallow water

* Then we will focus on trans critical with Carlos

* Then friction

* Then another system



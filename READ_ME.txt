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

* With analitycal B'(x) grid convergence for Lake at rest 
* With analitycal B'(x) grid convergence for super/subcritical
* add option to chose AB/AM outside and the steps outside
* add option to chose AB/AM or ABSW/AMSW outside
* fix size of sumHx 
* copy/past the GF with B'(x) numerical for AB/AM from 4 to 8
* With numerical  B'(x) show machine error Lake at rest 
* With numerical B'(x) grid convergence for super/subcritical
* perturbations of above and also super.sub critical
* MS for shallow water

* Then we will focus on trans critical with Carlos

* Then friction

* Then another system

Then we pass to shallow water

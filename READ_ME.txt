How to run:
* use python3 (or 2 if it works)
* pythonX main.py -c config/config_file.py -N=npts
* atm only upwind flux is completed with the gf part, only for burgers
* added time to a bunch of functions types (funH, euler, RKk, BCs,nm_upwing, odiintegrator etc) 
* added MMS case to burgers !be careful: you have to change S in eq_burgers otherwise 

* git ff conflict :   git config pull.rebase false 

DONE
* Grid convergence  with   MMS 
* Grid convergence on steady state:
   ** initialise with steady exact
   **  run finite time and compute error
   ** do this with WENOk-AMp -- error == p+1=?
* Perturbations for fun
* add option to chose AB/AM outside and the steps outside
* add option to chose AB/AM or ABSW/AMSW outside

* discontinuous H: first order ok, missing AMK and AB
* discontinuous H: first order ok, missing extrapolation
* discontinuous H: Carlo's stuff


To do 
* Latex the above



* Then we will focus on trans critical with Carlos
* add friction? 
* reversed flow case


* Euler: check mismatches just to be able to run
* Euler:  look into case of figure 4 in https://arxiv.org/pdf/2307.12089
* Euler:  compute jumps to linearize with discontinuous section






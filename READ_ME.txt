How to run:
* use python3 (or 2 if it works)
* pythonX main.py -c config/config_file.py -N=npts
* atm only upwind flux is completed with the gf part, only for burgers
* added time to a bunch of functions types (funH, euler, RKk, BCs,nm_upwing, odiintegrator etc) 
* added MMS case to burgers !be careful: you have to change S in eq_burgers otherwise 

* git ff conflict :   git config pull.rebase false 

DONE

To do 

* Euler: check mismatches just to be able to run
* Euler:  look into case of figure 4 in https://arxiv.org/pdf/2307.12089
* Euler:  compute jumps to linearize with discontinuous section






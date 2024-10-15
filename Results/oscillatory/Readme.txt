BURGERS steady oscillatory solution. 
S(U)=U^2
H = x+ 0.1*sin(100*x)

folders T0p9 and T1P0 we run the solution to time 0.9 and 1.0 respectively. We used the AM unless stated differently in the name of the file.
all the *txt files contain :
x, u_{in}, u_{final}
-----------------------------------------------------------------
folder pertrurbation_0P005
*.txt files :
x, u_{in}, u_{final} (i don't remeber the t_{final})

u_{in} is the initial solution exp(H)+ the GAUSSIAN pertrubation 10^{-3}

*_stationary.txt

x, u_{in}, u_{final} (i don't remeber the t_{final})

u_{in} is the initial solution exp(H)
u_{final} the convergent solution

CHECK THE plotinpython.py file
-----------------------------------------------------------------
folder pertrurbation_0P05

Same but with a bigger GAUSSIAN pertrubation
-----------------------------------------------------------------
folder perturbation_0P1 and 0P3

Same. We kept also a refernce solution reference_p.txt (with the pertrubation) using the upwind scheme

CHECK THE plotinpython.py file
-----------------------------------------------------------------
folder pertrurbation_DISC_t0p7

Dicontinious pertrubtion 

plotinpython.py : plots the difference of the stationary solution and the one with the propageted pertrubation 

-----------------------------------------------------------------


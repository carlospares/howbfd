# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import newton
from Funciones_salto_estacionario import phi
from equation import Equation
from nosteadyexc import NoSteadyError

class EulerEquationGRAV(Equation):
    """ 1D Euler  equation, vars [rho,q=rho u,  E],  E = rho e + rho u^2/2 = p/(gm-1)  + rho u^2/2

        rho_t +  q_x               = 0
        q_t + (q^2/rho + p)_x = -rho H_x
        E_t + (uE + pu/rho)_x = -rho u Hx
        H is the gravitational potential (equivalent of -g b for  SW)
    """
    
    # other function parameters
    g = 1.4 # heat coefficient ratio - perfect gas
    fcoeff = 0.0

    def F(self, U):
        """ Flux function """
        ret = np.empty(U.shape)
        rho = U[0,:]
        q   = U[1,:]
        E   = U[2,:]
        
        u = q/rho
        p = ( self.g - 1. )*( E - q*u*0.5 )
        
        ret[0,:] = q
        ret[1,:] = u*q +  p
        ret[2,:] = u*E + p*u
        return ret

    def Fp(self, U):
        """ Flux function """
        ret = np.empty(U.shape)
        rho = U[0,:]
        q   = U[1,:]
        E   = U[2,:]
        
        u = q/rho
        p = ( self.g - 1. )*( E - q*u*0.5 )
        a = np.sqrt(self.g*p/rho)
        
        
        l1 = u - a
        l2 = u
        l3 = u + a
        
        r11 = 1.
        r21 = u - a
        r31 = E/rho + p/rho - u*a
        
        r12 = 1.
        r22 = u
        r32 = u*u*0.5
        
        r13 = 1.
        r23 = u + a
        r33 = E/rho + p/rho + u*a
        
        
        l11 = 0.5*( 0.5*(self.g-1.)*u/a + 1. )*u/a
        l12 = -0.5*( (self.g-1.)*u/a + 1. )/a
        l13 = 0.5*(self.g-1.)/( a*a )
        
        l21 = 1. - 0.5*(self.g-1.)*u*u/(a*a)
        l22 = 0.5*(self.g-1.)*u/(a*a)
        l23 = -0.5*(self.g-1.)/(a*a)
        
        l31 = 0.5*( 0.5*(self.g-1.)*u/a - 1. )*u/a
        l32 = -0.5*( (self.g-1.)*u/a - 1. )/a
        l33 = l13
        
        l1p = l1*(l1 >= 0)
        l2p = l2*(l2 >= 0)
        l3p = l3*(l3 >= 0)
        
        l11p = l1p*l11
        l12p = l1p*l12
        l13p = l1p*l13
        
        l21p = l2p*l21
        l22p = l2p*l22
        l23p = l2p*l23
        
        l31p = l3p*l31
        l32p = l3p*l32
        l33p = l3p*l33
        
        vv1 = l11p*rho + l12p*q + l13p*E
        vv2 = l21p*rho + l22p*q + l23p*E
        vv3 = l31p*rho + l32p*q + l33p*E
        
        ret[0,:] = r11*vv1 + r12*vv2 + r13*vv3
        ret[1,:] = r21*vv1 + r22*vv2 + r23*vv3
        ret[2,:] = r31*vv1 + r32*vv2 + r33*vv3
        return ret

    def Fm(self, U):
        """ Flux function """
        ret = np.empty(U.shape)
        rho = U[0,:]
        q   = U[1,:]
        E   = U[2,:]
        
        u = q/rho
        p = ( self.g - 1. )*( E - q*u*0.5 )
        a = np.sqrt(self.g*p/rho)
        
        
        l1 = u - a
        l2 = u
        l3 = u + a
        
        r11 = 1.
        r21 = u - a
        r31 = E/rho + p/rho - u*a
        
        r12 = 1.
        r22 = u
        r32 = u*u*0.5
        
        r13 = 1.
        r23 = u + a
        r33 = E/rho + p/rho + u*a
        
        
        l11 = 0.5*( 0.5*(self.g-1.)*u/a + 1. )*u/a
        l12 = -0.5*( (self.g-1.)*u/a + 1. )/a
        l13 = 0.5*(self.g-1.)/( a*a )
        
        l21 = 1. - 0.5*(self.g-1.)*u*u/(a*a)
        l22 = 0.5*(self.g-1.)*u/(a*a)
        l23 = -0.5*(self.g-1.)/(a*a)
        
        l31 = 0.5*( 0.5*(self.g-1.)*u/a - 1. )*u/a
        l32 = -0.5*( (self.g-1.)*u/a - 1. )/a
        l33 = l13
        
        l1m = l1*(l1 <= 0)
        l2m = l2*(l2 <= 0)
        l3m = l3*(l3 <= 0)
        
        l11m = l1m*l11
        l12m = l1m*l12
        l13m = l1m*l13
        
        l21m = l2m*l21
        l22m = l2m*l22
        l23m = l2m*l23
        
        l31m = l3m*l31
        l32m = l3m*l32
        l33m = l3m*l33
        
        vv1 = l11m*rho + l12m*q + l13m*E
        vv2 = l21m*rho + l22m*q + l23m*E
        vv3 = l31m*rho + l32m*q + l33m*E
        
        ret[0,:] = r11*vv1 + r12*vv2 + r13*vv3
        ret[1,:] = r21*vv1 + r22*vv2 + r23*vv3
        ret[2,:] = r31*vv1 + r32*vv2 + r33*vv3
        return ret
    

    def dF(self, U):
        """ Derivative of flux function """
        DF = np.zeros((U.shape[1], self.dim(), self.dim()))
        for i in range(U.shape[1]):
            rho = U[0,i]
            q   = U[1,i]
            E   = U[2,i]
            u = q/rho
            k = 0.5*u*u
            gm1 = self.g-1.
            gm3 = self.g-3.
            p = (self.g-1.)*( E - rho*k )
            Hs = E/rho  + p/rho
            DF[i] = np.array([[0,1,0],[gm3*k,-gm3*u,gm1],[u*(gm1*k-Hs),Hs-2.*gm1*k,self.g*u]])
        return DF

    def eig_of_dF(self, U):
        """ Returns eigenvalues of dF(U), as
            numpy array shaped like U, with
            eig[0,i] < eig[1,i] < ... for all i """
        eig = np.zeros(U.shape)
        rho = U[0,:]
        q   = U[1,:]
        E   = U[2,:]
        
        u = q/rho
        p = ( self.g - 1. )*( E - q*u*0.5 )
        a = np.sqrt(self.g*p/rho)
        eig[0,:] = u - a
        eig[1,:] = u
        eig[2,:] = u + a
        return eig


    def S(self, U):
        """ Return S(U) """
        return np.array([ 0, -U[0], - U[1] ])
        
    def sigma(self, U):
        """ Return sigma(U) """
        return np.array([ 0, self.fcoeff*U[1], self.fcoeff*U[1]*U[1]/U[0]  ])

#    def upw_criterion(self, uStencil):
#        """ Returns a pair (l,r) with the velocity for upwind criterion at
#            left and right intercells, for cell at center of uStencil
#        """
#        print "[ERROR] Upwind only implemented for scalar equations!"
#        raise NotImplementedError

    def discH_jumpF(self, ui, uip1, i, dH, x, t):
        # depends on S
        #hbar   = 0.5*( ui[0] + uip1[0] )
        #hubar  = 0.5*( ui[1] + uip1[1] )
        #ratio  = hubar*hubar/( self.g*(ui[0]*uip1[0])*(ui[0]*uip1[0]))
        #num    = ratio*( hbar*hbar - ui[0]*uip1[0] )
        #den    = 1.0 - ratio*hbar
        #htilde = hbar + num/den
        #
        # ATTENTION: to be coded
        
        delta = dH
        return delta
        
    def Piplus(self, ui, uip1):
        rhol = ui[0]
        ql   = ui[1]
        ul   = ui[1]/ui[0]
        El   = ui[2]
        kl   = 0.5*ul*ul
        rhor = uip1[0]
        qr   = uip1[1]
        ur   = uip1[1]/uip1[0]
        Er   = uip1[2]
        kr   = 0.5*ur*ur
        
        rho = .5*(np.sqrt(rhol) + np.sqrt(rhor))
        rho = rho*rho
        u   = (np.sqrt(rhol)*ul + np.sqrt(rhor)*ur)/(np.sqrt(rhol) + np.sqrt(rhor))
        Hs  = np.sqrt(rhol)*( self.g*El - (self.g-1.)*kl ) + np.sqrt(rhor)*( self.g*Er - (self.g-1.)*kr )
        Hs  = Hs/(np.sqrt(rhol) + np.sqrt(rhor))
        k   = 0.5*u*u
        
        p = (self.g-1.)*( Hs - 0.5*u*u )/self.g
        a = np.sqrt(self.g*p/rho)
        
        l1 = u - a
        l2 = u
        l3 = u + a
        
        r11 = 1.
        r21 = u - a
        r31 = Hs - u*a
        
        r12 = 1.
        r22 = u
        r32 = u*u*0.5
        
        r13 = 1.
        r23 = u + a
        r33 = Hs + u*a
        
        
        l11 = 0.5*( 0.5*(self.g-1.)*u/a + 1. )*u/a
        l12 = -0.5*( (self.g-1.)*u/a + 1. )/a
        l13 = 0.5*(self.g-1.)/( a*a )
        
        l21 = 1. - 0.5*(self.g-1.)*u*u/(a*a)
        l22 = 0.5*(self.g-1.)*u/(a*a)
        l23 = -0.5*(self.g-1.)/(a*a)
        
        l31 = 0.5*( 0.5*(self.g-1.)*u/a - 1. )*u/a
        l32 = -0.5*( (self.g-1.)*u/a - 1. )/a
        l33 = l13
        
        l1p = .5*(1 + np.sign(l1))
        l2p = .5*(1 + np.sign(l2))
        l3p = .5*(1 + np.sign(l3))
        
        l11p = l1p*l11
        l12p = l1p*l12
        l13p = l1p*l13
        
        l21p = l2p*l21
        l22p = l2p*l22
        l23p = l2p*l23
        
        l31p = l3p*l31
        l32p = l3p*l32
        l33p = l3p*l33
        
        A11 = r11*l11p + r12*l21p + r13*l31p
        A12 = r11*l12p + r12*l22p + r13*l32p
        A13 = r11*l13p + r12*l23p + r13*l33p
                
        A21 = r21*l11p + r22*l21p + r23*l31p
        A22 = r21*l12p + r22*l22p + r23*l32p
        A23 = r21*l13p + r22*l23p + r23*l33p
        
        A31 = r31*l11p + r32*l21p + r33*l31p
        A32 = r31*l12p + r32*l22p + r33*l32p
        A33 = r31*l13p + r32*l23p + r33*l33p
        
        A = np.array([[A11, A12, A13], [A21, A22, A23],[A31,A32,A33]])
        return A
    
    def Piminus(self, ui, uip1):
        hl = ui[0]
        ql = ui[1]
        ul = ui[1]/ui[0]
        El = ui[2]
        hr = uip1[0]
        qr = uip1[1]
        ur = uip1[1]/uip1[0]
        Er = uip1[2]
        
        h = .5*(np.sqrt(hl) + np.sqrt(hr))
        h = h*h
        u = (np.sqrt(hl)*ul + np.sqrt(hr)*ur)/(np.sqrt(hl) + np.sqrt(hr))
        Hs = np.sqrt(hl)*( self.g*El - (self.g-1.)*0.5*ql*ql/(hl*hl) ) + np.sqrt(hr)*( self.g*Er - (self.g-1.)*0.5*qr*qr/(hr*hr) )
        Hs = Hs/(np.sqrt(hl) + np.sqrt(hr))
        
        p = (self.g-1.)*( Hs - 0.5*u*u )/self.g
        a = np.sqrt(self.g*p/h)
        
        
        l1 = u - a
        l2 = u
        l3 = u + a
        
        r11 = 1.
        r21 = u - a
        r31 = Hs - u*a
        
        r12 = 1.
        r22 = u
        r32 = u*u*0.5
        
        r13 = 1.
        r23 = u + a
        r33 = Hs + u*a
        
        
        l11 = 0.5*( 0.5*(self.g-1.)*u/a + 1. )*u/a
        l12 = -0.5*( (self.g-1.)*u/a + 1. )/a
        l13 = 0.5*(self.g-1.)/( a*a )
        
        l21 = 1. - 0.5*(self.g-1.)*u*u/(a*a)
        l22 = 0.5*(self.g-1.)*u/(a*a)
        l23 = -0.5*(self.g-1.)/(a*a)
        
        l31 = 0.5*( 0.5*(self.g-1.)*u/a - 1. )*u/a
        l32 = -0.5*( (self.g-1.)*u/a - 1. )/a
        l33 = l13
        
        l1m = .5*(1 - np.sign(l1))
        l2m = .5*(1 - np.sign(l2))
        l3m = .5*(1 - np.sign(l3))
        
        l11m = l1m*l11
        l12m = l1m*l12
        l13m = l1m*l13
        
        l21m = l2m*l21
        l22m = l2m*l22
        l23m = l2m*l23
        
        l31m = l3m*l31
        l32m = l3m*l32
        l33m = l3m*l33
        
        A11 = r11*l11m + r12*l21m + r13*l31m
        A12 = r11*l12m + r12*l22m + r13*l32m
        A13 = r11*l13m + r12*l23m + r13*l33m
                
        A21 = r21*l11m + r22*l21m + r23*l31m
        A22 = r21*l12m + r22*l22m + r23*l32m
        A23 = r21*l13m + r22*l23m + r23*l33m
        
        A31 = r31*l11m + r32*l21m + r33*l31m
        A32 = r31*l12m + r32*l22m + r33*l32m
        A33 = r31*l13m + r32*l23m + r33*l33m
        
        A = np.array([[A11, A12, A13], [A21, A22, A23],[A31,A32,A33]])
        return A
    
    def dim(self):
        """ Returns dimension of the problem: 1 for scalars """
        return 3

    def steady(self, H,x):
        U0 = np.ones((self.dim(), len(H)))
        """ Returns an arbitrary steady state for the equation.
            Input: 
                x: spatial coordinates
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work!
        """
# # Water at rest solution: [h, q](x) = [eta + H(x), 0]
        #arbitrary_eta = 1
        #U0 = np.zeros((2, len(x)))
        #U0[0,:] = arbitrary_eta + self.H(x)
        #return U0
# BUMP2
#        HConst = 0.
#        qConst = 2.5
#        hConst = 2.
        
#        HConst = 0.13
#        qConst = 1.
#        hConst = 1.

#BUMPS
        #----supercritical
#        HConst = 0.
#        qConst = 24.
#        hConst = 2.

        #----subcritical
        HConst = 0.
        qConst = 4.42
        hConst = 2.

#BUMPT
        #----transcritical with shock 
#        HConst = 0.
#        qConst = 0.18
#        hConst = 0.33

        #----transcritical without shock 
#        HConst = 0.
#        qConst = 1.53
#        hConst = 0.4057809453450358#0.66
        
#BUMPD  
#        HConst = -.5
#        qConst = 2.5
#        hConst = .5
#        U0[1,:] = qConst
#        hConst = 2.

        uConst = [hConst, qConst]
        return self.steady_constraint(HConst, uConst, H,x, U0)
    
    def steady_trans(self,H,x): 
        U0 = np.ones((self.dim(), len(H)))
        g = self.g
        HConst = -.5
        qConst = 2.5
        U0[1,:] = qConst
        hConst = np.abs(qConst)**(2/3.)/g**(1/3.)
        uConst = [hConst, qConst]
        return self.steady_constraint(HConst, uConst, H,x, U0)
#        print 'steady',index, H[index]
#        index = np.argmin(np.abs(H-HConst)) 
#        U0 = np.zeros((self.dim(), len(H)))
#        for j in range(len(H)):
#            # phi may fail to find a steady state even if it formally exists
#            # when flow is close to critical. This is a failsafe for that.
#            try: 
#                (hsuperc, hsubc) = phi(hConst, uConst, HConst, H[j])
#            except Exception:
#                raise NoSteadyError("Steady state exists but failed to find it. Too close to critical flow?\
#                                    (Hi={}, hi={}, ui={}), H-Hstar={}".format(HConst, hConst, uConst, H[j]-Const ))
#            # hstar = polyNewton[j] # Halley's method
#            U0[0,j] = hsuperc if j >= index else hsubc
#            U0[1,j] = qConst
#        return U0

    def twosteady(self,H,x):
        U0 = np.ones((self.dim(), len(H)))
        cond= np.where(x>=0)
        if len(cond[0]) == 0:
            ind = len(H)+1
        else:
            ind = cond[0][0]
        qConst = 1.
        if ind > 0:
            HConstL = 1.
            hConstL = 1.
            uConstL = [hConstL, qConst]
            U0[:,0:ind] = self.steady_constraint(HConstL,uConstL, H[0:ind], x[0:ind], U0[:,0:ind])
        if ind < len(H) + 1:
            HConstR = 11.
            hConstR = 9.
            uConstR = [hConstR, qConst]
            U0[:,ind:] = self.steady_constraint(HConstR,uConstR, H[ind:], x[ind:], U0[:,ind:])

        return U0

    def Froude(self, u, h):
        return abs(u) / np.sqrt(self.g*h)

    def Pi(self, V):
        """ Return projection of V into the space of non-conservative subsys.
            (i.e. Pi([4,2]) = [0,2]) """
        PV = np.zeros_like(V)
        PV[1] = V[1]
        return PV
    
    def critical_H(self, Hc, qc, hc):
        """ Let P be the polynomial for which, if (q_i*, h_i*) is a steady state solution,
            P(h_i*) = 0. H* is the value such that:
            H > H* : two roots exist (sub/supercritical)
            H = H* : one roots exists (critical)
            H < H* : no root exists (ie no steady state solution)
        """
        return Hc + 1.5*(qc*qc/self.g)**(1./3) - (0.5/self.g)*qc*qc/(hc*hc) - hc

    def steady_constraint(self, HConstr, uConstr, H,x, U0):
        """ Returns a steady state solution of the equation, u*, constrained
            to u*(xConstr) = uConstr
            Input:
                xConstr: double, x to fix the constraint
                uConstr: double or (dims,1) np array, u*(x) = uConstr for u* we look for
                x: values of x at which to evaluate u*
            Output:
                (nvars, len(x)) numpy array with the values
                If nvars = 1, this must still be a (1,len(x)) matrix;
                a len(x) array will not work! """
        Ustar = np.zeros((self.dim(), len(H)))
        (hi, qi, ui) = (uConstr[0], uConstr[1], uConstr[1]/uConstr[0])
        Hstar = self.critical_H(HConstr, qi, hi)
        if np.min(H) < Hstar - 1.e-6:
            print ('no existe solucion estacionaria')
            raise NoSteadyError("No steady state exists for constraint \
                                (Hi={}, hi={}, ui={}), H={}"\
                                .format(HConstr, hi, ui, [myH for myH in H if myH<Hstar] ))
        #polyNewton = self.solve_steady_poly_newton(qi, hi, H)
        #print self.find_steady_poly_roots(qi, hi, HConstr, H)

        Fr_i = self.Froude(ui, hi)
        for j in range(len(H)):
            # phi may fail to find a steady state even if it formally exists
            # when flow is close to critical. This is a failsafe for that.
            try: 
                (hsuperc, hsubc) = phi(hi, ui, HConstr, H[j], U0[0,j])
            except Exception:
                print ('no se ha encontrado la solucion estacionaria')
                raise NoSteadyError("Steady state exists but failed to find it. Too close to critical flow?\
                                    (Hi={}, hi={}, ui={}), H-Hstar={}".format(HConstr, hi, ui, H-Hstar ))
            # hstar = polyNewton[j] # Halley's method
            Ustar[0,j] = hsuperc if Fr_i > 1 else hsubc
#            Ustar[0,j] = hsuperc if x[j] > -1.2 else hsubc  # transcritical stationary solution with critical point at x = 0
            Ustar[1,j] = uConstr[1]
#        i = (U0.shape[1]-1)/2
#        if x[i]== -1.05:
#            print ''
#            print U0
#            print Ustar
        return Ustar

    # # TO DO: very slow. Delete me?
    # def find_steady_poly_roots(self, qi, hi, Hi, Hj):
        # #i = (len(Hj)-1)/2 # central point is the constraint
        # Ci = 0.5*qi*qi/hi/hi + self.g*(hi - Hi)
        # # hstar = np.empty_like(Hj)
        # # for j in range(len(Hj)):
        # #     hstar[j] = max(np.roots([self.g, -(Ci+self.g*Hj[j]), 0, 0.5*qi*qi]))
        # # return hstar
        # hstars = np.zeros((len(Hj), 3), dtype=complex)
        # for j in range(len(Hj)):
            # hstars[j] = np.roots([self.g, -(Ci+self.g*Hj[j]), 0, 0.5*qi*qi])
        # return hstars
    
    # def solve_steady_poly_newton(self, qi, hi, Hj):
        # """ find h*(x) at x_j so that (h*, qi) is a steady state.
            # Input:
                # qi: double, qi = q(x_i) = q*(x_i) = q(x_j) forall j in stencil
                # hi: double, hi = h(x_i) = h*(x_i)
                # Hj: numpy array with [H(x_j)] forall j in stencil
            # Output:
                # numpy array with [h*(x_j)] forall j in stencil
        # """
        # i = (len(Hj)-1)/2 # central point is the constraint
        # Ci = 0.5*qi*qi/hi/hi + self.g*(hi - Hj[i])
        # hstar = np.empty_like(Hj)
        # for j in range(len(Hj)):
            # hstar[j] = newton(steady_poly, hi, args=(qi, Ci, Hj[j], self.g), 
                              # fprime=steady_poly_prime, 
                              # fprime2=steady_poly_second) # for Halley's method
        # return hstar

    def prepare_plot(self,x,u,H,t):
        """ Plot x and u in whichever way is appropriate for the equation.
            This function will be called by io_manager.
            This function should produce a finished plot, including title,
            legend, labels and so on; io_manager will do plt.show() or savefig() 
            as required """
        plt.subplot(311)
        plt.title(t)
        plt.plot(x, -H, 'b', label='-H')
        plt.plot(x, u[0], 'g', label='$rho$')
#        plt.plot(x, u[0], 'r', label='h')
        plt.legend()
        plt.subplot(312)
        plt.plot(x, u[1]/u[0], label='u')
        #plt.plot(x, u[1], label='rho u')
        plt.legend()
        plt.subplot(313)
        plt.plot(x, u[2], label='E')
        plt.legend()

    def exact(self, x, t, H, params):
        return self.steady(H,x)
# def steady_poly(h, q, Ci, Hj, g):
    # """ Polynomial in h, such that P(h) = 0 iff h is a value at x_j """
    # return g*h*h*h - (Ci+g*Hj)*h*h + 0.5*q*q
# def steady_poly_prime(h, q, Ci, Hj, g):
    # """ Derivative in h of steady_poly """
    # return 3*g*h*h - 2*(Ci+g*Hj)*h
# def steady_poly_second(h, q, Ci, Hj, g):
    # """ Second derivative in h of steady_poly """
    # return 6*g*h - 2*(Ci+g*Hj)
    

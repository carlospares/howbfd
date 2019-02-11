from equation import Equation

class Flux:
	UPWIND = 0
	RUSANOV = 1

	def __init__(self, flux):
		self.numflux = flux

	def flux(self, fluxL, fluxR, eq, uL, uR):
		""" Computes the numerical flux
			Input:
				fluxL: flux reconst. at left of interface
				flurR: same at right of interface
				eq: object of class Equation
				uL: current value at cell to left of interface
				uR: same at right of interface
		"""
		if self.numflux == self.UPWIND:
			crit = eq.upw_criterion(uL, uR)
			return self.upwind(fluxL, fluxR, crit)
		elif self.numflux == self.RUSANOV:
			return self.rusanov(fluxL, fluxR, uL, uR, eq)

	def upwind(self, fluxL, fluxR, crit):
	    return fluxL if crit>=0 else fluxR

	def rusanov(self, fluxL, fluxR, uL, uR, eq):
		return 0.5*(fluxL + fluxR) - \
				max(abs(eq.dF(uL)), abs(eq.dF(uR)))*0.5*(uR-uL)
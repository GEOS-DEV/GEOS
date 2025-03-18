
class RateAndStateParameters:
    def __init__(self):
        self.amax = 0.016
        self.a0 = 0.004
        self.b = 0.01
        self.Vstar = 1.0e-6
        self.f = 0.6
        self.Rvw = 200
        self.Drs = 5.0e-4
        self.shearModulus = 32.04e9 
        self.cs = 3.464e3

    def a(self, r):
        if (r > self.Rvw ):
            return self.amax
        else:
            return self.a0    

    def shearImpedance(self):
        return self.shearModulus / (2*self.cs)     
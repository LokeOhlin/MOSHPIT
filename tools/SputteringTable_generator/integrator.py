import numpy as np
import sys

class dopri5_integrator():
    def __init__(self, function, *args, y_independent = False, rtol = 1e-6, atol = 1e-6, min_stepsize = 1e-10): 
        self.cs = np.array([1/5, 3/10, 4/5, 8/9, 1, 1])
        self.aij = np.zeros((6,7))
        
        self.aij[0,0] = 1/5
        self.aij[1,:2] = np.array([3/40      ,  9/40])
        self.aij[2,:3] = np.array([44/45     , -56/15     , 32/9])
        self.aij[3,:4] = np.array([19372/6561, -25360/2187, 64448/6561, -212/729])
        self.aij[4,:5] = np.array([9017/3168 , -355/33    , 46732/5247, 49/176  , -5103/18656])
        self.aij[5,:6] = np.array([35/384    , 0          , 500/1113  , 125/192 , -2187/6784 , 11/84])

        self.b2 = np.array([5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40])

        self.x = 0
        self.y = np.array([0])
        self.args = args
        self.function = function
        self.atol = atol
        self.rtol = rtol
        self.min_stepsize = min_stepsize
        self.y_independent = y_independent
    def set_ic(self, x, y):
        self.x = x
        self.y = np.atleast_1d(y).astype(np.float64)
        self.ys = np.zeros((7,) +self.y.shape, dtype = np.float64)
    def step(self, dx):
        self.ys[:] = 0
        self.ys[0] = self.function(self.x, self.y, *self.args)
        self.ys[1] = self.function(self.x + self.cs[0]*dx, self.y + dx*np.sum(self.ys*self.aij[0,:,None], axis = 0), *self.args)
        self.ys[2] = self.function(self.x + self.cs[1]*dx, self.y + dx*np.sum(self.ys*self.aij[1,:,None], axis = 0), *self.args)
        self.ys[3] = self.function(self.x + self.cs[2]*dx, self.y + dx*np.sum(self.ys*self.aij[2,:,None], axis = 0), *self.args)
        self.ys[4] = self.function(self.x + self.cs[3]*dx, self.y + dx*np.sum(self.ys*self.aij[3,:,None], axis = 0), *self.args)
        self.ys[5] = self.function(self.x + self.cs[4]*dx, self.y + dx*np.sum(self.ys*self.aij[4,:,None], axis = 0), *self.args)
        yfifth = self.y + dx * np.sum(self.ys*self.aij[5,:,None], axis = 0)
        self.ys[6] = self.function(self.x + self.cs[5]*dx, yfifth, *self.args)
        ysixth = self.y + dx * np.sum(self.ys*self.b2[:,None], axis = 0)  
        
        diff = np.abs(yfifth - ysixth)
        tol = self.atol + self.rtol*np.maximum(np.abs(yfifth), np.abs(ysixth))
        error = np.sqrt(1/self.y.size * np.sum((diff/tol)**2))
        if error == 0:
            dx_new = dx * 1.2
        else:
            dx_new = 0.8 * min(error**(-1/4), 1.5)*dx
        if np.any(diff > tol):
            return dx_new, False
        else:
            self.y = yfifth
            return dx_new, True

    def step_y_independent(self, dx):

        self.ys[:] = 0
        self.ys[0] = self.function(self.x, 0, *self.args)
        self.ys[1] = self.function(self.x + self.cs[0]*dx, 0, *self.args)
        self.ys[2] = self.function(self.x + self.cs[1]*dx, 0, *self.args)
        self.ys[3] = self.function(self.x + self.cs[2]*dx, 0, *self.args)
        self.ys[4] = self.function(self.x + self.cs[3]*dx, 0, *self.args)
        self.ys[5] = self.function(self.x + self.cs[4]*dx, 0, *self.args)
        yfifth = self.y + dx * np.sum(self.ys*self.aij[5,:,None], axis = 0)
        self.ys[6] = self.function(self.x + self.cs[5]*dx, 0, *self.args)
        ysixth = self.y + dx * np.sum(self.ys*self.b2[:,None], axis = 0)  
        
        diff = np.abs(yfifth - ysixth)
        tol = self.atol + self.rtol*np.maximum(np.abs(yfifth), np.abs(ysixth))
        error = np.sqrt(1/self.y.size * np.sum((diff/tol)**2))
        if error == 0:
            dx_new = dx * 1.2
        else:
            dx_new = 0.8 * min(error**(-1/4), 1.5)*dx
        if np.any(diff > tol):
            return dx_new, False
        else:
            self.y = yfifth
            return dx_new, True

    def integrate(self, xend, first_step = 1e-3):
        dx_tot = xend - self.x
        dx_sign = np.sign(dx_tot)
        # TODO: make this settable
        dx = dx_tot * first_step
        while dx_sign*self.x < dx_sign*xend:
            if(dx < self.min_stepsize*dx_tot):
                print("WARNING: dx has reached the minimum stepsize for x = %f with xend = %f"%(self.x, xend))
                sys.exit(0)
                break
            if self.y_independent:
                dx_new, success = self.step_y_independent(dx)
            else:
                dx_new, success = self.step(dx)
            if success:
                self.x += dx
            dx = dx_new
            if(dx_sign*dx > dx_sign*(xend - self.x)):
                dx = xend - self.x
    
        return self.y


if __name__ == "__main__":
    def sin(x,integral):
        return np.sin(x)

    integrator = dopri5_integrator(sin)
    integrator.set_ic(0, -1)
    from scipy.integrate import ode as ode
    import time

    sci = ode(sin).set_integrator("dopri5")
    sci.set_initial_value(-1,0)

    xs = np.linspace(0,2*np.pi)
    ys = np.zeros(xs.shape)
    ys2 = np.zeros(xs.shape)
    t0 = time.time()
    for i, x in enumerate(xs):
        ys[i] = integrator.integrate(x)
    
    t1 = time.time()
    for i, x in enumerate(xs):
        ys2[i] = sci.integrate(x)
    t2 = time.time()
    print(t1 - t0, t2 -t1)
    import matplotlib.pyplot as plt

    plt.plot(xs, ys)
    plt.plot(xs, ys2)
    plt.plot(xs, -np.cos(xs))
    plt.show()

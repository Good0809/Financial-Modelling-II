#importing numpy as scipy.stats for the log and norm.cdf,pdf respectively
import numpy as np
import scipy.stats as ss

'''For all of these they will have inputs: s,k,r,si,t.

s refers to the stock price
k refers to the strike of the option
r refers to the risk free rate
si refers to the implied/assumed volatility
t refers to the tenor (T = 1 means 1yr)

All functions will have a standard call. For example call theta is:

calltheta(s,k,r,si,t)

and put theta

puttheta(s,k,r,si,t)

For all functions the variables are in the same. exact. order.

errors will be returned in any of them are not floats/ints and if si or t is 0.
'''

#creating the "anonymous functions" for d1 and d2 which will help us in our later endeavors.
d1 = lambda s,k,r,si,t: (np.log(s/k) + (r + ((si**2)/2)*t))/(si*np.sqrt(t))
d2 = lambda s,k,r,si,t: d1(s,k,r,si,t) - si*np.sqrt(t) #d2 can call d1 for a simpler calculation.

#call price and put price using the BSM formulas from lecture
callprice = lambda s,k,r,si,t: s*ss.norm.cdf(d1(s,k,r,si,t)) -k*np.exp(-r*t)*ss.norm.cdf(d2(s,k,r,si,t))
putprice = lambda s,k,r,si,t: k*np.exp(-r*t)*ss.norm.cdf(-d2(s,k,r,si,t)) - s*ss.norm.cdf(-d1(s,k,r,si,t))

#once again defining values from lectures and using the functions defined earlier
deltacall = lambda s,k,r,si,t: ss.norm.cdf(d1(s,k,r,si,t))
deltaput = lambda s,k,r,si,t: ss.norm.cdf(d1(s,k,r,si,t)) - 1

#defining gamma and vega
gamma = lambda s,k,r,si,t: ss.norm.pdf(d1(s,k,r,si,t))/(s*si*np.sqrt(t))
vega = lambda s,k,r,si,t: s*ss.norm.pdf(d1(s,k,r,si,t))*np.sqrt(t)

#defining call theta and put theta
calltheta = lambda s,k,r,si,t: -(s*ss.norm.pdf(d1(s,k,r,si,t))*si)/(2*np.sqrt(t)) - r*k*np.exp(-r*t)*ss.norm.cdf(d2(s,k,r,si,t))
puttheta = lambda s,k,r,si,t: -(s*ss.norm.pdf(d1(s,k,r,si,t))*si)/(2*np.sqrt(t)) + r*k*np.exp(-r*t)*ss.norm.cdf(-d2(s,k,r,si,t))

#defining call rho and put rho
callrho = lambda s,k,r,si,t: k*t*np.exp(-r*t)*ss.norm.cdf(d2(s,k,r,si,t))
putrho = lambda s,k,r,si,t: -k*t*np.exp(-r*t)*ss.norm.cdf(-d2(s,k,r,si,t))


#all done, with 12 functions
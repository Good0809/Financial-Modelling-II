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
#importing all the lambda functions for the output following binomial tree
#creating the "anonymous functions" for d1 and d2 which will help us in our later endeavors.
d1 = lambda s,q,k,r,si,t: (np.log(s/k) + ((r-q) + ((si**2)/2)*t))/(si*np.sqrt(t))
d2 = lambda s,q,k,r,si,t: d1(s,q,k,r,si,t) - si*np.sqrt(t) #d2 can call d1 for a simpler calculation.

#call price and put price using the BSM formulas from lecture
callprice = lambda s,q,k,r,si,t: s*ss.norm.cdf(d1(s,q,k,r,si,t)) -k*np.exp(-(r-q)*t)*ss.norm.cdf(d2(s,q,k,r,si,t))
putprice = lambda s,q,k,r,si,t: k*np.exp(-(r-q)*t)*ss.norm.cdf(-d2(s,q,k,r,si,t)) - s*ss.norm.cdf(-d1(s,q,k,r,si,t))

#once again defining values from lectures and using the functions defined earlier
deltacall = lambda s,q,k,r,si,t: ss.norm.cdf(d1(s,q,k,r,si,t))
deltaput = lambda s,q,k,r,si,t: ss.norm.cdf(d1(s,q,k,r,si,t)) - 1

#defining gamma and vega
gamma = lambda s,q,k,r,si,t: ss.norm.pdf(d1(s,q,k,r,si,t))/(s*si*np.sqrt(t))
vega = lambda s,q,k,r,si,t: s*ss.norm.pdf(d1(s,q,k,r,si,t))*np.sqrt(t)

#defining call theta and put theta
calltheta = lambda s,q,k,r,si,t: -(s*ss.norm.pdf(d1(s,q,k,r,si,t))*si)/(2*np.sqrt(t)) - (r-q)*k*np.exp(-(r-q)*t)*ss.norm.cdf(d2(s,q,k,r,si,t))
puttheta = lambda s,q,k,r,si,t: -(s*ss.norm.pdf(d1(s,q,k,r,si,t))*si)/(2*np.sqrt(t)) + (r-q)*k*np.exp(-(r-q)*t)*ss.norm.cdf(-d2(s,q,k,r,si,t))

#defining call rho and put rho
callrho = lambda s,q,k,r,si,t: k*t*np.exp(-(r-q)*t)*ss.norm.cdf(d2(s,q,k,r,si,t))
putrho = lambda s,q,k,r,si,t: -k*t*np.exp(-(r-q)*t)*ss.norm.cdf(-d2(s,q,k,r,si,t))

def testgreek(s,q,k,r,si,t): #this is the function we will call to print the greeks
    print('deltacall is', deltacall(s,q,k,r,si,t))
    print('deltaput is', deltaput(s,q,k,r,si,t))
    
    print('gamma is', gamma(s,q,k,r,si,t))
    print('vega is', vega(s,q,k,r,si,t))
    
    print('call theta is', calltheta(s,q,k,r,si,t))
    print('put theta is', puttheta(s,q,k,r,si,t))
    
    print('call rho is', callrho(s,q,k,r,si,t))
    print('put rho is', putrho(s,q,k,r,si,t))

def stock_sim(s0,s,u,n,l = 0, m=0,):
    if m <= l and l<n+1: #this will iterate through for each row and generate the appropriate stock values. For example, s[0,0] = s0 and s[1,2] = s0*u**(4 - 1), so on and so forth
        s[l,m] = s0*u**(2*m - l) #eventually generating the tree as we know it
        m += 1
        return  stock_sim(s0,s,u,n,l,m)
    elif l < n+1:
        m=0
        l+=1
        return stock_sim(s0,s,u,n,l,m)
    else:
        return s
def step1(s0,k,r,sigma,T,n):
    s = np.zeros([n+1,n+1])
    j = n+1
    dt = T/n
    u = np.exp(sigma*np.sqrt(dt))

    stockprices = stock_sim(s0,s,u,n)

    return stockprices # this just calls the stock prices.... honestly might not be completely necessary but I dont wanna test it now haha
def pricefind(S,k,Acall, Ecall, Aput, Eput, p, r, q,dt,n, l = 0, m = 0):
    if m < n and l < n: #these are gonna find the option prices. There seems to be some rounding error on the american call but they otherwise they all have the appropriate values
        Acall[n - (l+1), m] = max(np.exp(-(r-q)*dt)*((1-p)*Acall[n - l, m] + p*Acall[n-l, m+1]), S[n - (l+1), m] - k)
        Aput[n - (l+1), m] = max(np.exp(-(r-q)*dt)*((1-p)*Aput[n - l, m] + p*Aput[n-l, m+1]), k - S[n - (l+1), m])
        
        Ecall[n - (l+1), m] = max(np.exp(-(r-q)*dt)*((1-p)*Ecall[n - l, m] + p*Ecall[n-l, m+1]), 0)
        Eput[n - (l+1), m] = max(np.exp(-(r-q)*dt)*((1-p)*Eput[n - l, m] + p*Eput[n-l, m+1]), 0)
        m+=1
        pricefind(S,k,Acall,Ecall,Aput,Eput,p,r,q,dt,n,l,m)
    if l < n:
        m = 0
        l +=1 #this is gonna reset and move rows
        return pricefind(S,k,Acall,Ecall,Aput,Eput,p,r,q,dt,n,l,m)
    else:
        return [Acall, Ecall, Aput, Eput]



def finished_product(s0,k,r,q,sigma,T,n):
    j = n+1
    dt = T/n 
    u = np.exp(sigma*np.sqrt(dt))
    d= 1/u
    p = (np.exp((r-q)*dt) -d)/(u-d) #getting all values we need for future functions
    
    stockprices = step1(s0,k,r,sigma,T,n)
    Acall = np.zeros([n+1,n+1])
    Ecall = np.zeros([n+1,n+1]) #making our zero arrays to be used later
    Aput = np.zeros([n+1,n+1])
    Eput = np.zeros([n+1,n+1])
    
    Call = np.triu(np.maximum(stockprices - k,0))
    Put = np.triu(np.maximum(k - stockprices,0)) #our end prices of options, should be the same for both the european and american counterparts
    
    Acall[:,-1] , Ecall[:,-1] = Call[:,-1], Call[:,-1] #assigning aforementioned end prices
    Aput[:,-1], Eput[:,-1] = Put[:,-1], Put[:,-1]
    
    finallist = pricefind(stockprices, k, Acall,Ecall,Aput,Eput,p,r,q,dt,n) #getting our end values for the prices
    
    amcall = finallist[0][0,0] #assigning values to make printing easier down below
    eucall = finallist[1][0,0]
    amput = finallist[2][0,0]
    euput = finallist[3][0,0]

    print('The American Call is ', amcall) 
    print('The European Call is ', eucall)
    print('The American Put is ', amput)
    print('The European Put is ', euput)
    
    testgreek(s0,q,k,r,sigma,T) #printing all the greeks
    return

finished_product(50,50,0.05,0,0.5,1,6)

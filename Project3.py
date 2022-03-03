#implied vol finder
import numpy as np
import scipy.stats as ss
from yahoo_fin import options
from scipy.optimize import fmin
import pandas as pd
import matplotlib.pyplot as plt
#copy over necessary greeks

d1 = lambda s,k,r,si,t: (np.log(s/k) + (r + ((si**2)/2)*t))/(si*np.sqrt(t))
d2 = lambda s,k,r,si,t: d1(s,k,r,si,t) - si*np.sqrt(t) #d2 can call d1 for a simpler calculation.
#call price and put price using the BSM formulas from lecture
callprice = lambda s,k,r,si,t: s*ss.norm.cdf(d1(s,k,r,si,t)) -k*np.exp(-r*t)*ss.norm.cdf(d2(s,k,r,si,t))
putprice = lambda s,k,r,si,t: k*np.exp(-r*t)*ss.norm.cdf(-d2(s,k,r,si,t)) - s*ss.norm.cdf(-d1(s,k,r,si,t))
vega = lambda s,k,r,si,t: s*ss.norm.pdf(d1(s,k,r,si,t))*np.sqrt(t)


#user should be able to provide:

#underlying price (S), strike price (k), risk-free rate (r), tenor (T), option price (op), and call/put type (boolean)
#newton method, bisection

def bisection(rng,s,k,r,T,op,Call,voldiff = [0,100], hf = 0): #range is going to be an np.linspace array
    #print(voldiff)
    if hf > 21:
        return voldiff[0]
    for x in range(len(rng)):
        if Call == True:
            #print(voldiff, 1)
            deltaxn = op - callprice(s,k,r,rng[x],T)
            if deltaxn < 0:
                deltaxn = deltaxn*-1
            #print(deltaxn)
            if voldiff[1] < deltaxn and np.abs(voldiff[1]) > 0.001: #doesnt work great for extremely small volatilities
                hf+=1
                newrng = np.linspace(rng[x - 1],rng[x],100)
                return bisection(newrng,s,k,r,T,op,Call,voldiff = voldiff,hf=hf)
            elif voldiff[1] < deltaxn and voldiff[1] <= 0.001:
                return voldiff[0]
            else:
                voldiff[0] = rng[x]
                voldiff[1] = deltaxn
        else:
            #print(voldiff, 2)
            deltaxn = op - putprice(s,k,r,rng[x],T)
            if deltaxn < 0:
                deltaxn = deltaxn*-1
            if voldiff[1] < deltaxn and np.abs(voldiff[1]) > 0.001:
                hf+=1
                newrng = np.linspace(rng[x-1],rng[x],100)
                return bisection(newrng,s,k,r,T,op,Call,voldiff = voldiff,hf=hf)
            elif voldiff[1] < deltaxn and voldiff[1] <= 0.001:
                return voldiff[0]
            else:
                voldiff[0] = rng[x]
                voldiff[1] = deltaxn
    hf+=1
    newrng = np.linspace(rng[-1], rng[-1]+2,100) #if nothing triggers above then this needs to occur cause deltaxn still isnt within 0.01%
    return bisection(newrng, s,k,r,T,op,Call, voldiff = voldiff)
                
def volfinder(s,k,r,T,op,Call):
    deltaxn = 10
    attvol = 0.5 #attempted vol, initial guess is 0.5
    while deltaxn > 0.00001: #fine with being  such a small number off
        if Call == True:
            deltaxn = (op - callprice(s,k,r,attvol, T)) / vega(s,k,r,attvol,T)
            attvol += deltaxn
            #print(deltaxn) #logic checkers
        else:
            deltaxn = (op - putprice(s,k,r,attvol, T)) / vega(s,k,r,attvol,T)
            attvol += deltaxn
            #print(deltaxn)
    newtonvol = attvol #will giv vol from Newtons Approx atm
    rang = np.linspace(0,2,100) #initial guess is between 0 and 200%, though I guess I could be cheeky and use a range with newtons approx in it haha
    #print(newtonvol)
    if Call == False:
        bivol = bisection(rang,s,k,r,T,op,False)
    else:
        bivol = bisection(rang,s,k,r,T,op,True)
    
    return newtonvol, bivol
#now we need to do the Gatheral SVI Skew Fit
from yahoo_fin import options
from scipy.optimize import fmin
import pandas as pd
import matplotlib.pyplot as plt
def GathSVI(tick,options, call): #takes an array of options prices given some put, then the usual stuff with Call once again being a boolean 
    def func(x):
        diff_list = np.array([])
        for y in range(options.shape[0]):
            if call == True:
                siv = lambda x: x[0] + np.abs(x[1])*(x[2]*(options['Call Strikes'][y] - x[3]) + np.sqrt((options['Call Strikes'][y] - x[3])**2 + x[4]**2))
                
                diffsq = (siv(x) - options['Call Vols'][y])**2
                diff_list = np.append(diff_list, diffsq)
            else:
                siv = lambda x: x[0] + np.abs(x[1])*(x[2]*(options['Put Strikes'][y] - x[3]) + np.sqrt((options['Put Strikes'][y] - x[3])**2 + x[4]**2))
                
                diffsq = (siv(x) - options['Put Vols'][y])**2
                diff_list = np.append(diff_list, diffsq)
        return diff_list.sum()
    if call == True:
        variables = fmin(func, [options['Call Vols'].median(),0.02,-0.04,options['Call Strikes'].median(),0.3])
    else:
        variables = fmin(func, [options['Put Vols'].median(),0.02,-0.5,options['Put Strikes'].median(),0.3])
    
    return variables #returns in order a,b,p,m,sigma

def getoptions(tick):
    net_options = options.get_options_chain(tick)
    #print(net_options)
    #net_options['calls']
    df = pd.DataFrame()
    df['Call Vols'] = (pd.to_numeric(net_options['calls']['Implied Volatility'].str[:-1]).div(100).mask(net_options['calls']['Implied Volatility'] == '%', 0)) #changing into numbers
    df['Put Vols'] = (pd.to_numeric(net_options['puts']['Implied Volatility'].str[:-1]).div(100).mask(net_options['puts']['Implied Volatility'] == '%', 0))
    df['Call Strikes'] = net_options['calls']['Strike']
    df['Put Strikes'] = net_options['puts']['Strike']
    
    return df


def finalfunc(s,k,r,T,op,call,ticker): #this is gonna test everything
    vols = volfinder(s,k,r,T,op,call)
    newtonvol, bivol = vols[0],vols[1]
    
    print('vol using newton for given option price is', newtonvol)
    print('vol using bisection for given option price is', bivol)
    
    options = getoptions(ticker)
    
    variables = GathSVI(ticker,options, call)
    
    print('Using the ticker and the options found from yahoo_fin the best parameters are:')
    
    print('a = ', variables[0])
    print('b = ', np.abs(variables[1]))
    print('rho = ', variables[2])
    print('m = ', variables[3])
    print('sigma = ', np.abs(variables[4]))
    
    print('Now the plot')
    hope = lambda x: variables[0] + np.abs(variables[1])*(variables[2]*(x - x[3]) + np.sqrt((x - variables[3])**2 + variables[4]**2))
    if call == True:
        plt.figure(figsize = (15,10))
        plt.title('Options from ' + ticker + ' with imp vol versus Gatherall SVI Fit')
        plt.scatter(options['Call Strikes'], options['Call Vols'], label ='imp vols')
        plt.scatter(options['Call Strikes'], hope(options['Call Strikes']), label = 'gatherall svi')
        plt.legend()
    else:
        plt.figure(figsize = (15,10))
        plt.title('Options from ' + ticker + ' with imp vol versus Gatherall SVI Fit')
        plt.scatter(options['Put Strikes'], options['Put Vols'], label ='imp vols')
        #x = np.linspace(200,900,options.shape[0])
        plt.scatter(options['Put Strikes'], hope(options['Put Strikes']), label = 'gatherall svi')
        #plt.scatter(x, hope(x), label = 'gatherall svi')
        plt.legend()
        
pprice = putprice(100,90,0.02,0.5,1)
cprice = callprice(100,90,0.02,0.5,1)
finalfunc(100,90,0.02,1,cprice,True,'nflx') #it can take a bit because of fmin sadly
#implied vol finder
import numpy as np
import scipy.stats as ss
import yfinance as yf
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

def bisection(rng,s,k,r,T,op,Call,voldiff = [0,100], hf = 0): #range is going to be an np.linspace array. voldiff has a high setting
    #initially as we need some initials and want something less than 100 haha

    #print(voldiff) #logic check
    if hf > 30: #only want it to iterate through so many times or else it's only going to go nowhere fast. As such after
        #30 runs either the volatility is astronomically high or we're so close it's looping
        #this is also to save processing power
        return voldiff[0]
    for x in range(len(rng)):
        if Call == True:
            #print(voldiff, 1)
            deltaxn = op - callprice(s,k,r,rng[x],T)
            if deltaxn < 0:
                deltaxn = deltaxn*-1 #need absolute deltaxn for the if statements
                #for some reason np.abs wasn't working great and would give the occasional error
                #which messed up the entire process

            #print(deltaxn)
            if voldiff[1] < deltaxn and np.abs(voldiff[1]) > 0.001: #doesnt work great for extremely small volatilities
                hf+=1 #add a count here since we're restarting basically
                newrng = np.linspace(rng[x - 1],rng[x],1000) #since the new deltaxn is larger than the last
                #this implies it's somwhere between the old and new
                #as such we refit our range to a better approximation
                return bisection(newrng,s,k,r,T,op,Call,voldiff = voldiff,hf=hf)
            elif voldiff[1] < deltaxn and voldiff[1] <= 0.001:
                return voldiff[0] #return our impvol
            else:
                voldiff[0] = rng[x] #rng[x] is now our best imp_vol bet
                voldiff[1] = deltaxn #check the distance since we'll use this for our previous if statements
        else: #just repeating everything but for puts

            #print(voldiff, 2) #logic check
            deltaxn = op - putprice(s,k,r,rng[x],T)
            if deltaxn < 0:
                deltaxn = deltaxn*-1
            if voldiff[1] < deltaxn and np.abs(voldiff[1]) > 0.001:
                hf+=1
                newrng = np.linspace(rng[x-1],rng[x],1000)
                return bisection(newrng,s,k,r,T,op,Call,voldiff = voldiff,hf=hf)
            elif voldiff[1] < deltaxn and voldiff[1] <= 0.001:
                return voldiff[0]
            else:
                voldiff[0] = rng[x]
                voldiff[1] = deltaxn
    hf+=1 #once again increasing since we're using incursion
    newrng = np.linspace(rng[-1], rng[-1]+2,1000) #if nothing triggers above then this
    #needs to occur cause deltaxn still isnt within 0.01%
    return bisection(newrng, s,k,r,T,op,Call, voldiff = voldiff) #recurring with new range (200% increase)
                
def volfinder(s,k,r,T,op,Call):
    deltaxn = 10
    attvol = 0.5 #attempted vol, initial guess is 0.5
    while deltaxn > 0.00001: #This is newtons method
        #ultimately fine with being  such a small number off
        if Call == True:
            deltaxn = (op - callprice(s,k,r,attvol, T)) / vega(s,k,r,attvol,T)
            attvol += deltaxn
            #print(deltaxn) #logic checkers
        else:
            deltaxn = (op - putprice(s,k,r,attvol, T)) / vega(s,k,r,attvol,T)
            attvol += deltaxn
            #print(deltaxn)
    newtonvol = attvol #will giv vol from Newtons Approx atm
    rang = np.linspace(0,2,1000) #initial guess is between 0 and 200%,
    #though I guess I could be cheeky and use a range with newtons approx in it haha
    #though that'd kinda defeat the purpose of this one by itself


    #print(newtonvol) #logic check
    if Call == False:
        bivol = bisection(rang,s,k,r,T,op,False) #refer to bisection
    else:
        bivol = bisection(rang,s,k,r,T,op,True)
    
    return newtonvol, bivol

#now we need to do the Gatheral SVI Skew Fit
def GathSVI(tick,options, call): #takes an array of options prices given some put, then the usual stuff with Call once again being a boolean 
    def func(x):
        diff_list = np.array([])#our diff_list, this will be used at the end to
        #sum our squared differences

        for y in range(options.shape[0]):
            if call == True: #call process

                siv = lambda x: x[0] + np.abs(x[1])*(x[2]*(options['Call Strikes'][y] - x[3]) + np.sqrt((options['Call Strikes'][y] - x[3])**2 + x[4]**2))
                #our siv function for calls, with puts we simply replace the options columns
                diffsq = (siv(x) - options['Call Vols'][y])**2 #our squared error

                diff_list = np.append(diff_list, diffsq) #append it to the array
            else:
                #repeat notes from call
                siv = lambda x: x[0] + np.abs(x[1])*(x[2]*(options['Put Strikes'][y] - x[3]) + np.sqrt((options['Put Strikes'][y] - x[3])**2 + x[4]**2))
                
                diffsq = (siv(x) - options['Put Vols'][y])**2
                diff_list = np.append(diff_list, diffsq)
        return diff_list.sum() #returning the sum, we're trying to minimize this value
    if call == True:
        maxstrike = options['Call Strikes'].max() #getting the max and minstrikes for our if statements
        #see more notes inside statement as it relates to what the values mean
        minstrike = options['Call Strikes'].min()
        volofmax = options[options['Call Strikes'] == maxstrike]['Call Vols'].max()
        #print(volofmax) #logic check
        volofmin = options[options['Call Strikes'] == minstrike]['Call Vols'].max()
        minvol = options['Call Vols'].min()
        if minvol < volofmin and minvol < volofmax: #indicates non-linear vol curve
            variables = fmin(func, [options['Call Vols'].min(),0.5,-0.5,options['Call Strikes'].min(),1])
        else:
            variables = fmin(func, [(options['Call Vols'].max() +options['Call Vols'].min())/2,0.0,-0.5,options['Call Strikes'].median(),1])
        #our initial guesses are based off the observations I noticed on Desmos after a conversation with Eric in class

        #alpha changes how high up the distribution starts. As such the best guess is the minimum implied volatility
        #adding one since sometimes it likes to drop low

        #beta changes if the curve is up or down (down if neg, up if positive). We know beta is always positive so it's
        #either higher positive or 0 in the case of linearity.

        #rho changes curvature of the left side. Good bet usually is negative as we need that side to be decreasing and not linear
        #as such -0.5 should be good as if theres an inverted vol curve thats a bigger problem

        #m changes where the critical point is. So, best bet is always where the min vols strike is in the case of curvature.
        # Otherwise the median point is pretty good

        #finally sigma flattens the curve the higher we get. As such I want the bet to be low initially since there's most likely
        #going to be a curve no matter how deep. Somewhere around 1 looks good
    else:
        #see notes in above if loop
        variables = fmin(func, [options['Put Vols'].min(),0.02,-0,options['Put Strikes'].median(),0.3])
        maxstrike = options['Put Strikes'].max()
        minstrike = options['Put Strikes'].min()
        volofmax = options[options['Put Strikes'] == maxstrike]['Put Vols'].max()
        volofmin = options[options['Put Strikes'] == minstrike]['Put Vols'].max()
        minvol = options['Put Vols'].min()
        if minvol < volofmin and minvol < volofmax: #indicates non-linear vol curve
            variables = fmin(func, [options['Put Vols'].min(),0.5,-0.5,options['Put Strikes'].min(),1])
        else:
            variables = fmin(func, [(options['Put Vols'].max() +options['Put Vols'].min())/2,0,-0.5,options['Put Strikes'].median(),1])
    return variables #returns in order a,b,p,m,sigma

def getoptions(tick,date=np.nan): #just gonna get our options using yfinance's option_chain. Including an optional date parameter.
    tickhist = yf.Ticker(tick)
    if type(date) == type(np.nan): 
        net_options = tickhist.option_chain()
    else: #if date is invalid this will throw an error
        net_options = tickhist.option_chain(date)
    #print(net_options)
    #net_options['calls']
    df = pd.DataFrame()
    df['Call Vols'] = net_options.calls.impliedVolatility #get everything going
    df['Put Vols'] = net_options.puts.impliedVolatility
    df['Call Strikes'] = net_options.calls.strike
    df['Put Strikes'] = net_options.puts.strike
    
    return df #return the dataframe


def finalfunc(s,k,r,T,op,call,ticker,date = np.nan): #this is gonna test everything we've done so far
    #also added the optional date parameter here.
    vols = volfinder(s,k,r,T,op,call) #get our newtonian and bisectional volatilities
    newtonvol, bivol = vols[0],vols[1] #define them here
    
    print('vol using newton for given option price is', newtonvol) #tell the use what they are
    print('vol using bisection for given option price is', bivol)
    
    options = getoptions(ticker,date) #getting our dataframe nice and ready
    
    variables = GathSVI(ticker,options, call) #finding the paramaters for the Gatheral SVI Skew Fit
    
    print('Using the ticker and the options found from yfinance the best parameters are:')
    
    print('a = ', variables[0]) #tell the user the parameters
    print('b = ', np.abs(variables[1]))
    print('rho = ', variables[2])
    print('m = ', variables[3])
    print('sigma = ', np.abs(variables[4]))
    
    print('Now the plot')
    #define the gatheral function
    #called hope cause I had hoped this would work. It does. I am happy. Transcendant even.
    hope = lambda x: variables[0] + np.abs(variables[1])*(variables[2]*(x - x[3]) + np.sqrt((x - variables[3])**2 + variables[4]**2))
    if call == True:
        plt.figure(figsize = (15,10)) #just plotting it out with scatters
        plt.title('Call Options from ' + ticker + ' with imp vol versus Gatherall SVI Fit')
        plt.scatter(options['Call Strikes'], options['Call Vols'], label ='imp vols')
        plt.scatter(options['Call Strikes'], hope(options['Call Strikes']), label = 'Gatheral SVI Skew Fit')
        plt.legend()
        plt.show()
    else:
        plt.figure(figsize = (15,10))#once again just plotting it out
        plt.title('Put Options from ' + ticker + ' with imp vol versus Gatherall SVI Fit')
        plt.scatter(options['Put Strikes'], options['Put Vols'], label ='imp vols')
        plt.scatter(options['Put Strikes'], hope(options['Put Strikes']), label = 'Gatheral SVI Skew Fit')
        plt.legend()
        plt.show()


#modify the following values if you please
StockPrice = 100 #underlying price
StrikePrice = 90
RF = 0.02 #risk free rate
Tenor = 1
Sigma = 0.5
pprice = putprice(StockPrice,StrikePrice,RF,Sigma,Tenor)
cprice = callprice(StockPrice,StrikePrice,RF,Sigma,Tenor)
theticker = 'nflx'
date = np.nan
def runit(call): #this will be the function that shows two examples, one of a put and one of a call
    if call == False:
        return finalfunc(StockPrice,StrikePrice,RF,Tenor,pprice,False,theticker,date) #it can take a bit because of fmin sadly
    else:
        return finalfunc(StockPrice,StrikePrice,RF,Tenor,cprice,True,theticker,date)

print(runit(True)) #change this to True of False. True is Call, False is Put.
#I've noticed that sometimes when you run this program at night you get some really, really weird options from yfinance so be forewarned
#best results I've gotten from yfinance were around 12-4

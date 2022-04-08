using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Project5
{
    class Program
    {
        static void Main(string[] args)
        {
            Boolean done = false;
            double initp = 100d;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the price of your stock?");
                    initp = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    initp = 100d;
                    Console.WriteLine("Invalid entry, input must be a number!");
                }
            }

            done = false;
            double stockvol = 0.2;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the vol of your stock?");
                    stockvol = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    stockvol = 0.2;
                    Console.WriteLine("Invalid entry, input must be a number!");
                }
            }

            done = false;
            double riskfree = 0.02;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the risk-free rate of your stock?");
                    riskfree = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    riskfree = 0.02;
                    Console.WriteLine("Invalid entry, input must be a number!");
                }
            }

            done = false;
            double Tenor = 1d;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the Tenor your option?");
                    Tenor = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    Tenor = 1d;
                    Console.WriteLine("Invalid entry, input must be a number!");
                }
            }

            done = false;
            int steps = 252;
            while (!done)
            {
                try
                {
                    Console.WriteLine("How many steps for your option? ");
                    steps = int.Parse(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    steps = 252;
                    Console.WriteLine("Invalid entry, input must be an integernumber!");
                }
            }

            done = false;
            int paths = 100;
            while (!done)
            {
                try
                {
                    Console.WriteLine("How many paths for your option?");
                    paths = int.Parse(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    paths = 100;
                    Console.WriteLine("Invalid entry, input must be an integer number!");
                }
            }
            done = false;
            double strike = 100;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the strike for your option?");
                    strike = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    strike = 100;
                    Console.WriteLine("Invalid entry, input must be an integer number!");
                }
            }

            done = false;
            Boolean call = true;
            while (!done)
            {
                try
                {
                    Console.WriteLine("Is your option a call or a put? (0 = put, 1 = call)");
                    int input = int.Parse(Console.ReadLine());
                    if (input == 1 || input == 0)
                    {
                        if (input == 1)
                        {
                            call = true;
                        }
                        else
                        {
                            call = false;
                        }
                        done = true;
                    }
                    else
                    {
                        Console.WriteLine("Input must be 0 or 1. Try again.");
                    }
                }
                catch
                {
                    strike = 100;
                    Console.WriteLine("Invalid entry, input must be an integer number 0 or 1!");
                }
            }

            Stock stock = new Stock() {Value = initp, Vol = stockvol};

            List<double[,]> stockrandsims = stock.StockSim(riskfree, Tenor, steps, paths);

            EuropeanOption option = new EuropeanOption() {Underlying = stock, Tenor = Tenor, Strike = strike, Call = call};

            double op_price = option.PriceViaMonte(stockrandsims[0], riskfree);
            double[,] step_payoffs = option.steppayoff(stockrandsims[0],riskfree);

            //Console.WriteLine(step_payoffs.GetLength(0));
            //Console.WriteLine(step_payoffs.GetLength(1));

            double StdErr = option.StandardError(step_payoffs);
            double optionDelta = option.optionsDelta(stockrandsims[1], riskfree);
            double optionGamma = option.optionsGamma(stockrandsims[1], riskfree);
            double optionTheta = option.optionsTheta(step_payoffs);
            double optionVega = option.optionsVega(stockrandsims[1], stockrandsims[0], riskfree);
            double optionRho = option.optionsRho(stockrandsims[1], stockrandsims[0], riskfree);


            Console.WriteLine("The option price via Monte Carlo is {0}", op_price);
            Console.WriteLine();
            Console.WriteLine("The StdError of the run for payoffs was {0}", StdErr);
            Console.WriteLine();
            Console.WriteLine("The option Delta is {0}", optionDelta);
            Console.WriteLine();
            Console.WriteLine("The option Gamma is {0}", optionGamma);
            Console.WriteLine();
            Console.WriteLine("The option Theta is {0}", optionTheta);
            Console.WriteLine();
            Console.WriteLine("The option Vega is {0}", optionVega);
            Console.WriteLine();
            Console.WriteLine("The option Rho is {0}", optionRho);

        }



    }
    class Gaussian // need a gaussian, might as well make a class that can be modified later for other methods
    {
        public static List<double> BoxMuller(Random x1, Random x2)
        {
            double y1, y2;
            y1 = x1.NextDouble();
            y2 = x2.NextDouble(); //get our uniform numbers between 0 and 1

            List<double> gaussians = new List<double>(); // our output list of gaussians

            double z1,z2;

            z1 = Math.Sqrt(-2*Math.Log(y1)) *Math.Cos(2*Math.PI*y2); // generate the gaussian using the formulas from slide 52 on lec 8
            z2 = Math.Sqrt(-2*Math.Log(y1)) *Math.Sin(2*Math.PI*y2);

            gaussians.Add(z1); // add them to the list
            gaussians.Add(z2);

            return gaussians; // return the list


        }
    }
    class Stock
    {
        public double Value {get;set;}
        public double Vol {get;set;} // every stock has a volatility. As such the volatility of the options will just be that of the stock
        
        public List<double[,]> StockSim(double r, double T, int steps, int paths) // dont want static since it uses the stock characteristics
        {
            double dt = T/Convert.ToDouble(steps);
            double[,] PricePaths = new double[paths,steps+1]; // each row should be a path and each column a step
            double[,] PathRandoms = new double[paths,steps+1];

            for (int x = 0; x < steps+1;x++)
            {
                for (int y = 0; y < paths; y++)
                {
                    if (x == 0)
                    {
                        PricePaths[y,x] = Value; // each initial step of each path is value and their random is 0
                        PathRandoms[y,x] = 0;
                    }
                    else
                    {
                        Random rand1 = new Random();
                        Random rand2 = new Random();
                        double Gauss = Gaussian.BoxMuller(rand1, rand2)[0];
                        PathRandoms[y,x] = Gauss;

                        //Console.WriteLine("Random from path {0}, step {1} is {2}", y, x, PathRandoms[y,x]); // confirm each gaussian is unique

                        PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2.0)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                    }
                }
            }

            List<double[,]> output = new List<double[,]>();

            output.Add(PricePaths);
            output.Add(PathRandoms);
            return output;
        }

    }
    class EuropeanOption
    {
        public Stock Underlying {get;set;} // this gets us the initial price and volatility
        public double Tenor {get;set;} // need the tenor of our option

        public double Strike {get;set;} // our strike

        public Boolean Call {get;set;} // is it a call or put?

        public double PriceViaMonte(double[,] simmedstock, double r) // r represents the risk free rate
        {
            int paths = simmedstock.GetLength(0);
            int steps = simmedstock.GetLength(1);
            List<double> payoffs = new List<double>();
            //int count = 0;
            for (int x = 0; x < paths; x++)
            {
                if(Call)
                {
                    payoffs.Add(Math.Max(simmedstock[x, steps-1] - Strike, 0d)*Math.Exp(-r * Tenor)); // last column is the last step. So we want to iterate through rows
                    //count++;
                    //Console.WriteLine(count);
                }
                else
                {
                    payoffs.Add(Math.Max(Strike - simmedstock[x,steps-1],0d)*Math.Exp(-r * Tenor));
                }
            }
            return payoffs.Average();
        }

        public double[,] steppayoff(double[,] simmedstock, double r) //array of payoffs which will help us with calculating Greeks
        {
            // r represents the risk-free rate
            int paths = simmedstock.GetLength(0);
            int steps = simmedstock.GetLength(1);
            double[,] oppayoffs = new double[paths, steps]; //discounted payoffs
            double dt = Tenor/steps;
            for (int x = 0; x < paths; x++)
            {
                for (int y = 0; y<steps; y++)
                {
                    if (Call)
                    {
                        oppayoffs[x,y] = Math.Max(simmedstock[x,y] - Strike, 0)*Math.Exp(-r*dt*y);
                    }
                    else
                    {
                        oppayoffs[x,y] = Math.Max(Strike - simmedstock[x,y], 0)*Math.Exp(-r*dt*y);
                    }

                }
            }
            return oppayoffs; //return the discounted payoffs
        }
        // the following Greeks formulas are from Lecture 4 slides 45,46 as well as lecture 10 slides 31-33
        public double optionsDelta(double[,] randoms, double r)
        {
            int steps = randoms.GetLength(1);
            int paths = randoms.GetLength(0);
            double dt = Tenor/Convert.ToDouble(steps);
            //double[,] deltas = new double[paths, steps-1];
            List<double> deltas = new List<double>();
            double v = Underlying.Vol;
            double changeS = (1.0/1000.0)*Underlying.Value;
            //int count = 0;
            double[,] upperstock = new double[paths,steps];
            double[,] lowerstock = new double[paths,steps];
            for (int x = 0; x < steps;x++)
            {
                for (int y = 0; y < paths; y++)
                {
                    if (x == 0)
                    {
                        //Console.WriteLine(0);
                        upperstock[y,x] = Underlying.Value + changeS; // each initial step of each path is value and their random is 0
                        lowerstock[y,x] = Underlying.Value - changeS;
                    }
                    else
                    {
                        //count+=1;
                        //Console.WriteLine(count);

                        upperstock[y,x] = (upperstock[y,x-1])*Math.Exp((r - (Math.Pow(v, 2))/2.0)*dt + (Underlying.Vol)*Math.Sqrt(dt)*randoms[y,x]);
                        lowerstock[y,x] = (lowerstock[y,x-1])*Math.Exp((r - (Math.Pow(v, 2))/2.0)*dt + (Underlying.Vol)*Math.Sqrt(dt)*randoms[y,x]);
                    }
                }
            }
            double[,] upperpayoffs = new double[paths,steps];
            double[,] lowerpayoffs = new double[paths,steps];
            for (int x = 0; x < paths; x++)
            {
                for (int y = 0; y<steps; y++)
                {
                    if (Call)
                    {
                        upperpayoffs[x,y] = Math.Max(upperstock[x,y] - Strike, 0)*Math.Exp(-(r)*dt*y);
                        lowerpayoffs[x,y] = Math.Max(lowerstock[x,y] - Strike,0)*Math.Exp(-(r)*dt*y);
                    }
                    else
                    {
                        upperpayoffs[x,y] = Math.Max(Strike - upperstock[x,y], 0)*Math.Exp(-(r)*dt*y);
                        lowerpayoffs[x,y] = Math.Max(Strike - lowerstock[x,y],0)*Math.Exp(-(r)*dt*y);
                    }

                }
            }
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 1; y <steps; y++)
                {
                    double k = upperpayoffs[x,y] -lowerpayoffs[x,y];
                    double j = 2.0*changeS;
                    //Console.WriteLine(k);
                    deltas.Add(k/j);
                    //Console.WriteLine(vegalist.Count());
                    //Console.WriteLine(vegalist[count]);
                    //count +=1;

                }
            }
            double avgdelta = deltas.Average();

            return avgdelta;
        }

        public double optionsGamma(double[,] randoms, double r)
        {
            int steps = randoms.GetLength(1);
            int paths = randoms.GetLength(0);
            //double[,] gammas = new double[paths, steps-2];
            double dt = Tenor/Convert.ToDouble(steps);
            //double[,] deltas = new double[paths, steps-1];
            List<double> gammas = new List<double>();
            double v = Underlying.Vol;
            double changeS = (1.0/100.0)*Underlying.Value;
            //int count = 0;
            double[,] upperstock = new double[paths,steps];
            double[,] regstock = new double[paths,steps];
            double[,] lowerstock = new double[paths,steps];
            for (int x = 0; x < steps;x++)
            {
                for (int y = 0; y < paths; y++)
                {
                    if (x == 0)
                    {
                        //Console.WriteLine(0);
                        upperstock[y,x] = Underlying.Value + changeS;
                        regstock[y,x] = Underlying.Value;
                        lowerstock[y,x] = Underlying.Value - changeS;
                    }
                    else
                    {
                        //count+=1;
                        //Console.WriteLine(count);

                        upperstock[y,x] = (upperstock[y,x-1])*Math.Exp((r - (Math.Pow(v, 2))/2.0)*dt + (Underlying.Vol)*Math.Sqrt(dt)*randoms[y,x]);
                        regstock[y,x] = (regstock[y,x-1])*Math.Exp((r - (Math.Pow(v, 2))/2.0)*dt + (Underlying.Vol)*Math.Sqrt(dt)*randoms[y,x]);
                        lowerstock[y,x] = (lowerstock[y,x-1])*Math.Exp((r - (Math.Pow(v, 2))/2.0)*dt + (Underlying.Vol)*Math.Sqrt(dt)*randoms[y,x]);
                    }
                }
            }
            double[,] upperpayoffs = new double[paths,steps];
            double[,] regpayoffs = new double[paths,steps];
            double[,] lowerpayoffs = new double[paths,steps];
            for (int x = 0; x < paths; x++)
            {
                for (int y = 0; y<steps; y++)
                {
                    if (Call)
                    {
                        upperpayoffs[x,y] = Math.Max(upperstock[x,y] - Strike, 0)*Math.Exp(-(r)*dt*y);
                        regpayoffs[x,y] = Math.Max(regstock[x,y] - Strike, 0)*Math.Exp(-(r)*dt*y);
                        lowerpayoffs[x,y] = Math.Max(lowerstock[x,y] - Strike,0)*Math.Exp(-(r)*dt*y);
                    }
                    else
                    {
                        upperpayoffs[x,y] = Math.Max(Strike - upperstock[x,y], 0)*Math.Exp(-(r)*dt*y);
                        regpayoffs[x,y] = Math.Max(Strike - regstock[x,y], 0)*Math.Exp(-(r)*dt*y);
                        lowerpayoffs[x,y] = Math.Max(Strike - lowerstock[x,y],0)*Math.Exp(-(r)*dt*y);
                    }

                }
            }
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 1; y <steps; y++)
                {
                    double k = upperpayoffs[x,y] - 2.0*regpayoffs[x,y] + lowerpayoffs[x,y];
                    double j = Math.Pow(changeS,2);
                    //Console.WriteLine(k);
                    gammas.Add(k/j);
                    //Console.WriteLine(vegalist.Count());
                    //Console.WriteLine(vegalist[count]);
                    //count +=1;

                }
            }
            double avggamma = gammas.Average();

            return avggamma;
        }

        public double optionsTheta(double[,] opvals)
        {
            int steps = opvals.GetLength(1);
            int paths = opvals.GetLength(0);
            double dsteps = Convert.ToDouble(steps);
            double dt = Tenor/dsteps;
            double v = Underlying.Vol;
            List<double> thetalist = new List<double>();
            //int count = 0;
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 0; y <steps-1; y++)
                {
                    thetalist.Add((opvals[x,y+1] - opvals[x,y])/(dt));
                    //Console.WriteLine(vegalist.Count());
                    //Console.WriteLine(vegalist[count]);
                    //count +=1;

                }
            }
            double averagetheta = thetalist.Average();

            return averagetheta*-1.0;
        }

        public double optionsVega(double[,] randoms, double[,] stockvals, double r)
        {
            int steps = randoms.GetLength(1);
            int paths = randoms.GetLength(0);
            double dt = Tenor/steps;
            List<double> vegalist = new List<double>();
            double v = Underlying.Vol;
            double deltavol = (1/100.0)*v;
            //int count = 0;
            double[,] upperstock = new double[paths,steps];
            double[,] lowerstock = new double[paths,steps];
            for (int x = 0; x < steps;x++)
            {
                for (int y = 0; y < paths; y++)
                {
                    if (x == 0)
                    {
                        //Console.WriteLine(0);
                        upperstock[y,x] = Underlying.Value; // each initial step of each path is value and their random is 0
                        lowerstock[y,x] = Underlying.Value;
                    }
                    else
                    {
                        //count+=1;
                        //Console.WriteLine(count);

                        upperstock[y,x] = upperstock[y,x-1]*Math.Exp((r - (Math.Pow(v + deltavol, 2))/2.0)*dt + (Underlying.Vol+deltavol)*Math.Sqrt(dt)*randoms[y,x]);
                        lowerstock[y,x] = lowerstock[y,x-1]*Math.Exp((r - (Math.Pow(v - deltavol, 2))/2.0)*dt + (Underlying.Vol-deltavol)*Math.Sqrt(dt)*randoms[y,x]);
                    }
                }
            }
            double[,] upperpayoffs = new double[paths,steps];
            double[,] lowerpayoffs = new double[paths,steps];
            for (int x = 0; x < paths; x++)
            {
                for (int y = 0; y<steps; y++)
                {
                    if (Call)
                    {
                        upperpayoffs[x,y] = Math.Max(upperstock[x,y] - Strike, 0)*Math.Exp(-r*dt*y);
                        lowerpayoffs[x,y] = Math.Max(lowerstock[x,y] - Strike,0)*Math.Exp(-r*dt*y);
                    }
                    else
                    {
                        upperpayoffs[x,y] = Math.Max(Strike - upperstock[x,y], 0)*Math.Exp(-r*dt*y);
                        lowerpayoffs[x,y] = Math.Max(Strike - lowerstock[x,y],0)*Math.Exp(-r*dt*y);
                    }

                }
            }
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 0; y <steps; y++)
                {
                    double k = upperpayoffs[x,y] -lowerpayoffs[x,y];
                    double j = 2.0*deltavol;
                    //Console.WriteLine(k);
                    vegalist.Add(k/j);
                    //Console.WriteLine(vegalist.Count());
                    //Console.WriteLine(vegalist[count]);
                    //count +=1;

                }
            }
            double averagevega = vegalist.Average();

            return averagevega;
        }

        public double optionsRho(double[,] randoms, double[,] stockvals, double r)
        {
            int steps = stockvals.GetLength(1);
            int paths = stockvals.GetLength(0);
            double dt = Tenor/steps;
            List<double> rholist = new List<double>();
            double v = Underlying.Vol;
            double deltar = (1/100.0)*r;
            //int count = 0;
            double[,] upperstock = new double[paths,steps];
            double[,] lowerstock = new double[paths,steps];
            for (int x = 0; x < steps;x++)
            {
                for (int y = 0; y < paths; y++)
                {
                    if (x == 0)
                    {
                        //Console.WriteLine(0);
                        upperstock[y,x] = Underlying.Value; // each initial step of each path is value and their random is 0
                        lowerstock[y,x] = Underlying.Value;
                    }
                    else
                    {
                        //count+=1;
                        //Console.WriteLine(count);

                        upperstock[y,x] = upperstock[y,x-1]*Math.Exp((r+deltar - (Math.Pow(Underlying.Vol, 2))/2.0)*dt + Underlying.Vol*Math.Sqrt(dt)*randoms[y,x]);
                        lowerstock[y,x] = lowerstock[y,x-1]*Math.Exp((r-deltar - (Math.Pow(Underlying.Vol, 2))/2.0)*dt + Underlying.Vol*Math.Sqrt(dt)*randoms[y,x]);
                    }
                }
            }
            double[,] upperpayoffs = new double[paths,steps];
            double[,] lowerpayoffs = new double[paths,steps];
            for (int x = 0; x < paths; x++)
            {
                for (int y = 0; y<steps; y++)
                {
                    if (Call == true)
                    {
                        upperpayoffs[x,y] = Math.Max(upperstock[x,y] - Strike, 0)*Math.Exp(-(r+deltar)*dt*y);
                        lowerpayoffs[x,y] = Math.Max(lowerstock[x,y] - Strike,0)*Math.Exp(-(r-deltar)*dt*y);
                    }
                    else
                    {
                        upperpayoffs[x,y] = Math.Max(Strike - upperstock[x,y], 0)*Math.Exp(-(r+deltar)*dt*y);
                        lowerpayoffs[x,y] = Math.Max(Strike - lowerstock[x,y],0)*Math.Exp(-(r-deltar)*dt*y);
                    }

                }
            }
            
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 0; y <steps; y++)
                {
                    //upperpath = Math.Max(stockvals[x,y]*Math.Exp(( (r+deltar) - (Math.Pow(v, 2))/2.0)*dt + (v)*Math.Sqrt(dt)*randoms[x,y]) - Strike,0);
                    //lowerpath = Math.Max(stockvals[x,y]*Math.Exp(( (r-deltar) - (Math.Pow(v, 2))/2.0)*dt + (v)*Math.Sqrt(dt)*randoms[x,y]) - Strike,0);
                    //Console.WriteLine(upperpath);
                    //Console.WriteLine(lowerpath);
                    //double k = upperpath - lowerpath;
                    double k = upperpayoffs[x,y] - lowerpayoffs[x,y];
                    //Console.WriteLine(k);
                    rholist.Add((k)/(2.0*deltar));
                    //Console.WriteLine(rholist.Count());
                    //Console.WriteLine(vegalist[count]);
                    //count +=1;

                }
            }
            double averagerho = rholist.Average(); // annualized

            return averagerho;
        }

        public double StandardError(double[,] oppayoffs)
        {
            int steps = oppayoffs.GetLength(1);
            int paths = oppayoffs.GetLength(0);
            List<double> pathpayoff = new List<double>();
            for (int x=0; x<paths;x++)
            {
                pathpayoff.Add(oppayoffs[x,steps-1]);
            }
            double avgpayoff = pathpayoff.Average();
            //Console.WriteLine("Average Payoff is {0}", avgpayoff);
            double standardeviation = 0d;
            double dpaths = Convert.ToDouble(paths);
            for(int x=0;x<paths;x++)
            {
                //double outscome = (1/(paths-1))*Math.Pow(pathpayoff[x] - avgpayoff,2);
                //Console.WriteLine(outscome);
                standardeviation += Convert.ToDouble((1.0/(dpaths-1.0)))*Convert.ToDouble(Math.Pow(pathpayoff[x] - avgpayoff,2));
            }
            double SE = Math.Sqrt(standardeviation/dpaths);

            return SE;
        }
    }
}
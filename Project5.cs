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
                    Console.WriteLine("How many steps for your option?");
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
            double optionDelta = option.optionsDelta(step_payoffs, stockrandsims[0]);
            double optionGamma = option.optionsGamma(step_payoffs,stockrandsims[0]);
            double optionTheta = option.optionsTheta(step_payoffs);
            double optionVega = option.optionsVega(stockrandsims[1], stockrandsims[0], riskfree);
            double optionRho = option.optionsRho(stockrandsims[1], stockrandsims[0], riskfree);


            Console.WriteLine("The option price via Monte Carlo is {0}", op_price);
            Console.WriteLine();
            Console.WriteLine("The StdError of the run was {0}", StdErr);
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

            //for (int x = 0; x < paths; x++)
            //{
                //Console.WriteLine(step_payoffs[x, steps-1]);
            //}

        }



    }
    class Gaussian // need a gaussian, might as well make a class that can be modified later for other methods
    {
        public static double sumtwelve(Random x1) // gonna use the sumtwelve method
        {
            List<double> gaussians = new List<double>();

            double k = 0d; //make zero a double using d as Chris showed me. Initial is 0

            for (int ii = 0; ii < 12; ii++)
            {
                k += x1.NextDouble(); // generate 12 numbers on the unit interval between 0 and 1 and that them together
            }

            k -= 6d; // once again change 6 to a double so that it can work with the operands 
            return k;
        }
    }
    class Stock
    {
        public double Value {get;set;}
        public double Vol {get;set;} // every stock has a volatility. As such the volatility of the options will just be that of the stock
        
        public List<double[,]> StockSim(double r, double T, int steps, int paths) // dont want static since it uses the stock characteristics
        {
            double dt = T/Convert.ToDouble(steps);
            double[,] PricePaths = new double[paths,steps]; // each row should be a path and each column a step
            double[,] PathRandoms = new double[paths,steps];

            for (int x = 0; x < steps;x++)
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
                        Random rand = new Random();
                        double Gauss = Gaussian.sumtwelve(rand);
                        PathRandoms[y,x] = Gauss;

                        PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
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
                    payoffs.Add(Math.Max(simmedstock[x, steps-1] - Strike, 0d)); // last column is the last step. So we want to iterate through rows
                    //count++;
                    //Console.WriteLine(count);
                }
                else
                {
                    payoffs.Add(Math.Max(Strike - simmedstock[x,steps-1],0d));
                }
            }
            return payoffs.Average()*Math.Exp(-r * Tenor);
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
        // the following Greeks formulas are from Lecture 4 slides 45,46
        public double optionsDelta(double[,] oppayoffs, double[,] stocks)
        {
            int steps = oppayoffs.GetLength(1);
            int paths = oppayoffs.GetLength(0);
            //double[,] deltas = new double[paths, steps-1];
            double delta = 0;

            for (int x = 0; x<paths;x++)
            {
                for (int y = 0; y<steps-1;y++)
                {
                    delta += (oppayoffs[x,y+1] - oppayoffs[x,y])/(stocks[x,y+1]- stocks[x,y]);
                }
            }
            int elements = (steps-1)*paths;
            double averagedelta = delta/elements; //average delta

            return averagedelta;
        }

        public double optionsGamma(double[,] oppayoffs, double[,] stocks)
        {
            int steps = oppayoffs.GetLength(1);
            int paths = oppayoffs.GetLength(0);
            //double[,] gammas = new double[paths, steps-2];
            double gamma = 0;

            for (int x = 0; x < paths ; x++)
            {
                for (int y = 0; y < steps-2; y++)
                {
                    gamma += Convert.ToDouble(((oppayoffs[x, y+2] - oppayoffs[x, y+1])/(stocks[x,y+2] - stocks[x,y+1]) - ((oppayoffs[x, y+1] - oppayoffs[x, y])/(stocks[x,y+1] - stocks[x,y])))/((1.0/2.0)*(stocks[x,y+2] - stocks[x,y])));
                    //Console.WriteLine(gamma);
                }
            }
            double dpaths = Convert.ToDouble(paths);
            double dsteps = Convert.ToDouble(steps);

            double averagegamma = (gamma)/(dpaths*(dsteps-2.0));
            return averagegamma;
        }

        public double optionsTheta(double[,] oppayoffs)
        {
            int steps = oppayoffs.GetLength(1);
            int paths = oppayoffs.GetLength(0);
            //double[,] thetas = new double[paths, steps-2];
            double theta = 0;
            double dt = Tenor/steps;

            for (int x = 0; x<paths;x++)
            {
                for (int y = 0; y<steps-2; y++)
                {
                    theta += (oppayoffs[x,y+2] - oppayoffs[x,y])/(2*dt);
                }
            }
            double averagetheta = theta/((steps-2)*paths);
            return averagetheta;
        }

        public double optionsVega(double[,] randoms, double[,] stockvals, double r)
        {
            int steps = randoms.GetLength(1);
            int paths = randoms.GetLength(0);
            double dt = Tenor/steps;
            double upperpath;
            double lowerpath;
            List<double> vegalist = new List<double>();
            double v = Underlying.Vol;
            double deltavol = (1/10.0)*v;
            //int count = 0;
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 1; y <steps; y++)
                {
                    upperpath = Math.Max(stockvals[x,y-1]*Math.Exp((r - (Math.Pow(v + deltavol, 2))/2)*dt + (v + deltavol)*Math.Sqrt(dt)*randoms[x,y]) - Strike,0);
                    lowerpath = Math.Max(stockvals[x,y-1]*Math.Exp((r - (Math.Pow(v - deltavol, 2))/2)*dt + (v - deltavol)*Math.Sqrt(dt)*randoms[x,y]) - Strike,0);
                    //Console.WriteLine(upperpath);
                    //Console.WriteLine(lowerpath);
                    double k = upperpath*1.0 -lowerpath*1.0;
                    double j = 2.0*deltavol;
                    //Console.WriteLine(k);
                    vegalist.Add(k*(1.0/j));
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
            int steps = randoms.GetLength(1);
            int paths = randoms.GetLength(0);
            double dt = Tenor/steps;
            double upperpath;
            double lowerpath;
            List<double> rholist = new List<double>();
            double v = Underlying.Vol;
            double deltar = (1/10.0)*r;
            //int count = 0;
            for (int x = 0; x < paths; x++)
            {
                //PricePaths[y,x] = PricePaths[y,x-1]*Math.Exp((r - (Math.Pow(Vol, 2))/2)*dt + Vol*Math.Sqrt(dt)*PathRandoms[y,x]);
                for (int y = 1; y <steps; y++)
                {
                    upperpath = Math.Max(stockvals[x,y-1]*Math.Exp((r - (Math.Pow(v + deltar, 2))/2)*dt + (v + deltar)*Math.Sqrt(dt)*randoms[x,y]) - Strike,0);
                    lowerpath = Math.Max(stockvals[x,y-1]*Math.Exp((r - (Math.Pow(v - deltar, 2))/2)*dt + (v - deltar)*Math.Sqrt(dt)*randoms[x,y]) - Strike,0);
                    //Console.WriteLine(upperpath);
                    //Console.WriteLine(lowerpath);
                    double k = upperpath*1.0 - lowerpath*1.0;
                    //Console.WriteLine(k);
                    rholist.Add((k)/(2.0*deltar));
                    //Console.WriteLine(rholist.Count());
                    //Console.WriteLine(vegalist[count]);
                    //count +=1;

                }
            }
            double averagerho = rholist.Average();

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
                standardeviation += Convert.ToDouble((1/(dpaths-1)))*Convert.ToDouble(Math.Pow(pathpayoff[x] - avgpayoff,2));
            }
            double SE = Math.Sqrt(standardeviation/dpaths);

            return SE;
        }
    }
}
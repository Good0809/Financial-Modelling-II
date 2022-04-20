using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Project6
{

    class Program
    {
        static void Main(string[] args)
        {
            //Console.WriteLine("check");
            Boolean done = false;
            double ival = 100d;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the price of your stock?");
                    ival = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    ival = 100d;
                    Console.WriteLine("Invalid entry, input must be a number!");
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
            double vola = 0.2;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the vol of your stock?");
                    vola = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    vola = 0.2;
                    Console.WriteLine("Invalid entry, input must be a number!");
                }
            }

            done = false;
            double rfree = 0.02;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the risk-free rate of your stock?");
                    rfree = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    rfree = 0.02;
                    Console.WriteLine("Invalid entry, input must be a number!");
                }
            }

            done = false;
            double T = 1d;
            while (!done)
            {
                try
                {
                    Console.WriteLine("What is the Tenor your option (in years)?");
                    T = Convert.ToDouble(Console.ReadLine());
                    done = true;
                }
                catch
                {
                    T = 1d;
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
            Boolean iscall = true;
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
                            iscall = true;
                        }
                        else
                        {
                            iscall = false;
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
                    iscall = true;
                    Console.WriteLine("Invalid entry, input must be an integer number 0 or 1!");
                }
            }
            done = false;
            Boolean isanti = true;
            while (!done)
            {
                try
                {
                    Console.WriteLine("Is this anti-corr sim? (0 = false, 1 = true)");
                    int input = int.Parse(Console.ReadLine());
                    if (input == 1 || input == 0)
                    {
                        if (input == 1)
                        {
                            isanti = true;
                        }
                        else
                        {
                            isanti = false;
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
                    isanti = true;
                    Console.WriteLine("Invalid entry, input must be an integer number 0 or 1!");
                }
            }
            Stock s = new Stock {Value = ival, Vol = vola};
            EuroOption e = new EuroOption {Underlying = s, Tenor = T, Strike = strike, Call = iscall};
            double[,] rands = Gaussian.GenerateNorms(paths,steps);
            List<double[,]> sa = s.simstock(paths,steps,rands, rfree,e.Tenor,isanti);

            double eval = e.PriceViaMonte(sa, rfree, isanti);
            double ds = s.Value/1000.0;
            double edelta = e.OpDelta(ds,rands,rfree,T,isanti);
            double egamma = e.OpGamma(ds,rands,rfree,T,isanti);
            double dvol = s.Vol/1000.0;
            double evega = e.OpVega(dvol,rands,rfree,T,isanti);
            double etheta = e.OpTheta(rands,rfree,T,isanti);
            double drho = rfree/1000.0;
            double erho = e.OpRho(drho, rands,rfree,T,isanti);
            string opname = "dumbdumb";
            //double stderror = e.StandardError(sa, rfree, isanti);
            if(iscall)
            {
                opname = "European Call";
            }
            else
            {
                opname = "European Put";
            }

            Console.WriteLine("Option price of {0} with initial val {1}, strike {2}, risk free {3}, and implied vol {4} with tenor {5}year(s) is {6} \n",opname, ival, strike, rfree, vola,T, eval);
            Console.WriteLine("It's Delta is approximately {0} \n",edelta);
            Console.WriteLine("It's Gamma is approximately {0} \n",egamma);
            Console.WriteLine("It's Vega is approximately {0} \n",evega);
            Console.WriteLine("It's Theta is approximately {0} \n", etheta);
            Console.WriteLine("It's Rho is approximately {0} \n", erho);
            //Console.WriteLine("The Standard Error of the run was {0}", stderror);


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

        public static double[,] GenerateNorms(int rows, int columns) //make the randoms outside of the functions to try and solve the greeks/stocks w/out issues of last project
        {
            double[,] rands = new double[rows,columns+1];
            for (int x = 0; x<rows; x++)
            {
                for (int y=0;y<columns+1;y++)
                {
                    Random x1, x2;
                    x1 = new Random();
                    x2 = new Random();
                    rands[x,y] = BoxMuller(x1,x2)[0];
                }
            }

            return rands;
        }
    }

    //the equations for stocks and greeks are pulled from lectures 10 and 12
    class Stock
    {
        public double Value {get;set;}
        public double Vol {get;set;}
        public List<double[,]> simstock(int rows, int columns, double[,] rands, double r, double Tenor, Boolean isanti, double ds = 0.0, double dv = 0.0, double dr = 0.0)
        {
            List<double[,]> endvals = new List<double[,]>();
            double dt = Tenor/Convert.ToDouble(columns);
            double[,] stockvals = new double[rows,columns+1];

            if (isanti)
            {
                double[,] antivals = new double[rows,columns+1];
                for(int x=0;x<columns+1;x++)
                {
                    for(int y=0;y<rows;y++)
                    {
                        if (x==0)
                        {
                            stockvals[y,x] = Value+ds;
                            antivals[y,x] = Value+ds;
                        }
                        else
                        {
                            stockvals[y,x] = stockvals[y,x-1]*Math.Exp((r+dr - (Math.Pow((Vol+dv), 2))/2.0)*(dt) + (Vol+dv)*Math.Sqrt(dt)*rands[y,x-1]);
                            antivals[y,x] = antivals[y,x-1]*Math.Exp((r+dr - (Math.Pow((Vol+dv), 2))/2.0)*(dt) - (Vol+dv)*Math.Sqrt(dt)*rands[y,x-1]);
                        }
                    }

                }
                endvals.Add(stockvals);
                endvals.Add(antivals);
                return endvals;
            }
            else
            {
                for(int x=0;x<columns+1;x++)
                {
                    for(int y=0;y<rows;y++)
                    {
                        if (x==0)
                        {
                            stockvals[y,x] = Value+ds;
                        }
                        else
                        {
                            //Console.WriteLine(y);
                            stockvals[y,x] = stockvals[y,x-1]*Math.Exp((r+dr - (Math.Pow(Vol+dv, 2))/2.0)*(dt) + (Vol+dv)*Math.Sqrt(dt)*rands[y,x-1]);
                        }
                    }

                }
                endvals.Add(stockvals);
                return endvals;
            }
        }
    }

    class EuroOption
    {
        public Stock Underlying {get;set;} //inherit the underlying
        public double Tenor {get;set;} // each european option has a tenor, strike, and type (call/put)
        public double Strike {get;set;}
        public Boolean Call {get;set;}
        // most of these are pulled from my last project and modified to fit anti-thetical runs. Modified as well to take randoms and use priveviamonte function w/efficency
        public double PriceViaMonte(List<double[,]> simmedstocks,double r, Boolean isanti, double dt = 0) // simmed stocks will either have just the regular stocks sims or that + antis
        {
            double price;
            int rows = simmedstocks[0].GetLength(0);
            int columns = simmedstocks[0].GetLength(1);
            //double dt = Tenor/columns;

            if (isanti)
            {
                double avgprice = 0;
                for(int x=0;x<rows;x++)
                {
                    if(Call)
                    {
                        avgprice += Math.Max(simmedstocks[0][x, columns-1] - Strike,0)*Math.Exp(-r*(Tenor+dt));
                        avgprice += Math.Max(simmedstocks[1][x, columns-1] - Strike,0)*Math.Exp(-r*(Tenor+dt));
                    }
                    else
                    {
                        avgprice += Math.Max(Strike - simmedstocks[0][x, columns-1],0)*Math.Exp(-r*(Tenor+dt));
                        avgprice += Math.Max(Strike - simmedstocks[1][x, columns-1],0)*Math.Exp(-r*(Tenor+dt));
                    }
                }

                price = avgprice/(2.0*rows);
                return price;
            }
            else
            {
                double avgprice = 0;
                for(int x=0;x<rows;x++)
                {
                    if(Call)
                    {
                        avgprice += Math.Max(simmedstocks[0][x, columns-1] - Strike,0)*Math.Exp(-r*(Tenor+dt));
                    }
                    else
                    {
                        avgprice += Math.Max(Strike - simmedstocks[0][x, columns-1],0)*Math.Exp(-r*(Tenor+dt));
                    }
                }

                price = avgprice/(rows);
                return price;
            } 

        }
        public double OpDelta(double dval, double[,] rands, double r, double Tenor, Boolean isanti)
        {
            int rows = rands.GetLength(0);
            int columns = rands.GetLength(1);

            List<double[,]> uppervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti, dval);
            List<double[,]> lowervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti, -dval);

            int listlen = uppervals.Count;

            double upper = PriceViaMonte(uppervals, r, isanti);
            double lower = PriceViaMonte(lowervals, r, isanti);

            double delta = (upper - lower)/(2.0*dval);
            return delta;
        }
        public double OpGamma(double dval, double[,] rands, double r, double Tenor, Boolean isanti)
        {
            int rows = rands.GetLength(0);
            int columns = rands.GetLength(1);

            List<double[,]> uppervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti, dval);
            List<double[,]> regvals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti, 0);
            List<double[,]> lowervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti, -dval);

            double upper = PriceViaMonte(uppervals, r, isanti);
            double reg = PriceViaMonte(regvals, r, isanti);
            double lower = PriceViaMonte(lowervals, r, isanti);

            double gamma = (upper - 2.0*reg + lower)/(Math.Pow(dval, 2));
            return gamma;
        }
        public double OpVega(double dval, double[,] rands, double r, double Tenor, Boolean isanti)
        {
            int rows = rands.GetLength(0);
            int columns = rands.GetLength(1);

            List<double[,]> uppervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti,0.0,dval,0.0);
            List<double[,]> lowervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti,0.0,-dval,0.0);

            double upper = PriceViaMonte(uppervals, r, isanti);
            double lower = PriceViaMonte(lowervals, r, isanti);

            double vega = (upper - lower)/(2.0*dval);
            return vega;
        }
        public double OpTheta(double[,] rands, double r, double Tenor, Boolean isanti)
        {
            int rows = rands.GetLength(0);
            int columns = rands.GetLength(1);
            double dt = Tenor/(columns);
            List<double> theta = new List<double>();
            int count = 0;
            for (int x = 0;x<columns;x++)
            {
                count +=1;
                List<double[,]> uppervals = Underlying.simstock(rows, x+1, rands, r, Tenor, isanti);
                List<double[,]> lowervals = Underlying.simstock(rows, x, rands, r, Tenor, isanti);
                
                double upper = PriceViaMonte(uppervals, r, isanti);
                double lower = PriceViaMonte(lowervals, r, isanti);


                theta.Add((upper - lower)/dt);
            }
            double theta_avg = (theta.Average())*-1.0;
            return theta_avg;
        }
        public double OpRho(double dval, double[,] rands, double r, double Tenor, Boolean isanti)
        {
            int rows = rands.GetLength(0);
            int columns = rands.GetLength(1);

            List<double[,]> uppervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti,0.0,0.0,dval);
            List<double[,]> lowervals = Underlying.simstock(rows, columns, rands, r, Tenor, isanti,0.0,0.0,-dval);

            double upper = PriceViaMonte(uppervals, r+dval, isanti);
            double lower = PriceViaMonte(lowervals, r-dval, isanti);

            double rho = (upper - lower)/(2.0*dval);
            return rho;
        }

        public double StandardError(List<double[,]> simmedstocks, double r, Boolean isanti)
        {
            int rows = simmedstocks[0].GetLength(0);
            int columns = simmedstocks[0].GetLength(1);
            double varsum = 0.0;
            if(isanti)
            {
                List<double[,]> antivals = new List<double[,]>();
                antivals.Add(simmedstocks[1]);
                for(int y=0;y<columns;y++)
                {
                    double column_mean_anti = PriceViaMonte(antivals, r, false);
                    double column_mean = PriceViaMonte(simmedstocks, r, false);
                    double temp1 = 0;
                    double temp2 = 0;
                    for(int x = 0; x < rows; x++)
                    {
                        temp1 += Math.Pow(simmedstocks[0][x,y] - column_mean,2);
                        temp2 += Math.Pow(simmedstocks[1][x,y] - column_mean_anti,2);
                        if(x == (rows-1))
                        {
                            varsum += (temp1/columns)/2.0 + (temp2/columns)/2.0;
                        }
                    }

                }
                double std = Math.Sqrt((1.0/2.0) * varsum * Math.Exp(-r*Tenor));
                double SE = Math.Sqrt(std/rows);
                return SE;
            }
            else
            {
                for(int y=0;y<columns;y++)
                {
                    double column_mean = PriceViaMonte(simmedstocks, r, false);
                    double temp1 = 0;
                    for(int x = 0; x < rows; x++)
                    {
                        temp1 += Math.Pow(simmedstocks[0][x,y] - column_mean,2);
                        if(x == (rows-1))
                        {
                            varsum += temp1;
                        }
                    }
                }
                double std = Math.Sqrt(varsum/columns);
                double SE = Math.Sqrt(std/rows);
                return SE;
            }

        }
    }



}

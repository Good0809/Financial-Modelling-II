using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace myapp
{
    class Program
    {
        static void Main(string [] args)
        {
            Random x1 = new Random(); //these will be our initial randoms
            Random x2 = new Random();
            Random x3 = new Random();

            Console.WriteLine("What would you like your correlation to be? (Must be between -1 and 1)");
            double p = -1.5; // need an initial value thats caught for the while

            while (p>1 || p<-1)
            {
                try
                {
                    p = Convert.ToDouble(Console.ReadLine()); // correlation must be between -1 and 1, as such we let the user know and have them renter
                    // otherwise, it doesnt enter any of the loops and goes through fine!
                    if (p>1 || p<-1)
                    {
                        Console.WriteLine("Invalid selection, please enter a number between -1 and 1.");
                    }
                }
                catch
                {
                    p = -1.5; // make sure it stays in the loop and lets the user know it must be a number
                    Console.WriteLine("Input must be a number. Please enter a number between -1 and 1.");
                }
            }

            double sumtwelvegauss = sumtwelve(x1); // this is our sum twelve gaussian outcome

            List<double> BoxMullerGauss, PolarRejectionGauss, CorrelatedGauss; // our lists for the other gaussian methods

            BoxMullerGauss = BoxMuller(x1,x2); // define the lists better
            PolarRejectionGauss = PolarRejection(x1,x2);
            CorrelatedGauss = Correlated(x1,x2,x3,p);

            //outputs 4 dayz

            Console.WriteLine("Gaussian from Sum Twelve is " + sumtwelvegauss);
            Console.WriteLine("Gaussians from Box Muller are z1 = " + BoxMullerGauss[0] +" , z2 = " + BoxMullerGauss[1]);
            Console.WriteLine("Gaussians from Polar Rejection are z1 = " + PolarRejectionGauss[0] +" , z2 = " + PolarRejectionGauss[1]);
            Console.WriteLine("Gaussians from correlated are z1 = " + CorrelatedGauss[0] +" , and z2 = " +CorrelatedGauss[1]);
        }

        static double sumtwelve(Random x1)
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

        static List<double> BoxMuller(Random x1, Random x2)
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

        static List<double> PolarRejection(Random x1, Random x2)
        {
            double y1,y2;
            y1 = x1.NextDouble(); // once again get our uniforms
            y2 = x2.NextDouble();

            double w = Math.Pow(y1,2) + Math.Pow(y2,2); // if initial is less the 1 it wont pass through the while
            while (w > 1) // if they're greater than one, it goes through the while loop until it isnt.
            {
                y1 = x1.NextDouble();
                y2 = x2.NextDouble();

                w = Math.Pow(y1,2) + Math.Pow(y2,2); 
            }

            //Console.WriteLine(w);

            double c = Math.Sqrt(-2*Math.Log(w)/w); // define our c and do the rest using formulas from page 53

            double z1,z2;
            z1 = c*y1;
            z2 = c*y2;

            List<double> gaussians = new List<double>();

            gaussians.Add(z1);
            gaussians.Add(z2);

            return gaussians; //return our list

        }

        static List<double> Correlated(Random x1, Random x2, Random x3, double p)
        {
            double y1,y2,y3;
            y1 = x1.NextDouble();
            y2 = x2.NextDouble();
            y3 = x3.NextDouble(); // now for correlation we need 3
            
            double z1,z2,z3; // z1 and z2 will double as our e1 and e2

            z1 = y1; // z1 and z2 are our e's as indicated on slide 55
            z2 = y2;

            z3 = p*z1 + Math.Sqrt(1-Math.Pow(p,2))*z2; // here we get our z3 using the correlated z1 and z2

            List<double> gaussians = new List<double>();

            gaussians.Add(z1); // add them to our list
            gaussians.Add(z3);

            return gaussians; // return the list

        }

    }

}

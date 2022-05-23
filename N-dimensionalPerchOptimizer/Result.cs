using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalPerchOptimizer
{
    public class Result
    {
        private Result(){ }

        public List<double> X;
        public List<double> X2;
        public List<double> X3;
        public List<double> U;
        public double fitness;

        private static Result Instance;
        public static Result GetInstance()
        {
            if (Instance == null)
            {
                Instance = new Result();
            }
            return Instance;
        }
    }
}

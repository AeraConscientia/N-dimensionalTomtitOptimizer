using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalPerchOptimizer
{
    public class Perch
    {
        public Vector coords;
        public double fitness = 0;

        public Perch(int dim)
        {
            coords = new Vector(dim);

            for (int i = 0; i < dim; i++)
                coords[i] = 0;
            fitness = 0;
        }
    }
}

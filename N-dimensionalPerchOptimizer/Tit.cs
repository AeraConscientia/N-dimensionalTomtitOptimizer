using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalPerchOptimizer
{
    public class Tit
    {
        public Vector coords;
        public double fitness = 0;
        public Tit best;
        public Tit local_best;
        public int dim;

        public Tit(int dim)
        {
            best = this;
            local_best = this;
            coords = new Vector(dim);

            for (int i = 0; i < dim; i++)
                coords[i] = 0;
            fitness = 0;
            this.dim = dim;
        }

        public Tit(Tit previousTit)
        {
            coords = new Vector(previousTit.coords);
            fitness = previousTit.fitness;
            best = previousTit.best;
            local_best = previousTit.local_best;
            dim = previousTit.dim;
        }
    }
}

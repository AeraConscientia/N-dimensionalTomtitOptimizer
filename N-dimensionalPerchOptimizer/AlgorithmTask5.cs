using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalTomtitOptimizer
{
    public class AlgorithmTask5 : AlgorithmTit
    {
        public AlgorithmTask5(double U1, double U2, double x0, double x1)
        {
            this.x0 = new List<double>(2) { x0, x1};   //Список одного элемента
            U = new List<Tuple<double, double>>();
            U.Add(new Tuple<double, double>(U1, U2));
        }

        public override void I(Tit tit, bool flag = false)
        {
            List<double> x1 = new List<double>();
            List<double> x2 = new List<double>();
            x1.Add(x0[0]); x2.Add(x0[1]);

            for (int i = 0; i < N_dim / 2; i++) // тут решилась проблема с количеством x
            {
                x1.Add(x1[i] + tit.coords[i]);
                x2.Add(2 * x1[i] + x2[i]);
            }

            double res1 = 0;
            double res2 = 0;
            double res3 = 0;
            for (int t = 0; t < N_dim / 2; t++) // N_dim или N_dim-1 ?
            {
                res1 += x1[t] * x1[t];
                res2 += x2[t] * x2[t];
                res3 += tit.coords[t] * tit.coords[t];
            }


            res1 += res2 + res3;
            tit.fitness = res1;
            if (flag == true)
            {
                Result result = Result.GetInstance();
                result.X = x1;
                result.X2 = x2;
                result.fitness = tit.fitness;
                result.U = new List<double>(N_dim / 2);

                for (int i = 0; i < N_dim / 2; i++)
                {
                    result.U.Add(tit.coords[i]);
                }
            }
        }

        public override Vector Levy()
        {
            Vector koefLevy = new Vector(N_dim);
            for (int i = 0; i < N_dim; i++)
            {
                if (i % 2 == 0)
                {
                    R1 = random.Next(Convert.ToInt32(0), Convert.ToInt32((U[0].Item2 - U[0].Item1) * 100)) / 100f; // (0, b1-a1)
                    thetta1 = R1 * 2 * Math.PI;
                    L1 = Math.Pow(R1 + 0.0001f, -1 / lambda);
                    koefLevy[i] = L1 * Math.Sin(thetta1);
                }
                else
                {
                    R2 = random.Next(Convert.ToInt32(0), Convert.ToInt32((U[0].Item2 - U[0].Item1) * 100)) / 100f; // (0, b2-a2)
                    thetta2 = R2 * 2 * Math.PI;
                    L2 = Math.Pow(R2 + 0.0001f, -1 / lambda);
                    koefLevy[i] = L2 * Math.Cos(thetta2);
                }
            }
            return koefLevy;
        }

        public override void GenerationAroundPool()
        {
            for (int i = 1; i < NP; i++)
            {
                for (int j = 0; j < N_dim; j++)
                {
                    double val = individuals[0].coords[j] + r * (-0.5 * (U[0].Item2 + U[0].Item1) + ((U[0].Item2 - U[0].Item1) * random.NextDouble() + U[0].Item1));

                    if (val < U[0].Item1)
                        val = (individuals[0].coords[j] - U[0].Item1) * random.NextDouble() + U[0].Item1;

                    if (val > U[0].Item2)
                        val = (U[0].Item2 - individuals[0].coords[j]) * random.NextDouble() + individuals[0].coords[j];

                    individuals[i].coords[j] = val;
                }

                I(individuals[i]);
            }
        }

        public override void GenerationAroundBest()
        {
            for (int i = 1; i < NP; i++)
            {
                for (int j = 0; j < N_dim; j++)
                {
                    double val = best.coords[j] + r * (-0.5 * (U[0].Item2 + U[0].Item1) + ((U[0].Item2 - U[0].Item1) * random.NextDouble() + U[0].Item1));

                    if (val < U[0].Item1)
                        val = (best.coords[j] - U[0].Item1) * random.NextDouble() + U[0].Item1;

                    if (val > U[0].Item2)
                        val = (U[0].Item2 - best.coords[j]) * random.NextDouble() + best.coords[j];

                    individuals[i].coords[j] = val;
                }
                I(individuals[i]);
            }
        }

        public override void BestCoordCorrect(Vector bestCoords, int k)
        {
            for (int j = 0; j < N_dim; j++)
            {
                if ((best.coords[j] < U[0].Item1) || (best.coords[j] > U[0].Item2))
                {
                    int NumTries = 0;
                    Vector tmp;
                    do
                    {
                        tmp = new Vector(Levy());
                        tmp[j] *= (alpha / (k + 1));
                        tmp[j] += bestCoords[j];
                        ++NumTries;
                    } while (((tmp[j] < U[0].Item1) || (tmp[j] > U[0].Item2)) && (NumTries <= 10));

                    best.coords[j] = (NumTries > 10) ? ((U[0].Item2 - U[0].Item1) * random.NextDouble() + U[0].Item1) : tmp[j];
                }
            }
        }

        public override void StohasticJump(Tit new_tit, double beta, int j)
        {
            if (beta < mu * h)
            {
                Vector Thetta = new Vector(individuals[j].coords.dim);
                for (int i = 0; i < individuals[j].coords.dim; i++)
                {
                    double delta_i = Math.Min(U[0].Item2 - new_tit.coords[i], new_tit.coords[i] - U[0].Item1);
                    Thetta[i] = random.NextDouble() * 2 * delta_i - delta_i;
                }

                new_tit.coords += Thetta;
                I(new_tit);
            }
        }

        public override void CheckBorders(Tit new_tit)
        {
            for (int i = 0; i < N_dim; i++)
            {
                if (new_tit.coords[i] < U[0].Item1)
                    new_tit.coords[i] = U[0].Item1;

                if (new_tit.coords[i] > U[0].Item2)
                    new_tit.coords[i] = U[0].Item2;
            }
        }

        public override void InitalPopulationGeneration()
        {
            for (int i = 0; i < NP; i++)
            {
                Tit tit = new Tit(N_dim);
                for (int j = 0; j < N_dim; j++)
                {
                    tit.coords[j] = U[0].Item1 + random.NextDouble() * (U[0].Item2 - U[0].Item1);
                }

                I(tit);
                individuals.Add(tit);
            }
        }
    }
}

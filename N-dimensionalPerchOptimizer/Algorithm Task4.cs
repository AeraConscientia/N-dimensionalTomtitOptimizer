using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalPerchOptimizer
{
    public class AlgorithmTask4 : AlgorithmPerch
    {
        public AlgorithmTask4(double U1, double U2, double x0)
        {
            this.x0 = new List<double>(1) { x0 };   //Список одного элемента
            U = new List<Tuple<double, double>>();
            U.Add(new Tuple<double, double>(U1, U2));
        }
        public override void FormingPopulation() // +
        {
            for (int i = 0; i < population; i++)
            {
                Perch perch = new Perch(N_dim);
                for (int j = 0; j < N_dim; j++)
                    perch.coords[j] = U[0].Item1 + rand.NextDouble() * (U[0].Item2 - U[0].Item1);

                I(perch);
                individuals.Add(perch);
            }
        }

        public override List<double> Levy() // *
        {
            List<double> koefLevy = new List<double>(N_dim);
            for (int i = 0; i < N_dim; i++)
            {
                if (i % 2 == 0)
                {
                    R1 = rand.Next(Convert.ToInt32(0), Convert.ToInt32((U[0].Item2 - U[0].Item1) * 100)) / 100f; // (0, b1-a1)
                    thetta1 = R1 * 2 * Math.PI;
                    L1 = Math.Pow(R1 + 0.0001f, -1 / lambda);
                    koefLevy.Add(L1 * Math.Sin(thetta1));
                }
                else
                {
                    R2 = rand.Next(Convert.ToInt32(0), Convert.ToInt32((U[0].Item2 - U[0].Item1) * 100)) / 100f; // (0, b2-a2)
                    thetta2 = R2 * 2 * Math.PI;
                    L2 = Math.Pow(R2 + 0.0001f, -1 / lambda);
                    koefLevy.Add(L2 * Math.Cos(thetta2));
                }
            }
            return koefLevy;

        }

        public override void I(Perch perch, bool flag = false)
        {
            List<double> x = new List<double>();
            x.Add(x0[0]);

            for (int i = 1; i < N_dim + 1; i++)
                x.Add(x[x.Count - 1] + perch.coords[i - 1]);

            double res = 0;
            for (int t = 0; t < N_dim; t++)
                res += ((x[t] * x[t]) + (perch.coords[t] * perch.coords[t]));

            perch.fitness = res;

            if (flag == true)
            {
                Result result = Result.GetInstance();
                result.X = x;
                result.fitness = perch.fitness;
                result.U = new List<double>(N_dim);
                for (int i = 0; i < N_dim; i++)
                {
                    result.U.Add(perch.coords[i]);
                }
            }
        }

        public override void MoveEPerchEFlock() // +
        {
            sigma = rand.NextDouble() * 0.4 + 0.1; // sigma [0.1,  0.5]

            for (int i = 0; i < NumFlocks; i++)
            {
                for (int j = 0; j < NumPerchInFlock; j++)
                {
                    int moveCount = (int)Math.Floor(sigma * NStep);

                    List<Perch> move = new List<Perch>();
                    for (int k = 0; k < moveCount; ++k)
                    {
                        Perch perch = new Perch(N_dim);
                        for (int l = 0; l < N_dim; l++)
                        {
                            double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                            if (tmp < U[0].Item1 || tmp > U[0].Item2)
                                tmp = flock[i, j].coords[l];
                            perch.coords[l] = tmp;
                        }
                        I(perch);
                        move.Add(perch);
                    }
                    move = move.OrderBy(s => s.fitness).ToList();
                    flock[i, j] = move[0];
                }
                Sort(flock, i);
            }
            SortFlocks();
        }

        protected override void BestFlockSwim() // +
        {
            sigma = rand.NextDouble() * 0.5 + 1; // sigma [1,  1.5]
            int i = 0;
            for (int j = 0; j < NumPerchInFlock; j++)
            {
                int moveCount = (int)Math.Floor(sigma * NStep);

                List<Perch> move = new List<Perch>();
                for (int k = 0; k < moveCount; ++k)
                {
                    Perch perch = new Perch(N_dim);
                    for (int l = 0; l < N_dim; l++) // если координаты не выйдут за область определения - меняем. Иначе - оставляем.
                    {
                        double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                        if (tmp < U[0].Item1 || tmp > U[0].Item2)
                            tmp = flock[i, j].coords[l];
                        perch.coords[l] = tmp;
                    }
                    I(perch);
                    move.Add(perch);
                }
                //Sort(move);
                move = move.OrderBy(s => s.fitness).ToList();
                flock[i, j] = move[0];
            }
            Sort(flock, i);
        }

        protected override void AverFlockSwim() // +
        {
            sigma = rand.NextDouble() / 20 + 0.6; // sigma [0.6,  0.8]

            for (int l = 1; l < NumFlocks - 1; l++) // если не изменяет память, передвижение лидеров средних стай
            {
                int moveCount = (int)Math.Floor(sigma * NStep);

                List<Perch> move = new List<Perch>();
                for (int k = 0; k < moveCount; ++k)
                {
                    Perch perch = new Perch(N_dim);
                    for (int m = 0; m < N_dim; m++)
                    {
                        double tmp = flock[l, 0].coords[m] + k * ((flock[0, 0].coords[m] - flock[l, 0].coords[m]) / (NStep));
                        if (tmp < U[0].Item1 || tmp > U[0].Item2)
                            tmp = flock[l, 0].coords[m];
                        perch.coords[m] = tmp; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    }
                    I(perch);
                    move.Add(perch);
                }

                //Sort(move);
                move = move.OrderBy(s => s.fitness).ToList();
                flock[l, 0] = move[0];

                for (int j = 0; j < NumPerchInFlock; j++)
                {
                    int moveCount1 = (int)Math.Floor(sigma * NStep);

                    List<Perch> move1 = new List<Perch>();
                    for (int k = 0; k < moveCount1; ++k)
                    {
                        Perch perch = new Perch(N_dim);
                        for (int p = 0; p < N_dim; p++)
                        {
                            double tmp = flock[l, j].coords[p] + k * ((flock[l, 0].coords[p] - flock[l, j].coords[p]) / (NStep));
                            if (tmp < U[0].Item1 || tmp > U[0].Item2)
                                tmp = flock[l, j].coords[p];
                            perch.coords[p] = tmp;
                        }
                        I(perch);
                        move1.Add(perch);
                    }
                    move1 = move1.OrderBy(s => s.fitness).ToList();
                    flock[l, j] = move1[0];
                }
                Sort(flock, l);
            }
            SortFlocks();
        }

        protected override void PoorFlockSwim() // +
        {
            List<double> PoorLeaderCoord = new List<double>();
            for (int p = 0; p < N_dim; p++)
            {
                PoorLeaderCoord.Add(flock[NumFlocks - 1, 0].coords[p]);
            }
            PoorLeaderSwim();
            //WrongCoord();
            for (int p = 0; p < N_dim / 2; p++)
            {
                int NumTries = 0;
                while (((flock[NumFlocks - 1, 0].coords[p] < U[0].Item1) || (flock[NumFlocks - 1, 0].coords[p] > U[0].Item2)))
                {
                    for (int w = 0; w < N_dim; w++)
                    {
                        flock[NumFlocks - 1, 0].coords[w] = PoorLeaderCoord[w];
                    }
                    PoorLeaderSwim();       //!
                    NumTries += 1;
                    if (NumTries == 10)
                        break;
                }

                bool ok = false;
                for (int d = 0; d < N_dim; ++d)
                {
                    if (flock[NumFlocks - 1, 0].coords[d] < U[0].Item1 || flock[NumFlocks - 1, 0].coords[d] > U[0].Item2)
                    {
                        ok = true;
                        break;
                    }
                }
                if (ok)
                    for (int pL = 0; pL < N_dim; pL++)
                    {
                        flock[NumFlocks - 1, 0].coords[pL] = U[0].Item1 + rand.NextDouble() * (U[0].Item2 - U[0].Item1);
                    }
            }

            sigma = rand.NextDouble() * 0.4 + 0.1; // sigma [0.1,  0.5]

            List<double> Min = new List<double>(N_dim);

            for (int j = 0; j < N_dim; j++)
                Min.Add(Math.Min((flock[NumFlocks - 1, 0].coords[j] - U[0].Item1), (U[0].Item2 - flock[NumFlocks - 1, 0].coords[j])));

            for (int j = 1; j < NumPerchInFlock; j++)
            {
                List<Tuple<double, double>> coords = new List<Tuple<double, double>>(N_dim);
                List<double> res = new List<double>(N_dim);

                for (int p = 0; p < N_dim; p++)
                {
                    double x1 = flock[NumFlocks - 1, 0].coords[p] - Min[p];
                    double x2 = flock[NumFlocks - 1, 0].coords[p] + Min[p];
                    coords.Add(new Tuple<double, double>(x1, x2));
                }

                Perch perch = new Perch(N_dim);
                for (int k = 0; k < N_dim; k++)
                {
                    res.Add(coords[k].Item1 + rand.NextDouble() * (coords[k].Item2 - coords[k].Item1));
                    //double tmp;
                    //do
                    //{
                    //    tmp = ((rand.NextDouble()) * 2 - 1) * (flock[NumFlocks - 1, 0].coords[k] - Min[k]);
                    //} while (tmp < U[0].Item1 || tmp > U[0].Item2);
                    perch.coords[k] = res[k]; //*
                }

                I(perch);
                flock[NumFlocks - 1, j] = perch; //*
            }

            int i = 1;

            for (int j = 0; j < NumPerchInFlock; j++) // всех окуней из худших двигаем к лидеру худшей стаи
            {
                int moveCount = (int)Math.Floor(sigma * NStep);

                List<Perch> move = new List<Perch>();
                for (int k = 0; k < moveCount; ++k)
                {
                    Perch perch = new Perch(N_dim);
                    for (int p = 0; p < N_dim; p++)
                    {
                        double tmp = flock[i, j].coords[p] + k * ((flock[i, 0].coords[p] - flock[i, j].coords[p]) / (NStep));
                        if (tmp < U[0].Item1 || tmp > U[0].Item2)
                            tmp = flock[i, j].coords[p];
                        perch.coords[p] = tmp;
                    }
                    I(perch);
                    move.Add(perch);
                }
                //Sort(move);
                move = move.OrderBy(s => s.fitness).ToList();
                flock[i, j] = move[0];
            }
            Sort(flock, i);
        }

    }
}

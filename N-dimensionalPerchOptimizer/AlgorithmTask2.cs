using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalPerchOptimizer
{
    public class AlgorithmTask2 : AlgorithmPerch
    {
        public AlgorithmTask2(double U11, double U12, double U21, double U22, double U31, double U32, double x0, double x1, double x2)
        {
            this.x0 = new List<double>(3) { x0, x1, x2 };   //Список одного элемента
            U = new List<Tuple<double, double>>();
            U.Add(new Tuple<double, double>(U11, U12));
            U.Add(new Tuple<double, double>(U21, U22));
            U.Add(new Tuple<double, double>(U31, U32));
        }
        public override void FormingPopulation() // +
        {
            for (int i = 0; i < population; i++)
            {
                Perch perch = new Perch(N_dim);
                for (int j = 0; j < N_dim/3; j++)
                    perch.coords[j] = U[0].Item1 + rand.NextDouble() * (U[0].Item2 - U[0].Item1);

                for (int j = N_dim/3; j < 2 * N_dim/3; j++)
                    perch.coords[j] = U[1].Item1 + rand.NextDouble() * (U[1].Item2 - U[1].Item1);

                for (int j = 2 * N_dim/3; j < N_dim; j++)
                    perch.coords[j] = U[2].Item1 + rand.NextDouble() * (U[2].Item2 - U[2].Item1);

                I(perch);
                individuals.Add(perch);
            }
        }

        public override List<double> Levy() // *
        {
            List<double> koefLevy = new List<double>(N_dim);
            for (int i = 0; i < N_dim/3; i++)
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

            for (int i = N_dim/3; i < 2*N_dim/3; i++)
            {
                if (i % 2 == 0)
                {
                    R1 = rand.Next(Convert.ToInt32(0), Convert.ToInt32((U[1].Item2 - U[1].Item1) * 100)) / 100f; // (0, b1-a1)
                    thetta1 = R1 * 2 * Math.PI;
                    L1 = Math.Pow(R1 + 0.0001f, -1 / lambda);
                    koefLevy.Add(L1 * Math.Sin(thetta1));
                }
                else
                {
                    R2 = rand.Next(Convert.ToInt32(0), Convert.ToInt32((U[1].Item2 - U[1].Item1) * 100)) / 100f; // (0, b2-a2)
                    thetta2 = R2 * 2 * Math.PI;
                    L2 = Math.Pow(R2 + 0.0001f, -1 / lambda);
                    koefLevy.Add(L2 * Math.Cos(thetta2));
                }
            }

            for (int i = 2*N_dim / 3; i < N_dim; i++)
            {
                if (i % 2 == 0)
                {
                    R1 = rand.Next(Convert.ToInt32(0), Convert.ToInt32((U[2].Item2 - U[2].Item1) * 100)) / 100f; // (0, b1-a1)
                    thetta1 = R1 * 2 * Math.PI;
                    L1 = Math.Pow(R1 + 0.0001f, -1 / lambda);
                    koefLevy.Add(L1 * Math.Sin(thetta1));
                }
                else
                {
                    R2 = rand.Next(Convert.ToInt32(0), Convert.ToInt32((U[2].Item2 - U[2].Item1) * 100)) / 100f; // (0, b2-a2)
                    thetta2 = R2 * 2 * Math.PI;
                    L2 = Math.Pow(R2 + 0.0001f, -1 / lambda);
                    koefLevy.Add(L2 * Math.Cos(thetta2));
                }
            }
            return koefLevy;
        }

        public override void I(Perch perch, bool flag = false)
        {
            //List<double> x = new List<double>(); // убрать х

            List<double> x1 = new List<double>();
            List<double> x2 = new List<double>();
            List<double> x3 = new List<double>();
            x1.Add(x0[0]); x2.Add(x0[1]); x3.Add(x0[2]);

            for (int i = 0; i < N_dim / 3; i++) // тут решилась проблема с количеством x
            {
                x1.Add((double)x1[i] / (1f + 0.01 * perch.coords[i] * (3 + perch.coords[N_dim / 3 + i])));
                x2.Add(((double)x2[i] + perch.coords[i] * x1[i+1]) / (1f + perch.coords[i] * (1f + perch.coords[N_dim / 3 + i])));
                x3.Add((double)x3[i] / (1f + 0.01 * perch.coords[N_dim / 3 + i] * (1 + perch.coords[2 * N_dim / 3 + i])));
            }

            double res1 = 0;
            double res2 = 0;
            for (int t = 0; t < N_dim/3; t++) // N_dim или N_dim-1 ?
            {
                res1 += x1[t] * x1[t] + x2[t] * x2[t] + 2 * perch.coords[2 * N_dim/3 + t] * perch.coords[2 * N_dim/3 + t];
                res2 += x3[t] * x3[t] + 2 * perch.coords[t] * perch.coords[t] + 2 * perch.coords[N_dim/3 + t] * perch.coords[N_dim/3 + t];
            }
                

            res1 *=  res2;
            res1 = Math.Sqrt(res1);
            res1 += x1[N_dim/3 - 1] * x1[N_dim/ 3 - 1] + x2[N_dim/ 3 - 1] * x2[N_dim/ 3 - 1] + x3[N_dim/ 3 - 1] * x3[N_dim/ 3 - 1];
            perch.fitness = res1;
            if (flag == true)
            {
                Result result = Result.GetInstance();
                result.X  = x1;
                result.X2 = x2;
                result.X3 = x3;
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
                        for (int l = 0; l < N_dim / 3; l++)
                        {
                            double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                            if (tmp < U[0].Item1 || tmp > U[0].Item2)
                                tmp = flock[i, j].coords[l];
                            perch.coords[l] = tmp;
                        }
                        for (int l = N_dim/3; l < 2*N_dim/3; l++)
                        {
                            double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                            if (tmp < U[1].Item1 || tmp > U[1].Item2)
                                tmp = flock[i, j].coords[l];
                            perch.coords[l] = tmp;
                        }
                        for (int l = 2*N_dim/3; l < N_dim; l++)
                        {
                            double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                            if (tmp < U[2].Item1 || tmp > U[2].Item2)
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
                    for (int l = 0; l < N_dim/3; l++) // если координаты не выйдут за область определения - меняем. Иначе - оставляем.
                    {
                        double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                        if (tmp < U[0].Item1 || tmp > U[0].Item2)
                            tmp = flock[i, j].coords[l];
                        perch.coords[l] = tmp;
                    }
                    for (int l = N_dim/3; l < 2*N_dim/3; l++) // если координаты не выйдут за область определения - меняем. Иначе - оставляем.
                    {
                        double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                        if (tmp < U[1].Item1 || tmp > U[1].Item2)
                            tmp = flock[i, j].coords[l];
                        perch.coords[l] = tmp;
                    }
                    for (int l = 2*N_dim/3; l < N_dim; l++) // если координаты не выйдут за область определения - меняем. Иначе - оставляем.
                    {
                        double tmp = flock[i, j].coords[l] + k * ((flock[i, 0].coords[l] - flock[i, j].coords[l]) / (NStep));
                        if (tmp < U[2].Item1 || tmp > U[2].Item2)
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

                    for (int m = 0; m < N_dim/3; m++)
                    {
                        double tmp = flock[l, 0].coords[m] + k * ((flock[0, 0].coords[m] - flock[l, 0].coords[m]) / (NStep));
                        if (tmp < U[0].Item1 || tmp > U[0].Item2)
                            tmp = flock[l, 0].coords[m];
                        perch.coords[m] = tmp;
                    }
                    for (int m = N_dim/3; m < 2*N_dim/3; m++)
                    {
                        double tmp = flock[l, 0].coords[m] + k * ((flock[0, 0].coords[m] - flock[l, 0].coords[m]) / (NStep));
                        if (tmp < U[1].Item1 || tmp > U[1].Item2)
                            tmp = flock[l, 0].coords[m];
                        perch.coords[m] = tmp;
                    }
                    for (int m = 2*N_dim/3; m < N_dim; m++)
                    {
                        double tmp = flock[l, 0].coords[m] + k * ((flock[0, 0].coords[m] - flock[l, 0].coords[m]) / (NStep));
                        if (tmp < U[2].Item1 || tmp > U[2].Item2)
                            tmp = flock[l, 0].coords[m];
                        perch.coords[m] = tmp;
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
                        for (int p = 0; p < N_dim/3; p++)
                        {
                            double tmp = flock[l, j].coords[p] + k * ((flock[l, 0].coords[p] - flock[l, j].coords[p]) / (NStep));
                            if (tmp < U[0].Item1 || tmp > U[0].Item2)
                                tmp = flock[l, j].coords[p];
                            perch.coords[p] = tmp;
                        }
                        for (int p = N_dim/3; p < 2*N_dim/3; p++)
                        {
                            double tmp = flock[l, j].coords[p] + k * ((flock[l, 0].coords[p] - flock[l, j].coords[p]) / (NStep));
                            if (tmp < U[1].Item1 || tmp > U[1].Item2)
                                tmp = flock[l, j].coords[p];
                            perch.coords[p] = tmp;
                        }
                        for (int p = 2*N_dim/3; p < N_dim; p++)
                        {
                            double tmp = flock[l, j].coords[p] + k * ((flock[l, 0].coords[p] - flock[l, j].coords[p]) / (NStep));
                            if (tmp < U[2].Item1 || tmp > U[2].Item2)
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

        /// <summary>Леви специально для задачи с 3 управлениями.</summary>
        private void PoorLeaderSwimTask2(int N_dimStart, int N_dimEnd)
        {
            List<double> koefLevy = Levy();
            for (int i = N_dimStart; i < N_dimEnd; i++)
                flock[NumFlocks - 1, 0].coords[i] = flock[NumFlocks - 1, 0].coords[i] + (alfa / currentIteration) * koefLevy[i];
        }

        private void CheckCoords(int i, int j, int N_dimStart, int N_dimStop, double ULeft, double URight, int moveCount)
        {
            for (int p = N_dimStart; p < N_dimStop; p++)
            {
                double tmp = flock[i, j].coords[p] + moveCount * ((flock[i, 0].coords[p] - flock[i, j].coords[p]) / (NStep));
                if (tmp < ULeft || tmp > URight)
                    tmp = flock[i, j].coords[p];
                perch.coords[p] = tmp;
            }
        }


        /// <summary>Проверка, есть ли выход за границы после Леви. Если есть - применяется равномерное распределение</summary>
        private void PoorFlockSwimCheckCoords(int N_dimStart, int N_dimStop, double ULeft, double URight)
        {
            bool ok = false;
            for (int d = N_dimStart; d < N_dimStop; ++d)
            {
                if (flock[NumFlocks - 1, 0].coords[d] < ULeft || flock[NumFlocks - 1, 0].coords[d] > URight)// ||flock[NumFlocks - 1, 0].coords[d]  U[0].Item1 || flock[NumFlocks - 1, 0].coords[d] > U[0].Item2)
                {
                    ok = true;
                    break;
                }
            }

            if (ok == true)
                for (int pL = N_dimStart; pL < N_dimStop; pL++)
                {
                    flock[NumFlocks - 1, 0].coords[pL] = ULeft + rand.NextDouble() * (URight - U[0].Item1);
                }
        }

        /// <summary>Попытка сделать через распределение Леви</summary>
        private void PoorFlockSwimLevy(int N_dimStart, int N_dimStop, double ULeft, double URight, List<double> PoorLeaderCoord)
        {
            for (int p = N_dimStart; p < N_dimStop; p++)
            {
                int NumTries = 0;
                while (((flock[NumFlocks - 1, 0].coords[p] < ULeft) || (flock[NumFlocks - 1, 0].coords[p] > URight)))// && (NumTries < 10))
                {
                    for (int w = N_dimStart; w < N_dimStop; w++)
                        flock[NumFlocks - 1, 0].coords[w] = PoorLeaderCoord[w];

                    PoorLeaderSwimTask2(N_dimStart, N_dimStop);

                    NumTries += 1;
                    if (NumTries == 10)
                        break;
                }
            }
        }

        protected override void PoorFlockSwim() // +
        {
            List<double> PoorLeaderCoord = new List<double>();
            for (int p = 0; p < N_dim; p++)
            {
                PoorLeaderCoord.Add(flock[NumFlocks - 1, 0].coords[p]);
            }

            PoorFlockSwimLevy(          0,          N_dim / 3,      U[0].Item1,     U[0].Item2, PoorLeaderCoord);
            PoorFlockSwimCheckCoords(   0,          N_dim / 3,      U[0].Item1,     U[0].Item2);
            //WrongLeader();
            //WrongCoord();
            PoorFlockSwimLevy(          N_dim / 3, 2 * N_dim / 3,   U[1].Item1,     U[1].Item2, PoorLeaderCoord);
            PoorFlockSwimCheckCoords(   N_dim / 3, 2 * N_dim / 3,   U[1].Item1,     U[1].Item2);

            //WrongLeader();
            //WrongCoord();
            PoorFlockSwimLevy(       2 * N_dim / 3,         N_dim,  U[2].Item1,     U[2].Item2, PoorLeaderCoord);
            PoorFlockSwimCheckCoords(2 * N_dim / 3,         N_dim,  U[2].Item1,     U[2].Item2);


            //WrongLeader();
            //WrongCoord();
            sigma = rand.NextDouble() * 0.4 + 0.1; // sigma [0.1,  0.5]

            List<double> Min = new List<double>(N_dim);

            for (int j = 0; j < N_dim/3; j++)
                Min.Add(Math.Min((flock[NumFlocks - 1, 0].coords[j] - U[0].Item1), (U[0].Item2 - flock[NumFlocks - 1, 0].coords[j])));
            for (int j = N_dim/3; j < 2*N_dim/3; j++)
                Min.Add(Math.Min((flock[NumFlocks - 1, 0].coords[j] - U[1].Item1), (U[1].Item2 - flock[NumFlocks - 1, 0].coords[j])));
            for (int j = 2*N_dim/3; j < N_dim; j++)
                Min.Add(Math.Min((flock[NumFlocks - 1, 0].coords[j] - U[2].Item1), (U[2].Item2 - flock[NumFlocks - 1, 0].coords[j])));

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

                    CheckCoords(i, j,               0,      N_dim / 3, U[0].Item1, U[0].Item2, k);
                    CheckCoords(i, j,       N_dim / 3,  2 * N_dim / 3, U[1].Item1, U[1].Item2, k);
                    CheckCoords(i, j,   2 * N_dim / 3,          N_dim, U[2].Item1, U[2].Item2, k);

                    I(perch);
                    move.Add(perch);
                }
                //Sort(move);
                move = move.OrderBy(s => s.fitness).ToList();
                flock[i, j] = move[0];
            }
            //WrongCoord();
            Sort(flock, i);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // часть для проверки результатов

        //private void WrongCoord()
        //{
        //    for (int i = 0; i < NumFlocks; i++)
        //    {
        //        for (int j = 0; j < NumPerchInFlock; j++)
        //        {
        //            for (int k = 0; k < N_dim; k++)
        //            {
        //                if (flock[i, j].coords[k] > 4)
        //                {
        //                    //int yyy = 0;
        //                    throw new Exception();
        //                }
        //            }
        //        }
        //    }
        //}

        //private void WrongLeader()
        //{
        //    for (int i = 0; i < N_dim; i++)
        //    {
        //        if (flock[NumFlocks - 1, 0].coords[i] > 4)
        //            throw new Exception();
        //    }
        //}
    }
    
}



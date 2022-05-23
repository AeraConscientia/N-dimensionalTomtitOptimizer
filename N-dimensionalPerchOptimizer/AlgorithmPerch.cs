using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization.Formatters;
using System.Text;

namespace N_dimensionalPerchOptimizer
{
    public abstract class AlgorithmPerch
    {
        public double R1, R2;
        public double thetta1, thetta2;
        public double L1, L2;

        protected Random rand = new Random();
        public Perch perch;
        public Perch best;
        public Perch result;

        public int N_dim = 0;
        public List<Tuple<double, double>> U;
        public List<double> x0;
        /// <summary>Дополнительные параметры системы</summary>
        public List<double> additionalParams;
        /// <summary>Номер выбранной задачи</summary>
        public int ExampleNum = 0;

        /// <summary>Начальные состояния системы</summary>
        public double Xn1, Xn2, Xn3;

        public List<double> Xn;
        public List<double> X;

        /// <summary>Размер популяции окуней </summary>
        public int population;

        /// <summary>Количество стай</summary>
        public int NumFlocks = 0;
        /// <summary>Количество окуней в стае</summary>
        public int NumPerchInFlock = 0;
        /// <summary>Количество шагов до окончания движения внутри стаи</summary>
        public int NStep = 0;
        /// <summary>Глубина продвижения внутри котла</summary>
        public double sigma = 0;

        /// <summary>Параметр распределения Леви</summary>
        public double lambda = 0;
        /// <summary>Величина шага</summary>
        public double alfa = 0;

        /// <summary>Число перекоммутаций</summary>
        public int PRmax = 0;
        /// <summary>Число шагов при перекоммутации</summary>
        public int deltapr = 0;

        /// <summary>Номер выбранной функции</summary>
        public int f;

        /// <summary>Максимальное число итераций</summary>
        public int MaxCount { get; set; }

        /// <summary>Текущая итерация </summary>
        public int currentIteration = 1;

        /// <summary>Массив средней приспособленности</summary>
        public List<double> averageFitness = new List<double>();

        /// <summary>Массив лучшей приспособленности</summary>
        public List<double> bestFitness = new List<double>();

        /// <summary>Популяция окуней</summary>
        public List<Perch> individuals = new List<Perch>();

        public List<Perch> Pool = new List<Perch>();

        public bool flagCreate = false;

        /// <summary>Все стаи с окунями</summary>
        public Perch[,] flock;

        public AlgorithmPerch() { }

        /// <summary>Сортировка ВСЕХ окуней</summary>
        protected void Sort(List<Perch> list) // не зашел
        {

            for (int i = 0; i < list.Count; i++)
                for (int j = 0; j < list.Count - i - 1; j++)
                    if (list[j].fitness > list[j + 1].fitness)
                    {
                        Perch tmp = list[j];
                        list[j] = list[j + 1];
                        list[j + 1] = tmp;
                    }
        }


        /// <summary>Сортировка всех окуней в своей стае</summary>
        protected void Sort(Perch[,] perches, int flockIndex) // *
        {
            for (int i = 0; i < NumPerchInFlock; i++)
                for (int j = 0; j < NumPerchInFlock - i - 1; j++)
                    if (perches[flockIndex, j].fitness > perches[flockIndex, j + 1].fitness)
                    {
                        Perch tmp = perches[flockIndex, j];
                        perches[flockIndex, j] = perches[flockIndex, j + 1];
                        perches[flockIndex, j + 1] = tmp;
                    }
        }

        protected void SortFlocks() // *
        {
            for (int i = 0; i < NumFlocks; i++)
            {
                for (int j = 0; j < NumFlocks - i - 1; j++)
                {
                    if (flock[j, 0].fitness > flock[j + 1, 0].fitness)
                    {
                        Perch[] tmp = new Perch[NumPerchInFlock];
                        for (int k = 0; k < NumPerchInFlock; k++)
                            tmp[k] = flock[j, k];
                        for (int k = 0; k < NumPerchInFlock; k++)
                            flock[j, k] = flock[j + 1, k];
                        for (int k = 0; k < NumPerchInFlock; k++)
                            flock[j + 1, k] = tmp[k];
                    }
                }
            }
        }

        /// <summary>Разбивка окуней на стаи</summary>
        public void MakeFlocks() // *
        {
            individuals = individuals.OrderBy(s => s.fitness).ToList();
            flock = new Perch[NumFlocks, NumPerchInFlock];

            for (int i = 0; i < individuals.Count; i++)
            {
                Math.DivRem(i, NumFlocks, out int tmp);
                flock[tmp, i / (NumFlocks)] = individuals[i];
            }
            //flock[0][i] - перебор окуней в первой стае
            //flock[-1][i] - перебор окуней в первой стае
        }

        /// <summary>Движения каждого окуня в каждой стае, создание котлов</summary>
        public abstract void MoveEPerchEFlock();

        protected abstract void BestFlockSwim();

        /// <summary>Движение средних окуней</summary>
        protected abstract void AverFlockSwim();

        protected abstract void PoorFlockSwim();

        public abstract void I(Perch perch, bool flag = false);
        /// <summary>Новые координаты лидера худшей стаи</summary>
        public void PoorLeaderSwim() // *
        {
            List<double> koefLevy = Levy();
            for (int i = 0; i < N_dim; i++)
                flock[NumFlocks - 1, 0].coords[i] = flock[NumFlocks - 1, 0].coords[i] + (alfa / currentIteration) * koefLevy[i];
        }

        /// <summary>Создание N-мерных коэф Леви</summary>
        /// <returns>Список коэффициентов для распределения Леви</returns>
        public abstract List<double> Levy();

        /// <summary>Начальное формирование популяции </summary>
        public abstract void FormingPopulation();

        public Perch StartAlg(int MaxCount,
            int NumFlocks, int NumPerchInFlock,
            int NStep,
            double lambda, double alfa,
            int PRmax, int deltapr, int N_dim)
        {

            this.MaxCount = MaxCount;

            this.NumFlocks = NumFlocks;
            this.NumPerchInFlock = NumPerchInFlock;
            population = NumFlocks * NumPerchInFlock;

            this.NStep   = NStep;
            this.lambda  = lambda;
            this.alfa    = alfa;
            this.PRmax   = PRmax;
            this.deltapr = deltapr;
            this.N_dim   = N_dim;

            perch   = new Perch(N_dim);
            best    = new Perch(N_dim);
            result  = new Perch(N_dim);

            FormingPopulation();

            Console.WriteLine("Выполняемые итерации:");
            for (int currentIteration = 1; currentIteration < MaxCount; currentIteration++)
            {
                Console.WriteLine(currentIteration);
                MakeFlocks();
                //WrongCoord();

                MoveEPerchEFlock();
                //WrongCoord();

                FlocksSwim();
                //WrongCoord();

                this.currentIteration++;
            }

            //WrongCoord();
            Recommutation(); // эта штука мне все портит.
            //WrongCoord();

            perch = Pool[0];
            I(perch, true);
            return perch;
        }

        public void FlocksSwim() // +
        {
            BestFlockSwim();
            //WrongCoord();

            PoorFlockSwim();
            //WrongCoord();

            AverFlockSwim();
            //WrongCoord();

            individuals = new List<Perch>();

            for (int i = 0; i < NumFlocks; i++)
                for (int j = 0; j < NumPerchInFlock; j++)
                    individuals.Add(flock[i, j]);
            Pool.Add(flock[0, 0]);
        }

        public void Recommutation() // *
        {
            int p, q, r;

            for (int pr = 0; pr < PRmax; pr++)
            {
                Perch perchResult = new Perch(N_dim);
                Perch Xp_pool = new Perch(N_dim);
                Perch Xq_pool = new Perch(N_dim);
                Perch Xr_pool = new Perch(N_dim);
                p = rand.Next(0, MaxCount - 1);
                q = rand.Next(0, MaxCount - 1);
                while (q == p)
                {
                    q = rand.Next(0, MaxCount - 1);
                }
                r = rand.Next(0, MaxCount - 1);
                while (r == p || r == q)
                {
                    r = rand.Next(0, MaxCount - 1);
                }

                Xp_pool = Pool[p];
                Xq_pool = Pool[q];
                Xr_pool = Pool[r];

                double min = Double.MaxValue; // попався
                for (int i = 0; i < deltapr - 1; i++)
                {
                    Perch perch = new Perch(N_dim);
                    for (int j = 0; j < N_dim; j++) // Очень логично, здесь начинаются проблемы с выходом за границы.
                                                    // идея правки... я пока не знаю.
                    {
                        double tmp = Xp_pool.coords[j] + i * (Xq_pool.coords[j] - Xp_pool.coords[j]) / deltapr; // TODO: а это что за бред?
                        perch.coords[j] = tmp;
                    }
                    I(perch);

                    if (perch.fitness < min)
                    {
                        perchResult = perch;
                        min = perch.fitness;
                    }
                }

                min = Double.MaxValue; // попався
                Perch inPool = new Perch(N_dim);
                for (int i = 0; i < deltapr - 1; i++)
                {
                    Perch perch = new Perch(N_dim);
                    for (int j = 0; j < N_dim; j++) // Очень логично, здесь начинаются проблемы с выходом за границы.
                        // идея правки... я пока не знаю.
                    {
                        double tmp = perchResult.coords[j] + i * (Xr_pool.coords[j] - perchResult.coords[j]) / deltapr;
                        //for (int k = 0; k < U.Count; k++)
                        //{
                        //    if(tmp U[k].Item1)
                        //}
                        perch.coords[j] = tmp;
                    }
                    I(perch);

                    if (perch.fitness < min)
                    {
                        inPool = perch;
                        min = perch.fitness;
                    }
                }
                Pool.Add(inPool);
            }

            //Sort(Pool);
            Pool = Pool.OrderBy(s => s.fitness).ToList();
            result = Pool[0];
        }

        //private void WrongCoord()
        //{
        //    for (int i = 0; i < NumFlocks; i++)
        //    {
        //        for (int j = 0; j < NumPerchInFlock; j++)
        //        {
        //            for (int k = 0; k < N_dim / 3; k++)
        //            {
        //                if (flock[i, j].coords[k] < 0 || flock[i, j].coords[k] > 4)
        //                {
        //                    //int yyy = 0;
        //                    throw new Exception();
        //                }
        //            }
        //        }
        //    }
        //}

        //private void WrongCoord()
        //{
        //    for (int i = 0; i < NumFlocks; i++)
        //    {
        //        for (int j = 0; j < NumPerchInFlock; j++)
        //        {
        //            if (flock[i, j].coords[0] < -1 || flock[i, j].coords[1] > 1 || flock[i, j].coords[0] > 1 || flock[i, j].coords[1] < -1)
        //            {
        //                //int yyy = 0;
        //                throw new Exception();
        //            }
        //
        //        }
        //    }
        //}
    }
}

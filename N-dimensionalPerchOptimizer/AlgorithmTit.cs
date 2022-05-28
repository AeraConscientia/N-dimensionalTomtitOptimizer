using System;
using System.Collections.Generic;
using System.Linq;

namespace N_dimensionalTomtitOptimizer
{
    public abstract class AlgorithmTit
    {
        public double R1, R2;
        public double thetta1, thetta2;
        public double L1, L2;

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


        uint I_CallsCount = 0;
        public Tit result;
        public int NP;
        public double alpha;
        public double gamma;
        public double lambda;
        public double eta;
        public double rho;
        public double c1, c2, c3;
        public double K;
        public double h;
        public double L;
        public double T;
        public double P;
        public double mu;
        public double eps;
        public Tit best;
        public List<Tit> individuals = new List<Tit>();            //Текущее множество синиц
        public List<Tit> search_tits = new List<Tit>();            //Массив лучших положений всех синиц после скачков (см. Шаг 2.6)
        public List<Tit> Pool = new List<Tit>();
        public List<Tit> memory;
        public double r;

        public Random random;

        public AlgorithmTit() { }

        public virtual void I(Tit tit, bool flag = false) 
        {
#if DEBUG
            I_CallsCount++;
#endif
        }
        public abstract void InitalPopulationGeneration();
        public abstract void GenerationAroundPool();
        public abstract void GenerationAroundBest();
        public abstract void BestCoordCorrect(Vector bestCoords, int k);
        /// <summary>Создание N-мерных коэф Леви</summary>
        /// <returns>Список коэффициентов для распределения Леви</returns>
        public abstract Vector Levy();
        public abstract void CheckBorders(Tit new_tit);
        public abstract void StohasticJump(Tit new_tit, double beta, int j);

        public Tit StartAlg(int NP, double alpha, double gamma, double lambda, double eta, double rho, double c1, double c2, double c3,
            double K, double h, double L, double P, double mu, double eps, int N_dim)
        {
            I_CallsCount = 0;       //Обнуляем число подсчетов

            //Шаг 1.1
            this.NP = NP;
            this.alpha = alpha;
            this.gamma = gamma;
            this.lambda = lambda;
            this.eta = eta;
            this.rho = rho;
            this.c1 = c1;
            this.c2 = c2;
            this.c3 = c3;
            this.K = K;
            this.h = h;
            this.L = L;
            this.P = P;
            this.mu = mu;
            this.eps = eps;
            this.N_dim = N_dim;
            r = 1;

            random = new Random();
            memory = new List<Tit>();

            InitalPopulationGeneration();       //Шаг 1.2

            for (int p = 0; p < P; p++)
            {
                //Диффузионный поиск со скачками
                int k = 0;                      //Шаг 2.1
                do
                {
                    individuals = individuals.OrderBy(t => t.fitness).ToList();     //Шаг 2.2
                    best = new Tit(individuals[0]);                                 //x^1,k

                    ProcessInfoAboutFlock();                                        //Шаг 2.3
                    SolveStohasticDiffEq();                                         //Шаг 2.4-2.6

                    if ((k >= K) || (Math.Pow(r, k) < eps))     //Шаг 2.7
                    {
                        Pool.Add(new Tit(memory.OrderBy(t => t.fitness).ToList()[0]));      //Шаг 2.7 K=k
                        break;
                    }

                    r *= gamma;
                    ++k;

                    // >>>>> Положение Вожака       Шаг 2.8
                    Vector bestCoords = new Vector(best.coords.dim);
                    for (int i = 0; i < best.coords.dim; i++)
                        bestCoords[i] = best.coords[i];

                    best.coords += (alpha / (k + 1)) * Levy();
                    BestCoordCorrect(bestCoords, k);
                    // <<<<<< Положение Вожака

                    I(best);
                    individuals[0] = new Tit(best);

                    //Шаг 2.9
                    GenerationAroundBest();

                } while (true);

                r = Math.Pow(eta, p);           //Шаг 3
                memory = new List<Tit>();

                individuals[0] = new Tit(Pool.OrderBy(t => t.fitness).ToList()[0]);
                I(individuals[0]);

                GenerationAroundPool();
            }

            result = new Tit(Pool.OrderBy(t => t.fitness).ToList()[0]);
            I(result, true);

#if DEBUG
            Console.WriteLine("------------------------------------------------");
            Console.WriteLine("Число подсчетов целевой функции: " + I_CallsCount);        //Вывод числа подсчета целевой функции
            Console.WriteLine("------------------------------------------------");
#endif
            return result;        //Шаг 4
        }

        public void ProcessInfoAboutFlock()
        {
            for (int j = 1; j < NP; j++)
            {
                if (individuals[j].fitness < individuals[j].best.fitness)
                    individuals[j].best = new Tit(individuals[j]);

                FindLocalBest(individuals[j]);
            }
        }

        public void SolveStohasticDiffEq()
        {
            search_tits = new List<Tit>();      

            for (int j = 1; j < NP; j++)
            {
                Tit current_tit = new Tit(individuals[j]);

                List<Tit> search = new List<Tit>(); //список всех скачков

                search.Add(individuals[j]);         //x^j,k (0)
                for (int l = 0; l < L; l++)
                {
                    Tit new_tit = new Tit(individuals[j].coords.dim);

                    double r1 = random.NextDouble();
                    double r2 = random.NextDouble();
                    double r3 = random.NextDouble();
                    double alpha1 = random.NextDouble();
                    double alpha2 = random.NextDouble();
                    double beta = random.NextDouble();
                    double ksi = Math.Sqrt(-2 * Math.Log(alpha1)) * Math.Cos(2 * Math.PI * alpha2);

                    Vector f = new Vector(c1 * r1 * (best.coords - current_tit.coords));
                    Vector sigma = new Vector(c2 * r2 * (individuals[j].best.coords - current_tit.coords) + c3 * r3 * (individuals[j].local_best.coords - current_tit.coords));

                    new_tit.coords = new Vector(current_tit.coords + h * f + Math.Sqrt(h) * sigma * ksi);
                    I(new_tit);

                    StohasticJump(new_tit, beta, j);

                    CheckBorders(new_tit);

                    I(new_tit);

                    search.Add(new_tit);   //x^j,k (l)
                    current_tit = new Tit(new_tit);
                }

                search = search.OrderBy(t => t.fitness).ToList();
                search_tits.Add(search[0]);
            }

            //>>>>>>>>Шаг 2.6
            search_tits.Add(best);
            search_tits = search_tits.OrderBy(t => t.fitness).ToList();
            memory.Add(search_tits[0]);
            best = search_tits[0];
            //<<<<<<<<
        }

        //Шаг 2.3
        public void FindLocalBest(Tit tit)
        {
            foreach (Tit item in individuals)
            {
                bool ok = true;
                for (int i = 0; i < item.coords.dim; i++)
                {
                    if (Math.Abs(item.coords[i] - tit.coords[i]) > rho)
                    {
                        ok = false;
                        break;
                    }
                }
                if (ok == true)
                    if (item.fitness < tit.local_best.fitness)
                        tit.local_best = new Tit(item);
            }
        }
    }
}

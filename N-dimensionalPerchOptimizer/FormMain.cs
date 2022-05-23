using System;
using System.Windows.Forms;
using System.IO;
using System.Diagnostics;
using System.Threading;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Drawing.Design;

namespace N_dimensionalPerchOptimizer
{
    public partial class FormMain : Form
    {
        /// <summary>Размерность задачи</summary>
        public int N_dim = 0;

        Graphics graphics = null;
        ErrorGraph errorGraph = null;
        public int task;


        private int MaxIteration = 0;
        private Perch resultBest;

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

        //private int population = 0;
        public int population = 0;

        public FormMain()
        {
            InitializeComponent();
            InitDataGridView();

            ToolStripMenuItem rereshProtocol = new ToolStripMenuItem("Очистить протокол");
            contextMenuStripProtocol.Items.AddRange(new[] { rereshProtocol});
            rereshProtocol.Click += buttonRefresh_Click;
            buttonProtocol.ContextMenuStrip = contextMenuStripProtocol;


            FileStream fs = new FileStream("protocol.txt", FileMode.Create, FileAccess.Write);
            fs.Close();
            fs = new FileStream("protocol.txt", FileMode.Append, FileAccess.Write); // деваааачки, какие костыли, я не могу T_T
            StreamWriter r = new StreamWriter(fs);
            r.Write(
                @"
| <><    <><    <><    <><    <><    <><    <><      ><>    ><>    ><>    ><>    ><>    ><>    ><> |
|                                                                                                  |
| <><                          Протокол применения метода стаи окуней                          ><> |
|              к задаче поиска оптимального управления и траектории дискретных систем              |
| <><                                                                                          ><> |
|__________________________________________________________________________________________________|
");
            r.Close();
            fs.Close();
        }

        /// <summary>Загрузка параметров аглоритма</summary>
        private void InitDataGridView()
        {
            dataGridView2.RowCount = 6;
            dataGridView2.Rows[0].Cells[0].Value = "Кол-во шагов до окончания движения";
            dataGridView2.Rows[0].Cells[1].Value = 100;

            dataGridView2.Rows[1].Cells[0].Value = "Максимальное количество итераций";
            dataGridView2.Rows[1].Cells[1].Value = 4;

            dataGridView2.Rows[2].Cells[0].Value = "Количество стай";
            dataGridView2.Rows[2].Cells[1].Value = 4;

            dataGridView2.Rows[3].Cells[0].Value = "Количество окуней в стае";
            dataGridView2.Rows[3].Cells[1].Value = 3;

            dataGridView2.Rows[4].Cells[0].Value = "Число перекоммутаций";
            dataGridView2.Rows[4].Cells[1].Value = 10;

            dataGridView2.Rows[5].Cells[0].Value = "Число шагов в перекоммутации";
            dataGridView2.Rows[5].Cells[1].Value = 5;


            dataGridView4.RowCount = 2;
            dataGridView4.Rows[0].Cells[0].Value = "Параметр распределения";
            dataGridView4.Rows[0].Cells[1].Value = (1.5).ToString();

            dataGridView4.Rows[1].Cells[0].Value = "Величина шага";
            dataGridView4.Rows[1].Cells[1].Value = (0.6).ToString();
        }

        private void LoadParams()
        {
            // Параметры окуней
            NumFlocks       = Convert.ToInt32(dataGridView2.Rows[2].Cells[1].Value);
            NumPerchInFlock = Convert.ToInt32(dataGridView2.Rows[3].Cells[1].Value);
            NStep           = Convert.ToInt32(dataGridView2.Rows[0].Cells[1].Value);
            MaxIteration    = Convert.ToInt32(dataGridView2.Rows[1].Cells[1].Value);
            PRmax           = Convert.ToInt32(dataGridView2.Rows[4].Cells[1].Value);
            deltapr         = Convert.ToInt32(dataGridView2.Rows[5].Cells[1].Value);

            // Для Леви
            lambda  = Convert.ToDouble(dataGridView4.Rows[0].Cells[1].Value);
            alfa    = Convert.ToDouble(dataGridView4.Rows[1].Cells[1].Value);
        }

        private async void button1_Click(object sender, EventArgs e)
        {
            task = tabControl2.SelectedIndex;

            buttonStartAlg.Enabled = false;
            buttonStartAlg.BackColor = System.Drawing.Color.Gainsboro;
            labelTimeStart.Text = (System.DateTime.Now.ToLongTimeString());

            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            LoadParams();
            AlgorithmPerch algPerch;

            object[] X;
            object[] X2;
            object[] X3;
            object[] U, U_2, U_3;

            switch (tabControl2.SelectedIndex) // Считывание параметров для задачи
            {
                case 0:
                    N_dim = Convert.ToInt32(numericUpDownN1.Value);
                    double U1_1 = Convert.ToDouble(textBoxU1_1.Text);       double U2_1 = Convert.ToDouble(textBoxU2_1.Text);

                    double x0_1 = Convert.ToDouble(textBoxX0_1.Text);

                    algPerch = new AlgorithmTask1(U1_1, U2_1, x0_1);

                    break;
                case 1:
                    N_dim = 3 * Convert.ToInt32(numericUpDownN2.Value);
                    double U11 = Convert.ToDouble(textBoxU11.Text);     double U12 = Convert.ToDouble(textBoxU12.Text);
                    double U21 = Convert.ToDouble(textBoxU21.Text);     double U22 = Convert.ToDouble(textBoxU22.Text);
                    double U31 = Convert.ToDouble(textBoxU31.Text);     double U32 = Convert.ToDouble(textBoxU32.Text);

                    double x00 = Convert.ToDouble(textBoxX11.Text);
                    double x11 = Convert.ToDouble(textBoxX22.Text);
                    double x22 = Convert.ToDouble(textBoxX33.Text);

                    algPerch = new AlgorithmTask2(U11, U12, U21, U22, U31, U32, x00, x11, x22);
                    break;
                case 2:
                    N_dim = Convert.ToInt32(numericUpDownN3.Value);
                    double U1_3 = Convert.ToDouble(textBoxU1_3.Text);   double U2_3 = Convert.ToDouble(textBoxU2_3.Text);

                    double x0_3 = Convert.ToDouble(textBoxX0_3.Text);

                    algPerch = new AlgorithmTask3(U1_3, U2_3, x0_3);
                    break;
                case 3:
                    N_dim = Convert.ToInt32(numericUpDownN4.Value);
                    double U1_4 = Convert.ToDouble(textBoxU1_4.Text);   double U2_4 = Convert.ToDouble(textBoxU2_4.Text);

                    double x0_4 = Convert.ToDouble(textBoxX0_4.Text);

                    algPerch = new AlgorithmTask4(U1_4, U2_4, x0_4);

                    break;
                case 4:
                    N_dim = 2 * Convert.ToInt32(numericUpDownN5.Value);
                    double U1_5 = Convert.ToDouble(textBoxU1_5.Text);   double U2_5 = Convert.ToDouble(textBoxU2_5.Text);
                
                    double x01_5 = Convert.ToDouble(textBoxX01_5.Text); double x02_5 = Convert.ToDouble(textBoxX02_5.Text);

                    algPerch = new AlgorithmTask5(U1_5, U2_5, x01_5, x02_5);
                
                    break;
                case 5:
                    N_dim = 2 * Convert.ToInt32(numericUpDownN6.Value);
                    double U1_6 = Convert.ToDouble(textBoxU1_6.Text);   double U2_6 = Convert.ToDouble(textBoxU2_6.Text);

                    double x01_6 = Convert.ToDouble(textBoxX01_6.Text); double x02_6 = Convert.ToDouble(textBoxX02_6.Text);

                    algPerch = new AlgorithmTask6(U1_6, U2_6, x01_6, x02_6);

                    break;
                case 6:
                    N_dim = 2 * Convert.ToInt32(numericUpDownN7.Value);
                    double U1_7 = Convert.ToDouble(textBoxU1_7.Text); double U2_7 = Convert.ToDouble(textBoxU2_7.Text);

                    double x01_7 = Convert.ToDouble(textBoxX01_7.Text); double x02_7 = Convert.ToDouble(textBoxX02_7.Text);

                    algPerch = new AlgorithmTask7(U1_7, U2_7, x01_7, x02_7);

                    break;
                case 7:
                    N_dim = 2 * Convert.ToInt32(numericUpDownN8.Value);
                    double U1_8 = Convert.ToDouble(textBoxU1_8.Text); double U2_8 = Convert.ToDouble(textBoxU2_8.Text);

                    double x01_8 = Convert.ToDouble(textBoxX01_8.Text); double x02_8 = Convert.ToDouble(textBoxX02_8.Text);

                    algPerch = new AlgorithmTask8(U1_8, U2_8, x01_8, x02_8);

                    break;
                case 8:
                    N_dim = Convert.ToInt32(numericUpDownN9.Value);
                    double U1_9 = Convert.ToDouble(textBoxU1_9.Text); double U2_9 = Convert.ToDouble(textBoxU2_9.Text);

                    double x0_9 = Convert.ToDouble(textBoxX0_9.Text);

                    double gamma_9 = Convert.ToDouble(textBoxGamma_9.Text);

                    algPerch = new AlgorithmTask9(U1_9, U2_9, x0_9, gamma_9);

                    break;
                default:
                    return;
            }
           // resultBest = 
           
            await Task.Run(() => algPerch.StartAlg(MaxIteration, NumFlocks, NumPerchInFlock, NStep, lambda, alfa, PRmax, deltapr, N_dim)); 

            Result result = Result.GetInstance();
            
            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                                            ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);

            FileStream fs = new FileStream("protocol.txt", FileMode.Append, FileAccess.Write);
            StreamWriter r = new StreamWriter(fs);
            r.Write(String.Format(
                @"
1. ПОСТАНОВКА ЗАДАЧИ
    Решаемая задача: Пример {0,3}", tabControl2.SelectedIndex + 1));
            switch (tabControl2.SelectedIndex) // первая часть записей протокола
            {
                case 0: // Пример 1
                case 2: // Пример 3
                case 3: // Пример 4
                case 8: // Пример 9
                    r.Write(String.Format(@"
    Количество точек разбиения (шагов): {0, 5}", N_dim));
                    switch (tabControl2.SelectedIndex)
                    {
                        case 0:
                            r.Write(String.Format(@"
    Ограничения на управление:          {1, 5:f1} <= u <= {2, 5:f1}
    Начальные условия:                           x ={3, 5:f1}", N_dim, Convert.ToDouble(textBoxU1_1.Text), Convert.ToDouble(textBoxU2_1.Text), Convert.ToDouble(textBoxX0_1.Text)));
                            break;
                        case 2:
                            r.Write(String.Format(@"
    Ограничения на управление:          {1, 5:f1} <= u <= {2, 5:f1}
    Начальные условия:                           x ={3, 5:f1}", N_dim, Convert.ToDouble(textBoxU1_3.Text), Convert.ToDouble(textBoxU2_3.Text), Convert.ToDouble(textBoxX0_3.Text)));
                            break;
                        case 3:
                            r.Write(String.Format(@"
    Ограничения на управление:          {1, 5:f1} <= u <= {2, 5:f1}
    Начальные условия:                           x ={3, 5:f1}", N_dim, Convert.ToDouble(textBoxU1_4.Text), Convert.ToDouble(textBoxU2_4.Text), Convert.ToDouble(textBoxX0_4.Text)));
                            break;
                        case 8:
                            r.Write(String.Format(@"
    Ограничения на управление:          {1, 5:f1} <= u <= {2, 5:f1}
    Начальные условия:                           x ={3, 5:f1}", N_dim, Convert.ToDouble(textBoxU1_9.Text), Convert.ToDouble(textBoxU2_9.Text), Convert.ToDouble(textBoxX0_9.Text)));
                            break;
                    }
                    
                    break;
                case 1: // Пример 2
                    r.Write(String.Format(
                                    @"
    Количество точек разбиения (шагов): {0, 5}
                                        {1, 5:f1} <= u1 <= {2, 5:f1}
    Ограничения на управление:          {3, 5:f1} <= u2 <= {4, 5:f1}
                                        {5, 5:f1} <= u3 <= {6, 5:f1}

    Начальные условия:                           x1 ={7, 5:f1}
                                                 x2 ={8, 5:f1}
                                                 x3 ={9, 5:f1}", 
                                    N_dim / 3, 
                                    Convert.ToDouble(textBoxU11.Text), Convert.ToDouble(textBoxU12.Text),
                                    Convert.ToDouble(textBoxU21.Text), Convert.ToDouble(textBoxU22.Text),
                                    Convert.ToDouble(textBoxU31.Text), Convert.ToDouble(textBoxU32.Text), 
                                    Convert.ToDouble(textBoxX11.Text), Convert.ToDouble(textBoxX22.Text), Convert.ToDouble(textBoxX33.Text)));
                    break;
                case 4: // Пример 5
                case 5: // Пример 6
                case 6: // Пример 7
                case 7: // Пример 8
                    r.Write(String.Format(@"
    Количество точек разбиения (шагов): {0, 5}", N_dim / 2));
                    switch (tabControl2.SelectedIndex)
                    {
                        case 4:
                            r.Write(String.Format(
                                    @"
    Ограничения на управление:          {0, 5:f1} <= u <= {1, 5:f1}

    Начальные условия:                           x1 ={2, 5:f1}
                                                 x2 ={3, 5:f1}",
                                    Convert.ToDouble(textBoxU1_5.Text), Convert.ToDouble(textBoxU1_5.Text),
                                    Convert.ToDouble(textBoxX01_5.Text), Convert.ToDouble(textBoxX02_5.Text)));
                            break;
                        case 5:
                            r.Write(String.Format(
                                    @"
    Ограничения на управление:          {0, 5:f1} <= u <= {1, 5:f1}

    Начальные условия:                           x1 ={2, 5:f1}
                                                 x2 ={3, 5:f1}",
                                    Convert.ToDouble(textBoxU1_6.Text), Convert.ToDouble(textBoxU1_6.Text),
                                    Convert.ToDouble(textBoxX01_6.Text), Convert.ToDouble(textBoxX02_6.Text)));
                            break;
                        case 6:
                            r.Write(String.Format(
                                    @"
    Ограничения на управление:          {0, 5:f1} <= u <= {1, 5:f1}

    Начальные условия:                           x1 ={2, 5:f1}
                                                 x2 ={3, 5:f1}",
                                    Convert.ToDouble(textBoxU1_7.Text), Convert.ToDouble(textBoxU2_7.Text),
                                    Convert.ToDouble(textBoxX01_7.Text), Convert.ToDouble(textBoxX02_7.Text)));
                            break;
                        case 7:
                            r.Write(String.Format(
                                    @"
    Ограничения на управление:          {0, 5:f1} <= u <= {1, 5:f1}

    Начальные условия:                           x1 ={2, 5:f1}
                                                 x2 ={3, 5:f1}",
                                    Convert.ToDouble(textBoxU1_8.Text), Convert.ToDouble(textBoxU1_8.Text),
                                    Convert.ToDouble(textBoxX01_8.Text), Convert.ToDouble(textBoxX02_8.Text)));
                            break;
                    }
                    
                    break;
            }
            r.Write(String.Format(@"
2. ПАРАМЕТРЫ МЕТОДА СТАИ ОКУНЕЙ
    Количество шагов до окончания
                        движения:       {0,5}
    Максимальное количество итераций:   {1,5}
    Количество стай:                    {2,5}
    Количество окуней в стае:           {3,5}
    Число перекоммутаций:               {4,5}
    Число шагов в перекоммутации:       {5,5}
    Параметр распределения Леви:        {6,5:f3}
    Величина шага:                      {7,5:f3}", NStep, MaxIteration, NumFlocks, NumPerchInFlock, PRmax, deltapr, lambda, alfa));

            r.Write(String.Format(@"
3. РЕЗУЛЬТАТЫ РАБОТЫ"));
            switch (tabControl2.SelectedIndex) // занесение в таблицу результатов
            {
                case 0:     // одномерный случай
                case 2:     // одномерный случай
                case 3:     // одномерный случай
                case 8:
                    X = new object[N_dim+1];
                    U = new object[N_dim];
                    for (int i = 0; i < N_dim; i++)
                    {
                        U[i] = result.U[i];
                    }
                    for (int i = 0; i < N_dim+1; i++)
                    {
                        X[i] = result.X[i];
                    }
                    dataGridViewX_separate.RowCount = 1;
                    dataGridViewX_separate.ColumnCount = N_dim + 1;
                    dataGridViewX_separate.Rows[0].SetValues(X);

                    dataGridViewU_separate.RowCount = 1;
                    dataGridViewU_separate.ColumnCount = N_dim;
                    
                    dataGridViewU_separate.Rows[0].SetValues(U);

                    
                    r.Write(String.Format(@"
    Оптимальное управление u*:
"));
                    for (int i = 0; i < N_dim; i++)     r.Write(String.Format(@" {0, 10:f5}", U[i]));    r.Write(String.Format("\r\n"));
                    r.Write(String.Format(@"
    Оптимальная траектория x*:
"));
                    for (int i = 0; i < N_dim+1; i++)   r.Write(String.Format(@" {0, 10:f5}", X[i]));    r.Write(String.Format("\r\n"));
                    break;
                case 4:     // двумерный случай
                case 5:
                case 6:
                case 7:
                    X   = new object[N_dim / 2 + 1];
                    X2  = new object[N_dim / 2 + 1];
                    U   = new object[N_dim / 2];
                    for (int i = 0; i < N_dim / 2 ; i++)
                    {
                        U[i] = result.U[i];
                    }
                    for (int i = 0; i < N_dim / 2 + 1; i++)
                    {
                        X[i] = result.X[i];
                        X2[i] = result.X2[i];
                    }
                    dataGridViewX_separate.Rows.Clear();
                    dataGridViewX_separate.RowCount = 2;
                    dataGridViewX_separate.ColumnCount = N_dim / 2 + 1;

                    dataGridViewX_separate.Rows[0].SetValues(X);
                    dataGridViewX_separate.Rows[1].SetValues(X2);

                    dataGridViewU_separate.Rows.Clear();
                    dataGridViewU_separate.RowCount = 1;
                    dataGridViewU_separate.ColumnCount = N_dim / 2;
                    dataGridViewU_separate.Rows[0].SetValues(U);

                    dataGridViewX_separate.Rows[0].DefaultCellStyle.Format = "n5";
                    dataGridViewX_separate.Rows[1].DefaultCellStyle.Format = "n5";
                    dataGridViewU_separate.Rows[0].DefaultCellStyle.Format = "n5";
                    
                    r.Write(String.Format(@"
    Оптимальное управление u*:
")); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 2; i++) r.Write(String.Format(@" {0, 10:f5}", U[i])); r.Write(String.Format("\r\n"));
                    r.Write(String.Format(@"
    Оптимальная траектория x*:
")); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 2 + 1; i++) r.Write(String.Format(@" {0, 10:f5}", X[i])); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 2 + 1; i++) r.Write(String.Format(@" {0, 10:f5}", X2[i])); r.Write(String.Format("\r\n"));
                    break;

                case 1:     // трехмерный случай
                    X = new object[N_dim + 3];
                    X2 = new object[N_dim + 3];
                    X3 = new object[N_dim + 3];
                    
                    U = new object[N_dim];
                    object[] U_0 = new object[N_dim / 3];
                    U_2 = new object[N_dim / 3];
                    U_3 = new object[N_dim / 3];
                    for (int i = 0; i < N_dim; i++)
                    {
                        U[i] = result.U[i];
                    }
                    Array.Copy(U, 0, U_0, 0, N_dim / 3);
                    Array.Copy(U,    N_dim/3, U_2, 0, N_dim/3);
                    Array.Copy(U, 2* N_dim/3, U_3, 0, N_dim/3);
                    for (int i = 0; i < N_dim/3 + 1; i++)
                    {
                        X[i] = result.X[i];
                        X2[i] = result.X2[i];
                        X3[i] = result.X3[i];
                    }
                    dataGridViewX_separate.Rows.Clear();
                    dataGridViewX_separate.RowCount = 3;
                    dataGridViewX_separate.ColumnCount = N_dim/3 + 1;

                    dataGridViewX_separate.Rows[0].SetValues(X);
                    dataGridViewX_separate.Rows[1].SetValues(X2);
                    dataGridViewX_separate.Rows[2].SetValues(X3);

                    dataGridViewU_separate.Rows.Clear();
                    dataGridViewU_separate.RowCount = 3;
                    dataGridViewU_separate.ColumnCount = N_dim / 3;
                    dataGridViewU_separate.Rows[0].SetValues(U_0);
                    dataGridViewU_separate.Rows[1].SetValues(U_2);
                    dataGridViewU_separate.Rows[2].SetValues(U_3);

                    dataGridViewX_separate.Rows[1].DefaultCellStyle.Format = "n5";
                    dataGridViewX_separate.Rows[2].DefaultCellStyle.Format = "n5";
                    dataGridViewU_separate.Rows[1].DefaultCellStyle.Format = "n5";
                    dataGridViewU_separate.Rows[2].DefaultCellStyle.Format = "n5";
                    
                    r.Write(String.Format(@"
    Оптимальное управление u*:
")); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 3; i++)     r.Write(String.Format(@" {0, 10:f5}", U_0[i])); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 3; i++)     r.Write(String.Format(@" {0, 10:f5}", U_2[i])); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 3; i++)     r.Write(String.Format(@" {0, 10:f5}", U_3[i])); r.Write(String.Format("\r\n"));
                    r.Write(String.Format(@"
    Оптимальная траектория x*:
")); r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 3 + 1; i++) r.Write(String.Format(@" {0, 10:f5}", X[i]));   r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 3 + 1; i++) r.Write(String.Format(@" {0, 10:f5}", X2[i]));  r.Write(String.Format("\r\n"));
                    for (int i = 0; i < N_dim / 3 + 1; i++) r.Write(String.Format(@" {0, 10:f5}", X3[i]));  r.Write(String.Format("\r\n"));
                    break;
                default:
                    return;
            }
            dataGridViewX_separate.Rows[0].DefaultCellStyle.Format = "n5";
            dataGridViewU_separate.Rows[0].DefaultCellStyle.Format = "n5";

            labelMinI.Text = result.fitness.ToString();

            r.Write(String.Format(@"
    Значение функционала: {0, 5:f1}", labelMinI.Text)); r.Write(String.Format("\r\n"));
            r.Write(String.Format(@"
Время расчета (ч:м:с): {0}", elapsedTime));
            r.Write(@"
------------------------------------------------------------------------
");
            r.Close();
            fs.Close();

            if (errorGraph == null && tabControl2.SelectedIndex == 6)
            {
                errorGraph = new ErrorGraph(6, N_dim);
            }
            if (tabControl2.SelectedIndex == 6 && errorGraph.IsDisposed)
            {
                errorGraph = new ErrorGraph(6, N_dim);
            }
            if (tabControl2.SelectedIndex == 6)
                errorGraph.UpdateErrorGraph(6, N_dim);
            //task = tabControl2.SelectedIndex;
            if (graphics == null)
                switch (tabControl2.SelectedIndex)
                {
                    case 0:
                    case 2:
                    case 3:
                    case 8:
                        graphics = new Graphics(1, tabControl2.SelectedIndex);
                        break;
                    case 1:
                        graphics = new Graphics(3, tabControl2.SelectedIndex);
                        break;
                    case 4:
                    case 5:
                    case 6:
                    case 7:
                        graphics = new Graphics(2, tabControl2.SelectedIndex);
                        break;
                    default:
                        break;
                }
            if (graphics.IsDisposed)
                switch (tabControl2.SelectedIndex)
                {
                    case 0:
                    case 2:
                    case 3:
                    case 8:
                        graphics = new Graphics(1, tabControl2.SelectedIndex);
                        break;
                    case 1:
                        graphics = new Graphics(3, tabControl2.SelectedIndex);
                        break;
                    case 4:
                    case 5:
                    case 6:
                    case 7:
                        graphics = new Graphics(2, tabControl2.SelectedIndex);
                        if (tabControl2.SelectedIndex == 6)
                            errorGraph = new ErrorGraph(6, N_dim);
                        break;
                    default:
                        break;
                }
            switch (tabControl2.SelectedIndex)
                {
                    case 0:
                    case 2:
                    case 3:
                    case 8:
                        graphics.UpdateGraph(1, tabControl2.SelectedIndex);
                        break;
                    case 1:
                        graphics.UpdateGraph(3, tabControl2.SelectedIndex);
                        break;
                    case 4:
                    case 5:
                    case 6:
                    case 7:
                        graphics.UpdateGraph(2, tabControl2.SelectedIndex);
                        //if (tabControl2.SelectedIndex == 6)
                    break;
                    default:
                        break;
                }
            
            graphics.Show();

            if (tabControl2.SelectedIndex == 6)
                if (errorGraph.IsDisposed == false && errorGraph != null)
                    errorGraph.Show();
            buttonStartAlg.BackColor = System.Drawing.Color.DarkSeaGreen;
            buttonStartAlg.Enabled = true;
            labelTimeStart.Text = "---";

        }

        /// <summary>Вызов протокола</summary>
        private void buttonProtocol_Click(object sender, EventArgs e)
        {
            Process.Start("protocol.txt");
        }

        /// <summary>Очистка протокола</summary>
        private void buttonRefresh_Click(object sender, EventArgs e)
        {
            FileStream fs = new FileStream("protocol.txt", FileMode.Create, FileAccess.Write);
            StreamWriter r = new StreamWriter(fs);
            r.Write(
                @"
| <><    <><    <><    <><    <><    <><    <><      ><>    ><>    ><>    ><>    ><>    ><>    ><> |
|                                                                                                  |
| <><                          Протокол применения метода стаи окуней                          ><> |
|              к задаче поиска оптимального управления и траектории дискретных систем              |
| <><                                                                                          ><> |
|__________________________________________________________________________________________________|
");
            r.Close();
            fs.Close();
        }

        private void CleanAll(object sender, EventArgs e)
        {
            dataGridViewX_separate.Rows.Clear();    dataGridViewU_separate.Rows.Clear();

            dataGridViewX_separate.RowCount     = 1;
            dataGridViewX_separate.ColumnCount  = 1;
            
            dataGridViewU_separate.RowCount     = 1;
            dataGridViewU_separate.ColumnCount  = 1;

            labelMinI.Text = "---";
            labelTimeStart.Text = "---";

        }

        private void FormMain_KeyUp(object sender, KeyEventArgs e)
        {
            if (e.KeyCode == Keys.Enter)
                buttonStartAlg.PerformClick();
        }
    }
}

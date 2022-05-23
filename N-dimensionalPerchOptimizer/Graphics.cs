using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace N_dimensionalPerchOptimizer
{
    public partial class Graphics : Form
    {
        int TaskNumber;
        int Dim;

        public Graphics(int dim, int task)
        {
            InitializeComponent();
            Result result = Result.GetInstance();

            TaskNumber = task + 1;
            Dim = dim;


            CleanGraphics();

            switch (dim)
            {
                case 1: // случай, когда рисуются только два графика (1,3 примеры)  - одномерный случай
                    for (int i = 0; i < result.X.Count; i++)
                        chartX_1.Series[0].Points.AddXY(i, result.X[i]);

                    for (int i = 0; i < result.U.Count; i++)
                        chartU_1.Series[0].Points.AddXY(i, result.U[i]);
                    break;
                case 2: // случай, когда рисуются три графика (4 пример)         - двумерный случай, 2 траектории, 1 управление
                    for (int i = 0; i < result.X.Count; i++)
                    {
                        chartX_1.Series[0].Points.AddXY(i, result.X[i]);
                        chartX_2.Series[0].Points.AddXY(i, result.X2[i]);
                    }

                    for (int i = 0; i < result.U.Count; i++)
                        chartU_1.Series[0].Points.AddXY(i, result.U[i]);
                    break;
                case 3: // случай, когда рисуются шесть графиков (2 пример)         - трехмерный случай
                    for (int i = 0; i < result.X.Count; i++)
                    {
                        chartX_1.Series[0].Points.AddXY(i, result.X[i]);
                        chartX_2.Series[0].Points.AddXY(i, result.X2[i]);
                        chartX_3.Series[0].Points.AddXY(i, result.X3[i]);
                    }

                    int N_dim = result.U.Count;
                    List<double> U_2 = new List<double>(N_dim / 3);
                    List<double> U_3 = new List<double>(N_dim / 3);

                    for (int i = 0; i < N_dim / 3; i++)
                    {
                        U_2.Add(result.U[i + N_dim/3]);
                        U_3.Add(result.U[i + 2 * N_dim / 3]);
                    }

                    for (int i = 0; i < result.U.Count / 3; i++)
                    {
                        chartU_1.Series[0].Points.AddXY(i, result.U[i]);
                        chartU_2.Series[0].Points.AddXY(i, U_2[i]);
                        chartU_3.Series[0].Points.AddXY(i, U_3[i]);
                    }
                        
                    break;
                default:
                    break;
            }
        }

        public void RoundChart()
        {
            //chartX_1.;
        }

        public void UpdateGraph(int dim, int task)
        {
            Result result = Result.GetInstance();
            CleanGraphics();
            TaskNumber = task + 1;
            Dim = dim;

            switch (dim)
            {
                case 1:
                    CleanTabControlHeader(1);
                    for (int i = 0; i < result.X.Count; i++)    chartX_1.Series[0].Points.AddXY(i, result.X[i]);

                    for (int i = 0; i < result.U.Count; i++)    chartU_1.Series[0].Points.AddXY(i, result.U[i]);
                    //chartX_1.ChartAreas[0].AxisX.Interval = 1;
                    //tabControlU.TabPages[1].Dispose();

                    break;
                case 2: // тут 2 траектории и 1 управление
                    CleanTabControlHeader(2);
                    for (int i = 0; i < result.X.Count; i++)
                    {
                        chartX_1.Series[0].Points.AddXY(i, result.X[i]);
                        chartX_2.Series[0].Points.AddXY(i, result.X2[i]);
                    }
                    for (int i = 0; i < result.U.Count; i++)    chartU_1.Series[0].Points.AddXY(i, result.U[i]);
                    break;
                case 3:
                    CleanTabControlHeader(3);
                    for (int i = 0; i < result.X.Count; i++)
                    {
                        chartX_1.Series[0].Points.AddXY(i, result.X[i]);
                        chartX_2.Series[0].Points.AddXY(i, result.X2[i]);
                        chartX_3.Series[0].Points.AddXY(i, result.X3[i]);
                    }
                    //tabControlU.TabPages[1].CreateControl();
                    int N_dim = result.U.Count;
                    List<double> U_2 = new List<double>(N_dim / 3);
                    List<double> U_3 = new List<double>(N_dim / 3);

                    for (int i = 0; i < N_dim / 3; i++)
                    {
                        U_2.Add(result.U[i + N_dim / 3]);
                        U_3.Add(result.U[i + 2 * N_dim / 3]);
                    }

                    for (int i = 0; i < U_2.Count; i++)
                    {
                        chartU_1.Series[0].Points.AddXY(i, result.U[i]);
                        chartU_2.Series[0].Points.AddXY(i, U_2[i]);
                        chartU_3.Series[0].Points.AddXY(i, U_3[i]);
                    }
                    break;
                default:
                    break;
            }
        }

        /// <summary>Очистка самих графиков</summary>
        private void CleanGraphics()
        {
            chartX_1.Series[0].Points.Clear(); chartX_2.Series[0].Points.Clear(); chartX_3.Series[0].Points.Clear();
            chartU_1.Series[0].Points.Clear(); chartU_2.Series[0].Points.Clear(); chartU_3.Series[0].Points.Clear();
        }

        /// <summary>Названия вкладок графиков</summary>
        /// <param name="dim">на самом деле, тип задачи</param>
        private void CleanTabControlHeader(int dim)
        {
            switch (dim)
            {
                case 1:
                    tabControlX.TabPages[0].Text = "X"; tabControlX.TabPages[1].Text = "----------"; tabControlX.TabPages[2].Text = "----------";
                    tabControlU.TabPages[0].Text = "U"; tabControlU.TabPages[1].Text = "----------"; tabControlU.TabPages[2].Text = "----------";
                    break;
                case 2:
                    tabControlX.TabPages[0].Text = "X1"; tabControlX.TabPages[1].Text = "X2"; tabControlX.TabPages[2].Text = "----------";
                    tabControlU.TabPages[0].Text = "U"; tabControlU.TabPages[1].Text = "----------"; tabControlU.TabPages[2].Text = "----------";
                    break;
                case 3:
                    tabControlX.TabPages[0].Text = "X1"; tabControlX.TabPages[1].Text = "X2"; tabControlX.TabPages[2].Text = "X3";
                    tabControlU.TabPages[0].Text = "U1"; tabControlU.TabPages[1].Text = "U2"; tabControlU.TabPages[2].Text = "U3";
                    break;
                default:
                    break;
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {

        }

        private void buttonError_Click(object sender, EventArgs e)
        {
           // if (errorGraph == null)
           //     switch (tabControl2.SelectedIndex)
           //     {
           //         case 6:
           //             errorGraph = new ErrorGraph(6, N_dim);
           //             break;
           //         default:
           //             break;
           //     }
           // if (errorGraph.IsDisposed)
           //     switch (tabControl2.SelectedIndex)
           //     {
           //         case 6:
           //             errorGraph = new ErrorGraph(6, N_dim);
           //             break;
           //         default:
           //             break;
           //     }
           // switch (tabControl2.SelectedIndex)
           // {
           //     case 6:
           //         errorGraph.UpdateErrorGraph(6, N_dim);
           //         break;
           //     default:
           //         break;
           // }
           //
           // if (errorGraph != null)
           //     errorGraph.Show();
        }

        /// <summary>Сохранение всех графиков в папку</summary>
        private void buttonSaveImg_Click(object sender, EventArgs e)
        {

            switch (Dim)
            {
                case 1:
                    chartX_1.SaveImage($"SavedTask{TaskNumber}_X1.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartU_1.SaveImage($"SavedTask{TaskNumber}_U1.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    break;
                case 2:
                    chartX_1.SaveImage($"SavedTask{TaskNumber}_X1.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartX_2.SaveImage($"SavedTask{TaskNumber}_X2.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartU_1.SaveImage($"SavedTask{TaskNumber}_U1.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    break;
                case 3:
                    chartX_1.SaveImage($"SavedTask{TaskNumber}_X1.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartX_2.SaveImage($"SavedTask{TaskNumber}_X2.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartX_3.SaveImage($"SavedTask{TaskNumber}_X3.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartU_1.SaveImage($"SavedTask{TaskNumber}_U1.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartU_2.SaveImage($"SavedTask{TaskNumber}_U2.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    chartU_3.SaveImage($"SavedTask{TaskNumber}_U3.png", System.Windows.Forms.DataVisualization.Charting.ChartImageFormat.Png);
                    break;
                default:
                    break;
            }
        }
    }
}

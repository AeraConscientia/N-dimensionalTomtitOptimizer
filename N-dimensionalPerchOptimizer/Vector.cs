using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace N_dimensionalPerchOptimizer
{
    public class Vector
    {
        public double[] vector;
        public int dim;

        public Vector(int dim)
        {
            vector = new double[dim];
            this.dim = dim;
        }

        public Vector(double x, double y)
        {
            vector = new double[2];
            vector[0] = x;
            vector[1] = y;
            dim = 2;
        }

        public Vector()
        {
        }

        public static Vector operator *(Vector vector, double value)
        {
            Vector tmp = new Vector(vector.dim);
            for (int i = 0; i < vector.dim; i++)
                tmp[i] = vector[i] * value;
            return tmp;
        }
        public static Vector operator *(double value, Vector vector)
        {
            return vector * value;
        }

        public static Vector operator *(Vector vector1, Vector vector2)
        {
            Vector tmp = new Vector(vector1.dim);
            for (int i = 0; i < vector1.dim; i++)
                tmp[i] = vector1[i] * vector2[i];
            return tmp;
        }

        public static Vector Abs(Vector vector)
        {
            Vector tmp = new Vector(vector.dim);
            for (int i = 0; i < vector.dim; i++)
                tmp[i] = Math.Abs(vector[i]);
            return tmp;
        }

        public static Vector operator -(Vector vec1, Vector vec2)
        {
            return vec1 + (-1 * vec2);
        }

        public static Vector operator +(Vector vec1, Vector vec2)
        {
            Vector tmp = new Vector(vec1.dim);
            for (int i = 0; i < vec1.dim; i++)
                tmp[i] = vec1[i] + vec2[i];
            return tmp;
        }

        public static Vector operator +(Vector vec, double val)
        {
            Vector tmp = new Vector(vec.dim);
            for (int i = 0; i < vec.dim; i++)
                tmp[i] = vec[i] + val;
            return tmp;
        }

        public static Vector operator +(double val, Vector vec)
        {
            return vec + val;
        }
        public static Vector operator /(Vector vec, double val)
        {
            return vec * (1.0 / val);
        }

        public double this[int index]
        {
            get { return vector[index]; }
            set { vector[index] = value; }
        }
    }
}

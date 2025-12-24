using System;
using System.IO;
using System.Collections.Generic;

// no significant .NET dependencies

namespace MatrixSingularValueDecomposition
{
  internal class SVDProgram
  {
    static void Main(string[] args)
    {
      Console.WriteLine("\nBegin SVD using scratch C# ");
      double[][] A = new double[6][];
      A[0] = new double[] { 1, 2, 3, 4, 5 };
      A[1] = new double[] { 0, -3, 5, -7, 9 };
      A[2] = new double[] { 2, 0, -2, 0, -2 };
      A[3] = new double[] { 4, -1, 5, 6, 1 };
      A[4] = new double[] { 3, 6, 8, 2, 2 };
      A[5] = new double[] { 5, -2, 4, -4, 3 };

      Console.WriteLine("\nSource (tall) matrix: ");
      MatShow(A, 1, 6);

      Console.WriteLine("\nComputing SVD ");
      double[][] U;
      double[] s;
      double[][] Vh;
      SVD.MatDecomp(A, out U, out s, out Vh);
      Console.WriteLine("Done ");

      Console.WriteLine("\nU = ");
      MatShow(U, 6, 11);

      Console.WriteLine("\ns = ");
      for (int i = 0; i < s.Length; ++i)
        Console.Write(s[i].ToString("F6").PadLeft(11));
      Console.WriteLine("");

      Console.WriteLine("\nVh = ");
      MatShow(Vh, 6, 11);

      double[][] B = new double[3][];
      B[0] = new double[] { 3, 1, 1, 0, 5 };
      B[1] = new double[] { -1, 3, 1, -2, 4 };
      B[2] = new double[] { 0, 2, 2, 1, -3 };

      Console.WriteLine("\n============================ ");

      Console.WriteLine("\nSource (wide) matrix: ");
      MatShow(B, 1, 6);

      Console.WriteLine("\nComputing SVD ");
      SVD.MatDecomp(B, out U, out s, out Vh);
      Console.WriteLine("Done ");

      Console.WriteLine("\nU = ");
      MatShow(U, 6, 11);

      Console.WriteLine("\ns = ");
      for (int i = 0; i < s.Length; ++i)
        Console.Write(s[i].ToString("F6").PadLeft(11));
      Console.WriteLine("");

      Console.WriteLine("\nVh = ");
      MatShow(Vh, 6, 11);

      Console.WriteLine("\nEnd demo ");
      Console.ReadLine();
    } // Main

    // helpers for Main()

    // ------------------------------------------------------

    static void MatShow(double[][] M, int dec, int wid)
    {
      for (int i = 0; i < M.Length; ++i)
      {
        for (int j = 0; j < M[0].Length; ++j)
        {
          double v = M[i][j];
          Console.Write(v.ToString("F" + dec).
            PadLeft(wid));
        }
        Console.WriteLine("");
      }
    }

    // ------------------------------------------------------

    static double[][] MatLoad(string fn, int[] usecols,
      char sep, string comment)
    {
      // ex: double[][] A =
      // MatLoad("data.txt", new int[] { 0,1,2,3 }, ',', "#");
      List<double[]> result =
        new List<double[]>();
      string line = "";
      FileStream ifs = new FileStream(fn, FileMode.Open);
      StreamReader sr = new StreamReader(ifs);
      while ((line = sr.ReadLine()) != null)
      {
        if (line.StartsWith(comment) == true)
          continue;
        string[] tokens = line.Split(sep);
        List<double> lst = new List<double>();
        for (int j = 0; j < usecols.Length; ++j)
          lst.Add(double.Parse(tokens[usecols[j]]));
        double[] row = lst.ToArray();
        result.Add(row);
      }
      sr.Close(); ifs.Close();
      return result.ToArray();
    }

  } // class Program

  // --------------------------------------------------------

  public class SVD
  {
    public static void MatDecomp(double[][] A,
      out double[][] U, out double[] s, out double[][] Vt,
      double tol = 1.0e-12)
    {
      int m = A.Length; int n = A[0].Length;
      if (m >= n)
      {
        double[][] Q; double[][] R;
        MatDecompQR(A, out Q, out R);

        double[][] RtR = MatProduct(MatTranspose(R), R);

        double[] evals; double[][] V;
        MatDecompEigen(RtR, out evals, out V);

        int[] sortIdxs = VecReversed(ArgSort(evals));
        evals = VecSortedUsingIdxs(evals, sortIdxs);
        V = MatSortedUsingIdxs(V, sortIdxs);

        double[] ss = VecSqrt(VecClip(evals, 0.0));
        int[] idxs = VecIndicesWhereNonNegative(ss, tol);
        double[] ssNonZero = VecRemoveZeros(ss, idxs);
        double[][] VNonZero = MatRemoveCols(V, idxs);

        double[][] tmp1 = MatProduct(R, VNonZero);
        double[][] tmp2 = MatDivVector(tmp1, ssNonZero);
        double[][] UU = MatProduct(Q, tmp2);

        U = UU;
        s = ssNonZero;
        Vt = MatTranspose(VNonZero);
      } // m >= n
      else
      {
        // use the transpose trick
        double[][] dummyU;
        double[] dummyS;
        double[][] dummyVt;
        MatDecomp(MatTranspose(A),
          out dummyU, out dummyS, out dummyVt);
        U = MatTranspose(dummyVt);
        s = dummyS;
        Vt = MatTranspose(dummyU);
      } // m < n

    } // MatDecomp()

    // ------------------------------------------------------
    // primary helpers: MatDecompQR() and MatDecompEigen()
    // ------------------------------------------------------

    private static void MatDecompQR(double[][] A,
      out double[][] Q, out double[][] R)
    {
      // reduced QR decomp using Householder reflections
      int m = A.Length; int n = A[0].Length;
      double[][] QQ = MatIdentity(m);  // working Q
      double[][] RR = MatCopy(A);  // working R

      for (int k = 0; k < n; ++k)
      {
        double[] x = MatColToEnd(RR, k); // RR[k:, k]
        double normx = VecNorm(x);
        if (Math.Abs(normx) < 1.0e-12) continue;

        double[] v = VecCopy(x);
        double sign;
        if (x[0] >= 0.0) sign = 1.0; else sign = -1.0;
        v[0] += sign * normx;
        double normv = VecNorm(v);
        for (int i = 0; i < v.Length; ++i)
          v[i] /= normv;

        // apply reflection to R
        double[][] tmp1 = 
          MatRowColToEnd(RR, k); // RR[k:, k:]
        double[] tmp2 = VecMatProduct(v, tmp1);
        double[][] tmp3 = VecOuter(v, tmp2);
        double[][] B = MatScalarMult(2.0, tmp3);
        MatSubtractCorner(RR, k, B);  // R[k:, k:] -= B

        // accumulate Q
        double[][] tmp4 = 
          MatExtractAllRowsColToEnd(QQ, k); // QQ[:, k:]
        double[] tmp5 = MatVecProduct(tmp4, v);
        double[][] tmp6 = VecOuter(tmp5, v);
        double[][] C = MatScalarMult(2.0, tmp6);
        MatSubtractRows(QQ, k, C); // QQ[:, k:] -= C
      } // for k

      // return parts of QQ and RR
      Q = MatExtractFirstCols(QQ, n); // Q[:, :n]
      R = MatExtractFirstRows(RR, n); // R[:n, :]
    } // MatDecompQR()

    // ------------------------------------------------------

    private static void MatDecompEigen(double[][] A,
      out double[] evals, out double[][] Evecs,
      double tol = 1.0e-12, int maxIter = 5000)
    {
      // eigen-decomp of symmetric matrix using QR
      int m = A.Length; int n = A[0].Length;
      double[][] AA = MatCopy(A);  // don't destroy A
      double[][] V = MatIdentity(n);

      int iter = 0;
      while (iter < maxIter)
      {
        double[][] Q;
        double[][] R;
        MatDecompQR(AA, out Q, out R);
        double[][] Anew = MatProduct(R, Q);
        V = MatProduct(V, Q);

        // # check convergence (off-diagonal norm)
        double[] tmp1 = MatDiag(Anew);
        double[][] tmp2 = MatCreateDiag(tmp1);
        double[][] tmp3 = MatSubtract(Anew, tmp2);
        double norm = MatNorm(tmp3);
        if (norm < tol)
        {
          AA = Anew; break;
        }
        AA = Anew;
        ++iter;
      } // iter

      if (iter == maxIter)
        Console.WriteLine("Warn: exceeded maxIter ");

      evals = MatDiag(AA);
      Evecs = V;
    } // MatDecompEigen()

    // ------------------------------------------------------
    // 32 secondary helpers
    // ------------------------------------------------------

    private static int[] ArgSort(double[] vec)
    {
      int n = vec.Length;
      int[] idxs = new int[n];
      for (int i = 0; i < n; ++i)
        idxs[i] = i;

      double[] dup = new double[n];
      for (int i = 0; i < n; ++i)
        dup[i] = vec[i]; // don't modify vec

      // special C# overload
      Array.Sort(dup, idxs);  // sort idxs based on dup vals
      return idxs;
    }

    // ------------------------------------------------------

    private static double[][] MatMake(int nr, int nc)
    {
      double[][] result = new double[nr][];
      for (int i = 0; i < nr; ++i)
        result[i] = new double[nc];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatIdentity(int n)
    {
      double[][] result = MatMake(n, n);
      for (int i = 0; i < n; ++i)
        result[i][i] = 1.0;
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatCopy(double[][] A)
    {
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, n);
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
          result[i][j] = A[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static double[] MatDiag(double[][] A)
    {
      // return diag elements of A as a vector
      int m = A.Length; int n = A[0].Length; // m == n
      double[] result = new double[m];
      for (int i = 0; i < m; ++i)
        result[i] = A[i][i];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatCreateDiag(double[] v)
    {
      // create a matrix using v as diagonal elements
      int n = v.Length;
      double[][] result = MatMake(n, n);
      for (int i = 0;i < n; ++i)
        result[i][i] = v[i];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatTranspose(double[][] A)
    {
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(n, m);  // note
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
          result[j][i] = A[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatSubtract(double[][] A,
      double[][] B)
    {
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, n);
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
          result[i][j] = A[i][j] - B[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static double[] MatColToEnd(double[][] A,
      int k)
    {
      // last part of col k, from [k,k] to end
      int m = A.Length; int n = A[0].Length;
      double[] result = new double[m - k];
      for (int i = 0; i < m - k; ++i)
        result[i] = A[i + k][k];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatRowColToEnd(double[][] A,
      int k)
    {
      // block of A, row k to end, col k to end
      // A[k:, k:]
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m-k, n-k);
      for (int i = 0; i < m - k; ++i)
        for (int j = 0; j < n - k; ++j)
          result[i][j] = A[i + k][j + k];
      return result;
    }

    // ------------------------------------------------------

    private static double[] MatVecProduct(double[][] A,
      double[] v)
    {
      // A * v
      int m = A.Length; int n = A[0].Length;
      double[] result = new double[m];
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
          result[i] += A[i][j] * v[j];
      return result;
    }

    private static double[] VecMatProduct(double[] v,
      double[][] A)
    {
      // v * A
      int m = A.Length; int n = A[0].Length;
      double[] result = new double[n];
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < m; ++i)
          result[j] += v[i] * A[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatScalarMult(double x,
      double[][] A)
    {
      // x * A element-wise
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, n);
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
          result[i][j] = x * A[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static void MatSubtractCorner(double[][] A,
      int k, double[][] C)
    {
      // subtract elements in C from lower right A
      // A[k:, k:] -= C
      int m = A.Length; int n = A[0].Length;
      for (int i = 0; i < m - k; ++i)
        for (int j = 0; j < n - k; ++j)
          A[i + k][j + k] -= C[i][j];

      return;  // modify A in-place 
    }

    // ------------------------------------------------------

    private static void MatSubtractRows(double[][] A,
      int k, double[][] C)
    {
      // A[:, k:] -= C
      int m = A.Length; int n = A[0].Length;
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n-k; ++j)
          A[i][j+k] -= C[i][j];

      return;  // modify A in-place
    }

    // ------------------------------------------------------

    private static double[][] 
      MatExtractAllRowsColToEnd(double[][] A, int k)
    {
      // A[:, k:]
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, n - k);
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n - k; ++j)
          result[i][j] = A[i][j + k];
      return result;
    }

    // -----------------------------------------------------

    private static double[][] 
      MatExtractFirstCols(double[][] A, int k)
    {
      // A[:, :n]
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, k);
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < k; ++j)
          result[i][j] = A[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] 
      MatExtractFirstRows(double[][] A, int k)
    {
      // A[:n, :]
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(k, n);
      for (int i = 0; i < k; ++i)
        for (int j = 0; j < n; ++j)
          result[i][j] = A[i][j];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatProduct(double[][] A,
        double[][] B)
    {
      int aRows = A.Length;
      int aCols = A[0].Length;
      int bRows = B.Length;
      int bCols = B[0].Length;
      if (aCols != bRows)
        throw new Exception("Non-conformable matrices");

      double[][] result = MatMake(aRows, bCols);

      for (int i = 0; i < aRows; ++i) // each row of A
        for (int j = 0; j < bCols; ++j) // each col of B
          for (int k = 0; k < aCols; ++k)
            result[i][j] += A[i][k] * B[k][j];

      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatRemoveCols(double[][] A,
      int[] idxs)
    {
      // remove the columns specified in idxs
      int m = A.Length; int n = A[0].Length;
      int countNonzeos = 0;
      for (int i = 0; i < idxs.Length; ++i)
        if (idxs[i] == 1) ++countNonzeos;

      double[][] result = MatMake(m, countNonzeos);
      int k = 0; // destination col
      for (int j = 0; j < n; ++j)
      {
        if (idxs[j] == 0) continue; // skip this col
        for (int i = 0; i < m; ++i)
        {
          result[i][k] = A[i][j];
        }
        ++k;
      }
      return result;
    }

    // ------------------------------------------------------

    private static double[][] MatDivVector(double[][] A,
      double[] v)
    {
      // v must be all non-zero
      // divide each row of A by v, element-wise
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, n);
      for (int i = 0; i < m; ++i)
        for (int j= 0; j < n; ++j)
          result[i][j] = A[i][j] / v[j];
      return result;
    }

    // ------------------------------------------------------

    private static double MatNorm(double[][] A)
    {
      // treat A as one big vector
      int m = A.Length; int n = A[0].Length;
      double sum = 0.0;
      for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
          sum += A[i][j] * A[i][j];
      return Math.Sqrt(sum); 
    }

    private static double VecNorm(double[] v)
    {
      // standard vector norm
      double sum = 0.0;
      for (int i = 0; i < v.Length; ++i)
        sum += v[i] * v[i];
      return Math.Sqrt(sum);
    }

    // ------------------------------------------------------

    private static double[][] 
      MatSortedUsingIdxs(double[][] A, int[] idxs)
    {
      // arrange columns of A using order in idxs
      int m = A.Length; int n = A[0].Length;
      double[][] result = MatMake(m, n);
      for (int i = 0; i < m; ++i)
      {
        for (int j = 0; j < n; ++j)
        {
          int srcCol = idxs[j];
          result[i][j] = A[i][srcCol];
        }

      }
      return result;
    }

    // ------------------------------------------------------

    private static double[] 
      VecSortedUsingIdxs(double[] v, int[] idxs)
    {
      // arrange values in v using order specified in idxs
      int n = v.Length;
      double[] result = new double[n];
      for (int i = 0; i < n; ++i)
        result[i] = v[idxs[i]];
      return result;
    }

    // ------------------------------------------------------

    private static double[] VecCopy(double[] v)
    {
      double[] result = new double[v.Length];
      for (int i = 0; i < v.Length; ++i)
        result[i] = v[i];
      return result;
    }

    // ------------------------------------------------------

    private static double[][] VecOuter(double[] v1,
      double[] v2)
    {
      // vector outer product
      double[][] result = MatMake(v1.Length, v2.Length);
      for (int i = 0; i < v1.Length; ++i)
        for (int j = 0; j < v2.Length; ++j)
          result[i][j] = v1[i] * v2[j];
      return result;
    }

    // ------------------------------------------------------

    private static int[] VecReversed(int[] v)
    {
      // return vector of v elements in reversed order
      int n = v.Length;
      int[] result = new int[n];
      for (int i = 0; i < n; ++i)
        result[i] = v[n-1-i];
      return result;
    }

    // ------------------------------------------------------

    private static double[] VecClip(double[] v, double minVal)
    {
      // used to make sure no negative values
      int n = v.Length;
      double[] result = new double[n];
      for (int i = 0; i < n; ++i)
      {
        if (v[i] < minVal)
          result[i] = minVal;
        else
          result[i] = v[i];
      }
      return result;
    }

    // ------------------------------------------------------

    private static double[] VecSqrt(double[] v)
    {
      // sqrt each element in v
      // assumes v has no negative values
      int n = v.Length;
      double[] result = new double[n];
      for (int i = 0; i < n; ++i)
        result[i] = Math.Sqrt(v[i]);
      return result;
    }

    // ------------------------------------------------------

    private static int[] 
      VecIndicesWhereNonNegative(double[] v, double tol)
    {
      // indices of non-negative values (within tolerance)
      int n = v.Length;
      int[] result = new int[n];
      for (int i = 0; i < n; ++i)
      {
        if (Math.Abs(v[i]) > tol)
          result[i] = 1;
      }
      return result;
    }

    // ------------------------------------------------------

    private static double[] VecRemoveZeros(double[] v,
      int[] idxs)
    {
      // if idx = 1, val is non-zero, if 0 val is near zero
      int n = v.Length;
      List<double> lst = new List<double>();
      for (int i = 0; i < n; ++i)
      {
        if (idxs[i] == 1)
          lst.Add(v[i]);
      }
      return lst.ToArray();
    }

    // ------------------------------------------------------

  } // class SVD

} // ns

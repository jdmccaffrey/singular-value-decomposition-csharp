# singular-value-decomposition-csharp
Function to compute singular value decomposition (SVD) of a matrix, implemented using C#

Function MatDecomp() is contained in a class SVD. Results are returned as out parameters. Example call:
```
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
double[][] U; double[] s; double[][] Vh;
SVD.MatDecomp(A, out U, out s, out Vh);
Console.WriteLine("Done ");

Console.WriteLine("\nU = "); MatShow(U, 6, 11);
Console.WriteLine("\ns = ");
for (int i = 0; i < s.Length; ++i)
  Console.Write(s[i].ToString("F6").PadLeft(11));
Console.WriteLine("");

Console.WriteLine("\nVh = ");  MatShow(Vh, 6, 11);
```
Output:
```
Begin SVD using scratch C#

Source (tall) matrix:
   1.0   2.0   3.0   4.0   5.0
   0.0  -3.0   5.0  -7.0   9.0
   2.0   0.0  -2.0   0.0  -2.0
   4.0  -1.0   5.0   6.0   1.0
   3.0   6.0   8.0   2.0   2.0
   5.0  -2.0   4.0  -4.0   3.0

Computing SVD
Done

U =
   0.319267   0.295047   0.358160   0.465731   0.652040
   0.625260  -0.614481   0.189087   0.174205  -0.138090
  -0.129989   0.032917  -0.342862  -0.160369   0.572690
   0.284056   0.494504  -0.490845   0.494308  -0.381366
   0.486033   0.495867   0.263713  -0.648653  -0.103790
   0.416300  -0.209426  -0.638701  -0.248874   0.267559

s =
  15.967661  12.793149   6.297366   5.706879   2.478679

Vh =
   0.296544   0.035215   0.708796  -0.130799   0.625557
   0.217254   0.416871   0.261754   0.803400  -0.255056
  -0.745282   0.555722  -0.030757   0.039094   0.365040
  -0.187161  -0.609726  -0.196995   0.579569   0.467438
   0.523820   0.379983  -0.623971   0.003540   0.438033
```

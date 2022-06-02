public class Solution {

    public static class Rational {

        private int num, denom;

        public Rational(double d) {
            String s = String.valueOf(d);
            int digitsDec = s.length() - 1 - s.indexOf('.');
            int denom = 1;
            for (int i = 0; i < digitsDec; i++) {
                d *= 10;
                denom *= 10;
            }

            int num = (int) Math.round(d);
            this.num = num;
            this.denom = denom;
        }

        public Rational(int num, int denom) {
            this.num = num;
            this.denom = denom;
        }

        public String toString() {
            return String.valueOf(num) + "/" + String.valueOf(denom);
        }

        public static void main(String[] args) {
            System.out.println(new Rational(1.5));
        }
    }

    static double[][] multiplyMatrix(double[][] A, double[][] B) {

        int row1 = A.length;
        int col1 = A[1].length;
        int row2 = B.length;
        int col2 = B[1].length;

        int i, j, k;

        /*
        // Check if multiplication is Possible
        if (row2 != col1) {

            System.out.println(
                    "\nMultiplication Not Possible");
            return;
        }
        */

        // Matrix to store the result
        // The product matrix will
        // be of size row1 x col2
        double[][] C = new double[row1][col2];

        // Multiply the two matrices
        for (i = 0; i < row1; i++) {
            for (j = 0; j < col2; j++) {
                for (k = 0; k < row2; k++)
                    C[i][j] += A[i][k] * B[k][j];
            }
        }

        return C;

    }

    // Method to carry out the partial-pivoting Gaussian elimination.
    // Here index[] stores pivoting order.

    public static void gaussian(double[][] a, int[] index)
    {
        int n = index.length;
        double[] c = new double[n];

        // Initialize the index
        for (int i=0; i<n; ++i)
            index[i] = i;

        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i)
        {
            double c1 = 0;
            for (int j=0; j<n; ++j)
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j)
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i)
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1)
                {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i)
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }

    public static double[][] invert(double[][] a) {
        int n = a.length;
        double[][] x = new double[n][n];
        double [][] b = new double[n][n];
        int[] index = new int[n];
        for (int i=0; i<n; ++i)
            b[i][i] = 1;

        // Transform the matrix into an upper triangle
        gaussian(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                            -= a[index[j]][i]*b[index[i]][k];

        // Perform backward substitutions
        for (int i=0; i<n; ++i)
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j)
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k)
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }
    public static int[] solution(int[][] m) {
        // Your code here
        return null;
    }
}
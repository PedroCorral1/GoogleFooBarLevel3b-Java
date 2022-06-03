import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class Solution {

    static int findOverallGCD(Integer[] resultList) {
        int result = resultList[0];

        for(int i=1; i < resultList.length; i++){

            result = findLCM(resultList[i], result);

        }

        return result;

    }

    static int findGCD(int a, int b) {
        if (a == 0)
            return b;

        return findGCD(b % a, a);
    }

    // function to calculate lcm of
    // two numbers.
    static int findLCM(int a, int b) {
        return (a * b) / findGCD(a, b);
    }

    static class Rational {

        int num, denom;

        Rational(double value) {

            this.num = toFractionPos(new BigDecimal(String.valueOf(value)))[0];
            this.denom = toFractionPos(new BigDecimal(String.valueOf(value)))[1];

        }

        static int[] toFractionPos(BigDecimal x) {
            String[] parts = x.toString().split("\\.");
            BigDecimal den = BigDecimal.TEN.pow(parts[1].length()); // denominator
            BigDecimal num = (new BigDecimal(parts[0]).multiply(den)).add(new BigDecimal(parts[1])); // numerator

           return reduceFraction(num.intValue(), den.intValue());

        }

        static int[] reduceFraction(int num, int den) {
            int gcd = BigInteger.valueOf(num).gcd(BigInteger.valueOf(den)).intValue(); // greatest
            // common
            // divisor
            int[] rf = { num / gcd, den / gcd };
            return rf;
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
        HashMap<Integer, Integer> rowTotals = new HashMap<>();
        int currKey = 0;

        for (int[] row : m) {
            int sum = 0;
            for (int val : row) {
                sum += val;
            }

            rowTotals.put(currKey, sum);
            currKey += 1;
        }

        HashMap<Integer, Integer[]> masterMatrix = new HashMap<>();

        for (int i = 0; i < m.length; i++) {

            masterMatrix.put(i, Arrays.stream(m[i]).boxed().toArray(Integer[]::new));

        }

        HashMap<Integer, int[]> iMatrix = new HashMap<>();
        HashMap<Integer, int[]> oMatrix = new HashMap<>();
        HashMap<Integer, Double[]> rMatrix = new HashMap<>();
        HashMap<Integer, Double[]> qMatrix = new HashMap<>();
        HashMap<Integer, int[]> iQMatrix = new HashMap<>();

        /*

        */

        //prep iMatrix
        int numAbsStates = 0;

        for (int key : rowTotals.keySet()) {

            if (rowTotals.get(key) == 0) {

                numAbsStates++;

            }

        }

        for (int row : rowTotals.keySet()) {
            if (rowTotals.get(row) == 0) {
                iMatrix.put(row, new int[numAbsStates]);
            }
        }

        int valCount = 0;
        for (int key : iMatrix.keySet()) {
            iMatrix.get(key)[valCount] = 1;
            valCount++;
        }

        //prep rMatrix
        ArrayList<Integer> rKeys = new ArrayList<>();
        ArrayList<Integer> qKeys = new ArrayList<>();

        for (Integer key : masterMatrix.keySet()) {

            if (iMatrix.containsKey(key)) {

                rKeys.add(key);

            } else {

                qKeys.add(key);

            }

        }

        for (Integer key : masterMatrix.keySet()) {

            Double[] keyValues = new Double[iMatrix.size()];

            if (!iMatrix.containsKey(key)) {

                for (int i = 0; i < rKeys.size(); i++) {

                    keyValues[i] = masterMatrix.get(key)[rKeys.get(i)].doubleValue() / rowTotals.get(key);

                }

                rMatrix.put(key, keyValues);

            }

        }


        //prepare Q Matrix

        for (Integer key : masterMatrix.keySet()) {

            Double[] keyValues = new Double[masterMatrix.size() - iMatrix.size()];

            if (!iMatrix.containsKey(key)) {

                for (int i = 0; i < qKeys.size(); i++) {

                    keyValues[i] = masterMatrix.get(key)[qKeys.get(i)].doubleValue() / rowTotals.get(key);

                }

                qMatrix.put(key, keyValues);

            }

        }

        //prepare iQ Matrix
        for (int row : qMatrix.keySet()) {
            iQMatrix.put(row, new int[qMatrix.size()]);
        }

        int valCountIQ = 0;
        for (int key : iQMatrix.keySet()) {
            iQMatrix.get(key)[valCountIQ] = 1;
            valCountIQ++;
        }

        //prepare F Matrix

        HashMap<Integer, Double[]> fNIMatrix = new HashMap<>();

        for (Integer iQVal : iQMatrix.keySet()) {

            int[] iQValsList = iQMatrix.get(iQVal);
            Double[] qValsList = qMatrix.get(iQVal);
            Double[] subResult = new Double[qValsList.length];

            for (int i = 0; i < iQValsList.length; i++) {

                subResult[i] = iQValsList[i] - qValsList[i];

            }

            fNIMatrix.put(iQVal, subResult);

        }
            //Convert fNIMatrix & rMatrix to double[][] to invert & multiply

            double[][] fMatrix = new double[fNIMatrix.size()][fNIMatrix.size()];
            double[][] rMatrixPost = new double[rMatrix.size()][rMatrix.size()];

            for (int i = 0; i < fNIMatrix.size(); i++) {

                fMatrix[i] = Arrays.stream(fNIMatrix.get(i)).mapToDouble(j -> j).toArray();

            }

            for (int k = 0; k < rMatrix.size(); k++) {

                rMatrixPost[k] = Arrays.stream(rMatrix.get(k)).mapToDouble(j -> j).toArray();

            }

            fMatrix = invert(fMatrix);

            double[][] FRMatrix = multiplyMatrix(fMatrix, rMatrixPost);

            double[] result = FRMatrix[0];

            ArrayList<Integer> resultNums = new ArrayList<>();
            ArrayList<Integer> resultDenoms = new ArrayList<>();

            for (double value : result) {

                Rational fractionValue = new Rational(value);

                resultNums.add(fractionValue.num);
                resultDenoms.add(fractionValue.denom);

            }

            int lcd = findOverallGCD(resultDenoms.toArray(new Integer[0]));


            for (int i = 0; i < resultNums.size(); i++) {

                int newNum = resultNums.get(i);

                newNum *= (lcd / resultDenoms.get(i));

                resultNums.set(i, newNum);
            }

            resultNums.add(lcd);


        return resultNums.stream().mapToInt(i -> i).toArray();
    }

    public static void main(String[] args) {

        int[][] testArray = {{0,1,0,0,0,1}, {4,0,0,3,2,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};

        System.out.println(Arrays.toString(solution(testArray)));

    }

}
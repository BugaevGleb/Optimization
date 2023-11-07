
template<typename T>
void readVector(vector<T> &v, int n) {
    for (int i = 0; i < n; ++i) {
        T t;
        cin >> t;
        v.push_back(t);
    }
}


double calculateDeterminantUsingREF(Matrix &a) {

    int n = a.getN();
    SquareMatrix res{a};
    int numberOfPermutation = 0;
    for (int column = 0; column < n; ++column) {
        for (int row = 0; row < n; ++row) {
            // Direct way
            if (column < row) {
                int rowNumberWithMaximumValue = res.findRowNumberWithMaximumAbsolutValueInColumnBelow(column);
                if (rowNumberWithMaximumValue != column) {
                    PermutationMatrix p{n, rowNumberWithMaximumValue, column};
                    res = p * res;
                    numberOfPermutation++;
                }
                if (res[row][column] != 0.0) {
                    EliminationMatrix e1{res, row, column};
                    EliminationMatrix e2{res, row, column};
                    res = e1 * res;
                }

            }

        }

    }
    double det = 1.0;
    for (int i = 0; i < res.getN(); ++i) {
        det = det * res[i][i];
    }
    if (numberOfPermutation % 2 != 0) det = det * (-1);
    return det;
}

bool checkStandardFormVectorB(vector<double> &b) {
    for (double &i: b) {
        if (i < 0) {
            return false;
        }
    }
    return true;
}

class Utils {
public:
    static Matrix *transpose(Matrix &A) {
        Matrix *res = new Matrix(A.getM(), A.getN());
        for (int i = 0; i < A.getN(); ++i) {
            for (int j = 0; j < A.getM(); ++j) {
                (*res)[j][i] = A[i][j];
            }
        }

        return res;
    }

    static SquareMatrix *calculateInverse(Matrix &a) {

        int n = a.getN();
        IdentityMatrix *augmentedMatrix = new IdentityMatrix{n};
        SquareMatrix res{a};
        for (int column = 0; column < n; ++column) {
            for (int row = 0; row < n; ++row) {
                // Direct way
                if (column < row) {
                    int rowNumberWithMaximumValue = res.findRowNumberWithMaximumAbsolutValueInColumnBelow(column);
                    if (rowNumberWithMaximumValue != column) {
                        PermutationMatrix p{n, rowNumberWithMaximumValue, column};
                        res = p * res;
                        *augmentedMatrix = p * *augmentedMatrix;
                    }
                    if (res[row][column] != 0.0) {
                        EliminationMatrix e1{res, row, column};
                        EliminationMatrix e2{res, row, column};
                        res = e1 * res;
                        *augmentedMatrix = e2 * (*augmentedMatrix);

                    }
                }

            }

        }

        for (int row = n - 1; row >= 0; --row) {
            for (int column = n - 1; column >= 0; --column) {
                // Direct way

                if (column < row) {


                    int rowNumberWithMaximumValue = res.findRowNumberWithMaximumAbsolutValueInColumnBelow(column);
                    if (rowNumberWithMaximumValue != column) {
                        PermutationMatrix p{n, rowNumberWithMaximumValue, column};
                        res = p * res;
                        *augmentedMatrix = p * *augmentedMatrix;
                    }
                    if (res[column][row] != 0.0) {
                        EliminationMatrix e1{res, column, row};
                        EliminationMatrix e2{res, column, row};
                        res = e1 * res;
                        *augmentedMatrix = e2 * *augmentedMatrix;
                    }
                }

            }

        }

        // Normalization
        for (int column = 0; column < n; ++column) {
            for (int row = 0; row < n; ++row) {
                (*augmentedMatrix)[column][row] = (*augmentedMatrix)[column][row] / res[column][column];

            }
            res[column][column] = 1.0;


        }

        return augmentedMatrix;
    }

    static Matrix *getMatrixOfOnes(int size) {
        Matrix *res = new Matrix(size, 1);
        for (int i = 0; i < size; ++i) {
            (*res)[i][0] = 1;
        }

        return res;
    }

    static double getNorm(Matrix &A) {
        double res = 0;
        for (int i = 0; i < A.getN(); ++i) {
            res += (A[i][0] * A[i][0]);
        }
        return sqrt(res);
    }

};

void solve(vector<double> c, Matrix &A, vector<double> &b, vector<double> x, double alpha, bool isMinimize) {
    ColumnVector C = ColumnVector(c);


    ColumnVector X = ColumnVector(x);
    ColumnVector v = ColumnVector(x);
    int i = 0;
    while (true) {

        if(i>10000){
            cout<<"The problem does not have solution!"<<endl;
            exit(0);
        }
        v = X;
        DiagonalMatrix D = DiagonalMatrix(X);

        Matrix AA = A * D;

        Matrix CC = D * C;
        Matrix *AA_t = Utils::transpose(AA);
        Matrix F = AA * (*AA_t);

        double det = calculateDeterminantUsingREF(F);

        if (isnan(det) || det==0) {
            cout << "The method is not applicable!" << endl;
           return;
        }
        for (int i = 0; i < F.getN(); ++i) {
            for (int j = 0; j < F.getM(); ++j) {
                if (isnan(F[i][j])) {
                    cout << "The method is not applicable!" << endl;
                    return;
                }
            }
        }

        Matrix *FI = Utils::calculateInverse(F);
        Matrix H = *AA_t * (*FI);
        Matrix Q = H * AA;
        IdentityMatrix I = IdentityMatrix(x.size());

        Matrix P = I - Q;

        Matrix cp = P * CC;

        double nu = abs(cp.getMinimalNegativeValue());

        Matrix *ones = Utils::getMatrixOfOnes(x.size());

        Matrix y = *ones + (cp * (alpha / nu));
        Matrix yy = D * y;
        X = yy;
        Matrix diff = yy - v;

        if (Utils::getNorm(diff) < 0.00001) {
            break;
        }

        i++;

        free(AA_t);
        free(ones);
        free(FI);

    }

    cout << "Result:\n";
    double functionResult;
    for (int j = 0; j < c.size(); ++j) {
        functionResult = functionResult + X[j] * c[j];
    }
    if(isMinimize) functionResult = functionResult*(-1);

    for (int j = 0; j < c.size(); ++j) {
        cout<<"x"<<(j+1)<<" = "<<X[j]<<endl;
    }

    cout << "Value of the function: " << functionResult << endl;
    cout << "Number of iterations: " << i << endl;

}


vector<double> checkLinearityAndConvertVector(vector<string> v) {
    vector<double> f;
    for (int i = 0; i < size(v); i++) {
        try {
            f.push_back(stod(v[i]));
        }
        catch (...) {
            cout << "The method is not applicable!\n";
            exit(0);
        }
    }
    return f;
}



class Matrix {
protected:
    int n;
    int m;
    double **matrix;


public:

    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        // Allocate memory for the matrix
        matrix = (double **) malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            matrix[i] = (double *) malloc(m * sizeof(double));
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = 0.0;
            }
        }
    }

    // Copy constructor
    Matrix(Matrix &a) : Matrix(a.getN(), a.getM()) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = a[i][j];
            }
        }
    }


    // Destructor
    ~Matrix() {
        for (int i = 0; i < n; ++i) {
            //  delete matrix[i];
        }
        // delete matrix;
    }

    int getN() { return n; }

    int getM() { return m; }

    double *operator[](int index) {
        return matrix[index];
    }

    vector<double> getColumn(int column) {
        vector<double> res(getN());
        for (int i = 0; i < getN(); ++i) {
            for (int j = 0; j < getM(); ++j) {
                if (j == column) {
                    res[i] = matrix[i][j];
                }
            }
        }
        return res;
    }

    // Return the index of row (below current row) with maximum absolute value in current column

    int findRowNumberWithMaximumAbsolutValueInColumnBelow(int column) {
        double maxValue = matrix[column][column];
        int resRow = column;
        for (int row = column + 1; row < n; ++row) {
            if (abs(matrix[row][column]) > abs(maxValue)) {
                maxValue = abs(matrix[row][column]);
                resRow = row;
            }
        }
        return resRow;
    }

    double getMinimalNegativeValue() {
        double res = DBL_MIN;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (matrix[i][j] < 0)
                    res = min(res, matrix[i][j]);
            }
        }
        return res;
    }

    friend istream &operator>>(istream &stream, Matrix &matrix1) {
        for (int i = 0; i < matrix1.getN(); ++i) {
            for (int j = 0; j < matrix1.getM(); ++j) {
                stream >> matrix1[i][j];
            }
        }
        return stream;

    }


    friend ostream &operator<<(ostream &stream, Matrix &matrix1) {

        for (int i = 0; i < matrix1.n; ++i) {
            for (int j = 0; j < matrix1.m; ++j) {
                if (fabs(matrix1[i][j]) < epsilon) matrix1[i][j] = 0.0;
                stream << matrix1[i][j] << " ";
            }
            stream << endl;
        }
        return stream;
    }

    Matrix operator-(Matrix &input) {
        Matrix resMatrix(n, m);
        if (this->n == input.n && this->m == input.m) {

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    resMatrix[i][j] = this->matrix[i][j] - input[i][j];
                }
            }
        } else {
            cout << "Error: the dimensional problem occurred. Please check the dimensions of the input" << endl;
            exit(1);

        }
        return resMatrix;
    }

    Matrix operator+(Matrix &input) {
        Matrix resMatrix(n, m);
        if (this->n == input.n && this->m == input.m) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    resMatrix[i][j] = input[i][j] + this->matrix[i][j];
                }
            }
        } else {
            cout << "Error: the dimensional problem occurred. Please check the dimensions of the input" << endl;
            exit(1);
        }
        return resMatrix;
    }

    Matrix &operator*(Matrix &secondMatrix) {
        auto resMatrix = new Matrix{n, secondMatrix.m};
        if (m == secondMatrix.n) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < secondMatrix.m; ++j) {
                    double res = 0.0;
                    for (int k = 0; k < m; ++k) {
                        res = res + matrix[i][k] * secondMatrix[k][j];

                    }
                    (*resMatrix)[i][j] = res;

                }
            }
        } else {
            cout << "Error: the dimensional problem occurred. Please check the dimensions of the input" << endl;
            exit(1);
        }
        return *resMatrix;
    }

    Matrix &operator*(double c) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = matrix[i][j] * c;

            }
        }
        return *this;
    }
};


class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix(Matrix &a) : Matrix(a) {

    }

    int getSize() {
        return n;
    }

    using Matrix::operator=;
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) matrix[i][j] = 1.0;
                else matrix[i][j] = 0.0;
            }
        }
    }

    using Matrix::operator=;

};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(SquareMatrix &a, int i, int j) : IdentityMatrix(a.getSize()) {
        matrix[i][j] = -a[i][j] / a[j][j];

    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int n, int row1, int row2) : IdentityMatrix(n) {

        swap(matrix[row1], matrix[row2]);

    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {

    }

    ColumnVector(Matrix &a) : Matrix(a) {

    }

    ColumnVector(vector<double> a) : Matrix(a.size(), 1) {

        for (int i = 0; i < a.size(); ++i) {
            matrix[i][0] = a[i];

        }
    }

    // Returns the length of a ColumnVector
    double computeNorm() {
        double res = 0.0;
        for (int i = 0; i < n; ++i) {
            res += (*this)[i] + (*this)[i];
        }
        return ::sqrt(res);
    }

    // Overloaded operators
    double &operator[](int index) {
        return matrix[index][0];
    }

    using Matrix::operator=;
    using Matrix::operator*;
    using Matrix::operator+;

};

class DiagonalMatrix : public Matrix {
public:
    DiagonalMatrix(ColumnVector &a) : Matrix(a.getN(), a.getN()) {
        int n = a.getN();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    matrix[i][j] = a[i];
                } else {
                    matrix[i][j] = 0;
                }
            }
        }
    }

};
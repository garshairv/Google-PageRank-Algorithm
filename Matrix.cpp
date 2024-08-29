//
// Created by Misam Ibrahimi on 2021-10-05.
//

#include "Matrix.hpp"
#define PRECISION 2
#define PERCENTAGE 100

constexpr double INITIAL_VALUE = 0.0;
constexpr int INITIAL_SIZE_OF_COL_ROWS = 1;
constexpr double TOLERANCE = 0.001;
constexpr double WALK_PROB = 0.85;


using namespace std;

Matrix::Matrix() : colSize(INITIAL_SIZE_OF_COL_ROWS), rowSize(INITIAL_SIZE_OF_COL_ROWS) {
    vector <double> vec;
    vec.push_back(INITIAL_VALUE);
    matrix.push_back(vec);
}

Matrix::Matrix(const int n) : colSize(n), rowSize(n) {
    if (n < 0) {
        throw invalid_argument( "Values must be greater than 0!" );
    }
    matrix.resize(n);
    for (int i = 0; i < n; i++) {
        matrix[i].resize(n);
        for (int j = 0; j < n; j++) {
            matrix[i][j] = INITIAL_VALUE;
        }
    }
}

Matrix::Matrix(int r, int c) : colSize(c), rowSize(r) {
    if (r < 0 || c < 0) {
        throw invalid_argument("Column and Row values must be greater than 0!");
    }
    matrix.resize(r, vector<double>(c,0));
    for (int i =0; i < rowSize; i++) {
        for (int j = 0; j < colSize; j++) {
          matrix[i][j] = INITIAL_VALUE;
        }
    }
}

Matrix::Matrix(vector<double> values) {
    if (sqrt(values.size()) != (floor(sqrt(values.size())))) {
        throw invalid_argument("Size of the vector is not an integer square root!");
    }
    double size = sqrt(values.size());
    colSize = floor(size);
    rowSize = floor(size);

    matrix.resize(size);
    int cnt = 0;
    for (int i = 0; i < size; i++) {
        matrix[i].resize(size);
        for (int j = 0; j < size; j++) {
            matrix[i][j] = values[cnt];
            cnt++;
        }
    }
}


Matrix::Matrix(const Matrix &otherMatrix) : colSize(otherMatrix.getColSize()), rowSize(otherMatrix.getRowSize()) {
    matrix = otherMatrix.matrix;
    rowSize = otherMatrix.getRowSize();
    colSize = otherMatrix.getColSize();
}

int Matrix::getColSize() const {return colSize;}
int Matrix::getRowSize() const {return rowSize;}

void Matrix::setValue(int r, int c, double val) {
    if (r < 0.0 || c < 0.0) {
        throw invalid_argument("Rows and columns must be greater than 0!");
    }
    if (c > getColSize() || r > getRowSize()) {
        throw invalid_argument("Rows and columns must not exceed the size of the matrix!");
    }
    matrix[r][c] = val;
}

double Matrix::getValue(int r, int c) const{
    if (r < 0.0 || c < 0.0) {
        throw invalid_argument("rows and columns must be greater than 0!");
    }
    if (r > getRowSize() || c > getColSize()) {
        throw invalid_argument("Rows and columns must not exceed the size of the matrix!");
    }
    return matrix[r][c];
}

void Matrix::clear() {
    for (int i = 0; i < rowSize; i++) {
        for (int j = 0; j < colSize; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

Matrix::~Matrix() = default;

ostream &operator<<(ostream &os, const Matrix &obj) {
    for (int i = 0; i < obj.getRowSize(); i++) {
        for (int j = 0; j < obj.getColSize(); j++) {
            os << obj.getValue(i, j) << '\t';
        }
        cout << " " << endl;
    }
    return os;
}

bool operator==(const Matrix &lhs, const Matrix &rhs) {
    if (!(lhs.getRowSize() == rhs.getRowSize() && lhs.getColSize() == rhs.getColSize())){
        return false;
    }

    for (int i = 0; i < lhs.getRowSize(); i++) {
        for (int j = 0; j < lhs.getColSize(); j++) {
            if (abs(lhs.getValue(i, j) - rhs.getValue(i, j)) > TOLERANCE) {
                return false;
            }
        }
    }
    return true;
}

bool operator!=(const Matrix &lhs, const Matrix &rhs) {
    if (!(lhs.getRowSize() == rhs.getRowSize() && lhs.getColSize() == rhs.getColSize())){
        return true;
    }

    for (int i = 0; i < lhs.getRowSize(); i++) {
        for (int j = 0; j < lhs.getColSize(); j++) {
            if (abs(lhs.getValue(i, j) - rhs.getValue(i, j)) > TOLERANCE) {
                return true;
            }
        }
    }
    return false;
}

Matrix &Matrix::operator++() {
    for (int i = 0; i < this->getRowSize(); ++i) {
        for (int j = 0; j < this->getColSize(); ++j) {
            this->setValue(i, j, (getValue(i, j) + 1));
        }
    }
    return *this;
}

Matrix Matrix::operator++(int) {
    Matrix temp = *this;
    operator++();
    return temp;
}

Matrix &Matrix::operator--() {
    for (int i = 0; i < this->getRowSize(); ++i) {
        for (int j = 0; j < this->getColSize(); ++j) {
            this->setValue(i, j, (getValue(i, j) - INITIAL_SIZE_OF_COL_ROWS));
        }
    }
    return *this;
}

Matrix Matrix::operator--(int) {
    Matrix temp = *this;
    operator--();
    return temp;
}

void swapMatrix(Matrix &first, Matrix &second) {
    swap(first.matrix, second.matrix);
    swap(first.colSize, second.colSize);
    swap(first.rowSize, second.rowSize);
}

Matrix &Matrix::operator=(Matrix other) {
    swapMatrix(*this, other);
    return *this;
}

Matrix &Matrix::operator+=(const Matrix &rhs) {
    if (this->colSize != rhs.colSize && this->rowSize != rhs.rowSize) {
        throw invalid_argument("Both matrices must be the same size to do '+='!");
    }
    for (int i = 0; i < this->rowSize; i++) {
        for (int j = 0; j < this->colSize; j++) {
            setValue(i, j, getValue(i, j) + rhs.getValue(i, j));
        }
    }
    return *this;
}

Matrix operator+(Matrix lhs, const Matrix &rhs) {
    lhs += rhs;
    return lhs;
}

Matrix &Matrix::operator-=(const Matrix &rhs) {
    if (this->colSize != rhs.colSize && this->rowSize != rhs.rowSize) {
        throw invalid_argument("Both matrices must be the same size to do '+='!");
    }
    for (int i = 0; i < this->rowSize; i++) {
        for (int j = 0; j < this->colSize; j++) {
            setValue(i, j, getValue(i, j) - rhs.getValue(i, j));
        }
    }
    return *this;
}

Matrix operator-(Matrix lhs, const Matrix &rhs) {
    lhs -= rhs;
    return lhs;
}

Matrix& Matrix::operator*=(const Matrix &rightMatrix) {
    if (colSize != rightMatrix.rowSize) {
        throw invalid_argument("First matrix's number of columns must be equal to second matrix's number of rows.");
    }
    Matrix resultMatrix(rowSize, rightMatrix.colSize);

    for(int i = 0; i < resultMatrix.rowSize; i++) {
        for (int j = 0; j < resultMatrix.colSize; j++) {
            for (int k = 0; k < rightMatrix.rowSize; k++) {
                resultMatrix.matrix[i][j] += matrix[i][k] * rightMatrix.matrix[k][j];
            }
        }
    }
    *this = resultMatrix;
    return *this;
}

Matrix operator*(Matrix lhs, const Matrix& rhs){
    lhs *= rhs;
    return lhs;
}

void generateImportanceMatrix(Matrix& values) {
    int currentCount= 0;
    for (int j = 0; j < values.getRowSize(); j ++) {
        for (int i = 0; i < values.getColSize(); i++) {
            currentCount += values.getValue(i,j);
        }
        if (currentCount > 0) {
            for (int i = 0; i < values.getRowSize(); i++) {
                values.setValue(i, j, (values.getValue(i, j) / currentCount));
            }
        } else {
            for (int i = 0; i < values.getRowSize(); i++) {
                values.setValue(i, j, (1.0 / values.getColSize()));
            }
        }
        currentCount = 0;
    }
}


void generateQMatrix(Matrix& values) {
    int currentCount = 0;
    for (int j = 0; j < values.getColSize(); j ++) {
        for (int i = 0; i < values.getRowSize(); i++) {
            currentCount += values.getValue(i,j);
            for (int i = 0; i < values.getRowSize(); i++) {
                values.setValue(i, j, (1.0 / values.getColSize()));
            }
        }
        currentCount = 0;
    }
}


void generateTransitionMatrix(Matrix& Q, Matrix& S) {
    for (int j = 0; j < Q.getColSize(); j++) {
        for (int i = 0; i < Q.getRowSize(); i++) {
            Q.setValue(i, j, (Q.getValue(i, j) * (1-WALK_PROB)));
        }
    }
    for (int j = 0; j < S.getColSize(); j++) {
        for (int i = 0; i < S.getRowSize(); i++) {
            S.setValue(i, j, (S.getValue(i, j) * WALK_PROB));
            S.setValue(i, j, (S.getValue(i, j) + Q.getValue(i, j)));
        }
    }
}

void markovProcess(Matrix& transitionMatrix) {
    Matrix rankMatrix(transitionMatrix.getRowSize(), 1);
    for (int row = 0; row < rankMatrix.getRowSize(); row++) {
        for (int col = 0; col < rankMatrix.getColSize(); col++) {
            rankMatrix.setValue(row, 0, 1);
        }
    }

    Matrix initialMatrix(rankMatrix);
    rankMatrix = transitionMatrix * rankMatrix;
    while (initialMatrix != rankMatrix) {
        initialMatrix = rankMatrix;
        rankMatrix = transitionMatrix * rankMatrix;
    }
    double sum = 0.0;
    for (int row = 0; row < rankMatrix.getRowSize(); row++) {
        sum = sum + rankMatrix.getValue(row, 0);
    }

    for (int row = 0; row < rankMatrix.getRowSize(); row++) {
        rankMatrix.setValue(row, 0, rankMatrix.getValue(row, 0) * PERCENTAGE / sum);
    }
    getOutput(rankMatrix);
}

void getOutput(Matrix &rankMatrix) {
    for (int i = 0; i < rankMatrix.getRowSize(); i++) {
        string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        cout << "Page " << alphabet[i] << ": " << fixed << setprecision(PRECISION) << rankMatrix.getValue(i, 0) << "%" << endl;
    }
}

vector<double> fileReader() {
    fstream myFile;
    myFile.open("../connectivity.txt");
    double y;
    int counter = 0;
    vector<double> values;

    if (!myFile.is_open()) {
        cerr << "File cannot be opened!" << endl;
        exit(0);
    }

    while (!myFile.eof()) {
        myFile >> y;
        values.push_back(y);
        counter++;
    }
    myFile.close();
    return values;
}
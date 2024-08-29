//
// Created by Misam Ibrahimi on 2021-10-05.
//

#ifndef LAB1TEMPLATE_MATRIX_HPP
#define LAB1TEMPLATE_MATRIX_HPP

#include <vector>
#include <ostream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

class Matrix {

protected:
    int colSize;
    int rowSize;
    vector<vector<double>> matrix;

public:

    /**
     * Default constructor for a Matrix.
     */
    Matrix();

    /**
     * Constructs an n x n matrix.
     * @param n side size
     */
    explicit Matrix(const int n);

    /**
     * Constructs an r x c matrix.
     * @param r rows
     * @param c columns
     */
    Matrix(int r, int c);

    /**
     * Constructs a matrix containing a vector of doubles
     */
    explicit Matrix(vector<double>); // works

    /**
     * Copies contents of a matrix to otherMatrix.
     * @param otherMatrix toBeCopied matrix
     */
    Matrix(const Matrix &otherMatrix); // works

    /**
     * Sets the matrix value at matrix[r][c]
     * @param r row
     * @param c column
     * @param val value
     */
    void setValue(int r, int c, double val); // works

     /**
      * Gets the matrix value at matrix[r][c]
      * @param r row
      * @param c column
      * @return double
      */
    double getValue(int r, int c) const;

    /**
     * Clears the Matrix.
     */
    void clear(); // works

    /**
     * Destructor for a Matrix.
     */
    ~Matrix(); // works

    /**
     * Overloaded insertion operator.
     * @param os stream
     * @param obj matrix
     * @return ostream obj
     */
    friend ostream& operator<<(ostream &os, const Matrix& obj);

    /**
     * Overloaded == operator.
     * @param lhs this matrix
     * @param rhs other matrix
     * @return boolean
     */
    friend bool operator==(const Matrix& lhs, const Matrix& rhs);

    /**
     * Overloaded != operator.
     * @param lhs this matrix
     * @param rhs  other matrix
     * @return boolean
     */
    friend bool operator!=(const Matrix& lhs, const Matrix& rhs);

    /**
     * Increments all values by 1. Prefix.
     * @return ++value
     */
    Matrix& operator++();

    /**
     * Increments all values by 1. Postfix.
     * @return value++
     */
    Matrix operator++(int);

    /**
     * Decrement all values by 1. Prefix.
     * @return --value
     */
    Matrix& operator--();

    /**
     * Decrements all values by 1. Postfix.
     * @return value--
     */
    Matrix operator--(int);

    /**
     * Overloaded = operator using copy and swap method.
     * @param first this
     * @param second other matrix
     */
    friend void swapMatrix(Matrix& first, Matrix& second);

    /**
     * Overloaded = operator.
     * @param other matrix
     * @return matrix
     */
    Matrix& operator=( Matrix other);

    /**
     * Overloaded += operator.
     * @param rhs matrix
     * @return matrix
     */
    Matrix& operator+=(const Matrix& rhs);

    /**
     * Overloaded + operator.
     * @param lhs this matrix
     * @param rhs other matrix
     * @return matrix
     */
    friend Matrix operator+(Matrix lhs, const Matrix& rhs);

    /**
     * Overloaded -= operator.
     * @param rhs other matrix
     * @return matrix
     */
    Matrix& operator-=(const Matrix& rhs);

    /**
     * Overloaded - operator.
     * @param lhs this matrix
     * @param rhs other matrix
     * @return matrix
     */
    friend Matrix operator-(Matrix lhs, const Matrix& rhs);

    /**
     * Overloaded *= operator
     * @param rhs
     * @return matrix
     */
    Matrix& operator*=(const Matrix& rhs);

    /**
     * Overloaded * operator.
     * @param lhs this
     * @param rhs other
     * @return matrix
     */
    friend Matrix operator*(Matrix lhs, const Matrix& rhs);


    /**
     * Gets the column size.
     * @return int
     */
    int getColSize() const;

    /**
     * Gets the row size.
     * @return int
     */
    int getRowSize() const;
};

/**
 * Generates the stochastic matrix.
 * @param values connectivity matrix
 */
void generateImportanceMatrix(Matrix& values);

/**
 * Generates the Q matrix.
 * @param values stochastic matrix.
 */
void generateQMatrix(Matrix& values);

/**
 * Generates the transition matrix.
 * @param Q Q matrix
 * @param S stchastic matrix
 */
void generateTransitionMatrix(Matrix& Q, Matrix& S);

/**
 * Does the markov process on the transition matrix.
 * @param TransitionMatrix transition matrix
 */
void markovProcess(Matrix& TransitionMatrix);

/**
 * Gets the report.
 * @param rankMatrix matrix with ranks
 */
void getOutput(Matrix &rankMatrix);

/**
 * Reads the connectivity.txt file.
 * @return double
 */
vector<double> fileReader();

#endif //LAB1TEMPLATE_MATRIX_HPP

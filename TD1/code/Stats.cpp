// Include the library that defines input/output streams and operations on them
#include <iostream>
// Include the library with basic mathematical functions, e.g. sqrt
#include <cmath>
// Need that for memset
#include <cstring>

/********************
 * Helper functions *
 ********************/

// Read the data matrix from the standard input
void readMatrix (double **matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new double[columns];
        for (int j = 0; j < columns; j++)
            std::cin >> matrix[i][j];
    }
}

// Print an array on the standard output
void printArray (double values[], int length)
{
    for (int i = 0; i < length - 1; i++)
        std::cout << values[i] << " ";

    std::cout << values[length - 1];
}

// Prepare an empty matrix
double **prepareMatrix(int rows, int columns) {
    double **matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[columns];
        for (int j = 0; j < columns; j++)
            matrix[i][j] = 0;
    }

    return matrix;
}

/*****************************************************
 * Exercise 3: Mean, variance and standard deviation *
 *****************************************************/

/**
 * This function computes the mean of the given array of values
 *
 * @param values the array with the values
 * @param length the length of the array 'values'
 * @return the mean of the values in the array
 */
double computeMean (double values[], int length)
{
    double sum = 0;

    for (int i = 0; i < length; i++) {
        sum += values[i];
    }

    return sum / length;
}

/**
 * This function computes the variance of the given array of values
 *
 * @param values the array with the values
 * @param length the length of the array 'values'
 * @return the variance of the values in the array
 */
double computeVariance (double values[], int length)
{
    double mean = computeMean(values, length);
    double variance = 0;

    for (int i = 0; i < length; i++) {
        variance += values[i] * values[i];
    }

    variance /= length;
    variance -= mean * mean;

    return variance;
}


/**
 * This function computes the unbiased sample variance of the given
 * array of values
 *
 * @param values the array with the values
 * @param length the length of the array values
 * @return the variance of the values in the array
 */
double computeSampleVariance (double values[], int length)
{
    // TODO
    // Note that you can call the functions defined above
    double variance = computeVariance(values, length);

    double r = length;
    if(length != 1) r /= (length - 1);
    variance *= r;
    // Do not forget to replace this return by a correct one!
    return variance;
}

/**
 * This function computes the standard deviation of the given
 * array of values
 *
 * @param values the array with the values
 * @param length the length of the array values
 * @return the variance of the values in the array
 */
double computeStandardDeviation (double values[], int length)
{
    double standardDeviation = sqrt(computeVariance(values, length));

    return standardDeviation;
}

/**
 * This function computes the unbiased sample standard deviation
 * of the given array of values
 *
 * @param values the array with the values
 * @param length the length of the array values
 * @return the variance of the values in the array
 */
double computeSampleStandardDeviation (double values[], int length)
{
    double sampleStandardDeviation = sqrt(computeSampleVariance(values, length));

    return sampleStandardDeviation;
}

/************* Additional tests **************/

void runStudentTests3 (double values[], int length)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests3() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 3

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

/*************************************
 * Exercise 4: Working with matrices *
 *************************************/

/**
 * This function prints a rectangular matrix on the standard output,
 * placing each row on a separate line.
 *
 * @param matrix the matrix to print
 * @param rows the number of rows in the matrix
 * @param columns the number of columns
 */
void printMatrix (double **matrix, int rows, int columns)
{
    // TODO
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            std::cout<<matrix[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
}

/**
 * This function extracts one row from a data matrix
 *
 * @param matrix the matrix with the data
 * @param columns the number of columns in the matrix
 * @param index the index of the row to extract
 * @param row the array where the extracted values are to be placed
 */
void getRow (double **matrix, int columns, int index, double row[])
{
    // Clean out completely the row
    std::memset(row, 0, columns*sizeof(double));

    for (int i = 0; i < columns; i++) {
        row[i] = matrix[index][i];
    }

}

/**
 * This function extracts one column from a data matrix
 *
 * @param matrix the matrix with the data
 * @param rows the number of rows in the matrix
 * @param index the index of the column to extract
 * @param column the array where the extracted values are to be placed
 */
void getColumn (double **matrix, int rows, int index, double column[])
{
    // Clean out completely the column
    std::memset(column, 0, rows*sizeof(double));

    for (int i = 0; i < rows; i++) {
        column[i] = matrix[i][index];
    }

}

/************* Additional tests **************/

void runStudentTests4 (double **matrix, int rows, int columns)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests4() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 4

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

/****************************************************
 * Exercise 5: Covariance and correlation of arrays *
 ****************************************************/

/**
 * This function computes the covariance of two vectors of data of the same length
 * @param values1 the first vector
 * @param values2 the second vector
 * @param length the length of the two vectors
 * @return the covariance of the two vectors
 */
double computeCovariance(double values1[], double values2[], int length)
{
    if(length == 0)
        return 0;

    double covariance = 0;

    for (int i = 0; i < length; i++) {
        covariance += values1[i] * values2[i];
    }
    covariance/=(length);

    covariance -= computeMean(values1, length) * computeMean(values2, length);

    return covariance;
}

/**
 * This function computes the correlation of two vectors of data of the same length
 *
 * @param values1 the first vector
 * @param values2 the second vector
 * @param length the length of the two vectors
 * @return the correlation of the two vectors
 */
double computeCorrelation(double values1[], double values2[], int length)
{
    double correlation = computeCovariance(values1, values2, length) / ( computeStandardDeviation(values1, length) * computeStandardDeviation(values2, length)) ;
    return correlation;
}

/************* Additional tests **************/

void runStudentTests5 (double values1[], double values2[], int length)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests5() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 5

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

/***************************************************
 * Exercise 6: Covariance and correlation matrices *
 ***************************************************/

/**
 * This function computes the covariance matrix of the matrix provided as argument
 *
 * @param dataMatrix the input matrix
 * @param rows the number of rows in the matrix
 * @param columns the number of columns in the matrix
 * @return the covariance matrix
 */
double **computeCovarianceMatrix (double **dataMatrix, int rows, int columns)
{
    // Initialise a square matrix
    double **matrix = prepareMatrix(columns, columns);

    for (int i = 0; i < columns; i++) {
            double *X = new double[rows],*Y = new double[rows];
            getColumn (dataMatrix, rows, i, X);
            for (int j = 0; j < columns; j++) {
                getColumn (dataMatrix, rows, j, Y);
                matrix[i][j] = computeCovariance(X, Y, rows);
            }
    }
    return matrix;

}

/**
 * This function computes the correlation matrix of the matrix provided as argument
 *
 * @param dataMatrix the input matrix
 * @param rows the number of rows in the matrix
 * @param columns the number of columns in the matrix
 * @return the correlation matrix
 */
double **computeCorrelationMatrix (double **dataMatrix, int rows, int columns)
{
    // Initialise a square matrix
    double **matrix = prepareMatrix(columns, columns);
    for (int i = 0; i < columns; i++) {
        double *X = new double[rows],*Y = new double[rows];
        getColumn (dataMatrix, rows,i , X);
        for (int j = 0; j < columns; j++) {
            getColumn (dataMatrix, rows, j, Y);
            matrix[i][j] = computeCorrelation(X, Y, rows);
        }
    }
    return matrix;
}

/************* Additional tests **************/

void runStudentTests6 (double **matrix, int rows, int columns)
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests6() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 6

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}


/**
 * @file                Matrix.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               Offers useful transformations and functions on a 2D
 *                      array of floats that represent an mxn matrix. Can be
 *                      instantiated using m, n, float** or more simply by 
 *                      providing Vects in the form Matrix m({&v1,...,&vN});
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "VectND.h"
#include "VWrap.h"

namespace TPP_VMath
{
    #pragma region Matrix_Declaration
    //========================================================================//
    // START OF MATRIX CLASS DECLARATION
    //========================================================================//

    class Matrix
    {
    public:

        /**
         * @brief       Default constructor for a matrix is the 0 matrix
         * @param       void
         */
        Matrix();

        /**
         * @brief       Constructor for a matrix that inits with a single Vect
         *              reference to a child of Vect
         * @param       vect a child of the Vect class
         */
        Matrix(Vect&);

        /**
         * @brief       Constructor for a matrix that inits with an array of
         *              polymorphic Vect children. If all the provided Vects
         *              don't have the same dimensions, a TPP_VMath_Exception
         *              is thrown.
         * @param       vects an array of Vect* pointing to Vect children
         */
        Matrix(std::vector<Vect*>);

        /**
         * @brief       Constructor takes a VWrap object and unpacks the Vects
         *              out of it in order to init a matrix.  
         * @param wrapped 
         */
        Matrix(VWrap&);

        /**
         * @brief       Overloaded constructer accepts a float** in the form
         *              f[rows][cols]
         * @param       fMatrix a 2D array of floats
         * @param       m the number of rows in the matrix
         * @param       n the number of columns in the matrix
         */ 
        Matrix(int, int, float**);

        /**
         * @brief       Matrix copy constructor
         * @param       original a const Matrix reference
         */ 
        Matrix(const Matrix&);

        /**
         * @brief       Destructor deletes the entries of the matrix
         */ 
        ~Matrix();

        /**
         * @brief       Adds a Vect with matching dimensions as the final column
         *              to the matrix. If the dimensions don't match a
         *              TPP_VMath_Exception is thrown.
         * @param       vect a polymorphic pointer to a child vect
         */
        void            append_vect_to_matrix(Vect*);

        /**
         * @brief       Reduce this matrix to reduced echelon form
         * @param       void
         */ 
        void            reduce_matrix();

        /**
         * @brief       Produce a new reduced echelon matrix from this matrix
         * @param       void
         * @return      this matrix in RE form
         */ 
        Matrix          get_reduced();

        /**
         * @brief       Ax = b. Faulty dimensions throws TPP_VMath_Exception
         * @param       vect a polymorphic pointer to a child of Vect
         * @return      a pointer to b where b = Ax
         */ 
        Vect*           get_matrix_vector_product(Vect&);

        /**
         * @brief       Ax = b. Faulty dimensions throws TPP_VMath_Exception
         * @param       matrix a Matrix reference
         * @param       vect a polymorphic pointer to a child of Vect
         * @return      a pointer to b where b = Ax
         */ 
        static Vect*    get_matrix_vector_product(Matrix&, Vect&);

        /**
         * @brief       AB = C. Faulty dimensions throws TPP_VMath_Exception
         * @param       B a matrix whose m dimension must match this n dimension
         * @return      a matrix that is the product of AB
         */
        Matrix          mat_mul(Matrix&); 

        /**
         * @brief       AB = C. Faulty dimensions throws TPP_VMath_Exception
         * @param       A a matrix
         * @param       B a matrix
         * @return      a matrix that is the product of AB
         */
        static Matrix   mat_mul(Matrix&, Matrix&); 

        /**
         * @brief       Identity matrix is an nxn diagonal matrix in RE form
         * @param       n the number of rows and columns in the resulting matrix
         * @return      A newly created nxn identity matrix
         */
        static Matrix   get_identity_matrix_of_size(int);

        /**
         * @brief       Transpose the rows and columns of a copy of this matrix
         * @param       original the matrix to be copied
         * @return      original copy transposed
         */
        static Matrix   transpose_matrix(Matrix&);

        /**
         * @brief       Access the nth column of the matrix as a VectND.
         * @param       n the column to be returned
         * @return      a Vect of this->n dimensions
         */
        VectND          operator[](const int);

        /**
         * @brief       get the ith, jth value of the entries in the matrix
         * @param       i the ith row
         * @param       j the jth col
         * @return      entries[i][j] 
         */
        float           operator()(const int, const int);

        /**
         * @brief       Console representation of coefficient matrix for testing
         * @param       void
         */ 
        void            print();

    private:

        ///             The number of rows in the matrix
        unsigned int    m;

        ///             The number of columns in the matrix
        unsigned int    n;

        ///             2D array using pointer notation representing a matrix
        float**         entries;

        /**
         * @brief       Swaps the pointers for A[rowA] and A[rowB]
         * @param       rowA a value representing a row in A
         * @param       rowB a value representing a row in A
        */
        void            row_interchange(int rowA, int rowB);

        /**
         * @brief       Generate a zero vector in Rm
         * @param       m the number of entries in the vector
         * @return      a new instance of a Vect child
         */
        Vect*           get_b(const int&);
    };

    //========================================================================//
    // END OF MATRIX CLASS DECLARATION
    //========================================================================//
    #pragma endregion Matrix_Declaration
}

#endif
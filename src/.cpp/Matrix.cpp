/**
 * @file                Matrix.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/Matrix.h"
#include "../.h/VWrap.h"
#include "../.h/TPP_VMath_Exception.h"
#include "../.h/Vect.h"
#include "../.h/Vect1D.h"
#include "../.h/Vect2D.h"
#include "../.h/Vect3D.h"
#include "../.h/Vect4D.h"
#include "../.h/VectND.h"
#include <iostream>
#include <iomanip>

namespace TPP_VMath
{
    #pragma region Matrix_Implementation
    //========================================================================//
    // START OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

    Matrix::Matrix()
    {
        m = n = 0;
        entries = nullptr;
    }

    Matrix::Matrix(Vect& vect)
    {
        m = vect.get_dimension();
        n = 1;
        entries = new float*[m];

        for (int i = 0; i < m; ++i)
        {
            entries[i] = new float[1];
            entries[i][0] = vect.get_components()[i];
        }
    }

    Matrix::Matrix(std::vector<Vect*> vects)
    {
        n = vects.size();

        if (n == 0)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "Matrix(std::vector<Vect*>) ";
            error += "Argument must have size >= 1";
            throw TPP_VMath_Exception(error, 0x3a7);
        }

        m = vects[0]->get_dimension();
        entries = new float*[m];

        for (int i = 0; i < m; ++i)
        {
            entries[i] = new float[n];

            for (int j = 0; j < n; ++j)
            {
                Vect* temp = vects[j];
                if (temp->get_dimension() != m)
                {
                    std::string error = "TPP_VMath_Exception: ";
                    error += "Matrix(std::vector<Vect*>) ";
                    error += "All entries must have matching dimensions";
                    throw TPP_VMath_Exception(error, 0x3a7);
                }

                entries[i][j] = temp->get_components()[i];
            }
        }

    }

    Matrix::Matrix(VWrap& wrapped)
    {
        // wrapped contains a single wrapped Vect
        if (wrapped.get_vect())
        {
            Vect* vect = wrapped.get_vect();
            m = vect->get_dimension();
            n = 1;
            entries = new float*[m];

            for (int i = 0; i < m; ++i)
            {
                entries[i] = new float[1];
                entries[i][0] = vect->get_components()[i];
            }
        }

        // wrapped contains multiple wrapped Vects
        else 
        {
            Vect** vpp = wrapped.get_vpp();
            m = vpp[0]->get_dimension();
            n = wrapped.get_n();
            entries = new float*[m];

            for (int i = 0; i < m; ++i)
            {
                entries[i] = new float[n];

                for (int j = 0; j < n; ++j)
                {
                    Vect*  vect = vpp[j];

                    if (vect->get_dimension() != m)
                    {
                        std::string error = "TPP_VMath_Exception: ";
                        error += "Matrix(VWrap&) ";
                        error += "All entries must have matching dimensions";
                        throw TPP_VMath_Exception(error, 0x3a7);
                    }

                    float* comp = vect->get_components();
                    entries[i][j] = comp[i];
                }
            }
        }
    }

    Matrix::Matrix(int m, int n, float** fMatrix)
    {
        this->m = m;
        this->n = n;
        entries = fMatrix;
    }

    Matrix::Matrix(const Matrix& original)
    {
        this->m = original.m;
        this->n = original.n;
        entries = new float*[m];
        for (int i = 0; i < m; ++i)
        {
            entries[i] = new float[n];
            for (int j = 0; j < n; ++j)
            {
                entries[i][j] = original.entries[i][j];
            }
        }
    }

    Matrix::~Matrix()
    {
        for (int i = 0; i < m; ++i)
        {
            if (entries[i])
                delete [] entries[i];
        }

        if (entries)
            delete [] entries;
    }

    void    Matrix::append_vect_to_matrix(Vect* vect)
    {
        // matrix already contains entries
        if (entries)
        {
            if (vect->get_dimension() != m)
            {
                std::string error = "TPP_VMath_Exception: ";
                error += "append_vect_to_matrix(Vect*) requires all Vects to ";
                error += "have the same dimensions";
                throw TPP_VMath_Exception(
                    error,
                    0x3a7
                );
            }

            // increase the column count of the matrix
            ++n;
            float** p = new float*[m];

            // include all the current entries in the matrix
            for (int i = 0; i < m; ++i)
            {
                p[i] = new float[n];

                for (int j = 0; j < n - 1; ++j)
                {
                    p[i][j] = entries[i][j];
                }
            }

            // deallocate the old memory resource
            for (int i = 0; i < m; ++i)
            {
                delete [] entries[i];
            }
            delete [] entries;

            // add the components of vect to the remaining entries of p
            for (int i = 0; i < m; ++i)
            {
                p[i][n - 1] = vect->get_components()[i];
            }

            // entries = new entries
            entries = p;
        }

        // first entry in the matrix
        else
        {
            // 1 entry in the matrix
            ++n;

            // define the number of rows
            m = vect->get_dimension();

            // establish the matrix in memory
            entries = new float*[m];

            for (int i = 0; i < m; ++i)
            {
                entries[i] = new float[1];
                entries[i][0] = vect->get_components()[i];
            }
        }
    }

    void    Matrix::reduce_matrix()
    {
        int pivotCol = 0, pivotRow = 0;
        bool notInReducedEchelonForm = true;

        while (notInReducedEchelonForm)
        {
            // sf RE form achieved, the loop will break
            notInReducedEchelonForm = false;
            
            // search the row for the first non zero entry as our pivot
            while (pivotCol < n && entries[pivotCol][pivotRow] == 0)
            {
                ++pivotCol;
            }

            // if the row is all zeros, interchange with the first non-zero row
            if (pivotCol == n) {
                // TODO
                // TODO
                // TODO
                // TODO
                // TODO
                // TODO
            }

        }
    }

    Matrix  Matrix::get_reduced()
    {
        Matrix copy(*this);
        copy.reduce_matrix();
        return copy;
    }

    Vect*   Matrix::get_matrix_vector_product(Vect& vect)
    {
        unsigned int n = vect.get_dimension();

        // the dimensions of the vector must match the number of columns in A
        if (n != this->n) 
        {
            throw TPP_VMath_Exception(
                "x is not in the range of the transformation Ax",
                0x3a7
            );
        }

        // The matrix vector product
        Vect* b = get_b(this->m);

        float* x = vect.get_components();

        // Ax = b
        for (int i = 0; i < m; ++i) 
        {
            int sum = 0;

            for (int j = 0; j < n; ++j) 
            {
                sum += entries[i][j] * x[j];
            }

            b->get_components()[i] = sum;
        }

        return b;
    }

    Vect*   Matrix::get_matrix_vector_product(Matrix& matrix, Vect& vect)
    {
        unsigned int n = vect.get_dimension();

        // the dimensions of the vector must match the number of columns in A
        if (n != matrix.n) 
        {
            throw TPP_VMath_Exception(
                "x is not in the range of the transformation Ax",
                0x3a7
            );
        }

        // The matrix vector product
        Vect* b = matrix.get_b(matrix.m);

        float* x = vect.get_components();

        // Ax = b
        for (int i = 0; i < matrix.m; ++i) 
        {
            int sum = 0;

            for (int j = 0; j < n; ++j) 
            {
                sum += matrix.entries[i][j] * x[j];
            }

            b->get_components()[i] = sum;
        }

        return b;
    }

    Vect*   Matrix::get_b(const int& m)
    {
        switch (m)
        {
            case  1: return new Vect1D;
            case  2: return new Vect2D;
            case  3: return new Vect3D;
            case  4: return new Vect4D;
            default: return new VectND(m); // VECTND WHEN CREATED TODO TODO
        }
    }

    Matrix  Matrix::mat_mul(Matrix& B)
    {
        if (this->n != B.m)
        {
            throw TPP_VMath_Exception("The product of AB is undefined", 0xAB);
        }

        const int mRows = this->m;
        const int nCols = B.n;
        float** f = new float*[mRows];

        for (int i = 0; i < mRows; ++i)
        {
            f[i] = new float[nCols];
        }

        for (int i = 0; i < nCols; ++i)
        {
            for (int j = 0; j < mRows; ++j)
            {
                f[j][i] = 0;

                for (int k = 0; k < this->n; ++k)
                {
                    f[j][i] += this->entries[j][k] * B.entries[k][i];
                }
            }
        }
        
        return Matrix(mRows, nCols, f);
    }

    Matrix  Matrix::mat_mul(Matrix& A, Matrix& B)
    {
        if (A.n != B.m)
        {
            throw TPP_VMath_Exception("The product of AB is undefined", 0xAB);
        }

        const int mRows = A.m;
        const int nCols = B.n;
        float** f = new float*[mRows];

        for (int i = 0; i < mRows; ++i)
        {
            f[i] = new float[nCols];
        }

        for (int i = 0; i < nCols; ++i)
        {
            for (int j = 0; j < mRows; ++j)
            {
                f[j][i] = 0;

                for (int k = 0; k < A.n; ++k)
                {
                    f[j][i] += A.entries[j][k] * B.entries[k][i];
                }
            }
        }
        
        return Matrix(mRows, nCols, f);
    }

    Matrix  Matrix::get_identity_matrix_of_size(int n)
    {
        float** fMatrix = new float*[n];
        
        for (int i = 0; i < n; ++i)
        {
            float* row = new float[n];

            for (int j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    row[j] = 0;
                }
                else
                {
                    row[j] = 1;
                }
            }
            fMatrix[i] = row;
        }

        return Matrix(n, n, fMatrix);
    }

    void    Matrix::print()
    {
        for (int i = 0; i < n; ++i)
        {
            std::cout << std::right
                      << std::setw(10)
                      << "Col "
                      << i + 1;
        }
        std::cout << std::endl;
        for (int i = 0; i < (n*11) + 4; ++i)
        {
            std::cout << "=";
        }
        std::cout << std::endl;
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                std::cout << std::setprecision(2) 
                          << std::fixed
                          << std::right
                          << std::setw(10)
                          << entries[i][j]
                          << " ";
            }

            std::cout << "    |Row " << i + 1 << "\n" << std::endl;
        }
        std::cout << std::endl;
    }

    Matrix  Matrix::transpose_matrix(Matrix& original)
    {
        Matrix result;
        
        for (int i = 0; i < original.m; ++i)
        {
            Vect*  temp = original.get_b(original.n);
            std::vector<float> v;

            for (int j = 0; j < original.n; ++j)
            {
                v.push_back(original.entries[i][j]);
            }

            temp->set_components(v);

            result.append_vect_to_matrix(temp);

            delete temp;
        }

        return result;
    }

    VectND  Matrix::operator[](const int n)
    {
        float f[m];

        for (int i = 0; i < m; ++i)
        {
            f[i] = (entries[i][n]);
        }

        return VectND(m, f);
    }
 
    float   Matrix::operator()(const int i, const int j)
    {
        return entries[i][j];
    }

    //========================================================================//
    // END OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Matrix_Implementation
}
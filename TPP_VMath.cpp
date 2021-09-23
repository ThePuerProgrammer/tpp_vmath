/**
 * @file        TPP_VMath.h
 * @author      Jesse Rankins
 * @since       2021
 * @version     0.1
 */ 
#pragma region Includes
#include "TPP_VMath.h"
#include <cmath>
#include <iostream>
#include <iomanip>  // for printing the matrix to the console
#pragma endregion Includes

namespace TPP_VMath
{

    #pragma region Vect_Implementation
    //========================================================================//
    // START OF VECT CLASS IMPLEMENTATION
    //========================================================================//

    Vect::~Vect()
    {   }

    float   Vect::get_magnitude()
    {
        float sumOfSquares = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfSquares += components[i] * components[i];
        }

        return std::sqrt(sumOfSquares);
    }

    float   Vect::get_magnitude(Vect& vect)
    {
        float sumOfSquares = 0;
        float* components = vect.get_components(); 
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            sumOfSquares += components[i] * components[i];
        }

        return std::sqrt(sumOfSquares);
    }

    void    Vect::scale_by(int c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    void    Vect::scale_by(int c, Vect& vect)
    {
        float* components = vect.get_components();
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    void    Vect::scale_by(float c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    void    Vect::scale_by(float c, Vect& vect)
    {
        float* components = vect.get_components();
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            components[i] *= c;
        }
    }

    float   Vect::dot_product(const Vect& that)
    {
        if (this->dimension != that.get_dimension())
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "The dot product of a vector in R";
            error += std::to_string(this->dimension) + " ";
            error += "and a vector in R";
            error += std::to_string(that.dimension) + " ";
            error += "is undefined.";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        float sumOfProducts = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += this->components[i] * that.get_components()[i];
        }

        return sumOfProducts;
    }

    float   Vect::dot_product(const Vect& a, const Vect& b)
    {
        int dimension = a.get_dimension();

        if (dimension != b.get_dimension())
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "The dot product of a vector in R";
            error += std::to_string(a.dimension) + " ";
            error += "and a vector in R";
            error += std::to_string(b.dimension) + " ";
            error += "is undefined.";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        float sumOfProducts = 0;
        float* acomponents = a.get_components();
        float* bcomponents = b.get_components();

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += acomponents[i] * bcomponents[i];
        }

        return sumOfProducts;
    }

    void    Vect::normalize()
    {
        float magnitude = this->get_magnitude();

        if (magnitude == 0)
        {
            std::string error = "TPP_VMath_Exception: Divide by zero ";
            error += "while attempting to normalize Vect";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        float inverseMagnitude = 1 / magnitude;

        this->scale_by(inverseMagnitude);
    }

    void    Vect::normalize(Vect& vect)
    {
        float magnitude = vect.get_magnitude();

        if (magnitude == 0)
        {
            std::string error = "TPP_VMath_Exception: Divide by zero ";
            error += "while attempting to normalize Vect";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        float inverseMagnitude = 1 / magnitude;

        vect.scale_by(inverseMagnitude);
    }

    void    Vect::vect_addition(Vect& vect)
    {
        if (this->dimension != vect.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->components[i] += vect.components[i];
        }
    }

    void    Vect::vect_addition(Vect& a, Vect& b)
    {
        if (a.dimension != b.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        for (int i = 0; i < a.dimension; ++i)
        {
            a.components[i] += b.components[i];
        }
    }

    void    Vect::vect_subtraction(Vect& vect)
    {
        if (this->dimension != vect.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->components[i] -= vect.components[i];
        }
    }

    void    Vect::vect_subtraction(Vect& a, Vect& b)
    {
        if (a.dimension != b.dimension)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "vect_addition(Vect&) requires matching dimensions\n";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        for (int i = 0; i < a.dimension; ++i)
        {
            a.components[i] -= b.components[i];
        }
    }

    void    Vect::set_components(std::vector<float> v)
    {
        if (dimension != 0 && v.size() != dimension)
        {
            std::string error = "TPP_VMath_Exception: Provided entries of ";
            error += "set_components(std::vector<float>) do not match the ";
            error += "dimensions of the Vect from which is was called.";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        if (dimension == 0)
        {
            dimension = v.size();
            components = new float[v.size()];
        }

        for (int i = 0; i < dimension; ++i)
        {
            components[i] = v[i];
        }
    }

    float*  Vect::operator[](int i)
    {
        return &components[i];
    }

    void    Vect::operator*=(float c)
    {
        this->scale_by(c);
    }

    void    Vect::operator/=(float c)
    {
        if (c == 0)
        {
            std::string error = "TPP_VMath_Exception: ";
            error += "operator/=(float c) divide by zero.";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        this->scale_by(1/c);
    }

    void    Vect::operator+=(Vect& right)
    {
        this->vect_addition(right);
    }

    void    Vect::operator-=(Vect& right)
    {
        this->vect_subtraction(right);
    }

    Vect*   Vect::operator+(Vect& right)
    {
        Vect* result = get_vect_of_valid_dimensions(this);
        result->vect_addition(right);
        return result;
    }

    Vect*   Vect::operator-(Vect& right)
    {
        Vect* result = get_vect_of_valid_dimensions(this);
        result->vect_subtraction(right);
        return result;
    }

    Vect&   Vect::operator=(const Vect& right)
    {
        if (this->dimension != right.get_dimension())
        {
            std::string error = "TPP_VMath_Exception: Dimensions of ";
            error += "operator=(const Vect& right) do not match the ";
            error += "dimensions of the left side Vect.";
            throw TPP_VMath_Exception(
                error,
                0x4ec7
            );
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->components[i] = right.get_components()[i];
        }

        return *this;
    }

    Vect*   Vect::get_vect_of_valid_dimensions(Vect* vect)
    {
        float* c = vect->get_components();
        int    m = vect->get_dimension();

        switch (m)
        {
            case  1: return new Vect1D(c);
            case  2: return new Vect2D(c);
            case  3: return new Vect3D(c);
            case  4: return new Vect4D(c);
            default: return new VectND(m, c);
        }
    }

    float*  Vect::get_components() const
    {
        return components;
    }

    unsigned int Vect::get_dimension() const
    {
        return dimension;
    }

    void    Vect::print()
    {
        std::cout << "---\n";
        for (int i = 0; i < dimension; ++i)
        {
            std::cout << std::setprecision(10) << std::fixed
                      << "|" << components[i] << "|\n";
        }
        std::cout << "---\n";
    }

    //========================================================================//
    // END OF VECT CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect_Implementation

    #pragma region VectND_Implementation
    //========================================================================//
    // START OF VECTND CLASS IMPLEMENTATION
    //========================================================================//

    VectND::VectND()
    {   
        dimension = 0;
        components = nullptr;
    }

    VectND::VectND(int n)
    {
        dimension = n;
        components = new float[n];

        for (int i = 0; i < n; ++i)
        {
            components[i] = 0;
        }
    }

    VectND::VectND(int n, float* components)
    {
        dimension = n;
        this->components = new float[n];

        for (int i = 0; i < n; ++i)
        {
            this->components[i] = components[i];
        }
    }

    VectND::VectND(std::vector<float> v)
    {
        if (v.size() > 0)
        {
            dimension = v.size();
            components = new float[v.size()];

            for (int i = 0; i < dimension; ++i)
            {
                components[i] = v[i];
            }
        }
        else
        {
            dimension = 0;
            components = nullptr;
        }
    }

    VectND::~VectND()
    {
        if (components) 
        {
            delete [] components;
            components = nullptr;
        }
    }

    //========================================================================//
    // END OF VECTND CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion VectND_Implementation

    #pragma region Vect1D_Implementation
    //========================================================================//
    // START OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//

    Vect1D::Vect1D() : VectND()
    {   }

    Vect1D::Vect1D(std::array<float, 1> a) : VectND(1, a.begin())
    {   }

    Vect1D::Vect1D(const Vect1D& original) : VectND(1, original.components)
    {   }

    Vect1D::Vect1D(float* components) : VectND(1, components)
    {   }

    //========================================================================//
    // END OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect1D_Implementation

    #pragma region Vect2D_Implementation
    //========================================================================//~~
    // START OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//

    Vect2D::Vect2D() : Vect2D({0, 0})
    {   }

    Vect2D::Vect2D(std::array<float, 2> a) : VectND(2, a.begin())
    {   }

    Vect2D::Vect2D(const Vect2D& original) : VectND(2, original.components)
    {   }

    Vect2D::Vect2D(float* components) : VectND(2, components)
    {   }

    //========================================================================//
    // END OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect2D_Implementation

    #pragma region Vect3D_Implementation
    //========================================================================//
    // START OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

    Vect3D::Vect3D() : Vect3D({0, 0, 0}) 
    {   }

    Vect3D::Vect3D(std::array<float, 3> a) : VectND(3, a.begin())
    {   }

    Vect3D::Vect3D(const Vect3D& original) : VectND(3, original.components)
    {   }

    Vect3D::Vect3D(float* components) : VectND(3, components)
    {   }

    Vect3D  Vect3D::cross_product(Vect3D& that)
    {
        return Vect3D
        ({
           this->components[1] * that.components[2] - 
           this->components[2] * that.components[1],

           this->components[2] * that.components[0] - 
           this->components[0] * that.components[2],

           this->components[0] * that.components[1] - 
           this->components[1] * that.components[0]
        });
    }

    Vect3D  Vect3D::cross_product(Vect3D& a, Vect3D& b)
    {
        return Vect3D
        ({
           a.components[1] * b.components[2] - 
           a.components[2] * b.components[1],

           a.components[2] * b.components[0] - 
           a.components[0] * b.components[2],

           a.components[0] * b.components[1] - 
           a.components[1] * b.components[0]
        });
    }

    //========================================================================//
    // END OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect3D_Implementation

    #pragma region Vect4D_Implementation
    //========================================================================//
    // START OF VECT4D CLASS IMPLEMENTATION
    //========================================================================//

    Vect4D::Vect4D() : Vect4D({0, 0, 0, 0})
    {   }

    Vect4D::Vect4D(std::array<float, 4> a) : VectND(4, a.begin())
    {   }

    Vect4D::Vect4D(const Vect4D& original) : VectND(4, original.components)
    {   }

    Vect4D::Vect4D(float* components) : VectND(4, components)
    {   }

    //========================================================================//
    // END OF VECT4D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect4D_Implementation

    #pragma region VWrap_Implementation
    //========================================================================//
    // START OF VWRAP CLASS IMPLEMENTATION
    //========================================================================//

    VWrap::VWrap()
    {
        n = 0;
        this->vect = nullptr;
        this->vpp = nullptr;
    }

    VWrap::VWrap(Vect* vect)
    {
        n = 1;
        this->vect = vect;
        this->vpp = nullptr;
    }

    VWrap::VWrap(int n, Vect** vpp)
    {
        this->n = n;
        this->vect = nullptr;
        this->vpp = vpp;
    }

    VWrap::~VWrap()
    {
        if (this->vect)
        {
            delete vect;
        }
        else if (this->vpp)
        {
            for (int i = 0; i < n; ++i)
            {
                delete vpp[i];
            }
            delete [] vpp;
        }
    }

    void    VWrap::wrap(Vect* vect)
    {
        n = 1;
        if (!this->vect)
            this->vect = vect;
    }

    void    VWrap::wrap(int n, Vect** vpp)
    {
        this->n = n;
        if (!this->vpp)
            this->vpp = vpp;
    }

    void    VWrap::unwrap()
    {
        if (this->vect)
        {
            delete this->vect;
            this->vect = nullptr;
        }
        else if (this->vpp)
        {
            for (int i = 0; i < n; ++i)
            {
                delete this->vpp[i];
                this->vpp[i] = nullptr;
            }
            delete [] this->vpp;
            this->vpp = nullptr;
        }
    }

    Vect*   VWrap::get_vect()
    {
        return vect;
    }

    Vect**  VWrap::get_vpp()
    {
        return vpp;
    }

    int     VWrap::get_n()
    {
        return n;
    }

    //========================================================================//
    // END OF VWRAP CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion VWrap_Implementation

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

    #pragma region TPP_VMath_Exception_Implementation
    //========================================================================//
    // START OF TPP_VMATH_EXCEPTION CLASS IMPLEMENTATION
    //========================================================================//

    TPP_VMath_Exception::TPP_VMath_Exception(std::string msg, int num)
    : errCode(num), errMsg(msg), std::runtime_error(msg)
    {   }

    int     TPP_VMath_Exception::get_error_code() const
    {
        return errCode;
    }

    //========================================================================//
    // END OF TPP_VMATH_EXCEPTION CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion TPP_VMath_Exception_Implementation

}
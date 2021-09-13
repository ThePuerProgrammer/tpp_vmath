/**
 * @file        TPP_VMath.h
 * @author      Jesse Rankins
 * @since       2021
 * @version     0.1
 */ 
#include "TPP_VMath.h"
#include <cmath>
#include <iostream>
#include <iomanip> // for printing the matrix to the console

namespace TPP_VMath
{

    //========================================================================//
    // START OF VECT CLASS IMPLEMENTATION
    //========================================================================//

    Vect::~Vect()
    {   }

    float Vect::get_magnitude()
    {
        float sumOfSquares = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfSquares += coordinates[i] * coordinates[i];
        }

        return sqrtf(sumOfSquares);
    }

    void Vect::scale_by(int c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            coordinates[i] *= c;
        }
    }

    void Vect::scale_by(float c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            coordinates[i] *= c;
        }
    }

    float* Vect::operator[](int i)
    {
        return &coordinates[i];
    }

    float* Vect::get_coordinates() const
    {
        return coordinates;
    }

    unsigned int Vect::get_dimension() const
    {
        return dimension;
    }

    //========================================================================//
    // END OF VECT CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF VWRAPPER CLASS IMPLEMENTATION
    //========================================================================//

    VWrapper::VWrapper()
    {
        this->vect = nullptr;
    }

    VWrapper::VWrapper(Vect* vect)
    {
        this->vect = vect;
    }

    VWrapper::~VWrapper()
    {
        if (this->vect)
            delete vect;
    }

    void VWrapper::wrap(Vect* vect)
    {
        if (!this->vect)
            this->vect = vect;
    }

    void VWrapper::unwrap()
    {
        if (this->vect)
        {
            delete this->vect;
            this->vect = nullptr;
        }
    }

    Vect* VWrapper::get_vect()
    {
        return vect;
    }

    //========================================================================//
    // END OF VWRAPPER CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//

    Vect1D::Vect1D()
    {
        dimension = 1;
        coordinates = new float[1];
        coordinates[0] = 0;
    }

    Vect1D::Vect1D(float x)
    {
        dimension = 1;
        coordinates = new float[1];
        coordinates[0] = x;
    }

    Vect1D::Vect1D(const Vect1D& original)
    {
        dimension = 1;
        coordinates = new float[1];
        
        for (int i = 0; i < 1; ++i)
        {
            coordinates[i] = original.get_coordinates()[i];
        }
    }

    Vect1D::Vect1D(float* coordinates)
    {
        dimension = 1;
        this->coordinates = coordinates;
    }

    Vect1D::~Vect1D()
    {
        delete [] coordinates;
    }

    Vect1D& Vect1D::operator=(const Vect1D& right)
    {
        coordinates[0] = right.coordinates[0];

        return *this;
    }

    void Vect1D::set_coordinates(float x)
    {
        coordinates[0] = x;
    }

    float Vect1D::dot_product(const Vect1D& that)
    {
        return this->coordinates[0] * that.get_coordinates()[0];
    }

    //========================================================================//
    // END OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

    Vect3D::Vect3D() : Vect3D(0, 0, 0) 
    {   }

    Vect3D::Vect3D(const float x, const float y, const float z)
    {
        dimension = 3;
        coordinates = new float[3];
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    Vect3D::Vect3D(const Vect3D& original)
    {
        dimension = 3;
        coordinates = new float[3];
        
        for (int i = 0; i < 3; ++i)
        {
            coordinates[i] = original.get_coordinates()[i];
        }
    }

    Vect3D::Vect3D(const float* coordinates)
    {
        this->dimension = 3;
        this->coordinates = new float[3];

        for (int i = 0; i < 3; ++i)
        {
            this->coordinates[i] = coordinates[i];
        }
    }

    Vect3D::~Vect3D()
    {
        delete [] coordinates;
        coordinates = nullptr;
    }

    Vect3D& Vect3D::operator=(const Vect3D& right)
    {
        for (int i = 0; i < 3; ++i)
        {
            coordinates[i] = right.coordinates[i];
        }

        return *this;
    }

    void Vect3D::set_coordinates(float x, float y, float z)
    {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    float Vect3D::dot_product(const Vect3D& that)
    {
        float sumOfProducts = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += this->coordinates[i] * that.get_coordinates()[i];
        }

        return sumOfProducts;
    }

    //========================================================================//
    // END OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF VECT4D CLASS IMPLEMENTATION
    //========================================================================//

    Vect4D::Vect4D() : Vect4D(0, 0, 0, 0)
    {   }

    Vect4D::Vect4D(float x, float y, float z, float w)
    {
        dimension = 4;
        coordinates = new float[4];
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
        coordinates[3] = w;
    }

    Vect4D::Vect4D(const Vect4D& original)
    {
        dimension = 4;
        coordinates = new float[4];

        for (int i = 0; i < 4; ++i)
        {
            coordinates[i] = original.get_coordinates()[i];
        }
    }

    Vect4D::Vect4D(const float* coordinates)
    {
        dimension = 4;
        this->coordinates = new float[4];

        for (int i = 0; i < 4; ++i)
        {
            this->coordinates[i] = coordinates[i];
        }
    }

    Vect4D::~Vect4D()
    {
        delete [] coordinates;
    }

    Vect4D& Vect4D::operator=(const Vect4D& right)
    {
        for (int i = 0; i < 4; ++i)
        {
            coordinates[i] = right.coordinates[i];
        }

        return *this;
    }

    void Vect4D::set_coordinates(float x, float y, float z, float w)
    {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
        coordinates[3] = w;
    }

    float Vect4D::dot_product(const Vect4D& that)
    {
        float sumOfProducts = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += this->coordinates[i] * that.get_coordinates()[i];
        }

        return sumOfProducts;
    }

    //========================================================================//
    // END OF VECT4D CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF VSET CLASS IMPLEMENTATION
    //========================================================================//

    VSet::VSet()
    {
        m = n = 0;
        setOfVectors = nullptr;
    }

    VSet::VSet(const Vect3D& v3D)
    {
        n = 1;
        m = v3D.get_dimension();
        setOfVectors = new Vect*[1];
        setOfVectors[0] = new Vect3D(v3D.get_coordinates());
    }

    VSet::VSet(const int& n, Vect3D a[])
    {
        this->n = n;
        m = 3;
        setOfVectors = new Vect*[n];
        for (int i = 0; i < n; ++i)
        {
            setOfVectors[i] = new Vect3D(a[i].get_coordinates());
        }
    }

    VSet::VSet(const Vect4D& v4D)
    {
        n = 1;
        m = v4D.get_dimension();
        setOfVectors = new Vect*[1];
        setOfVectors[0] = new Vect4D(v4D.get_coordinates());
    }

    VSet::VSet(const VSet& original)
    {
        this->n = original.n;
        this->m = original.m;
        setOfVectors = new Vect*[n];

        for (int i = 0; i < n; ++i)
        {
            Vect* v = original.setOfVectors[i];
            setOfVectors[i] = add_vect_of_valid_dimensions(m, v);
        }
    }

    VSet::~VSet()
    {
        for (int i = 0; i < n; ++i)
        {
            delete setOfVectors[i];
        }
        delete [] setOfVectors;
    }

    unsigned int VSet::get_n() const
    {
        return n;
    }

    unsigned int VSet::get_m() const
    {
        return m;
    }

    Vect** VSet::get_set_of_vectors() const
    {
        return setOfVectors;
    }

    void VSet::add_vect_to_set(Vect* vect)
    {
        // set already contains entries
        if (setOfVectors)
        {
            if (vect->get_dimension() != m)
            {
                std::cerr << "The dimensions for vect must match the "
                          << "established VSet dimensions\n";
            }

            // increase the size of the set
            ++n;
            Vect** p = new Vect*[n];

            // add all in the current set to the new set
            for (int i = 0; i < n - 1; ++i)
            {
                p[i] = setOfVectors[i];
            }

            // deallocate the old memory resource
            delete [] setOfVectors;

            // add the new vect of m dimensions to the new set p
            p[n - 1] = add_vect_of_valid_dimensions(m, vect);

            // set = new set
            setOfVectors = p;
        }

        // first entry in the set
        else
        {
            // 1 entry in the set
            ++n;

            // set the dimension of the set
            m = vect->get_dimension();

            // establish the set in memory
            setOfVectors = new Vect*[1];

            // add the new vect of m dimension as the first entry in the set
            setOfVectors[0] = add_vect_of_valid_dimensions(m, vect);
        }
    }

    Vect* VSet::add_vect_of_valid_dimensions(int m, Vect* vect)
    {
        float* c = vect->get_coordinates();
        switch (m)
        {
            case  1: return new Vect1D(c);
            case  2: // return new Vect2D;
            case  3: return new Vect3D(c);
            case  4: return new Vect4D(c);
            default: return new Vect3D(c); // VECTND WHEN CREATED TODO TODO
        }
    }

    //========================================================================//
    // END OF VSET CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

    Matrix::Matrix(const VSet& set)
    {
        this->m = set.get_m();
        this->n = set.get_n();
        A = new float*[m];

        for (int i = 0; i < m; ++i)
        {
            A[i] = new float[n];

            for (int j = 0; j < n; ++j)
            {
                Vect* jthVect  = set.get_set_of_vectors()[j];
                float ithEntry = jthVect->get_coordinates()[i];
                A[i][j] = ithEntry;
            }
        }
    }

    Matrix::Matrix(int m, int n, float** fMatrix)
    {
        this->m = m;
        this->n = n;
        A = fMatrix;
    }

    Matrix::Matrix(const Matrix& original)
    {
        this->m = original.m;
        this->n = original.n;
        A = new float*[m];
        for (int i = 0; i < m; ++i)
        {
            A[i] = new float[n];
            for (int j = 0; j < n; ++j)
            {
                A[i][j] = original.A[i][j];
            }
        }
    }

    Matrix::~Matrix()
    {
        for (int i = 0; i < m; ++i)
        {
            delete [] A[i];
        }
        delete [] A;
    }

    void Matrix::print_matrix()
    {
        std::cout << "    ";
        for (int i = 0; i < n; ++i)
        {
            std::cout << std::right
                      << std::setw(10)
                      << "Col "
                      << i + 1;
        }
        std::cout << std::endl;
        for (int i = 0; i < m; ++i)
        {
            std::cout << "Row " << i + 1;
            for (int j = 0; j < n; ++j)
            {
                std::cout << std::setprecision(2) 
                          << std::fixed
                          << std::right
                          << std::setw(10)
                          << A[i][j]
                          << " ";
            }
            std::cout << std::endl;
        }
    }

    void Matrix::reduce_matrix()
    {
        int pivotCol = 0, pivotRow = 0;
        bool notInReducedEchelonForm = true;

        while (notInReducedEchelonForm)
        {
            // sf RE form achieved, the loop will break
            notInReducedEchelonForm = false;
            
            // search the row for the first non zero entry as our pivot
            while (pivotCol < n && A[pivotCol][pivotRow] == 0)
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

    Vect*   Matrix::get_matrix_vector_product(Vect* vect)
    {
        unsigned int n = vect->get_dimension();

        // the dimensions of the vector must match the number of columns in A
        if (n != this->n) 
        {
            std::cerr << "x is not in the range of the transformation Ax\n";
            return nullptr;
        }

        // The matrix vector product
        Vect* b = get_b(this->m);

        float* x = vect->get_coordinates();

        // Ax = b
        for (int i = 0; i < m; i++) 
        {
            int sum = 0;

            for (int j = 0; j < n; j++) 
            {
                sum += A[i][j] * x[j];
            }

            b->get_coordinates()[i] = sum;
        }

        return b;
    }

    Vect* Matrix::get_b(const int& m)
    {
        switch (m)
        {
            case  1: return new Vect1D;
            case  2: // return new Vect2D;
            case  3: return new Vect3D;
            case  4: return new Vect4D;
            default: return new Vect3D; // VECTND WHEN CREATED TODO TODO
        }
    }

    Matrix Matrix::get_identity_matrix_of_size(int n)
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

        Matrix result(n, n, fMatrix);
        return result;
    }

    //========================================================================//
    // END OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

}
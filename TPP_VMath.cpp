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
    { }

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
        if (m == 1)
        {
            // Vect1D
        }
        else if (m == 2)
        {
            // Vect2D
        }
        else if (m == 3) 
        {
            for (int i = 0; i < n; ++i)
            {
                setOfVectors[i] = 
                    new Vect3D(original.setOfVectors[i]->get_coordinates());
            }
        }
        else if (m == 4)
        {
            for (int i = 0; i < n; ++i)
            {
                setOfVectors[i] = 
                    new Vect4D(original.setOfVectors[i]->get_coordinates());
            }
        }
        else
        {
            // VectND
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
        try
        {
            // set already contains entries
            if (setOfVectors)
            {
                if (vect->get_dimension() != m)
                {
                    throw std::exception();
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

                // add the new vector
                // if set of Vect3D
                if (m == 3)
                {
                    p[n - 1] = new Vect3D(vect->get_coordinates());
                }
                else if (m == 4)
                {
                    p[n - 1] = new Vect4D(vect->get_coordinates());
                }

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

                // if set of Vect3D
                if (m == 3) 
                {
                    setOfVectors[0] = new Vect3D(vect->get_coordinates());
                }
                else if (m == 4)
                {
                    setOfVectors[0] = new Vect4D(vect->get_coordinates());
                }
            }
        }

        // the set dimensions didn't match
        catch(const std::exception& e)
        {
            std::cerr << "The dimensions for vect must match the current "
                      << "VSet dimensions\n";
            std::cerr << e.what() << '\n';
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
        Vect* b;

        if (n == 3)
        {
            b = new Vect3D();
        }

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
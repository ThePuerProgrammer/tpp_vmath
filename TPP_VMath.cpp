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
#include <cfloat>   // for -FLT_MAX
#pragma endregion Includes

namespace TPP_VMath
{

    #pragma region Vect_Implementation
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

    float Vect::get_magnitude(Vect& vect)
    {
        float sumOfSquares = 0;
        float* coordinates = vect.get_coordinates(); 
        int      dimension = vect.get_dimension();

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

    void Vect::scale_by(int c, Vect& vect)
    {
        float* coordinates = vect.get_coordinates();
        int      dimension = vect.get_dimension();

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

    void Vect::scale_by(float c, Vect& vect)
    {
        float* coordinates = vect.get_coordinates();
        int      dimension = vect.get_dimension();

        for (int i = 0; i < dimension; ++i)
        {
            coordinates[i] *= c;
        }
    }

    float Vect::dot_product(const Vect& that)
    {
        if (this->dimension != that.get_dimension())
        {
            std::cerr << "Dimensions do not match. Cannot evaulate dot product";
            return -FLT_MAX;
        }

        float sumOfProducts = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += this->coordinates[i] * that.get_coordinates()[i];
        }

        return sumOfProducts;
    }

    float Vect::dot_product(const Vect& a, const Vect& b)
    {
        int dimension = a.get_dimension();

        if (dimension != b.get_dimension())
        {
            std::cerr << "Dimensions do not match. Cannot evaulate dot product";
            return -FLT_MAX;
        }

        float sumOfProducts = 0;
        float* aCoordinates = a.get_coordinates();
        float* bCoordinates = b.get_coordinates();

        for (int i = 0; i < dimension; ++i)
        {
            sumOfProducts += aCoordinates[i] * bCoordinates[i];
        }

        return sumOfProducts;
    }

    void Vect::normalize_vect()
    {
        float magnitude = this->get_magnitude();

        if (magnitude == 0)
        {
            std::cerr << "Cannot normalize the zero vector. Divide by zero "
                      << "error\n";
            return;
        }

        float inverseMagnitude = 1 / magnitude;

        this->scale_by(inverseMagnitude);
    }

    void Vect::normalize_vect(Vect& vect)
    {
        float magnitude = vect.get_magnitude();

        if (magnitude == 0)
        {
            std::cerr << "Cannot normalize the zero vector. Divide by zero "
                      << "error\n";
            return;
        }

        float inverseMagnitude = 1 / magnitude;

        vect.scale_by(inverseMagnitude);
    }

    void Vect::set_coordinates(std::vector<float> v)
    {
        if (dimension != 0 && v.size() != dimension)
        {
            std::cerr << "Dimensions do not match. Cannot assign coordinates\n";
            return;
        }

        if (dimension == 0)
        {
            dimension = v.size();
            coordinates = new float[v.size()];
        }

        for (int i = 0; i < dimension; ++i)
        {
            coordinates[i] = v[i];
        }
    }

    float* Vect::operator[](int i)
    {
        return &coordinates[i];
    }


    Vect& Vect::operator=(const Vect& right)
    {
        if (this->dimension != right.get_dimension())
        {
            std::cerr << "Dimensions do not match. Cannot assign coordinates\n";
            return *this;
        }

        for (int i = 0; i < this->dimension; ++i)
        {
            this->coordinates[i] = right.get_coordinates()[i];
        }

        return *this;
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
    #pragma endregion Vect_Implementation

    #pragma region VectND_Implementation
    //========================================================================//
    // START OF VECTND CLASS IMPLEMENTATION
    //========================================================================//

    VectND::VectND()
    {   
        dimension = 0;
        coordinates = nullptr;
    }

    VectND::VectND(int n)
    {
        dimension = n;
        coordinates = new float[n];

        for (int i = 0; i < n; ++i)
        {
            coordinates[i] = 0;
        }
    }

    VectND::VectND(int n, float* coordinates)
    {
        dimension = n;
        this->coordinates = new float[n];

        for (int i = 0; i < n; ++i)
        {
            this->coordinates[i] = coordinates[i];
        }
    }

    VectND::VectND(std::vector<float> v)
    {
        if (v.size() > 0)
        {
            dimension = v.size();
            coordinates = new float[v.size()];

            for (int i = 0; i < dimension; ++i)
            {
                coordinates[i] = v[i];
            }
        }
        else
        {
            dimension = 0;
            coordinates = nullptr;
        }
    }

    VectND::~VectND()
    {
        if (coordinates) 
        {
            delete [] coordinates;
            coordinates = nullptr;
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

    Vect1D::Vect1D(const Vect1D& original) : VectND(1, original.coordinates)
    {   }

    Vect1D::Vect1D(float* coordinates) : VectND(1, coordinates)
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

    Vect2D::Vect2D(const Vect2D& original) : VectND(2, original.coordinates)
    {   }

    Vect2D::Vect2D(float* coordinates) : VectND(2, coordinates)
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

    Vect3D::Vect3D(const Vect3D& original) : VectND(3, original.coordinates)
    {   }

    Vect3D::Vect3D(float* coordinates) : VectND(3, coordinates)
    {   }

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

    Vect4D::Vect4D(const Vect4D& original) : VectND(4, original.coordinates)
    {   }

    Vect4D::Vect4D(float* coordinates) : VectND(4, coordinates)
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

    void VWrap::wrap(Vect* vect)
    {
        n = 1;
        if (!this->vect)
            this->vect = vect;
    }

    void VWrap::wrap(int n, Vect** vpp)
    {
        this->n = n;
        if (!this->vpp)
            this->vpp = vpp;
    }

    void VWrap::unwrap()
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

    Vect* VWrap::get_vect()
    {
        return vect;
    }

    Vect** VWrap::get_vpp()
    {
        return vpp;
    }

    int VWrap::get_n()
    {
        return n;
    }

    //========================================================================//
    // END OF VWRAP CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion VWrap_Implementation

    #pragma region VSet_Implementation
    //========================================================================//
    // START OF VSET CLASS IMPLEMENTATION
    //========================================================================//

    VSet::VSet()
    {
        m = n = 0;
        setOfVectors = nullptr;
    }

    VSet::VSet(Vect& vect)
    {
        n = 1;
        m = vect.get_dimension();
        setOfVectors = new Vect*[1];
        setOfVectors[0] = add_vect_of_valid_dimensions(m, &vect);
    }

    VSet::VSet(int n, Vect** a)
    {
        this->n = n;
        m = a[0]->get_dimension();
        setOfVectors = a;
    }

    VSet::VSet(VWrap& wrapper)
    {
        this->n = wrapper.get_n();
        Vect** temp = wrapper.get_vpp();
        m = temp[0]->get_dimension();
        setOfVectors = new Vect*[n];
        
        for (int i = 0; i < n; ++i)
        {
            setOfVectors[i] = add_vect_of_valid_dimensions(m, temp[i]);
        }
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
                return;
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
            case  2: return new Vect2D(c);
            case  3: return new Vect3D(c);
            case  4: return new Vect4D(c);
            default: return new VectND(m, c);
        }
    }

    //========================================================================//
    // END OF VSET CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion VSet_Implementation

    #pragma region Matrix_Implementation
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
                          << A[i][j]
                          << " ";
            }

            std::cout << "    ||Row " << i + 1 << std::endl;
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

    Vect*   Matrix::get_matrix_vector_product(Vect& vect)
    {
        unsigned int n = vect.get_dimension();

        // the dimensions of the vector must match the number of columns in A
        if (n != this->n) 
        {
            std::cerr << "x is not in the range of the transformation Ax\n";
            return nullptr;
        }

        // The matrix vector product
        Vect* b = get_b(this->m);

        float* x = vect.get_coordinates();

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

    Vect*   Matrix::get_matrix_vector_product(Matrix& matrix, Vect& vect)
    {
        unsigned int n = vect.get_dimension();

        // the dimensions of the vector must match the number of columns in A
        if (n != matrix.n) 
        {
            std::cerr << "x is not in the range of the transformation Ax\n";
            return nullptr;
        }

        // The matrix vector product
        Vect* b = matrix.get_b(matrix.m);

        float* x = vect.get_coordinates();

        // Ax = b
        for (int i = 0; i < matrix.m; i++) 
        {
            int sum = 0;

            for (int j = 0; j < n; j++) 
            {
                sum += matrix.A[i][j] * x[j];
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
            case  2: return new Vect2D;
            case  3: return new Vect3D;
            case  4: return new Vect4D;
            default: return new VectND(m); // VECTND WHEN CREATED TODO TODO
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
    #pragma endregion Matrix_Implementation
}
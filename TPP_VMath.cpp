// TPP_VMath.cpp
// written by Jesse Rankins 2021
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

    // Returns sqrtf(v₁^2+...+vᵢ^2)
    float Vect::get_magnitude()
    {
        float sumOfSquares = 0;

        for (int i = 0; i < dimension; ++i)
        {
            sumOfSquares += coordinates[i] * coordinates[i];
        }

        return sqrtf(sumOfSquares);
    }

    // Scales each entry in the vector by integer scalar value
    void Vect::scale_by(int c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            coordinates[i] *= c;
        }
    }

    // Scales each entry in the vector by float scalar value
    void Vect::scale_by(float c)
    {
        for (int i = 0; i < dimension; ++i)
        {
            coordinates[i] *= c;
        }
    }

    // Overloaded [] operator returns the ith coordinate in the vector
    float* Vect::operator[](int i)
    {
        return &coordinates[i];
    }

    // Returns the coordinate array as a pointer to coordinates[0]
    float* Vect::get_coordinates() const
    {
        return coordinates;
    }

    // Returns the dimension of the vector (the number of entries)
    int Vect::get_dimension() const
    {
        return dimension;
    }

    //========================================================================//
    // END OF VECT CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

    // Default constructor for a vector in R3 that inits to the zero vector
    Vect3D::Vect3D() : Vect3D(0, 0, 0)
    { }

    // Overloaded constructor for a vector in R3 that accepts x,y,z coordinates
    Vect3D::Vect3D(const float x, const float y, const float z)
    {
        dimension = 3;
        coordinates = new float[3];

        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    // Copy constructor
    Vect3D::Vect3D(const Vect3D& original)
    {
        dimension = 3;
        coordinates = new float[3];
        
        for (int i = 0; i < 3; ++i)
        {
            coordinates[i] = original.get_coordinates()[i];
        }
    }

    /*
    Constructor useful for polymorphic copies of Vect* v = new Vect3D();
    Error prone!! Call should ALWAYS be wrapped with guard in the form
    if (v->get_dimension() == 3)
    {
        Vect3D(v->get_components());
    }
    */
    Vect3D::Vect3D(const float* coordinates)
    {
        this->dimension = 3;
        this->coordinates = new float[3];
        for (int i = 0; i < 3; ++i)
        {
            this->coordinates[i] = coordinates[i];
        }
    }

    // Destructor deletes coordinate array
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

    // Manual assignment of x, y, z vector entries
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
    // START OF VSET CLASS IMPLEMENTATION
    //========================================================================//

    // Default constructor is a null set
    VSet::VSet()
    {
        m = n = 0;
        setOfVectors = nullptr;
    }

    // Overloaded constructor for one Vect3D
    VSet::VSet(const Vect3D& v3D)
    {
        n = 1;
        m = v3D.get_dimension();
        setOfVectors = new Vect*[1];
        setOfVectors[0] = new Vect3D(v3D.get_coordinates());
    }

    // Overloaded constructor for an array of Vect3D
    VSet::VSet(int n, Vect3D** a)
    {
        this->n = n;
        m = 3;
        setOfVectors = new Vect*[n];
        for (int i = 0; i < n; ++i)
        {
            setOfVectors[i] = new Vect3D(a[i]->get_coordinates());
        }
    }

    VSet::VSet(const VSet& original)
    {
        this->n = original.n;
        this->m = original.m;
        setOfVectors = new Vect*[n];
        for (int i = 0; i < n; ++i)
        {
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
                setOfVectors[i] = 
                    new Vect3D(original.setOfVectors[i]->get_coordinates());
            }
            else if (m == 4)
            {
                // Vect4D
            }
            else
            {
                // VectND
            }
        }
    }

    // Destructor deletes the set
    VSet::~VSet()
    {
        std::cout << n << std::endl;
        for (int i = 0; i < n; ++i)
        {
            delete setOfVectors[i];
            setOfVectors[i] = nullptr;
        }
        delete [] setOfVectors;
    }

    // Returns the number of vectors in the set
    int VSet::get_n()
    {
        return n;
    }

    // Returns the number of entries in each vector
    int VSet::get_m()
    {
        return m;
    }

    // Returns the set as an array of Vects
    Vect** VSet::get_set_of_vectors()
    {
        return setOfVectors;
    }

    /*
    Add a new vect of *matching dimension* to the set
    If the vector dimension doesn't == m, a runtime exception is thrown
    It is best to wrap the call to this funtion in a guard in the form
    if (Vect*->get_dimension() == Set.get_m())
    {
        add_vect_to_set(Vect*);
    }
    This does not apply if the set is currently null
    */
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
                    p[n - 1] = new Vect3D(vect->get_coordinates());

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
                    setOfVectors[0] = new Vect3D(vect->get_coordinates());
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

    // Constructor builds an mxn matrix A from a set of n Vects
    Matrix::Matrix(VSet& set)
    {
        this->m = set.get_m();
        this->n = set.get_n();
        A = new float*[n];
        for (int i = 0; i < n; ++i)
        {
            A[i] = new float[m];
            for (int j = 0; j < m; ++j)
            {
                A[i][j] = set.get_set_of_vectors()[i]->get_coordinates()[j]; 
            }
        }
    }

    // Matrix copy constructor
    Matrix::Matrix(const Matrix& original)
    {
        this->m = original.m;
        this->n = original.n;
        A = new float*[n];
        for (int i = 0; i < n; ++i)
        {
            A[i] = new float[m];
            for (int j = 0; j < m; ++j)
            {
                A[i][j] = original.A[i][j];
            }
        }
    }

    // Destructor deletes A
    Matrix::~Matrix()
    {
        for (int i = 0; i < n; ++i)
        {
            delete [] A[i];
            A[i] = nullptr;
        }
        delete [] A;
        A = nullptr;
    }

    // Console representation of coefficient matrix for testing
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
                          << A[j][i]
                          << " ";
            }
            std::cout << std::endl;
        }
    }



    // Reduce this matrix to reduced echelon form
    void Matrix::reduce_matrix()
    {
        bool notInReducedEchelonForm = true;

        while (notInReducedEchelonForm)
        {
            // If RE form achieved, the loop will break
            notInReducedEchelonForm = false;

            // STEP 1. TODO

        }
    }

    // Produce a new reduced echelon matrix from this matrix
    Matrix  Matrix::get_reduced()
    {
        Matrix copy(*this);
        copy.reduce_matrix();
        return copy;
    }

    //========================================================================//
    // END OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

}
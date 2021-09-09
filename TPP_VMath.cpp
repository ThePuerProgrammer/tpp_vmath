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
    float* Vect::get_coordinates()
    {
        return coordinates;
    }

    // Returns the dimension of the vector (the number of entries)
    int Vect::get_dimension()
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
    Vect3D::Vect3D(const Vect3D& oldVect3D)
    {
        dimension = 3;
        coordinates = new float[3];
        
        for (int i = 0; i < 3; ++i)
        {
            coordinates[i] = oldVect3D.coordinates[i];
        }
    }

    // Destructor deletes coordinate array
    Vect3D::~Vect3D()
    {
        delete [] coordinates;
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

    // Intended only to keep Vect as a pure virtual abstract class
    void Vect3D::virtualizer()
    { }
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
    VSet::VSet(Vect3D v3D)
    {
        n = 1;
        m = v3D.get_dimension();
        setOfVectors = new Vect*[1];
        setOfVectors[0] = &v3D;
    }

    // Overloaded constructor for an array of Vect3D
    VSet::VSet(int n, Vect3D** a)
    {
        this->n = n;
        m = 3;
        setOfVectors = new Vect*[n];
        for (int i = 0; i < n; ++i)
        {
            setOfVectors[i] = a[i];
        }
    }

    // Destructor deletes the set
    VSet::~VSet()
    {
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
                p[n - 1] = vect;

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
                setOfVectors[0] = vect;
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
    Matrix::Matrix(VSet set)
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

    // Destructor deletes A
    Matrix::~Matrix()
    {
        for (int i = 0; i < n; ++i)
        {
            delete [] A[i];
        }
        delete [] A;
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

    //========================================================================//
    // END OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

}
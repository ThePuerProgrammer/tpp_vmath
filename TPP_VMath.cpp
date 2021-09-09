// TPP_VMath.cpp
// written by Jesse Rankins 2021
#include "TPP_VMath.h"
#include <cmath>

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

    void Vect3D::virtualizer()
    { }
    //========================================================================//
    // END OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

    //========================================================================//
    // START OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

    // Constructor builds an mxn matrix A from a set of n Vects
    Matrix::Matrix(int m, int n, Vect **set)
    {
        this->m = m;
        this->n = n;
        A = new float*[n];
        for (int i = 0; i < n; ++i)
        {
            A[i] = new float[m];
            for (int j = 0; j < m; ++j)
            {
                A[i][j] = set[i]->get_coordinates()[j]; 
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

    //========================================================================//
    // END OF MATRIX CLASS IMPLEMENTATION
    //========================================================================//

}
// TPP_VMath.cpp
// written by Jesse Rankins 2021
#include "TPP_VMath.h"
#include <cmath>

namespace TPP_VMath
{
    //========================================================================//
    // VECT CLASS
    //========================================================================//

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
    // VECT3D CLASS
    //========================================================================//

    // Default constructor for a vector in R3 that inits to the zero vector
    Vect3D::Vect3D()
    {
        Vect3D(0,0,0);
    }

    // Overloaded constructor for a vector in R3 that accepts x,y,z coordinates
    Vect3D::Vect3D(float x, float y, float z)
    {
        coordinates = new float[3];
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    // Destructor deletes coordinate array
    Vect3D::~Vect3D()
    {
        delete [] coordinates;
    }

    // Manual assignment of x, y, z vector entries
    void Vect3D::set_coordinates(float x, float y, float z)
    {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    // Returns sqrtf(x^2 + y^2 + z^2)
    float Vect3D::get_magnitude()
    {
        float xSquared = coordinates[0] * coordinates[0];
        float ySquared = coordinates[1] * coordinates[1];
        float zSquared = coordinates[2] * coordinates[2];
        float sumOfSquares = xSquared + ySquared + zSquared;
        return sqrtf(sumOfSquares);
    }

    // Scales each entry in the vector by integer scalar value
    void Vect3D::scale_by(int c)
    {
        coordinates[0] *= c;
        coordinates[1] *= c;
        coordinates[2] *= c;
    }

    // Scales each entry in the vector by float scalar value
    void Vect3D::scale_by(float c)
    {
        coordinates[0] *= c;
        coordinates[1] *= c;
        coordinates[2] *= c;
    }

    //========================================================================//
    // END OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

}
// TPP_VMath.cpp
// written by Jesse Rankins 2021
#include "TPP_VMath.h"
#include <cmath>

namespace TPP_VMath
{

    // VECTOR CLASS
    //========================================================================//
    float Vector::operator[](int i)
    {
        return coordinates[i];
    }

    float* Vector::get_coordinates()
    {
        return coordinates;
    }
    //========================================================================//

    // VECTOR3D CLASS
    //========================================================================//
    Vector3D::Vector3D()
    {
        // default: initialize to the zero vector in R3;
        Vector3D(0,0,0);
    }

    Vector3D::Vector3D(float x, float y, float z)
    {
        coordinates = new float(3);
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    Vector3D::~Vector3D()
    {
        delete [] coordinates;
    }

    void Vector3D::set_coordinates(float x, float y, float z)
    {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    float Vector3D::get_magnitude()
    {
        float xSquared = coordinates[0] * coordinates[0];
        float ySquared = coordinates[1] * coordinates[1];
        float zSquared = coordinates[2] * coordinates[2];
        float sumOfSquares = xSquared + ySquared + zSquared;
        return sqrtf(sumOfSquares);
    }

    void Vector3D::scale_by(int c)
    {
        coordinates[0] *= c;
        coordinates[1] *= c;
        coordinates[2] *= c;
    }

    void Vector3D::scale_by(float c)
    {
        coordinates[0] *= c;
        coordinates[1] *= c;
        coordinates[2] *= c;
    }
    //========================================================================//

}
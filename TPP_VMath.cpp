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

    float* Vector::getCoordinates()
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

    void Vector3D::setCoordinates(float x, float y, float z)
    {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    }

    float Vector3D::getMagnitude()
    {
        float xSquared = coordinates[0] * coordinates[0];
        float ySquared = coordinates[1] * coordinates[1];
        float zSquared = coordinates[2] * coordinates[2];
        float sumOfSquares = xSquared + ySquared + zSquared;
        return sqrtf(sumOfSquares);
    }
    //========================================================================//

}
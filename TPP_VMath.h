// TPP_VMath.h
// written by Jesse Rankins 2021
#pragma once

namespace TPP_VMath
{
    //========================================================================//
    // VECT CLASS
    //========================================================================//

    class Vect
    {
    public:
        virtual float   get_magnitude() = 0;
        virtual void    scale_by(int) = 0;
        virtual void    scale_by(float) = 0;

        // Overloaded [] operator returns the ith coordinate in the vector
        float*          operator[](int);

        // Returns the coordinate array as a pointer to coordinates[0]
        float*          get_coordinates();
        virtual void    set_coordinates(float, float, float) = 0;

    protected:
        float* coordinates;
    };

    //========================================================================//
    // END OF VECT CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // VECT CLASS
    //========================================================================//

    class Vect3D : protected Vect
    {
    public:
        // Default constructor for a vector in R3 that inits to the zero vector
        Vect3D();

        /*
        Overloaded constructor for a vector in R3 that accepts x,y,z coordinates
        */
        Vect3D(float, float, float);

        // Destructor deletes coordinate array
        ~Vect3D();

        // Returns sqrtf(x^2 + y^2 + z^2)
        float   get_magnitude();

        // Scales each entry in the vector by integer scalar value
        void    scale_by(int);

        // Scales each entry in the vector by float scalar value
        void    scale_by(float);

        // Manual assignment of x, y, z vector entries
        void    set_coordinates(float, float, float);
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//
};
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
        // Returns sqrtf(v₁^2+...+vᵢ^2)
        float           get_magnitude();

        // Scales each entry in the vector by integer scalar value
        void            scale_by(int);

        // Scales each entry in the vector by float scalar value
        void            scale_by(float);

        // Returns pointer to the ith coordinate in the vector
        float*          operator[](int);

        // Returns the coordinate array as a pointer to coordinates[0]
        float*          get_coordinates();

        // Returns the dimension of the vector
        int             get_dimension();

        // Temporary virtualizer to force abstraction
        virtual void    virtualizer() = 0;

    protected:
        int             dimension;
        float*          coordinates;
    };

    //========================================================================//
    // END OF VECT CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // VECT CLASS
    //========================================================================//

    class Vect3D : public Vect
    {
    public:
        // Default constructor for a vector in R3 that inits to the zero vector
        Vect3D();

        /*
        Overloaded constructor for a vector in R3 that accepts x,y,z coordinates
        */
        Vect3D(const float, const float, const float);

        // Copy constructor
        Vect3D(const Vect3D&);

        // Destructor deletes coordinate array
        ~Vect3D();

        // Copies coordinates from the right side Vect3D to the left side Vect3D
        Vect3D& operator=(const Vect3D&);

        // Manual assignment of x, y, z vector entries
        void    set_coordinates(float, float, float);

        // Unimplemented virtualizer forces abstraction
        void    virtualizer();
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//
};
// TPP_VMath.h
// written by Jesse Rankins 2021
#pragma once

namespace TPP_VMath
{
    //========================================================================//
    // START OF VECT CLASS DECLARATION
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
        // The number of entries in coordinates
        int             dimension;

        // The array representing the vector
        float*          coordinates;
    };

    //========================================================================//
    // END OF VECT CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT3D CLASS DECLARATION
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

        // Unimplemented. Do NOT call
        void    virtualizer();
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF MATRIX CLASS DECLARATION
    //========================================================================//

    class Matrix
    {
    public:
    
        // Overloaded constructor accepts a set of Vects and generates a matrix
        Matrix(int, int, Vect**);

        // Destructor deletes the matrix
        ~Matrix();

    private:
        int     m;
        int     n;
        float   **A;
    };

    //========================================================================//
    // END OF MATRIX CLASS DECLARATION
    //========================================================================//
};
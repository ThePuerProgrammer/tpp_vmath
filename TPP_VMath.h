// TPP_VMath.h
// written by Jesse Rankins 2021
#pragma once

namespace TPP_VMath
{
    class Vector
    {
    public:
        // PUBLIC FUNCTIONS
        //--------------------------------------------------------------------//
        virtual float   get_magnitude() = 0;
        virtual void    scale_by(int) = 0;
        virtual void    scale_by(float) = 0;
        //--------------------------------------------------------------------//

        // OP OVERLOADS
        //--------------------------------------------------------------------//
        float           operator[](int i);
        //--------------------------------------------------------------------//

        // COORDINATES GETTERS/SETTERS
        //--------------------------------------------------------------------//
        float*          get_coordinates();
        virtual void    set_coordinates(float, float, float) = 0;
        //--------------------------------------------------------------------//

    protected:
        // xyz coordinates
        float* coordinates;
    };

    class Vector3D : protected Vector
    {
    public:
        // CONSTRUCTORS/DESTRUCTOR
        //--------------------------------------------------------------------//
        Vector3D();
        Vector3D(float, float, float);

        ~Vector3D();
        //--------------------------------------------------------------------//

        // PUBLIC FUNCTIONS
        //--------------------------------------------------------------------//
        float   get_magnitude();
        void    scale_by(int);
        void    scale_by(float);
        //--------------------------------------------------------------------//

        // COORDINATES GETTERS/SETTERS
        //--------------------------------------------------------------------//
        void  set_coordinates(float, float, float);
        //--------------------------------------------------------------------//
    };
};
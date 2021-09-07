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
        virtual float   getMagnitude() = 0;
        //--------------------------------------------------------------------//

        // OP OVERLOADS
        //--------------------------------------------------------------------//
        float           operator[](int i);
        //--------------------------------------------------------------------//

        // COORDINATES GETTERS/SETTERS
        //--------------------------------------------------------------------//
        float*          getCoordinates();
        virtual void    setCoordinates(float, float, float) = 0;
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
        float getMagnitude();
        //--------------------------------------------------------------------//

        // COORDINATES GETTERS/SETTERS
        //--------------------------------------------------------------------//
        void  setCoordinates(float, float, float);
        //--------------------------------------------------------------------//
    };
};
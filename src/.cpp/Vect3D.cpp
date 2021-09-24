/**
 * @file                Vect3D.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/VectND.h"
#include "../.h/Vect3D.h"

namespace TPP_VMath
{
    #pragma region Vect3D_Implementation
    //========================================================================//
    // START OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//

    Vect3D::Vect3D() : Vect3D({0, 0, 0}) 
    {   }

    Vect3D::Vect3D(std::array<float, 3> a) : VectND(3, a.begin())
    {   }

    Vect3D::Vect3D(const Vect3D& original) : VectND(3, original.components)
    {   }

    Vect3D::Vect3D(float* components) : VectND(3, components)
    {   }

    Vect3D  Vect3D::cross_product(Vect3D& that)
    {
        return Vect3D
        ({
           this->components[1] * that.components[2] - 
           this->components[2] * that.components[1],

           this->components[2] * that.components[0] - 
           this->components[0] * that.components[2],

           this->components[0] * that.components[1] - 
           this->components[1] * that.components[0]
        });
    }

    Vect3D  Vect3D::cross_product(Vect3D& a, Vect3D& b)
    {
        return Vect3D
        ({
           a.components[1] * b.components[2] - 
           a.components[2] * b.components[1],

           a.components[2] * b.components[0] - 
           a.components[0] * b.components[2],

           a.components[0] * b.components[1] - 
           a.components[1] * b.components[0]
        });
    }

    //========================================================================//
    // END OF VECT3D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect3D_Implementation
}
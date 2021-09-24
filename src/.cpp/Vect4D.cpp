/**
 * @file                Vect4D.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/VectND.h"
#include "../.h/Vect4D.h"

namespace TPP_VMath
{
    #pragma region Vect4D_Implementation
    //========================================================================//
    // START OF VECT4D CLASS IMPLEMENTATION
    //========================================================================//

    Vect4D::Vect4D() : Vect4D({0, 0, 0, 0})
    {   }

    Vect4D::Vect4D(std::array<float, 4> a) : VectND(4, a.begin())
    {   }

    Vect4D::Vect4D(const Vect4D& original) : VectND(4, original.components)
    {   }

    Vect4D::Vect4D(float* components) : VectND(4, components)
    {   }

    //========================================================================//
    // END OF VECT4D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect4D_Implementation
}
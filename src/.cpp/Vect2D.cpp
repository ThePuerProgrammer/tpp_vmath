/**
 * @file                Vect2D.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/VectND.h"
#include "../.h/Vect2D.h"

namespace TPP_VMath
{
    #pragma region Vect2D_Implementation
    //========================================================================//
    // START OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//

    Vect2D::Vect2D() : Vect2D({0, 0})
    {   }

    Vect2D::Vect2D(std::array<float, 2> a) : VectND(2, a.begin())
    {   }

    Vect2D::Vect2D(const Vect2D& original) : VectND(2, original.components)
    {   }

    Vect2D::Vect2D(float* components) : VectND(2, components)
    {   }

    //========================================================================//
    // END OF VECT1D CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion Vect2D_Implementation
}
/**
 * @file                VectND.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/VectND.h"

namespace TPP_VMath
{
    #pragma region VectND_Implementation
    //========================================================================//
    // START OF VECTND CLASS IMPLEMENTATION
    //========================================================================//

    VectND::VectND()
    {   
        dimension = 0;
        components = nullptr;
    }

    VectND::VectND(int n)
    {
        dimension = n;
        components = new float[n];

        for (int i = 0; i < n; ++i)
        {
            components[i] = 0;
        }
    }

    VectND::VectND(int n, float* components)
    {
        dimension = n;
        this->components = new float[n];

        for (int i = 0; i < n; ++i)
        {
            this->components[i] = components[i];
        }
    }

    VectND::VectND(std::vector<float> v)
    {
        if (v.size() > 0)
        {
            dimension = v.size();
            components = new float[v.size()];

            for (int i = 0; i < dimension; ++i)
            {
                components[i] = v[i];
            }
        }
        else
        {
            dimension = 0;
            components = nullptr;
        }
    }

    VectND::~VectND()
    {
        if (components) 
        {
            delete [] components;
            components = nullptr;
        }
    }

    //========================================================================//
    // END OF VECTND CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion VectND_Implementation
}
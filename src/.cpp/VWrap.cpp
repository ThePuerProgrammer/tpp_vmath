/**
 * @file                VWrap.cpp
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#include "../.h/Vect.h"
#include "../.h/VWrap.h"

namespace TPP_VMath
{
    #pragma region VWrap_Implementation
    //========================================================================//
    // START OF VWRAP CLASS IMPLEMENTATION
    //========================================================================//

    VWrap::VWrap()
    {
        n = 0;
        this->vect = nullptr;
        this->vpp = nullptr;
    }

    VWrap::VWrap(Vect* vect)
    {
        n = 1;
        this->vect = vect;
        this->vpp = nullptr;
    }

    VWrap::VWrap(int n, Vect** vpp)
    {
        this->n = n;
        this->vect = nullptr;
        this->vpp = vpp;
    }

    VWrap::~VWrap()
    {
        if (this->vect)
        {
            delete vect;
        }
        else if (this->vpp)
        {
            for (int i = 0; i < n; ++i)
            {
                delete vpp[i];
            }
            delete [] vpp;
        }
    }

    void    VWrap::wrap(Vect* vect)
    {
        n = 1;
        if (!this->vect)
            this->vect = vect;
    }

    void    VWrap::wrap(int n, Vect** vpp)
    {
        this->n = n;
        if (!this->vpp)
            this->vpp = vpp;
    }

    void    VWrap::unwrap()
    {
        if (this->vect)
        {
            delete this->vect;
            this->vect = nullptr;
        }
        else if (this->vpp)
        {
            for (int i = 0; i < n; ++i)
            {
                delete this->vpp[i];
                this->vpp[i] = nullptr;
            }
            delete [] this->vpp;
            this->vpp = nullptr;
        }
    }

    Vect*   VWrap::get_vect()
    {
        return vect;
    }

    Vect**  VWrap::get_vpp()
    {
        return vpp;
    }

    int     VWrap::get_n()
    {
        return n;
    }

    //========================================================================//
    // END OF VWRAP CLASS IMPLEMENTATION
    //========================================================================//
    #pragma endregion VWrap_Implementation
}
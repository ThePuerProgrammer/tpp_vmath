/**
 * @file                VWrap.h
 * @author              Jesse Rankins (https://www.github.com/thepuerprogrammer)
 * @brief               This class should be used to wrap any newly allocated 
 *                      Vect memory such as Vect* p = new Vect3D or any function
 *                      that returns a Vect* to a newly allocated resource.
 * 
 *                      Wrapping a new Vect allows the destructor to handle 
 *                      deallocation of heap resources. If two VWraps are 
 *                      used to wrap the same Vect*, this will result in 
 *                      undefined behavior. Therefore, a single address should 
 *                      only ever be wrapped once.
 * @version             0.1
 * @date                2021-09-23
 * @copyright           Copyright (c) 2021
 */

#ifndef VWRAP_H
#define VWRAP_H

#include "Vect.h"

namespace TPP_VMath
{
    #pragma region VWrap_Declaration
    //========================================================================//
    // START OF VWRAP CLASS DECLARATION
    //========================================================================//

    class VWrap
    {
    public:

        /**
         * @brief       Set vect to nullptr and wrap as a seperate step
         * @param       void
         */ 
        VWrap();

        /**
         * @brief       Wrap upon construction
         * @param       vect a polymorphic pointer to a Vect child
         */ 
        VWrap(Vect*);

        /**
         * @brief       Wrap array upon construction
         * @param       n the number of wrapped entries
         * @param       vpp a polymorphic pointer to an array Vect children
         */ 
        VWrap(int, Vect**);

        /**
         * @brief       Destructor handles deallocation of vect
         */ 
        ~VWrap();

        /**
         * @brief       If vect is null, wrap a Vect*
         * @param       vect a polymorphic pointer to a Vect child
         */  
        void            wrap(Vect*);

        /**
         * @brief       If vpp is null, wrap a Vect**
         * @param       vpp a polymorphic pointer to an array of Vect children
         */  
        void            wrap(int, Vect**);

        /**
         * @brief       Explicit deallocation of vect
         * @param       void
         */ 
        void            unwrap();

        /**
         * @brief       Since address is wrapped, the pointer can be assigned
         * @param       void
         * @return      the wrapped pointer
         */ 
        Vect*           get_vect();

        /**
         * @brief       Since address is wrapped, the pointer can be assigned
         * @param       void
         * @return      the wrapped vpp
         */ 
        Vect**          get_vpp();

        /**
         * @return      the number of wrapped entries
         */ 
        int             get_n();

    protected:

        ///             The number of wrapped entries
        int             n;

        ///             A wrapped Vect*
        Vect*           vect;

        ///             A wrapped Vect**
        Vect**          vpp;
    };

    //========================================================================//
    // END OF VWRAP CLASS DECLARATION
    //========================================================================//
    #pragma endregion VWrap_Declaration
}

#endif
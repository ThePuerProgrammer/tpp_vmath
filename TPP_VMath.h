/**
 * @file        TPP_VMath.h
 * @author      Jesse Rankins
 * @since       2021
 * @version     0.1
 */ 
#pragma once

#include <array>
#include <vector>

namespace TPP_VMath
{
    #pragma region Vect_Declaration
    //========================================================================//
    // START OF VECT CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           Parent class for all vectors in the library. Children 
     *                  that inherit from Vect should be thought of as having a 
     *                  vertical orientation when considering their relationship
     *                  to VSets and Matrices. 
     * @file            TPP_VMath.h
     */ 
    class Vect
    {
    public:

        /**
         * @brief       The magnitude of a Vect is its distance from the origin
         * @param       void
         * @return      sqrtf(v₁^2+...+vᵢ^2) as a float
         */
        float           get_magnitude();

        /**
         * @brief       The magnitude of a Vect is its distance from the origin
         * @param       vect a reference to a child of Vect 
         * @return      sqrtf(v₁^2+...+vᵢ^2) as a float
         */
        static float    get_magnitude(Vect&);

        /**
         * @brief       Scales each entry in the vector by integer scalar value
         * @param       c an integer scalar value
         */ 
        void            scale_by(int);

        /**
         * @brief       Scales each entry in the vector by integer scalar value
         * @param       c an integer scalar value
         * @param       vect a reference to a child of Vect
         */ 
        static void     scale_by(int, Vect&);

        /**
         * @brief       Scales each entry in the vector by float scalar value
         * @param       c a floating point scalar value
         */ 
        void            scale_by(float);

        /**
         * @brief       Scales each entry in the vector by float scalar value
         * @param       c an float scalar value
         * @param       vect a reference to a child of Vect
         */ 
        static void     scale_by(float, Vect&);

        /**
         * @brief       If the Vects have unmatching dimensions an error will
         *              print to the console and the function will return 
         *              -FLT_MAX
         * @param       that a const Vect reference
         * @return      sum of the entrywise products of two Vects
         */ 
        float           dot_product(const Vect&);

        /**
         * @brief       If the Vects have unmatching dimensions an error will
         *              print to the console and the function will return 
         *              -FLT_MAX
         * @param       a a const Vect reference
         * @param       b a const Vect reference
         * @return      sum of the entrywise products of two Vects
         */ 
        static float    dot_product(const Vect&, const Vect&);

        /**
         * @brief       Normalization scales the vector by the inverse of its
         *              magnitude, which has the effect of setting the magnitude
         *              of the vector to 1.
         * @param       void
         */
        void            normalize_vect(); 

        /**
         * @brief       Normalization scales the vector by the inverse of its
         *              magnitude, which has the effect of setting the magnitude
         *              of the vector to 1.
         * @param       vect a reference to a child of Vect
         */
        static void     normalize_vect(Vect&); 

        /**
         * @brief       Manual assignment of entires in the vector
         * @param       v the provided array of floats
         */ 
        void            set_coordinates(std::vector<float>);

        /**
         * @param       i the ith coordinate
         * @return      a pointer to the ith coordinate in the vector
         */
        float*          operator[](int);

        /**
         * @param       void
         * @return      the coordinate array as a pointer to coordinates[0]
         */
        float*          get_coordinates() const;

        /**
         * @param       void
         * @return      the dimension of the vector (the number of entries)
         */ 
        unsigned int    get_dimension() const;

        /**
         * @brief       Copies coordinates from right vector to left vector
         * @param       right a const Vect reference to the right side of =
         */ 
        Vect&           operator=(const Vect&);

        /**
         * @brief       Pure virtual destructor
         */ 
        virtual         ~Vect() = 0;

    protected:

        ///             The number of entries in coordinates
        unsigned int    dimension;

        ///             The array representing the vector
        float*          coordinates;
    };

    //========================================================================//
    // END OF VECT CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect_Declaration

    #pragma region VectND_Declaration
    //========================================================================//
    // START OF VECTND CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           A floating point vector in RN (variable length) that
     *                  inherits from Vect
     * @file            TPP_VMath.h
     */
    class VectND : public Vect
    {
    public:

        VectND();

        /**
         * @brief       Constructor for a vector in RN where N is provided that 
         *              defaults to the zero vector
         * @param       n the dimension of the vector
         */ 
        VectND(int);

        /**
         * @brief       Constructor for a vector in RN where N is provided that
         *              initializes to the provided array of floats
         * @param       n the dimension of the vector
         * @param       coordinates the provided array of floats
         */
        VectND(int, float*);

        /**
         * @brief       Constructor for a vector in RN where N is the size of a
         *              std::vector<float>
         * @param       v a std::vector of floats
         */
        VectND(std::vector<float>); 

        /**
         * @brief       Destructor deletes the coordinate array
         */
        ~VectND();
    }; 

    //========================================================================//
    // END OF VECTND CLASS DECLARATION
    //========================================================================//
    #pragma endregion VectND_Declaration

    #pragma region Vect1D_Declaration
    //========================================================================//
    // START OF VECT1D CLASS DECLARATION
    //========================================================================//
    /**
     * @brief           A one dimensional vector can be thought of as a scalar 
     *                  but there are some unique edge cases where this is a
     *                  very interesting and useful property. A floating point
     *                  vector in R1 that inherits from VectND
     * @file            TPP_VMath.h
     */
    class Vect1D : public VectND
    {
    public:

        /**
         * @brief       constructor for a vector in R1 that inits to zero
         * @param       void
         */ 
        Vect1D();

        /**
         * @brief       constructor for a vector in R1 that inits to (x)
         * @param       a an array of floats representing the value(s) x
         */ 
        Vect1D(std::array<float, 1>);

        /**
         * @brief       Copy constructor for a vector in R1
         * @param       original a const Vect1D reference
         */
        Vect1D(const Vect1D&); 

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect1D(); 
         *              Error prone!! Call should ALWAYS be wrapped with guard 
         *              in the form if (v->get_dimension() == 1)
         * @param       coordinates the entries in the 1D vector
         */
        Vect1D(float*);
    }; 
    //========================================================================//
    // END OF VECT1D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect1D_Declaration

    #pragma region Vect2D_Declaration
    //========================================================================//
    // START OF VECT2D CLASS DECLARATION
    //========================================================================//
    /**
     * @brief           A floating point vector in R2 that inherits from Vect
     * @file            TPP_VMath.h
     */
    class Vect2D : public VectND
    {
    public:

        /**
         * @brief       constructor for a vector in R2 that inits to zero
         * @param       void
         */ 
        Vect2D();

        /**
         * @brief       constructor for a vector in R2 that inits to (x,y)
         * @param       a an array of floats representing the value(s) x,y
         */ 
        Vect2D(std::array<float, 2>);

        /**
         * @brief       Copy constructor for a vector in R2
         * @param       original a const Vect2D reference
         */
        Vect2D(const Vect2D&); 

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect2D(); 
         *              Error prone!! Call should ALWAYS be wrapped with guard 
         *              in the form if (v->get_dimension() == 2)
         * @param       coordinates the entries in the 2D vector
         */
        Vect2D(float*);
    }; 
    //========================================================================//
    // END OF VECT2D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect2D_Declaration

    #pragma region Vect3D_Declaration
    //========================================================================//
    // START OF VECT3D CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           A floating point vector in R3 that inherits from Vect
     * @file            TPP_VMath.h
     */ 
    class Vect3D : public VectND
    {
    public:
    
        /**
         * @brief       Constructor for a vector in R3 that inits to zero
         * @param       void
         */ 
        Vect3D();

        /**
         * @brief       Constructor for a vector in R3 that inits to (x,y,z)
         * @param       a an array of floats representing the value(s) x,y,z
         */
        Vect3D(std::array<float, 3>);

        /**
         * @brief       Copy Constructor for a vector in R3
         * @param       original a const Vect3D reference
         */ 
        Vect3D(const Vect3D&);

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect3D(); 
         *              Error prone!! Call should ALWAYS be wrapped with guard 
         *              in the form if (v->get_dimension() == 3)
         * @param       coordinates the entries in the 3D vector
         */
        Vect3D(float*);
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect3D_Declaration

    #pragma region Vect4D_Declaration
    //========================================================================//
    // START OF VECT4D CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           A floating point vector in R4 that inherits from Vect
     * @file            TPP_VMath.h
     */ 
    class Vect4D : public VectND
    {
    public:

        /**
         * @brief       A vector in R4 that defaults to zero
         */
        Vect4D(); 

        /**
         * @brief       A vector in R4 that inits to (x, y, z, w)
         * @param       a an array of floats representing the value(s) x,y,z,w
         */
        Vect4D(std::array<float, 4>);

        /**
         * @brief       Copy Constructor for a vector in R4
         * @param       original a const Vect4D reference
         */
        Vect4D(const Vect4D&);

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect4D(); 
         *              Error prone!! Call should ALWAYS be wrapped with guard 
         *              in the form if (v->get_dimension() == 4)
         * @param       coordinates the entries in the 4D vector
         */
        Vect4D(float*); 
    };
    //========================================================================//
    // END OF VECT4D CLASS DECLARATION
    //========================================================================//
    #pragma endregion Vect4D_Declaration

    #pragma region VWrap_Declaration
    //========================================================================//
    // START OF VWRAP CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           This class should be used to wrap any newly allocated 
     *                  Vect memory such as Vect* p = new Vect3D[3] or any 
     *                  function that returns a Vect* to a newly allocated 
     *                  resource, such as the function 
     * 
     *                  Vect* Matrix::get_matrix_vector_product() 
     * 
     *                  Wrapping a new Vect allows the destructor to handle 
     *                  deallocation of heap resources. If two VWraps are 
     *                  used to wrap the same Vect*, this will result in 
     *                  undefined behavior. Therefore, a single address should 
     *                  only ever be wrapped once.
     * @file            TPP_VMath.h
    */
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

    #pragma region VSet_Declaration
    //========================================================================//
    // START OF VSET CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           A set of vectors with matching dimensions
     * @file            TPP_VMath.h
     */ 
    class VSet
    {
    public:

        /**
         * @brief       Default constructor is a null set
         * @param       void
         */ 
        VSet();

        /**
         * @brief       Constructor for a set containing a single Vect child
         * @param       v a Vect child reference
         */
        VSet(Vect&);

        /**
         * @brief       Constructor using an array of Vect where n is the number
         *              of elements in a. In order to avoid errors, a should be 
         *              instatntiated with a const int that is also used as the 
         *              argument for n. This function should not be used to 
         *              create a set from a wrapped Vect**. Instead, a
         *              constructor has been defined to accept a VWrap.
         * @param       n the number of elements in a
         * @param       a a polymorphic pointer to an array of Vect children
         */
        VSet(int, Vect**); 

        /**
         * @brief       Constructor using a VWrap that can either contain a
         *              single wrapped Vect* or an array of Wrapped Vect*. This
         *              function creates a copy of every wrapped Vect.
         * @param       wrapper a VWrap containing one or more Vects
         */
        VSet(VWrap&); 

        /**
         * @brief       Copy constructor
         * @param       original a const VSet reference
         */ 
        VSet(const VSet&);

        /**
         * @brief       Destructor deletes the set
         */ 
        ~VSet();

        /**
         * @brief       The number of vectors in the set is always non-negative
         * @param       void
         * @return      The number of vectors in the set
         */ 
        unsigned int    get_n() const;
        
        /**
         * @brief       The number of entries is always non-negative
         * @param       void
         * @return      The number of entries in each vector
         */ 
        unsigned int    get_m() const;

        /** 
         * @brief       Returns the set as an array of Vects
         * @param       void
         * @return      an array of Vect* that point to children of Vect
         */ 
        Vect**          get_set_of_vectors() const;

        /**
        * @brief        Add a new vect of *matching dimensions* to the set. If 
        *               the vector dimension isn't m, a runtime exception is 
        *               thrown. It is best to wrap the call to this funtion in 
        *               a guard ensuring that the dimensions are a match. This 
        *               does not apply if the set is currently null.
        * @param        vect a polymorphic pointer to a child of Vect
        */
        void            add_vect_to_set(Vect*);

    private:

        ///             The number of entries in each vector in the set
        unsigned int    m;

        ///             The number of vectors in the set
        unsigned int    n;

        ///             The set of vectors as an array in pointer notation
        Vect**          setOfVectors;

        /**
         * @brief       Called by add_vect_to_set() to reduce redundancy
         * @param       m the dimension of the vect
         * @param       vect the original vect passed to add_vect_to_set()
         * @return      a Vect* to a new child vect of m dimensions
         */ 
        Vect*           add_vect_of_valid_dimensions(int, Vect*);
    };

    //========================================================================//
    // END OF VSET CLASS DECLARATION
    //========================================================================//
    #pragma endregion VSet_Declaration

    #pragma region Matrix_Declaration
    //========================================================================//
    // START OF MATRIX CLASS DECLARATION
    //========================================================================//

    class Matrix
    {
    public:

        /**
         * @brief       Overloaded constructor accepts a set of Vects and 
         *              generates a matrix
         * @param       set a const VSet reference as a set of vectors with 
         *              matching dimensions
         */ 
        Matrix(const VSet&);

        /**
         * @brief       Overloaded constructer accepts a float** in the form
         *              f[rows][cols]
         * @param       fMatrix a 2D array of floats
         * @param       m the number of rows in the matrix
         * @param       n the number of columns in the matrix
         */ 
        Matrix(int, int, float**);

        /**
         * @brief       Matrix copy constructor
         * @param       original a const Matrix reference
         */ 
        Matrix(const Matrix&);

        /**
         * @brief       Destructor deletes the matrix A
         */ 
        ~Matrix();

        /**
         * @brief       Reduce this matrix to reduced echelon form
         * @param       void
         */ 
        void            reduce_matrix();

        /**
         * @brief       Produce a new reduced echelon matrix from this matrix
         * @param       void
         * @return      this matrix in RE form
         */ 
        Matrix          get_reduced();

        /**
         * @brief       Ax = b
         * @param       vect a polymorphic pointer to a child of Vect
         * @return      a pointer to b where b = Ax
         */ 
        Vect*           get_matrix_vector_product(Vect*);

        /**
         * @brief       Console representation of coefficient matrix for testing
         * @param       void
         */ 
        void            print_matrix();

        /**
         * @brief       Identity matrix is an nxn diagonal matrix in RE form
         * @param       n the number of rows and columns in the resulting matrix
         * @return      A newly created nxn identity matrix
         */
        static Matrix   get_identity_matrix_of_size(int);

    private:

        ///             The number of rows in A
        unsigned int    m;

        ///             The number of columns in A
        unsigned int    n;

        ///             2D array using pointer notation representing a matrix
        float**         A;

        /**
         * @brief       Swaps the pointers for A[rowA] and A[rowB]
         * @param       rowA a value representing a row in A
         * @param       rowB a value representing a row in A
        */
        void            row_interchange(int rowA, int rowB);

        /**
         * @brief       Generate a zero vector in Rm
         * @param       m the number of entries in the vector
         * @return      a new instance of a Vect child
         */
        Vect*           get_b(const int&);
    };

    //========================================================================//
    // END OF MATRIX CLASS DECLARATION
    //========================================================================//
    #pragma endregion Matrix_Declaration
};
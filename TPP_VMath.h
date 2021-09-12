/**
 * @file        TPP_VMath.h
 * @author      Jesse Rankins
 * @since       2021
 * @version     0.1
 */ 
#pragma once

namespace TPP_VMath
{

    //========================================================================//
    // START OF VECT CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           Parent class for all vectors in the library
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
         * @brief       Scales each entry in the vector by integer scalar value
         * @param       c an integer scalar value
         */ 
        void            scale_by(int);

        /**
         * @brief       Scales each entry in the vector by float scalar value
         * @param       c a floating point scalar value
         */ 
        void            scale_by(float);

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

    //========================================================================//
    // START OF VWRAPPER CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           This class should be used to wrap any function that 
     *                  returns a Vect* to a newly allocated resource, such as 
     *                  the function 
     * 
     *                  Vect* Matrix::get_matrix_vector_product() 
     * 
     *                  Wrapping a new Vect allows the destructor to handle 
     *                  deallocation of heap resources. If two VWrappers are 
     *                  used to wrap the same Vect*, this will result in 
     *                  undefined behavior. Therefore, a single address should 
     *                  only ever be wrapped once.
     * @file            TPP_VMath.h
    */
    class VWrapper
    {
    public:

        /**
         * @brief       Set vect to nullptr and wrap as a seperate step
         * @param       void
         */ 
        VWrapper();

        /**
         * @brief       Wrap upon construction
         * @param       vect a polymorphic pointer to a Vect child
         */ 
        VWrapper(Vect*);

        /**
         * @brief       Destructor handles deallocation of vect
         */ 
        ~VWrapper();

        /**
         * @brief       If vect is null, wrap a Vect*
         * @param       vect a polymorphic pointer to a Vect child
         */  
        void            wrap(Vect*);

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

    protected:

        ///             The wrapped Vect*
        Vect*           vect;
    };

    //========================================================================//
    // END OF VWRAPPER CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT1D CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECT1D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT2D CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECT2D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT3D CLASS DECLARATION
    //========================================================================//

    /**
     * @brief           A vector in R3 that inherits from Vect
     * @file            TPP_VMath.h
     */ 
    class Vect3D : public Vect
    {
    public:
    
        /**
         * @brief       Constructor for a vector in R3 that inits to zero
         * @param       void
         */ 
        Vect3D();

        /**
         * @brief       Constructor for a vector in R3 that inits to (x,y,z)
         * @param       x the 1st entry in the Vect3D
         * @param       y the 2nd entry in the Vect3D
         * @param       z the 3rd entry in the Vect3D
         */
        Vect3D(const float, const float, const float);

        /**
         * @brief       Copy Constructor
         * @param       original a const Vect3D reference
         */ 
        Vect3D(const Vect3D&);

        /**
         * @brief       Constructor useful for polymorphic copies of 
         *              Vect* v = new Vect3D(); 
         *              Error prone!! Call should ALWAYS be wrapped with guard 
         *              in the form if (v->get_dimension() == 3)
         * @param       coordinates the float pointer representing the entries
         *              in the vector
         */
        Vect3D(const float*);

        /**
         * @brief       Destructor deletes coordinate array
         */ 
        ~Vect3D();

        /**
         * @brief       Copies coordinates from right vector to left vector
         * @param       right a const Vect3D reference to the right side of =
         */ 
        Vect3D&         operator=(const Vect3D&);

        /**
         * @brief       Manual assignment of x, y, z vector entries
         * @param       x the 1st entry in the Vect3D
         * @param       y the 2nd entry in the Vect3D
         * @param       z the 3rd entry in the Vect3D
         */ 
        void            set_coordinates(float, float, float);

        /**
         * @param       that a const Vect3D
         * @return      sum of the entrywise products of two Vect3Ds
         */ 
        float           dot_product(const Vect3D&);
    };

    //========================================================================//
    // END OF VECT3D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECT4D CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECT4D CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF VECTND CLASS DECLARATION
    //========================================================================//
    // TODO
    // TODO
    // TODO
    // TODO
    //========================================================================//
    // END OF VECTND CLASS DECLARATION
    //========================================================================//

    //========================================================================//
    // START OF SET CLASS DECLARATION
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
         * @brief       Overloaded constructor for one Vect3D
         * @param       v3D a const Vect3D reference is a vector in R3
         */ 
        VSet(const Vect3D&);

        /** 
         * @brief       Overloaded VSet constructor using an array of Vect3D 
         *              where n is the number of elements in a. In order to 
         *              avoid errors, a should be instantiated with a const int 
         *              that is also used as the argument for n.
         * @param       n the number of elements in a
         * @param       a an array of Vect3D representing a set
        */
        VSet(const int&, Vect3D[]);

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
         * @brief       Returns the number of vectors in the set
         * @param       void
         * @return      the number of vectors in the set is always non-negative
         */ 
        unsigned int    get_n() const;
        
        /**
         * @brief       Returns the number of entries in each vector
         * @param       void
         * @return      the number of entries is always non-negative
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
    };

    //========================================================================//
    // END OF SET CLASS DECLARATION
    //========================================================================//

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
    };

    //========================================================================//
    // END OF MATRIX CLASS DECLARATION
    //========================================================================//

};
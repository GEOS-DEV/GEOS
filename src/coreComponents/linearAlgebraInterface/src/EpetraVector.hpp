/**
 * @file EpetraVector.hpp
 */

#ifndef EPETRAVECTOR_HPP_
#define EPETRAVECTOR_HPP_

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * \class EpetraVector
 * \brief This class creates and provides basic support for the Epetra_Vector
 *        vector object type used in Trilinos.
 */

using integer     = int;
using localIndex  = int;

class EpetraVector
{
public:
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty vector constructor.
   *
   * Create an empty (distributed) vector.
   */
  EpetraVector();

  /**
   * @brief Copy constructor.
   */
  EpetraVector( EpetraVector const & vector );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraVector() = default;
  //@}

  //! @name Create Methods
  //@{we

  /**
   * @brief Construct vector from array.
   *
   * Create a vector from an Epetra_Map and an array of values.
   */
  void create( const Epetra_Map &map,
               real64     *values );

  /**
   * @brief Construct vector from array.
   *
   * Create a vector from an Epetra_Map and an array of values.
   */
  void create( const globalIndex  size,
               real64      *values );

  /**
   * @brief Construct vector from std::vector.
   *
   * Create a vector from an std vector.
   */
  void create( std::vector<real64> &vec );
  //@}

  //! @name Linear Algebra Methods
  //@{

  /**
   * @brief Multiply all elements by scalingFactor.
   */
  void scale( real64 const scalingFactor );

  /**
   * @brief Dot product with the vector vec.
   */
  void dot( EpetraVector const &vec,
            real64 *dst );

  /**
   * @brief Update (name to be changed) vector as this = alpha*vec + beta*this.
   */
  void update( real64 const alpha,
               EpetraVector const &vec,
               real64 const beta );

  /**
   * @brief 1-norm of the vector.
   */
  void norm1( real64 &result ) const;

  /**
   * @brief 2-norm of the vector.
   */
  void norm2( real64 &result ) const;

  /**
   * @brief Infinity-norm of the vector.
   */
  void normInf( real64 &result ) const;

  //@}

  //! @name Accessor Methods
  //@{

  /**
   * @brief Returns the global of the vector.
   */
  globalIndex globalSize() const;

  /**
   * @brief Returns the local of the vector.
   */
  localIndex localSize() const;

  /**
   * @brief Returns a const pointer to the underlying Epetra_Vector.
   */
  const Epetra_Vector* getPointer() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra_Vector.
   */
  Epetra_Vector* getPointer();

  //@}

  //! @name I/O Methods
  //@{

  /**
   * @brief Print the vector in Trilinos format to the terminal.
   */
  void print() const;

  //@}

protected:
  // Unique pointer to underlying Epetra_Vector type.
  std::unique_ptr<Epetra_Vector> m_vector = nullptr;
};

}

#endif /* EPETRAVECTOR_HPP_ */

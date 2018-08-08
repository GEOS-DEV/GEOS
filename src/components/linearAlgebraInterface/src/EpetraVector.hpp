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
  void create(const Epetra_Map &map,
      real64     *values);

  /**
   * @brief Construct vector from array.
   *
   * Create a vector from an Epetra_Map and an array of values.
   */
  void create(const globalIndex  size,
                    real64      *values);

  /**
   * @brief Construct vector from std::vector.
   *
   * Create a vector from an std vector.
   */
  void create(std::vector<real64> &vec);
  //@}


  //! @name Accessor Methods
  //@{

  /**
   * @brief Returns the global of the vector.
   */
  globalIndex globalSize();

  /**
   * @brief Returns the local of the vector.
   */
  localIndex localSize();

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
  void print();

  //@}

protected:
  // Pointer to underlying Epetra_Vector type.
  Epetra_Vector *vector;
};

}

#endif /* EPETRAVECTOR_HPP_ */

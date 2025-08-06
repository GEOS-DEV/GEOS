/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TableLayout.hpp
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP

#include "common/DataTypes.hpp"
#include "TableTypes.hpp"
#include <variant>
#include "common/logger/Logger.hpp"


namespace geos
{

/**
 * @brief Class for setup the table layout
 */
class TableLayout
{

public:

  /// default value for m_maxColumnWidth when it is not set
  static constexpr size_t noColumnMaxWidth = std::numeric_limits< size_t >::max();

  /// Type of aligment for a column
  enum Alignment { right, left, center };

  /// default value for columns header cells alignement
  static constexpr Alignment defaultHeaderAlignment = Alignment::center;

  /// default value for data cells alignement
  static constexpr Alignment defaultValueAlignment = Alignment::right;

  /// Space to apply between all data and border
  enum MarginValue : integer
  {
    tiny = 0,
    small = 1,
    medium = 2,
    large = 3
  };

  /**
   * @brief Enumeration for table sections.
   */
  enum Section { header, values };

  /**
   * @brief Structure to set up values m_alignment for each colum.
   */
  struct ColumnAlignement
  {
    /// Alignment for column name. By default aligned to center
    Alignment headerAlignment = defaultHeaderAlignment;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment = defaultValueAlignment;
  };

  /**
   * @brief View on cell data with information to display it in a table (content, type, alignment, ...).
   * @note the source text must not be freeed/moved as the CellLayout will *only be a view on the text data*.
   */
  class CellLayout
  {
public:
    /// The type of the cell (Header,Value, Merge, ...).
    CellType m_cellType;
    /// The alignment of the cell (left, center, right).
    Alignment m_alignment;

    /**
     * @brief Constructor to initialize a Cell with a default settings. Use prepareLayout() when setup.
     */
    CellLayout();

    /**
     * @brief Constructor to initialize an empty Cell of a given type.
     * @param cellType The type of the cell.
     */
    CellLayout( CellType cellType );

    /**
     * @brief Constructor to fully initialize a cell with given celltype, text and alignment.
     * m_cellWidth will be initialized aDter
     * @param cellType The type of the cell.
     * @param alignment The alignment of the cell (left, right, or center).
     */
    CellLayout( CellType cellType, TableLayout::Alignment alignment );

    /**
     * @return The width of the cell, which must be constrained by the content lines length.
     */
    size_t getWidth() const
    { return m_cellWidth; }

    /**
     * @brief Set the width of the cell, which must be constrained by the content lines length.
     * @param cellWidth the new width to consider for this cell.
     */
    void setWidth( size_t cellWidth )
    { m_cellWidth = cellWidth; }

    /**
     * @return The view on each cell line.
     */
    stdVector< string_view > const & getLines() const
    { return m_lines; }

    /**
     * @return get the height of the cell (its number of lines).
     */
    size_t getHeight() const
    { return m_lines.size(); }

    /**
     * @return True if the cell has no text.
     */
    size_t isEmpty() const
    { return m_lines.empty() || m_lines[0].empty(); }

    /**
     * @brief Set the data view to the given string_view & precompute display settings.
     * @param value The view on the full cell text, with '\n' for manual line breaks. Must not be deallocated!
     *              `getLines()` will then contain each lines, and `m_cellWidth`, the maximum line width.
     * @param maxLineWidth The maximum allowed line width. Use `noColumnMaxWidth` to disable.
     */
    void prepareLayout( string_view value, size_t maxLineWidth );

private:
    /// The width of the cell, which must be constrained by the content lines length.
    size_t m_cellWidth;
    /// vector containing each cell content, separated by lines.
    stdVector< string_view > m_lines;
  };

  /**
   * @brief Represents a cell in a table with ownership of its text data.
   * @note Unlike CellLayout which is just a view, Cell owns its text data (m_text).
   *       For safety & performance reasons, Cell is non-copyable/movable after layout has been prepared.
   */
  class Cell
  {
public:
    /// The view & display setting on m_text.
    CellLayout m_layout;

    /**
     * @brief Constructor to initialize a Cell with a default settings. Use prepareLayout() after setup.
     */
    Cell();

    /**
     * @brief Constructor to partially initialize a cell with display settings. Use prepareLayout() after setup.
     * @param cellType The type of the cell.
     * @param alignment The alignment of the cell (left, right, or center).
     */
    Cell( CellType cellType, TableLayout::Alignment alignment );

    /**
     * @brief Constructor to partially initialize a cell with all settings. Use prepareLayout() after setup.
     * @param cellType The type of the cell.
     * @param alignment The alignment of the cell (left, right, or center).
     * @param value The text to set in the cell (stored in m_text).
     */
    Cell( CellType cellType, TableLayout::Alignment alignment, string_view value );

    /**
     * @brief Copy data, or throw an error if the layout has already been prepared (which means instance
     *        will reference potencially outdated reference, and we do not want to compute the layout twice)
     * @param other The source data.
     */
    Cell( Cell const & other );

    /**
     * @brief Move data, or throw an error if the layout has already been prepared (which means instance
     *        will reference potencially outdated reference, and we do not want to compute the layout twice)
     * @param other The source data.
     */
    Cell( Cell && other );

    /**
     * @brief Copy data, or throw an error if the layout has already been prepared (which means instance
     *        will reference potencially outdated reference, and we do not want to compute the layout twice)
     * @param other The source data.
     * @return The instance reference.
     */
    Cell & operator=( Cell const & other );

    /**
     * @brief Move data, or throw an error if the layout has already been prepared (which means instance
     *        will reference potencially outdated reference, and we do not want to compute the layout twice)
     * @param other The source data.
     * @return The instance reference.
     */
    Cell & operator=( Cell && other );

    /**
     * @return The full cell text.
     */
    string_view getText() const
    { return m_text; }

    /**
     * @brief Set the full cell text.
     * @param text The full cell text, with '\n' for manual line breaks.
     */
    void setText( string_view text );

    /**
     * @brief Precompute m_layout display settings and link it with m_text.
     * @param maxLineWidth The maximum allowed line width. Use `noColumnMaxWidth` to disable.
     */
    void prepareLayout( size_t maxLineWidth );

private:
    /// The text data of the cell (potencially multiline).
    string m_text;
  };

  /**
   * @class Column
   * @brief Class representing a column in a table layout.
   */
  class Column
  {
public:
    /// Alias for the list of columns.
    using ColumnsList = stdVector< Column >;

    /// The header cell
    Cell m_header;
    /// A vector containing all sub-columns in the column.
    ColumnsList m_subColumns;
    /// struct containing m_alignment for the column (header and values)
    ColumnAlignement m_alignment;

    /**
     * @brief Construct a default column with no parameter (must be configurated).
     */
    Column();

    /**
     * @brief Construct a default column with minimal parameters.
     * @param name The name of the Column.
     */
    explicit Column( string_view name ):
      Column( name, ColumnAlignement() )
    {}

    /**
     * @brief Construct a default column with minimal parameters.
     * @param name The name of the Column.
     * @param alignment The alignment setting of the column header and values.
     */
    Column( string_view name, ColumnAlignement alignment );

    /**
     * @brief Get the parent column.
     * @return Pointer to the parent column, or `nullptr` if no parent is set.
     */
    Column * getParent()
    { return m_parent; }

    /**
     * @brief Get the parent column.
     * @return Pointer to the parent column, or `nullptr` if no parent is set.
     */
    Column const * getParent() const
    { return m_parent; }

    /**
     * @brief Set the parent column.
     * @param parent Pointer to the parent column to set.
     */
    void setParent( Column * parent )
    { m_parent = parent; }

    /**
     * @return Pointer to the next column that has the same parent or `nullptr` if no next column exists.
     */
    Column * getNext()
    { return m_next; }

    /**
     * @return Pointer to the next column that has the same parent or `nullptr` if no next column exists.
     */
    Column const * getNext() const
    { return m_next; }

    /**
     * @param nextCell The next column in the table layout that has the same parent.
     */
    void setNext( Column * nextCell )
    {  m_next = nextCell; }

    /**
     * @brief Sets the name of the column.
     * @param name The name to set for the column.
     * @return The current column object.
     */
    Column & setName( string_view name );

    /**
     * @brief Set the column and its children visibility.
     * @param visible True to make the column visible.
     * @return The current column .
     */
    Column & setVisibility( bool visible );

    /**
     * @brief Adds multiple sub-columns to the column.
     * @param subCol A list of sub-column names to add.
     * @return The current column object
     */
    Column & addSubColumns( std::initializer_list< Column > subCol );

    /**
     * @brief Adds multiple sub-columns to the column.
     * @param subColNames A list of sub-column names to add.
     * @return The current column object
     */
    Column & addSubColumns( std::initializer_list< string > subColNames );

    /**
     * @brief Adds multiple sub-columns to the column.
     * @param subColNames A list of sub-column names to add.
     * @return The current column object
     */
    Column & addSubColumns( stdVector< string > const & subColNames );

    /**
     * @brief Adds a single sub-column to the column.
     * @param subColName The name of the sub-column to add.
     * @return The current column object.
     */
    Column & addSubColumn( string_view subColName );

    /**
     * @brief Adds a single sub-column to the column.
     * @param subCol The sub-column to add.
     * @return The current column object.
     */
    Column & addSubColumn( Column const & subCol );

    /**
     * @brief Sets the header alignment for the column.
     * @param headerAlignment The alignment to set for the column header (left, right, or center).
     * @return The current column object
     */
    Column & setHeaderAlignment( Alignment headerAlignment );

    /**
     * @brief Sets the values alignment for the column.
     * @param valueAlignment The alignment to set for the column values (left, right, or center).
     * @return The current column object
     */
    Column & setValuesAlignment( Alignment valueAlignment );

    /**
     * @brief Checks if the column has any child columns.
     * @return bool True if the column has child columns, otherwise false.
     */
    bool hasChild() const
    { return !this->m_subColumns.empty(); }

    /**
     * @brief Checks if the column has a parent column.
     * @return bool True if the column has a parent, otherwise false.
     */
    bool hasParent() const
    { return this->m_parent != nullptr; }

    /**
     * @return bool True if the column has a neightboor to its right that has the same parent.
     */
    bool hasNext() const
    { return this->m_next != nullptr; }

    /**
     * @return True if the column and its children are visible.
     */
    bool isVisible() const
    { return m_header.m_layout.m_cellType!=CellType::Hidden; }

private:
    /// Pointer to the parent cell (if any).
    Column * m_parent = nullptr;
    /// Pointer to the next cell (if any).
    Column * m_next = nullptr;
  };

  /**
   * @brief Iterator to loop over all columns, starting by the deepest sub columns,
   * then to their parents, then to their siblings.
   */
  class DeepFirstIterator
  {
public:
    ///alias for column
    using ColumnType = Column const;

    /**
     * @brief Construct a new Leaf Iterator object
     * @param columnPtr The first deepest column of vector
     * @param idxLayer the layer associated with the column
     */
    DeepFirstIterator( ColumnType * columnPtr, size_t idxLayer ):
      m_currentColumn( columnPtr ), m_currentLayer( idxLayer )
    {}

    /**
     * @brief Copy assignment operator
     * @param[in] columnPtr Coulmn  to copy
     * @return Leaf iterator
     */
    DeepFirstIterator & operator=( ColumnType * columnPtr )
    {
      this->m_currentColumn = columnPtr;
      return *this;
    }

    /**
     * @brief Prefix ++ overload
     * @return Leaf iterator
     */
    DeepFirstIterator & operator++();

    /**
     * @brief Postfix ++ overload
     * @return Leaf iterator
     */
    DeepFirstIterator operator++( int );

    /**
     * @brief Dereference operator.
     * @return Reference to the current Column object pointed to by the iterator.
     */
    ColumnType & operator*() const
    { return *m_currentColumn; }

    /**
     * @return Pointer to the current Column object pointed to by the iterator.
     */
    ColumnType * getPtr() const
    { return m_currentColumn; }

    /**
     * @brief Arrow operator.
     * @return Pointer to the current Column object.
     */
    ColumnType * operator->() const
    { return m_currentColumn; }

    /**
     * @brief Equality comparison operator.
     * @param a The first iterator.
     * @param b The second iterator.
     * @return True if both iterators point to the same column; false otherwise.
     */
    friend bool operator== ( DeepFirstIterator const & a, DeepFirstIterator const & b )
    { return a.m_currentColumn == b.m_currentColumn; };
    /**
     * @brief Inequality comparison operator.
     * @param a The first iterator.
     * @param b The second iterator.
     * @return True if the iterators point to different columns; false otherwise.
     */
    friend bool operator!= ( DeepFirstIterator const & a, DeepFirstIterator const & b )
    { return a.m_currentColumn != b.m_currentColumn; };

    /**
     * @brief Gets the current layer (depth) of the iterator.
     * @return The current layer (depth) of the iterator.
     */
    size_t getCurrentLayer() const
    { return m_currentLayer; }

private:
    /// Pointer to the current column
    ColumnType * m_currentColumn;
    /// The current depth of the iterator
    size_t m_currentLayer;
  };

  /**
   * @return Return an itarator pointing on the first leaf of the first columns vector
   * Example on 2 column with Column A : 2 layer and Column B : 3 layers
   * A.A -> A-B -> A-C -> A -> B-A-A -> B-A-B -> B-A -> B-B-A -> B-B-B -> B-B -> B
   */
  DeepFirstIterator beginDeepFirst() const;

  /**
   * @return Return a end itarator
   * This iterator is initialized with a null pointer
   * representing a position after the last valid element
   */
  DeepFirstIterator endDeepFirst() const
  { return DeepFirstIterator( nullptr, 0 ); }

  /// Alias for an initializer list of variants that can contain either a string or a layout column.
  using TableLayoutArgs = std::initializer_list< std::variant< string_view, TableLayout::Column > >;

  /// Alias for the list of columns.
  using ColumnsList = Column::ColumnsList;


  TableLayout()
  {
    setMargin( MarginValue::medium );
  }

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param columns A vector containing all column initialized
   */
  TableLayout( string_view title,
               stdVector< TableLayout::Column > const & columns )
  {
    setMargin( MarginValue::medium );
    setTitle( title );
    addColumns( columns );
  }

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param args An initializer_list containing string / column
   */
  TableLayout( string_view title,
               TableLayoutArgs args )
  {
    setMargin( MarginValue::medium );
    setTitle( title );
    processArguments( args );
  }

  /**
   * @brief Construct a new Table Layout object
   * @param args An initializer_list containing string / column
   */

  TableLayout( TableLayoutArgs args )
  {
    setMargin( MarginValue::medium );
    processArguments( args );
  }

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param args An initializer_list containing string / column
   */
  TableLayout( string_view title,
               stdVector< string > const & args )
  {
    setMargin( MarginValue::medium );
    setTitle( title );
    addColumns( args );
  }

  /**
   * @return The columns list
   */
  ColumnsList const & getColumns() const
  { return m_tableColumns; }

  /**
   * @return The columns list
   */
  ColumnsList & getColumns()
  { return m_tableColumns; }

  /**
   * @return The table name. Returned as a for multiline support.
   */
  CellLayout const & getTitleLayout() const
  { return m_tableTitleLayout; }

  /**
   * @return The table name. Returned as a for multiline support.
   */
  CellLayout & getTitleLayout()
  { return m_tableTitleLayout; }

  /**
   * @return The table name. Returned as a for multiline support.
   */
  string_view getTitleStr() const
  { return m_tableTitleStr; }

  /**
   * @param title The table title
   * @return The tableLayout reference
   */
  TableLayout & setTitle( string_view title );

  /**
   * @brief Remove the return line at the end & begenning of the table
   * @param value Value to desactivate or not wrapLine at the end
   * @return The tableLayout reference
   */
  TableLayout & enableLineBreak( bool value );

  /**
   * @brief Set the minimal margin width between cell content and borders.
   * @param marginValue The margin value
   * @return The tableLayout reference
   */
  TableLayout & setMargin( MarginValue marginValue );

  /**
   * @brief Set the maximal width for each column
   * @param width The max column width
   * @return The tableLayout reference
   */
  TableLayout & setMaxColumnWidth( size_t width );

  /**
   * @brief check if a column max width has been set
   * @return Truef a column max width has been set, otherwise false
   */
  bool isMaxColumnWidthSet()
  { return m_maxColumnWidth != noColumnMaxWidth; }

  /**
   * @return check if the line break at the end & beginning is activated
   */
  bool isLineBreakEnabled() const;

  /**
   * @return The number of spaces at each table sides
   */
  integer const & getBorderMargin() const
  { return m_borderMargin; }

  /**
   * @return The number of character between two columns (spaces + the separacting character).
   */
  integer const & getColumnMargin() const
  { return m_columnMargin; }

  /**
   * @return The number of margin spaces around contents.
   */
  integer const & getMarginValue() const
  { return m_marginValue; }

  /**
   * @return The margin title
   */
  size_t const & getMaxColumnWidth() const
  { return m_maxColumnWidth; }

  /**
   * @brief Create and add columns to the columns vector given a string vector
   * @param columnNames The columns name
   */
  void addColumns( stdVector< TableLayout::Column > const & columnNames );

  /**
   * @brief Create and add columns to the columns vector given a string vector
   * @param columns The columns list
   */
  void addColumns( stdVector< string > const & columns );

  /**
   * @brief Create and add a column to the columns vector given a string
   * @param columnName The column name
   */
  void addColumn( string_view columnName );

  /**
   * @brief Create and add a column to the columns vector given a Column
   * @param column Vector containing addition information on the column
   */
  void addColumn( TableLayout::Column const & column );

protected:

  /**
   * @brief Add a column to the table given an initializer_list of string & Column
   * @param args An initializer_list containing string / column
   */
  void processArguments( TableLayoutArgs args )
  {
    for( auto const & arg : args )
    {
      std::visit( [this]( auto const & value ) {
        addColumn( value );
      }, arg );
    }
  }

  /**
   * @tparam Ts The remaining arguments
   * @param args The remaining arguments to be processed
   */
  template< typename ... Ts >
  void processArguments( Ts &... args )
  {
    addColumns( args ... );
  }

  /// Columns settings hierarchy
  ColumnsList m_tableColumns;

  /// Indicate if we have a line break a the beginning of the table
  bool m_lineBreakAtBegin = true;

  /// Table title text
  string m_tableTitleStr;

  /// Table title cell layout settings
  CellLayout m_tableTitleLayout = CellLayout( CellType::Header, Alignment::center );

  /// Max width for each column
  size_t m_maxColumnWidth = noColumnMaxWidth;


  /// The number of spaces at each table sides
  integer m_borderMargin;

  /// The number of character between two columns (spaces + the separacting character).
  integer m_columnMargin;

  /// The number of margin spaces around contents.
  integer m_marginValue;

};

/**
 * @brief Variation of the TableLayout to store precomputed layout information, ready to be formatted.
 */
class PreparedTableLayout : public TableLayout
{
public:

  /**
   * @brief Construct a default Table Formatter without layout specification (to only insert data in it,
   * without any column / title). Feature is not tested.
   */
  PreparedTableLayout();

  /**
   * @brief Precompute various information for formatting from a configurated TableLayout:
   *        - parent-child relationships between columns and sub-columns,
   *        - layout elements size,
   *        - line wrapping.
   *        For now, called automatically at TableFormatter construction.
   * @note If an error happen while this process, it must output the table name and the error
   *       message in a GEOS_WARNING().
   * @param other The table layout configuration.
   */
  PreparedTableLayout( TableLayout const & other );

  /**
   * @brief As prepared CellLayout & Column types have internal pointers, we cannot copy this class.
   */
  PreparedTableLayout( PreparedTableLayout const & ) = delete;

  /**
   * @brief as prepared CellLayout & Column types have internal pointers, we cannot move this class
   *        (SSO breaks string<-string_view move).
   */
  PreparedTableLayout( PreparedTableLayout && ) = delete;

  /**
   * @return The count of column layers
   */
  size_t getColumnLayersCount() const
  { return m_columnLayersCount; }

  /**
   * @return The number of visible columns that does not contain child (useful to know the maximum number of
   *         column to show in a given row).
   */
  size_t getLowermostColumnsCount() const
  { return m_lowermostColumnCount; }

private:

  size_t m_columnLayersCount;
  size_t m_lowermostColumnCount;

  /**
   * @brief Recursive part of column layout preparation, see constructor documentation.
   * @param columns The list of columns to prepare.
   */
  void prepareLayoutRecusive( stdVector< TableLayout::Column > & columns, size_t level );

};

}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */

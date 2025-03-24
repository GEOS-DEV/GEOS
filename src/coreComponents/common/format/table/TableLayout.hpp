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
    Alignment headerAlignment = Alignment::center;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment = Alignment::right;
  };

  /**
   * @struct CellLayout
   * @brief View on cell data grouping the cell information to display it in a table (content, type, alignment, ...).
   * @note the source text must not be freeed/moved as the CellLayout will *only be a view on the text data*.
   */
  struct CellLayout
  {
    /// Maximum length of the data in the cell.
    size_t m_cellWidth;
    /// The type of the cell (Header,Value, Merge, ...).
    CellType m_cellType;
    /// The alignment of the cell (left, center, right).
    Alignment m_alignment;
    /// vector containing each cell content, separated by lines.
    std::vector< string_view > m_lines;

    /**
     * @brief Constructor to initialize a Cell with a default settings.
     */
    CellLayout();

    /**
     * @brief Constructor to initialize an empty Cell of a given type.
     * @param cellType The type of the cell.
     */
    CellLayout( CellType cellType );

    /**
     * @brief Constructor to fully initialize a cell given with given celltype, text and alignment.
     * m_cellWidth will be initialized aDter
     * @param cellType The type of the cell.
     * @param alignment The alignment of the cell (left, right, or center).
     */
    CellLayout( CellType cellType, TableLayout::Alignment alignment );

    /**
     * @param inputText The view on the text of the cell. `m_lines` will contain each separated lines, and
     *                  `m_cellWidth`, the maximum line width. Called automatically by PreparedTableLayout.
     * @param maxLineWidth The maximum allowed line width. Use `noColumnMaxWidth` to disable.
     */
    void prepareLayout( string_view value, size_t maxLineWidth );
  };

  /**
   * @class Column
   * @brief Class representing a column in a table layout.
   */
  class Column
  {
public:
    /// Alias for the list of columns.
    using ColumnsList = std::vector< Column >;

    // The text of the header.
    string m_headerStr;
    /// The header cell layout (view on m_headerStr).
    CellLayout m_headerLayout;
    /// A vector containing all sub-columns in the column.
    ColumnsList m_subColumns;
    /// struct containing m_alignment for the column (header and values)
    ColumnAlignement m_alignment;

    /**
     * @brief Default constructor.
     * Initializes a column with default values.
     */
    Column();

    /**
     * @brief Move constructor. Ignore any input pointer (m_next, m_parent).
     */
    explicit Column( string_view name ):
      Column( name, ColumnAlignement() )
    {}

    /**
     * @brief Move constructor. Ignore any input pointer (m_next, m_parent).
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
    Column & addSubColumns( std::vector< string > const & subColNames );

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
    { return m_headerLayout.m_cellType!=CellType::Hidden; }

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
               std::vector< TableLayout::Column > const & columns )
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
               std::vector< string > const & args )
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
   * @return The border margin,
   * number of spaces at each table sides
   */
  integer const & getBorderMargin() const
  { return m_borderMargin; }

  /**
   * @return The column margin,
   * numbers of spaces separating both left and right side from a vertical line
   */
  integer const & getColumnMargin() const
  { return m_columnMargin; }

  /**
   * @return The table margin value
   */
  integer const & getMarginValue() const
  { return m_marginValue; }

  /**
   * @return The margin title
   */
  integer const & getMarginTitle() const
  { return m_titleMargin; }

  /**
   * @return The margin title
   */
  size_t const & getMaxColumnWidth() const
  { return m_maxColumnWidth; }

  /**
   * @brief Create and add columns to the columns vector given a string vector
   * @param columnNames The columns name
   */
  void addColumns( std::vector< TableLayout::Column > const & columnNames );

  /**
   * @brief Create and add columns to the columns vector given a string vector
   * @param columns The columns list
   */
  void addColumns( std::vector< string > const & columns );

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

  // Indicate if we have a line break a the beginning of the table
  bool m_lineBreakAtBegin = true;

  /// Table title text
  string m_tableTitleStr;

  /// Table title cell layout settings
  CellLayout m_tableTitleLayout = CellLayout( CellType::Header, Alignment::center );

  /// Max width for each column
  size_t m_maxColumnWidth = noColumnMaxWidth;


  integer m_borderMargin;
  integer m_columnMargin;
  integer m_marginValue;
  integer m_titleMargin = 2;

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
  void prepareLayoutRecusive( std::vector< TableLayout::Column > & columns, size_t level );

};

}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */

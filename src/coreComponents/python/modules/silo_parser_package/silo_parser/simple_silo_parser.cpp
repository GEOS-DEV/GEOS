


#include <Python.h>
#include <silo.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <map>


static PyObject * parse_file(PyObject *self, PyObject *args, PyObject *keywords)
{
  //---------------------------------------------------------------------
  // Parse inputs
  //---------------------------------------------------------------------
  std::string ghostRankString("ghostRank");

  // Command line arguments
  const char *file_name = "plot_000000";
  const char *field_type = "Fracture_ElementFields";
  PyObject *field_names;
  const char *parallel_folder="";

  static char *keywordlist[] = {"file_name", "field_type", "field_names", "parallel_folder", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywords, "ssOs", keywordlist, &file_name,
                                   &field_type, &field_names, &parallel_folder))
  {
    std::cout << "Error: inputs were not specified correctly!  Correct argments:" << std::endl;
    std::cout << "parse_file(file_name, field_type, field_names)" << std::endl;
    return NULL;
  }

  // Check the field_names input
  std::vector<std::string> field_names_requested;
  if (PyList_Check(field_names))
  {
    size_t field_names_size = PyList_Size(field_names);
    if (field_names_size > 0)
    {
      field_names_requested.resize(field_names_size);
      for (size_t ii=0; ii<field_names_size; ++ii)
      {
        PyObject *list_val = PyList_GetItem(field_names, ii);
        if (PyUnicode_Check(list_val))
        {
          field_names_requested[ii] = PyUnicode_AsUTF8(list_val);
        }
        else
        {
          std::cout << "Error: field_names[" << ii << "] is not a string!" << std::endl;
        }
      }
    }
    else
    {
      std::cout << "Error: field_names is empty" << std::endl;
    }
  }
  else
  {
    std::cout << "Error: field_names must be a list containing the field names to be exported!" << std::endl;
  }

  if (field_names_requested.size() == 0)
  {
    return NULL;
  }

  // Add location string based that matches the field type
  if (strcmp(field_type, "Fracture_ElementFields") == 0)
  {
    field_names_requested.push_back("elementCenter");
  }
  else
  {
    field_names_requested.push_back("ReferencePosition");
  }
  

  //---------------------------------------------------------------------
  // Check for available field values, find files to pars
  //---------------------------------------------------------------------
  std::vector<std::string> field_names_available;
  std::map< std::string, std::vector<std::string> > active_map;

  // Root of file
  int db_result = 0;
  DBfile *db = DBOpen(file_name, DB_UNKNOWN, DB_READ);
  int current_cycle = *((int *)DBGetVar(db, "cycle"));
  double current_time = *((double *)DBGetVar(db, "dtime"));

  // Field Dir
  db_result = DBSetDir(db, field_type);
  if (db_result == 0)
  {
    // Check the available keys
    for (size_t ii=0; ii<field_names_requested.size(); ++ii)
    {
      // Grab the field
      DBmultivar* var = DBGetMultivar(db, field_names_requested[ii].c_str());
      if (var != NULL)
      {
        field_names_available.push_back(field_names_requested[ii]);
      }
    }

    // Use ghostRank to check data locations
    DBmultivar* var = DBGetMultivar(db, ghostRankString.c_str());
    if (var != NULL)
    {
      for (int jj=0; jj<var->nvars; ++jj)
      {
        if (strcmp(var->varnames[jj], "EMPTY") != 0)
        {
          // Split the var name into file / domain
          std::stringstream ss_a(var->varnames[jj]);
          std::string segment;
          std::vector<std::string> varsplit_a;
          while(std::getline(ss_a, segment, ':'))
          {
             varsplit_a.push_back(segment);
          }
          
          std::stringstream ss_b(varsplit_a[1]);
          std::vector<std::string> varsplit_b;
          while(std::getline(ss_b, segment, '/'))
          {
             varsplit_b.push_back(segment);
          }
          
          // Add values to map
          std::map< std::string, std::vector<std::string> >::iterator it = active_map.find(varsplit_a[0]);
          if ( it == active_map.end())
          {
            std::vector<std::string> new_domain_list;
            new_domain_list.push_back(varsplit_b[1]);
            active_map.insert(std::pair< std::string, std::vector<std::string> >(varsplit_a[0], new_domain_list));
          }
          else
          {
            it->second.push_back(varsplit_b[1]);
          }
        }
      }
    }
  }
  DBClose(db);

  


  //---------------------------------------------------------------------
  // Extract Data
  //---------------------------------------------------------------------
  std::map< std::string, std::vector< std::vector< int > > > field_values_ints;
  std::map< std::string, std::vector< std::vector< double > > > field_values_floats;
  int init_maps = 1;
  int max_size = 0;
  int current_size = 0;

  // Walk through all of the files, calculate max size
  std::map< std::string, std::vector<std::string> >::iterator it;
  for (it = active_map.begin(); it != active_map.end(); ++it)
  {
    db = DBOpen(it->first.c_str(), DB_UNKNOWN, DB_READ);

    for (int ii=0; ii<it->second.size(); ++ii)
    {
      // Open domain
      db_result = DBSetDir(db, it->second[ii].c_str());

      if (db_result == 0)
      {
        // Open field type
        db_result = DBSetDir(db, field_type);

        if (db_result == 0)
        {
          // Count the number of elements
          DBucdvar *var_ghost = DBGetUcdvar(db, ghostRankString.c_str());
          if (var_ghost != NULL)
          {
            max_size += var_ghost->nels;
          }

          // Go back to the domain
          db_result = DBSetDir(db, "..");
        }
        else
        {
          std::cout << "Could not open field directory!" << std::endl;
        }

        // Go back to the root
        db_result = DBSetDir(db, "..");
      }
      else
      {
        std::cout << "Could not open domain directory!" << std::endl;
      }
    }
    DBClose(db);
  }


  // Iterate over the file map again, extracting data 
  for (it = active_map.begin(); it != active_map.end(); ++it)
  {
    // std::cout << "Opening file: " << it->first << std::endl;
    db = DBOpen(it->first.c_str(), DB_UNKNOWN, DB_READ);

    for (int ii=0; ii<it->second.size(); ++ii)
    {
      // Open domain
      // std::cout << "  " << it->second[ii] << std::endl;
      db_result = DBSetDir(db, it->second[ii].c_str());

      if (db_result == 0)
      {
        // Open field type
        db_result = DBSetDir(db, field_type);

        if (db_result == 0)
        {
          // Setup the value maps
          if (init_maps == 1)
          {
            // std::cout << "Map setup" << std::endl;

            for (int jj=0; jj<field_names_available.size(); jj++)
            {
              DBucdvar *var_test = DBGetUcdvar(db, field_names_available[jj].c_str());
              if (var_test != NULL)
              {
                int ndim = var_test->nvals;
                if (var_test->datatype == DB_DOUBLE)
                {
                  std::vector< double > va(max_size, 0.0);
                  std::vector< std::vector< double > > vb(ndim, va);
                  field_values_floats.insert(std::pair< std::string, std::vector< std::vector< double > > >(field_names_available[jj], vb));
                }
                else if (var_test->datatype == DB_INT)
                {
                  std::vector< int > va(max_size, 0);
                  std::vector< std::vector< int > > vb(ndim, va);
                  field_values_ints.insert(std::pair< std::string, std::vector< std::vector< int > > >(field_names_available[jj], vb));
                }
              }
            }

            init_maps = 0;
          }

          // Build the list of valid elements from the ghostRank
          std::vector< int > ghost_map;
          DBucdvar *var_ghost = DBGetUcdvar(db, ghostRankString.c_str());
          if (var_ghost != NULL)
          {
            // std::cout << "Ghost rank setup" << std::endl;
            int *var_ghost_sub = (int *)(var_ghost->vals[0]);
            for (int jj=0; jj<var_ghost->nels; ++jj)
            {
              if (var_ghost_sub[jj] < 0)
              {
                ghost_map.push_back(jj);
              }
            }
          }
          int sub_size = ghost_map.size();
          
          // Extract float data
          // std::cout << "Extracting float data" << std::endl;
          std::map< std::string, std::vector< std::vector< double > > >::iterator it_double;
          for (it_double = field_values_floats.begin(); it_double != field_values_floats.end(); ++it_double)
          {
            // std::cout << "  " << it_double->first << std::endl;
            DBucdvar *var_double = DBGetUcdvar(db, it_double->first.c_str());
            for (int jj=0; jj<it_double->second.size(); ++jj)
            {
              double *var_double_sub = (double *)(var_double->vals[jj]);
              int old_size = it_double->second[jj].size();
              it_double->second[jj].resize(old_size + sub_size);
              
              for (int kk=0; kk<sub_size; ++kk)
              {
                it_double->second[jj][kk+current_size] = var_double_sub[ghost_map[kk]];
              }
            }
          }

          // Extract int data
          std::map< std::string, std::vector< std::vector< int > > >::iterator it_int;
          for (it_int = field_values_ints.begin(); it_int != field_values_ints.end(); ++it_int)
          {
            // std::cout << it_int->first << std::endl;
            DBucdvar *var_int = DBGetUcdvar(db, it_int->first.c_str());
            for (int jj=0; jj<it_int->second.size(); ++jj)
            {
              int *var_int_sub = (int *)(var_int->vals[jj]);
              int old_size = it_int->second[jj].size();
              it_int->second[jj].resize(old_size + sub_size);

              for (int kk=0; kk<sub_size; ++kk)
              {
                it_int->second[jj][kk+current_size] = var_int_sub[ghost_map[kk]];
              }
            }
          }

          // Go back to the domain
          current_size += sub_size;
          db_result = DBSetDir(db, "..");
        }

        // Go back to the root
        db_result = DBSetDir(db, "..");
      }
    }

    DBClose(db);
  }


  //---------------------------------------------------------------------
  // Build Outputs
  //---------------------------------------------------------------------
  // Setup the dictionary that will be returned
  PyObject *plot_dict = PyDict_New();
  PyDict_SetItem(plot_dict, PyUnicode_FromString("time"), PyFloat_FromDouble(current_time));
  PyDict_SetItem(plot_dict, PyUnicode_FromString("cycle"), PyLong_FromLong(current_cycle));

  // Write float data
  std::map< std::string, std::vector< std::vector< double > > >::iterator it_double;
  for (it_double = field_values_floats.begin(); it_double != field_values_floats.end(); ++it_double)
  {
    // Construct as list of lists
    PyObject *current_list = PyList_New(it_double->second.size());
    for(size_t ii=0; ii<it_double->second.size(); ++ii)
    {
      PyObject *sub_list = PyList_New(it_double->second[ii].size());
      for (size_t jj=0; jj<it_double->second[ii].size(); ++jj)
      {
        PyList_SetItem(sub_list, jj, PyFloat_FromDouble(it_double->second[ii][jj]));
      }

      PyList_SetItem(current_list, ii, sub_list);
    }      

    // Add to dict
    PyDict_SetItem(plot_dict, PyUnicode_FromString(it_double->first.c_str()), current_list);
  }

  // Write int data
  std::map< std::string, std::vector< std::vector< int > > >::iterator it_ints;
  for (it_ints = field_values_ints.begin(); it_ints != field_values_ints.end(); ++it_ints)
  {
    // Construct as list of lists
    PyObject *current_list = PyList_New(it_ints->second.size());
    for(size_t ii=0; ii<it_ints->second.size(); ++ii)
    {
      PyObject *sub_list = PyList_New(it_ints->second[ii].size());
      for (size_t jj=0; jj<it_ints->second[ii].size(); ++jj)
      {
        PyList_SetItem(sub_list, jj, PyLong_FromLong((long)it_ints->second[ii][jj]));
      }
      PyList_SetItem(current_list, ii, sub_list);
    }      

    // Add to the dict
    PyDict_SetItem(plot_dict, PyUnicode_FromString(it_ints->first.c_str()), current_list);
  }

  return plot_dict;
}


// Method definition
static PyMethodDef mod_methods[] = {{"parse_file", (PyCFunction)parse_file, METH_VARARGS|METH_KEYWORDS, "Method to parse a single silo file"},
                                {NULL, NULL, 0, NULL}};

// Module definition
static struct PyModuleDef mod_definition = { 
    PyModuleDef_HEAD_INIT,
    "simple_silo_parser",
    "Simple SILO parser",
    -1,
    mod_methods
};


PyMODINIT_FUNC PyInit_simple_silo_parser(void)
{
  return PyModule_Create(&mod_definition);
}








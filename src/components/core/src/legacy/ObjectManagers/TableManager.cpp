/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TableManager.cpp
 * @author walsh24
 * @date Sep 12, 2012
 */

#include "TableManager.h"
#include "mpi.h"

void TableManager::ReadXML(TICPP::HierarchicalDataNode* TablesNode)
{
  int initialized;
  int numtasks, rank;

  MPI_Initialized(&initialized);
  if (!initialized)
    MPI_Init(NULL, NULL);
//  MPI_Status Stat;
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (TICPP::HierarchicalDataNode* TableNode = TablesNode->Next(true) ; TableNode ;
       TableNode = TablesNode->Next())
  {
    const std::string label(TableNode->Heading());

    unsigned int dim = 0;
    localIndex nComponents = 0;
    if(streq("Table1D",label))
    {
      dim = 1;
      nComponents = 1;
    }
    else if(streq("Table2D",label))
    {
      dim = 2;
      nComponents = 1;
    }
    else if(streq("Table3D",label))
    {
      dim = 3;
      nComponents = 1;
    }
    else if(streq("Table4D",label))
    {
      dim = 4;
      nComponents = 1;
    }
    else if(streq("VectorField1D",label))
    {
      dim = 1;
      nComponents = 3;
    }
    else if(streq("VectorField2D",label))
    {
      dim = 2;
      nComponents = 3;
    }
    else if(streq("VectorField3D",label))
    {
      dim = 3;
      nComponents = 3;
    }
    else if(streq("VectorField4D",label))
    {
      dim = 4;
      nComponents = 3;
    }
    else
      throw GPException("TableManager::ReadXML - Table heading not recognized");

    //create table
    const std::string tableName = TableNode->GetAttributeString("name");
    if(tableName.empty())
      throw GPException("Table specified without a name");

    //read grid data
    array<array<real64> > x(dim);
    array<array<real64> > bufX(dim);
    gArray1d bufXLength(dim);

    if (rank == 0)
    {
      for(localIndex i = 0 ; i < dim ; i++)
      {
        //read values
        {
          const std::string fname = i == 0 ? "x_file" : (i == 1 ? "y_file" : (i == 2 ? "z_file" : "t_file") );
          std::string ticksFile = TableNode->GetAttributeStringOrDefault(fname, "");
          if (!ticksFile.empty())
            dlmreadVector(ticksFile, x(i), ' ');
          else
          {
            if(dim == 1)
              x(i) = TableNode->GetAttributeVector<realT>("coord", ",");
            else if( dim==3 )
            {
              if( i==0 )
                x(i) = TableNode->GetAttributeVector<realT>("xcoord", ",");
              else if( i==1 )
                x(i) = TableNode->GetAttributeVector<realT>("ycoord", ",");
              else if( i==2 )
                x(i) = TableNode->GetAttributeVector<realT>("zcoord", ",");
            }
          }
        }

        //read units
        {
          const std::string name = i == 0 ? "x_units" : (i == 1 ? "y_units" : (i == 2 ? "z_units" : "t_units") );
          const realT units = TableNode->GetAttributeOrDefault<realT>(name, 1.0);
          if (!isEqual(units, 1.0, 0.0))
            x(i) *= units;
        }
        bufXLength[i] = x(i).size();
      }
    }

    // Optional offset for time tables
    if(dim == 1)
    {
      realT offsetTime = TableNode->GetAttributeOrDefault<realT>("offset_time", 0.0);
      x(0) += offsetTime;
    }

    //Broadcast
    {
      //bufXLength
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(bufXLength.data(), dim, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    }

    //Update
    for (localIndex i = 0 ; i < dim ; i++)
      bufX[i].resize(bufXLength[i]);

    //Pack up
    if (rank == 0)
    {
      for (localIndex i = 0 ; i < dim ; i++)
        for (globalIndex j = 0 ; j < bufX[i].size() ; j++)
          bufX[i][j] = x(i)[j];
    }

    //Broadcast
    {
      //bufX
      MPI_Barrier(MPI_COMM_WORLD);
      for (localIndex i = 0 ; i < dim ; i++)
        MPI_Bcast(bufX[i].data(), bufX[i].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    //Update
    if (rank != 0)
    {
      for (localIndex i = 0 ; i < dim ; i++)
        for (globalIndex j = 0 ; j < bufX[i].size() ; j++)
          x(i).push_back(bufX[i][j]);
    }

    //read value data
    array<real64> values;
    array<real64> bufValues;
    globalIndex bufValuesLength = 0;

    if (rank == 0)
    {
      std::cout << tableName << std::endl;

      //read values
      {
        std::string ticksFile = TableNode->GetAttributeStringOrDefault("value_file", "");
        if (!ticksFile.empty())
          dlmreadVector(ticksFile, values, ' ');
        else// if(dim == 1)
        {
          values = TableNode->GetAttributeVector<realT>("value", ",");
        }

        if(values.size()==0 && ticksFile.empty())
        {
          if(dim == 3)
          {
            ticksFile = TableNode->GetAttributeStringOrDefault("voxel_file", "");
            if (!ticksFile.empty())
              ReadVoxelFile(ticksFile, nComponents, values);
          }
          else if(dim == 4)
          {
            ticksFile = TableNode->GetAttributeStringOrDefault("time_voxel_file", "");
            if (!ticksFile.empty())
              ReadTimeVoxelFile(ticksFile, nComponents, values);
            else if(x[0].size() == 0)
            {
              ticksFile = TableNode->GetAttributeStringOrDefault("nuft_voxel_file", "");
              if (!ticksFile.empty())
              {
                ReadNUFTFile(ticksFile, x, values);
                nComponents = 1;
              }
            }
          }
        }
      }

      //read units and apply
      {
        const realT units = TableNode->GetAttributeOrDefault<realT>("value_units", 1.0);
        if (!isEqual(units, 1.0, 0.0))
          values *= units;
      }

      bufValuesLength = values.size();
    }

    //Broadcast
    {
      //bufXLength
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&bufValuesLength, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    }

    //Update
    bufValues.resize(bufValuesLength);

    //Pack up
    if (rank == 0)
    {
      for (globalIndex i = 0 ; i < bufValues.size() ; i++)
        bufValues[i] = values[i];
    }

    //Broadcast
    {
      //bufValues
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(bufValues.data(), bufValues.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    //Update
    if (rank != 0)
    {
      for (globalIndex i = 0 ; i < bufValues.size() ; i++)
        values.push_back(bufValues[i]);
    }

    //check that values were set
    if(values.size()==0)
    {
      throw GPException("TableManager::ReadXML : failed to set values");
    }

    // set interpolation order
    TableInterpolation::Order interpOrder;
    {
      const int interp = TableNode->GetAttributeOrDefault<int>("interpolation", 1); // linear
                                                                                    // interpolation
                                                                                    // by
                                                                                    // default
      switch(interp)
      {
      case 0:
        interpOrder = TableInterpolation::zeroth;
        break;
      case 1:
        interpOrder = TableInterpolation::linear;
        break;
      default:
        throw GPException("TableManager::ReadXML : only 0th and 1st order table interpolation currently supported.");

      }
    }



    if(nComponents == 1)
    {
      switch(dim)
      {
      case 1:
        NewTable<1>(tableName, x, values,interpOrder);
        break;
      case 2:
        NewTable<2>(tableName, x, values,interpOrder);
        break;
      case 3:
        NewTable<3>(tableName, x, values,interpOrder);
        break;
      case 4:
        NewTable<4>(tableName, x, values,interpOrder);
        break;
      }
    }
    else if(nComponents == 3)
    {
      array<R1Tensor> values1(values.size() / 3);
      {
        localIndex i= 0, j = 0;
        for(array<real64>::const_iterator it = values.begin() ; it != values.end() ; ++it)
        {
          values1[i][j++] = *it;
          if(j==4)
          {
            ++i;
            j = 0;
          }
        }
      }
      switch(dim)
      {
      case 1:
        NewVectorField<1>(tableName, x, values1);
        break;
      case 2:
        NewVectorField<2>(tableName, x, values1);
        break;
      case 3:
        NewVectorField<3>(tableName, x, values1);
        break;
      case 4:
        NewVectorField<4>(tableName, x, values1);
        break;
      }
    }
    else
    {
      throw GPException("TableManager::ReadXML : can only support data structures that are rank 0 or 1 tensors");
    }
  } //end for
}

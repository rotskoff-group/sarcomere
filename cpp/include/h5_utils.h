#ifndef H5_UTILS_H
#define H5_UTILS_H

#include <iostream>
#include <H5Cpp.h>
#include <vector>
#include <string>
#include <cstdio>       // For std::remove
#include "components.h" // This header must define Filament, Myosin, and utils::MoleculeConnection

// For convenience, use the vec type from the utils namespace.
using vec = utils::vec;

//---------------------------------------------------------------------
// Function Prototypes
//---------------------------------------------------------------------

// Serialize and pad actinIndicesPerActin.
// It returns a 2D vector where each row corresponds to the connections of an actin,
// padded to max_bonds (using -1 for missing bonds).
std::vector<std::vector<int>> serializeActinIndicesPerActin(
    utils::MoleculeConnection& actinIndicesPerActin, int n_actins, int& max_bonds);

// Flatten a vector of 2D points (vec) into a 1D array of doubles.
std::vector<double> flatten_2d_array(std::vector<vec> array);

// Create an empty (extendable) dataset within the specified group of the file.
void create_empty_dataset(H5::H5File& file, const std::string& groupName,
                          const std::string& datasetName,
                          const std::vector<hsize_t>& initialDims,
                          const std::vector<hsize_t>& maxDims,
                          const std::vector<hsize_t>& chunkDims);

// Append new data to an existing dataset within a given group.
void append_to_dataset(H5::Group& group, const std::string& datasetName,
                       const std::vector<double> newData,
                       const std::vector<hsize_t>& newDims);

// Create a new HDF5 file and set up datasets for actin and myosin.
void create_file(std::string& filename, Filament& actin, Myosin& myosin);

// Append simulation data (actin, myosin, connection indices) to the file.
void append_to_file(std::string& filename, Filament& actin, Myosin& myosin, 
    std::vector<double>& flatActinBonds, 
    std::vector<double>& flatMyosinBonds);

// Load data from a dataset into a 1D vector of doubles.
// The dataset dimensions are returned in the dims vector.
std::vector<double> load_from_dataset(H5::Group& group, const std::string& datasetName,
                                      std::vector<hsize_t>& dims);

// Load actin and myosin data (and energy, if needed) from the file.
void load_from_file(std::string& filename, Filament& actin, Myosin& myosin, 
                    std::vector<std::vector<int>>& actin_actin_bonds, int& n_frames);

#endif // H5_UTILS_H

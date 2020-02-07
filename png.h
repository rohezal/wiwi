#ifndef WW_PNG_HEADER
#define WW_PNG_HEADER

#include <iostream>
#include <omp.h>
#include <png++/png.hpp>
#include <tiff.h>
#include <tiffio.h>
#include <string>

//load the tiff into a vector
std::vector<std::vector<float> > loadTiff(const std::string& _name);

//converts loadtiff  vector into an array
std::pair<size_t, size_t> loadTiffArray(const std::string& _name, float **array);

//saves a float array into a file in the tiff format
void saveTiffArray(const std::string& _name, float* array, size_t width, size_t height);
#endif

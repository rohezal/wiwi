#ifndef WW_PNG_HEADER
#define WW_PNG_HEADER

#include <iostream>
#include <omp.h>
#include <png++/png.hpp>
#include <tiff.h>
#include <tiffio.h>
#include <string>

std::vector<std::vector<float> > loadTiff(const std::string& _name);
std::pair<size_t, size_t> loadTiffArray(const std::string& _name, float **array);
void saveTiffArray(const std::string& _name, float* array, size_t width, size_t height);
#endif

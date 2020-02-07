#include <iostream>
#include <omp.h>
#include <png++/png.hpp>
#include "png.h"
//#include <cmath>

//how many pixels are inside of a circle in a certain
int numberOfPixelsInCirlce(int border, int radius_squarded)
{
    int border_half = border/2;
    int counter = 0;

    for(int a = -border_half; a <= border_half; a++)
    {
        for(int b = -border_half; b <= border_half; b++)
        {
            counter += (a*a+b*b) < radius_squarded;
        }
    }

    return counter;
}

float getValuesInsideCircle(float* array, size_t x, size_t y, int size_x, int size_y, int _radius_squared, int _border_half)
{
    const int border_half = _border_half;
    int counter = 0;
    float result = 0;
    const float radius_squarded = _radius_squared;

    for(int a = -border_half; a <= border_half; a++)
    {
        for(int b = -border_half; b <= border_half; b++)
        {
            const int current_position_y = y+a;
            const int current_position_x = x+b;
            bool valid_x = (current_position_x >= 0) && (current_position_x < size_x);
            bool valid_y = (current_position_y >= 0) && (current_position_y < size_y);

            if(valid_x &&  valid_y && (a*a+b*b) < radius_squarded)
            {
                counter++;
                result += array[current_position_y*size_x+current_position_x];
                //result = 30;
            }
        }
    }

    if(counter > 0)
    {
        //std::cout << "Counter = " << counter << " | result: " << result;
        return result / counter;
    }
    else
    {
        return 0;
    }
}

using namespace std;

int main()
{

    float *tifdata;
    std::pair<size_t, size_t> tifsize = loadTiffArray("brightness_win1.tif",&tifdata); //load the tiff data into an array
    const int tsize_x = tifsize.second;
    const int tsize_y = tifsize.first;

    const int size_x = tsize_x;
    const int size_y = tsize_y;

    double start_time;

    std::cout << "Width: " << tsize_x << " | Height: " << tsize_y << std::endl;

    float* image = tifdata;
    float* new_image = new float[size_x*size_y];
    float* new_image2 = new float[size_x*size_y];
    //float* new_image2 = combined_array;



    const int border = 16;
    const int border_half = border/2;
    const int border_squared = border*border;
    const float circle_radius_squared = (border/2)*(border/2);
    const int number_of_pixels_in_circle = numberOfPixelsInCirlce(border,circle_radius_squared); //calculate this
    std::cout <<  "number_of_pixels_in_circle: " << number_of_pixels_in_circle << std::endl;


    #pragma omp parallel for
    for(int y = 0; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            new_image[y*size_x+x] = 0;
            new_image2[y*size_x+x] = 0;
        }
    }

    //border regions
    start_time = omp_get_wtime();

    #pragma omp parallel for
    for(int y = 0; y < border_half; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            new_image[y*size_x+x] = getValuesInsideCircle(image, x, y, size_x, size_y, circle_radius_squared, border_half);
            //new_image[y*size_x+x] = 20;
        }
    }

    #pragma omp parallel for
    for(int y = size_y-border_half; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            new_image[y*size_x+x] = getValuesInsideCircle(image, x, y, size_x, size_y, circle_radius_squared, border_half);
        }
    }

    #pragma omp parallel for
    for(int y = border_half; y < size_y-border_half; y++)
    {
        for(int x = 0; x < border_half; x++)
        {
            new_image[y*size_x+x] = getValuesInsideCircle(image, x, y, size_x, size_y, circle_radius_squared, border_half);
        }
    }

    #pragma omp parallel for
    for(int y = border_half; y < size_y-border_half; y++)
    {
        for(int x = size_x-border_half; x < size_x; x++)
        {
            new_image[y*size_x+x] = getValuesInsideCircle(image, x, y, size_x, size_y, circle_radius_squared, border_half);
        }
    }

    std::cout << "Runtime borders: " << omp_get_wtime()-start_time << std::endl;

    //end border regions


    //average on CPU
    start_time = omp_get_wtime();

    //central area

    #pragma omp parallel for
    for(int y = border_half; y < size_y-border_half; y++)
    {
        for(int x = border_half; x < size_x-border_half; x++)
        {
            for(int a = -border_half; a <= border_half; a++)
            {
                //#pragma omp for simd
                for(int b = -border_half; b <= border_half; b++)
                {
                    new_image[y*size_x+x] += (a*a+b*b < circle_radius_squared) * image[(y+a)*size_x+(x+b) ];
                }
            }
            new_image[y*size_x+x] /= number_of_pixels_in_circle;
        }
    }
    std::cout << "Runtime classic: " << omp_get_wtime()-start_time << std::endl;

    saveTiffArray("filtered_tiff.tif",new_image,size_x,size_y);


    return 0;
}

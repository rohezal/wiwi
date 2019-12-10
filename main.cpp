#include <iostream>
#include <omp.h>
#include <png++/png.hpp>
#include "png.h"
//#include <cmath>

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
    //const int size_x = 4096;
    //const int size_y = 4096;
    //const float half_x = size_x*0.5f;
    //const float half_y = size_y*0.5f;

    float *tifdata;
    std::pair<size_t, size_t> tifsize = loadTiffArray("brightness_win1.tif",&tifdata);
    const int tsize_x = tifsize.second;
    const int tsize_y = tifsize.first;

    std::cout << "Width: " << tsize_x << " | Height: " << tsize_y << std::endl;

    saveTiffArray("output.tif", tifdata,tsize_x, tsize_y);

    png::image<png::gray_pixel_16> test("output_png.png");

    std::cout << "Before: fore loop. ysize:" << tsize_y << " | xsize:" << tsize_x << std::endl;
    for(int y = 0; y < tsize_y; y++)
    {
        for(int x = 0; x < tsize_x; x++)
        {
            test[y][x] = tifdata[y*tsize_x+x]*255;
            //test[y][x] = dataxx[y][x]*255;
        }
    }

    std::cout << "Before: test.write(output_png.png)" << std::endl;
    test.write("output_png.png");
    //exit(0);


    png::image<png::gray_pixel_16> light("light.png");
    png::image<png::gray_pixel_16> frequency("frequency.png");
    png::image<png::gray_pixel_16> combined("light.png");

    const int size_x = light.get_width();
    const int size_y = light.get_height();
    const int size_bytes = size_x*size_y;

    float* light_array = new float[size_x*size_y];
    float* frequency_array = new float[size_x*size_y];
    float* combined_array = new float[size_x*size_y];

    float* image = combined_array;
    float* new_image = new float[size_x*size_y];
    float* new_image2 = new float[size_x*size_y];
    //float* new_image2 = combined_array;

    image = tifdata;

    float maximum_frequency = 0;


    //create a random image

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
            /*
            const float val_li =  light[y][x];
            const float val_fre = frequency[y][x];
            light_array[y*size_x+x] = val_li;
            frequency_array[y*size_x+x] = val_fre;
            image[y*size_x+x] = val_li;
            */

            new_image[y*size_x+x] = 0;
            new_image2[y*size_x+x] = 0;
        }
    }


    #pragma omp parallel for
    for(int y = 0; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            combined[y][x] = image[y*size_x+x];

        }
    }

    combined.write("combined_input.png");



    /*
    // Show a debug cirlce
    for(int y = border_half; y < border_half+1; y++)
    {
        for(int x = border_half; x < border_half+1; x++)
        {
            for(int a = -border_half; a < border_half; a++)
            {
                for(int b = -border_half; b < border_half; b++)
                {
                    combined[y+a][x+b] = (a*a+b*b) < circle_radius_squared;
                    combined[y+a][x+b] *= 255;
                }
            }
        }
    }



    combined.write("combined.png");
    exit(0);
    */


    //average on GPU
    double start_time = omp_get_wtime();
    /*
    #pragma omp target data map(alloc:new_image) map(image)
    {
        #pragma omp target teams distribute parallel for collapse(2)
        for(int y = border_half; y < size_y-border_half; y++)
        {
            for(int x = border_half; x < size_x-border_half; x++)
            {
                for(int a = -border_half; a < border_half; a++)
                {
                    for(int b = -border_half; b < border_half; b++)
                    {
                        //new_image[y][x] += image[y+border_half][x+border_half];
                        new_image[y*size_x+x] += (a*a+b*b < circle_radius_squared) ? image[(y+a)*size_x+(x+b) ] : 0;
                    }
                }

                new_image[y*size_x+x] /= number_of_pixels_in_circle;
            }
        }
    }
    */

    std::cout << "Runtime gpu: " << omp_get_wtime()-start_time << std::endl;


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


    //average on CPU
    start_time = omp_get_wtime();
    //#pragma omp parallel for

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

    //=====================================================

    /*
    #pragma omp parallel for
    for(int y = 0; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            combined[y][x] = new_image[y*size_x+x];

        }
    }

    combined.write("combined_output.png");
    */

    return 0;
}

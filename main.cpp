#include <iostream>
#include <omp.h>
#include <png++/png.hpp>

using namespace std;

int main()
{
    //const int size_x = 4096;
    //const int size_y = 4096;
    //const float half_x = size_x*0.5f;
    //const float half_y = size_y*0.5f;


    png::image<png::gray_pixel> light("light.png");
    png::image<png::gray_pixel_16> frequency("frequency.png");
    png::image<png::gray_pixel> combined("light.png");

    const int size_x = light.get_width();
    const int size_y = light.get_height();
    const int size_bytes = size_x*size_y;

    float* light_array = new float[size_x*size_y];
    float* frequency_array = new float[size_x*size_y];
    float* combined_array = new float[size_x*size_y];

    float* image = combined_array;
    float* new_image = new float[size_x*size_y];
    //float* new_image2 = combined_array;

    float maximum_frequency = 0;


    #pragma omp parallel for
    for(int y = 0; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            const float val_li =  light[y][x];
            const float val_fre = frequency[y][x];
            light_array[y*size_x+x] = val_li;
            frequency_array[y*size_x+x] = val_fre;
        }
    }

    float max_light = 0;

    #pragma omp parallel for reduction(max : maximum_frequency)
    for(int i=0;i<size_bytes; i++)
    {
        if(frequency_array[i] > maximum_frequency)
        {
            maximum_frequency = frequency_array[i];
        }
    }

    std::cout << "maximum_frequency: " << maximum_frequency << std::endl;
    maximum_frequency = 1/maximum_frequency;
    std::cout << "maximum_frequency: " << maximum_frequency << std::endl;

    std::cout << "max_light: " << max_light << std::endl;


    #pragma omp parallel for
    for(int y = 0; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            const float val_li =  light_array[y*size_x+x];
            const float val_fre = frequency_array[y*size_x+x];
            const float val_comb = val_li*val_fre*maximum_frequency; //normalizing the frequency
            combined_array[y*size_x+x] = val_comb;
            //combined[y][x] = val_comb;

        }
    }

    //combined.write("combined.png");


    //create a random image

    const int border = 16;
    const int border_half = border/2;
    const int border_squared = border*border;
    const float circle_radius_squared = (border/2)*(border/2);

	int A[1] = {-1};
	#pragma omp target
	{
	  A[0] = omp_is_initial_device();
	}

	if (!A[0])
	{
		std::cout << "Able to use offloading!" << std::endl;
	}

    std::cout << "Devices: " << omp_get_num_devices() << std::endl;



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

                new_image[y*size_x+x] /= border_squared;
            }
        }
    }
    */

    std::cout << "Runtime gpu: " << omp_get_wtime()-start_time << std::endl;

    //average on CPU
    start_time = omp_get_wtime();
    //#pragma omp parallel for
    #pragma omp parallel for
    for(int y = border_half; y < size_y-border_half; y++)
    {
        for(int x = border_half; x < size_x-border_half; x++)
        {
            for(int a = -border_half; a < border_half; a++)
            {
                //#pragma omp for simd
                for(int b = -border_half; b < border_half; b++)
                {
                    new_image[y*size_x+x] += (a*a+b*b < circle_radius_squared) ? image[(y+a)*size_x+(x+b) ] : 0;
                }
            }
            new_image[y*size_x+x] /= border_squared;
        }
    }
    std::cout << "Runtime classic: " << omp_get_wtime()-start_time << std::endl;



    cout << "Hello World!" << endl;
    return 0;
}

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

    png::image< png::gray_pixel> light("light.png");
    png::image< png::gray_pixel> frequency("frequency.png");
    png::image< png::gray_pixel> combined("light.png");

    const int size_x = light.get_width();
    const int size_y = light.get_height();

    float* image = new float[size_x*size_y];
    float* new_image = new float[size_x*size_y];
    float* new_image2 = new float[size_x*size_y];



    exit(0);

    const int border = 32;
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

    //create a random image
    #pragma omp parallel for
    for(int y = 0; y < size_y; y++)
    {
        for(int x = 0; x < size_x; x++)
        {
            //image[y*size_x+x] = 1/(1+((x-half_x)*(x-half_x)))*1/(1+((y-half_y)*(y-half_y)));
            image[y*size_x+x] = 1;
        }
    }


    //average on GPU
    double start_time = omp_get_wtime();
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
                        new_image2[(y+a)*size_x+(x+b) ] = (a*a+b*b < circle_radius_squared) ? image[(y+a)*size_x+(x+b) ] : 0;
                    }
                }

                new_image[y*size_x+x] /= border_squared;
            }
        }
    }

    std::cout << "Runtime gpu: " << omp_get_wtime()-start_time << std::endl;

    //average on CPU
    start_time = omp_get_wtime();
    #pragma omp parallel for
    for(int y = border_half; y < size_y-border_half; y++)
    {
        for(int x = border_half; x < size_x-border_half; x++)
        {
            for(int a = -border_half; a < border_half; a++)
            {
                for(int b = -border_half; b < border_half; b++)
                {
                    //new_image[y*size_x+x] += image[(y+a)*size_x+(x+b) ];
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

#include "png.h"

std::vector<std::vector<float> > loadTiff(const std::string& _name)
{
    std::vector<std::vector<float> > _pixels;
    TIFF* tif = TIFFOpen(_name.c_str(), "r");
    if (tif)
    {
        uint32 imagelength;
        tsize_t scanline;
        tdata_t buf;
        uint32 row;
        uint32 col;
        uint16 depth;
        uint32 width;

        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);

        scanline = TIFFScanlineSize(tif);

        depth = scanline / width;

        uint32_t* temp = new uint32_t[scanline/depth];




        std::cout << "length: " << imagelength << " | scanlinesize: " << scanline  << " | depth: " << depth << std::endl;

        buf = _TIFFmalloc(scanline);
        for (row = 0; row < imagelength; row++)
        {
            _pixels.push_back(std::vector<float>() );
            TIFFReadScanline(tif, buf, row);
            /*
            for (col = 0; col < scanline; col++)
            {
                const uint32_t value = ((uint32_t*)buf)[col];
                _pixels.back().push_back(value);
            }
            */
            _TIFFmemcpy(temp, buf, scanline);

            for (col = 0; col < scanline/depth; col++)
            {
                //_pixels.back().push_back(reverse_bytes(temp[col]));
                float* dirty_hack = (float*) &temp[col];
                _pixels.back().push_back(*dirty_hack);
                if(col == scanline/depth/2)
                {
                    //std::cout << *dirty_hack << std::endl;
                }
            }


        }
        _TIFFfree(buf);
        TIFFClose(tif);
        delete[] temp;
    }
    return _pixels;
}

std::pair<size_t, size_t> loadTiffArray(const std::string &_name, float *array)
{
    std::vector<std::vector<float> > image = loadTiff(_name);

    array = new float[image.size() * image.front().size()];

    int size_x = image.front().size();
    int size_y = image.size();

    for(size_t y = 0; y < size_y;y++)
    {
        for(size_t x = 0; x < size_x; x++)
        {
            array[y*size_x+x] = image[y][x];
        }
    }
    return std::pair<size_t, size_t> (size_y,size_x);
}

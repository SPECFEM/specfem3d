
/*
 *****************************************************************************
 *
 * Description:
 *
 *    Applies a Gaussian filter to an ascii point set.
 *
 *****************************************************************************
 *
 * Author: Joe Henry
 *
 * Copyright 2002 Delta Search Labs
 */

/* DK DK define threshold at which we truncate the basement surface */
/* DK DK check that threshold is the same as in constants.h for mesher */
#define Z_THRESHOLD -4700.

/* DK DK call with

  xfilter inputfile outputfile rows cols kernelSize sigma

  with   kernelSize = 10 gridpoints in each direction
  and    sigma = 2 meters standard deviation on vertical Z data

  the input file is an ASCII file with x y z on each line
  separated by simple white spaces (no commas etc.)

  compile code with gcc using " g++ -O -o xfilter filter_gauss_joe.cpp "

  on the SGI compile with " CC -o xfilter filter_gauss_joe.cpp -lm "

  i.e. for basement of LA basin, call with:

  xfilter reggridbase2_raw_data.dat reggridbase2_filtered.dat  161 144 10 2

 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

class FGauss
{
    public:

        /*!
         *   Constructor.
         *
         *   \param int    width  - Filter width, must be odd.
         *   \param double sigma  - Exponential multiplier (defaults to 1).
         *   \param double aspect - Physical aspect ratio (dy/dx) of data.
         */

        FGauss(int width, double sigma = 1.0, double aspect = 1.0);

        /*!
         *   Destructor.
         */

        ~FGauss();

        /*!
         *   Filter a single sample.
         *
         *   \param int    x  - X coordinate of sample.
         *   \param int    y  - Y coordinate of sample.
         *   \param float *data - The array to filter.
         *   \param int    rows - The array size in y.
         *   \param int    cols - The array size in x.
         *
         *   \return int - The filtered sample.
         */

        double filter(int x, int y, float *data, int rows, int cols);

        /*!
         *   Filter the whole thing.
         *
         *   \param float *dataIn   - The array to filter.
         *   \param float *dataOut  - The array to filter.
         *   \param int    rows     - The array size in y.
         *   \param int    cols     - The array size in x.
         */

        void filter(float *dataIn, float *dataOut, int rows, int cols);

    private:

        double filter(float *origin, int sizex); //!< filter a sample.

        double *_kernel;   //!< filter kernel.
        double  _norm;     //!< normalization factor
        int     _width;    //!< kernel width.
        int     _w2;       //!< width / 2.
};

int
main(int argc, char **argv)
{
    int row, col, rows, cols;
    FILE *fp;

    float sigma      = 2.0;
    int   kernelSize = 10;

    if (argc < 7 ||
        sscanf(argv[3], "%d", &rows) != 1 ||
        sscanf(argv[4], "%d", &cols) != 1 ||
        sscanf(argv[5], "%d", &kernelSize) != 1 ||
        sscanf(argv[6], "%f", &sigma) != 1)
    {
        fprintf(stderr,
               "\nUsage:\n%s inputfile outputfile rows cols kernelSize sigma\n",
                argv[0]);

        return 1;
      }

    printf("Reading %s\n", argv[1]);

    if (!(fp = fopen(argv[1], "r")))
    {
        fprintf(stderr, "\n%s not found.\n\n", argv[1]);
        return 1;
    }

    double meanAspect = 1.0;

    int idummy;
    float zvalue,z1dummy,z2dummy;

    float *x    = new float[rows * cols];
    float *y    = new float[rows * cols];
    float *z    = new float[rows * cols];
    float *fltz = new float[rows * cols];
    float *xp, *yp, *zp;

    int *t = new int[rows * cols];
    int *tp;

    xp = x; yp = y; zp = z; tp = t;

    for (row = 0; row < rows; ++row)
    {
        for (col = 0; col < cols; ++col)
        {

/* DK DK dummy variables read from input file as well */
            if (fscanf(fp, "%d %g %g %g %g %g\n",&idummy,xp,yp,&zvalue,&z1dummy,&z2dummy) != 6)
            {
                fprintf(stderr,
                        "\nRead error on %s at %d, %d.\n\n",
                        argv[1], col, row);
                fclose(fp);
                delete[] x;
                delete[] y;
                delete[] z;
                delete[] fltz;
                delete[] t;
                return 1;
            }

/* DK DK apply threshold and store new z value */
         if (zvalue > Z_THRESHOLD) { zvalue = Z_THRESHOLD; }
         *zp = zvalue;

            ++xp; ++yp; ++zp; ++tp;
        }
    }

    fclose(fp);


    // Allocate a temporary buffer and filter the data into it.

    printf("Filtering %d rows, %d cols, width = %d, sigma = %f\n",
           rows, cols, kernelSize, sigma);

    FGauss gauss(kernelSize, sigma, meanAspect);

    gauss.filter(z, fltz, rows, cols);

    memcpy(z, fltz, (rows * cols) * sizeof(float));


    // Write out the filtered data

    printf("Writing %s\n", argv[2]);

    if (!(fp = fopen(argv[2], "w")))
    {
        fprintf(stderr, "\nCan't open %s for write.\n\n", argv[2]);
        delete[] x;
        delete[] y;
        delete[] z;
        delete[] fltz;
        delete[] t;
        return 1;
    }

    xp = x; yp = y; zp = z; tp = t;

    for (row = 0; row < rows; ++row)
    {
        for (col = 0; col < cols; ++col)
        {
/* DK DK UGLY            if (fprintf(fp, "%f %f %f\n", *xp, *yp, *zp) < 0) */
            if (fprintf(fp, "1 %f %f %f 0 0\n", *xp, *yp, *zp) < 0)
            {
                fprintf(stderr,
                        "\nWrite error on %s at %d, %d.\n\n",
                        argv[2], col, row);
                fclose(fp);
                delete[] x;
                delete[] y;
                delete[] z;
                delete[] fltz;
                delete[] t;
                return 1;
            }

            ++xp; ++yp; ++zp; ++tp;
        }
    }

    fclose(fp);

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] fltz;
    delete[] t;

    return 0;
}

static bool verbose = false;

FGauss::FGauss(int width, double sigma, double aspect)
{
    // Make sure it's odd

    if (!(width % 2))
    {
        ++width;
    }

    _width    = width;
    _w2       = width / 2;
    _kernel   = new double[_width * _width];
    _norm     = 0.0;

    double ec  = 1.0 / (sigma * sigma * width);
    double asq = aspect * aspect;

    // Build the kernel

    for (int i = -_w2; i <= _w2; ++i)
    {
        for (int j = -_w2; j <= _w2; ++j)
        {
            double t = (double)width * exp(-ec * (i * i + j * j * asq));

            _kernel[(i + _w2) * _width + j + _w2] = t;

            _norm += t;
        }
    }

    if (verbose)
    {
        printf("\n");

        for (int i = 0; i < _width; ++i)
        {
            for (int j = 0; j < _width; ++j)
            {
                printf("%g ", (float)_kernel[i * _width + j]);
            }
            printf("\n");
        }
        printf("normfactor = %f\n", (float)_norm);
    }
}

FGauss::~FGauss()
{
    if (_kernel)
    {
        delete[] _kernel;
    }
}

double
FGauss::filter(int x, int y, float *data, int rows, int cols)
{
    if (data)
    {
        if (x >= _w2 && x <= cols - _w2 &&
            y >= _w2 && y <= rows - _w2)
        {
            return filter(data + y * rows + x, cols);
        }
        else if (x && x < cols && y && y < rows)
        {
            return (double)data[y * rows + x];
        }
    }

    return 0;
}

void
FGauss::filter(float *dataIn, float *dataOut, int rows, int cols)
{
    if (dataIn && dataOut)
    {
        // Leave unfiltered data at the border

        memcpy(dataOut, dataIn, cols * rows * sizeof(float));

        printf("\n");

        float *dip    = dataIn + _w2 * cols + _w2;
        float *dop    = dataOut + _w2 * cols + _w2;
        float *endCol = dataIn + _w2 * cols + cols - _w2;
        float *endRow = dataIn + (rows - _w2) * cols;

        int nRows = rows - _width;

        for (int row = 0; dip < endRow; row += 100)
        {
            while (dip < endCol)
            {
                *dop++ = (float)filter(dip++, cols);
            }

            dip    += _width - 1;
            dop    += _width - 1;
            endCol += cols;

            printf("\r%d%%", row / nRows);
            fflush(stdout);
        }
        printf("\n");
    }
}

double
FGauss::filter(float *origin, int cols)
{
    double sample, sum;

    double *kp     = _kernel;
    double *endCol = kp + _width;
    double *endRow = kp + _width * _width;
    float  *dp     = origin - (_w2 * cols) - _w2;

    sum = 0.0;

    while (kp < endRow)
    {
        while (kp < endCol)
        {
            sample  = (double)(*dp++);
            sum    += sample * *kp++;
        }
        dp     += cols - _width;
        endCol += _width;
    }

    return sum / _norm;
}


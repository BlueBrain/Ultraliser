#pragma once 

#include <cmath>
#include <cstdio>
#include <cstdlib>
namespace Ultraliser
{


/**
 * @brief The Mask class
 *  Class  used  to  store a mask and its  rotations (90, 180 and 270 degrees)
 */
class Mask
{
public:

    /**
     * @brief The AXIS enum
     */
    enum AXIS { X, Y, Z };

    /**
     * @brief The DIRECTION enum
     */
    enum DIRECTION { U, D, E, W, N, S, UNSPECIFIED };

    Mask();
    ~Mask();

    void printMask();
    void printMasks();

    void generateRotations();					 // Generates  the rotated  masks
    void printDirection();

    void setDirection(char d);

    /**
     * @brief updateDirection
     * @param newDirection
     */
    void updateDirection(DIRECTION newDirection);

    // Tests if one of the four masks  (0, 90, 180 and 270) matchs with the neighborhood of the
    // voxel p  in the volume Vol
    bool matches(const int8_t *subVolume);


    // Generates the four masks from the mask "umask" in up direction
    void set_mask_from_u(int8_t umask[]);



private:
    int8_t *_mask0;
    int8_t *_mask90;
    int8_t *_mask180;
    int8_t *_mask270;

    char direction;

    DIRECTION _direction;

    // Static functions to match and rotate the vectors
    static bool matchf(int ***Vol, int *vec);
    bool matchf1d(int*vol, int *vec);
    bool matches(const int8_t *vol, int8_t *vec);



    static void rotate(int8_t vector[], int8_t *vector_rot, char axis);
};

}

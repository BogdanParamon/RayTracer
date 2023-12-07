#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    float x = floor(texCoord[0] * image.width), y = floor((1 - texCoord[1]) * image.height);

    x = std::clamp(x, 0.f, (float)image.width);
    y = std::clamp(y, 0.f, (float)image.height);

    return image.pixels[x + y * image.width];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    if (texCoord.y < 0.0f || texCoord.y > 1.0f || texCoord.x < 0.0f || texCoord.x > 1.0f) {
        return glm::vec3(0.0f);
    }

    float x = texCoord[0] * image.width, y = (1 - texCoord[1]) * image.height;

    x = std::clamp(x - 0.5f, 0.f, (float)image.width - 1);
    y = std::clamp(y - 0.5f, 0.f, (float)image.height - 1);

    float t1 = x - floor(x);
    float t2 = y - floor(y);
    glm::vec3 c1, c2, c3, c4;

    c1 = image.pixels[floor(y) * image.width + floor(x)];
    c2 = image.pixels[ceil(y) * image.width + ceil(x)];
    c3 = image.pixels[floor(y) * image.width + ceil(x)];
    c4 = image.pixels[ceil(y) * image.width + floor(x)];

    glm::vec3 result = c1 * ((1 - t2) * (1 - t1)) + c2 * (t1 * t2) + c3 * ((1 - t2) * t1) + c4 * (t2 * (1 - t1));

    return result;
}
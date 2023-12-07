#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "scene.h"
#include "texture.h"
#include <framework/trackball.h>
#include <iostream>
#include <intersect.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    float lensSize = (float) features.extra.lensSize;
    float focalLength = (float) features.extra.focalLength;
    int numSamples = features.extra.numSamplesDOF;

    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            // Generate a ray from each pixel to the focal point
            glm::vec2 position = (glm::vec2(x, y) + 0.5f) / glm::vec2(screen.resolution()) * 2.f - 1.f;
            Ray cameraToFocus = camera.generateRay(position);

            glm::vec3 color = {};
            glm::vec3 focalPoint = cameraToFocus.origin + glm::normalize(cameraToFocus.direction) * focalLength;

            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };

            for (int h = 0; h < numSamples; h++) {
                // Generate rays from random positions inside the pixel to the focal point
                Ray ray;
                ray.t = cameraToFocus.t;
                ray.direction = glm::normalize(-cameraToFocus.origin + focalPoint);

                float random1 = state.sampler.next_1d() * (state.sampler.next_1d() > 0.5f ? 1 : -1);
                float random2 = state.sampler.next_1d() * (state.sampler.next_1d() > 0.5f ? 1 : -1);
                ray.origin = camera.position() + random1 * glm::normalize(camera.up()) * 0.5f * lensSize + random2 * glm::normalize(camera.left()) * 0.5f * lensSize;

                color += renderRay(state, ray);
            }

            color /= (float)numSamples;
            screen.setPixel(x, y, color);            
        }
    }
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {

        glm::vec3 minCoord = glm::vec3(-1.0f, -1.0f, -1.0f);
        glm::vec3 maxCoord = glm::vec3(1.0f, 1.0f, 1.0f);

        AxisAlignedBox aabb = { .lower = minCoord, .upper = maxCoord };
        Ray copy;
        copy.direction = glm::normalize(ray.direction);

        intersectRayWithShape(aabb, copy);
        glm::vec3 intersection = copy.origin + copy.direction * copy.t;
        const auto& [x, y, z] = std::tie(intersection.x, intersection.y, intersection.z);

        int index;

        float axis = std::max(fabs(x), std::max(fabs(y), fabs(z)));
        bool positiveX = x > 0.0f;
        bool positiveY = y > 0.0f;
        bool positiveZ = z > 0.0f;

        float projU, projV;

        if (axis == fabs(x) && positiveX) {
            projU = -z;
            projV = y;
            index = 0;
        }
        if (axis == fabs(x) && !positiveX) {
            projU = z;
            projV = y;
            index = 1;
        }
        if (axis == fabs(y) && positiveY) {
            projU = x;
            projV = -z;
            index = 2;
        }
        if (axis == fabs(y) && !positiveY) {
            projU = x;
            projV = z;
            index = 3;
        }
        if (axis == fabs(z) && positiveZ) {
            projU = x;
            projV = y;
            index = 4;
        }
        if (axis == fabs(z) && !positiveZ) {
            projU = -x;
            projV = y;
            index = 5;
        }

        projU = (projU / axis + 1.0f) / 2.0f;
        projV = (projV / axis + 1.0f) / 2.0f;

        return sampleTextureNearest(state.scene.images[index], { projU, projV });

    } else {
        return glm::vec3(0.0f);
    }
}

// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    std::sort(primitives.begin(), primitives.end(), [axis](const Primitive& a, const Primitive& b) {
        return computePrimitiveCentroid(a)[axis] < computePrimitiveCentroid(b)[axis];
    });

    AxisAlignedBox leftAABB;
    AxisAlignedBox rightAABB;
    int partition = 0;
    size_t currentBestPartition;
    float cost = FLT_MAX;

    std::vector<std::vector<Primitive>> bins(16);
    float k0 = aabb.lower[axis];
    float k1 = (16.0f * (1 - FLT_EPSILON)) / (aabb.upper[axis] - aabb.lower[axis]);

    for (Primitive p : primitives) {
        int binIndex = k1 * (computePrimitiveCentroid(p)[axis] - k0);
        bins[binIndex].push_back(p);
    }
    size_t size = primitives.size();

    while (partition < 15) {
        std::vector<Primitive> leftPrimitives;
        std::vector<Primitive> rightPrimitives;
        for (int i = 0; i <= 15; i++) {
            if (i <= partition)
                leftPrimitives.insert(leftPrimitives.end(), bins[i].begin(), bins[i].end());
            else
                rightPrimitives.insert(rightPrimitives.end(), bins[i].begin(), bins[i].end());
        }

        if (rightPrimitives.size() > 0 && leftPrimitives.size() > 0) {
            leftAABB = computeSpanAABB(leftPrimitives);
            rightAABB = computeSpanAABB(rightPrimitives);

            glm::vec3 diagonalLeftAABB = leftAABB.upper - leftAABB.lower;
            float areaOfLeftAABB = 2.0f * (diagonalLeftAABB.x * diagonalLeftAABB.y + diagonalLeftAABB.y * diagonalLeftAABB.z + diagonalLeftAABB.x * diagonalLeftAABB.z);

            glm::vec3 diagonalRightAABB = rightAABB.upper - rightAABB.lower;
            float areaOfRightAABB = 2.0f * (diagonalRightAABB.x * diagonalRightAABB.y + diagonalRightAABB.y * diagonalRightAABB.z + diagonalRightAABB.x * diagonalRightAABB.z);

            float currentCost = areaOfLeftAABB * leftPrimitives.size() + areaOfRightAABB * rightPrimitives.size();
            if (currentCost < cost) {
                cost = currentCost;
                currentBestPartition = leftPrimitives.size();
            }

        }
        partition++;
    }

    return currentBestPartition;
}
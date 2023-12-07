#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.

    //Compute the position and color using linear interpolation.
    position    =   (1 - sample) * light.endpoint0  +   sample * light.endpoint1;
    color       =   (1 - sample) * light.color0     +   sample * light.color1;
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.

    //Compute the position using bilinear interpolation
    position    =   (1 - sample.x)  *   (1 - sample.y)  *   (light.v0)
                +   (sample.x)      *   (1 - sample.y)  *   (light.v0+light.edge01)
                +   (1 - sample.x)  *   (sample.y)      *   (light.v0+light.edge02)
                +   (sample.x)      *   (sample.y)      *   (light.v0+light.edge01+light.edge02);
    
    //Compute the color using bilinear interpolation
    color       =   (1 - sample.x)  *   (1 - sample.y)  *   (light.color0)
                +   (sample.x)      *   (1 - sample.y)  *   (light.color1)
                +   (1 - sample.x)  *   (sample.y)      *   (light.color2)
                +   (sample.x)      *   (sample.y)      *   (light.color3);
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer
        // TODO: implement this function; currently, the light simply passes through

        // Set up a testing ray at the light source towards the intersection
        Ray r { lightPosition, glm::normalize(ray.origin + ray.t * ray.direction - lightPosition) };

        // Set up the maximal depth parameter (i.e. the maximum depth in the lightview before intersecting with the current surface)
        double t_max = glm::length(ray.origin + ray.t * ray.direction - lightPosition);

        // Check that the surface is facing the light
        //if (glm::dot(r.direction, hitInfo.normal) < 0) {
        //    return false;
        //}
        
        
        // Determine whether there obstructing surfaces between the intersection and the light position
        HitInfo temp;
        state.bvh.intersect(state, r, temp);

        drawRay(Ray { ray.origin + ray.t * ray.direction, glm::normalize(lightPosition - (ray.origin + ray.t * ray.direction)) }, lightColor);


        // Determine whether the obstruction occurs before the surface
        // Compute using an epsilon as a rounding correction. (arbitrarily chosen in this case).
        if (t_max - r.t < 0.0001f) {
            return true;
        }

        return false;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: implement this function; currently, the light simply passes through
    
    // Set up a testing ray at the light position towards the intersection
    Ray r {lightPosition, glm::normalize(ray.origin + ray.t * ray.direction - lightPosition) };

    // Set up the maximal depth parameter (i.e. the maximum depth in the lightview before intersecting with the current surface)
    double t_max = glm::length(ray.origin + ray.t * ray.direction - lightPosition);

    // Check that the surface is facing the light
    // if (glm::dot(r.direction, hitInfo.normal) < 0) {
    //    return glm::vec3(0.0);
    //}

    // Declare a new light color vector
    glm::vec3 color = lightColor;

    // Compute iteratively the effect of obstructing surfaces
    HitInfo obstructingHitInfo;

    
    while (state.bvh.intersect(state, r, obstructingHitInfo) && t_max - r.t > 0.0001) {

        
        // Determine whether the obstruction occurs before the surface
        // Compute using an epsilon as a rounding correction. (arbitrarily chosen in this case).
        color *= obstructingHitInfo.material.kd * obstructingHitInfo.material.transparency;

        r = Ray { r.direction * r.t + r.origin + r.direction * 0.000001f, r.direction };
    }

    drawRay(Ray { ray.origin + ray.t * ray.direction, lightPosition - (ray.origin + ray.t * ray.direction) }, color);


    return color;
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: modify this function to incorporate visibility corerctly
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;
    return computeShading(state, v, l, visibilityOfLightSample(state,light.position,light.color,ray,hitInfo), hitInfo);
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model

    glm::vec3 color;
    glm::vec3 position;
    glm::vec3 Lo;
    glm::vec3 shade { 0, 0, 0 };

    for (int i = 0; i < numSamples; i++) {
        // Sample the segment light
        sampleSegmentLight(state.sampler.next_1d(), light, position, color);

        // Compute the visible light
        Lo = visibilityOfLightSample(state, position, color, ray, hitInfo);

        // Evaluate the Phong Model
        shade += computePhongModel(state, -ray.direction, position - (ray.direction * ray.t + ray.origin), Lo, hitInfo);
    }

    shade = shade / (float)numSamples;
    return shade;
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection 
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model
    

    glm::vec3 color;
    glm::vec3 position;
    glm::vec3 Lo;
    glm::vec3 shade { 0, 0, 0 };


   

    for (int i = 0; i < numSamples; i++){
        // Sample the parallelogram light
        sampleParallelogramLight(state.sampler.next_2d(), light, position, color);

        // Compute the visible light
        Lo = visibilityOfLightSample(state, position, color, ray, hitInfo);

        // Evaluate the Phong Model
        shade += computePhongModel(state, -ray.direction, position - (ray.direction * ray.t + ray.origin), Lo, hitInfo);
    }

    shade = shade / (float) numSamples;
    return shade;
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}
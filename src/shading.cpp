#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded linear gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            { 0.1f, glm::vec3(215.f / 256.f, 210.f / 256.f, 203.f / 256.f) },
            { 0.22f, glm::vec3(250.f / 256.f, 250.f / 256.f, 240.f / 256.f) },
            { 0.5f, glm::vec3(145.f / 256.f, 170.f / 256.f, 175.f / 256.f) },
            { 0.78f, glm::vec3(255.f / 256.f, 250.f / 256.f, 205.f / 256.f) },
            { 0.9f, glm::vec3(170.f / 256.f, 170.f / 256.f, 170.f / 256.f) },
        }
    };

    if (state.features.enableShading) {
        switch (state.features.shadingModel) {
            case ShadingModel::Lambertian:
                return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::Phong:
                return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::BlinnPhong:
                return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::LinearGradient:
                return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 result = {}, newLight = glm::normalize(lightDirection), newNormal = glm::normalize(hitInfo.normal);

    float cosAngle = glm::dot(newLight, newNormal);
    if (cosAngle < 0) {
        return glm::vec3(0.0f);
    }

    return lightColor * sampleMaterialKd(state, hitInfo) * cosAngle;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 specular = {}, newLight = glm::normalize(lightDirection), newNormal = glm::normalize(hitInfo.normal);
    glm::vec3 diffuse = computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);

    float cosAngle = glm::dot(newLight, newNormal);
    if (cosAngle < 0) {
        return glm::vec3(0.0f);
    }

    glm::vec3 reflection = glm::reflect(-newLight, newNormal);
    float cosCamera = glm::dot(glm::normalize(cameraDirection), glm::normalize(reflection));

    if (cosCamera < 0) {
        specular = glm::vec3(0.0f);
    } else {
        specular = lightColor * hitInfo.material.ks * float(pow(cosCamera, hitInfo.material.shininess));
    }

    glm::vec3 result = diffuse + specular;
    return result;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    glm::vec3 specular = {}, newLight = glm::normalize(lightDirection), newNormal = glm::normalize(hitInfo.normal);
    glm::vec3 diffuse = computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);

    float cosAngle = glm::dot(newLight, newNormal);
    if (cosAngle < 0) {
        return glm::vec3(0.0f);
    }

    glm::vec3 halfway = glm::normalize(newLight + cameraDirection);
    float cosCamera = glm::dot(halfway, newNormal);
    specular = lightColor * hitInfo.material.ks * float(pow(cosCamera, hitInfo.material.shininess));

    glm::vec3 result = diffuse + specular;

    return result;
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    float min = -2, max = 2, indexMin = -1, indexMax = -1;
    
    for (int i = 0, size = components.size(); i < size; i++) {
        if (ti >= components[i].t && min < components[i].t) {
            min = components[i].t;
            indexMin = i;
        } else if (ti <= components[i].t && max > components[i].t) {
            max = components[i].t;
            indexMax = i;
        }
    }

    if (ti == min || indexMax == -1) {
        return components[indexMin].color;
    }
    if (ti == max || indexMin == -1) {
        return components[indexMax].color;
    }

    float t1 = components[indexMin].t, t2 = components[indexMax].t;
    glm::vec3 result = components[indexMin].color * ((max - ti) / (max - min)) + components[indexMax].color * ((ti - min) / (max - min));

    return result;
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);

    if (cos_theta < 0) {
        return glm::vec3(0.0f);
    }

    glm::vec3 result = lightColor * cos_theta * gradient.sample(cos_theta);

    return result;
}
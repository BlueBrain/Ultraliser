#version 450 core
#extension GL_ARB_bindless_texture : enable

// Object space position
layout (location = 0) in vec3 inLocalPos;
// Camera space position
layout (location = 1) in vec3 inCamSpacePos;

layout (location = 0) out vec4 outColor;

// Object space light direction
uniform vec3 lightDirection;
//const vec3 lightDirection = vec3(-1.0, 1.0, 1.0);
// Object space volume max and min bounds
uniform vec3 localMinBound;
uniform vec3 localMaxBound;

uniform mat4 invModelView;

layout(bindless_sampler) uniform sampler3D volumeTexture;

#define MARCH_SAMPLES 500
#define SAMPLE_HIT_EPSILON 0.01
#define MAX_FIX_ITERATIONS 100

/**
 * Returns true wether a normalized position is in range 0,0,0 
 * to 1,1,1
 */
bool intersectVolume(in vec3 normalizedPos)
{
  return normalizedPos.x >= 0.0 && normalizedPos.x <= 1.0
      && normalizedPos.y >= 0.0 && normalizedPos.y <= 1.0
      && normalizedPos.z >= 0.0 && normalizedPos.z <= 1.0;
}

/**
 * Attempts to intersect the volume by testing a ray with a given
 * origin and a direction. Stores the result position (if intersected)
 * in the out paramter hitPos
 */
bool intersectLocalVolume(in vec3 rayO, in vec3 rayDir, out vec3 hitPos)
{
  float tmin, tmax, tymin, tymax, tzmin, tzmax; 
 
  const float minBoundX = rayDir.x < 0.0? localMaxBound.x : localMinBound.x;
  const float maxBoundX = rayDir.x < 0.0? localMinBound.x : localMaxBound.x;
  const float minBoundY = rayDir.y < 0.0? localMaxBound.y : localMinBound.y;
  const float maxBoundY = rayDir.y < 0.0? localMinBound.y : localMaxBound.y;
  const float minBoundZ = rayDir.z < 0.0? localMaxBound.z : localMinBound.z;
  const float maxBoundZ = rayDir.z < 0.0? localMinBound.z : localMaxBound.z;

  const vec3 invDir = 1.0 / rayDir;

  tmin = (minBoundX - rayO.x) * invDir.x; 
  tmax = (maxBoundX - rayO.x) * invDir.x; 
  tymin = (minBoundY - rayO.y) * invDir.y; 
  tymax = (maxBoundY - rayO.y) * invDir.y; 

  if ((tmin > tymax) || (tymin > tmax)) 
      return false; 
    
  tmin = tymin > tmin? tymin : tmin;
  tmax = tymax < tmax? tymax : tmax;

  tzmin = (minBoundZ - rayO.z) * invDir.z; 
  tzmax = (maxBoundZ - rayO.z) * invDir.z; 

  if ((tmin > tzmax) || (tzmin > tmax)) 
      return false;

  tmin = tzmin > tmin? tzmin : tmin;
  tmax = tzmax < tmax? tzmax : tmax;

  hitPos = rayO + rayDir * tmax;

  return tmax >= 0.0; 
}

/**
 * Normalizes a position inside the volume within the
 * range 0,0,0 to 1,1,1
 */
vec3 normalizeVolumePos(in vec3 localPos)
{
  return (localPos - localMinBound) / (localMaxBound - localMinBound);
}

/**
 * Given a volume entry point, find the exiting point to
 * compute the most appropiate ray marching distance
 */
float findAppropiateDistance(in vec3 startPoint, in vec3 dir)
{
  vec3 localEndPoint;
  intersectLocalVolume(startPoint, dir, localEndPoint);

  return length(localEndPoint - startPoint);
}

/**
 * Binary adjustement to set the hit point, whithin the volume, to the
 * closest point of the volume data surface. Defined margin is
 * SAMPLE_HIT_EPSILON
 */
vec3 fixHitPoint(in vec3 hitPoint, in vec3 dir, in float stepSize)
{
  // Make sure the input positions can sample volume data at LOD 0
  // Otherwise, adjust them
  bool realHit = false;
  vec3 end = hitPoint;
  float tempStepSize = stepSize / 10.0;
  while(!realHit)
  {
    vec3 temp = normalizeVolumePos(end);
    realHit = intersectVolume(temp) && textureLod(volumeTexture, temp, 0.0).r > 0.0;
    end = !realHit? end + dir * tempStepSize : end;
  }
  vec3 start = end - dir * stepSize;

  // Binary search for the real boundary
  bool goodPos = false;
  int i = 0;
  while(!goodPos && i < MAX_FIX_ITERATIONS)
  {
    vec3 middle = (start + end) / 2.0;
    const vec3 normalizedMiddle = normalizeVolumePos(middle);
    if(intersectVolume(normalizedMiddle))
    {
      const float val = textureLod(volumeTexture, normalizedMiddle, 0.0).r;
      end = val >= 0.0? middle : end;
      start = val < 0.0? middle : start;
    }
    else
    {
      start = middle;
    }

    goodPos = length(end - start) <= SAMPLE_HIT_EPSILON;
    i++;
  }

  // The closest OUTSIDE the surface to avoid lighting artifacts
  return start;
}

/*
 * Main raymarch algorithm.
 * Returns the sampled value from the volume data and sets
 * the out parameter hitPos to the local space hit position.
 * Returns 0.0 if the volume wasn't hit
 */
float raymarch(in vec3 pos, 
               in vec3 dir,
               out vec3 hitPos)
{
  const float stepSize = findAppropiateDistance(pos, dir) / float(MARCH_SAMPLES);

  int i = 0;
  while(i < MARCH_SAMPLES)
  {
    hitPos = pos + dir * stepSize * i;
    vec3 normHitPos = normalizeVolumePos(hitPos);

    if(intersectVolume(normHitPos))
    {
      const float volSample = texture(volumeTexture, normHitPos).r;//textureLod(volumeTexture, normHitPos, 5.0).r;
      if(volSample > 0.0)
      {
        hitPos = fixHitPoint(hitPos, dir, stepSize);
        return volSample;
      }
    }
    
    i++;
  }

  return 0.0;
}

/**
 * Peforms a raymarch from a valid hit position within the volume
 * towards the light direction to compute shadows
 */
vec4 shadeSample(in vec3 volumeHitPos)
{
  vec3 dummy;
  const float shadeSample = raymarch(volumeHitPos, normalize(lightDirection), dummy);
  return vec4(vec3(1,0,0) * (shadeSample > 0.0? 0.5 : 1.0), 1.0);
}

/**
 * Main
 */
void main( void )
{
  const vec3 localSpaceRayDir = normalize((invModelView * vec4(normalize(inCamSpacePos), 0.0)).xyz);
  
  vec3 rayHitPos;
  float s = raymarch(inLocalPos, localSpaceRayDir, rayHitPos);
  outColor = s > 0.0? shadeSample(rayHitPos) : vec4(0.0);
  //outColor = vec4(1,0,0,1);
}
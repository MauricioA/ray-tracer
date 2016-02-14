#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <random>
#include "ray.h"

// g++ test.cc -o test -std=c++11 -O3 -Wall
// ./test > test.ppm
// convert test.ppm test.png

// TODO paralelize
// TODO some unitVector() calls might be redundant
// TODO stop bouncing if color is close to black

vec3 blend(const vec3& u, const vec3& v, float t) {
  return (1-t) * u + t * v;
}

vec3 background(const Ray& r) {
  vec3 unit_dir = r.B.unitVector();
  float t = 0.5 + unit_dir.y() * 0.5;
  return blend(vec3(1, 1, 1), vec3(0.5, 0.7, 1.0), t);
}

vec3 color(const Ray& r, const Hitable& world, int depth) {
  HitRecord record;
  vec3 att;
  Ray r_out;
  if (world.hit(r, 0.001, FLT_MAX, record)) {
    if (depth > 0 && record.mat->scatter(r, record, att, r_out)) {
      vec3 rcol = color(r_out, world, depth - 1);
      return vec3(att.r() * rcol.r(), att.g() * rcol.g(), att.b() * rcol.b());
    } else {
      return vec3(0, 0, 0);
    }
  } else {
    return background(r);
  }
}

inline uint8_t gicol(float f) {
  return int(255.99 * sqrt(f));
}

HitableList randomScene() {
  const int limit = 11;
  HitableList world;
  world.add(new Sphere(vec3(0, -1000, 0), 1000, 
            new Lambertian(vec3(0.5, 0.5, 0.5)))); // floor
  for (int a = -limit; a < limit; ++a) {
    for (int b = -limit; b < limit; ++b) {
      float mat = dis(gen);
      vec3 center(a + 0.9 * dis(gen), 0.2, b + 0.9 * dis(gen));
      if ((center - vec3(4, 0.2, 0)).norm2() > 0.9) {
        if (mat < 0.8) {
          Lambertian* lamb = new Lambertian(vec3(
              dis(gen)*dis(gen), dis(gen)*dis(gen), dis(gen)*dis(gen)));
          world.add(new Sphere(center, 0.2, lamb));
        } else if (mat < 0.95) {
          Metal* met = new Metal(vec3(0.5 * (1+dis(gen)), 0.5 * (1+dis(gen)),
                                      0.5 * (1+dis(gen))), 0.5 * (dis(gen)));
          world.add(new Sphere(center, 0.2, met));
        } else {
          world.add(new Sphere(center, 0.2, new Dielectric(1.5)));
        }
      }
    }
  }

  world
    .add(new Sphere(vec3(0, 1, 0), 1, new Dielectric(1.5)))
    .add(new Sphere(vec3(-4, 1, 0), 1, new Lambertian(vec3(0.4, 0.2, 1))))
    .add(new Sphere(vec3(4, 1, 0), 1, new Metal(vec3(0.7, 0.6, 0.5), 0)));

  return world;
}

int main() {
  const int nx = 800;
  const int ny = 400;
  const float nxf = float(nx);
  const float nyf = float(ny);
  const int nsamples = 300;
  const int max_depth = 50;
  gen = std::mt19937(dv());
  dis = std::uniform_real_distribution<>(0, 1);
  
  vec3 loofkfrom(15, 2, 4);
  vec3 lookat(4, 1, 0);
  float dist_focus = (loofkfrom - lookat).norm2();
  float aperture = 0.1;
  Camera cam(loofkfrom, lookat, vec3(0, 1, 0), 20,
             nxf / nyf, aperture, dist_focus);  
  HitableList world = randomScene();

  printf("P3\n%d %d\n255\n", nx, ny);
  for (int j = ny-1; j >= 0; --j) {
    for (int i = 0; i < nx; ++i) {
      vec3 col(0, 0, 0);
      for (int s = 0; s < nsamples; ++s) {
        float u = (i + dis(gen)) / nxf;
        float v = (j + dis(gen)) / nyf;
        col += (1.0 / nsamples) * color(cam.getRay(u, v), world, max_depth);
      }
      printf("%d %d %d\n", gicol(col.r()), gicol(col.g()), gicol(col.b()));
    }
    fprintf(stderr, ".");
  }
  fprintf(stderr, "\n");
}

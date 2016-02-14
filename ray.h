#pragma once

#include <cassert>
#include <math.h>
#include <memory>
#include <random>
#include <vector>

std::random_device dv;
std::mt19937 gen;
std::uniform_real_distribution<> dis;

class vec3 {
 public:
  vec3() {};

  vec3(float e0, float e1, float e2) {
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
  }

  inline float x() const { return e[0]; }
  inline float y() const { return e[1]; }
  inline float z() const { return e[2]; }
  inline float r() const { return e[0]; }
  inline float g() const { return e[1]; }
  inline float b() const { return e[2]; }

  inline float& operator[](int i) { return e[i]; }
  inline float operator[](int i) const { return e[i]; }

  inline vec3 operator+(const vec3& v) const {
    return vec3(e[0] + v[0], e[1] + v[1], e[2] + v[2]);
  }

  inline vec3 operator-(const vec3& v) const {
    return vec3(e[0] - v[0], e[1] - v[1], e[2] - v[2]);
  }

  // dot product
  inline float operator*(const vec3& v) const {
    return e[0] * v[0] + e[1] * v[1] + e[2] * v[2];
  }

  inline vec3 operator*(float f) const {
    return vec3(e[0]*f, e[1]*f, e[2]*f);
  }

  inline void operator+=(const vec3& v) {
    e[0] += v[0];
    e[1] += v[1];
    e[2] += v[2];
  }

  inline void operator*=(float f) {
    e[0] *= f;
    e[1] *= f;
    e[2] *= f;
  }

  inline vec3 operator+() const {
    return vec3(e[0], e[1], e[2]);
  }

  inline vec3 operator-() const {
    return vec3(-e[0], -e[1], -e[2]);
  }

  inline float norm2() const {
    return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  }

  inline vec3 unitVector() const {
    float norm = norm2();
    return vec3(e[0] / norm, e[1] / norm, e[2] / norm);
  }

 private:
  float e[3];
};

inline vec3 operator*(float t, const vec3& B) {
  return vec3(B[0]*t, B[1]*t, B[2]*t);
}

inline vec3 cross(const vec3& u, const vec3& v) {
  return vec3(
    +u[1]*v[2] - u[2]*v[1],
    -u[0]*v[2] + u[2]*v[0],
    +u[0]*v[1] - u[1]*v[0]
  );
}

class Ray {
 public:
  Ray() {};

  Ray(const vec3& A, const vec3& B) : A(A), B(B) {};

  inline vec3 getPoint(float t) const {
    return A + t*B;
  }

  vec3 A, B;
};

class Material;

struct HitRecord {
  HitRecord() {};
  float t;
  vec3 P;
  vec3 N;
  Material* mat;
};

class Material {
 public:
  // If the ray was reflected it indicates how much it gets
  // attenuated and where it reflects, and returns true.
  virtual bool scatter(const Ray& r_in, const HitRecord& record, 
                       vec3& attenuation, Ray& r_out) const = 0;
};

vec3 randomPointSphere() {
  vec3 p;
  do {
    p = vec3(dis(gen)*2-1, dis(gen)*2-1, dis(gen)*2-1);
  } while (p * p >= 1);
  return p;
}

class Lambertian : public Material {
 public:
  Lambertian(vec3(albedo)) : albedo(albedo) {};

  // Always reflect, but attenuates
  bool scatter(const Ray& r_in, const HitRecord& record, 
               vec3& attenuation, Ray& r_out) const override;

 private:
  vec3 albedo;
};

bool Lambertian::scatter(const Ray& r_in, const HitRecord& record,
                         vec3& attenuation, Ray& r_out) const {
  vec3 target = record.P + record.N + randomPointSphere();
  r_out.A = record.P;
  r_out.B = target - record.P;
  attenuation = albedo;
  return true;
}

inline vec3 reflect(const vec3& v, const vec3& n) {
  return v - 2 * (v * n) * n;
}

class Metal : public Material {
 public:
  Metal(const vec3& albedo, float f) : 
        albedo(albedo), f(f) {};

  bool scatter(const Ray& r_in, const HitRecord& record, 
               vec3& attenuation, Ray& r_out) const override;

 private:
  vec3 albedo;
  float f;
};

bool Metal::scatter(const Ray& r_in, const HitRecord& record, 
                    vec3& attenuation, Ray& r_out) const {
  attenuation = albedo;
  r_out.A = record.P;
  vec3 v = r_in.B.unitVector();
  r_out.B = reflect(v, record.N);
  r_out.B += f * randomPointSphere();
  return r_out.B * record.N > 0;
}

class Dielectric : public Material {
 public:
  Dielectric(float ref) : ref(ref) {};

  bool scatter(const Ray& r_in, const HitRecord& record, 
               vec3& attenuation, Ray& r_out) const override;

 private:
  inline float schlick(const vec3& normal, const vec3& incident,
                       float ni, float no) const;

  inline bool refract(const vec3& v, const vec3& normal, float ni, float no,
                      vec3& refracted) const;

  float ref;
};

bool Dielectric::scatter(const Ray& r_in, const HitRecord& record,
                         vec3& attenuation, Ray& r_out) const {
  attenuation = vec3(1, 1, 1);
  vec3 out_normal;
  vec3 refracted;
  float ni, no;
  bool should_reflect;
  if ((r_in.B * record.N) > 0) {
    // From inside the sphere to outside
    out_normal = -record.N; 
    ni = ref;
    no = 1.0;
  } else {
    // From outside to inside
    out_normal = +record.N;
    ni = 1.0;
    no = ref;
  }
  if (refract(r_in.B, out_normal, ni, no, refracted)) {
    should_reflect = dis(gen) < schlick(record.N, r_in.B, ni, no);
  } else {
    should_reflect = true;
  }
  if (should_reflect) {
    r_out = Ray(record.P, reflect(r_in.B, record.N));
  } else {
    r_out = Ray(record.P, refracted);
  }
  return true;
}

float Dielectric::schlick(const vec3& normal, const vec3& incident,
                          float ni, float no) const {
  vec3 normal_u = normal.unitVector();
  vec3 incident_u = incident.unitVector();
  float r0 = (ni - no) / (ni + no);
  r0 *=r0;
  float cosX = -(normal_u * incident_u);
  if (ni > no) {
    float n = ni / no;
    float sinT2 = n * n * (1.0 - cosX * cosX);
    if (sinT2 > 1) return 1;
    cosX = sqrt(1 - sinT2);
  }
  float x = 1 - cosX;
  return r0 + (1 - r0) * x * x * x * x * x;
}

bool Dielectric::refract(const vec3& v, const vec3& normal, float ni, float no,
                         vec3& refracted) const {
  vec3 uv = v.unitVector();
  float dt = uv * normal;
  float n = ni / no;
  float disc = 1 - n * n * (1 - dt * dt);
  if (disc > 0) {
    refracted = n * (uv - normal * dt) - normal * sqrt(disc);
    return true;
  } else {
    return false;
  }
}

class Hitable {
 public:
  virtual bool hit(const Ray& r, float t_min, float t_max, 
                   HitRecord& rec) const = 0;
};

class Sphere : public Hitable {
 public:
  Sphere(const vec3& C, float R, Material *mat) : C(C), R(R), mat(mat) {};

  bool hit(const Ray& r, float t_min, float t_max, 
           HitRecord& rec) const override;

  vec3 C;
  float R;
  std::unique_ptr<Material> mat;
};

bool Sphere::hit(const Ray& r, float t_min, float t_max, 
                 HitRecord& rec) const {
  vec3 amC = r.A - C;
  float a2 = 2 * r.B * r.B;
  float b = 2 * amC * r.B;
  float c = amC * amC - R*R;
  float sub = b*b - 2*a2*c;
  for (float mult : {-1, 1}) {
    float t = (-b + mult * sqrt(sub)) / a2;
    if (t >= t_min && t < t_max) {
      rec.t = t;
      rec.P = r.getPoint(t);
      rec.N = (rec.P - C).unitVector();
      rec.mat = mat.get();
      return true;
    }
  }
  return false;
}

class HitableList : public Hitable {
 public:
  HitableList() {};
  
  HitableList& add(Hitable *h) {
    list.push_back(std::unique_ptr<Hitable>(h));
    return *this;
  }

  bool hit(const Ray& r, float t_min, float t_max, 
           HitRecord& rec) const override;

 private:
  std::vector<std::unique_ptr<Hitable>> list;
};

bool HitableList::hit(const Ray& r, float t_min, float t_max, 
                      HitRecord& rec) const {
  bool hitted = false;
  float closest_t = t_max;
  for (const auto& h : list) {
    HitRecord r_temp;
    if (h->hit(r, t_min, closest_t, r_temp)) {
      hitted = true;
      rec = r_temp;
      closest_t = rec.t;
    }
  }
  return hitted;
}

class Camera {
 public:
  Camera(const vec3& lookfrom, const vec3& lookat, const vec3& vup, 
         float vfov, float aspect, float apperture, float focus_dist);

  Ray getRay(float s, float t) const;

 private:
  vec3 lower_left, horizontal, vertical, origin, u, v, w;
  float lens_radius;
};

Camera::Camera(const vec3& lookfrom, const vec3& lookat, const vec3& vup, 
               float vfov, float aspect, float apperture, float focus_dist) {
  lens_radius = apperture / 2;
  float theta = vfov * M_PI / 180;
  float half_h = tan(theta / 2);
  float half_w = aspect * half_h;
  origin = lookfrom;
  w = (lookfrom - lookat).unitVector();
  u = cross(vup, w).unitVector();
  v = cross(w, u);
  lower_left = origin - half_w * focus_dist * u - half_h * focus_dist * v 
               - focus_dist * w;
  horizontal = 2 * half_w * focus_dist * u;
  vertical = 2 * half_h * focus_dist * v;
}

Ray Camera::getRay(float s, float t) const {
  vec3 rd = lens_radius * randomPointSphere();
  vec3 offset = u * rd.x() + v * rd.y();
  return Ray(origin + offset, 
             lower_left + s * horizontal + t * vertical - origin - offset);
}

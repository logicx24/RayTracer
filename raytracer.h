// raytracer.h
//

#ifndef LZZ_raytracer_h
#define LZZ_raytracer_h
#include <stdio.h>
#include <math.h>
#include "lodepng.h"

#include <vector>
#include <iostream>
#define LZZ_INLINE inline
class Vector
{
public:
  float x;
  float y;
  float z;
  float magnitude;
  Vector (float _x, float _y, float _z);
  Vector * operator + (Vector * vec);
  Vector * operator - (Vector * vec);
  float operator * (Vector * vec);
  Vector * operator * (float scale);
  void normalize ();
};
class Point
{
public:
  float x;
  float y;
  float z;
  Point (float _x, float _y, float _z);
  Vector * operator - (Point * p);
};
class Ray
{
public:
  Point * pos;
  Vector * dir;
  float t_min;
  float t_max;
  Ray (Point * p, Vector * v);
};
class Matrix
{
public:
  float ((mat) [4]) [4];
  Matrix (int kind, float theta);
  Matrix (int kind);
};
class Color
{
public:
  float r;
  float g;
  float b;
  Color (float _r = 0.0f, float _g = 0.0f, float _b = 0.0f);
};
class Sample
{
public:
  int x;
  int y;
  Sample (int _x, int _y);
};
class Sampler
{
public:
  int width;
  int height;
  int x;
  int y;
  Sampler (int _width, int _height);
  bool generateSample (Sample * sample);
};
class Film
{
public:
  int width;
  int height;
  Color * film;
  Film (int _width, int _height);
  void writeToFilm (Sample * sample, Color * color);
  void writeFile ();
};
class Scene
{
public:
  Point * eye;
  Point * ul;
  Point * ur;
  Point * ll;
  Point * lr;
  int width;
  int height;
  Film * film;
  Scene (Point * _eye, Point * _ul, Point * _ur, Point * _ll, Point * _lr, Film * _film, int _width, int _height);
  void writeToFilm ();
  void renderLoop ();
};
class Sphere
{
public:
  Point * center;
  float radius;
  Sphere (Point * _center, float _radius);
  bool intersects (Ray * ray);
  Point intersectsAt (Ray * ray);
};
void testFilm ();
int main (int argc, char * (argv) []);
#undef LZZ_INLINE
#endif

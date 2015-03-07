#include <stdio.h>
#include <math.h>
#include "lodepng.h"

#include <vector>
#include <iostream>


#define MAT_ROTATE_X 101
#define MAT_ROTATE_Y 102
#define MAT_ROTATE_Z 103

using namespace std;

class Vector {
public:
	float x, y, z, magnitude;

	Vector(float _x, float _y, float _z) {
		x = _x;
		y = _y;
		z = _z;
		magnitude = sqrt(x * x + y * y + z * z);
	}

	Vector* operator+(Vector *vec) {
		return new Vector(x + vec->x, y + vec->y, z + vec->z);
	}

	Vector* operator-(Vector *vec) {
		return new Vector(x - vec->x, y - vec->y, z - vec->z);
	}

	float operator*(Vector *vec) {
		return x * vec->x + y * vec->y + z * vec->z;
	}

	Vector* operator*(float scale) {
		return new Vector(x * scale, y * scale, z * scale);
	}

	void normalize() {
		x = x/magnitude;
		y = y/magnitude;
		z = z/magnitude;
	}
};

class Point {
public:
	float x, y, z;

	Point(float _x, float _y, float _z) {
		x = _x;
		y = _y;
		z = _z;
	}

	Vector* operator-(Point *p) {
		return new Vector(x - p->x, y - p->y, z - p->z);
	}
};

class Ray {
public:
	Point *pos;
	Vector *dir;
	float t_min, t_max;
	
	Ray(Point *p, Vector *v) {
		// Possibly make copies?
		pos = p;
		dir = v;
		t_min = 1;
		t_max = 10000000.0f;
	}
};

class Matrix {
public:
	float mat[4][4];

	Matrix(int kind, float theta) {
		for (int x = 0; x < 4; x++) {
			for (int y = 0; y < 4; y++) {
				mat[y][x] = 0.0f;
			}
		}
		switch (kind) {
			case MAT_ROTATE_X:
				mat[0][0] = 1.0f;
				mat[1][1] = cos(theta);
				mat[2][1] = -sin(theta);
				mat[2][1] = sin(theta);
				mat[2][2] = cos(theta);
				mat[3][3] = 1.0f;
				break;
			case MAT_ROTATE_Y:
				mat[0][0] = cos(theta);
				mat[0][2] = -sin(theta);
				mat[1][1] = 1.0f;
				mat[2][0] = sin(theta);
				mat[2][1] = cos(theta);
				mat[3][3] = 1.0f;
				break;
			case MAT_ROTATE_Z:
				mat[0][0] = cos(theta);
				mat[0][2] = -sin(theta);
				mat[1][0] = sin(theta);
				mat[1][1] = cos(theta);
				mat[2][2] = 1.0f;
				mat[3][3] = 1.0f;
				break;
		}
	}

	Matrix(int kind) {

	}
};

class Color {
public:
    float r, g, b;
    Color(float _r=0.0f, float _g=0.0f, float _b=0.0f) {
    	r = _r;
    	g = _g;
    	b = _b;
    }
};

class Sample {
public:
    int x, y;
    Sample(int _x, int _y) {
    	x = _x;
    	y = _y;
    }
};

class Sampler {
public:
    int width, height;
    int x, y;

    Sampler(int _width, int _height) {
    	width = _width;
    	height = _height;
        x = 0;
        y = 0;
    }

    bool generateSample(Sample *sample) {
        if (y  == height) {
        	return false;
        }
    	sample->x = x;
    	sample->y = y;
        if (x == width - 1) {
        	x = 0;
            y++;
        } else {
        	x++;
        }
        return true;
    }

};

class Film {
public:
    int width, height;
    Color *film;
    Film(int _width, int _height) {
    	width = _width;
    	height = _height;
    	film = new Color[width * height];
    }

    void writeToFilm(Sample *sample, Color *color) {
        film[sample->y * width + sample->x].r = color->r;
        film[sample->y * width + sample->x].g = color->g;
        film[sample->y * width + sample->x].b = color->b;
    }

    void writeFile() {
        vector<unsigned char> image;
        int area = width * height;
        
        for (int i = 0; i < area; i++) {
        	image.push_back(0);
        	image.push_back(0);
        	image.push_back(0);
        	image.push_back(255);
        }
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
            	int position = 4 * (y * width + x);
            	image.at(position)     = film[y * width + x].r;
            	image.at(position + 1) = film[y * width + x].g;
            	image.at(position + 2) = film[y * width + x].b;
            	image.at(position + 3) = 255; 
            }
        }
        unsigned error = lodepng::encode("img.png", image, width, height);
        if (error) {
        	cout << "lodepng error" << endl;
        }
    }
};

class Sphere {
public:
	Point *center;
	float radius;

	Sphere(Point *_center, float _radius) {
		center = _center;
		radius = _radius;
	}

	bool intersects(Ray *ray) {
		// (d * (e - c)) - (d * d)((e-c) * (e-c))
        float discriminant = pow(ray->dir * (ray->pos - center), 2) - 
					 (ray->dir * ray->dir)*(((ray->pos - center)*(ray->pos - center))
					 - pow(radius, 2));
		return discriminant >= 0;
	}

	Point intersectsAt(Ray *ray) {
		float discriminant = pow(ray->dir * (ray->pos - center), 2) - 
					 (ray->dir * ray->dir)*(((ray->pos - center)*(ray->pos - center))
					 - pow(radius, 2));
		float soln1 = (((ray->dir) * -1.0f) * (ray->pos - center) + sqrt(discriminant))/(ray->dir * ray->dir);
		float soln2 = (((ray->dir) * -1.0f)* (ray->pos - center) - sqrt(discriminant))/(ray->dir * ray->dir);
		//Figure out which solution is correct.
		return new Point(ray->pos->x + soln1*ray->dir->x, ray->pos->y + soln1*ray->dir->y, ray->pos->z + soln1*ray->dir->z);

	}

};

class Scene {
public:
	Point *eye, *ul, *ur, *ll, *lr;
	int width, height;
	Film *film;

	Scene(Point *_eye, Point *_ul, Point *_ur, Point *_ll, Point *_lr, 
		  Film *_film, int _width, int _height) {
		eye = _eye;
		ul = _ul;
		ur = _ur;
		ll = _ll;
		lr = _lr;
		film = _film;
		width = _width;
		height = _height;
	}
    
    void testSphere() {
    	Sampler *sampler = new Sampler(width, height);
    	Sample *sample = new Sample(0, 0);
    	Color *black = new Color(0, 0, 0);
    	Color *red = new Color(255, 0, 0);
    	Point *center = new Point(50, 50, 10);
    	Sphere *sphere = new Sphere(center, 10);
        while (sampler->generateSample(sample)) {
        	Vector *testRayDir = new Vector(sample->x - eye->x, sample->y - eye->y, -eye->z);
        	Point *testRayPoint = new Point(eye->x, eye->y, eye->z);
        	Ray *testRay = new Ray(testRayPoint, testRayDir);
            if (sphere->intersects(testRay)) {
            	film->writeToFilm(sample, red);
            	//cout << "RED" << endl;
            } else {
            	film->writeToFilm(sample, black);
            	//cout << "BLACK" << endl;
            }
            delete testRayDir;
            delete testRayPoint;
            delete testRay;
        }
        film->writeFile();
    }

    void renderLoop() {
        
    }

};



void testFilm() {
	int height = 100;
	int width = 100;
	Film *film = new Film(width, height);
	Sampler *sampler = new Sampler(width, height);
	Sample *sample = new Sample(0, 0);
	Color *color = new Color(0, 0, 0);
    while (sampler->generateSample(sample)) {
    	//cout << "x: " << sample->x << endl;
    	//cout << "y: " << sample->y << endl;
        film->writeToFilm(sample, color);
    }
    film->writeFile();
}


int main(int argc, char* argv[]) {
	printf("Ray Tracer1!\n");
	Point *eye = new Point(50, 50, -10);
	Point *ul = new Point (0, 99, 0);
	Point *ur = new Point (99, 99, 0);
	Point *ll = new Point (0, 0, 0);
	Point *lr = new Point (99, 0, 0);
	Film *film = new Film(100, 100);
	Scene *scene = new Scene(eye, ul, ur, ll, lr, film, 100, 100);
	scene->testSphere();
	return 0;
}


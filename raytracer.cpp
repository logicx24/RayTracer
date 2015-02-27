#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAT_ROTATE_X 101
#define MAT_ROTATE_Y 102
#define MAT_ROTATE_Z 103

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


class Film {
public:
    int width, height;
    Color **film;
    Film(int _width, int _height) {
    	width = _width;
    	height = _height;
    	film = new Color*[width * height];
    }

    void commit(Sample& sample, Color& color) {

    }
};

class Scene {
public:
	Point *eye, *ul, *ur, *ll, *lr;
	int width, height;

	Scene(Point *_eye, Point *_ul, Point *_ur, Point *_ll, Point *_lr, 
		  int _width, int _height) {
		eye = _eye;
		ul = _ul;
		ur = _ur;
		ll = _ll;
		lr = _lr;
		width = _width;
		height = _height;
	}
};


int main(int argc, char* argv[]) {
	printf("Ray Tracer!\n");
	return 0;
}

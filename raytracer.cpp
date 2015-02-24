#include <stdio.h>
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
		switch (kind) {
			case MAT_ROTATE_X:
				float _mat[4] = {{1, 0, 		   0,  			0}, 
					   {0, cos(theta), -sin(theta), 0},
					   {0, sin(theta), cos(theta),  0},
					   {0, 0, 		   0, 			1}};
				break;
			case MAT_ROTATE_Y:
				float _mat[4] = {{cos(theta), 0, -sin(theta), 0}, 
					   {0, 			1, 0, 			0},
					   {sin(theta), 0, cos(theta),  0},
					   {0, 			0, 0, 			1}};
				break;
			case MAT_ROTATE_Z:
				float _mat[4][4] = {{cos(theta), -sin(theta), 0, 0}, 
					   {sin(theta), cos(theta),  0, 0},
					   {0, 			0, 			 1, 0},
					   {0, 			0, 			 0, 1}};
				break;
		}
		mat = _mat;
	}

	Matrix(int kind) {

	}
};


int main(int argc, char* argv[]) {
	printf("Ray Tracer!\n");
	return 0;
}
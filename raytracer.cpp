#include <stdio.h>
#include <math.h>

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


int main(int argc, char* argv[]) {
	printf("Ray Tracer!\n");
	return 0;
}
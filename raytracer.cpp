#include <stdio.h>
#include <math.h>
#include "lodepng.h"
#include "algebra3.h"

#include <vector>
#include <iostream>


using namespace std;

class Color {
public:
    float r, g, b;
    Color(float _r=0.0f, float _g=0.0f, float _b=0.0f) {
    	r = _r;
    	g = _g;
    	b = _b;
    }
};


class Ray {
public:
	vec3 pos;
	vec3 vec;
	float t_min, t_max;
	
	Ray(vec3 &p, vec3 &v) {
		// Possibly make copies?
		pos = p;
		vec = v;
		t_min = 1;
		t_max = 10000000.0f;
	}
};



class Sample {
public:
    int x, y;

    Sample(int _x = 0, int _y = 0) {
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
        if (y == height) {
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

    Film(){}

    Film(int _width, int _height) {
    	width = _width;
    	height = _height;
    	film = new Color[width * height];
    }

    void writeToFilm(Sample &sample, Color &color) {
        film[sample.y * width + sample.x].r = color.r;
        film[sample.y * width + sample.x].g = color.g;
        film[sample.y * width + sample.x].b = color.b;
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
	vec3 center;
	float radius;

	Sphere(vec3 &_center, float _radius) {
		center = _center;
		radius = _radius;
	}

	bool intersects(Ray &ray) {
		// (d * (e - c)) - (d * d)((e-c) * (e-c))
		vec3 p_c = ray.pos - center;
        //p_c.normalize();
        float discriminant = pow(ray.vec * p_c, 2) - 
        ((ray.vec * ray.vec) * ((p_c * p_c) - pow(radius, 2)));
        cout << "d: " << discriminant << endl;
		return discriminant >=  0;
	}
};

class Camera {
public:

	vec3 eye, ul, ur, ll, lr;
    int width, height;

	Camera(){}

	Camera(vec3 &_eye, vec3 &_ul, vec3 &_ur, vec3 &_ll, vec3 &_lr, 
	       int _width, int _height) {
	    eye = _eye;
		ul = _ul;
		ur = _ur;
		ll = _ll;
		lr = _lr;
		width = _width;
		height = _height;
	}

	Ray generateRay(Sample &sample) {
        float u = sample.x / float(width);
        float v = sample.y / float(height);
        // cout << sample.x << " : " << sample.y << endl;
        // cout << width << " : " << height << endl;
        // cout << u << " : " << v << endl;
        vec3 c = (ll*(1-u) + lr*u)*(1-v) + (ul*(1-u) + ur*u)*v;
        //cout << "1. " << c[VX] << " : " << c[VY] << " : " << c[VZ] << endl;
        vec3 testRayVec = c - eye;
        testRayVec = testRayVec.normalize();
        //cout << "2. " << c[VX] << " : " << c[VY] << " : " << c[VZ] << endl;
	    return Ray(eye, testRayVec); 
	}
};


class Scene {
public:
	Camera cam;
	Film film;
	vector<Sphere> objects;
	int width, height;
	
	Scene(Camera &_cam, Film &_film, int _width, int _height)  {
		cam = _cam;
		film = _film;
		width = _width;
		height = _height;
	}

	void addObject(Sphere &s) {
 		objects.push_back(s);
	}
    

    void renderLoop() {
    	Sampler sampler = Sampler(width, height);
    	Sample sample = Sample();
    	// Ray testRay = Ray(rayPos, rayVec);
    	while (sampler.generateSample(&sample)) {
    		//cout << "x: " << sample.x << " y: " << sample.y << endl;
    		Ray testRay = cam.generateRay(sample);
    		Color c = traceRay(testRay);
    		//cout << "r: " << c.r << " g: " << c.g << " b: " << c.b << endl;
    		film.writeToFilm(sample, c);
    	}

        film.writeFile();
    }

    Color traceRay(Ray &ray) {
    	Color black = Color(0, 0, 0);
    	Color red = Color(255, 0, 0);
    	//cout << objects.size() << endl;
		for (int i = 0; i < objects.size(); i++) {
			if (objects[i].intersects(ray)) {
				return Color(255, 0, 0);
				// Recurse from objects[i]->intersectsAt(testRay) and bounce dir!
				// Shade!

			} else {
				return Color(0, 0, 0);
			}
		}
		return Color(0, 0, 0);
    }

};




int main(int argc, char* argv[]) {
	printf("Ray Tracer!\n");
	int width = 100;
	int height = 100;
	vec3 eye = vec3(0, 0, 0);
	vec3 ul = vec3(-1, 1, -1);
	vec3 ur = vec3(1, 1, -1);
	vec3 ll = vec3(-1, -1, -1);
	vec3 lr = vec3(1, -1, -1);
	Film film = Film(width, height);
    Camera cam = Camera(eye, ul, ur, ll, lr, width, height);
	Scene scene = Scene(cam, film, width, height);


	vec3 sphereCenter = vec3(0, 0, -2);
	Sphere sphere = Sphere(sphereCenter, 1);
	scene.addObject(sphere);


	scene.renderLoop();
	return 0;
}


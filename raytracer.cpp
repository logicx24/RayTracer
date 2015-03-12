#include <stdio.h>
#include <math.h>
#include "lodepng.h"
#include "algebra3.h"

#include <vector>
#include <iostream>


using namespace std;
void printVec3(vec3 &vec) {
	cout << "x: " << vec[VX] << " y: " << vec[VY] << " Z: " << vec[VZ] << endl;
}

float distance(vec3 &p1, vec3 &p2) {
	return sqrt(p1[VX]*p2[VX] + p1[VY]*p2[VY] + p1[VZ]*p2[VZ]);
}

vec3 multiplyComponents(vec3 v1, vec3 v2) {
	return vec3(v1[VX] * v2[VX], v1[VY] * v2[VY], v1[VZ] * v2[VZ]);
}

class Material {
public:
    vec3 ka, kd, ks, kr;
    float sp;
    
    Material() {
    	ka = vec3(0, 0, 0);
    	kd = vec3(0, 0, 0);
    	ks = vec3(0, 0, 0);
    	kr = 0.0f;
    	sp = 0.0f;
    }

    Material(vec3 &_ka, vec3 &_kd, vec3 &_ks, vec3 &_kr, float _sp) {
    	ka = _ka;
    	kd = _kd;
    	ks = _ks;
    	kr = _kr;
    	sp = _sp;
    }
};
class Ray {
public:
	vec3 pos;
	vec3 vec;
	float t_min, t_max;
	
	Ray(vec3 &p, vec3 &v) {
		pos = p;
		vec = v;
		t_min = 0.001f;
		t_max = 10000000.0f;
	}

    vec3 pointAt(float t) {
    	return pos + t * vec;
    }
};

class Light {
public:
    vec3 color;
    vec3 pos;
    bool pointLight;

    Light(vec3 &_color, vec3 &_pos, int falloff, bool _pointLight) {
    	color = _color;
    	pos = _pos;
    	pointLight = _pointLight;
    }

    vec3 lightVec(vec3 &point) {
        if (pointLight) {
        	return -(point - pos);
        } else {
        	return -pos;
        }
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
    vector<vec3> film;

    Film(){}

    Film(int _width, int _height) {
    	width = _width;
    	height = _height;
        for (int i = 0; i < width * height; i++) {
        	film.push_back(vec3(0, 0, 0));
        }
    }

    void writeToFilm(Sample &sample, vec3 &color) {
    	if (color[RED] > 255) {color[RED] = 255;}
    	if (color[GREEN] > 255) {color[GREEN] = 255;}
    	if (color[BLUE] > 255) {color[BLUE] = 255;}
        film.at(sample.y * width + sample.x)[RED] = color[RED];
        film.at(sample.y * width + sample.x)[GREEN] = color[GREEN];
        film.at(sample.y * width + sample.x)[BLUE] = color[BLUE];
    }

    void writeFile(int pixWidth, int pixHeight) {
        if (pixWidth % width != 0 || pixHeight % height != 0) {
        	cerr << "Image will not fit nicely into output dimensions." << endl;
        	exit(0);
        }
        int widthScale = pixWidth / width;
        int heightScale = pixHeight / height;
        vector<unsigned char> image;
        int area = pixWidth * pixHeight;
        
        for (int i = 0; i < area; i++) {
        	image.push_back(0);
        	image.push_back(0);
        	image.push_back(0);
        	image.push_back(255);
        }
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
            	int position = 4 * (y * width + x);
            	image.at(position)     = film.at(y * width + x)[RED];
            	image.at(position + 1) = film.at(y * width + x)[GREEN];
            	image.at(position + 2) = film.at(y * width + x)[BLUE];
            	image.at(position + 3) = 255; 
            }
        }
        unsigned error = lodepng::encode("img.png", image, width, height);
        if (error) {
        	cout << "lodepng error" << endl;
        }
    }
};

class Polygon {
public:
    Material mat;
    Polygon(){};
    virtual bool intersects(Ray &ray) = 0;
    virtual float intersection(Ray &ray, float eps) = 0;
    virtual vec3 normalAt(vec3 &point) = 0;
};

class Sphere : public Polygon {
public:
	vec3 center;
	float radius;

	Sphere(vec3 &_center, float _radius, Material &_mat) {
		center = _center;
		radius = _radius;
		mat = _mat;
	}

	bool intersects(Ray &ray) {
		vec3 p_c = ray.pos - center;
        float discriminant = pow(ray.vec * p_c, 2) - 
        ((ray.vec * ray.vec) * ((p_c * p_c) - pow(radius, 2)));
		return discriminant >= 0.0f;
	}
    
	float intersection(Ray &ray, float eps) {
		vec3 p_c = ray.pos - center;
        vec3 d = ray.vec;
        float discriminant = pow(d * p_c, 2) - ((d * d) * ((p_c * p_c) - pow(radius, 2)));
        if (discriminant < 0) {
        	return MAXFLOAT;
        }
		float t1 = (-d*(p_c) + sqrt(discriminant)) / (d * d);
		float t2 = (-d*(p_c) - sqrt(discriminant)) / (d * d);
        if (t1 < eps && t2 < eps) {
        	return MAXFLOAT;
        }
	    if (t1 > 0) {
        	return t2;
        } else {
        	return t1;
        }
	}

    vec3 normalAt(vec3 &point) {
    	vec3 normal = (point - center) / radius;
    	//printVec3(normal);
    	return normal;
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
        vec3 c = (ll*(1-u) + lr*u)*(1-v) + (ul*(1-u) + ur*u)*v;
        vec3 testRayVec = c - eye;
        testRayVec = testRayVec.normalize();
	    return Ray(eye, testRayVec); 
	}
};


class Scene {
public:
	Camera cam;
	Film film;
	vector<Polygon*> objects;
	vector<Light*> lights;

	int width, height;
	int maxDepth;
	float eps;

	Scene(Camera &_cam, Film &_film, float _eps, int _width, int _height)  {
		cam = _cam;
		film = _film;
		width = _width;
		height = _height;
		maxDepth = 4;
		eps = _eps;
	}

	void addObject(Polygon *p) {
 		objects.push_back(p);
	}

    void addLight(Light *l) {
    	lights.push_back(l);
    }

    void writeFile(int pixWidth, int pixHeight) {
    	film.writeFile(pixWidth, pixHeight);
    }

    void render() {
    	Sampler sampler = Sampler(width, height);
    	Sample sample = Sample();
    	while (sampler.generateSample(&sample)) {
    		Ray testRay = cam.generateRay(sample);
    		vec3 c = traceRay(testRay, maxDepth);
    		film.writeToFilm(sample, c);
    	}

    }

    vec3 traceRay(Ray &ray, int depth) {
    	float min_t = MAXFLOAT;
    	int closeIndex = 0;
		for (int i = 0; i < objects.size(); i++) {
			float t = objects[i]->intersection(ray, eps);
            if (t != MAXFLOAT) {
                if (t < min_t) {
                	min_t = t;
                	closeIndex = i;
                }
		    }
		}
		if (min_t != MAXFLOAT) {
			Polygon *hitObject = objects[closeIndex];
		    vec3 intersection = ray.pointAt(min_t);
            vec3 normal = hitObject->normalAt(intersection);
            vec3 reflectionVec = (ray.vec - 2*(ray.vec * normal)*normal);
            reflectionVec.normalize();
            Ray reflected = Ray(intersection, reflectionVec);
            if (depth == 0) {
            	return phongShade(ray, hitObject, normal, intersection);
            } else {
        	    return phongShade(ray, hitObject, normal, intersection)
        	       + multiplyComponents(traceRay(reflected, depth - 1), hitObject->mat.kr);
		    }
		} else {
		    return vec3(0, 0, 0);

        }
    }

    bool traceShadowRay(Ray &ray, Polygon *source) {
        for (int i = 0; i < objects.size(); i++) {
            if (objects[i]->intersects(ray) && objects[i] != source) {
            	return true;
            }
        }
        return false;
    }

    vec3 phongShade(Ray &ray, Polygon *polygon, vec3 &normal, vec3 &point) {
        vec3 total_intensity = vec3(0, 0, 0);
        for (int i = 0; i < lights.size(); i++) {
        	vec3 intensity = lights[i]->color;
        	vec3 shadowRayVec = (lights[i]->pos - point).normalize();
        	Ray shadowRay = Ray(point, shadowRayVec);
            if (traceShadowRay(shadowRay, polygon)) {
            	continue;
            }
        	vec3 l = lights[i]->lightVec(point);
            vec3 v = -ray.vec;
            vec3 r = -l + ((2*(l * normal)) * normal);
            l.normalize();
        	v.normalize();
        	r.normalize();
            vec3 coeffs = polygon->mat.ka + polygon->mat.kd * MAX(l * normal, 0.0f) 
                        + polygon->mat.ks * pow(MAX(r * v, 0.0f), polygon->mat.sp);
            total_intensity += multiplyComponents(intensity, coeffs);
        }
        return total_intensity; 
    }

};




int main(int argc, char* argv[]) {
	printf("Ray Tracer!\n");
	int width = 1000;
	int height = 1000;
	float eps = 0.01f;
	vec3 eye = vec3(0, 0, 0);
	vec3 ul = vec3(-1, 1, -1);
	vec3 ur = vec3(1, 1, -1);
	vec3 ll = vec3(-1, -1, -1);
	vec3 lr = vec3(1, -1, -1);
	Film film = Film(width, height);
    Camera cam = Camera(eye, ul, ur, ll, lr, width, height);
	Scene scene = Scene(cam, film, eps, width, height);
    

    vec3 ka = vec3(0.05, 0.05, 0.05);
    vec3 kd = vec3(1, 1, 1);
    vec3 ks = vec3(1, 1, 1);
    vec3 kr = vec3(0.94, 0.94, 0.94);
    float sp = 50.0f;

    vec3 ka1 = vec3(0.5, 0.04, 0.15);
    vec3 kd1 = vec3(1, 1, 1);
    vec3 ks1 = vec3(1, 1, 1);
    vec3 kr1 = vec3(0.88, 0.88, 0.88);
    float sp1 = 45.0f;

    Material sphereMat = Material(ka, kd, ks, kr, sp);
    Material sphereMat1 = Material(ka1, kd1, ks1, kr1, sp1);
	
	vec3 sphereCenter1 = vec3(-5, -5, -10);
	Sphere sphere1 = Sphere(sphereCenter1, 3.0f, sphereMat);
	scene.addObject(&sphere1);
	
	vec3 sphereCenter2 = vec3(0, 2, -5);
	Sphere sphere2 = Sphere(sphereCenter2, 2.0f, sphereMat1);
	scene.addObject(&sphere2);

	vec3 sphereCenter3 = vec3(5, -5, -10);
	Sphere sphere3 = Sphere(sphereCenter3, 2.0f, sphereMat);
	scene.addObject(&sphere3);

    vec3 lightColor1 = vec3(122, 130, 20);
    vec3 lightPos1 = vec3(1, 1, -1);
	Light dl1 = Light(lightColor1, lightPos1, 0, false);
	scene.addLight(&dl1);

    vec3 lightColor2 = vec3(40, 130, 250);
    vec3 lightPos2 = vec3(1, 1, 1);
	Light pl2 = Light(lightColor2, lightPos2, 0, true);
	//scene.addLight(&pl2);

	scene.render();
	scene.writeFile(width, height);
	
	return 0;
}


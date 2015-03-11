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
class Color {
public:
    float r, g, b;
    Color(float _r=0.0f, float _g=0.0f, float _b=0.0f) {
    	r = _r;
    	g = _g;
    	b = _b;
    }

    Color add(Color color) {
    	return Color(r + color.r, g + color.g, b + color.b);
    }
    Color scale(float scale) {
    	return Color(scale*r, scale*g, scale*b);
    }
};

class Material {
public:
    vec3 ka, kd, ks;
    float sp, kr;
    
    Material() {
    	ka = vec3(0, 0, 0);
    	kd = vec3(0, 0, 0);
    	ks = vec3(0, 0, 0);
    	kr = 0.0f;
    	sp = 0.0f;
    }

    Material(vec3 &_ka, vec3 &_kd, vec3 &_ks, float _kr, float _sp) {
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
		// Possibly make copies?
		pos = p;
		vec = v;
		t_min = 1;
		t_max = 10000000.0f;
	}
};

class Light {
public:
    Color color;
    vec3 pos;
    bool pointLight;

    Light(Color &_color, vec3 &_pos, bool _pointLight) {
    	color = _color;
    	pos = _pos;
    	pointLight = _pointLight;
    }

    vec3 lightVec(vec3 &point) {
        if (pointLight) {
        	return (point - pos);
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
    Color *film;

    Film(){}

    Film(int _width, int _height) {
    	width = _width;
    	height = _height;
    	film = new Color[width * height];
    }

    void writeToFilm(Sample &sample, Color &color) {
    	if (color.r > 255) {color.r = 255;}
    	if (color.g > 255) {color.g = 255;}
    	if (color.b > 255) {color.b = 255;}
        film[sample.y * width + sample.x].r = color.r;
        film[sample.y * width + sample.x].g = color.g;
        film[sample.y * width + sample.x].b = color.b;
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

class Polygon {
public:
    Material mat;
    Polygon(){};
    virtual bool intersects(Ray &ray) = 0;
    virtual float intersectsAt(Ray &ray) = 0;
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
		// discriminant = (d * (e - c)) - (d * d)((e-c) * (e-c))
		// d = ray.pos
		vec3 p_c = ray.pos - center;
        float discriminant = pow(ray.vec * p_c, 2) - 
        ((ray.vec * ray.vec) * ((p_c * p_c) - pow(radius, 2)));
		return discriminant >= 0.0f;
	}
    
	float intersectsAt(Ray &ray) {
		vec3 p_c = ray.pos - center;
        vec3 d = ray.vec;
        float discriminant = pow(d * p_c, 2) - ((d * d) * ((p_c * p_c) - pow(radius, 2)));
		float t1 = (-d*(p_c) + sqrt(discriminant)) / (d * d);
		float t2 = (-d*(p_c) - sqrt(discriminant)) / (d * d);
		/*
		vec3 p1 = ray.pos + ray.vec * t1;
		vec3 p2 = ray.pos + ray.vec * t2;
		float d1 = distance(ray.pos, p1);
		float d2 = distance(ray.pos, p2);
        if (d1 < d2) {
        	return ray.pos + p1;
        } else {
        	return ray.pos + p2;
        }*/
        return MIN(t1, t2);
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
	int maxDepth = 5;

	Scene(Camera &_cam, Film &_film, int _width, int _height)  {
		cam = _cam;
		film = _film;
		width = _width;
		height = _height;
	}

	void addObject(Polygon *p) {
 		objects.push_back(p);
	}

    void addLight(Light *l) {
    	lights.push_back(l);
    }

    void render() {
    	Sampler sampler = Sampler(width, height);
    	Sample sample = Sample();
    	while (sampler.generateSample(&sample)) {
    		Ray testRay = cam.generateRay(sample);
    		Color c = traceRay(testRay);
    		film.writeToFilm(sample, c);
    	}

    }

    void writeFile(int pixWidth, int pixHeight) {
    	film.writeFile(pixWidth, pixHeight);
    }

    Color traceRay(Ray &ray) {
    	float min_t = MAXFLOAT;
    	int closeIndex = 0;
		for (int i = 0; i < objects.size(); i++) {
			if (objects[i]->intersects(ray)) {
		        float t = objects[i]->intersectsAt(ray);
                if (t < min_t) {
                	min_t = t;
                	closeIndex = i;
                }
		    }
		}
		if (min_t != MAXFLOAT) {
		    vec3 intersection = ray.pos + min_t * ray.vec;
            vec3 normal = objects[closeIndex]->normalAt(intersection);
            vec3 reflectionVec = (ray.vec - 2*(ray.vec * normal)*normal);
            Ray reflected = Ray(intersection, reflectionVec);
        	return phongShade(ray, objects[closeIndex], normal, intersection).add(traceReflection(reflected, maxDepth, objects[closeIndex]).scale(objects[closeIndex]->mat.kr));
		} else {
		    return Color(0, 0, 0);

        }
    }

    Color traceReflection(Ray &ray, int depth, Polygon *original) {
    	if (depth == 0) {
    		return Color(0, 0, 0);
    	}
    	float min_t = MAXFLOAT;
    	int closeIndex = 0;
		for (int i = 0; i < objects.size(); i++) {
			if (objects[i]->intersects(ray) && objects[i] != original) {
		        float t = objects[i]->intersectsAt(ray);
                if (t < min_t) {
                	min_t = t;
                	closeIndex = i;
                }
		    }
		}
		if (min_t != MAXFLOAT) {
		    vec3 intersection = ray.pos + min_t * ray.vec;
            vec3 normal = objects[closeIndex]->normalAt(intersection);
            vec3 reflectionVec = (ray.vec - 2*(ray.vec * normal)*normal);
            Ray reflected = Ray(intersection, reflectionVec);
        	return phongShade(ray, objects[closeIndex], normal, intersection).add(traceReflection(reflected, depth-1, objects[closeIndex]).scale(objects[closeIndex]->mat.kr));
		} else {
		    return Color(0, 0, 0);
        }
    }

    bool traceShadowRay(Ray &ray, Polygon *p) {
    	//cout << objects.size() << endl;
        for (int i = 0; i < objects.size(); i++) {
            if (objects[i]->intersects(ray) && objects[i] != p) {
            	//cout << "HALLO" << endl;
            	return true;
            }
        }
        return false;
    }

    Color phongShade(Ray &ray, Polygon *polygon, vec3 &normal, vec3 &point) {
        Color total_intensity = Color(0, 0, 0);
        for (int i = 0; i < lights.size(); i++) {
        	Color intensity = lights[i]->color;
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
            vec3 coeffs = polygon->mat.ka + polygon->mat.kd * MAX(l * normal, -0.0f) 
                        + polygon->mat.ks * pow(MAX(r * v, 0.0f), polygon->mat.sp);
            total_intensity.r += intensity.r * coeffs[RED];
            total_intensity.g += intensity.g * coeffs[GREEN];
            total_intensity.b += intensity.b * coeffs[BLUE];
        }
        return total_intensity; 
    }

};




int main(int argc, char* argv[]) {
	printf("Ray Tracer!\n");
	int width = 1000;
	int height = 1000;
	vec3 eye = vec3(0, 0, 0);
	vec3 ul = vec3(-1, 1, -1);
	vec3 ur = vec3(1, 1, -1);
	vec3 ll = vec3(-1, -1, -1);
	vec3 lr = vec3(1, -1, -1);
	Film film = Film(width, height);
    Camera cam = Camera(eye, ul, ur, ll, lr, width, height);
	Scene scene = Scene(cam, film, width, height);
    

    vec3 ka = vec3(0.05, 0.05, 0.05);
    vec3 kd = vec3(1, 1, 1);
    vec3 ks = vec3(1, 1, 1);
    float kr = 0.94f;
    float sp = 50.0f;

    vec3 ka1 = vec3(0.5, 0.04, 0.15);
    vec3 kd1 = vec3(1, 1, 1);
    vec3 ks1 = vec3(1, 1, 1);
    float kr1 = 0.88f;
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

    Color lightColor1 = Color(122, 130, 20);
    vec3 lightPos1 = vec3(0, 0, -1);
	Light dl1 = Light(lightColor1, lightPos1, false);
	scene.addLight(&dl1);

	scene.render();
	scene.writeFile(width, height);
	
	return 0;
}


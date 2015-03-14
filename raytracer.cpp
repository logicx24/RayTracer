#include <stdio.h>
#include <math.h>
#include "lodepng.h"
#include "algebra3.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

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
		t_max = MAXFLOAT;
	}

    vec3 pointAt(float t) {
    	return pos + t * vec;
    }
};

class Light {
public:
    vec3 color;
    vec3 pos;
    int falloff;
    bool pointLight;

    Light(vec3 &_color, vec3 &_pos, int _falloff, bool _pointLight) {
    	color = _color;
    	pos = _pos;
    	falloff = _falloff;
    	pointLight = _pointLight;
    }

    vec3 lightVec(vec3 &point) {
        if (pointLight) {
        	return (pos - point);
        } else {
        	return -pos;
        }
    }

    Ray lightRay(vec3 &point) {
        vec3 rayVec = lightVec(point);
        rayVec.normalize();
        return Ray(point, rayVec);
    }

    vec3 intensityAt(vec3 &point) {
        if (falloff == 0) {
        	return color;
        } else {
        	float d = distance(pos, point);
            return color / pow(d, falloff);
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

            	image.at(position)     = film.at((pixHeight - y - 1) * width + x)[RED];
            	image.at(position + 1) = film.at((pixHeight - y - 1) * width + x)[GREEN];
            	image.at(position + 2) = film.at((pixHeight - y - 1) * width + x)[BLUE];
            	image.at(position + 3) = 255; 
            }
        }
        unsigned error = lodepng::encode("img.png", image,  width, height);
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
    virtual vec3 normal(vec3 &point) = 0;
};

class Triangle : public Polygon {
public:
    vec3 vertex1;
    vec3 vertex2;
    vec3 vertex3;

    vec3 edge1;
    vec3 edge2;

    Triangle(vec3 &_v1, vec3 &_v2, vec3 &_v3, Material &_mat) {
        vertex1 = _v1;
        vertex2 = _v2;
        vertex3 = _v3;
        mat = _mat;

        edge1 = vertex2 - vertex1;
        edge2 = vertex3 - vertex1;
    }

    bool intersects(Ray &ray) {
        // vec3 normal = edge1 ^ edge2; 
        // vec3 w0 = ray.pos - vertex1;
        // float a = -(normal * w0);
        // //float b = w0 * ray.vec;
        // float b = normal * ray.vec;

        // if (fabs(b) < 0.1f) {
        //     if (a == 0.0f) {
        //         return true;
        //     }
        //     return false;
        // }

        // float r = a/b;

        // if (r < 0.0f) {
        //     return false;
        // }

        // vec3 intersect = ray.pointAt(r);
        // float uu, uv, vv, wu, wv, D;

        // uu = edge1 * edge1;
        // uv = edge1 * edge2;
        // vv = edge2 * edge2;
        // vec3 w = intersect - vertex1;
        // wu = w * edge1;
        // wv = w * edge2;
        // D = uv * uv - uu * vv;

        // float s, t;

        // s = (uv * wv - vv * wu) / D;
        // if (s < 0.0f || s > 1.0f) {
        //     return false;
        // }
        // t = (uv * wu - uu * wv) / D;
        // if (t < 0.0f || (s + t) > 1.0f) {
        //     return false;
        // }

        // if (r < 0.1f) {
        //     return false;
        // }

        // return true;
        return intersection(ray, 0.1f) != MAXFLOAT;

    }

    float intersection(Ray &ray, float eps) {
        vec3 normal = edge1 ^ edge2; 
        vec3 w0 = ray.pos - vertex1;
        float a = -(normal * w0);
        float b = normal * ray.vec;
        float r = a/b;

        if (fabs(b) < eps) {
            if (a == 0.0f) {
                return r;
            }
            return MAXFLOAT;
        }

        if (r < 0.0f) {
            return MAXFLOAT;
        }

        vec3 intersect = ray.pointAt(r);
        float uu, uv, vv, wu, wv, D;

        uu = edge1 * edge1;
        uv = edge1 * edge2;
        vv = edge2 * edge2;
        vec3 w = intersect - vertex1;
        wu = w * edge1;
        wv = w * edge2;
        D = uv * uv - uu * vv;

        float s, t;

        s = (uv * wv - vv * wu) / D;
        if (s < 0.0f || s > 1.0f) {
            return MAXFLOAT;
        }
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0f || (s + t) > 1.0f) {
            return MAXFLOAT;
        }

        if (r < eps) {
            return MAXFLOAT;
        }
        return r;
    }

    vec3 normal(vec3 &point) {
        vec3 normal = (edge1 ^ edge2);
        return normal;
    }
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
        // vec3 p_c = ray.pos - center;
        // vec3 d = ray.vec;
        // float discriminant = pow(d * p_c, 2) - ((d * d) * ((p_c * p_c) - pow(radius, 2)));
		//return discriminant >= 0 && intersection(ray, 0.1f) != MAXFLOAT;
	   return intersection(ray, 0.1f) != MAXFLOAT;
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
        vec3 p1 = ray.pointAt(t1);
        vec3 p2 = ray.pointAt(t2);
        float d1 = distance(ray.pos, p1);
        float d2 = distance(ray.pos, p2);
	    if (d1 < d2) {
        	return t1;
        } else {
        	return t2;
        }
	}

    vec3 normal(vec3 &point) {
    	vec3 normal = (point - center).normalize();
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
        float u = float(sample.x) / float(width);
        float v = float(sample.y) / float(height);
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

	int width, height, maxDepth;
	float eps;

    Scene(){}

	Scene(Camera &_cam, Film &_film, float _eps, int _width, int _height)  {
		cam = _cam;
		film = _film;
		width = _width;
		height = _height;
		maxDepth = 1;
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
        //cout << "rendering bitch" << endl;
    	Sampler sampler = Sampler(width, height);
    	Sample sample = Sample();
    	while (sampler.generateSample(&sample)) {
    		Ray testRay = cam.generateRay(sample);
    		vec3 color = traceRay(testRay, maxDepth);
    		film.writeToFilm(sample, color);
    	}

    }

    vec3 traceRay(Ray &ray, int depth) {
        //cout << "tracing yeah" << endl;
    	float min_t = MAXFLOAT;
    	int closeIndex = 0;
		for (int i = 0; i < objects.size(); i++) {
			float t = objects[i]->intersection(ray, eps);
            if (t != MAXFLOAT && t < min_t) {
                min_t = t;
                closeIndex = i;
		    }
		}
		if (min_t != MAXFLOAT) {
			Polygon *hitObject = objects[closeIndex];
		    vec3 intersection = ray.pointAt(min_t);
            vec3 normal = hitObject->normal(intersection);
            normal.normalize();
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
            if (objects[i]->intersects(ray)) {
                if (objects[i] == source) {
                    continue;
                }
                else if (objects[i] != source) {
                	return true;
                }
            }
        }
        return false;
    }

    vec3 phongShade(Ray &ray, Polygon *polygon, vec3 &normal, vec3 &point) {
        vec3 total_intensity = vec3(0, 0, 0);
        for (int i = 0; i < lights.size(); i++) {
        	Light *light = lights[i];
        	vec3 intensity = light->intensityAt(point);
        	Ray shadowRay = light->lightRay(point);
            if (traceShadowRay(shadowRay, polygon)) {
                continue;
            }
        	vec3 l = light->lightVec(point);
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
    // float eps = 0.05f;
    // string curline;
    // string argument;
    // Camera pcam;
    // Scene scene = Scene();
    // vec3 eye, ll, lr, ul, ur;
    // vec3 ka = vec3(0,0,0);
    // vec3 kd = vec3(0,0,0);
    // vec3 ks = vec3(0,0,0);
    // vec3 kr = vec3(0,0,0);
    // vector<Sphere> spheres;
    // vector<Triangle> triangles;
    // float sp = 0.0f;
    // int index = 0;
    // Material curMat = Material(ka,kd,ks, kr,sp);
    // int width = 1000;
    // int height = 1000;
    // Film film = Film(width, height);
    // ifstream myfile (argv[1]);
    // if (myfile.is_open()){
    //     while (getline (myfile, curline)) {
    //         cout << curline << '\n';
    //         istringstream iss;
    //         iss.str(curline);
    //         iss >> argument;

    //         if (curline.empty()){
    //             continue;
    //         }
    //         if (strcmp(argument.c_str(), "cam") == 0){

    //             float args[15];

    //             string temp;
    //             index = 0;
    //             while (iss >> temp){
    //                 cout << temp << endl;
    //                 args[index] = atof(temp.c_str());
    //                 index++;
    //             }

    //             if (index != 15){
    //                 cerr << "incorrect number of arguments for camera"<<endl;
    //             }

    //             eye = vec3(args[0], args[1], args[2]);
    //             ll = vec3(args[3], args[4], args[5]);
    //             lr = vec3(args[6], args[7], args[8]);
    //             ul = vec3(args[9], args[10], args[11]);
    //             ur = vec3(args[12], args[13], args[14]);

    //             pcam = Camera(eye, ul, ur, ll, lr, width, height);
    //             scene = Scene(pcam, film, eps, width, height);
    //         }
    //         else if (strcmp(argument.c_str(), "mat") == 0){
    //             float args[13];
    //             string temp;
    //             int index = 0;
    //             while (iss >> temp){
    //                 cout << index << endl;
    //                 args[index] = atof(temp.c_str());
    //                 index++;
    //             }

    //             if (index != 13){
    //                 cout << index << endl;
    //                 cerr << "incorrect number of arguments for material" << endl;
    //             }
    //             ka = vec3(args[0],args[1],args[2]);
    //             kd = vec3(args[3],args[4],args[5]);
    //             ks = vec3(args[6],args[7],args[8]);
    //             sp = args[9];
    //             kr = vec3(args[10],args[11],args[12]);
    //             curMat = Material(ka,kd,ks, kr,sp);
    //         }
    //         else if (strcmp(argument.c_str(), "sph") == 0){
    //             float args[4];
    //             string temp;
    //             int index = 0;
    //             while (iss >> temp){
    //                 cout << index << endl;
    //                 args[index] = atof(temp.c_str());
    //                 index++;
    //             }

    //             if (index != 4){
    //                 cout << index << endl;
    //                 cerr << "incorrect number of arguments for sphere" << endl;
    //             }
    //             vec3 sphereCenter = vec3(args[0],args[1],args[2]);
    //             //should be okay? won't be changing same sphere next time?
    //             Sphere sphere = Sphere(sphereCenter, args[3], curMat);
    //             spheres.push_back(sphere);
    //             //scene.addObject(&sphere);
    //         }

    //         else if (strcmp(argument.c_str(), "tri") == 0) {
    //             float args[9];
    //             string temp;
    //             int index = 0;
    //             while (iss >> temp){
    //                 cout << index << endl;
    //                 args[index] = atof(temp.c_str());
    //                 index++;
    //             }

    //             if (index != 9){
    //                 cout << index << endl;
    //                 cerr << "incorrect number of arguments for triangle" << endl;
    //             }
    //             vec3 v1 = vec3(args[0],args[1],args[2]);
    //             vec3 v2 = vec3(args[3],args[4],args[5]);
    //             vec3 v3 = vec3(args[6],args[7],args[8]);

    //             //should be okay? won't be changing same sphere next time?
    //             Triangle triangle = Triangle(v1, v2, v3, curMat);
    //             triangles.push_back(triangle);
    //         }

    //         else if (strcmp(argument.c_str(), "ltp") == 0){
    //             float args[6];
    //             int falloff = 0;
    //             string temp;
    //             int index = 0;
    //             while (iss >> temp){
    //                 cout << index << endl;
    //                 if (index == 6){
    //                     falloff = atoi(temp.c_str());
    //                 }
    //                 else{
    //                     args[index] = atof(temp.c_str());
    //                 }
    //                 index++;
    //             }
    //             if (index > 7 || index < 6){
    //                 cout << index << endl;
    //                 cerr << "incorrect number of arguments for point light" << endl;
    //             }
    //             vec3 color = vec3(args[3],args[4],args[5]);
    //             vec3 pos = vec3(args[0],args[1],args[2]);
    //             Light pl = Light(color, pos, falloff, true);
    //             scene.addLight(&pl);


    //         }
    //         else if (strcmp(argument.c_str(), "ltd") == 0){
    //             float args[6];
    //             string temp;
    //             int index = 0;
    //             while (iss >> temp){
    //                 cout << index << endl;
    //                args[index] = atof(temp.c_str());
    //                 index++;
    //             }
    //             if (index != 6){
    //                 cout << index << endl;
    //                 cerr << "incorrect number of arguments for directional light" << endl;
    //             }
    //             vec3 color = vec3(args[3],args[4],args[5]);
    //             vec3 pos = vec3(args[0],args[1],args[2]);
    //             Light dl = Light(color, pos, 0, false);
    //             scene.addLight(&dl);
    //         }

    //     }

    //     for (int j = 0; j < spheres.size(); j++){
    //         scene.addObject(&spheres[j]);
    //     }
    //     for (int j = 0; j < triangles.size(); j++){
    //         scene.addObject(&triangles[j]);
    //     }
    //     scene.render();
    //     scene.writeFile(width, height);
    //     myfile.close();

    // }
    // else {
    //     cerr << "No file! Or bad input path!" << endl;
    //     exit(1);
    // }
    // printf("Ray Tracer!\n");
    int width = 500;
    int height = 500;
    float eps = 0.05f;
    vec3 eye = vec3(0, 0, 0);
    vec3 ul = vec3(-1, 1, -3);
    vec3 ur = vec3(1, 1, -3);
    vec3 ll = vec3(-1, -1, -3);
    vec3 lr = vec3(1, -1, -3);
    Film film = Film(width, height);
    Camera cam = Camera(eye, ul, ur, ll, lr, width, height);
    Scene scene = Scene(cam, film, eps, width, height);

    vec3 ka = vec3(0.1, 0.1, 0.1);
    vec3 kd = vec3(1, 0, 1);
    vec3 ks = vec3(1, 1, 1);
    vec3 kr = vec3(0, 0, 0);
    float sp = 50.0f;

    vec3 ka1 = vec3(0.1, 0.1, 0.1);
    vec3 kd1 = vec3(1, 1, 0);
    vec3 ks1 = vec3(1, 1, 1);
    vec3 kr1 = vec3(0, 0, 0);
    float sp1 = 50.0f;

    vec3 ka2 = vec3(0.1, 0.1, 0.1);
    vec3 kd2 = vec3(0, 1, 1);
    vec3 ks2 = vec3(1, 1, 1);
    vec3 kr2 = vec3(0, 0, 0);
    float sp2 = 50.0f;

    vec3 ka3 = vec3(0.1, 0.1, 0.1);
    vec3 kd3 = vec3(0.1, 0.1, 0.1);
    vec3 ks3 = vec3(1, 1, 1);
    vec3 kr3 = vec3(1, 1, 1);
    float sp3 = 50.0f;

    Material sphereMat = Material(ka, kd, ks, kr, sp);
    Material sphereMat1 = Material(ka1, kd1, ks1, kr1, sp1);
    Material sphereMat2 = Material(ka2, kd2, ks2, kr2, sp2);
    Material triMat = Material(ka3, kd3, ks3, kr3, sp3);

    vec3 sphereCenter1 = vec3(0, 0, -20);
    Sphere sphere1 = Sphere(sphereCenter1, 3.0f, sphereMat);
    scene.addObject(&sphere1);
    
    vec3 sphereCenter2 = vec3(-2, 2, -15);
    Sphere sphere2 = Sphere(sphereCenter2, 1.0f, sphereMat1);
    scene.addObject(&sphere2);

    vec3 sphereCenter3 = vec3(-2, -2, -15);
    Sphere sphere3 = Sphere(sphereCenter3, 1.0f, sphereMat2);
    scene.addObject(&sphere3);

    vec3 v1 = vec3(5, 5, -17);
    vec3 v2 = vec3(1, 4, -20);
    vec3 v3 = vec3(6, -1, -20);
    Triangle triangle = Triangle(v1, v2, v3, triMat);
    scene.addObject(&triangle);

    vec3 lightColor1 = vec3(255, 255, 255);
    vec3 lightPos1 = vec3(0.57735027f, -0.57735027f, -0.57735027f);
    Light dl1 = Light(lightColor1, lightPos1, 0, false);
    scene.addLight(&dl1);

    vec3 lightColor2 = vec3(0, 0, 255);
    vec3 lightPos2 = vec3(0.57735027f, 0.57735027f, -0.57735027f);
    Light pl2 = Light(lightColor2, lightPos2, 0, false);
    scene.addLight(&pl2);

    scene.render();
    scene.writeFile(width, height);
        
	return 0;
}


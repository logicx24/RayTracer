/*********************************************
 ____             _                           
|  _ \ __ _ _   _| |_ _ __ __ _  ___ ___ _ __ 
| |_) / _` | | | | __| '__/ _` |/ __/ _ \ '__|
|  _ < (_| | |_| | |_| | | (_| | (_|  __/ |   
|_| \_\__,_|\__, |\__|_|  \__,_|\___\___|_|   
            |___/    

Aakash Japi
Jonathan Tau
Gabriel Intrator
*********************************************/


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
	return sqrt(pow(p1[VX] - p2[VX], 2) + pow(p1[VY] - p2[VY], 2) + pow(p1[VZ] - p2[VZ], 2));
}

vec3 multiplyComponents(vec3 v1, vec3 v2) {
	return vec3(v1[VX] * v2[VX], v1[VY] * v2[VY], v1[VZ] * v2[VZ]);
}

float magnitude(vec3 v1) {
	return sqrt(v1[VX] * v1[VX] + v1[VY] * v1[VY] + v1[VZ] * v1[VZ]);
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
	vec3 pos, vec;
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
    vec3 color, pos;
    int falloff;
    bool pointLight;

    Light(){}

    Light(vec3 &_color, vec3 &_pos, int _falloff, bool _pointLight) {
    	color = _color;
    	pos = _pos;
    	falloff = _falloff;
    	pointLight = _pointLight;
    }

    Light(vec3 &ambientColor) {
    	color = ambientColor;
    	pos = vec3(0, 0, 0);
    	falloff = 0;
    	pointLight = false;
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
    int width, height, x, y;

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

    int percentage() {
        return int((float(y * width + x) / float(width * height)) * 100);
    }

};

vec3 floatToByte(vec3 &fcolor) {
	return vec3((int)(fcolor[RED] * 255.0f), 
	            (int)(fcolor[GREEN] * 255.0f), 
	            (int)(fcolor[BLUE] * 255.0f));
}

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

    void writeToFilm(Sample &sample, vec3 &fColor) {
    	vec3 color = floatToByte(fColor);
    	if (color[RED] > 255) {color[RED] = 255;}
    	if (color[GREEN] > 255) {color[GREEN] = 255;}
    	if (color[BLUE] > 255) {color[BLUE] = 255;}
        film.at(sample.y * width + sample.x)[RED] = color[RED];
        film.at(sample.y * width + sample.x)[GREEN] = color[GREEN];
        film.at(sample.y * width + sample.x)[BLUE] = color[BLUE];
    }

    void writeFile(string fileName, int pixWidth, int pixHeight) {
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
        unsigned error = lodepng::encode(fileName, image,  width, height);
        if (error) {
        	cout << "lodepng error." << endl;
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
    vec3 vertex1, vertex2, vertex3, edge1, edge2;

    Triangle(vec3 &_v1, vec3 &_v2, vec3 &_v3, Material &_mat) {
        vertex1 = _v1;
        vertex2 = _v2;
        vertex3 = _v3;
        mat = _mat;
        edge1 = vertex3 - vertex1;
        edge2 = vertex2 - vertex1;
    }

    bool intersects(Ray &ray) {
        return intersection(ray, 0.1f) != MAXFLOAT;

    }

    float intersection(Ray &ray, float eps) {
        vec3 normal = (edge2 ^ edge1).normalize(); 
        vec3 w0 = ray.pos - vertex1;
        float a = -(normal * w0);
        float b = normal * ray.vec;
        float r = a/b;

        if (fabs(b) < eps) {
            if (a == 0.0f) {
                return r;
            } else {
            	return MAXFLOAT;
            }
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
        vec3 normal = (edge2 ^ edge1);
        return normal;
    }
};

class Sphere : public Polygon {
public:
	vec3 center;
	float radius;

    Sphere(){}

	Sphere(vec3 &_center, float _radius, Material &_mat) {
		center = _center;
		radius = _radius;
		mat = _mat;
	}

	bool intersects(Ray &ray) {
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

class Ellipse : public Sphere {
public:

    mat4 transform;
    mat4 inverse;
    mat4 inverseTranspose;

    Ellipse(vec3 &_center, float _radius, Material &_mat, mat4 &_transform) {
        center = _center;
        radius = _radius;
        mat = _mat;
        transform = _transform;
        inverse = transform.inverse();
        inverseTranspose = inverse.transpose();
    }

    bool intersects(Ray &ray) {
        return intersection(ray, 0.1f) != MAXFLOAT;
    }

    float intersection(Ray &ray, float eps) {
        vec3 transing = inverse * ray.pos;
        vec4 tmp = vec4(ray.vec[VX], ray.vec[VY], ray.vec[VZ], 0.0f);
        tmp = inverse * tmp;
        vec3 transing2 = vec3(tmp[VX], tmp[VY], tmp[VZ]);
        //cout <<  "t1" << endl;
        //printVec3(transing);
        //cout << "t2" << endl;
        //printVec3(transing2);
        Ray ray2 = Ray(transing, transing2);
        return Sphere::intersection(ray2, eps);
    }

    vec3 normal(vec3 &point) {
        vec3 normal = (point - center);
        vec4 temp = vec4(normal[VX], normal[VY], normal[VZ], 0.0f);
        //cout << "temp: " << temp << endl;
        //cout << "iT: " << endl << inverseTranspose << endl;
        vec4 normal1 = inverseTranspose * temp;
        //cout << "normal: " << normal1 << endl;
        normal = vec3(normal1[VX], normal1[VY], normal1[VZ]);
        //cout << "3D normal: " << normal << endl;
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
	    return Ray(eye, testRayVec);
	}
};


class Scene {
public:
	Camera cam;
	Film film;
	vector<Polygon*> objects;
	vector<Light*> lights;
	Light ambientLight;

	int width, height, maxDepth;
	float eps;

    Scene(){}

	Scene(Camera &_cam, Film &_film, float _eps, int _width, int _height)  {
		cam = _cam;
		film = _film;
		vec3 defaultAmbientColor = vec3(0, 0, 0);
		Light ambientLight = Light(defaultAmbientColor);
		width = _width;
		height = _height;
		maxDepth = 3;
		eps = _eps;
	}

	void addObject(Polygon *p) {
 		objects.push_back(p);
	}

    void addLight(Light *l) {
    	lights.push_back(l);
    }

    void addAmbientLight(Light &l) {
    	ambientLight = l;
    }

    void writeFile(string fileName, int pixWidth, int pixHeight) {
    	film.writeFile(fileName, pixWidth, pixHeight);
    }

    void render() {
        Sampler sampler = Sampler(width, height);
        Sample sample = Sample();
        int count = 0;
        while (sampler.generateSample(&sample)) {
            Ray testRay = cam.generateRay(sample);
            vec3 color = traceRay(testRay, maxDepth);
            film.writeToFilm(sample, color);
            if (count > 100000) {
                cout << sampler.percentage() << "\% completed." << endl;
                count = 0;
            }
            count++;
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
            vec3 reflectionVec = ray.vec - (2*(ray.vec * normal))*normal; 
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
        total_intensity += multiplyComponents(ambientLight.color, polygon->mat.ka);
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
            vec3 coeffs = polygon->mat.kd * MAX(l * normal, 0.0f)
                        + polygon->mat.ks * pow(MAX(r * v, 0.0f), polygon->mat.sp);
            total_intensity += multiplyComponents(intensity, coeffs);
        }
        return total_intensity;
    }

};

void parseArgs(float args[], istringstream *iss, int numArgs, string errMessage) {
    string temp;
    int index = 0;
    while (*iss >> temp){
        args[index] = atof(temp.c_str());
        index++;
    }

    if (index != numArgs){
        cerr << errMessage << endl;
    }
}

void parseArgs(int args[], istringstream *iss, int numArgs, string errMessage) {
    string temp;
    int index = 0;
    while (*iss >> temp){
        args[index] = atoi(temp.c_str());
        index++;
    }

    if (index != numArgs){
        cerr << errMessage << endl;
    }

}

bool argumentMatches(string argument, const char *compare) {
    return (strcmp(argument.c_str(), compare) == 0);
}

void parseObj(string fileName, vector<Polygon*> *sceneObjects, Material &curMat) {
	cout << "Parsing obj file: " << fileName << "..." <<  endl;

	vector<vec3> verts;
	verts.push_back(vec3(0, 0, 0)); // Kill index 0, since vertices are numbered 1 ... n
	ifstream objFile(fileName);
	string currentLine;
	string argument;
    if (!objFile.is_open()) {
    	cerr << "Error opening obj file" << endl;
    	return;
    }
    while (getline(objFile, currentLine)) {
    	istringstream iss;
    	iss.str(currentLine);
    	iss >> argument;
        if (currentLine.empty()) {
        	continue;
        }
        if (argumentMatches(argument, "#")) {
        	continue;
        } else if (argumentMatches(argument, "v")) {
        	//cout << "found vertice" << endl;
        	float args[3];
            parseArgs(args, &iss, 3, "Incorrect number of arguments for vertex in obj file.");
            verts.push_back(vec3(args[0], args[1], args[2]));
        } else if (argumentMatches(argument, "f")) {
        	int args[3];
        	parseArgs(args, &iss, 3, "More than 3 vertices specified for face, only using first 3.");
            Triangle *t = new Triangle(verts[args[0]], verts[args[1]], verts[args[2]], curMat);
            sceneObjects->push_back(t);
        }
    }
    objFile.close();
    cout << "Finished parsing obj file." << endl;
}

Scene* parseSceneFile(string fileName, float eps, int width, int height) {
    float sp;
    vec3 eye, ll, lr, ul, ur, ka, kd, ks, kr;

    mat4 currTrans = identity3D();
    
    Film film = Film(width, height);
    Camera pcam;
    Scene *scene;
    
    Material curMat;
    vector<Polygon*> sceneObjects;
    vector<Light*> sceneLights;
    vector<vec3> vertices;

    int index = 0;
    string argument;
    string curline;
    ifstream myfile (fileName);
    if (myfile.is_open()){
        while (getline (myfile, curline)) {
            //cout << curline << '\n';
            istringstream iss;
            iss.str(curline);
            iss >> argument;

            if (curline.empty()){
                continue;
            }
            if (argumentMatches(argument, "cam")) {
                float args[15];
                parseArgs(args, &iss, 15, "Incorrect number of arguments for camera");
                eye = vec3(args[0], args[1], args[2]);
                ll = vec3(args[3], args[4], args[5]);
                lr = vec3(args[6], args[7], args[8]);
                ul = vec3(args[9], args[10], args[11]);
                ur = vec3(args[12], args[13], args[14]);

                pcam = Camera(eye, ul, ur, ll, lr, width, height);
                scene = new Scene(pcam, film, eps, width, height);
            
            } else if (argumentMatches(argument, "mat")) {
                float args[13];
                parseArgs(args, &iss, 13, "Incorrect number of arguments for material.");
                ka = vec3(args[0],args[1],args[2]);
                kd = vec3(args[3],args[4],args[5]);
                ks = vec3(args[6],args[7],args[8]);
                sp = args[9];
                kr = vec3(args[10],args[11],args[12]);
                curMat = Material(ka,kd,ks, kr,sp);
            
            } else if (argumentMatches(argument, "sph")) {
                float args[4];
                parseArgs(args, &iss, 4, "Incorrect number of arguments for sphere.");
                vec3 sphereCenter = vec3(args[0],args[1],args[2]);
                if (currTrans == identity3D()) {
                    sceneObjects.push_back(new Sphere(sphereCenter, args[3], curMat));
                } else {
                    sceneObjects.push_back(new Ellipse(sphereCenter, args[3], curMat, currTrans));
                }
            
            } else if (argumentMatches(argument, "tri")) {
                float args[9];
                parseArgs(args, &iss, 9, "Incorrect number of arguments for triangle.");
                vec3 v1 = vec3(args[0],args[1],args[2]);
                vec3 v2 = vec3(args[3],args[4],args[5]);
                vec3 v3 = vec3(args[6],args[7],args[8]);

                sceneObjects.push_back(new Triangle(v1, v2, v3, curMat));
            
            } else if (argumentMatches(argument, "ltp")) {
                float args[6];
                int falloff = 0;
                string temp;
                int index = 0;
                while (iss >> temp){
                    if (index == 6){
                        falloff = atoi(temp.c_str());
                    } else {
                        args[index] = atof(temp.c_str());
                    }
                    index++;
                }
                if (index > 7 || index < 6){
                    cerr << "Incorrect number of arguments for point light." << endl;
                }
                vec3 color = vec3(args[3],args[4],args[5]);
                vec3 pos = vec3(args[0],args[1],args[2]);
                sceneLights.push_back(new Light(color, pos, falloff, true));

            } else if (argumentMatches(argument, "ltd")) {
                float args[6];
                parseArgs(args, &iss, 6, "Incorrect number of arguments for direct light.");
                vec3 color = vec3(args[3],args[4],args[5]);
                vec3 pos = vec3(args[0],args[1],args[2]);
                sceneLights.push_back(new Light(color, pos, 0, false));
            } else if (argumentMatches(argument, "lta")) {
                float args[3];
                parseArgs(args, &iss, 3, "Incorrect number of arguments for ambient light.");
                vec3 ambientColor = vec3(args[0], args[1], args[2]);
                Light ambientLight = Light(ambientColor);
                scene->addAmbientLight(ambientLight);
            } else if (argumentMatches(argument, "obj")) {
                string temp;
                iss >> temp;
                parseObj(temp.c_str(), &sceneObjects, curMat);
            } else if (argumentMatches(argument, "xfs")) {
                float args[3];
                parseArgs(args, &iss, 3, "Incorrect number of arguments for scaling transformation.");
                vec3 transformation = vec3(args[0], args[1], args[2]);
                currTrans = currTrans * scaling3D(transformation); 
            } else if (argumentMatches(argument, "xft")) {
                float args[3];
                parseArgs(args, &iss, 3, "Incorrect number of arguments for translation transformation.");
                vec3 transformation = vec3(args[0], args[1], args[2]);
                currTrans = currTrans * translation3D(transformation); 
            } else if (argumentMatches(argument, "xfr")) {
                float args[3];
                parseArgs(args, &iss, 3, "Incorrect number of arguments for rotation transformation.");
                vec3 transformation = vec3(args[0], args[1], args[2]);
                currTrans = currTrans * rotation3D(transformation, magnitude(transformation)); 
            } else if (argumentMatches(argument, "xfz")) {
                currTrans = identity3D();
            }
        } 

        for (int i = 0; i < sceneObjects.size(); i++){
            scene->addObject(sceneObjects[i]);
        }
        for (int i = 0; i < sceneLights.size(); i++) {
            scene->addLight(sceneLights[i]);
        }
        
        myfile.close();

    } else {
        cerr << "No scene file found." << endl;
        exit(0);
    }
    return scene;
}

int main(int argc, char* argv[]) {
   
    cout << "Ray Tracer!" << endl;
     
    string outFile = "output.png";

    float eps = 0.05f;
    int width = 1000;
    int height = 1000;
   
    Scene *scene = parseSceneFile(argv[1], eps, width, height);

    cout << "Rendering Scene..." << endl;

    scene->render();
    
    cout << "Scene render completed." << endl;
    cout << "Writing to file..." << endl;
    
    if (argc > 2) { outFile = argv[2]; }
    scene->writeFile(outFile, width, height);
    
    cout << "Image saved." << endl;
        
	return 0;
}


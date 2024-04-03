#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static std :: default_random_engine engine[8] ;
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) { return x * x;}

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	const Vector& operator+=(const Vector& b){
		data[0] += b[0];
		data[1] += b[1];
		data[2] += b[2];
		return (*this);
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0]* b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector random_cos(const Vector& N){
	int tid = omp_get_thread_num();
	double r1 = uniform(engine[tid]);
	double r2 = uniform(engine[tid]);
	double x = cos(2. * M_PI * r1) * sqrt(1. - r2);
	double y = sin(2. * M_PI * r1) * sqrt(1. - r2);
	double z = sqrt(r2);

	Vector T1;
	if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])){
		T1 = Vector(0,-N[2], N[1]);
	}else {
		if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])){
			T1 = Vector(-N[2],0, N[0]);
		} else {
			T1 = Vector(-N[1],N[0], 0);
		}}

	T1.normalize();
	Vector T2 = cross(N, T1);

	return x * T1 + y * T2 + z*N;
}

class Ray {
public:
	Ray(const Vector& O, const Vector& u): O(O), u(u) {};

	Vector O, u;
};

class Sphere {
public:
	Sphere(const Vector& C, double R, const Vector& albedo, bool isMirror = false, bool isTransparent = false) : C(C), albedo(albedo), R(R), isMirror(isMirror), isTransparent(isTransparent) {};

	bool intersect(const Ray& r, Vector& P, Vector& N, double &t) const {
		double delta = sqr(dot(r.u, r.O-C)) - ((r.O-C).norm2() - R*R);
		if (delta < 0) return false;

		// Compute the two ray sphere intersection
		double t1 = dot(r.u,C-r.O) - sqrt(delta);
		double t2 = dot(r.u,C-r.O) + sqrt(delta);

		if (t2 < 0) return false;
		
		if (t1 > 0) {
			t = t1;
		} else {
			t = t2;
		}
		P = r.O + t * r.u;
		N = P - C;
		N.normalize();
		return true;
	}

	Vector C, albedo;
	double R;
	bool isMirror, isTransparent;
};

class Scene{
public:
	bool intersect(const Ray& r, Vector& P, Vector& N, int& objectID, double& bestt ) const {
		bestt = 1E9;
		double t;
		Vector Ptmp, Ntmp;
		bool has_inter = false;

		for (int i=0; i< objects.size(); i++){
			if (objects[i].intersect(r, Ptmp, Ntmp, t)) {
				if (t < bestt){
					bestt = t;
					P = Ptmp;
					N = Ntmp;
					objectID = i;
					has_inter = true;
				}
			}
		}
		return has_inter;
	}

	void addSphere(const Sphere &s) {
		objects.push_back(s);
	}

	Vector getColor(const Ray& r, int bounce){
		Vector color(0,0,0);
		if (bounce < 0) return color;

		double t;
		int objectID;
		Vector P,N;
		bool inter = intersect(r, P, N, objectID, t);
			
		if (inter){

			if (objects[objectID].isMirror) {
				Ray reflected(P+0.001 * N, r.u - 2*dot(r.u, N) * N); // direction of reflection
				return getColor(reflected, bounce-1); // recur on surfaces to continue mirroring
			}
			if (objects[objectID].isTransparent) { // computations for transparent surfaces following lecture 1
				double n1 = 1;
				double n2 = 1.5;

				Vector correctN = N;
				if (dot(N, r.u) > 0){
					correctN = -N; 
					std::swap(n1, n2);
				}

				Vector Tt = n1/ n2 * (r.u - dot(r.u, correctN) * correctN);
				double d = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, correctN)));
				if (d < 0) {
					Ray reflected(P+0.001 * correctN, r.u - 2 * dot(r.u, correctN) * correctN);
					return getColor(reflected, bounce-1);
				}

				Vector Tn = -sqrt(d) * correctN;
				Vector T = Tn + Tt;

				Ray refracted(P - 0.001 * correctN, T);
				return getColor(refracted, bounce-1);
			} // here add fresnel transparent

			// direct light
			Vector wlight = L - P;
			double dlight2 = wlight.norm2();
			wlight.normalize();
			double tshadow;
			Vector Pshadow, Nshadow;
			int objectshadow;
			Ray rShadow(P + 0.001 * N, wlight);

			double l = I/(4 * M_PI * dlight2) * std::max(0.0,dot(N, wlight)); // color with direct light
			color = l * objects[objectID].albedo / M_PI; // add the material in color
			if (intersect(rShadow, Pshadow, Nshadow, objectshadow, tshadow)){
				if (sqr(tshadow) < dlight2){
					color = Vector(0,0,0);
				}
			}

			// add indirect light
			Vector wi = random_cos(N);
			Ray indirectRay(P + 0.001 * N, wi); 
			return objects[objectID].albedo * getColor(indirectRay, bounce - 1) + color;
			/* double diffuse_probability = rho_d/(rho_d+rho_s) ; // we should use some color←-average of rho d and rho s
			if (uniform(engine) < diffuse_probability) { 
				// we sample the diffuse lobe
			Ray randomRay = . . . ; // randomly sample ray using random cos
			Lo += albedo/diffuse_probability * getColor(randomRay , ray_depth=1); }
			else{
			Ray randomRay = . . . ; // randomly sample ray using random pow and mirroring ←- of ray. direction
			if (dot(randomRay.direction, N) < 0) return Vector(0., 0., 0.); // make sure←- we sampled the upper hemisphere
			Vector brdf_indirect = rho_s * (alpha+8)/(8*M_PI)* PhongSpecularLobe (...) ; ←- // just the specular part of the Blinn=Phong model
			double pdf_pow = . . . ; // the pdf associated with our function random pow ←- with the reflection
			Lo += brdf_indirect* std::max(dot(N, randomRay.direction), 0.)/((1=-diffuse_probability)*pdf_pow) * getColor(randomRay , ray_depth=1) ;*/
		}
		return color;
	}

	std::vector<Sphere> objects;
	Vector L;
	double I;
};

int main() {
	int W = 512;
	int H = 512;
	int nrays = 64;

	double fov = 60 * M_PI / 180;
	double z = -W/(2*tan(fov/2));

	for (int i = 0; i<8; i++){
		engine[i].seed(i);
	}

	Scene s;
	// Spheres
	s.addSphere(Sphere(Vector(0,0,0), 10, Vector(1, 0.5, 0.3), false, false));
	s.addSphere(Sphere(Vector(-10,0,10), 10, Vector(1, 0.5, 0.3), true, false));
	s.addSphere(Sphere(Vector(10,0,-10), 10, Vector(1, 0.5, 0.3), false, true));
	// Ceiling and floor
	s.addSphere(Sphere(Vector(1000,0,0), 900, Vector(0.7, 0.4, 0.2))); 
	s.addSphere(Sphere(Vector(-1000,0,0), 900, Vector(0.5, 0.9, 0.1)));
	// Walls
	s.addSphere(Sphere(Vector(0,1000,0), 940, Vector(0.2, 0.4, 0.7))); // Change to red
	s.addSphere(Sphere(Vector(0,-1000,0), 990, Vector(0.2, 0.7, 0.4))); // Change to blue 
	s.addSphere(Sphere(Vector(0,0,1000), 940, Vector(0.5, 0.4, 0.5))); // change to pink
	s.addSphere(Sphere(Vector(0,0,-1000), 940, Vector(0.9, 0.3, 0.1))); // change to green

	Vector C(0,0,55);
	s.L = Vector(-10, 20, 40); // Light ray
	s.I = 2E10;

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		int tid = omp_get_thread_num();
		for (int j = 0; j < W; j++) {

			Vector color;
			for (int k = 0; k < nrays; k++){

				// box Muller
				double u1 = uniform(engine[tid]);
				double u2 = uniform(engine[tid]);
				double r1 = cos(2*M_PI*u1)*sqrt(-2*log(u2))*0.4;
				double r2 = sin(2*M_PI*u1)*sqrt(-2*log(u2))*0.4;

				// vector of the view u
				Vector u(j-W/2+0.5+r1, H/2-i-0.5+r2, z);
				u.normalize();
				Ray r(C, u);
				color += s.getColor(r, 5)/nrays;		
			}

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::pow(color[0], 0.45)); // including gamma correction
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::pow(color[1], 0.45)); 
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Ray(const Vector& O, const Vector& u): O(O), u(u) {};

	Vector O, u;
};

class Sphere {
public:
	Sphere(const Vector& C, double R, const Vector& albedo, bool isMirror = false, bool isTransparent = false) : C(C), R(R), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {};

	bool intersect(const Ray& r, Vector& P, Vector& N, double &t) const {
		double delta = sqr(dot(r.u, r.O-C)) - ((r.O-C).norm2() - R*R);
		if (delta < 0) return false;

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
				Ray reflected(P+0.001 * N, r.u - 2*dot(r.u, N) * N);
				return getColor(reflected, bounce-1);
			}
			if (objects[objectID].isTransparent) {
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
			}

			Vector wlight = L - P;
			double dlight2 = wlight.norm2();
			wlight.normalize();
			double tshadow;
			Vector Pshadow, Nshadow;
			int objectshadow;
			Ray rShadow(P + 0.001 * N, wlight);

			double l = I/(4 * M_PI * dlight2) * std::max(0.0,dot(N, wlight));
			color = l * objects[objectID].albedo / M_PI;
			if (intersect(rShadow, Pshadow, Nshadow, objectshadow, tshadow)){
				if (sqr(tshadow) < dlight2){
					color = Vector(0,0,0);
				}
			}
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

	double fov = 60 * M_PI /100;
	double z = -W/(2*tan(fov/2));

	Scene s;
	s.addSphere(Sphere(Vector(0,0,0), 10, Vector(1, 0.5, 0.3), true, false));
	s.addSphere(Sphere(Vector(-20,0,0), 10, Vector(1, 0.5, 0.3), false, false));
	s.addSphere(Sphere(Vector(20,0,0), 10, Vector(1, 0.5, 0.3), false, true));

	s.addSphere(Sphere(Vector(0,1000,0), 960, Vector(0.2, 0.4, 0.7)));
	s.addSphere(Sphere(Vector(0,-1000,0), 990, Vector(0.2, 0.7, 0.4)));
	s.addSphere(Sphere(Vector(1000,0,0), 900, Vector(0.7, 0.4, 0.2))); 
	s.addSphere(Sphere(Vector(-1000,0,0), 900, Vector(0.5, 0.9, 0.1)));
	s.addSphere(Sphere(Vector(0,0,1000), 900, Vector(0.5, 0.4, 0.5)));
	s.addSphere(Sphere(Vector(0,0,-1000), 900, Vector(0.9, 0.3, 0.1)));

	Vector C(0,0,55);
	s.L = Vector(-10, 20, 40);

	s.I = 2E10;

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector u(j-W/2+0.5, H/2-i-0.5, z);
			u.normalize();
			Ray r(C, u);

			Vector color = s.getColor(r, 5);

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::pow(color[0], 0.45));
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::pow(color[1], 0.45));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}
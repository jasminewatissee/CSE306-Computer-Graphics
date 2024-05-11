#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

static std::default_random_engine engine[8] ;
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
	Ray(const Vector& Origin, const Vector& direction): Origin(Origin), direction(direction) {};

	Vector Origin, direction;
};

class Geometry{
public:
	Geometry(const Vector& albedo, bool isMirror, bool isTransparent, bool isInverted) : albedo(albedo), isMirror(isMirror), isTransparent(isTransparent), isInverted(isInverted) {}

	virtual bool intersect(const Ray& ray, Vector& P, Vector& N, double &t) const = 0;

	Vector albedo;
	bool isMirror, isTransparent, isInverted;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {}
	
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class BoundingBox {
public:
	Vector m, M;

	bool intersect(const Ray& ray) const {
		Vector invr(1./ray.direction[0], 1. / ray.direction[1], 1./ray.direction[2]);
		double t0x = (m[0] - ray.Origin[0]) * invr[0];
		double t1x = (M[0] - ray.Origin[0]) * invr[0];
		double tminx = std::min(t0x, t1x);
		double tmaxx = std::max(t0x, t1x);

		double t0y = (m[1] - ray.Origin[1]) * invr[1];
		double t1y = (M[1] - ray.Origin[1]) * invr[1];
		double tminy = std::min(t0y, t1y);
		double tmaxy = std::max(t0y, t1y);

		double t0z = (m[2] - ray.Origin[2]) * invr[2];
		double t1z = (M[2] - ray.Origin[2]) * invr[2];
		double tminz = std::min(t0z, t1z);
		double tmaxz = std::max(t0z, t1z);

		double tmax = std::min(tmaxx, std::min(tmaxy, tmaxz));
		double tmin = std::max(tminx, std::max(tminy, tminz));
		if (tmax > tmin && tmax > 0) {
			return true;
		}
		return false;

	}
};

class BVH {
public:
	BVH *left, *right;
	int start, end;
	BoundingBox bbox;
};

class TriangleMesh : public Geometry {
public:
	TriangleMesh(const Vector& albedo, bool isMirror = false, bool isTransparent = false, bool isInverted = false) : ::Geometry(albedo, isMirror, isTransparent, isInverted) {};
	~TriangleMesh() {};

	void translate_and_scale(const Vector& translate, double scale) {
		for (int i=0; i < vertices.size(); i++) {
			vertices[i] = scale * vertices[i] + translate;
		}
	}

	BoundingBox compute_bbox(int triangle_start, int triangle_end) {
		BoundingBox bbox;

		bbox.m = Vector(1E9, 1E9, 1E9);
		bbox.M = Vector(-1E9, -1E9, -1E9);

		for (int i = triangle_start; i < triangle_end; i++){
			for (int j = 0; j < 3; j++){
				bbox.m[j] = std::min(bbox.m[j], vertices[indices[i].vtxi][j]);
				bbox.M[j] = std::max(bbox.M[j], vertices[indices[i].vtxi][j]);
				bbox.m[j] = std::min(bbox.m[j], vertices[indices[i].vtxj][j]); 
				bbox.M[j] = std::max(bbox.M[j], vertices[indices[i].vtxj][j]);
				bbox.m[j] = std::min(bbox.m[j], vertices[indices[i].vtxk][j]);
				bbox.M[j] = std::max(bbox.M[j], vertices[indices[i].vtxk][j]);
			}
		}
		return bbox;
	}

	bool intersect(const Ray& ray, Vector& P, Vector& N, double &t) const {
		bool hasInter = false;
		t = 1E9;

		if (!bvh.bbox.intersect(ray)) return false;

		// Without BVH
		/* for (int i=0; i<indices.size(); i++){
			const Vector &A = vertices[indices[i].vtxi];
			const Vector &B = vertices[indices[i].vtxj];
			const Vector &C = vertices[indices[i].vtxk];

			Vector e1 = B - A;
			Vector e2 = C - A;
			Vector localN = cross(e1, e2);
			Vector OA = A - ray.Origin;
			double invUdotN = 1. / (dot(ray.direction, localN));
			double local_t = dot(OA, localN) * invUdotN;
			if (local_t < 0) continue;

			Vector OAcrossU = cross(OA, ray.direction);
			double beta = dot(e2, OAcrossU) * invUdotN;
			if (beta < 0) continue;
			if (beta > 1) continue;

			double gamma = -dot(e1, OAcrossU) * invUdotN;
			if (gamma < 0) continue;

			double alpha = 1 - beta - gamma;
			if (alpha < 0) continue;

			hasInter = true;
			if (local_t < t) {
				t = local_t;
				P = A + beta * e1 + gamma * e2;
				N = localN;
			}
		}
		if (hasInter){
			N.normalize();
			if (isInverted){
				N = (-1.)*N;
			}
		}
		return hasInter; */
		
		//std::list<BVH*> l;
		const BVH* l[32];
		int list_size = 0;
		//l.push_back(&bvh);
		l[list_size] = &bvh;
		list_size++;

		while (list_size != 0){
			//const BVH* cur = l.front();
			const BVH* cur = l[list_size-1];
			//l.pop_front();
			list_size--;

			if (cur-> left){
				if (cur->left-> bbox.intersect(ray)){
					//l.push_back(cur->left);
					l[list_size] = cur -> left;
					list_size++;
				}
				if (cur->right-> bbox.intersect(ray)){
					//l.push_back(cur->right);
					l[list_size] = cur -> right;
					list_size++;
				}
			} else{
					for (int i=0; i<indices.size(); i++){
						const Vector &A = vertices[indices[i].vtxi];
						const Vector &B = vertices[indices[i].vtxj];
						const Vector &C = vertices[indices[i].vtxk];
						Vector e1 = B - A;
						Vector e2 = C - A;
						Vector localN = cross(e1, e2);
						Vector OA = A - ray.Origin;
						double invUdotN = 1. / (dot(ray.direction, localN));
						double local_t = dot(OA, localN) * invUdotN;
						if (local_t < 0) continue;

						Vector OAcrossU = cross(OA, ray.direction);
						double beta = dot(e2, OAcrossU) * invUdotN;
						if (beta < 0) continue;
						if (beta > 1) continue;

						double gamma = -dot(e1, OAcrossU) * invUdotN;
						if (gamma < 0) continue;
						if (gamma > 1) continue;

						double alpha = 1 - beta - gamma;
						if (alpha < 0) continue;

						hasInter = true;
						if (local_t < t) {
							t = local_t;
							P = A + beta * e1 + gamma * e2;
							N = localN;
						}
					}
					if (hasInter){
						N.normalize();
						if (isInverted){
							N = (-1.)*N;
						}
					}
					return hasInter;
			}
		}
		return false;
	}

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	void buildBVH(BVH* curNode, int triStart, int triEnd){
		curNode -> start = triStart;
		curNode -> end = triEnd;
		curNode -> left = NULL;
		curNode -> right = NULL;
		curNode -> bbox = compute_bbox(triStart, triEnd);

		Vector diag = curNode -> bbox.M - curNode -> bbox.m;
		int longestAxis;
		if (diag[0] >= diag[1] && diag[0] >= diag[2]){
			longestAxis = 0;
		} else if (diag[1] >= diag[0] && diag[1] >= diag[2]) {
			longestAxis = 1;
		} else {
			longestAxis = 2;
		}

		double middleSplit = curNode -> bbox.m[longestAxis] + diag[longestAxis] * 0.5;

		int pivot = triStart;
		for (int i= triStart; i < triEnd; i++){
			double middleTriangle = (vertices[indices[i].vtxi][longestAxis] + vertices[indices[i].vtxj][longestAxis] + vertices[indices[i].vtxk][longestAxis]) / 3;
			if (middleTriangle < middleSplit) {
				std::swap(indices[i], indices[pivot]);
				pivot++;
			}
		}

		if (pivot-triStart <= 1 || triEnd - pivot <= 1 || triEnd-triStart < 5){
			return;
		}

		curNode -> left = new BVH;
		curNode -> right = new BVH;

		buildBVH(curNode->left, triStart, pivot);
		buildBVH(curNode-> right, pivot, triEnd);
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	BVH bvh;
};

class Sphere : public Geometry {
public:
	Sphere(const Vector& Center, double Radius, const Vector& albedo, bool isMirror = false, bool isTransparent = false, bool isInverted = false) : ::Geometry(albedo, isMirror, isTransparent, isInverted), Center(Center), Radius(Radius) {};

	bool intersect(const Ray& ray, Vector& P, Vector& N, double &t) const {
		double delta = sqr(dot(ray.direction, ray.Origin - Center)) - ((ray.Origin - Center).norm2() - Radius * Radius);
		if (delta < 0) return false;

		// Compute the two ray sphere intersection
		double t1 = dot(ray.direction, Center - ray.Origin) - sqrt(delta);
		double t2 = dot(ray.direction, Center - ray.Origin) + sqrt(delta);

		if (t2 < 0) return false;
		
		if (t1 > 0) {
			t = t1;
		} else {
			t = t2;
		}
		P = ray.Origin + t * ray.direction;
		N = P - Center;
		N.normalize();

		if (isInverted){
			N = -N;
		}
		return true;
	}

	Vector Center;
	double Radius;
};

class Scene{
public:
	bool intersect(const Ray& ray, Vector& P, Vector& N, int& objectID, double& best_t ) const {
		best_t = 1E9;
		double t;
		Vector Ptmp, Ntmp;
		bool has_inter = false;

		for (int i=0; i< objects.size(); i++){
			if (objects[i]->intersect(ray, Ptmp, Ntmp, t)) {
				if (t < best_t){
					best_t = t;
					P = Ptmp;
					N = Ntmp;
					objectID = i;
					has_inter = true;
				}
			}
		}
		return has_inter;
	}

	void addSphere(Sphere* s) {
		objects.push_back(s);
	}

	void addMesh(TriangleMesh* m) {
		objects.push_back(m);
	}

	Vector getColor(const Ray& r, int bounce){
		Vector color(0,0,0);
		if (bounce < 0) return color;

		double t;
		int objectID;
		Vector P,N;
		bool inter = intersect(r, P, N, objectID, t);
			
		if (inter){

			if (objects[objectID]->isMirror) { // computations for mirror surfaces (reflection)
				Ray reflected(P + 0.001 * N, r.direction - 2 * dot(r.direction, N) * N); // ray of reflection
				return getColor(reflected, bounce-1); // recur on surfaces to continue mirroring
			}
			if (objects[objectID]->isTransparent) { // computations for transparent surfaces following lecture 1 (refraction)
				double n1 = 1;
				double n2 = 1.5;

				Vector correctN = N;
				if (dot(N, r.direction) > 0){
					correctN = -N; 
					std::swap(n1, n2);
				}

				Vector Tt = n1/ n2 * (r.direction - dot(r.direction, correctN) * correctN);
				double d = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.direction, correctN)));
				if (d < 0) {
					Ray reflected(P + 0.001 * correctN, r.direction - 2 * dot(r.direction, correctN) * correctN);
					return getColor(reflected, bounce-1);
				}

				Vector Tn = -sqrt(d) * correctN;
				Vector T = Tn + Tt;

				// Fresnel
				double k0 = sqr(n1 - n2) / sqr(n1 + n2);
				double R = k0 + (1 - k0)* pow(1 - std::abs(dot(correctN, r.direction)), 5); // reflection coefficient
				double u = (double)rand()/(double)RAND_MAX;
				if (u < R){
					Ray reflected(P + 0.001 * correctN, r.direction - 2 * dot(r.direction, correctN) * correctN);
					return getColor(reflected, bounce-1);
				}

				Ray refracted(P - 0.001 * correctN, T);
				return getColor(refracted, bounce-1);
			} 

			// direct light
			Vector wlight = L - P;
			double dlight2 = wlight.norm2();
			wlight.normalize();
			double tshadow;
			Vector Pshadow, Nshadow;
			int objectshadow;
			Ray rShadow(P + 0.001 * N, wlight);

			double l = I/(4 * M_PI * dlight2) * std::max(0.0,dot(N, wlight)); // color with direct light
			color = l * objects[objectID]->albedo / M_PI; // add the material in color
			if (intersect(rShadow, Pshadow, Nshadow, objectshadow, tshadow)){
				if (sqr(tshadow) < dlight2){
					color = Vector(0,0,0); 
				}
			}

			// add indirect light
			Vector wi = random_cos(N);
			Ray indirectRay(P + 0.001 * N, wi); 
			return objects[objectID]->albedo * getColor(indirectRay, bounce - 1) + color;
		}
		return color;
	}

	std::vector<Geometry*> objects; // list of objects in the scene
	Vector L; // light source
	double I; // light intensity
};

int main() {
	int W = 512;
	int H = 512;
	int nrays = 32;
	int bounces = 3;

	double fov = 60 * M_PI / 180;
	double z = -W / (2 * tan(fov / 2));

	for (int i = 0; i<8; i++){
		engine[i].seed(i);
	}

	Scene scene;
	// Spheres
	// transparent
	scene.addSphere(new Sphere(Vector(10,12,0), 5, Vector(1, 0.5, 0.3), false, true)); 
	// transparent hollow (does not work right now, but it works in my soft shadow script (without the soft shadow))
	//s.addSphere(new Sphere(Vector(-15,10,0), 5, Vector(1, 0.5, 0.3), false, true));
	//s.addSphere(new Sphere(Vector(-15,10,0), 4.5, Vector(1, 0.5, 0.3), false, true, true));
	// mirror (reflective)
	scene.addSphere(new Sphere(Vector(13,-5,20), 5, Vector(1, 0.5, 0.3), true, false));
	// Ceiling and floor
	scene.addSphere(new Sphere(Vector(0,1000,0), 940, Vector(1, 0, 0))); // Ceiling
	scene.addSphere(new Sphere(Vector(0,-1000,0), 990, Vector(0, 0, 1))); // Floor
	// Walls
	scene.addSphere(new Sphere(Vector(1000,0,0), 920, Vector(0.1, 0.5, 0.8))); // right
	scene.addSphere(new Sphere(Vector(-1000,0,0), 920, Vector(0.5, 0.9, 0.1))); // left
	scene.addSphere(new Sphere(Vector(0,0,1000), 940, Vector(0.7, 0.2, 0.5))); // back
	scene.addSphere(new Sphere(Vector(0,0,-1000), 940, Vector(0, 1, 0))); // front

	// Objects
	// cat
	TriangleMesh *cat = new TriangleMesh(Vector(0.3, 0.2, 0.25));
	cat->readOBJ("objects/cat.obj");
	cat->translate_and_scale(Vector(-15,-10,0), 0.5);
	cat->buildBVH(&cat->bvh, 0, cat->indices.size());
	scene.addMesh(cat);
	// stool
	TriangleMesh *stool = new TriangleMesh(Vector(0.6, 0.4, 0.2));
	stool->readOBJ("objects/stool.obj");
	stool->translate_and_scale(Vector(10,-2,0), 3);
	stool->buildBVH(&stool->bvh, 0, stool->indices.size());
	scene.addMesh(stool);

	Vector Cam(0,0,55); // camera point of view
	scene.L = Vector(-10, 20, 40); // Light source
	scene.I = 3E10; // Light intensity

	std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		int tid = omp_get_thread_num();
		for (int j = 0; j < W; j++) {

			Vector pixelColor;
			for (int k = 0; k < nrays; k++){

				// box Muller
				double u1 = uniform(engine[tid]);
				double u2 = uniform(engine[tid]);
				double r1 = cos(2 * M_PI * u1) * sqrt(-2 * log(u2)) * 0.4;
				double r2 = sin(2 * M_PI * u1) * sqrt(-2 * log(u2)) * 0.4;

				// vector of the direction of the view u from the camera
				Vector u(j - W / 2 + 0.5 + r1, H / 2 - i - 0.5 + r2, z);
				u.normalize();
				Ray r(Cam, u);
				pixelColor += scene.getColor(r, bounces)/nrays;		
			}

			image[(i * W + j) * 3 + 0] = std::min(255.0, std::pow(pixelColor[0], 0.45)); // including gamma correction
			image[(i * W + j) * 3 + 1] = std::min(255.0, std::pow(pixelColor[1], 0.45)); 
			image[(i * W + j) * 3 + 2] = std::min(255.0, std::pow(pixelColor[2], 0.45));
		}
	}
	stbi_write_png("images/image.png", W, H, 3, &image[0], 0);

	return 0;
}
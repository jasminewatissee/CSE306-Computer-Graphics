#include <vector>
#include <string>
#include <iostream>

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


// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {  
public:
	std::vector<Vector> vertices;
};	

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
	FILE* f = fopen(filename.c_str(), "w+"); 
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i=0; i<polygons.size(); i++) {
	    fprintf(f, "<g>\n");
	    fprintf(f, "<polygon points = \""); 
	    for (int j = 0; j < polygons[i].vertices.size(); j++) {
		    fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
	    }
	    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
    }
	fprintf(f, "</svg>\n");
	fclose(f);
}

//Sutherland Hodgman to clip cell by bissector of UV
Polygon clip_by_bissector(const Polygon &cell, const Vector &Xi, const Vector &Xj){
	Polygon result;
	int N = cell.vertices.size();
	for (int i=0; i< N; i++){
		const Vector& A = cell.vertices[i==0?(N-1): i-1];
		const Vector& B = cell.vertices[i];
		Vector M = (Xi + Xj) * 0.5;
		double t = dot(M-A, Xj - Xi) / dot(B-A, Xj - Xi);
		Vector P = A + t* (B-A);

		if ((B-Xi).norm2() <= (B-Xj).norm2()){ // B is inside
			if ((A-Xi).norm2() > (A-Xj).norm2()){ // A is outside
				result.vertices.push_back(P);
			}
			result.vertices.push_back(B);
		} else if ((A-Xi).norm2() <= (A-Xj).norm2()){ // A is inside
			result.vertices.push_back(P);
		}
	}
	return result;
}

class Voronoi{
public:
	Voronoi(){}

	void compute(){
		Polygon square;
		square.vertices.push_back(Vector(0,0,0)); // anticlockwise
		square.vertices.push_back(Vector(1,0,0));
		square.vertices.push_back(Vector(1,1,0));
		square.vertices.push_back(Vector(0,1,0));

		cells.resize(points.size());

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++){
			// cell i
			Polygon cell = square;
			for (int j = 0; j < points.size(); j++){
				if (i ==j) continue;
				cell = clip_by_bissector(cell, points[i], points[j]);
			}
			cells[i] = cell;
		}
	}

	std::vector<Vector> points;
	std::vector<Polygon> cells;
};

int main(){
	Voronoi vor;
	int N = 1000;
	vor.points.resize(N);
	for (int i = 0; i< N; i++){
		vor.points[i] = Vector(rand()/double(RAND_MAX),rand()/double(RAND_MAX),0);
	}
	vor.compute();
	
	save_svg(vor.cells, "test2.svg");

	return 0;
}
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <iostream>
#include <sstream>

#include <omp.h>
#include "lbfgs.h"

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
	double area(){
		if (vertices.size() <= 2) return 0;

		double result = 0;
		int n = vertices.size();
		for (int i = 0; i < n; i++){
			int j = (i+1)%n;
			result += vertices[i][0]*vertices[j][1] - vertices[j][0]*vertices[i][1];
		}
		result = std::abs(result)*0.5;
		return result;
	}

	double sumSquareDistance(const Vector& P){
		if (vertices.size() <= 2) return 0;
		if (area() == 0) return 0;

		double result = 0;
		int n = vertices.size();
		for (int i = 1; i < n-1; i++){
			Vector c[3] = {vertices[0], vertices[i], vertices[i+1]};
			double areaC = std::abs((c[2][1]-c[0][1])*(c[1][0]-c[0][0]) - (c[2][0]-c[0][0])*(c[1][1]-c[0][1]))*0.5;
			double s = 0;
			for (int k = 0; k < 3; k++){
				for (int l= k; l< 3; l++){
					s += dot(c[k]-P, c[l]-P);
				}
			}
			result += areaC/6 * s;
		}
		return result;
	}

	Vector centroid(){
		if (vertices.size() <= 2) return Vector(0,0,0);

		Vector result(0, 0, 0);
		int n = vertices.size();
		for (int i=0; i < n; i++){
			int j = (i+1) % n;
			double xy = vertices[i][0]* vertices[j][1] - vertices[j][0] * vertices[i][1]; 
			result += Vector((vertices[i][0] + vertices[j][0])*xy, (vertices[i][1] + vertices[j][1])*xy, 0);
		}
		double a = area();
		result = result / (6*a); 
		return result;
	}

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

int sgn(double x){
	if (x>0) return 1;
	if (x<0) return -1;
	return 0;
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int frameid = 0) {
		int W = 1000, H = 1000;
		int N = cells.size();
		std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < cells.size(); i++) {

			double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
			for (int j = 0; j < cells[i].vertices.size(); j++) {
				bminx = std::min(bminx, cells[i].vertices[j][0]);
				bminy = std::min(bminy, cells[i].vertices[j][1]);
				bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
				bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
			}
			bminx = std::min(W-1., std::max(0., W * bminx));
			bminy = std::min(H-1., std::max(0., H * bminy));
			bmaxx = std::max(W-1., std::max(0., W * bmaxx));
			bmaxy = std::max(H-1., std::max(0., H * bmaxy));

			for (int y = bminy; y < bmaxy; y++) {
				for (int x = bminx; x < bmaxx; x++) {
					int prevSign = 0;
					bool isInside = true;
					double mindistEdge = 1E9;
					for (int j = 0; j < cells[i].vertices.size(); j++) {
						double x0 = cells[i].vertices[j][0] * W;
						double y0 = cells[i].vertices[j][1] * H;
						double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
						double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
						double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
						int sign = sgn(det);
						if (prevSign == 0) prevSign = sign; else
							if (sign == 0) sign = prevSign; else
							if (sign != prevSign) {
								isInside = false;
								break;
							}
						prevSign = sign;
						double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
						double distEdge = std::abs(det)/ edgeLen;
						double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
						if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
						mindistEdge = std::min(mindistEdge, distEdge);
					}
					if (isInside) {
						if (i < N) {   // the N first particles may represent fluid, displayed in blue
							image[((H - y - 1)*W + x) * 3] = 0;
							image[((H - y - 1)*W + x) * 3 + 1] = 0;
							image[((H - y - 1)*W + x) * 3 + 2] = 255;
						}
						if (mindistEdge <= 2) {
							image[((H - y - 1)*W + x) * 3] = 0;
							image[((H - y - 1)*W + x) * 3 + 1] = 0;
							image[((H - y - 1)*W + x) * 3 + 2] = 0;
						}

					}
					
				}
			}
		}
		std::ostringstream os;
		os << filename << frameid << ".png";
		stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
	}

Polygon clip_by_line(const Polygon &cell, const Vector &u, const Vector &v){
	Vector normal(v[1]-u[1], u[0]-v[0],0);
	
	Polygon result;
	int N = cell.vertices.size();
	for (int i = 0; i < N; i++){
		const Vector& A = cell.vertices[i==0 ? (N-1): i-1];
		const Vector& B = cell.vertices[i];
		Vector P = A + (dot(u-A, normal)/ dot(B-A, normal)) * (B-A);

		if (dot(B-u, normal) <= 0){ // B is inside
			if (dot(A-u, normal) > 0){ // A is outside
				result.vertices.push_back(P);
			}
			result.vertices.push_back(B);
		} else if (dot(A-u, normal) <= 0){ // A is inside
			result.vertices.push_back(P);
		}
	}
	return result;
}

//Sutherland Hodgman to clip cell by bissector of UV
Polygon clip_by_bisector(const Polygon &cell, const Vector &Xi, const Vector &Xj, double& wi, double& wj){
	Vector M = (Xi + Xj) / 2;
	Vector offset = (wi - wj) / (2. * (Xi - Xj).norm2()) * (Xj - Xi);
	Vector Mprime = offset + M;

	Polygon result;
	int N = cell.vertices.size();
	for (int i = 0; i < N; i++){
		const Vector& A = cell.vertices[i==0? (N-1): i-1];
		const Vector& B = cell.vertices[i];
		double t = dot(Mprime - A, Xj - Xi) / dot(B - A, Xj - Xi);
		Vector P = A + t * (B-A);

		if ((B - Xi).norm2() - wi <= (B - Xj).norm2() - wj){ // B is inside
			if ((A - Xi).norm2() - wi > (A - Xj).norm2() - wj){ // A is outside
				result.vertices.push_back(P);
			}
			result.vertices.push_back(B);
		} else if ((A - Xi).norm2() - wi <= (A - Xj).norm2() - wj){ // A is inside
			result.vertices.push_back(P);
		}
	}
	return result;
}

class Voronoi{
public:
	Voronoi(){
		// Create a unit circle
		for (int i=0; i<NCircles; i++){
			double theta = i*2*M_PI / NCircles;
			circle[i] = Vector(cos(theta), sin(theta), 0);
		}
	}

	void compute(){
		Polygon square;
		square.vertices.push_back(Vector(0,0,0)); // anticlockwise
		square.vertices.push_back(Vector(1,0,0));
		square.vertices.push_back(Vector(1,1,0));
		square.vertices.push_back(Vector(0,1,0));

		cells.resize(points.size());

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++){
			Polygon cell = square;
			Polygon result;
			result.vertices.reserve(20);
			for (int j = 0; j < points.size(); j++){
				if (i == j) continue;
				result.vertices.clear();
				result = clip_by_bisector(cell, points[i], points[j], weights[i], weights[j]);
				std::swap(result, cell);
			}
			cells[i] = cell;
		}
	}

	void compute_circles(){
		Polygon square;
		square.vertices.push_back(Vector(0,0,0)); // anticlockwise
		square.vertices.push_back(Vector(1,0,0));
		square.vertices.push_back(Vector(1,1,0));
		square.vertices.push_back(Vector(0,1,0));

		cells.resize(points.size());

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++){
			Polygon cell = square;
			Polygon result, result2;
			result.vertices.reserve(20);
			for (int j = 0; j < points.size(); j++){
				if (i == j) continue;
				result.vertices.clear();
				result = clip_by_bisector(cell, points[i], points[j], weights[i], weights[j]);
				std::swap(result, cell);
			}
			for (int j = 0; j < NCircles; j++){
				result.vertices.clear();
				double radius = sqrt(weights[i] - w_air);
				Vector u = circle[j] * radius + points[i];
				Vector v = circle[(j+1)% NCircles]*radius + points[i];
				result = clip_by_line(cell, u, v);
				std::swap(result, cell);
			}
			cells[i] = cell;
		}
	}

	std::vector<Vector> points;
	std::vector<Polygon> cells;
	std::vector<double> weights;
	double w_air;
	static const int NCircles = 200;
	Vector circle[NCircles];
};

class SemiDiscreteOT {
public:

	void optimize(){
		int N = diagram.points.size();
		diagram.weights.resize(N, 1.);
		
		double objectivefct = -1;
		lbfgs_parameter_t param;
		lbfgs_parameter_init(&param);
		param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;

		std::vector<double> optimized_weights(N, 1.);
		int ret = lbfgs(N, &optimized_weights[0], &objectivefct, _evaluate, _progress, this, &param);
		std::cout << "L-BFGS optimization terminated with status code = " << ret<< std::endl;

		diagram.weights = optimized_weights;
		diagram.compute();
	}

	void optimize_fluid(){
		int N = diagram.points.size();
		diagram.weights.resize(N, 1.);
		diagram.w_air = 0;
		
		double objectivefct = -1;
		lbfgs_parameter_t param;
		lbfgs_parameter_init(&param);
		param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;

		std::vector<double> optimized_weights(N+1, 1.);
		optimized_weights[N] = 0.999999;
		int ret = lbfgs(N+1, &optimized_weights[0], &objectivefct, _evaluate_fluid, _progress_fluid, this, &param);
		std::cout << "L-BFGS optimization terminated with status code = " << ret<< std::endl;

		for(int i = 0; i < N; i++){
			diagram.weights[i] = optimized_weights[i];
		}
		diagram.w_air = optimized_weights[N];
		diagram.compute_circles();
	}

	Voronoi diagram;

protected:
	static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
		for (int i = 0; i < n; i++){
			diagram.weights[i] = x[i];
		}
		diagram.compute();

        lbfgsfloatval_t fx = 0.0;
		double sumSqrDists = 0;
		double sumAreaW = 0;
		double sumLambdaW = 0;

		double lambda = 1. / n;

        for (int i = 0; i < n;i++) {
			double area = diagram.cells[i].area();
			g[i] = area - lambda;

			sumAreaW += area * diagram.weights[i];
			sumLambdaW +=  lambda * diagram.weights[i];
			sumSqrDists += diagram.cells[i].sumSquareDistance(diagram.points[i]);
        }
		fx = sumAreaW - sumSqrDists - sumLambdaW;
		return fx;
    }

	static lbfgsfloatval_t _evaluate_fluid(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->evaluate_fluid(x, g, n, step);
    }

    lbfgsfloatval_t evaluate_fluid(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
		for (int i = 0; i < n-1; i++){
			diagram.weights[i] = x[i];
		}
		diagram.w_air = x[n-1];
		diagram.compute_circles();

        lbfgsfloatval_t fx = 0.0;
		double sumSqrDists = 0;
		double sumAreaW = 0;
		double sumLambdaW = 0;
		double sumArea = 0;

		double fluid = 0.4;
		double lambda = fluid / (n-1);
		double desired_air = 0.6;

        for (int i = 0; i < n-1; i++) {
			double area = diagram.cells[i].area();
			g[i] = area - lambda;

			sumArea += area;
			sumAreaW += area * diagram.weights[i];
			sumLambdaW +=  lambda * diagram.weights[i];
			sumSqrDists += diagram.cells[i].sumSquareDistance(diagram.points[i]);
        }

		double estimated_air = 1 - sumArea;
		g[n-1] = estimated_air - desired_air;

		fx = sumAreaW - sumSqrDists - sumLambdaW + diagram.w_air * g[n-1];
		return fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
		for (int i = 0; i < n; i++){
			diagram.weights[i] = x[i];
		}
		diagram.compute();
		
		double maxDiff = 0;
		for(int i = 0; i < n; i++){
			maxDiff = std::max(maxDiff, std::abs(diagram.cells[i].area() - (1. / n)));
		}

        printf("Iteration %d:\n", k);
        printf("  maxDiff = %f, fx = %f, \n", maxDiff * n, fx);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

static int _progress_fluid(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<SemiDiscreteOT*>(instance)->progress_fluid(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress_fluid(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
		for (int i = 0; i < n-1; i++){
			diagram.weights[i] = x[i];
		}
		diagram.w_air = x[n-1];
		diagram.compute_circles();
		
		double maxDiff = 0;
		for(int i = 0; i< n-1; i++){
			maxDiff = std::max(maxDiff, std::abs(diagram.cells[i].area() - (0.4 / (n-1))));
		}

        printf("Iteration %d:\n", k);
        printf("  maxDiff = %f, fx = %f, \n", maxDiff * (n-1), fx);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }
};

class Fluid {
public:
	Fluid() : epsilon(0.004), mass(200), gravity(Vector(0,-9.81, 0)) {}

	void time_step(double dt){
		ot.diagram.points = positions;
		ot.optimize_fluid();

		for (int i=0; i < positions.size(); i++){
			Vector ForceSpring = 1./sqr(epsilon)*(ot.diagram.cells[i].centroid() - positions[i]);
			Vector ForceGravity = mass * gravity;

			velocity[i] += dt/mass * (ForceSpring + ForceGravity);
			positions[i] += dt * velocity[i];
			positions[i][0] = std::min(1-1E-9, std::max(1E-9, positions[i][0]));
			positions[i][1] = std::min(1-1E-9, std::max(1E-9, positions[i][1]));
		}
	}

	void simulate(){
		int N = 100;
		positions.resize(N);
		velocity.resize(N, Vector(0,0,0));

		for (int i=0; i< N; i++){
			positions[i] = Vector(rand()/double(RAND_MAX),rand()/double(RAND_MAX),0);
		}

		for (int i=0; i<500; i++){
			time_step(0.005);
			save_frame(ot.diagram.cells, "results/animations/anim", i);
		}
	}
	
	SemiDiscreteOT ot;
	std::vector<Vector> positions;
	std::vector<Vector> velocity;
	double epsilon, mass;
	Vector gravity;

};

int main(int argc, char* argv[]){

	if (argc != 2){
		std::cerr << "Invalid command, please give your command as:" << std::endl;
        std::cerr << "./geo <opt>" << std::endl;
        std::cerr << "opt being:" << std::endl;
        std::cerr << "0 for a Simple Voronoi Diagram" << std::endl;
		std::cerr << "1 for a SemiDiscreteOT Diagram" << std::endl;
		std::cerr << "2 for a SemiDiscreteOT Diagram with circles" << std::endl;
		std::cerr << "3 for a Fluid Simulation" << std::endl;
        return 1;
	}

	auto start = std::chrono::high_resolution_clock::now();

	if (std::stoi(argv[1]) == 0){
		Voronoi vor;
		int N = 1000;
		vor.points.resize(N);
		vor.weights.resize(N, 1.);
		for (int i = 0; i< N; i++){
			vor.points[i] = Vector(rand()/double(RAND_MAX),rand()/double(RAND_MAX),0);
		}
		vor.compute();
		
		save_svg(vor.cells, "results/voronoi.svg");
	} else if (std::stoi(argv[1]) == 1){
		SemiDiscreteOT ot;
		int N = 200;
		ot.diagram.points.resize(N);
		for (int i = 0; i< N; i++){
			ot.diagram.points[i] = Vector(rand()/double(RAND_MAX),rand()/double(RAND_MAX),0);
		}
		ot.optimize();
		
		save_svg(ot.diagram.cells, "results/OT.svg");
	} else if (std::stoi(argv[1]) == 2){
		SemiDiscreteOT ot;
		int N = 200;
		ot.diagram.points.resize(N);
		for (int i = 0; i< N; i++){
			ot.diagram.points[i] = Vector(rand()/double(RAND_MAX),rand()/double(RAND_MAX),0);
		}
		ot.optimize_fluid();
		
		save_svg(ot.diagram.cells, "results/OT_circles.svg");
	} else if (std::stoi(argv[1]) == 3){
		Fluid fluid;
		fluid.simulate();
	} else {
		std::cerr << "Wrong option, please use 0 for Voronoi, 1 for OT or 2 for OT with circles and 3 for FluidSimulator" << std::endl;
        return 1;
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << "Execution time (seconds):" << std::endl;
	std::cout << duration.count()/1000/1000 << std::endl;

	return 0;
}
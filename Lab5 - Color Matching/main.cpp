#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
std::default_random_engine generator;
std::normal_distribution<double> N01(0.0, 1.0);

void sliced(double* src, double* tgt, int Npix){
	std::vector<std::pair<double, int> > projectionSrc(Npix);
	std::vector<double> projectionTgt(Npix);

	for (int iter=0; iter < 50; iter++){
		double xLine = N01(generator);
		double yLine = N01(generator);
		double zLine = N01(generator);
		double norm = sqrt(xLine * xLine + yLine * yLine + zLine * zLine);
		xLine /= norm;
		yLine /= norm;
		zLine /= norm;
		for (int i = 0; i < Npix; i++){
			projectionSrc[i] = std::pair<double, int>(src[i*3 + 0] * xLine + src[i*3 + 1] * yLine + src[i*3 + 2] * zLine, i);
			projectionTgt[i] = tgt[i*3 + 0] * xLine + tgt[i*3 + 1] * yLine + tgt[i*3 + 2] * zLine;
		}
		std::sort(projectionSrc.begin(), projectionSrc.end());
		std::sort(projectionTgt.begin(), projectionTgt.end());

		for (int i = 0; i < Npix; i++){
			src[projectionSrc[i].second * 3 + 0] += (projectionTgt[i] - projectionSrc[i].first) * xLine;
			src[projectionSrc[i].second * 3 + 1] += (projectionTgt[i] - projectionSrc[i].first) * yLine;
			src[projectionSrc[i].second * 3 + 2] += (projectionTgt[i] - projectionSrc[i].first) * zLine;
		}
	}
}

int main() {

	int W, H, C;
	int W2, H2, C2;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("imgA.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	unsigned char *target = stbi_load("redim.jpg",
                                 &W2,
                                 &H2,
                                 &C2,
                                 STBI_rgb);						 
	std::vector<double> image_double(W*H*3);
	std::vector<double> target_double(W2*H2*3);
	for (int i=0; i<W*H*3; i++)
		image_double[i] = image[i];
	for (int i=0; i<W2*H2*3; i++)
		target_double[i] = target[i];
	
	sliced(&image_double[0], &target_double[0], W * H);

	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < W2*H2*3; i++){
		image_result[i] = (unsigned char) std::max(0. , std::min(255., image_double[i]));
	}

	stbi_write_png("image.png", W, H, 3, &image_result[0], 0);

	return 0;
}
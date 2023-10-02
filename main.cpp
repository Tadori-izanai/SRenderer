#include <vector>
#include <cmath>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const TGAColor blue  = TGAColor(0,   0,   255, 255);


// Draws a red point
void demo0() {
	TGAImage image(100, 100, TGAImage::RGB);
	image.set(52, 41, red);
	image.set(41, 52, red);
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
}

void lineNaive(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c) {
	for (float t = 0; t < 1; t += .1) {
		int x = x0 + (x1 - x0) * t;
		int y = y0 + (y1 - y0) * t;
		img.set(x, y, c);
	}
}

void line1Helper(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c, bool isSteep) {
	if (x0 > x1) {
		line1Helper(img, x1, y1, x0, y0, c, isSteep);
	}
	for (int x = x0; x <= x1; x += 1) {
		float t = (float) (x - x0) / (x1 - x0);
		int y = y0 + (y1 - y0) * t;
		if (isSteep) {
			img.set(y, x, c);
		} else {
			img.set(x, y, c);
		}
	}
}

void line1(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c) {
	if (std::abs(x1 - x0) < std::abs(y1 - y0)) {
		line1Helper(img, y0, x0, y1, x1, c, true);
	} else {
		line1Helper(img, x0, y0, x1, y1, c, false);
	}
}

void line2(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c) {
	bool isSteep = false;
	if (std::abs(x1 - x0) < std::abs(y1 - y0)) {
		isSteep = true;
		std::swap(x0, y0);
		std::swap(x1, y1);
	}
	if (x1 < x0) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}

	int dx = x1 - x0;
	int dy = y1 - y0;
	// float derror = std::abs(dy / float(dx));
	int derror2 = std::abs(dy) * 2;		// consider comparing `error * (2 * dx)` with `0.5 * (2 * dx)`
	int yStep = (dy > 0) ? 1 : -1;
	
	// float error = 0;
	int error2 = 0;
	int y = y0;
	for (int x = x0; x <= x1; x += 1) {
		if (isSteep) {
			img.set(y, x, c);
		} else {
			img.set(x, y, c);
		}
		// error += derror;
		// if (error > .5) {
		// 	y += yStep;
		// 	error -= 1;
		// }
		error2 += derror2;
		if (error2 > dx) {
			y += yStep;
			error2 -= dx * 2;
		}
	}
}

void line(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c) {
	bool isSteep = false;
	if (std::abs(x1 - x0) < std::abs(y1 - y0)) {
		isSteep = true;
		std::swap(x0, y0);
		std::swap(x1, y1);
	}
	if (x1 < x0) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}

	int dx = x1 - x0;
	int dy = y1 - y0;
	int derror2 = std::abs(dy) * 2;
	int yStep = (dy > 0) ? 1 : -1;
	
	int error2 = 0;
	int y = y0;
	if (isSteep) {
		for (int x = x0; x <= x1; x += 1) {
			img.set(y, x, c);
			error2 += derror2;
			if (error2 > dx) {
				y += yStep;
				error2 -= dx * 2;
			}
		}
	} else {
		for (int x = x0; x <= x1; x += 1) {
			img.set(x, y, c);
			error2 += derror2;
			if (error2 > dx) {
				y += yStep;
				error2 -= dx * 2;
			}
		}
	}
}

void line(TGAImage &img, const Vec2i &v0, const Vec2i &v1, const TGAColor &c) {
	line(img, v0.x, v0.y, v1.x, v1.y, c);
}

// Draws two lines
void demo1() {
	TGAImage image(100, 100, TGAImage::RGB);
	
	line(image, 11, 15, 77, 33, white);
	line(image, 33, 77, 15, 11, white);
	
	image.flip_vertically();
	image.write_tga_file("output.tga");
}

// Wireframe rendering
void demo2() {
	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("./obj/african_head.obj");
	
	int n = model.nfaces();
	for (int i = 0; i < n; i += 1) {
		std::vector<int> f = model.face(i);
		for (int j = 0; j < 3; j += 1) {
			Vec3f v0 = model.vert(f[j]);
			Vec3f v1 = model.vert(f[(j + 1) % 3]);
			int x0 = (v0.x + 1.0f) * 0.5f * width;
			int y0 = (v0.y + 1.0f) * 0.5f * height;
			int x1 = (v1.x + 1.0f) * 0.5f * width;
			int y1 = (v1.y + 1.0f) * 0.5f * height;
			line(canvas, x0, y0, x1, y1, white);
		}
	}
	canvas.flip_vertically();
	canvas.write_tga_file("output/wireframe.tga");
}

void triangleBorder(TGAImage &img, const Vec2i &v0, const Vec2i &v1, const Vec2i &v2, const TGAColor &c) {
	line(img, v0, v1, c);
	line(img, v1, v2, c);
	line(img, v2, v0, c);
}

void sortTriangleVertices(Vec2i &v0, Vec2i &v1, Vec2i &v2) {
	if (v1.y < v0.y) {
		std::swap(v0, v1);
	}
	if (v2.y < v1.y) {
		std::swap(v1, v2);
	}
	if (v1.y < v0.y) {
		std::swap(v0, v1);
	}
}

inline void horizontalLine(TGAImage &img, int y, int l, int r, const TGAColor &c) {
	if (l > r) {
		horizontalLine(img, y, r, l, c);
		return;
	}
	for (int x = l; x <= r; x += 1) {
		img.set(x, y, c);
	}
}

void triangle(TGAImage &img, Vec2i v0, Vec2i v1, Vec2i v2, const TGAColor &c) {
	if (v0.y == v1.y && v0.y == v1.y) {
		return;
	}
	sortTriangleVertices(v0, v1, v2);

	int totalHeight = v2.y - v0.y;
	int segmentHeight = v1.y - v0.y;
	int dx0 = v2.x - v0.x;
	int dx1 = v1.x - v0.x;
	for (int y = v0.y; y < v1.y; y += 1) {
		int x1 = v0.x + dx0 * float(y - v0.y) / totalHeight;
		int x2 = v0.x + dx1 * float(y - v0.y) / segmentHeight;
		horizontalLine(img, y, x1, x2, c);
	}
	segmentHeight = v2.y - v1.y;
	dx1 = v2.x - v1.x;
	for (int y = v1.y; y <= v2.y; y += 1) {
		int x1 = v0.x + dx0 * float(y - v0.y) / totalHeight;
		int x2 = v1.x + dx1 * float(y - v1.y) / segmentHeight;
		horizontalLine(img, y, x1, x2, c);
	}
}

Vec3f barycentric(const Vec2i pts[], const Vec2i &p) {
	Vec2i ab = pts[1] - pts[0];
	Vec2i ac = pts[2] - pts[0];
	Vec2i pa = pts[0] - p;
	Vec3i crossProd = Vec3i(ab.x, ac.x, pa.x) ^ Vec3i(ab.y, ac.y, pa.y);
	if (std::abs(crossProd.z) < 1) {		// `abs(u[2])` < 1 means `u[2]` is 0
		return Vec3f(-1, 1, 1);		// that means triangle is degenerate,
									// in this case return something with negative coordinates
	}

	float u = (float) crossProd.x / crossProd.z;
	float v = (float) crossProd.y / crossProd.z;
	return Vec3f(1.0f - u - v, u, v);
}

Vec3f barycentric(const Vec3f pts[], const Vec2i &p) {
	Vec2i ptsProjected[3];
	for (int i = 0; i < 3; i += 1) {
		ptsProjected[i].x = pts[i].x;
		ptsProjected[i].y = pts[i].y;
	}
	return barycentric(ptsProjected, p);
}

void triangle(TGAImage &img, const Vec2i pts[], const TGAColor &c) {
	const float eps = 0.001;
	int xLeft = std::numeric_limits<int>::max();
	int xRight = std::numeric_limits<int>::min();
	int yBottom = std::numeric_limits<int>::max();
	int yTop = std::numeric_limits<int>::min();
	for (int i = 0; i < 3; i += 1) {
		xLeft = std::min(xLeft, pts[i].x);
		xRight = std::max(xRight, pts[i].x);
		yBottom = std::min(yBottom, pts[i].y);
		yTop = std::max(yTop, pts[i].y);
	}

	for (int y = yBottom; y <= yTop; y += 1) {
		for (int x = xLeft; x <= xRight; x += 1) {
			Vec3f bc = barycentric(pts, Vec2i(x, y));
			if (bc.x < -eps || bc.y < -eps || bc.z < -eps) {
				continue;
			}
			img.set(x, y, c);
		}
	}
}

/** 
 * Draws a triangle considering the z-buffer.
 * The x and y coordinates of pts[i] is in the screen space ([0, width) x [0, height)),
 * and the z coordinate of pts[i] is in the world space. (i = 0, 1, 2)
 */
void triangle(TGAImage &img, float *zBuffer, const Vec3f pts[], const TGAColor &c) {
	const float eps = 0.001;
	int xLeft = std::numeric_limits<int>::max();
	int xRight = std::numeric_limits<int>::min();
	int yBottom = std::numeric_limits<int>::max();
	int yTop = std::numeric_limits<int>::min();
	for (int i = 0; i < 3; i += 1) {
		xLeft = std::min(xLeft, (int) pts[i].x);
		xRight = std::max(xRight, (int) pts[i].x);
		yBottom = std::min(yBottom, (int) pts[i].y);
		yTop = std::max(yTop, (int) pts[i].y);
	}
	for (int y = yBottom; y <= yTop; y += 1) {
		for (int x = xLeft; x <= xRight; x += 1) {
			Vec3f bc = barycentric(pts, Vec2i(x, y));
			if (bc.x < -eps || bc.y < -eps || bc.z < -eps) {
				continue;
			}
			
			float zWorldCoord = 0;
			for (int i = 0; i < 3; i += 1) {
				zWorldCoord += bc[i] * pts[i].z;
			}
			if (zWorldCoord > zBuffer[x + img.get_width() * y]) {
				zBuffer[x + img.get_width() * y] = zWorldCoord;
				img.set(x, y, c);
			}
		}
	}
}

// Draws triangles
void demo3() {
	TGAImage image(200, 200, TGAImage::RGB);
	Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
	Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
	Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
	
	// triangleBorder(image, t0[0], t0[1], t0[2], blue); 
	// triangleBorder(image, t1[0], t1[1], t1[2], blue); 
	// triangleBorder(image, t2[0], t2[1], t2[2], blue);
	triangle(image, t0[0], t0[1], t0[2], red); 
	triangle(image, t1[0], t1[1], t1[2], white); 
	triangle(image, t2[0], t2[1], t2[2], green);

	Vec2i pts[3] = {Vec2i(10,10), Vec2i(100, 30), Vec2i(190, 160)}; 
    triangle(image, pts, blue); 
	
	image.flip_vertically();
	image.write_tga_file("output.tga");
}

// Fills with random colors
void demo4() {
	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("./obj/african_head.obj");
	
	int n = model.nfaces();
	for (int i = 0; i < n; i += 1) {
		std::vector<int> face = model.face(i);	// vertex indices
		Vec2i screenCoords[3];
		for (int j = 0; j < 3; j += 1) {
			Vec3f v = model.vert(face[j]);
			int x = (v.x + 1.f) * .5f * width;
			int y = (v.y + 1.f) * .5f * height;
			screenCoords[j] = Vec2i(x, y);
		}
		TGAColor randColor(rand() % 255, rand() % 255, rand() % 255, 255);
		triangle(canvas, screenCoords, randColor);
	}
	canvas.flip_vertically();
	canvas.write_tga_file("output/random_colors.tga");
}

// Do back-face culling
void demo5() {
	Vec3f lightDir(0, 0, -1);

	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("./obj/african_head.obj");
	
	int n = model.nfaces();
	for (int i = 0; i < n; i += 1) {
		std::vector<int> face = model.face(i);	// vertex indices
		Vec2i screenCoords[3];
		Vec3f worldCoords[3];
		for (int j = 0; j < 3; j += 1) {
			Vec3f v = model.vert(face[j]);
			int x = (v.x + 1.f) * .5f * width;
			int y = (v.y + 1.f) * .5f * height;
			screenCoords[j] = Vec2i(x, y);
			worldCoords[j] = v;
		}
		Vec3f n = (worldCoords[1] - worldCoords[0]) ^ (worldCoords[2] - worldCoords[0]);
		n.normalize();
		float intensity = -(lightDir * n);
		if (intensity > 0) {
			triangle(canvas, screenCoords, white * intensity);
		}
	}

	canvas.flip_vertically();
	canvas.write_tga_file("output/backface.tga");
}

/**
 * Converts the world coordinate to the screen.
 * Lets (x, y) converts from (-1, 1) x (-1, 1) to (0, width) x (0, height)
 * Remains z constant.
*/
Vec3f worldToScreen(const Vec3f &worldCoord, int width, int height) {
	int x = (worldCoord.x + 1.f) * .5f * width;
	int y = (worldCoord.y + 1.f) * .5f * height;
	return Vec3f(x, y, worldCoord.z);
}

// Uses z-buffer
void demoe6() {
	Vec3f lightDir(0, 0, -1);

	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("./obj/african_head.obj");

	int n = model.nfaces();
	float *zBuffer = new float[width * height];
	std::fill(zBuffer, zBuffer + width * height, -1);

	for (int i = 0; i < n; i += 1) {
		std::vector<int> face = model.face(i);
		Vec3f screenCoords[3];
		Vec3f worldCoords[3];
		for (int j = 0; j < 3; j += 1) {
			Vec3f v = model.vert(face[j]);
			screenCoords[j] = worldToScreen(v, width, height);
			worldCoords[j] = v;
		}
		Vec3f n = (worldCoords[1] - worldCoords[0]) ^ (worldCoords[2] - worldCoords[0]);
		n.normalize();
		float intensity = -(lightDir * n);
		triangle(canvas, zBuffer ,screenCoords, white * intensity);
	}	

	canvas.flip_vertically();
	canvas.write_tga_file("output/zbuffer.tga");	
}

// z-buffer with texture
void demo7() {
	
}

void test();
int main(int argc, char** argv) {
	demoe6();
	return 0;
}

void test() {
	Vec3f v(1, 2, 3);
	std::cout << v[1] << std::endl;
}

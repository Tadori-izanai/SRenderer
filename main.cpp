#include <vector>
#include <cmath>
#include <iostream>

#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "gl.h"

// Draws a RED point
void demo0() {
	TGAImage image(100, 100, TGAImage::RGB);
	image.set(52, 41, RED);
	image.set(41, 52, RED);
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
}

// Draws two lines
void demo1() {
	TGAImage image(100, 100, TGAImage::RGB);
	
	line(image, 11, 15, 77, 33, WHITE);
	line(image, 33, 77, 15, 11, WHITE);
	
	image.flip_vertically();
	image.write_tga_file("output.tga");
}

// Wireframe rendering
void demo2() {
	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("../obj/african_head.obj");
	
	int n = model.nfaces();
	for (int i = 0; i < n; i += 1) {
		std::vector<int> f = model.facePositions(i);
		for (int j = 0; j < 3; j += 1) {
			Vec3f v0 = model.vert(f[j]);
			Vec3f v1 = model.vert(f[(j + 1) % 3]);
			int x0 = (v0.x + 1.0f) * 0.5f * width;
			int y0 = (v0.y + 1.0f) * 0.5f * height;
			int x1 = (v1.x + 1.0f) * 0.5f * width;
			int y1 = (v1.y + 1.0f) * 0.5f * height;
			line(canvas, x0, y0, x1, y1, WHITE);
		}
	}
	canvas.flip_vertically();
	canvas.write_tga_file("output/wireframe.tga");
}

// Draws triangles
void demo3() {
	TGAImage image(200, 200, TGAImage::RGB);
	Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
	Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
	Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
	
	// triangleBorder(image, t0[0], t0[1], t0[2], BLUE);
	// triangleBorder(image, t1[0], t1[1], t1[2], BLUE);
	// triangleBorder(image, t2[0], t2[1], t2[2], BLUE);
	triangle(image, t0[0], t0[1], t0[2], RED);
	triangle(image, t1[0], t1[1], t1[2], WHITE);
	triangle(image, t2[0], t2[1], t2[2], GREEN);

	Vec2i pts[3] = {Vec2i(10,10), Vec2i(100, 30), Vec2i(190, 160)}; 
    triangle(image, pts, BLUE);
	
	image.flip_vertically();
	image.write_tga_file("output.tga");
}

// Fills with random colors
void demo4() {
	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("../obj/african_head.obj");
	
	int n = model.nfaces();
	for (int i = 0; i < n; i += 1) {
		std::vector<int> face = model.facePositions(i);	// vertex indices
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

// Do back-facePositions culling
void demo5() {
	Vec3f lightDir(0, 0, -1);

	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("../obj/african_head.obj");
	
	int n = model.nfaces();
	for (int i = 0; i < n; i += 1) {
		std::vector<int> face = model.facePositions(i);	// vertex indices
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
			triangle(canvas, screenCoords, WHITE * intensity);
		}
	}

	canvas.flip_vertically();
	canvas.write_tga_file("output/backface.tga");
}

// Uses z-buffer
void demo6() {
	Vec3f lightDir(0, 0, -1);

	int width = 800;
	int height = 800;
	TGAImage canvas(width, height, TGAImage::RGB);
	Model model("../obj/african_head.obj");

	int n = model.nfaces();
	float *zBuffer = new float[width * height];
	std::fill(zBuffer, zBuffer + width * height, -1);

	for (int i = 0; i < n; i += 1) {
		std::vector<int> face = model.facePositions(i);
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
		triangle(canvas, zBuffer, screenCoords, WHITE * intensity);
	}	

	canvas.flip_vertically();
	canvas.write_tga_file("output/zbuffer.tga");	
}

// z-buffer with texture
void demo7() {
    Vec3f lightDir(0, 0, -1);

    int width = 800;
    int height = 800;
    TGAImage canvas(width, height, TGAImage::RGB);
    TGAImage texture;
    if (!texture.read_tga_file("../obj/african_head_diffuse.tga")) {
        std::cout << "Failed to read ../obj/african_head_diffuse.tga" << std::endl;
        return;
    }
    Model model("../obj/african_head.obj");

    int n = model.nfaces();
    float *zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, -1);

    for (int i = 0; i < n; i += 1) {
        std::vector<int> facePos = model.facePositions(i);
        std::vector<int> faceTex = model.faceTextures(i);
        std::vector<int> faceNor = model.faceNormals(i);

        Vec3f screenCoords[3];
        Vec2f textureCoords[3];
        Vec3f normals[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = worldToScreen(model.vert(facePos[j]), width, height);
            textureCoords[j] = model.uv(faceTex[j]);
            normals[j] = model.normal(faceNor[j]);
        }

        triangle(canvas, texture, zBuffer, lightDir, screenCoords, textureCoords, normals);
    }

    canvas.flip_vertically();
    canvas.write_tga_file("output/zbuffer_with_tex.tga");
}

// Uses projection and viewport matrices
void demo8() {
    Vec3f lightDir(0, 0, -1);
    Vec3f camera(0, 0, 3);
    int width = 800;
    int height = 800;
    Matrix projection = getProjection(camera.z);
    Matrix viewport = getViewport(width, height);
    Matrix projectionForNormals = projection.inverse().transpose();

//    std::cout << projection << std::endl;
//    std::cout << viewport << std::endl;
//    std::cout << projectionForNormals << std::endl;

    TGAImage canvas(width, height, TGAImage::RGB);
    TGAImage texture;
    if (!texture.read_tga_file("../obj/african_head_diffuse.tga")) {
        std::cout << "Failed to read ../obj/african_head_diffuse.tga" << std::endl;
        return;
    }
    Model model("../obj/african_head.obj");

    int n = model.nfaces();
    float *zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, -1);

    for (int i = 0; i < n; i += 1) {
        std::vector<int> facePos = model.facePositions(i);
        std::vector<int> faceTex = model.faceTextures(i);
        std::vector<int> faceNor = model.faceNormals(i);

        Vec3f screenCoords[3];
        Vec2f textureCoords[3];
        Vec3f normals[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = mat2vec(
                    viewport * projection * vec2mat(model.vert(facePos[j]))
                    );
            textureCoords[j] = model.uv(faceTex[j]);
//            normals[j] = model.normal(faceNor[j]);
            normals[j] = vecTrans(projectionForNormals, model.normal(faceNor[j]));  // actually no change
        }

        triangle(canvas, texture, zBuffer, lightDir, screenCoords, textureCoords, normals);
    }

    canvas.flip_vertically();
    canvas.write_tga_file("output/zbuffer_with_tex_proj.tga");
}

// do mvp and viewport
void demo9() {
    const int width = 800;
    const int height = 800;
    const int depth = 255;

    TGAImage canvas(width, height, TGAImage::RGB);
    TGAImage texture;
    if (!texture.read_tga_file("../obj/african_head_diffuse.tga")) {
        std::cout << "Failed to read ../obj/african_head_diffuse.tga" << std::endl;
        return;
    }
    Model model("../obj/african_head.obj");

    Vec3f lightDir = Vec3f(-1, 1, -1).normalize();
    Vec3f eye(1, 1, 3);
    Vec3f center(0, 0, 0);
    Vec3f up(0, 1, 0);

    Matrix modelView = lookAt(eye, center, up);
    Matrix projection = getProjection((eye - center).norm());
    Matrix viewport = getViewport(width * 3 / 4, height * 3 / 4, depth, width / 8, height / 8);
    Matrix totTrans = viewport * projection * modelView;
    std::cout << totTrans << std::endl;

    int n = model.nfaces();
    float *zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, -1);

    for (int i = 0; i < n; i += 1) {
        std::vector<int> facePos = model.facePositions(i);
        std::vector<int> faceTex = model.faceTextures(i);
        std::vector<int> faceNor = model.faceNormals(i);

        Vec3f screenCoords[3];
        Vec2f textureCoords[3];
        Vec3f normals[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = mat2vec(
                    totTrans * vec2mat(model.vert(facePos[j]))
            );
            textureCoords[j] = model.uv(faceTex[j]);
            normals[j] = model.normal(faceNor[j]);
        }
        triangle(canvas, texture, zBuffer, lightDir, screenCoords, textureCoords, normals);
    }

    canvas.flip_vertically();
    canvas.write_tga_file("output/mvp.tga");
}

void test0() {
    Vec3f v(1, 2, 3);
    Vec4f vv(v, Vec4<float>::POINT);
    std::cout << v << vv << std::endl;
    std::cout << std::numeric_limits<float>::min() << std::endl;
    std::cout << std::numeric_limits<float>::lowest() << std::endl;
    std::cout << std::numeric_limits<float>::max() << std::endl;

    Matrix M = Matrix::identity(4);
    std::cout << M.multiply(v, Vec4<float>::VECTOR) << std::endl;
    std::cout << M.multiply(v, Vec4<float>::POINT) << std::endl;
}

void test() {
//    std::cout << int(1.1) << std::endl;
//    std::cout << int(1.5) << std::endl;
//    std::cout << int(1.6) << std::endl;

    Vec3f lightDir(0, 0, -1);
    int width = 800;
    int height = 800;
    TGAImage canvas(width, height, TGAImage::RGB);
    TGAImage texture;
    if (!texture.read_tga_file("../obj/african_head_diffuse.tga")) {
        std::cout << "Failed to read ../obj/african_head_diffuse.tga" << std::endl;
        return;
    }
    Mesh model("../obj/african_head.obj");

    int n = model.nFaces();
    float *zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, -1);

    for (int i = 0; i < n; i += 1) {
        Vec3f screenCoords[3];
        Vec2f textureCoords[3];
        Vec3f normals[3];
        for (int j = 0; j < 3; j += 1) {
            Vertex curr = model.getVertex(i, j);
            screenCoords[j] = worldToScreen(curr.position, width, height);
            textureCoords[j] = curr.uv;
            normals[j] = curr.normal;
        }

        triangle(canvas, texture, zBuffer, lightDir, screenCoords, textureCoords, normals);
    }

    canvas.flip_vertically();
    canvas.write_tga_file("output/test.tga");
}


//////////////////////////////////////////////////////////////////////////////////////////////////

Mesh *model = nullptr;
const int width = 800;
const int height = 800;
const int depth = 255;

Vec3f lightDir(-1, -1, -1);
Vec3f eyePos(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);


class PhongShader : public IShader {
    Matrix uniformMVP;          // projection * view * model
    Matrix uniformMVPIT;        // inverse transpose of (projection * view * model)
    Vec2f varyingUVs[3];

public:
    explicit PhongShader(Matrix mvp) : uniformMVP(mvp) {
        uniformMVPIT = mvp.inverse().transpose();
    }

    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;

        varyingUVs[nthVert] = curr.uv;
        return gl_Vertex;
    }

    bool fragment(Vec3f bary, TGAColor &color) override {
        Vec2f uv = varyingUVs[0] * bary[0] + varyingUVs[1] * bary[1] + varyingUVs[2] * bary[2];

        Vec3f n = model->getNormal(uv);
        n = uniformMVPIT.multiply(n, Vec4<float>::VECTOR).normalize();
        Vec3f l = -uniformMVP.multiply(lightDir, Vec4<float>::VECTOR).normalize();
        Vec3f r = (n * (n * l * 2.f) - l).normalize();

        int ambient = 5;
        float diffuse = std::max(l * n, 0.f);
        float specular = std::pow(std::max(r.z, 0.f), model->getSpecular(uv));
        TGAColor c = model->getDiffuse(uv);
        for (int i = 0; i < 3; i += 1) {
            color[i] = std::min(ambient + int(c[i] * (diffuse + specular * .6f)), 255);
        }
        return false;
    }
};

void demo() {
    model = new Mesh("../obj/african_head.obj", "../obj/african_head_diffuse.tga",
                     "../obj/african_head_nm.tga", "../obj/african_head_spec.tga");
    lightDir.normalize();

    modelView = lookAt(eyePos, center, up);
    projection = getProjection((eyePos - center).norm());
    viewport = getViewport(width * 3 / 4, height * 3 / 4, depth, width / 8, height / 8);

    TGAImage canvas(width, height, TGAImage::RGB);
    auto *zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, std::numeric_limits<float>::lowest());

    PhongShader shader(projection * modelView);
    for (int i = 0; i < model->nFaces(); i += 1) {
        Vec4f screenCoords[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = shader.vertex(i, j);
        }
        triangle(screenCoords, shader, canvas, zBuffer);
    }

    canvas.flip_vertically();
    canvas.write_tga_file("../output/output.tga");
    delete[] zBuffer;
    delete model;
}

int main(int argc, char** argv) {
	demo();
	return 0;
}

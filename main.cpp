#include <utility>
#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>

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

void test() {
//    Vec3f v(1, 2, 3);
//    Vec4f vv(v, Vec4<float>::POINT);
//    std::cout << v << vv << std::endl;
//    std::cout << std::numeric_limits<float>::min() << std::endl;
//    std::cout << std::numeric_limits<float>::lowest() << std::endl;
//    std::cout << std::numeric_limits<float>::max() << std::endl;

    auto points = getRandomPointsOnHemisphere(1000);
    for (const auto &p : points) {
        std::cout << p;
    }

    Vec2f v(3, 4);
    std::cout << v << std::endl;
    std::cout << v.normalize() << std::endl;
    std::cout << v << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Mesh *model = nullptr;
const int width = 800;
const int height = 800;
const int depth = 255;

Vec3f lightDir(-1, -1, 0);
Vec3f eyePos(1, 1, 4);
//Vec3f eyePos(1, -1, 4);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

float *shadowDepth;

class DepthShader : public IShader {
    Matrix uniformMVP;
    Vec3f varyingPositions[3];

public:
    explicit DepthShader(Matrix mvp) : uniformMVP(std::move(mvp)) {}

    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;
        varyingPositions[nthVert] = gl_Vertex.proj3();
        return gl_Vertex;
    }

    bool fragment(Vec3f bary, TGAColor &color) override {
        float z = 0;
        for (int i = 0; i < 3; i += 1) {
            z += varyingPositions[i].z * bary[i];
        }
        color = WHITE * (z / depth);
        return false;
    }
};

class DepthShaderLite : public IShader {
    Matrix uniformMVP;

public:
    explicit DepthShaderLite(Matrix mvp) : uniformMVP(std::move(mvp)) {}

    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;
        return gl_Vertex;
    }

    bool fragment(Vec3f bary, TGAColor &color) override {
        return false;
    }
};

class OcclusionShader: public IShader {
    Matrix uniformMVP;
    Matrix uniformMVPLight;
    Vec3f varyingPositions[3];
    Vec2f varyingUVs[3];
    TGAImage *occ;
    float *depthBuffer;

public:
    explicit OcclusionShader(Matrix mvp, Matrix mvpLight, TGAImage *occ, float *depthBuffer)
        : uniformMVP(std::move(mvp)), uniformMVPLight(std::move(mvpLight)), occ(occ), depthBuffer(depthBuffer) {}

    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;

        varyingPositions[nthVert] = curr.position;
        varyingUVs[nthVert] = curr.uv;
        return gl_Vertex;
    }

    bool fragment(Vec3f bary, TGAColor &color) override {
        const float eps = 11.11;
        Vec3f positionInSM = (viewport * uniformMVPLight).multiply(
            varyingPositions[0] * bary[0] + varyingPositions[1] * bary[1] + varyingPositions[2] * bary[2],
            Vec4<float>::POINT
        );
        if (positionInSM.z > depthBuffer[int(positionInSM.x) + width * int(positionInSM.y)] - eps) {
            Vec2f uv = varyingUVs[0] * bary[0] + varyingUVs[1] * bary[1] + varyingUVs[2] * bary[2];
            occ->set(occ->get_width() * uv.u, occ->get_height() * uv.v, WHITE);
        }
        return false;
    }
};

class DiffuseShader : public IShader {
    Matrix uniformMVP;
    Vec2f varyingUVs[3];
public:
    explicit DiffuseShader(Matrix mvp) : uniformMVP(std::move(mvp)) {}
    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;
        varyingUVs[nthVert] = curr.uv;
        return gl_Vertex;
    }
    bool fragment(Vec3f bary, TGAColor &color) override {
        Vec2f uv = varyingUVs[0] * bary[0] + varyingUVs[1] * bary[1] + varyingUVs[2] * bary[2];
        color = model->getDiffuse(uv);
        return false;
    }
};

class Shader : public IShader {
    Matrix uniformMVP;          // projection * view * model
    Matrix uniformMVPIT;        // inverse transpose of (projection * view * model)
    Matrix uniformMVPLight;
    Vec2f varyingUVs[3];
    Vec3f varyingPositions[3];

public:
    explicit Shader(Matrix mvp, Matrix mvpLight) : uniformMVP(mvp), uniformMVPLight(mvpLight) {
        uniformMVPIT = mvp.inverse().transpose();
    }

    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;

        varyingPositions[nthVert] = curr.position;
        varyingUVs[nthVert] = curr.uv;
        return gl_Vertex;
    }

    bool fragment(Vec3f bary, TGAColor &color) override {
        const float eps = 10;
        Vec3f positionInSM = (viewport * uniformMVPLight).multiply(
            varyingPositions[0] * bary[0] + varyingPositions[1] * bary[1] + varyingPositions[2] * bary[2],
            Vec4<float>::POINT
        );
        float shadow = .4f;
        if (positionInSM.z > shadowDepth[int(positionInSM.x) + width * int(positionInSM.y)] - eps) {
            shadow = 1.f;
        }

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
            color[i] = std::min(ambient + int(c[i] * (diffuse + specular * .6f) * shadow), 255);
        }
        return false;
    }
};

/** used for tangent space normal mapping */
class TGShader : public IShader {
    Matrix uniformMVP;          // projection * view * model
    Matrix uniformMVPIT;        // inverse transpose of (projection * view * model)
    Vec3f varyingPositions[3];
    Vec2f varyingUVs[3];
    Vec3f varyingNormals[3];

public:
    explicit TGShader(Matrix mvp) : uniformMVP(mvp) {
        uniformMVPIT = mvp.inverse().transpose();
    }

    Vec4f vertex(int iFace, int nthVert) override {
        Vertex curr = model->getVertex(iFace, nthVert);
        varyingPositions[nthVert] = uniformMVP.multiply(curr.position, Vec4<float>::POINT);
        varyingUVs[nthVert] = curr.uv;
        varyingNormals[nthVert] = uniformMVPIT.multiply(curr.normal, Vec4<float>::VECTOR);

        Vec4f gl_Vertex(curr.position, Vec4<float>::POINT);     // gl_Vertex
        gl_Vertex = viewport * uniformMVP * gl_Vertex;
        return gl_Vertex;
    }

    bool fragment(Vec3f bary, TGAColor &color) override {
        Vec2f uv = varyingUVs[0] * bary[0] + varyingUVs[1] * bary[1] + varyingUVs[2] * bary[2];
        Vec3f N = (
            varyingNormals[0] * bary[0] + varyingNormals[1] * bary[1] + varyingNormals[2] * bary[2]
        ).normalize();

        Matrix A(3, 3);
        A.assignRow(0, varyingPositions[1] - varyingPositions[0]);
        A.assignRow(1, varyingPositions[2] - varyingPositions[0]);
        A.assignRow(2, N);
        Matrix invA = A.inverse();
        Vec2f deltaUV1 = varyingUVs[1] - varyingUVs[0];
        Vec2f deltaUV2 = varyingUVs[2] - varyingUVs[0];
        Vec3f T = (invA * Vec3f(deltaUV1.u, deltaUV2.u, 0.f)).normalize();
        Vec3f B = (invA * Vec3f(deltaUV1.v, deltaUV2.v, 0.f)).normalize();

        Vec3f tbn = model->getNormal(uv);
        Vec3f normal = (T * tbn[0] + B * tbn[1] + N * tbn[2]).normalize();
        Vec3f l = -uniformMVP.multiply(lightDir, Vec4<float>::VECTOR).normalize();
        float intensity = std::max(l * normal, 0.f);

        color = model->getDiffuse(uv) * intensity;
        return false;
    }
};

/** Gets the shadow map, and returns the MVP on light (the viewport matrix is the same) */
Matrix getShadowMap(float *depthBuffer) {
    Matrix MV = lookAt(-lightDir, center, up);
    Matrix P = Matrix::identity(4);        // non perspective

    TGAImage shadowFrame(width, height, TGAImage::RGB);
//    DepthShader ds(P * MV);
    DepthShaderLite ds(P * MV);
    for (int i = 0; i < model->nFaces(); i += 1) {
        Vec4f screenCoords[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = ds.vertex(i, j);
        }
        triangle(screenCoords, ds, shadowFrame, depthBuffer);
    }
//    shadowFrame.flip_vertically();
//    shadowFrame.write_tga_file("../output/shadow_map.tga");
    return P * MV;
}

void precomputeOcclusionMap() {
    const int N = 1000;

    model = new Mesh("../obj/diablo3_pose/diablo3_pose.obj", "../obj/diablo3_pose/diablo3_pose_diffuse.tga");
    viewport = getViewport(width * 3 / 4, height * 3 / 4, depth, width / 8, height / 8);

    std::vector<Vec3f> lightPositions = getRandomPointsOnHemisphere(N);

    Vec2i textureSize = model->textureSize();
    TGAImage occlusionMaps[N];
    for (int i = 0; i < N; i += 1) {
        occlusionMaps[i] = TGAImage(textureSize[0], textureSize[1], TGAImage::RGB);
    }

#pragma omp parallel for shared(lightPositions, textureSize, center, up, model,\
projection, modelView, std::cout, N, occlusionMaps) default(none)
    for (int k = 0; k < N; k += 1) {
        const auto &l = lightPositions[k] * 2.f;

        // pass 1: SM
        auto depthBuffer = new float[height * width];
        std::fill(depthBuffer, depthBuffer + width * depth, std::numeric_limits<float>::lowest());
        Matrix MV = lookAt(l, center, up);
        Matrix P = Matrix::identity(4);

        TGAImage shadowFrame(width, height, TGAImage::RGB);
        DepthShaderLite ds(P * MV);
        for (int i = 0; i < model->nFaces(); i += 1) {
            Vec4f screenCoords[3];
            for (int j = 0; j < 3; j += 1) {
                screenCoords[j] = ds.vertex(i, j);
            }
            triangle(screenCoords, ds, shadowFrame, depthBuffer);
        }

        // pass 2: get occlusion map
        TGAImage frame(width, height, TGAImage::RGB);
        auto *zBuffer = new float[width * height];
        std::fill(zBuffer, zBuffer + width * height, std::numeric_limits<float>::lowest());
        modelView = lookAt(l, center, up);
        projection = Matrix::identity(4);

        OcclusionShader shader(projection * modelView, P *  MV, &occlusionMaps[k], depthBuffer);
        for (int i = 0; i < model->nFaces(); i += 1) {
            Vec4f screenCoords[3];
            for (int j = 0; j < 3; j += 1) {
                screenCoords[j] = shader.vertex(i, j);
            }
            triangle(screenCoords, shader, frame, zBuffer);
        }

        std::cout << "Done with sample point #" << k << ", total: " << N << std::endl;
        delete[] zBuffer;
        delete[] depthBuffer;
    }

    TGAImage occlusionMap(textureSize[0], textureSize[1], TGAImage::RGB);
#pragma omp parallel for shared(textureSize, occlusionMaps, occlusionMap) default(none)
//    for (int i = 0; i < textureSize[0]; i += 1) {
//        for (int j = 0; j < textureSize[1]; j += 1) {
//            int c = 0;
//            for (int k = 0; k < N; k += 1) {
//                c += occlusionMaps[k].get(i, j)[0];
//            }
//            c /= N;
//            occlusionMap.set(i, textureSize[1] - 1 - j, TGAColor(c, c, c, 255));
//        }
//    }
    for (int i=0; i<1024; i++) {
        for (int j=0; j<1024; j++) {
            for (int iter = 1; iter <= N; iter += 1) {
                float tmp = occlusionMap.get(i,j)[0];
                float c = (tmp*(iter-1)+occlusionMaps[iter].get(i,1023-j)[0])/(float)iter+.5f;
                occlusionMap.set(i, j, TGAColor(c,c,c,255));
            }
        }
    }

    occlusionMap.flip_vertically();
    occlusionMap.write_tga_file("../output/diablo3_pose_diffuse_ao.tga");
//    occlusionMaps[0].write_tga_file("../output/diablo3_pose_diffuse_aoooo.tga");
    delete model;
}

float getMaxSlope(const float *zBuffer, Vec2f p, Vec2f dir) {
    float t = 1.f;
    float z0 = zBuffer[int(p.x) + width * int(p.y)];

    Vec2f curr = p + dir * t;
    float maxSlope = 0.f;       // the tangent of the slope angle
    while (curr.x >= 0 && curr.y >= 0 && curr.x < width && curr.y < height) {
        float tanTheta = (zBuffer[int(curr.x) + width * int(curr.y)] - z0) / t;
        maxSlope = std::max(maxSlope, tanTheta);
        t += 1.f;
        curr = p + dir * t;
    }

    return maxSlope;
}

void demoSSAO() {
    eyePos = Vec3f(1.2, -.8, 3);
    model = new Mesh("../obj/diablo3_pose/diablo3_pose.obj");

    modelView = lookAt(eyePos, center, up);
    projection = getProjection((eyePos - center).norm());
    viewport = getViewport(width * 3 / 4, height * 3 / 4, depth, width / 8, height / 8);

    TGAImage frame(width, height, TGAImage::RGB);
    auto zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, std::numeric_limits<float>::lowest());

    DepthShaderLite ds(projection * modelView);
    for (int i = 0; i < model->nFaces(); i += 1) {
        Vec4f screenCoords[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = ds.vertex(i, j);
        }
        triangle(screenCoords, ds, frame, zBuffer);
    }

    for (int x = 0; x < width; x += 1) {
        for (int y = 0; y < height; y += 1) {
            if (zBuffer[x + y * width] < 0) {
                continue;
            }

            Vec2f curr((float) x + .5f, (float) y + .5f);
            float totSolidAngle = 0.f;

            for (int dx = -1; dx <= 1; dx += 1) {
                for (int dy = -1; dy <= 1; dy += 1) {
                    if (dx == 0 && dy == 0) {
                        continue;
                    }
                    Vec2f dir = Vec2f((float) dx, (float) dy).normalize();
                    float t = getMaxSlope(zBuffer, curr, dir);
                    totSolidAngle += M_PI / 4 * (1 - t / std::sqrt(1 + t * t));     // pi / 4 * (1 - cos(theta))
                }
            }
            float c = std::min(1.f, (float) std::pow(totSolidAngle / (2 * M_PI), 1.5f) * 1.5f);
            frame.set(x, y, WHITE * c);
        }
    }

    frame.flip_vertically();
    frame.write_tga_file("../output/output.tga");
    delete[] zBuffer;
    delete model;
}

void demo() {
//    model = new Mesh("../obj/african_head.obj", "../obj/african_head_diffuse.tga",
//                     "../obj/african_head_nm.tga", "../obj/african_head_spec.tga");
////                     "../obj/african_head_nm_tangent.tga", "../obj/african_head_spec.tga");
    model = new Mesh("../obj/diablo3_pose/diablo3_pose.obj", "../obj/diablo3_pose/diablo3_pose_diffuse.tga",
                     "../obj/diablo3_pose/diablo3_pose_nm.tga", "../obj/diablo3_pose/diablo3_pose_spec.tga");
    lightDir.normalize();

    modelView = lookAt(eyePos, center, up);
    projection = getProjection((eyePos - center).norm());
    viewport = getViewport(width * 3 / 4, height * 3 / 4, depth, width / 8, height / 8);

    shadowDepth = new float[width * height];
    std::fill(shadowDepth, shadowDepth + width * depth, std::numeric_limits<float>::lowest());
    Matrix lightMVP = getShadowMap(shadowDepth);

    TGAImage frame(width, height, TGAImage::RGB);
    auto *zBuffer = new float[width * height];
    std::fill(zBuffer, zBuffer + width * height, std::numeric_limits<float>::lowest());

    Shader shader(projection * modelView, lightMVP);
//    DiffuseShader shader(projection * modelView);
//    TGShader shader(projection * modelView);
    for (int i = 0; i < model->nFaces(); i += 1) {
        Vec4f screenCoords[3];
        for (int j = 0; j < 3; j += 1) {
            screenCoords[j] = shader.vertex(i, j);
        }
        triangle(screenCoords, shader, frame, zBuffer);
    }

    frame.flip_vertically();
    frame.write_tga_file("../output/output.tga");
    delete[] shadowDepth;
    delete[] zBuffer;
    delete model;
}

int main(int argc, char** argv) {
//    test();
//	demo();
    demoSSAO();
//	precomputeOcclusionMap();
	return 0;
}

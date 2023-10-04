#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<std::vector<int> > faces_pos;
	//
	std::vector<Vec2f> textures_;
	std::vector<Vec3f> normals_;
    std::vector<std::vector<int> > faces_tex;
    std::vector<std::vector<int> > faces_nor;
public:
	Model(const char *filename);
	~Model();
	int nverts();
	int nfaces();
	Vec3f vert(int i);
	std::vector<int> facePositions(int idx);
	//
	Vec2f uv(int i);
	Vec3f normal(int i);
    std::vector<int> faceTextures(int idx);      // texture indices of the idx-th face's vertices
    std::vector<int> faceNormals(int idx);   // normal indices of the idx-th face's vertices
};

//////////////////////////////////////////////////////////////////////////////////////////////////

struct Vertex {
    Vec3f position;
    Vec2f uv;
    Vec3f normal;
    Vertex(Vec3f p, Vec2f t, Vec3f n) : position(p), uv(t), normal(n) {}
};

class Mesh {
    std::vector<std::vector<size_t> > faces;
    std::vector<Vertex> vertices;
    TGAImage diffuseMap;
    TGAImage normalMap;
    TGAImage specularMap;
public:
    Mesh(const char *filename);
    Mesh(const char *filename, const char *textureFile);
    Mesh(const char *filename, const char *textureFile, const char *normalFile);
    Mesh(const char *filename, const char *textureFile, const char *normalFile, const char *specularFile);
    ~Mesh() {}

//    size_t nVertices() { return vertices.size(); }
    size_t nFaces() { return faces.size(); }
//    Vertex vert(size_t idx) const { return vertices[idx]; }
//    std::vector<size_t> face(int idx) { return faces[idx]; }
    Vertex getVertex(int iFace, int iVert) { return vertices[faces[iFace][iVert]]; }

    TGAColor getDiffuse(Vec2f uv) { return diffuseMap.get(uv); }
    Vec3f getNormal(Vec2f uv);
};

#endif //__MODEL_H__

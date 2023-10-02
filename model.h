#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"

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

#endif //__MODEL_H__

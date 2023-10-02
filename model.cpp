#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_pos(), textures_(), normals_() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 2, "vt")) {
            float zeroTrash;
            Vec2f v;
            iss >> trash >> trash >> v.x >> v.y >> zeroTrash;
            textures_.push_back(Vec2f(v.x, 1.f - v.y));
        } else if (!line.compare(0, 2, "vn")) {
            Vec3f v;
            iss >> trash >> trash >> v.x >> v.y >> v.z;
            normals_.push_back(v);
        }
//        else if (!line.compare(0, 2, "f ")) {
//            std::vector<int> f;
//            int itrash, idx;
//            iss >> trash;
//            while (iss >> idx >> trash >> itrash >> trash >> itrash) {
//                idx--; // in wavefront obj all indices start at 1, not zero
//                f.push_back(idx);
//            }
//            faces_pos.push_back(f);
//        }
        else if (!line.compare(0, 2, "f ")) {
            std::vector<int> fPositions;
            std::vector<int> fTextures;
            std::vector<int> fNormals;
            int idxPosition, idxTexture, idxNormal;
            iss >> trash;
            while (iss >> idxPosition >> trash >> idxTexture >> trash >> idxNormal) {
                fPositions.push_back(--idxPosition);
                fTextures.push_back(--idxTexture);
                fNormals.push_back(--idxNormal);
            }
            faces_pos.push_back(fPositions);
            faces_tex.push_back(fTextures);
            faces_nor.push_back(fNormals);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# " << faces_pos.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_pos.size();
}

std::vector<int> Model::facePositions(int idx) {
    return faces_pos[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

//
Vec2f Model::uv(int i) {
    return textures_[i];
}

Vec3f Model::normal(int i) {
    return normals_[i];
}

std::vector<int> Model::faceNormals(int i) {
    return faces_nor[i];
}

std::vector<int> Model::faceTextures(int i) {
    return faces_tex[i];
}




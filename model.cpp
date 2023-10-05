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
        } else if (!line.compare(0, 2, "f ")) {
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


//////////////////////////////////////////////////////////////////////////////////////////////////

Mesh::Mesh(const char *filename) {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) {
        std::cout << "failed to read " << filename << std::endl;
        std::exit(0);
    }

    std::vector<Vec3f> positions;
    std::vector<Vec2f> uvs;
    std::vector<Vec3f> normals;

    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());

        char trash;
        if (!line.compare(0, 2, "v ")) {
            Vec3f v;
            iss >> trash >> v.x >> v.y >> v.z;
            positions.push_back(v);

        } else if (!line.compare(0, 2, "vt")) {
            float zeroTrash;
            Vec2f v;
            iss >> trash >> trash >> v.x >> v.y >> zeroTrash;
            uvs.push_back(Vec2f(v.x, 1.f - v.y));

        } else if (!line.compare(0, 2, "vn")) {
            Vec3f v;
            iss >> trash >> trash >> v.x >> v.y >> v.z;
            normals.push_back(v);

        } else if (!line.compare(0, 2, "f ")) {
            int idxPosition, idxUV, idxNormal;
            std::vector<size_t> indices;
            iss >> trash;
            while (iss >> idxPosition >> trash >> idxUV >> trash >> idxNormal) {
                Vertex vtx(positions[--idxPosition], uvs[--idxUV], normals[--idxNormal]);
                indices.push_back(vertices.size());
                vertices.push_back(vtx);
            }
            faces.push_back(indices);
        }
    }
    std::cerr << "# v# " << vertices.size() << " f# " << faces.size() << std::endl;
}

Mesh::Mesh(const char *filename, const char *textureFile) : Mesh(filename) {
    if (!diffuseMap.read_tga_file(textureFile)) {
        std::cout << "Failed to read " << textureFile << std::endl;
        exit(0);
    }
}

Mesh::Mesh(const char *filename, const char *textureFile, const char *normalFile)
    : Mesh(filename, textureFile)
{
    if (!normalMap.read_tga_file(normalFile)) {
        std::cout << "Failed to read " << textureFile << std::endl;
        exit(0);
    }
}

Mesh::Mesh(const char *filename, const char *textureFile, const char *normalFile, const char *specularFile)
    : Mesh(filename, textureFile, normalFile)
{
    if (!specularMap.read_tga_file(specularFile)) {
        std::cout << "Failed to read " << textureFile << std::endl;
        exit(0);
    }
}

Vec3f Mesh::getNormal(Vec2f uv) {
    TGAColor c = normalMap.get(uv);
    Vec3f res;
    res.x = (float) c.r / 255.f * 2.f - 1.f;
    res.y = (float) c.g / 255.f * 2.f - 1.f;
    res.z = (float) c.b / 255.f * 2.f - 1.f;
    return res.normalize();
}

float Mesh::getSpecular(Vec2f uv) {
    return specularMap.get(uv)[0];
}

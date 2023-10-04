//
// Created by lxl on 2023/10/3.
//

#ifndef __GL_H__
#define __GL_H__

#include "tgaimage.h"
#include "geometry.h"

const TGAColor WHITE = TGAColor(255, 255, 255, 255);
const TGAColor RED   = TGAColor(255, 0, 0, 255);
const TGAColor GREEN = TGAColor(0, 255, 0, 255);
const TGAColor BLUE  = TGAColor(0, 0, 255, 255);

void lineNaive(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c);

void line1Helper(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c, bool isSteep);

void line1(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c);

void line2(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c);

void line(TGAImage &img, int x0, int y0, int x1, int y1, const TGAColor &c);

void line(TGAImage &img, const Vec2i &v0, const Vec2i &v1, const TGAColor &c);

void triangleBorder(TGAImage &img, const Vec2i &v0, const Vec2i &v1, const Vec2i &v2, const TGAColor &c);

void sortTriangleVertices(Vec2i &v0, Vec2i &v1, Vec2i &v2);

inline void horizontalLine(TGAImage &img, int y, int l, int r, const TGAColor &c);

void triangle(TGAImage &img, Vec2i v0, Vec2i v1, Vec2i v2, const TGAColor &c);

Vec3f barycentric(const Vec2i pts[], const Vec2i &p);

Vec3f barycentric(const Vec3f pts[], const Vec2i &p);

void triangle(TGAImage &img, const Vec2i pts[], const TGAColor &c);

/**
 * Draws a triangle considering the z-buffer.
 * The x and y coordinates of pts[i] is in the screen space ([0, width) x [0, height)),
 * and the z coordinate of pts[i] is in the world space. (i = 0, 1, 2)
 */
void triangle(TGAImage &img, float *zBuffer, const Vec3f pts[], const TGAColor &c);

/**
 * Converts the world coordinate to the screen.
 * Lets (x, y) converts from (-1, 1) x (-1, 1) to (0, width) x (0, height)
 * Remains z constant.
*/
Vec3f worldToScreen(const Vec3f &worldCoord, int width, int height);

void triangle(TGAImage &img, const TGAImage &texture, float *zBuffer, const Vec3f &light,
              const Vec3f positions[], const Vec2f uvs[], const Vec3f normals[]);

// for a point, not a vector
Matrix vec2mat(Vec3f v);

Vec3f mat2vec(Matrix m);

Vec3f vecTrans(Matrix m, Vec3f v);

//////////////////////////////////////////////////////////////////////////////////////////////////

extern Matrix modelView;
extern Matrix projection;
extern Matrix viewport;

Matrix lookAt(Vec3f eyePos, Vec3f center, Vec3f up);

/** zCamera is the distance between eye and center */
Matrix getProjection(float zCamera);

Matrix getViewport(int width, int height);

Matrix getViewport(int width, int height, int depth);

// [-1,1]*[-1,1]*[-1,1] is mapped onto the screen cube [x,x+w]*[y,y+h]*[0,d]
Matrix getViewport(int width, int height, int depth, int x, int y);

/** Shader interface */
class IShader {
protected:
    virtual ~IShader() {}

public:
    /** Transforms the coordinates of the vertices,
     * and prepares data (store in varying variables) for the fragment shader. */
    virtual Vec4f vertex(int iFace, int nthVert) = 0;

    /** Determines the color of the current pixel,
     * and we can discard current pixel by returning true */
    virtual bool fragment(Vec3f bary, TGAColor &color) = 0;
};

/**
 * Renders a triangle on the canvas.
 * @param pts the coordinates of the triangle's vertices (size of 3),
 *            usually calculated via shader.vertex()
 * @param shader the shader used on the triangle
 */
void triangle(Vec4f pts[], IShader &shader, TGAImage &canvas, float *zBuffer);

#endif //__GL_H__

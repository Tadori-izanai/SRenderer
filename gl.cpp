//
// Created by lxl on 2023/10/3.
//

#include "gl.h"

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

Vec3f worldToScreen(const Vec3f &worldCoord, int width, int height) {
    int x = (worldCoord.x + 1.f) * .5f * width;
    int y = (worldCoord.y + 1.f) * .5f * height;
    return Vec3f(x, y, worldCoord.z);
}

void triangle(TGAImage &img, const TGAImage &texture, float *zBuffer, const Vec3f &light,
              const Vec3f positions[], const Vec2f uvs[], const Vec3f normals[])
{
    const float eps = 0.001;
    int xLeft = std::numeric_limits<int>::max();
    int xRight = std::numeric_limits<int>::min();
    int yBottom = std::numeric_limits<int>::max();
    int yTop = std::numeric_limits<int>::min();
    for (int i = 0; i < 3; i += 1) {
        xLeft = std::min(xLeft, (int) positions[i].x);
        xRight = std::max(xRight, (int) positions[i].x);
        yBottom = std::min(yBottom, (int) positions[i].y);
        yTop = std::max(yTop, (int) positions[i].y);
    }
    xLeft = std::max(0, xLeft);
    xRight = std::min(img.get_width() - 1, xRight);
    yBottom = std::max(0, yBottom);
    yTop = std::min(img.get_height() - 1, yTop);

    for (int y = yBottom; y <= yTop; y += 1) {
        for (int x = xLeft; x <= xRight; x += 1) {
            Vec3f bc = barycentric(positions, Vec2i(x, y));
            if (bc.x < -eps || bc.y < -eps || bc.z < -eps) {
                continue;
            }

            float zWorldCoord = 0;
            Vec3f nn;
            Vec2f uv;
            for (int i = 0; i < 3; i += 1) {
                zWorldCoord += positions[i].z * bc[i];
                nn = nn + normals[i] * bc[i];
                uv = uv + uvs[i] * bc[i];
            }
            if (zWorldCoord > zBuffer[x + img.get_width() * y]) {
                zBuffer[x + img.get_width() * y] = zWorldCoord;

                nn.normalize();
                float intensity = -(light * nn);
                TGAColor c = texture.get(uv) * intensity;
                img.set(x, y, c);
            }
        }
    }
}

Matrix vec2mat(Vec3f v) {
    Matrix m(4, 4);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

Vec3f mat2vec(Matrix m) {
    return Vec3f(m[0][0], m[1][0], m[2][0]) / m[3][0];
}

Vec3f vecTrans(Matrix m, Vec3f v) {
    float elem[3] = {0};
    for (int i = 0; i < 3; i += 1) {
        for (int j = 0; j < 3; j += 1) {
            elem[i] += m[i][j] * v[j];
        }
    }
    return Vec3f(elem[0], elem[1], elem[2]);
}


//////////////////////////////////////////////////////////////////////////////////////////////////

Matrix modelView;
Matrix projection;
Matrix viewport;

Matrix lookAt(Vec3f eyePos, Vec3f center, Vec3f up) {
    Matrix translation = Matrix::identity(4);
    for (int i = 0; i < 3; i += 1) {
        translation[i][3] = -center[i];
    }

    Vec3f g = (center - eyePos).normalize();    // along -z'
    Vec3f r = (g ^ up).normalize();             // along x'
    Vec3f t = (r ^ g).normalize();              // along y'
    Matrix rotation = Matrix::identity(4);
    for (int j = 0; j < 3; j += 1) {
        rotation[0][j] = r[j];
        rotation[1][j] = t[j];
        rotation[2][j] = -g[j];
    }

    return rotation * translation;
}

Matrix getProjection(float zCamera) {
    Matrix m = Matrix::identity(4);
    m[3][2] = -1 / zCamera;
    return m;
}

Matrix getViewport(int width, int height) {
    Matrix m = Matrix::identity(4);
    m[0][0] = width / 2.f;
    m[0][3] = width / 2.f;
    m[1][1] = height / 2.f;
    m[1][3] = height / 2.f;
    return m;
}

Matrix getViewport(int width, int height, int depth) {
    Matrix m = Matrix::identity(4);
    m[0][0] = width / 2.f;
    m[0][3] = width / 2.f;
    m[1][1] = height / 2.f;
    m[1][3] = height / 2.f;
    m[2][2] = depth / 2.f;
    m[2][3] = depth / 2.f;
    return m;
}

Matrix getViewport(int width, int height, int depth, int x, int y) {
    Matrix m = Matrix::identity(4);
    m[0][0] = width / 2.f;
    m[0][3] = width / 2.f + x;
    m[1][1] = height / 2.f;
    m[1][3] = height / 2.f + y;
    m[2][2] = depth / 2.f;
    m[2][3] = depth / 2.f;
    return m;
}

Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
    float eps = 1e-2;
    Vec2f AB(B - A);
    Vec2f AC(C - A);
    Vec2f PA(A - P);
    Vec3f crossProd = Vec3f(AB.x, AC.x, PA.x) ^ Vec3f(AB.y, AC.y, PA.y);
    if (std::abs(crossProd.z) < eps) {
        return {-1, 1, 1};
    }
    float u = crossProd.x / crossProd.z;
    float v = crossProd.y / crossProd.z;
    return {1.0f - u - v, u, v};
}

/** Returns Vec4f(xLeft, xRight, yBottom, yTop) */
Vec4f getBBox(int width, int height, const Vec3f vertices[]) {
    float xLeft = std::numeric_limits<float>::max();
    float xRight = std::numeric_limits<float>::lowest();
    float yBottom = std::numeric_limits<float>::max();
    float yTop = std::numeric_limits<float>::lowest();
    for (int i = 0; i < 3; i += 1) {
        xLeft = std::min(xLeft, vertices[i].x);
        xRight = std::max(xRight, vertices[i].x);
        yBottom = std::min(yBottom, vertices[i].y);
        yTop = std::max(yTop, vertices[i].y);
    }
    xLeft = std::max(0.f, xLeft);
    xRight = std::min((float) width - 1, xRight);
    yBottom = std::max(0.f, yBottom);
    yTop = std::min((float) height - 1, yTop);

    return {xLeft, xRight, yBottom, yTop};
}

void triangle(Vec4f pts[], IShader &shader, TGAImage &canvas, float *zBuffer) {
    const float eps = 0.001;
    Vec3f vertCoords[3];
    for (int i = 0; i < 3; i += 1) {
        vertCoords[i] = pts[i].proj3();
    }

    Vec4f bbox = getBBox(canvas.get_width(), canvas.get_height(), vertCoords);
    for (int y = bbox[2]; y <= bbox[3]; y += 1) {
        for (int x = bbox[0]; x <= bbox[1]; x += 1) {
            Vec2f curr(x + .5f, y + .5f);
            Vec3f bc = barycentric(vertCoords[0].xy(), vertCoords[1].xy(), vertCoords[2].xy(), curr);
            if (bc.x < -eps || bc.y < -eps || bc.z < -eps) {
                continue;
            }
            float zCurr = bc[0] * vertCoords[0].z + bc[1] * vertCoords[1].z + bc[2] * vertCoords[2].z;
            TGAColor color;
            if (zCurr > zBuffer[x + canvas.get_width() * y] && !shader.fragment(bc, color)) {
                zBuffer[x + canvas.get_width() * y] = zCurr;
                canvas.set(x, y, color);
            }
        }
    }
}


#include <iostream>
#include <fstream>
#include <vector>
#include "struct.h"
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
using namespace rapidxml;

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

int ww, wh;
float alpha, beta, rad;

std::vector<std::string> renders;

std::vector<float> transf;
vertice pos, l, u, per;

std::vector<GLuint> buffers;
std::vector<int> bufferSizes;

struct CR_Curve {
    std::vector<vertice> points;
    float time;
    bool align = false;
};
std::vector<CR_Curve> curves;

void spherical2Cartesian() {
    pos.x = rad * cosf(beta) * sinf(alpha);
    pos.y = rad * sinf(beta);
    pos.z = rad * cosf(beta) * cosf(alpha);
}

void changeSize(int w, int h) {
    if (h == 0) h = 1;
    float ratio = w * 1.0f / h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(per.x, ratio, per.y, per.z);
    glMatrixMode(GL_MODELVIEW);
}

void multMatrixVector(float *m, float *v, float *res) {
	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}
}

float multVectorVector(float *a, float *b){
    float r;
    r = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    return r;
}

void buildRotMatrix(float *x, float *y, float *z, float *m) {
	m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
	m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
	m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
	m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
}

void cross(float *a, float *b, float *res) {
	res[0] = a[1]*b[2] - a[2]*b[1];
	res[1] = a[2]*b[0] - a[0]*b[2];
	res[2] = a[0]*b[1] - a[1]*b[0];
}

float length(float *v) {
	float res = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	return res;
}

void normalize(float *a) {
	float l = length(a);
	a[0] = a[0]/l;
	a[1] = a[1]/l;
	a[2] = a[2]/l;
}

void getCatmullRomPoint(float t, vertice p0, vertice p1, vertice p2, vertice p3, vertice *pos, float *deriv) {
	// catmull-rom matrix
	float m[4][4] = {{-0.5f,  1.5f, -1.5f,  0.5f},
					 { 1.0f, -2.5f,  2.0f, -0.5f},
					 {-0.5f,  0.0f,  0.5f,  0.0f},
					 { 0.0f,  1.0f,  0.0f,  0.0f}};
	
	float tt[4] = {t*t*t, t*t, t, 1};
	float td[4] = {3*t*t, 2*t, 1, 0};
	float aux[4];

    float px[4] = {p0.x, p1.x, p2.x, p3.x};
    multMatrixVector((float*)m, px, aux);
    deriv[0] = multVectorVector(td, aux);
    pos->x = multVectorVector(tt, aux);

    float py[4] = {p0.y, p1.y, p2.y, p3.y};
    multMatrixVector((float*)m, py, aux);
    deriv[1] = multVectorVector(td, aux);
    pos->y = multVectorVector(tt, aux);

    float pz[4] = {p0.z, p1.z, p2.z, p3.z};
    multMatrixVector((float*)m, pz, aux);
    deriv[2] = multVectorVector(td, aux);
    pos->z = multVectorVector(tt, aux);
}

void getGlobalCatmullRomPoint(float gt, int i_curve, vertice *pos, float *deriv) {
    int pointCount = curves[i_curve].points.size();
    float t = gt * pointCount;
    int index = (int)t;
    t = t - index;

    int indices[4];
    indices[0] = (index + pointCount - 1) % pointCount;
    indices[1] = (indices[0] + 1) % pointCount;
    indices[2] = (indices[1] + 1) % pointCount;
    indices[3] = (indices[2] + 1) % pointCount;

    getCatmullRomPoint(t,
        curves[i_curve].points[indices[0]],
        curves[i_curve].points[indices[1]],
        curves[i_curve].points[indices[2]],
        curves[i_curve].points[indices[3]],
        pos, deriv);
}

void renderCatmullRomCurve(int index) {
    // draw curve using line segments with GL_LINE_LOOP
        float deriv[3];
        vertice pos;
        float tess = 100.0f;
        glBegin(GL_LINE_LOOP);
        for (int i = 0; i < tess; i++) {
            getGlobalCatmullRomPoint(i / tess, index, &pos, deriv);
            glVertex3f(pos.x, pos.y, pos.z);
        }
        glEnd();
    }

void renderScene(void) {
    if (renders.empty()) {
        std::cerr << "Nenhum ficheiro 3D encontrado!" << std::endl;
        return;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt(pos.x, pos.y, pos.z,
              l.x, l.y, l.z,
              u.x, u.y, u.z);

    glBegin(GL_LINES);
    // X axis
    glColor3f(1.0f, 0.0f, 0.0f); //R
    glVertex3f(-100.0f, 0.0f, 0.0f);
    glVertex3f(100.0f, 0.0f, 0.0f);
    // Y Axis
    glColor3f(0.0f, 1.0f, 0.0f); //G
    glVertex3f(0.0f, -100.0f, 0.0f);
    glVertex3f(0.0f, 100.0f, 0.0f);
    // Z Axis
    glColor3f(0.0f, 0.0f, 1.0f); //B
    glVertex3f(0.0f, 0.0f, -100.0f);
    glVertex3f(0.0f, 0.0f, 100.0f);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(1.0f, 1.0f, 1.0f);

    int i_transf = 0;
    int i_vbo = 0;
    int i_curve = 0;
    
    for (const auto& render : renders) {
        if (render == "translate") {
            glTranslatef(transf[i_transf], transf[i_transf + 1], transf[i_transf + 2]);
            i_transf += 3;
        }
        else if (render == "rotate") {
            glRotatef(transf[i_transf], transf[i_transf + 1], transf[i_transf + 2], transf[i_transf + 3]);
            i_transf += 4;
        }
        else if (render == "scale") {
            glScalef(transf[i_transf], transf[i_transf + 1], transf[i_transf + 2]);
            i_transf += 3;
        }
        else if (render == "catmull") {
            float et = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
            float time = curves[i_curve].time;
            float gt = fmod(et / time, 1.0f); 
                        
            renderCatmullRomCurve(i_curve);
            
            static float vecY[3] = { 0, 1, 0 };
            float vecX[3], vecZ[3], mat[16];
            vertice pos;
            getGlobalCatmullRomPoint(gt, i_curve, &pos, vecX);
            glTranslatef(pos.x, pos.y, pos.z); //position
            
            if (curves[i_curve].align){
                cross(vecX, vecY, vecZ);
                cross(vecZ, vecX, vecY);
                normalize(vecX);
                normalize(vecY);
                normalize(vecZ);
                buildRotMatrix(vecX, vecY, vecZ, mat); 
                glMultMatrixf(mat); //rotation
            }
            i_curve++;
        }
        else if (render == "rotateTime") {
            float et = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
            float time = transf[i_transf];
            float ang = (et * 360.0f) / time;
            glRotatef(ang, transf[i_transf + 1], transf[i_transf + 2], transf[i_transf + 3]);
            i_transf += 4;
        }
        else if (render == "push") {
            glPushMatrix();
        }
        else if (render == "pop") {
            glPopMatrix();
        }
        else {
            glBindBuffer(GL_ARRAY_BUFFER, buffers[i_vbo]);
            glEnableClientState(GL_VERTEX_ARRAY);
            glVertexPointer(3, GL_FLOAT, 0, 0);
            glDrawArrays(GL_TRIANGLES, 0, bufferSizes[i_vbo]);
            glDisableClientState(GL_VERTEX_ARRAY);
            i_vbo++;
        }
    }
    
    glutSwapBuffers();
}

void processKeys(unsigned char c, int xx, int yy) {
    switch (c) {
    case 'a':
        alpha -= 0.1f;
        break;
    case 'd':
        alpha += 0.1f;
        break;
    case 'e':
        rad += 0.1f;
        break;
    case 'q':
        rad -= 0.1f;
        if (rad < 0.1f)
            rad = 0.1f;
        break;
    case 's':
        beta -= 0.1f;
        if (beta < -1.5f)
            beta = -1.5f;
        break;
    case 'w':
        beta += 0.1f;
        if (beta > 1.5f)
            beta = 1.5f;
        break;
    default:
        break;
    }
    spherical2Cartesian();
    glutPostRedisplay();
}

void processSpecialKeys(int key, int xx, int yy) {
    switch (key) {

    case GLUT_KEY_RIGHT:
        alpha -= 0.1; break;

    case GLUT_KEY_LEFT:
        alpha += 0.1; break;

    case GLUT_KEY_UP:
        beta += 0.1f;
        if (beta > 1.5f)
            beta = 1.5f;
        break;

    case GLUT_KEY_DOWN:
        beta -= 0.1f;
        if (beta < -1.5f)
            beta = -1.5f;
        break;

    case GLUT_KEY_PAGE_DOWN:
        rad -= 0.1f;
        if (rad < 0.1f)
            rad = 0.1f;
        break;

    case GLUT_KEY_PAGE_UP: 
        rad += 0.1f; 
        break;
    }
    spherical2Cartesian();
    glutPostRedisplay();
}

void processGroup(xml_node<>* group) {
    if (!group) return;

    renders.push_back("push");

    xml_node<>* transform = group->first_node("transform");
    if (transform) {

        for (xml_node<>* t = transform->first_node(); t; t = t->next_sibling()) {
            std::string name = t->name();

            if (name == "translate" && t->first_attribute("time")) {
                renders.push_back("catmull");
                CR_Curve curve;
                curve.time = std::stof(t->first_attribute("time")->value());
                curve.align = t->first_attribute("align") &&
                    std::string(t->first_attribute("align")->value()) == "true";

                for (xml_node<>* p = t->first_node("point"); p; p = p->next_sibling("point")) {
                    vertice v;
                    v.x = std::stof(p->first_attribute("x")->value());
                    v.y = std::stof(p->first_attribute("y")->value());
                    v.z = std::stof(p->first_attribute("z")->value());
                    curve.points.push_back(v);
                }

                curves.push_back(curve);
            }
            else if (name == "rotate" && t->first_attribute("time")) {
                renders.push_back("rotateTime");
                transf.push_back(std::stof(t->first_attribute("time")->value()));
                transf.push_back(std::stof(t->first_attribute("x")->value()));
                transf.push_back(std::stof(t->first_attribute("y")->value()));
                transf.push_back(std::stof(t->first_attribute("z")->value()));

            }
            else if (name == "translate" || name == "rotate" || name == "scale") {
                renders.push_back(name);
                if (name == "rotate" && t->first_attribute("angle")) {
                    transf.push_back(std::stof(t->first_attribute("angle")->value()));
                }
                if (t->first_attribute("x") && t->first_attribute("y") && t->first_attribute("z")) {
                    transf.push_back(std::stof(t->first_attribute("x")->value()));
                    transf.push_back(std::stof(t->first_attribute("y")->value()));
                    transf.push_back(std::stof(t->first_attribute("z")->value()));
                }
            }
        }
    }

    xml_node<>* models = group->first_node("models");
    if (models) {
        for (xml_node<>* model = models->first_node("model"); model; model = model->next_sibling("model")) {
            std::string filename = model->first_attribute("file")->value();
            renders.push_back(filename);

            std::ifstream fig(filename);
            std::vector<float> points;
            vertice aux;

            while (fig >> aux.x >> aux.y >> aux.z) {
                points.push_back(aux.x);
                points.push_back(aux.y);
                points.push_back(aux.z);
            }
            fig.close();

            GLuint vbo;
            glGenBuffers(1, &vbo);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(float), points.data(), GL_STATIC_DRAW);

            buffers.push_back(vbo);
            bufferSizes.push_back(points.size() / 3);
        }
    }

    for (xml_node<>* subgroup = group->first_node("group"); subgroup; subgroup = subgroup->next_sibling("group")) {
        processGroup(subgroup);
    }

    renders.push_back("pop");
}

int main(int argc, char** argv) {
    
    if (argc < 2) {
        std::cout << "Numero invalido de argumentos" << std::endl;
        return 1;
    }
    
    file<> xmlFile(argv[argc - 1]); //XML
    xml_document<> doc;
    doc.parse<0>(xmlFile.data());
    
    xml_node<>* win = doc.first_node("world")->first_node("window");
    if (win) {
        ww = std::stoi(win->first_attribute("width")->value());
        wh = std::stoi(win->first_attribute("height")->value());
    }
    
    xml_node<>* camera = doc.first_node("world")->first_node("camera");
    if (camera) {
        xml_node<>* position = camera->first_node("position");
        if (position) {
            pos.x = std::stof(position->first_attribute("x")->value());
            pos.y = std::stof(position->first_attribute("y")->value());
            pos.z = std::stof(position->first_attribute("z")->value());
            
            rad = sqrtf(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
            beta = acosf(pos.y / rad);
            alpha = atanf(pos.x / pos.z);
        }
        
        xml_node<>* look = camera->first_node("lookAt");
        if (look) {
            l.x = std::stof(look->first_attribute("x")->value());
            l.y = std::stof(look->first_attribute("y")->value());
            l.z = std::stof(look->first_attribute("z")->value());
        }
        
        xml_node<>* up = camera->first_node("up");
        if (up) {
            u.x = std::stof(up->first_attribute("x")->value());
            u.y = std::stof(up->first_attribute("y")->value());
            u.z = std::stof(up->first_attribute("z")->value());
        }
        
        xml_node<>* proj = camera->first_node("projection");
        if (proj) {
            per.x = std::stof(proj->first_attribute("fov")->value());
            per.y = std::stof(proj->first_attribute("near")->value());
            per.z = std::stof(proj->first_attribute("far")->value());
        }
    }
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(ww, wh);
    glutCreateWindow("CG@DI");

    glutReshapeFunc(changeSize);
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);

    glutKeyboardFunc(processKeys);
    glutSpecialFunc(processSpecialKeys);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    
    glewInit();

    xml_node<>* rootGroup = doc.first_node("world")->first_node("group");
    processGroup(rootGroup);

    doc.clear();    

    glutMainLoop();

    return 0;
}
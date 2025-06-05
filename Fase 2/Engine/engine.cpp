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
#include <GL/glut.h>
#endif

std::vector<std::string> renders;
std::vector<float> transf;
vertice pos, l, u, per;
int ww, wh;

float alpha, beta, rad;

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

    int i = 0;
    // Itera sobre todos os arquivos e desenha os modelos
    for (const auto& render : renders) {
        if (render == "translate") {
            glTranslatef(transf[i], transf[i + 1], transf[i + 2]);
            i += 3;
        }
        else if (render == "scale") {
            glScalef(transf[i], transf[i + 1], transf[i + 2]);
            i += 3;
        }
        else if (render == "rotate") {
            glRotatef(transf[i], transf[i + 1], transf[i + 2], transf[i + 3]);
            i += 4;
        }
        else if (render == "push") {
            glPushMatrix();
        }
        else if (render == "pop") {
            glPopMatrix();
        }
        else {
            std::ifstream fig(render);
            if (!fig) {
                std::cerr << "Erro ao abrir o ficheiro .3d: " << render << std::endl;
                continue;
            }

            glBegin(GL_TRIANGLES);
            vertice aux;
            while (fig >> aux.x >> aux.y >> aux.z) {
                glVertex3f(aux.x, aux.y, aux.z);
            }
            glEnd();

            fig.close();
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
            renders.push_back(name);
            if (name == "rotate") {
                float angle = std::stof(t->first_attribute("angle")->value());
                transf.push_back(angle);
            }
            transf.push_back(std::stof(t->first_attribute("x")->value()));
            transf.push_back(std::stof(t->first_attribute("y")->value()));
            transf.push_back(std::stof(t->first_attribute("z")->value()));
        }
    }

    xml_node<>* models = group->first_node("models");
    if (models) {
        for (xml_node<>* model = models->first_node("model"); model; model = model->next_sibling("model")) {
            renders.push_back(model->first_attribute("file")->value());
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
            alpha = asinf(pos.x / pos.y);
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

    xml_node<>* rootGroup = doc.first_node("world")->first_node("group");
    processGroup(rootGroup);

    doc.clear();

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(ww, wh);
    glutCreateWindow("CG@DI");

    glutReshapeFunc(changeSize);
    glutDisplayFunc(renderScene);

    glutKeyboardFunc(processKeys);
    glutSpecialFunc(processSpecialKeys);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glutMainLoop();

    //renders.clear();
    //transf.clear();

    return 0;
}
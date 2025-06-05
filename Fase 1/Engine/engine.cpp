#include <iostream>
#include <fstream>
#include <vector>
#include "struct.h"
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
using namespace rapidxml;

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

std::vector<std::string> files3d;
vertice pos, l, u, per;
int ww, wh;

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
    if (files3d.empty()) {
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

    // Itera sobre todos os arquivos e desenha os modelos
    for (const auto& file3d : files3d) {
        std::ifstream fig(file3d);
        if (!fig) {
            std::cerr << "Erro ao abrir o ficheiro .3d: " << file3d << std::endl;
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

    glutSwapBuffers();
}

int main(int argc, char** argv) {

    if (argc < 2) {
        std::cout << "Numero invalido de argumentos" << std::endl;
        return 1;
    }

    file<> xmlFile(argv[argc-1]); //XML
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

    xml_node<>* models = doc.first_node("world")->first_node("group")->first_node("models");
    if (models) {
        for (xml_node<>* model = models->first_node("model"); model; model = model->next_sibling("model")) {
            if(model){
                files3d.push_back(model->first_attribute("file")->value());
            }
        }
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(ww, wh);
    glutCreateWindow("CG@DI");

    glutReshapeFunc(changeSize);
    glutDisplayFunc(renderScene);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glutMainLoop();

    return 0;
}

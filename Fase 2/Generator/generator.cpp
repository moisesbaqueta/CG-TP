#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#define _USE_MATH_DEFINES
#include <math.h>

#include "struct.h" 


void drawPlane(float length, int divs, std::string file) {
    vertice v, v1, v2, v3;
    v.x = -(length / 2);
    v.y = 0.0f;
    v.z = -(length / 2);
    float inc = length / divs;

    std::ofstream plane(file);
    if (!plane) {
        std::cerr << "Erro ao abrir o arquivo: " << file << std::endl;
        return;
    }

    for (int i = 0; i < divs; i++) {
        for (int j = 0; j < divs; j++) {
            v1.x = v.x;
            v1.y = 0.0f;
            v1.z = v.z + inc;

            v2.x = v.x + inc;
            v2.y = 0.0f;
            v2.z = v.z;

            v3.x = v.x + inc;
            v3.y = 0.0f;
            v3.z = v.z + inc;

            plane << v.x << " " << v.y << " " << v.z << "\n";
            plane << v1.x << " " << v1.y << " " << v1.z << "\n";
            plane << v2.x << " " << v2.y << " " << v2.z << "\n";
            plane << v2.x << " " << v2.y << " " << v2.z << "\n";
            plane << v1.x << " " << v1.y << " " << v1.z << "\n";
            plane << v3.x << " " << v3.y << " " << v3.z << "\n";

            v.z += inc;
        }
        v.z = -(length / 2);
        v.x += inc;
    }
    plane.close();
}

void drawBox(float dimension, int divisions, std::string file) {
    float halfDim = dimension / 2.0f;
    float step = dimension / divisions;
    float x0, x1, y0, y1, z0, z1;

    vertice v0, v1, v2, v3;
    std::ofstream box(file);
    if (!box) {
        std::cerr << "Erro ao abrir o arquivo: " << file << std::endl;
        return;
    }

    for (int i = 0; i < divisions; i++) {
                // Faces Frontal e Traseira
        for (int j = 0; j < divisions; j++) {
            x0 = -halfDim + i * step;
            x1 = -halfDim + (i + 1) * step;
            y0 = -halfDim + j * step;
            y1 = -halfDim + (j + 1) * step;

            v0.x = x0;
            v0.y = y0;
            v0.z = halfDim;

            v1.x = x1;
            v1.y = y0;
            v1.z = halfDim;

            v2.x = x1;
            v2.y = y1;
            v2.z = halfDim;

            v3.x = x0;
            v3.y = y1;
            v3.z = halfDim;

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
            box << v3.x << " " << v3.y << " " << v3.z << "\n";

            v0.z = -halfDim;
            v1.z = -halfDim;
            v2.z = -halfDim;
            v3.z = -halfDim;

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v3.x << " " << v3.y << " " << v3.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
        }
        // Faces Superior e Inferior
        for (int j = 0; j < divisions; j++) {
            x0 = -halfDim + i * step;
            x1 = -halfDim + (i + 1) * step;
            z0 = -halfDim + j * step;
            z1 = -halfDim + (j + 1) * step;

            v0.x = x0;
            v0.y = halfDim;
            v0.z = z0;

            v1.x = x1;
            v1.y = halfDim;
            v1.z = z1;

            v2.x = x1;
            v2.y = halfDim;
            v2.z = z0;

            v3.x = x0;
            v3.y = halfDim;
            v3.z = z1;

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v3.x << " " << v3.y << " " << v3.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";

            v0.y = -halfDim;
            v1.y = -halfDim;
            v2.y = -halfDim;
            v3.y = -halfDim;

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";
            box << v3.x << " " << v3.y << " " << v3.z << "\n";
        }

        // Faces Direita e Esquerda
        for (int j = 0; j < divisions; j++) {
            y0 = -halfDim + i * step;
            y1 = -halfDim + (i + 1) * step;
            z0 = -halfDim + j * step;
            z1 = -halfDim + (j + 1) * step;

            v0.x = halfDim;
            v0.y = y0;
            v0.z = z0;

            v1.x = halfDim;
            v1.y = y1;
            v1.z = z0;

            v2.x = halfDim;
            v2.y = y1;
            v2.z = z1;

            v3.x = halfDim;
            v3.y = y0;
            v3.z = z1;

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
            box << v3.x << " " << v3.y << " " << v3.z << "\n";

            v0.x = -halfDim;
            v1.x = -halfDim;
            v2.x = -halfDim;
            v3.x = -halfDim;

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
            box << v1.x << " " << v1.y << " " << v1.z << "\n";

            box << v0.x << " " << v0.y << " " << v0.z << "\n";
            box << v3.x << " " << v3.y << " " << v3.z << "\n";
            box << v2.x << " " << v2.y << " " << v2.z << "\n";
        }

    }

    box.close();
}

void drawSphere(float radius, int slices, int stacks, std::string file) {
    float angx = 0.0f;
    float incx = 2 * M_PI / slices;
    float angy = M_PI / 2;
    float incy = M_PI / stacks;

    vertice v0, v1, v2, v3;
    std::ofstream sphere(file);
    if (!sphere) {
        std::cerr << "Erro ao abrir o arquivo: " << file << std::endl;
        return;
    }

    for (int i = 0; i < slices; i++) {
        for (int j = 1; j <= stacks; j++) {
            if (j == stacks) {
                v0.x = 0.0f;
                v0.y = radius * sinf(angy - incy);
                v0.z = 0.0f;

                v1.x = radius * cosf(angy) * sinf(angx + incx);
                v1.y = radius * sinf(angy);
                v1.z = radius * cosf(angy) * cosf(angx + incx);

                v2.x = radius * cosf(angy) * sinf(angx);
                v2.y = radius * sinf(angy);
                v2.z = radius * cosf(angy) * cosf(angx);

                sphere << v0.x << " " << v0.y << " " << v0.z << "\n";
                sphere << v1.x << " " << v1.y << " " << v1.z << "\n";
                sphere << v2.x << " " << v2.y << " " << v2.z << "\n";
            }
            else {
                v0.x = radius * cosf(angy) * sinf(angx);
                v0.y = radius * sinf(angy);
                v0.z = radius * cosf(angy) * cosf(angx);

                v1.x = radius * cosf(angy - incy) * sinf(angx);
                v1.y = radius * sinf(angy - incy);
                v1.z = radius * cosf(angy - incy) * cosf(angx);

                v2.x = radius * cosf(angy - incy) * sinf(angx + incx);
                v2.y = radius * sinf(angy - incy);
                v2.z = radius * cosf(angy - incy) * cosf(angx + incx);

                v3.x = radius * cosf(angy) * sinf(angx + incx);
                v3.y = radius * sinf(angy);
                v3.z = radius * cosf(angy) * cosf(angx + incx);

                sphere << v0.x << " " << v0.y << " " << v0.z << "\n";
                sphere << v1.x << " " << v1.y << " " << v1.z << "\n";
                sphere << v2.x << " " << v2.y << " " << v2.z << "\n";

                sphere << v2.x << " " << v2.y << " " << v2.z << "\n";
                sphere << v3.x << " " << v3.y << " " << v3.z << "\n";
                sphere << v0.x << " " << v0.y << " " << v0.z << "\n";
            }
            angy -= incy;
        }
        angy = M_PI / 2;
        angx += incx;
    }
    sphere.close();
}

void drawCone(float radius, float height, int slices, int stacks, std::string file) {
    float angx = 0.0f;
    float incx = 2 * M_PI / slices;
    float r = radius;
    float yy = 0.0f;
    float incy = height / stacks;

    vertice v0, v1, v2, v3;
    std::ofstream cone(file);
    if (!cone) {
        std::cerr << "Erro ao abrir o arquivo: " << file << std::endl;
        return;
    }

    for (int i = 0; i < slices; i++) {
        v0.x = 0.0f;
        v0.y = yy;
        v0.z = 0.0f;

        v1.x = r * sinf(angx + incx);
        v1.y = yy;
        v1.z = r* cosf(angx + incx);

        v2.x = r * sinf(angx);
        v2.y = yy;
        v2.z = r* cosf(angx);

        cone << v0.x << " " << v0.y << " " << v0.z << "\n";
        cone << v1.x << " " << v1.y << " " << v1.z << "\n";
        cone << v2.x << " " << v2.y << " " << v2.z << "\n";

        for (int j = 1; j <= stacks; j++) {
            float r1 = radius * (height - (yy + incy)) / height;
            if (j == stacks) {
                v0.x = 0.0f;
                v0.y = yy + incy;
                v0.z = 0.0f;

                v1.x = r * sinf(angx);
                v1.y = yy;
                v1.z = r* cosf(angx);

                v2.x = r * sinf(angx + incx);
                v2.y = yy;
                v2.z = r* cosf(angx + incx);

                cone << v0.x << " " << v0.y << " " << v0.z << "\n";
                cone << v1.x << " " << v1.y << " " << v1.z << "\n";
                cone << v2.x << " " << v2.y << " " << v2.z << "\n";
            }
            else {
                v0.x = r1 * sinf(angx);
                v0.y = yy + incy;
                v0.z = r1* cosf(angx);

                v1.x = r * sinf(angx);
                v1.y = yy;
                v1.z = r* cosf(angx);

                v2.x = r * sinf(angx + incx);
                v2.y = yy;
                v2.z = r* cosf(angx + incx);

                v3.x = r1 * sinf(angx + incx);
                v3.y = yy + incy;
                v3.z = r1* cosf(angx + incx);

                cone << v0.x << " " << v0.y << " " << v0.z << "\n";
                cone << v1.x << " " << v1.y << " " << v1.z << "\n";
                cone << v2.x << " " << v2.y << " " << v2.z << "\n";

                cone << v2.x << " " << v2.y << " " << v2.z << "\n";
                cone << v3.x << " " << v3.y << " " << v3.z << "\n";
                cone << v0.x << " " << v0.y << " " << v0.z << "\n";
            }
            yy += incy;
            r = r1;
        }
        yy = 0.0f;
        r = radius;
        angx += incx;
    }

    cone.close();
}

void drawTorus(float radiusI, float radiusE, int sides, int rings, std::string file){
    vertice v0, v1, v2, v3;
    std::ofstream torus(file);
    if (!torus) {
        std::cerr << "Erro ao abrir o arquivo: " << file << std::endl;
        return;
    }

    float angi = 0.0f, ange = 0.0f, inci, ince;
    inci = (2.0f * M_PI) / sides;
    ince = (2.0f * M_PI) / rings;

    for(int i = 0; i < rings; i++){
        for(int j = 0; j < sides; j++){
            v0.x = (radiusE + radiusI * cosf(angi)) * cosf(ange);
            v0.y = radiusI * sinf(angi);
            v0.z = (radiusE + radiusI * cosf(angi)) * sinf(ange);

            v1.x = (radiusE + radiusI * cosf(angi)) * cosf(ange+ince);
            v1.y = radiusI * sinf(angi);
            v1.z = (radiusE + radiusI * cosf(angi)) * sinf(ange+ince);  

            v2.x = (radiusE + radiusI * cosf(angi+inci)) * cosf(ange);
            v2.y = radiusI * sinf(angi+inci);
            v2.z = (radiusE + radiusI * cosf(angi+inci)) * sinf(ange); 

            v3.x = (radiusE + radiusI * cosf(angi+inci)) * cosf(ange+ince);
            v3.y = radiusI * sinf(angi+inci);
            v3.z = (radiusE + radiusI * cosf(angi+inci)) * sinf(ange+ince);

            torus << v0.x << " " << v0.y << " " << v0.z << "\n";
            torus << v1.x << " " << v1.y << " " << v1.z << "\n";
            torus << v2.x << " " << v2.y << " " << v2.z << "\n";

            torus << v2.x << " " << v2.y << " " << v2.z << "\n";
            torus << v1.x << " " << v1.y << " " << v1.z << "\n";
            torus << v3.x << " " << v3.y << " " << v3.z << "\n";

            angi += inci;
        }
        ange += ince;
    }
    torus.close();
}

void drawCylinder(float radius, float height, int slices, std::string file){
    std::ofstream cylinder(file);
    if (!cylinder) {
        std::cerr << "Erro ao abrir o arquivo: " << file << std::endl;
        return;
    }
    
    float ang = 0.0;
	float inc = 2 * M_PI / slices;

    vertice v0, v1, v2;
    v0.x = 0.0f;
	v0.y = height / 2;
	v0.z = 0.0f;

    v1.y = height / 2; v2.y = height / 2;
	for (int i = 0; i < slices; i++) {
        v1.x = radius * sinf(ang);
        v1.z = radius * cosf(ang);
        
        v2.x = radius * sinf(ang + inc);
        v2.z = radius * cosf(ang + inc);

		cylinder << v0.x << " " << v0.y << " " << v0.z << "\n";
        cylinder << v1.x << " " << v1.y << " " << v1.z << "\n";
        cylinder << v2.x << " " << v2.y << " " << v2.z << "\n";

		cylinder << v0.x << " " << -v0.y << " " << v0.z << "\n";
        cylinder << v2.x << " " << -v2.y << " " << v2.z << "\n";
        cylinder << v1.x << " " << -v1.y << " " << v1.z << "\n";

		cylinder << v1.x << " " << v1.y << " " << v1.z << "\n";
        cylinder << v1.x << " " << -v1.y << " " << v1.z << "\n";
		cylinder << v2.x << " " << v2.y << " " << v2.z << "\n";

		cylinder << v2.x << " " << v2.y << " " << v2.z << "\n";
		cylinder << v1.x << " " << -v1.y << " " << v1.z << "\n";
		cylinder << v2.x << " " << -v2.y << " " << v2.z << "\n";

		ang += inc;
	}

    cylinder.close();
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cout << "Numero invalido de argumentos" << std::endl;
        return 1;
    }

    if (strcmp(argv[1], "plane") == 0) {
        if (argc < 5) return 1;
        float a = std::stof(argv[2]);
        int b = std::stoi(argv[3]);
        drawPlane(a, b, argv[argc-1]);
    }
    else if (strcmp(argv[1], "box") == 0) {
        if (argc < 5) return 1;
        float a = std::stof(argv[2]);
        int b = std::stoi(argv[3]);
        drawBox(a, b, argv[argc-1]);
    }
    else if (strcmp(argv[1], "sphere") == 0) {
        if (argc < 6) return 1;
        float a = std::stof(argv[2]);
        int b = std::stoi(argv[3]);
        int c = std::stoi(argv[4]);
        drawSphere(a, b, c, argv[argc-1]);
    }
    else if (strcmp(argv[1], "cone") == 0) {
        if (argc < 7) return 1;
        float a = std::stof(argv[2]);
        float b = std::stof(argv[3]);
        int c = std::stoi(argv[4]);
        int d = std::stoi(argv[5]);
        drawCone(a, b, c, d, argv[argc-1]);
    }
    else if(strcmp(argv[1], "torus") == 0){
        if (argc < 7) return 1;
        float a = std::stof(argv[2]);
        float b = std::stof(argv[3]);
        int c = std::stoi(argv[4]);
        int d = std::stoi(argv[5]);
        drawTorus(a, b, c, d, argv[argc-1]);
    }
    else if (strcmp(argv[1], "cylinder") == 0){
        if (argc < 6) return 1;
        float a = std::stof(argv[2]);
        int b = std::stoi(argv[3]);
        int c = std::stoi(argv[4]);
        drawCylinder(a, b, c, argv[argc-1]);
    }
    else {
        std::cout << "Erro: Figura Invalida" << std::endl;
        return 1;
    }

    std::cout << "Ficheiro pronto." << std::endl;

    return 0;
}
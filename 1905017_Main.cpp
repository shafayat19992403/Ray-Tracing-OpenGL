#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include "1905017_Header_old.h"
#include "bitmap_image.hpp"

#include <iostream>
#include <memory>
#include <functional>
#include <unordered_map>
#include <string>

using namespace std;

#define pi (2 * acos(0.0))
#define HEIGHT 0
#define ANGLE 1
#define WIDTH 1

bool isGrid = false;
bool isAxis = false;

double angle;
int recursionLevel;
int numberOfObjects;
int numberOfPointLights;
int numberOfSpotLights;
double cameraDimension[2];
double imageDimension[2];
double window[2] = {500, 500};
double viewAngle = 80;

bitmap_image image;
vector<Object *> objects;
vector<PointLight *> pointLights;
vector<SpotLight *> spotLights;

ofstream out("output_amr.txt");

int imageCount = 1;

Point positionOfCamera, up_vec, right_vec, look_vec;
int numberOfSegments;

void drawAxes(int size)
{
    if (isAxis)
    {
        glBegin(GL_LINES);
        {
            glColor3f(1, 0, 0);
            glVertex3f(size, 0, 0);
            glVertex3f(-size, 0, 0);

            glColor3f(0, 1, 0);
            glVertex3f(0, -size, 0);
            glVertex3f(0, size, 0);

            glColor3f(0, 0, 1);
            glVertex3f(0, 0, size);
            glVertex3f(0, 0, -size);
        }
        glEnd();
    }
}

void rotate(Point &vec, Point &axis, double angle)
{
    Point tempOne = vec.multiply(cos(angle));
    Point tempTwo = (axis.cross(vec)).multiply(sin(angle));
    Point res = tempOne.add(tempTwo);
    vec = res;
}

void image_reset()
{
    int i = 0;
    while (i < imageDimension[WIDTH])
    {
        int j = 0;
        while (j < imageDimension[HEIGHT])
        {
            image.set_pixel(i, j, 0, 0, 0);
            j++;
        }
        i++;
    }
}

void capture()
{
    cout << "Started capturing\n";

    image_reset();

    double planeDistance = (window[HEIGHT] / 2) / tan((viewAngle / 2) * (pi / 180));

    Point lookAdjustment = look_vec.multiply(planeDistance);
    Point rightAdjustment = right_vec.multiply(window[WIDTH] / 2.0);
    Point upAdjustment = up_vec.multiply(window[HEIGHT] / 2.0);

    Point topLeft = positionOfCamera.add(lookAdjustment);
    topLeft = topLeft.subtract(rightAdjustment);
    topLeft = topLeft.add(upAdjustment);

    double du = window[WIDTH] / imageDimension[WIDTH] * 1.0;
    double dv = window[HEIGHT] / imageDimension[HEIGHT] * 1.0;

    rightAdjustment = right_vec.multiply(du / 2.0);
    upAdjustment = up_vec.multiply(dv / 2.0);
    topLeft = topLeft.add(rightAdjustment);
    topLeft = topLeft.subtract(upAdjustment);

    int nearIndex = -1;
    double t, tMin;

    int width = imageDimension[WIDTH];
    int height = imageDimension[HEIGHT];
    for (int i = 0; i < width; i++)
    {
        // cout << "i: " << i << endl;
        for (int j = 0; j < height; j++)
        {
            Point pixel = topLeft.add(right_vec.multiply(i * du)).subtract(up_vec.multiply(j * dv));

            Ray ray = Ray(positionOfCamera, pixel.subtract(positionOfCamera));
            Color color = Color(0, 0, 0);

            tMin = -1;

            nearIndex = -1;

            int kk = 0;
            while (kk < (int)objects.size())
            {
                t = objects[kk]->intersect(ray, color, 0);
                if (t > 0)
                {
                    bool shouldUpdate = (nearIndex == -1 || t < tMin);
                    tMin = shouldUpdate ? t : tMin;
                    nearIndex = shouldUpdate ? kk : nearIndex;
                }
                kk++;
            }

            if (nearIndex != -1)
            {
                // color = Color(0, 0, 0);

                // a
                // if (objects[nearIndex] == NULL)
                // {
                //     continue;
                // }

                double t = objects[nearIndex]->intersect(ray, color, 1);

                color.red = std::max(0.0, std::min(color.red, 1.0));
                color.green = std::max(0.0, std::min(color.green, 1.0));
                color.blue = std::max(0.0, std::min(color.blue, 1.0));

                out << color.red << " " << color.green << " " << color.blue << endl;
                image.set_pixel(i, j, color.red * 255, color.green * 255, color.blue * 255);
            }
        }
    }

    cout << "outside loop" << endl;

    image.save_image("output" + to_string(imageCount) + ".bmp");
    imageCount++;
    cout << "Finished capturing\n";
}
void rotateLeftRight(int pos)
{
    rotate(right_vec, up_vec, pos * pi / 180);
    rotate(look_vec, up_vec, pos * pi / 180);
}

void rotateUpDown(int pos)
{
    rotate(up_vec, right_vec, pos * pi / 180);
    rotate(look_vec, right_vec, pos * pi / 180);
}

void tilt(int pos)
{
    rotate(right_vec, look_vec, pos * pi / 180);
    rotate(up_vec, look_vec, pos * pi / 180);
}

void keyboardListener(unsigned char key, int x, int y)
{
    if (key == '0')
    {
        capture();
    }
    else if (key == '1')
    {

        rotateLeftRight(1);
    }
    else if (key == '2')
    {

        rotateLeftRight(-1);
    }
    else if (key == '3')
    {

        rotateUpDown(1);
    }
    else if (key == '4')
    {

        rotateUpDown(-1);
    }
    else if (key == '5')
    {

        tilt(1);
    }
    else if (key == '6')
    {

        tilt(-1);
    }
}

void specialKeyListener(int key, int x, int y)
{
    if (key == GLUT_KEY_DOWN)
    {
        positionOfCamera = positionOfCamera.subtract(look_vec.multiply(3));
    }
    else if (key == GLUT_KEY_UP)
    {
        positionOfCamera = positionOfCamera.add(look_vec.multiply(3));
    }
    else if (key == GLUT_KEY_RIGHT)
    {
        positionOfCamera = positionOfCamera.add(right_vec.multiply(3));
    }
    else if (key == GLUT_KEY_LEFT)
    {
        positionOfCamera = positionOfCamera.subtract(right_vec.multiply(3));
    }
    else if (key == GLUT_KEY_PAGE_UP)
    {
        positionOfCamera = positionOfCamera.add(up_vec.multiply(3));
    }
    else if (key == GLUT_KEY_PAGE_DOWN)
    {
        positionOfCamera = positionOfCamera.subtract(up_vec.multiply(3));
    }
}

void drawObjects()
{
    // cout << "inside drawObjects\n";

    int i = 0;
    while (i < objects.size())
    {
        // cout << "Inside object loop\n";
        Object *obj = objects[i];
        obj->draw();
        i++;
        // cout << "draw done of object " << i << "\n";
    }
}

void drawPointLights()
{
    int i = 0;
    while (i < pointLights.size())
    {
        PointLight *pointLightObj = pointLights[i];
        pointLightObj->draw();
        i++;
    }
}

void drawSpotLights()
{
    int i = 0;
    while (i < spotLights.size())
    {
        SpotLight *spotLightObj = spotLights[i];
        spotLightObj->draw();
        i++;
    }
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(positionOfCamera.values.x, positionOfCamera.values.y, positionOfCamera.values.z,
              positionOfCamera.values.x + look_vec.values.x, positionOfCamera.values.y + look_vec.values.y, positionOfCamera.values.z + look_vec.values.z,
              up_vec.values.x, up_vec.values.y, up_vec.values.z);
    glMatrixMode(GL_MODELVIEW);
    // glLoadIdentity();
    // glMatrixMode(GL_MODELVIEW);
    // glLoadIdentity();
    drawAxes(100);
    // drawGrid();
    drawObjects();
    drawPointLights();
    drawSpotLights();
    // glFlush();
    glutSwapBuffers();
}
void createObject(ifstream &file, string objectType)
{
    if (objectType == "sphere")
    {
        Sphere *sphere = new Sphere();
        file >> *(sphere);
        objects.push_back(sphere);
    }
    else if (objectType == "triangle")
    {
        Triangle *triangle = new Triangle();
        file >> *(triangle);
        objects.push_back(triangle);
    }
    else if (objectType == "general")
    {
        GeneralObject *general = new GeneralObject();
        file >> *(general);
        objects.push_back(general);
    }
}

void animate()
{
    glutPostRedisplay();
}

void loadData()
{

    ifstream file("scene.txt");
    file >> recursionLevel >> imageDimension[HEIGHT];
    cout << "re:" << recursionLevel << endl;

    imageDimension[WIDTH] = imageDimension[HEIGHT];

    file >> numberOfObjects;
    string objectType;
    int i = 0;
    while (i < numberOfObjects)
    {

        file >> objectType;

        createObject(file, objectType);
        i++;
    }

    cout << "finished loading objects\n";

    file >> numberOfPointLights;
    i = 0;
    while (i < numberOfPointLights)
    {
        PointLight *pointLight = new PointLight();
        file >> *(pointLight);
        pointLights.push_back(pointLight);
        i++;
    }

    cout << "finished loading point lights\n";

    file >> numberOfSpotLights;
    i = 0;
    while (i < numberOfSpotLights)
    {
        SpotLight *spotLight = new SpotLight();
        file >> *(spotLight);
        spotLights.push_back(spotLight);
        i++;
    }

    cout << "finished loading spot lights\n";

    // Floor *floor = new Floor(400, 10);
    Floor *floor = new Floor(1000, 20);
    floor->color = Color(0.5, 0.5, 0.5);

    vector<double> coEfficients;
    coEfficients.push_back(0.4);
    coEfficients.push_back(0.2);
    coEfficients.push_back(0.2);
    coEfficients.push_back(0.2);

    floor->coeff = coEfficients;
    objects.push_back(floor);

    cout << "done loading data\n";
    cout << pointLights.size() << " " << spotLights.size() << endl;

    for (int i = 0; i < spotLights.size(); i++)
    {
        spotLights[i]->print();
    }
}

// void printData()
// {
//     cout << recursionLevel << endl;
//     cout << imageDimension[HEIGHT] << endl;

//     cout << numberOfObjects << endl;
//     for (int i = 0; i < numberOfObjects; i++)
//     {
//         // cout << objects[i]->type << endl;
//         objects[i]->print();
//     }

//     for (int i = 0; i < numberOfPointLights; i++)
//     {
//         pointLights[i]->print();
//     }

//     for (int i = 0; i < numberOfSpotLights; i++)
//     {
//         spotLights[i]->print();
//     }
// }

void init()
{
    positionOfCamera = Point(200, 0, 10);
    // positionOfCamera = Point(100, 100, -200);
    up_vec = Point(0, 0, 1);

    double neg = -1 / sqrt(2);
    double pos = 1 / sqrt(2);

    right_vec = Point(neg, pos, 0);
    look_vec = Point(neg, neg, 0);

    isGrid = true;
    isAxis = true;

    cameraDimension[HEIGHT] = 150;
    cameraDimension[ANGLE] = 1.0;
    angle = 0;
    numberOfSegments = 36;

    loadData();
    // printData();

    image = bitmap_image(imageDimension[WIDTH], imageDimension[HEIGHT]);

    glClearColor(0, 0, 0, 0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(80, 1, 1, 1000.0);
}

void clearData()
{
    int i = 0;
    while (i < objects.size())
    {
        delete objects[i];
        i++;
    }
    i = 0;
    while (i < pointLights.size())
    {
        delete pointLights[i];
        i++;
    }
    i = 0;
    while (i < spotLights.size())
    {
        delete spotLights[i];
        i++;
    }
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("Ray Tracing");
    init();
    glEnable(GL_DEPTH_TEST);
    glutDisplayFunc(display);
    glutIdleFunc(animate);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();
    clearData();
    return 0;
}
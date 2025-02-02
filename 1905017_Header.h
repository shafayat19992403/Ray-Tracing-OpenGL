#include "bitmap_image.hpp"
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <iostream>
#include <map>
#include <vector>
#include <iomanip>

using namespace std;

#define pi (2 * acos(0.0))
#define HEIGHT 0
#define WIDTH 1
#define LENGTH 2

extern bitmap_image image;
class Object;
extern vector<Object *> objects;
struct Values
{
public:
    double x, y, z, n;
    Values(double x, double y, double z, double n)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->n = n;
    }
    Values()
    {
        x = 0;
        y = 0;
        z = 0;
        n = 1.0;
    }
};

class Point
{
public:
    Values values;
    Point()
    {
        values = Values();
    }

    Point(double x, double y, double z)
    {
        values = Values(x, y, z, 1.0);
    }

    Point(double x, double y, double z, double n)
    {
        values = Values(x, y, z, n);
    }

    Point(const Point &p)
    {
        values.x = p.values.x;
        values.y = p.values.y;
        values.z = p.values.z;
        values.n = p.values.n;
    }

    Point add(Point p)
    {
        Point temp;
        temp.values.x = values.x + p.values.x;
        temp.values.y = values.y + p.values.y;
        temp.values.z = values.z + p.values.z;
        return temp;
    }

    Point subtract(Point p)
    {
        Point temp;
        temp.values.x = values.x - p.values.x;
        temp.values.y = values.y - p.values.y;
        temp.values.z = values.z - p.values.z;
        return temp;
    }

    Point multiply(double n)
    {
        Point temp;
        temp.values.x = values.x * n;
        temp.values.y = values.y * n;
        temp.values.z = values.z * n;
        return temp;
    }

    Point divide(double n)
    {
        Point temp;
        temp.values.x = values.x / n;
        temp.values.y = values.y / n;
        temp.values.z = values.z / n;
        return temp;
    }

    Point cross(Point p)
    {
        Point temp;
        temp.values.x = values.y * p.values.z - values.z * p.values.y;
        temp.values.y = values.z * p.values.x - values.x * p.values.z;
        temp.values.z = values.x * p.values.y - values.y * p.values.x;
        return temp;
    }

    double dot(Point p)
    {
        return values.x * p.values.x + values.y * p.values.y + values.z * p.values.z;
    }

    double magnitude()
    {
        return sqrt(values.x * values.x + values.y * values.y + values.z * values.z);
    }

    Point normalize()
    {
        Point temp;
        double mag = magnitude();
        temp.values.x = values.x / mag;
        temp.values.y = values.y / mag;
        temp.values.z = values.z / mag;
        return temp;
    }

    friend istream &operator>>(istream &is, Point &p)
    {
        double x, y, z;
        is >> x >> y >> z;

        p.values.x = x;
        p.values.y = y;
        p.values.z = z;

        return is;
    }
};

class Color
{
public:
    double red, green, blue;
    Color()
    {
        red = 0;
        green = 0;
        blue = 0;
    }

    Color(double red, double green, double blue)
    {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }

    Color(const Color &c)
    {
        red = c.red;
        green = c.green;
        blue = c.blue;
    }

    Color multiply(double n)
    {
        Color temp;
        temp.red = red * n;
        temp.green = green * n;
        temp.blue = blue * n;
        return temp;
    }

    Color add(Color c)
    {
        Color temp;
        temp.red = red + c.red;
        temp.green = green + c.green;
        temp.blue = blue + c.blue;
        return temp;
    }

    Color subtract(Color c)
    {
        Color temp;
        temp.red = red - c.red;
        temp.green = green - c.green;
        temp.blue = blue - c.blue;
        return temp;
    }

    Color divide(double n)
    {
        Color temp;
        temp.red = red / n;
        temp.green = green / n;
        temp.blue = blue / n;
        return temp;
    }

    Color mulColor(Color c)
    {
        Color temp;
        temp.red = red * c.red;
        temp.green = green * c.green;
        temp.blue = blue * c.blue;
        return temp;
    }

    friend istream &operator>>(istream &is, Color &c)
    {
        double r, g, b;
        is >> r >> g >> b;
        c.red = r;
        c.green = g;
        c.blue = b;
        return is;
    }
};

class PointLight
{
public:
    Point position;
    Color color;

    void draw()
    {
        glPointSize(5);
        glBegin(GL_POINTS);
        double r, g, b, x, y, z;
        r = color.red;
        g = color.green;
        b = color.blue;
        x = position.values.x;
        y = position.values.y;
        z = position.values.z;
        glColor3f(r, g, b);
        glVertex3f(x, y, z);
        glEnd();
    }

    friend istream &operator>>(istream &is, PointLight &pl)
    {
        is >> pl.position;
        is >> pl.color;
        return is;
    }

    void print()
    {
        cout << "Position: " << position.values.x << " " << position.values.y << " " << position.values.z << endl;
        cout << "Color: " << color.red << " " << color.green << " " << color.blue << endl;
    }
};

class SpotLight
{
public:
    PointLight pointLight;
    double cutoff;
    Point direction;

    void draw()
    {

        glPointSize(15);
        glBegin(GL_POINTS);
        glColor3f(pointLight.color.red, pointLight.color.green, pointLight.color.blue);
        double x, y, z;
        x = pointLight.position.values.x;
        y = pointLight.position.values.y;
        z = pointLight.position.values.z;
        glVertex3f(x, y, z);
        glEnd();
    }

    friend istream &operator>>(istream &is, SpotLight &sl)
    {
        is >> sl.pointLight >> sl.direction >> sl.cutoff;
        return is;
    }

    void print()
    {
        cout << "Position: " << pointLight.position.values.x << " " << pointLight.position.values.y << " " << pointLight.position.values.z << endl;
        cout << "Color: " << pointLight.color.red << " " << pointLight.color.green << " " << pointLight.color.blue << endl;
        cout << "Direction: " << direction.values.x << " " << direction.values.y << " " << direction.values.z << endl;
        cout << "Cutoff: " << cutoff << endl;
    }
};

class Ray
{
public:
    Point origin;
    Point direction;

    Ray()
    {
        origin = Point();
        direction = Point();
    }

    Ray(Point origin, Point direction)
    {
        this->origin = origin;
        direction = direction.normalize();
        this->direction = direction;
    }

    Ray(const Ray &r)
    {
        origin = r.origin;
        direction = r.direction;
    }

    friend istream &operator>>(istream &is, Ray &r)
    {
        is >> r.origin >> r.direction;
        return is;
    }
};

extern vector<PointLight *> pointLights;
extern vector<SpotLight *> spotLights;
extern int recursionLevel;

class Object
{
public:
    Point ref_point;
    Color color;
    vector<double> dimentions;
    vector<double> coeff;
    int shine;

    Object()
    {
        ref_point = Point(0, 0, 0);
        color = Color(0, 0, 0);
        shine = 0;
        for (int i = 0; i < 4; i++)
        {
            coeff.push_back(0);
        }

        for (int i = 0; i < 3; i++)
        {
            dimentions.push_back(0);
        }
    }

    virtual Ray getNormal(Point point, Ray i) = 0;

    virtual Color getColorAt(Point point)
    {
        // return Color(this->color.red, this->color.green, this->color.blue);
        Color colorTemp;
        colorTemp.red = color.red;
        colorTemp.green = color.green;
        colorTemp.blue = color.blue;
        return colorTemp;
    }

    virtual void draw() = 0;

    virtual double findIntersectionDistance(Ray ray) = 0;

    virtual double intersect(Ray ray, Color &color, int level)
    {
        // out<<ray<<endl<<color<<endl<<level<<endl;
        double t = findIntersectionDistance(ray);

        bool isDistanceNeg = (t < 0);
        bool isNoLevel = (level == 0);
        if (isDistanceNeg)
            return -1;

        if (isNoLevel)
            return t;

        // origin + t*direction

        Point intersectPoint = ray.origin.add(ray.direction.multiply(t));

        Color colorOfIntersectPoint = getColorAt(intersectPoint);

        // out << "Color of intersect point: " << colorOfIntersectPoint << endl;
        //   apply ambiance color
        // color = colorOfIntersectPoint.multiply(coeff[0]);
        // cout << "Color of intersect point: " << color.red << " " << color.green << " " << color.blue << endl;
        Color tempColor = colorOfIntersectPoint.multiply(coeff[0]);
        color.red = tempColor.red;
        color.green = tempColor.green;
        color.blue = tempColor.blue;
        // cout << "Color after applying co0: " << color.red << " " << color.green << " " << color.blue << endl;

        int *ptr = new int;

        for (int i = 0; i < pointLights.size(); i++)
        {

            Point lightPosition = pointLights[i]->position;
            Point lightDirection = intersectPoint.subtract(lightPosition).normalize();

            Ray lightRay = Ray(lightPosition, lightDirection);
            Ray normal = getNormal(intersectPoint, lightRay);

            double dot = lightRay.direction.dot(normal.direction);
            Point scaledNormal = normal.direction.multiply(2 * dot);
            Point reflectionDirection = lightRay.direction.subtract(scaledNormal);
            Ray reflection = Ray(intersectPoint, reflectionDirection);

            Point tempForT2 = intersectPoint.subtract(lightPosition);
            double t2 = tempForT2.magnitude();

            bool isT2Enough = (t2 < 1e-5);
            if (isT2Enough)
                continue;

            bool obscured = false;

            for (Object *obj : objects)
            {
                double t3 = obj->findIntersectionDistance(lightRay);

                bool isT3Positive = (t3 > 0);
                bool isT3LessThanT2 = (t3 + 1e-5 < t2);

                if (isT3Positive)
                {
                    if (isT3LessThanT2)
                    {
                        obscured = true;
                        break;
                    }
                }
            }

            if (!obscured)
            {
                // cout << "HEre in pointLight loop" << i << endl;
                double val = max(0.0, -lightRay.direction.dot(normal.direction));
                // cout << "Val: " << val << endl;
                //  Ray reflection = Ray(intersectPoint, lightRay.direction.subtract(normal.direction.multiply(2 * lightRay.direction.dot(normal.direction))));
                double dotProduct = lightRay.direction.dot(normal.direction);
                Point scaledNormal = normal.direction.multiply(2 * dotProduct);
                Point reflectionDirection = lightRay.direction.subtract(scaledNormal);
                Ray reflection = Ray(intersectPoint, reflectionDirection);

                double phong = max(0.0, -1 * ray.direction.dot(reflection.direction));
                phong = pow(phong, shine);
                // cout << "Phong: " << phong << endl;

                Color temp = pointLights[i]->color.multiply(coeff[1] * val);
                temp = temp.mulColor(colorOfIntersectPoint);
                // color = color.add(temp);
                color.red += temp.red;
                color.green += temp.green;
                color.blue += temp.blue;

                temp = pointLights[i]->color.multiply(coeff[2] * phong);
                temp = temp.mulColor(colorOfIntersectPoint);
                // color = color.add(temp);
                color.red += temp.red;
                color.green += temp.green;
                color.blue += temp.blue;

                // cout << "Color after pointLight: " << color.red << " " << color.green << " " << color.blue << endl;
            }
        }

        for (int i = 0; i < spotLights.size(); i++)
        {

            if (spotLights[i] == NULL && ptr == NULL)
                continue;

            Point lightPosition = spotLights[i]->pointLight.position;
            Point lightDirection = intersectPoint.subtract(lightPosition).normalize();
            // lightDirection = lightDirection.normalize();

            double dot = lightDirection.dot(spotLights[i]->direction);
            double denominator = lightDirection.magnitude() * spotLights[i]->direction.magnitude();
            // double angle = acos(dot / (lightDirection.magnitude() * spotLights[i]->direction.magnitude())) * (180 / pi);
            double angle = acos(dot / denominator) * (180 / pi);

            bool isAngleEnough = (fabs(angle) < spotLights[i]->cutoff);
            if (isAngleEnough)
            {
                // cout << "Color before spotLight: " << color.red << " " << color.green << " " << color.blue << endl;
                Ray lightRay = Ray(lightPosition, lightDirection);
                Ray normal = getNormal(intersectPoint, lightRay);

                double dotProduct = lightRay.direction.dot(normal.direction);
                Point scaledNormal = normal.direction.multiply(2 * dotProduct);
                Point reflectionDirection = lightRay.direction.subtract(scaledNormal);
                Ray reflection = Ray(intersectPoint, reflectionDirection);

                Point tempForT2 = intersectPoint.subtract(lightPosition);
                double t2 = tempForT2.magnitude();

                bool isT2Enough = (t2 > 1e-5);
                if (!isT2Enough)
                    continue;

                bool obscured = false;

                for (Object *obj : objects)
                {
                    double t3 = obj->findIntersectionDistance(lightRay);

                    bool isT3Positive = (t3 > 0);
                    bool isT3LessThanT2 = (t3 + 1e-5 < t2);

                    if (isT3Positive)
                    {
                        if (isT3LessThanT2)
                        {
                            obscured = true;
                            break;
                        }
                    }
                }

                if (!obscured)
                {

                    double val = max(0.0, -lightRay.direction.dot(normal.direction));
                    // cout << "Val: " << val << endl;
                    //  Ray reflection = Ray(intersectPoint, lightRay.direction.subtract(normal.direction.multiply(2 * lightRay.direction.dot(normal.direction))));
                    double dotProduct = lightRay.direction.dot(normal.direction);
                    Point scaledNormal = normal.direction.multiply(2 * dotProduct);
                    Point reflectionDirection = lightRay.direction.subtract(scaledNormal);
                    Ray reflection = Ray(intersectPoint, reflectionDirection);

                    double phong = max(0.0, -1 * ray.direction.dot(reflection.direction));
                    phong = pow(phong, shine);
                    // cout << "Phong: " << phong << endl;

                    Color temp = spotLights[i]->pointLight.color.multiply(coeff[1] * val);
                    temp = temp.mulColor(colorOfIntersectPoint);
                    // color = color.add(temp);
                    color.red += temp.red;
                    color.green += temp.green;
                    color.blue += temp.blue;

                    temp = spotLights[i]->pointLight.color.multiply(coeff[2] * phong);
                    temp = temp.mulColor(colorOfIntersectPoint);
                    // color = color.add(temp);
                    color.red += temp.red;
                    color.green += temp.green;
                    color.blue += temp.blue;

                    // cout << "Color after spotLight: " << color.red << " " << color.green << " " << color.blue << endl;
                }
            }
        }

        if (level < recursionLevel)
        {
            Ray normal = getNormal(intersectPoint, ray);

            // getReflectedRay
            double dotProduct = ray.direction.dot(normal.direction);
            Point vectorToSubtract = normal.direction.multiply(2 * dotProduct);
            Point reflectedDirection = ray.direction.subtract(vectorToSubtract);
            Ray reflectionRay = Ray(intersectPoint, reflectedDirection);

            reflectionRay.origin = reflectionRay.origin.add(reflectionRay.direction.multiply(1e-5));

            int nearestObjIndex = -1;
            double t = -1;
            double minT = 1e9;

            int j = 0;
            while (j < (int)objects.size())
            {
                t = objects[j]->intersect(reflectionRay, color, 0);

                bool isTinRange = (t > 0 && t < minT);
                if (isTinRange)
                {
                    minT = t;
                    nearestObjIndex = j;
                }
                j++;
            }

            if (nearestObjIndex != -1)
            {
                Color temp(0, 0, 0);
                double t = objects[nearestObjIndex]->intersect(reflectionRay, temp, level + 1);

                // color = color.add(temp.multiply(coeff[3]));
                Color temp2 = temp.multiply(coeff[3]);
                color.red += temp2.red;
                color.green += temp2.green;
                color.blue += temp2.blue;
            }
        }
        return t;
    }

    virtual ~Object()
    {
        dimentions.clear();
        coeff.clear();
        dimentions.shrink_to_fit();
        coeff.shrink_to_fit();
    };

    void
    processPointLight(PointLight &pointLight, Point intersectPoint, Ray ray, Color &color, Color colorOfIntersectPoint)
    {
    }

    void processSpotLight(SpotLight &spotLight, Point intersectPoint, Ray ray, Color &color, Color colorOfIntersectPoint)
    {
    }
};

double det(double arr[3][3])
{
    double ans = 0;
    for (int i = 0; i < 3; i++)
    {
        ans += arr[0][i] * (arr[1][(i + 1) % 3] * arr[2][(i + 2) % 3] - arr[1][(i + 2) % 3] * arr[2][(i + 1) % 3]);
    }
    return ans;
}

class GeneralObject : public Object
{
public:
    Point center;
    double a, b, c, d, e, f, g, h, i, j;
    Point normal;
    GeneralObject(){};

    virtual Ray getNormal(Point point, Ray i_ray)
    {
        Point normal = Point(2 * a * point.values.x + d * point.values.y + e * point.values.z + g, 2 * b * point.values.y + d * point.values.x + f * point.values.z + h, 2 * c * point.values.z + e * point.values.x + f * point.values.y + i);
        return Ray(point, normal);
    }

    bool isOkay(Point point)
    {
        bool lengthCheck = (dimentions[LENGTH] > 1e-5);
        bool widthCheck = (dimentions[WIDTH] > 1e-5);
        bool heightCheck = (dimentions[HEIGHT] > 1e-5);

        if ((lengthCheck && point.values.x < ref_point.values.x) || (lengthCheck && point.values.x > ref_point.values.x + dimentions[LENGTH]))
            return false;

        if ((widthCheck && point.values.y < ref_point.values.y) || (widthCheck && point.values.y > ref_point.values.y + dimentions[WIDTH]))
            return false;

        if ((heightCheck && point.values.z < ref_point.values.z) || (heightCheck && point.values.z > ref_point.values.z + dimentions[HEIGHT]))
            return false;

        return true;
    }

    virtual void draw()
    {
        return;
    }

    virtual double findIntersectionDistance(Ray ray)
    {
        Point _0 = Point(ray.origin.values.x, ray.origin.values.y, ray.origin.values.z);
        Point _1 = Point(ray.direction.values.x, ray.direction.values.y, ray.direction.values.z);

        Point C = Point(
            a * _1.values.x * _1.values.x + b * _1.values.y * _1.values.y + c * _1.values.z * _1.values.z + d * _1.values.x * _1.values.y + e * _1.values.x * _1.values.z + f * _1.values.z * _1.values.y,
            2 * a * _0.values.x * _1.values.x + 2 * b * _0.values.y * _1.values.y + 2 * c * _0.values.z * _1.values.z + d * (_0.values.x * _1.values.y + _0.values.y * _1.values.x) + e * (_0.values.x * _1.values.z + _0.values.z * _1.values.x) + f * (_0.values.z * _1.values.y + _0.values.y * _1.values.z) + g * _1.values.x + h * _1.values.y + i * _1.values.z,
            a * _0.values.x * _0.values.x + b * _0.values.y * _0.values.y + c * _0.values.z * _0.values.z + d * _0.values.x * _0.values.y + e * _0.values.x * _0.values.z + f * _0.values.z * _0.values.y + g * _0.values.x + h * _0.values.y + i * _0.values.z + j);

        double discriminant = C.values.y * C.values.y - 4 * C.values.x * C.values.z;
        if (discriminant < 0)
            return -1;

        bool isCLess = (fabs(C.values.x) < 1e-5);
        if (isCLess)
            return -1;

        // pair<double, double> t = make_pair((-C.values.y - sqrt(discriminant)) / (2 * C.values.x), (-C.values.y + sqrt(discriminant)) / (2 * C.values.x));
        pair<double, double> t((-C.values.y - sqrt(discriminant)) / (2 * C.values.x),
                               (-C.values.y + sqrt(discriminant)) / (2 * C.values.x));

        bool isAllSolutionNeg = (t.first < 0 && t.second < 0);
        if (isAllSolutionNeg)
            return -1;

        if (t.first > t.second)
            swap(t.first, t.second);

        bool isTFirstPos = (t.first > 0);
        bool isTSecondPos = (t.second > 0);

        Point intersectPointFirst = ray.origin.add(ray.direction.multiply(t.first));
        Point intersectPointSecond = ray.origin.add(ray.direction.multiply(t.second));

        if (isOkay(intersectPointFirst))
        {
            return t.first;
        }

        if (isOkay(intersectPointSecond))
        {
            return t.second;
        }

        return -1;
    }

    friend istream &operator>>(istream &input, GeneralObject &generalObject)
    {
        input >> generalObject.a >> generalObject.b >> generalObject.c >> generalObject.d >> generalObject.e >> generalObject.f >> generalObject.g >> generalObject.h >> generalObject.i >> generalObject.j;
        input >> generalObject.ref_point.values.x >> generalObject.ref_point.values.y >> generalObject.ref_point.values.z;
        input >> generalObject.dimentions[LENGTH] >> generalObject.dimentions[WIDTH] >> generalObject.dimentions[HEIGHT];
        input >> generalObject.color.red >> generalObject.color.green >> generalObject.color.blue;

        input >> generalObject.coeff[0] >> generalObject.coeff[1] >> generalObject.coeff[2] >> generalObject.coeff[3];
        input >> generalObject.shine;
        return input;
    }
};

class Triangle : public Object
{
public:
    Point a, b, c;

    Triangle()
    {
        a = Point(0, 0, 0);
        b = Point(0, 0, 0);
        c = Point(0, 0, 0);
    }
    Triangle(Point A, Point B, Point C)
    {
        a = A;
        b = B;
        c = C;
    }

    virtual Ray getNormal(Point point, Ray i_ray)
    {
        Point normal = (b.subtract(a)).cross(c.subtract(a));
        normal = normal.normalize();

        double dotProduct = normal.dot(i_ray.direction);

        if (dotProduct > 0)
            return Ray(point, normal);
        else
            return Ray(point, normal.multiply(-1));
    }

    virtual void draw()
    {
        glColor3f(color.red, color.green, color.blue);
        glBegin(GL_TRIANGLES);
        glVertex3f(a.values.x, a.values.y, a.values.z);
        glVertex3f(b.values.x, b.values.y, b.values.z);
        glVertex3f(c.values.x, c.values.y, c.values.z);
        glEnd();
    }

    virtual double findIntersectionDistance(Ray ray)
    {
        Point tempOne = a.subtract(ray.origin);
        Point tempTwo = a.subtract(c);
        Point tempThree = ray.direction;

        double betaMatrix[3][3] = {
            {tempOne.values.x, tempTwo.values.x, tempThree.values.x},
            {tempOne.values.y, tempTwo.values.y, tempThree.values.y},
            {tempOne.values.z, tempTwo.values.z, tempThree.values.z}};

        tempOne = a.subtract(b);
        tempTwo = a.subtract(ray.origin);
        tempThree = ray.direction;

        double gammaMatrix[3][3] = {
            {tempOne.values.x, tempTwo.values.x, tempThree.values.x},
            {tempOne.values.y, tempTwo.values.y, tempThree.values.y},
            {tempOne.values.z, tempTwo.values.z, tempThree.values.z}};

        tempOne = a.subtract(b);
        tempTwo = a.subtract(c);
        tempThree = a.subtract(ray.origin);

        double tMatrix[3][3] = {
            {tempOne.values.x, tempTwo.values.x, tempThree.values.x},
            {tempOne.values.y, tempTwo.values.y, tempThree.values.y},
            {tempOne.values.z, tempTwo.values.z, tempThree.values.z}};

        tempOne = a.subtract(b);
        tempTwo = a.subtract(c);
        tempThree = ray.direction;

        double aMatrix[3][3] = {
            {tempOne.values.x, tempTwo.values.x, tempThree.values.x},
            {tempOne.values.y, tempTwo.values.y, tempThree.values.y},
            {tempOne.values.z, tempTwo.values.z, tempThree.values.z}};

        double detA = det(aMatrix);
        double detBeta = det(betaMatrix);
        double detGamma = det(gammaMatrix);
        double detT = det(tMatrix);

        if (detA == 0)
            return -1;

        double beta = detBeta / detA;
        double gamma = detGamma / detA;
        double t = detT / detA;

        bool betaAndGammaInRange = (beta > 0 && gamma > 0 && beta + gamma <= 1);
        bool tInRange = (t > 0);

        if (betaAndGammaInRange && tInRange)
        {
            return t;
        }
        else
        {
            return -1;
        }
    }

    friend istream &operator>>(istream &input, Triangle &triangle)
    {
        input >> triangle.a.values.x >> triangle.a.values.y >> triangle.a.values.z;
        input >> triangle.b.values.x >> triangle.b.values.y >> triangle.b.values.z;
        input >> triangle.c.values.x >> triangle.c.values.y >> triangle.c.values.z;
        input >> triangle.color.red >> triangle.color.green >> triangle.color.blue;
        input >> triangle.coeff[0] >> triangle.coeff[1] >> triangle.coeff[2] >> triangle.coeff[3];
        input >> triangle.shine;
        return input;
    }
};

class Sphere : public Object
{
public:
    Sphere() {}
    Sphere(Point center, double radius)
    {
        this->ref_point = center;
        this->dimentions[LENGTH] = radius;
    }

    virtual Ray getNormal(Point point, Ray i_ray)
    {
        Point normal = point.subtract(ref_point);

        // eta check kora lagbe
        normal = normal.normalize();
        return Ray(point, normal);
    }
    Point points[100][100];

    void drawHemisphere(int pos, int i, int j)
    {
        // glVertex3f(points[i][j].values.x, points[i][j].values.y, points[i][j].values.z);
        // glVertex3f(points[i][j + 1].values.x, points[i][j + 1].values.y, points[i][j + 1].values.z);
        // glVertex3f(points[i + 1][j + 1].values.x, points[i + 1][j + 1].values.y, points[i + 1][j + 1].values.z);
        // glVertex3f(points[i + 1][j].values.x, points[i + 1][j].values.y, points[i + 1][j].values.z);

        glVertex3f(points[i][j].values.x, points[i][j].values.y, points[i][j].values.z * pos);
        glVertex3f(points[i][j + 1].values.x, points[i][j + 1].values.y, points[i][j + 1].values.z * pos);
        glVertex3f(points[i + 1][j + 1].values.x, points[i + 1][j + 1].values.y, points[i + 1][j + 1].values.z * pos);
        glVertex3f(points[i + 1][j].values.x, points[i + 1][j].values.y, points[i + 1][j].values.z * pos);
    }
    virtual void draw()
    {
        int stacks = 30;
        int slices = 20;

        int i, j;
        double h, r;
        i = 0;
        while (i <= stacks)
        {
            double Theta = ((double)i / (double)stacks) * (pi / 2);
            h = dimentions[LENGTH] * sin(Theta);
            r = dimentions[LENGTH] * cos(Theta);

            j = 0;
            while (j <= slices)
            {
                double Phi = ((double)j / (double)slices) * 2 * pi;
                points[i][j].values.x = r * cos(Phi);
                points[i][j].values.y = r * sin(Phi);
                points[i][j].values.z = h;
                j++;
            }
            i++;
        }

        i = 0;
        while (i < stacks)
        {
            glPushMatrix();
            glTranslatef(ref_point.values.x, ref_point.values.y, ref_point.values.z);
            glColor3f(color.red, color.green, color.blue);

            j = 0;
            while (j < slices)
            {
                glBegin(GL_QUADS);
                // upper hemisphere
                drawHemisphere(1, i, j);

                drawHemisphere(-1, i, j);
                glEnd();
                j++;
            }
            glPopMatrix();
            i++;
        }
    }

    virtual double findIntersectionDistance(Ray ray)
    {
        ray.origin = ray.origin.subtract(ref_point);

        double a = 1;
        double b = 2 * ray.direction.dot(ray.origin);
        double c = ray.origin.dot(ray.origin) - dimentions[LENGTH] * dimentions[LENGTH];

        double discriminant = b * b - 4 * a * c;
        double t = -1;

        if (discriminant < 0)
        {
            t = -1;
        }
        else
        {
            // pair<double, double> roots = make_pair((-b - sqrt(discriminant)) / (2 * a), (-b + sqrt(discriminant)) / (2 * a));
            pair<double, double> roots((-b - sqrt(discriminant)) / (2 * a), (-b + sqrt(discriminant)) / (2 * a));

            if (roots.first > roots.second)
            {
                swap(roots.first, roots.second);
            }

            (roots.first > 0) ? t = roots.first : t = roots.second;
            (t > 0) ? t = t : t = -1;
        }
        return t;
    }

    friend istream &operator>>(istream &input, Sphere &sphere)
    {
        input >> sphere.ref_point.values.x >> sphere.ref_point.values.y >> sphere.ref_point.values.z;
        input >> sphere.dimentions[LENGTH];
        input >> sphere.color.red >> sphere.color.green >> sphere.color.blue;
        input >> sphere.coeff[0] >> sphere.coeff[1] >> sphere.coeff[2] >> sphere.coeff[3];
        input >> sphere.shine;
        return input;
    }
};

class Floor : public Object
{
public:
    int tileNumber;
    Floor()
    {
        tileNumber = 1;
    }

    Floor(int fw, int tw)
    {
        tileNumber = fw / tw;
        ref_point = Point(-fw / 2, -fw / 2, 0);
        dimentions[LENGTH] = tw;
    }

    virtual Color getColorAt(Point point)
    {
        // pair<int, int> tile = make_pair((point.values.x - ref_point.values.x) / dimentions[LENGTH], (point.values.y - ref_point.values.y) / dimentions[LENGTH]);
        pair<int, int> tile = make_pair(
            static_cast<int>((point.values.x - ref_point.values.x) / dimentions[LENGTH]),
            static_cast<int>((point.values.y - ref_point.values.y) / dimentions[LENGTH]));

        bool isTileXInRange = (tile.first >= 0 && tile.first < tileNumber);
        bool isTileYInRange = (tile.second >= 0 && tile.second < tileNumber);

        if (!isTileXInRange || !isTileYInRange)
        {
            return Color(0, 0, 0);
        }

        if ((tile.first + tile.second) % 2 == 0)
        {
            return Color(1, 1, 1);
        }
        else
        {
            return Color(0, 0, 0);
        }
    }

    virtual Ray getNormal(Point point, Ray i_ray)
    {
        return Ray(point, Point(0, 0, i_ray.direction.values.z > 0 ? 1 : -1));
    }

    virtual void draw()
    {

        int i = 0;
        while (i < tileNumber)
        {
            float xOrigin = ref_point.values.x + i * dimentions[LENGTH];
            float xNextOrigin = ref_point.values.x + (i + 1) * dimentions[LENGTH];

            int j = 0;
            while (j < tileNumber)
            {
                float yOrigin = ref_point.values.y + j * dimentions[LENGTH];
                float yNextOrigin = ref_point.values.y + (j + 1) * dimentions[LENGTH];

                // Use the ternary operator for color selection
                ((i + j) % 2 == 0) ? glColor3f(1, 1, 1) : glColor3f(0, 0, 0); // White : Black

                glBegin(GL_QUADS);
                {
                    glVertex3f(xOrigin, yOrigin, 0);
                    glVertex3f(xNextOrigin, yOrigin, 0);
                    glVertex3f(xNextOrigin, yNextOrigin, 0);
                    glVertex3f(xOrigin, yNextOrigin, 0);
                }
                glEnd();

                j++;
            }
            i++;
        }
    }

    virtual double findIntersectionDistance(Ray ray)
    {
        Point normal = Point(0, 0, 1);
        double dotProduct = normal.dot(ray.direction) * 100;

        if (round(dotProduct) == 0)
        {
            return -1;
        }

        dotProduct /= 100;

        double t = -1 * (ray.origin.dot(normal)) / dotProduct;

        Point p = ray.origin.add(ray.direction.multiply(t));

        bool firstCond = (p.values.x <= ref_point.values.x);
        bool secondCond = (p.values.x >= abs(ref_point.values.x));
        bool thirdCond = (p.values.y <= ref_point.values.y);
        bool fourthCond = (p.values.y >= abs(ref_point.values.y));

        if (firstCond || secondCond && thirdCond && fourthCond)
        {
            return -1;
        }
        return t;
    }
};
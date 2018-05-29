/*
 * Figure3D.h
 *
 *  Created on: Mar 7, 2018
 *      Author: uauser
 */

#ifndef FIGURE3D_H_
#define FIGURE3D_H_

#include <vector>
#include <fstream>
#include <stack>

#include "vector3d.h"
#include "lines2D.h"
#include "l_parser.h"

using namespace std;

class Color;

class Face {
private:
	vector<int> pointIndexes; //waardes refereren naar de indexen van de punten in in point in een Figure class
public:
	const vector<int>& getPointIndexes() const;

	Face(vector<int> pIndexes);

	Face();

	void addPoint(const int pointIndex);
};


class Figure3D {
private:
	vector<Matrix> transformationsToApply;

public:
	vector<Vector3D> points;
	vector<Face*> faces;
	Color* ambientReflection;
	Color* diffusereflection;

	Figure3D(vector<Vector3D>& pointsVector, vector<Face*>& faceVector, Color*& ambientRefl, Color*& diffuseRefl);

	Figure3D(vector<Vector3D>& pointsVector, vector<Face*>& faceVector, Color*& ambientRefl);

	void scaleFigure(const double scale);

	void rotateX(const double angle);

	void rotateY(const double angle);

	void rotateZ(const double angle);

	void translate(const Vector3D& vector);

	void applyTransformations();

	vector<Face*> triangulateFace(Face*& face);
};

class Figures3D {
private:
	vector<Figure3D*> figures;
	double eyeTheta;
	double eyePhi;
	double eyeR;

	void eyeToPolar(const Vector3D eye);

	Point2D* doPointProjection(const double d, const Vector3D& point);

	string replaceRule(string currentRule, unsigned int currentIteration, LParser::LSystem3D& lSystem);
public:
	Matrix eyepointTrans(const Vector3D& eyepoint);

	void applyTransformations(Matrix& eyepointMatrix);

	Lines2D doProjection();

	const vector<Figure3D*>& getFigures() const {
		return figures;
	}

	void addFigure(Figure3D* figure);

	Figure3D* createTetrahedron(Color*& col);

	Figure3D* createCube(Color*& col);

	Figure3D* createOctahedron(Color*& col);

	Figure3D* createIcosahedron(Color*& col);

	Figure3D* createDodecahedron(Color*& col);

	Figure3D* createSphere(const int n, Color*& col);

	Figure3D* createCone(const int n, const double h, Color*& col);

	Figure3D* createCylinder(const int n, const double h, Color*& col);

	Figure3D* createTorus(const double R, const double r, const int n, const int m, Color*& col);

	Figure3D* create3DLSystem(string inputfile, Color*& col);

	void deleteFigure(Figure3D*& figureToDelete);
};

#endif /* FIGURE3D_H_ */

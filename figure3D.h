/*
 * Figure3D.h
 *
 *  Created on: Mar 7, 2018
 *      Author: uauser
 */

#ifndef FIGURE3D_H_
#define FIGURE3D_H_

#include <vector>

#include "vector3d.h"
#include "lines2D.h"

using namespace std;

class Color;

class Face {
private:
	vector<int> pointIndexes; //waardes refereren naar de indexen van de punten in in point in een Figure class
public:
	const vector<int>& getPointIndexes() const;

	Face(vector<int> pIndexes);
};


class Figure3D {
private:
	vector<Matrix> transformationsToApply;

public:
	vector<Vector3D> points;
	vector<Face*> faces;
	Color* color;

	Figure3D(vector<Vector3D>& pointsVector, vector<Face*>& faceVector, Color*& color);

	void scaleFigure(const double scale);

	void rotateX(const double angle);

	void rotateY(const double angle);

	void rotateZ(const double angle);

	void translate(const Vector3D& vector);

	void applyTransformations();
};

class Figures3D {
private:
	vector<Figure3D*> figures;
	double eyeTheta;
	double eyePhi;
	double eyeR;

	void eyeToPolar(const Vector3D eye);

	Point2D* doPointProjection(const double d, const Vector3D& point);
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
};

#endif /* FIGURE3D_H_ */

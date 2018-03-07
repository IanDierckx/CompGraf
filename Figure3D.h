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

using namespace std;

class Face {
private:
	vector<int> pointIndexes; //waardes refereren naar de indexen van de punten in in point in een Figure class
public:
	const vector<int>& getPointIndexes() const;

	Face(vector<int> pIndexes);
};


class Figure3D {
private:
	vector<Vector3D> points;
	vector<Face*> faces;
	vector<Matrix> transformationsToApply;
public:
	Figure3D(vector<Vector3D> pointsVector, vector<Face*> faceVector);

	const vector<Face*>& getFaces() const;

	const vector<Vector3D>& getPoints() const;

	void scaleFigure(const double scale);

	void rotateX(const double angle);

	void rotateY(const double angle);

	void rotateZ(const double angle);

	void translate(const Vector3D& vector);

	void applyTransformations();
};

class Figures3D {
public:
	vector<Figure3D*> figures;
};

#endif /* FIGURE3D_H_ */

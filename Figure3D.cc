/*
 * Figure3D.cc
 *
 *  Created on: Mar 7, 2018
 *      Author: uauser
 */

#include "Figure3D.h"
#include <cmath>


const vector<int>& Face::getPointIndexes() const {
	return pointIndexes;
}

Face::Face(vector<int> pIndexes) {
	pointIndexes = pIndexes;
}

const vector<Face*>& Figure3D::getFaces() const {
	return faces;
}

Figure3D::Figure3D(vector<Vector3D> pointsVector, vector<Face*> faceVector) {
	points = pointsVector;
	faces = faceVector;
	transformationsToApply = {};
}

const vector<Vector3D>& Figure3D::getPoints() const {
	return points;
}

void Figure3D::scaleFigure(const double scale) {
	Matrix scaleMatrx;
	scaleMatrx(1,1) = scale;
	scaleMatrx(2,2) = scale;
	scaleMatrx(3,3) = scale;
	transformationsToApply.push_back(scaleMatrx);
}

void Figure3D::rotateY(const double angle) {
	Matrix rotYMatrx;
	double anglerad = angle * (M_PI/180);
	rotYMatrx(1,1) = cos(anglerad);
	rotYMatrx(3,1) = sin(anglerad);
	rotYMatrx(1,3) = -sin(anglerad);
	rotYMatrx(3,3) = cos(angle);
	transformationsToApply.push_back(rotYMatrx);
}

void Figure3D::rotateX(const double angle) {
	Matrix rotXMatrx;
	double anglerad = angle * (M_PI/180);
	rotXMatrx(2,2) = cos(anglerad);
	rotXMatrx(2,3) = sin(anglerad);
	rotXMatrx(3,2) = -sin(anglerad);
	rotXMatrx(3,3) = cos(anglerad);
	transformationsToApply.push_back(rotXMatrx);
}

void Figure3D::rotateZ(const double angle) {
	Matrix rotZMatrx;
	double anglerad = angle * (M_PI/180);
	rotZMatrx(1,1) = cos(anglerad);
	rotZMatrx(1,2) = sin(anglerad);
	rotZMatrx(2,1) = -sin(anglerad);
	rotZMatrx(2,2) = cos(anglerad);
	transformationsToApply.push_back(rotZMatrx);
}

void Figure3D::translate(const Vector3D& vector) {
	if (vector.is_vector()) {
		Matrix transMatrx;
		transMatrx(4,1) = vector.x;
		transMatrx(4,2) = vector.y;
		transMatrx(4,3) = vector.z;
		transformationsToApply.push_back(transMatrx);
	} else {
		cerr << "Trying to translate using a point instead of a vector" << endl;
	}
}

void Figure3D::applyTransformations() {
	Matrix transformatieMatrx = transformationsToApply[0];
	if (transformationsToApply.size() > 1) {
		for (unsigned int index = 0; index < transformationsToApply.size(); ++index) {
			transformatieMatrx *= transformationsToApply[index];
		}
	}
	for (unsigned int puntIndex = 0; puntIndex < points.size(); ++puntIndex) {
		points[puntIndex] *= transformatieMatrx;
	}
	transformationsToApply = {};
}

/*
 * Figure3D.cc
 *
 *  Created on: Mar 7, 2018
 *      Author: uauser
 */

#include <cmath>
#include "figure3D.h"


const vector<int>& Face::getPointIndexes() const {
	return pointIndexes;
}

Face::Face(vector<int> pIndexes) {
	pointIndexes = pIndexes;
}

Figure3D::Figure3D(vector<Vector3D> pointsVector, vector<Face*> faceVector, Color*& col) {
	points = pointsVector;
	faces = faceVector;
	color = col;
	transformationsToApply = {};
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

Matrix Figures3D::eyepointTrans(const Vector3D& eyepoint) {
	eyeToPolar(eyepoint);
	Matrix eyePointTransMatrx;
	eyePointTransMatrx(1,1) = -sin(eyeTheta);
	eyePointTransMatrx(1,2) = -cos(eyeTheta)*cos(eyePhi);
	eyePointTransMatrx(1,3) = cos(eyeTheta)*sin(eyePhi);
	eyePointTransMatrx(2,1) = cos(eyeTheta);
	eyePointTransMatrx(2,2) = -sin(eyeTheta)*cos(eyePhi);
	eyePointTransMatrx(2,3) = sin(eyeTheta)*sin(eyePhi);
	eyePointTransMatrx(3,2) = sin(eyePhi);
	eyePointTransMatrx(3,3) = cos(eyePhi);
	eyePointTransMatrx(4,3) = -eyeR;
	return eyePointTransMatrx;
}

void Figures3D::eyeToPolar(const Vector3D eye) {
	eyeR = sqrt((pow(eye.x,2)+pow(eye.y,2)+pow(eye.z,2)));
	eyeTheta = atan2(eye.y,eye.x);
	eyePhi = acos((eye.z/eyeR));
}

void Figures3D::applyTransformations(Matrix eyepointMatrix) {
	for (auto figure:figures) {
		figure->applyTransformations();
		for (unsigned int puntIndex = 0; puntIndex < figure->points.size(); ++puntIndex) {
			figure->points[puntIndex] *= eyepointMatrix;
		}
	}
}

Point2D* Figures3D::doPointProjection(const double d, const Vector3D& point) {
	double x = (d*point.x)/-point.z;
	double y = (d*point.y)/-point.z;
	Point2D* projectedPoint = new Point2D(x,y);
	return projectedPoint;
}

Lines2D Figures3D::doProjection() {
	vector<Line2D*> lines;
	for (auto figure:figures) {
		for (auto face:figure->faces) {
			unsigned int currentIndex = 0;
			while (currentIndex < face->getPointIndexes().size()-1) {
				Point2D* p1 = doPointProjection(1,figure->points[face->getPointIndexes()[currentIndex]]);
				unsigned int index2ePunt = currentIndex+1;
				while (index2ePunt < face->getPointIndexes().size()) {
					Point2D* p2 = doPointProjection(1,figure->points[face->getPointIndexes()[index2ePunt]]);
					Line2D* lijn = new Line2D(p1,p2,figure->color);
					lines.push_back(lijn);
					index2ePunt += 1;
				}
				currentIndex += 1;
			}
		}
	}
	Lines2D lijnen = Lines2D(lines);
	return lijnen;
}

void Figures3D::addFigure(Figure3D* figure) {
	figures.push_back(figure);
}

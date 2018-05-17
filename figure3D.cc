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

Figure3D::Figure3D(vector<Vector3D>& pointsVector, vector<Face*>& faceVector, Color*& col) {
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
	rotYMatrx(3,3) = cos(anglerad);
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
		for (unsigned int index = 1; index < transformationsToApply.size(); ++index) {
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

void Figures3D::applyTransformations(Matrix& eyepointMatrix) {
	for (auto& figure:figures) {
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
	for (auto& figure:figures) {
		for (auto& face:figure->faces) {
			unsigned int currentIndex = 0;
			while (currentIndex < face->getPointIndexes().size()-1) {
				Point2D* p1 = doPointProjection(1,figure->points[face->getPointIndexes()[currentIndex]]);
				unsigned int index2ePunt = currentIndex+1;
				Point2D* p2 = doPointProjection(1,figure->points[face->getPointIndexes()[index2ePunt]]);
				Line2D* lijn = new Line2D(p1,p2,figure->color);
				lines.push_back(lijn);
				currentIndex += 1;
			}
			Point2D* p1 = doPointProjection(1,figure->points[face->getPointIndexes()[currentIndex]]);
			unsigned int index2ePunt = currentIndex+1;
			Point2D* p2 = doPointProjection(1,figure->points[face->getPointIndexes()[0]]);
			Line2D* lijn = new Line2D(p1,p2,figure->color);
			lines.push_back(lijn);
		}
	}
	Lines2D lijnen = Lines2D(lines);
	return lijnen;
}

void Figures3D::addFigure(Figure3D* figure) {
	figures.push_back(figure);
}

Figure3D* Figures3D::createTetrahedron(Color*& col) {
	Vector3D p1 = Vector3D::point(1.0,-1.0,-1.0);
	Vector3D p2 = Vector3D::point(-1.0,1.0,-1.0);
	Vector3D p3 = Vector3D::point(1.0,1.0,1.0);
	Vector3D p4 = Vector3D::point(-1.0,-1.0,1.0);

	vector<Vector3D> points = {p1,p2,p3,p4};

	Face* face1 = new Face({0,1,2});
	Face* face2 = new Face({1,3,2});
	Face* face3 = new Face({0,3,1});
	Face* face4 = new Face({0,2,3});

	vector<Face*> faces = {face1,face2,face3,face4};

	Figure3D* tetrahedron = new Figure3D(points,faces,col);
	return tetrahedron;
}

Figure3D* Figures3D::createCube(Color*& col) {
	Vector3D p1 = Vector3D::point(1.0,-1.0,-1.0);
	Vector3D p2 = Vector3D::point(-1.0,1.0,-1.0);
	Vector3D p3 = Vector3D::point(1.0,1.0,1.0);
	Vector3D p4 = Vector3D::point(-1.0,-1.0,1.0);
	Vector3D p5 = Vector3D::point(1.0,1.0,-1.0);
	Vector3D p6 = Vector3D::point(-1.0,-1.0,-1.0);
	Vector3D p7 = Vector3D::point(1.0,-1.0,1.0);
	Vector3D p8 = Vector3D::point(-1.0,1.0,1.0);

	vector<Vector3D> points = {p1,p2,p3,p4,p5,p6,p7,p8};

	Face* face1 = new Face({0,4,2,6});
	Face* face2 = new Face({4,1,7,2});
	Face* face3 = new Face({1,5,3,7});
	Face* face4 = new Face({5,0,6,3});
	Face* face5 = new Face({6,2,7,3});
	Face* face6 = new Face({0,5,1,4});

	vector<Face*> faces = {face1,face2,face3,face4,face5,face6};

	Figure3D* cube = new Figure3D(points,faces,col);
	return cube;
}

Figure3D* Figures3D::createOctahedron(Color*& col) {
	Vector3D p1 = Vector3D::point(1.0,0.0,0.0);
	Vector3D p2 = Vector3D::point(0.0,1.0,0.0);
	Vector3D p3 = Vector3D::point(-1.0,0.0,0.0);
	Vector3D p4 = Vector3D::point(0.0,-1.0,0.0);
	Vector3D p5 = Vector3D::point(0.0,0.0,-1.0);
	Vector3D p6 = Vector3D::point(0.0,0.0,1.0);

	vector<Vector3D> points = {p1,p2,p3,p4,p5,p6};

	Face* face1 = new Face({0,1,5});
	Face* face2 = new Face({1,2,5});
	Face* face3 = new Face({2,3,5});
	Face* face4 = new Face({3,0,5});
	Face* face5 = new Face({1,0,4});
	Face* face6 = new Face({2,1,4});
	Face* face7 = new Face({3,2,4});
	Face* face8 = new Face({0,3,4});

	vector<Face*> faces = {face1,face2,face3,face4,face5,face6,face7,face8};

	Figure3D* octahedron = new Figure3D(points,faces,col);
	return octahedron;
}

Figure3D* Figures3D::createIcosahedron(Color*& col) {
	vector<Vector3D> points;
	for (int index = 1; index <= 12; ++index) {
		if (index == 1) {
			Vector3D point = Vector3D::point(0.0,0.0,sqrt(5)/2.0);
			points.push_back(point);
		} else if (index <= 6) {
			Vector3D point  = Vector3D::point(cos((index-2)*2*M_PI/5),sin((index-2)*2*M_PI/5),0.5);
			points.push_back(point);
		} else if (index <= 11) {
			Vector3D point = Vector3D::point(cos(M_PI/5+(index-7)*2*M_PI/5),sin(M_PI/5+(index-7)*2*M_PI/5),-0.5);
			points.push_back(point);
		} else {
			Vector3D point = Vector3D::point(0.0,0.0,-sqrt(5)/2);
			points.push_back(point);
		}
	}

	Face* face1 = new Face({0,1,2});
	Face* face2 = new Face({0,2,3});
	Face* face3 = new Face({0,3,4});
	Face* face4 = new Face({0,4,5});
	Face* face5 = new Face({0,5,1});
	Face* face6 = new Face({1,6,2});
	Face* face7 = new Face({2,6,7});
	Face* face8 = new Face({2,7,3});
	Face* face9 = new Face({3,7,8});
	Face* face10 = new Face({3,8,4});
	Face* face11 = new Face({4,8,9});
	Face* face12 = new Face({4,9,5});
	Face* face13 = new Face({5,9,10});
	Face* face14 = new Face({5,10,1});
	Face* face15 = new Face({1,10,6});
	Face* face16 = new Face({11,7,6});
	Face* face17 = new Face({11,8,7});
	Face* face18 = new Face({11,9,8});
	Face* face19 = new Face({11,10,9});
	Face* face20 = new Face({11,6,10});

	vector<Face*> faces = {face1,face2,face3,face4,face5,face6,face7,face8,face9,face10,face11,face12,face13,face14,
		face15,face16,face17,face18,face19,face20};

	Figure3D* icosahedron = new Figure3D(points,faces,col);
	return icosahedron;
}

Figure3D* Figures3D::createDodecahedron(Color*& col) {
	Figure3D* icosahedron = createIcosahedron(col);
	vector<Vector3D> points;
	for (auto face:icosahedron->faces) {
		Vector3D p1 = icosahedron->points[face->getPointIndexes()[0]];
		Vector3D p2 = icosahedron->points[face->getPointIndexes()[1]];
		Vector3D p3 = icosahedron->points[face->getPointIndexes()[2]];
		double midX = (p1.x+p2.x+p3.x)/3;
		double midY = (p1.y+p2.y+p3.y)/3;
		double midZ = (p1.z+p2.z+p3.z)/3;
		Vector3D newPoint = Vector3D::point(midX,midY,midZ);
		points.push_back(newPoint);
		delete face;
	}
	delete icosahedron;

	Face* face1 = new Face({0,1,2,3,4});
	Face* face2 = new Face({0,5,6,7,1});
	Face* face3 = new Face({1,7,8,9,2});
	Face* face4 = new Face({2,9,10,11,3});
	Face* face5 = new Face({3,11,12,13,4});
	Face* face6 = new Face({4,13,14,5,0});
	Face* face7 = new Face({19,18,17,16,15});
	Face* face8 = new Face({19,14,13,12,18});
	Face* face9 = new Face({18,12,11,10,17});
	Face* face10 = new Face({17,10,9,8,16});
	Face* face11 = new Face({16,8,7,6,15});
	Face* face12 = new Face({15,6,5,14,19});

	vector<Face*> faces = {face1,face2,face3,face4,face5,face6,face7,face8,face9,face10,face11,face12};

	Figure3D* dodecahedron = new Figure3D(points,faces,col);
	return dodecahedron;
}

Figure3D* Figures3D::createSphere(const int n, Color*& col) {
	Figure3D* icosahedron = createIcosahedron(col);
	vector<Face*> newFaces = icosahedron->faces;
	vector<Vector3D> newPoints = icosahedron->points;
	for (int times = 1; times <= n; ++times) {
		vector<Face*> newNewFaces;
		for (auto& face:newFaces) {
			Vector3D p1 = newPoints[face->getPointIndexes()[0]];
			Vector3D p2 = newPoints[face->getPointIndexes()[1]];
			Vector3D p3 = newPoints[face->getPointIndexes()[2]];
			Vector3D p4 = Vector3D::point((p1.x+p3.x)/2,(p1.y+p3.y)/2,(p1.z+p3.z)/2);
			Vector3D p5 = Vector3D::point((p1.x+p2.x)/2,(p1.y+p2.y)/2,(p1.z+p2.z)/2);
			Vector3D p6 = Vector3D::point((p2.x+p3.x)/2,(p2.y+p3.y)/2,(p2.z+p3.z)/2);

			int p1Index = face->getPointIndexes()[0];
			int p2Index = face->getPointIndexes()[1];
			int p3Index = face->getPointIndexes()[2];
			int p4Index;
			int p5Index;
			int p6Index;
			bool p4Found = false;
			bool p5Found = false;
			bool p6Found = false;
			for (unsigned int pointIndex = 0; pointIndex < newPoints.size(); ++pointIndex) {
				if ((newPoints[pointIndex].x == p4.x) and (newPoints[pointIndex].y == p4.y) and (newPoints[pointIndex].z == p4.z)) {
					p4Found = true;
					p4Index = pointIndex;
				} else if ((newPoints[pointIndex].x == p5.x) and (newPoints[pointIndex].y == p5.y) and (newPoints[pointIndex].z == p5.z)) {
					p5Found = true;
					p5Index = pointIndex;
				} else if ((newPoints[pointIndex].x == p6.x) and (newPoints[pointIndex].y == p6.y) and (newPoints[pointIndex].z == p6.z)) {
					p6Found = true;
					p6Index = pointIndex;
				}
			}
			if (not p4Found) {
				newPoints.push_back(p4);
				p4Index = newPoints.size()-1;
			}
			if (not p5Found) {
				newPoints.push_back(p5);
				p5Index = newPoints.size()-1;
			}
			if (not p6Found) {
				newPoints.push_back(p6);
				p6Index = newPoints.size()-1;
			}

			Face* newFace1 = new Face({p1Index,p5Index,p4Index});
			newNewFaces.push_back(newFace1);
			Face* newFace2 = new Face({p4Index,p5Index,p6Index});
			newNewFaces.push_back(newFace2);
			Face* newFace3 = new Face({p5Index,p2Index,p6Index});
			newNewFaces.push_back(newFace3);
			Face* newFace4 = new Face({p4Index,p6Index,p3Index});
			newNewFaces.push_back(newFace4);

			delete face;
		}
		newFaces = newNewFaces;
	}
	delete icosahedron;
	for (auto& point:newPoints) {
		point.normalise();
	}
	Figure3D* sphere = new Figure3D(newPoints,newFaces,col);
	return sphere;
}

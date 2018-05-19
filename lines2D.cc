/*
 * Lines2D.cpp
 *
 *  Created on: Feb 20, 2018
 *      Author: uauser
 */

#include "lines2D.h"


using namespace std;

class Color;
class Point2D;
class Line2D;

Color::Color(double r, double g, double b) {
	red=r;
	green=g;
	blue=b;
}

Point2D::Point2D(double inx, double iny) {
	x = inx;
	y = iny;
}

Line2D::Line2D(Point2D*& point1, Point2D*& point2, Color*& col) {
	p1 = point1;
	p2 = point2;
	color = col;
	z1 = 0;
	z2 = 0;
}

void Lines2D::addLine(Line2D* line2D) {
	this->Lines.push_back(line2D);
}


Lines2D::Lines2D(vector<Line2D*> lines) {
	Lines = lines;
}

vector<double> Lines2D::getMinMax() {
	double maxX;
	double maxY;
	double minX;
	double minY;
	vector<Line2D*> lines = this->getLines();
	vector<double> xen;
	vector<double> yen;
	for (unsigned int index = 0; index < lines.size(); index++) {
		xen.push_back(lines[index]->p1->x);
		xen.push_back(lines[index]->p2->x);
		yen.push_back(lines[index]->p1->y);
		yen.push_back(lines[index]->p2->y);
	}
	maxX = *max_element(xen.begin(), xen.end());
	maxY = *max_element(yen.begin(), yen.end());
	minX = *min_element(xen.begin(), xen.end());
	minY = *min_element(yen.begin(), yen.end());
	return {maxX,maxY,minX,minY};
}

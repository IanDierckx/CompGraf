/*
 * Lines2D.cpp
 *
 *  Created on: Feb 20, 2018
 *      Author: uauser
 */

#include "Lines2D.h"
#include <algorithm>

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
}

void Lines2D::addLine(Line2D* line2D) {
	this->Lines.push_back(line2D);
}


Lines2D::Lines2D(vector<Line2D*> lines) {
	Lines = lines;
}

img::EasyImage Lines2D::drawLines(const int size) {

	vector<double> maxmin = this->getMinMax();
	double maxX = maxmin[0];
	double maxY = maxmin[1];
	double minX = maxmin[2];
	double minY = maxmin[3];
	double xrange = maxX-minX;
	double yrange = maxY-minY;
	double imageX = size*(xrange/max(xrange,yrange));
	double imageY = size*(yrange/max(xrange,yrange));

	img::EasyImage img = img::EasyImage(imageX, imageY);

	double schaalfactor = 0.95*(imageX/xrange);

	double DCx = schaalfactor*(minX+maxX)/2;
	double DCy = schaalfactor*(minY+maxY)/2;
	double dx = (imageX/2)-DCx;
	double dy = (imageY/2)-DCy;

	for(auto line:this->Lines) {
		double newP1X = line->p1->x*schaalfactor;
		double newP2X = line->p2->x*schaalfactor;
		double newP1Y = line->p1->y*schaalfactor;
		double newP2Y = line->p2->y*schaalfactor;

		newP1X += dx;
		newP2X += dx;
		newP1Y += dy;
		newP2Y += dy;

		unsigned int red = static_cast<unsigned int>(rint(line->color->red*255));
		unsigned int green = static_cast<unsigned int>(rint(line->color->green*255));
		unsigned int blue = static_cast<unsigned int>(rint(line->color->blue*255));

		img::Color lijnkleur = img::Color(red,green, blue);

		img.draw_line(static_cast<unsigned int>(floor(newP1X)),static_cast<unsigned int>(floor(newP1Y)),
				static_cast<unsigned int>(floor(newP2X)), static_cast<unsigned int>(floor(newP2Y)),
				lijnkleur);
	}
	return img;
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

Color::Color() {
}

Point2D::Point2D() {
}

Line2D::Line2D() {
}

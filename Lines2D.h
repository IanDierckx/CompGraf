/*
 * Lines2D.h
 *
 *  Created on: Feb 20, 2018
 *      Author: uauser
 */

#ifndef LINES2D_H_
#define LINES2D_H_

#include <iostream>
#include <vector>
#include "easy_image.h"


using namespace std;


class Color {
public:
	double red;
	double green;
	double blue;

	Color(double r, double g, double b);

	Color();
};

class Point2D {
public:
	double x;
	double y;

	Point2D(double inx, double iny);

	Point2D();
};

class Line2D {
public:
	Point2D* p1;
	Point2D* p2;
	Color* color;

	Line2D(Point2D*& point1, Point2D*& point2, Color*& col);
	Line2D();
};

class Lines2D {
private:
	vector<Line2D*> Lines;


public:
	const vector<Line2D*>& getLines() const {
		return Lines;
	}

	void addLine(Line2D*);

	Lines2D(vector<Line2D*> lines);

	vector<double> getMinMax();
};

#endif /* LINES2D_H_ */

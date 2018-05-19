/*
 * ZBuffer.h
 *
 *  Created on: May 18, 2018
 *      Author: ian
 */

#ifndef ZBUFFER_H_
#define ZBUFFER_H_

#include <vector>
#include <limits>

using namespace std;

class ZBuffer {
private:
	vector<vector<double>> buffer;
public:
	ZBuffer(const int width, const int height);

	void setBufferValue(const int x, const int y, double zValue);

	double getBufferValue(const int x, const int y);
};

#endif /* ZBUFFER_H_ */

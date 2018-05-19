/*
 * ZBuffer.cc
 *
 *  Created on: May 18, 2018
 *  Author: ian
 */

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) {
	buffer = vector<vector<double>>(height,vector<double>(width,std::numeric_limits<double>::infinity()));
}

void ZBuffer::setBufferValue(const int x, const int y, double zValue) {
	buffer[y][x] = zValue;
}


double ZBuffer::getBufferValue(const int x, const int y) {
	return buffer[y][x];
}

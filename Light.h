/*
 * Light.h
 *
 *  Created on: May 24, 2018
 *      Author: ian
 */

#ifndef LIGHT_H_
#define LIGHT_H_

#include "figure3D.h"

class Light {
public:
	Color* ambientLight;
	Color* diffuseLight;
	Vector3D ldVector;

	bool infinite;

	Light(Color* ambLight,Color* diffLight, Vector3D ld, bool inf);

	~Light();
};


class Lights3D {
private:
	vector<Light*> lights;
public:
	const vector<Light*>& getLights() const;

	void setLights(const vector<Light*>& lights);

	void addLight(Light*& lightToAdd);
};
#endif /* LIGHT_H_ */

/*
 * Light.cc
 *
 *  Created on: May 24, 2018
 *      Author: ian
 */

#include "Light.h"

Light::Light(Color* ambLight,Color* diffLight, Vector3D ld, bool inf) {
	ambientLight = ambLight;
	diffuseLight = diffLight;
	ldVector = ld;
	infinite = inf;
}

const vector<Light*>& Lights3D::getLights() const {
	return lights;
}

void Lights3D::setLights(const vector<Light*>& lights) {
	this->lights = lights;
}

void Lights3D::addLight(Light*& lightToAdd) {
	lights.push_back(lightToAdd);
}

Light::~Light() {
	delete ambientLight;
	delete diffuseLight;
}

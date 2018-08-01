#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stack>
#include <assert.h>
#include "figure3D.h"
#include "l_parser.h"
#include "lines2D.h"
#include "vector3d.h"
#include "ZBuffer.h"
#include "Light.h"
#include "easy_image.h"
#include "ini_configuration.h"

using namespace std;

inline int roundToInt(double d) {
	return static_cast<int>(round(d));
}

double clamp(double toClamp) {
	if (toClamp > 1) {
		toClamp = 1;
	} else if (toClamp < 0) {
		toClamp = 0;
	}
	return toClamp;
}

img::EasyImage generate_ColorRect(const ini::Configuration &configuration)
{
	const unsigned int width = configuration["ImageProperties"]["width"].as_int_or_die();
	const unsigned int heigth = configuration["ImageProperties"]["height"].as_int_or_die();
	img::EasyImage image(width,heigth);
	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < heigth; j++) {
			image(i,j).red = i;
			image(i,j).green = j;
			image(i,j).blue = (i+j)%256;
		}
	}
	return image;
}

img::EasyImage generateBlocks(const ini::Configuration &configuration) {
	const unsigned int width = configuration["ImageProperties"]["width"].as_int_or_die();
	const unsigned int heigth = configuration["ImageProperties"]["height"].as_int_or_die();
	img::EasyImage image(width,heigth);
	const unsigned int Nx = configuration["BlockProperties"]["nrXBlocks"].as_int_or_die();
	const unsigned int Ny = configuration["BlockProperties"]["nrYBlocks"].as_int_or_die();
	double blockWidth = width/Nx;
	double blockHeight = heigth/Ny;
	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < heigth; j++) {
			int blockX = i/blockWidth;
			int blockY = j/blockHeight;
			int som = blockX+blockY;
			std::vector<double> wit;
			std::vector<double> zwart;
			if (configuration["BlockProperties"]["invertColors"].as_bool_or_default(false)) {
				wit = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
				zwart = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
			} else {
				wit = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
				zwart = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
			}
			if (som%2 == 0) {
				image(i,j).red = roundToInt(wit[0]*255);
				image(i,j).green = roundToInt(wit[1]*255);
				image(i,j).blue = roundToInt(wit[2]*255);
			} else {
				image(i,j).red = roundToInt(zwart[0]*255);
				image(i,j).green = roundToInt(zwart[1]*255);
				image(i,j).blue = roundToInt(zwart[2]*255);
			}
		}
	}
	return image;
}

img::EasyImage QuarterCircle(const ini::Configuration &configuration) {
	const unsigned int width = configuration["ImageProperties"]["width"].as_int_or_die();
	const unsigned int heigth = configuration["ImageProperties"]["height"].as_int_or_die();
	const unsigned int aantalLijnen = configuration["LineProperties"]["nrLines"].as_int_or_die();
	std::vector<double> kleurLijn = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();
	std::vector<double> backgroundcolor = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
	img::EasyImage image(width,heigth);
	img::Color lijnKLeur;
	lijnKLeur.red = kleurLijn[0]*255;
	lijnKLeur.green = kleurLijn[1]*255;
	lijnKLeur.blue = kleurLijn[2]*255;
	img::Color achtergrond;
	achtergrond.red = backgroundcolor[0]*255;
	achtergrond.green = backgroundcolor[1]*255;
	achtergrond.blue = backgroundcolor[2]*255;
	for(unsigned int i = 0; i < width; i++) {
		for(unsigned int j = 0; j < heigth; j++) {
			image(i,j).red = achtergrond.red;
			image(i,j).green = achtergrond.green;
			image(i,j).blue = achtergrond.blue;
		}
	}
	const unsigned int verticaleAfstand = heigth/(aantalLijnen-1);
	const unsigned int horizontaleAfstand = width/(aantalLijnen-1);
	for(unsigned int i = 0; i < aantalLijnen; i++) {
		image.draw_line(horizontaleAfstand*i,heigth-1,0,verticaleAfstand*i,lijnKLeur);
	}
	return image;
}

void draw_zbuff_line(ZBuffer*& buffer, img::EasyImage& img, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, const img::Color color) {
	assert(x0 < img.get_width() && y0 < img.get_height());
	assert(x1 < img.get_width() && y1 < img.get_height());
	if (x0 == x1)
	{
		//special case for x0 == x1
		double a = std::max(y0, y1)-std::min(y0, y1);
		int count = a;
		for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++)
		{
			unsigned int xi = x0;
			unsigned int yi = i;
			double zValue = ((count/a)/z0)+((1-(count/a))/z1);
			if (buffer->getBufferValue(x0,i) > zValue) {
				(img)(xi, yi) = color;

			}
			count--;
		}
	}
	else if (y0 == y1)
	{
		//special case for y0 == y1
		double a = std::max(x0, x1)-std::min(x0, x1);
		int count = a;
		for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++)
		{
			unsigned int xi = i;
			unsigned int yi = y0;
			double zValue = ((count/a)/z0)+((1-(count/a))/z1);
			if (buffer->getBufferValue(x0,i) > zValue) {
				(img)(xi, yi) = color;
				buffer->setBufferValue(i, y0,zValue);
			}
			count--;
		}
	}
	else
	{
		if (x0 > x1)
		{
			//flip points if x1>x0: we want x0 to have the lowest value
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
		if (-1.0 <= m && m <= 1.0)
		{
			double a = (x1 - x0);
			double count = a;
			for (unsigned int i = 0; i <= (x1 - x0); i++)
			{
				unsigned int xi = x0 + i;
				unsigned int yi = (unsigned int) round(y0 + m * i);
				double zValue = ((count/a)/z0)+((1-(count/a))/z1);
				if (buffer->getBufferValue(xi, yi) > zValue) {
					(img)(xi, yi) = color;
					buffer->setBufferValue(xi, yi,zValue);
				}
				count--;
			}
		}
		else if (m > 1.0)
		{
			double a = (y0 - y1);
			int count = a;
			for (unsigned int i = 0; i <= (y1 - y0); i++)
			{
				unsigned int xi = (unsigned int) round(x0 + (i / m));
				unsigned int yi = y0 + i;
				double zValue = ((count/a)/z0)+((1-(count/a))/z1);
				if (buffer->getBufferValue(xi, yi) > zValue) {
					(img)(xi, yi) = color;
					buffer->setBufferValue(xi, yi,zValue);
				}
				count--;
			}
		}
		else if (m < -1.0)
		{
			double a = (y0 - y1);
			int count = a;
			for (unsigned int i = 0; i <= (y0 - y1); i++)
			{
				unsigned int xi = (unsigned int) round(x0 - (i / m));
				unsigned int yi = y0 - i;
				double zValue = ((count/a)/z0)+((1-(count/a))/z1);
				if (buffer->getBufferValue(xi, yi) > zValue) {
					(img)(xi, yi) = color;
					buffer->setBufferValue(xi, yi,zValue);
				}
				count--;
			}
		}
	}
}

void drawZBuffTriangle(ZBuffer*& buffer,img::EasyImage& img, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d, double dx, double dy, Color*& ambientReflection, Color*& diffuseReflection, Lights3D& lights, Matrix eyePointMatrix) {
	Point2D* projectedA = new Point2D((d*A.x)/(-A.z)+dx,(d*A.y)/(-A.z)+dy);
	Point2D* projectedB = new Point2D((d*B.x)/(-B.z)+dx,(d*B.y)/(-B.z)+dy);
	Point2D* projectedC = new Point2D((d*C.x)/(-C.z)+dx,(d*C.y)/(-C.z)+dy);

	int yMin = roundToInt(min(min(projectedA->y,projectedB->y),projectedC->y)+0.5);
	int yMax = roundToInt(max(max(projectedA->y,projectedB->y),projectedC->y)-0.5);

	vector<vector<Point2D*>> pointsInTriangle;

	double xG = (projectedA->x+projectedB->x+projectedC->x)/3;
	double yG = (projectedA->y+projectedB->y+projectedC->y)/3;

	double inverseZG = (1/(3*A.z))+(1/(3*B.z))+(1/(3*C.z));

	Vector3D u = B-A;
	Vector3D v = C-A;

	Vector3D w = Vector3D::cross(u,v);

	double k = (w.x*A.x)+(w.y*A.y)+(w.z*A.z);

	double dzdx = w.x/(-d*k);
	double dzdy = w.y/(-d*k);

	Vector3D normal = Vector3D::normalise(w);


	Color tempColor = Color(0,0,0);

	vector<Light*> nonInfiniteLights;

	for (auto& light:lights.getLights()) {
		if (light->infinite) {
			Vector3D l = -(light->ldVector/light->ldVector.length());
			double cosAlph = Vector3D::dot(l,normal);
			if (cosAlph >= 0) {
				tempColor.red += (light->diffuseLight->red*diffuseReflection->red)*cosAlph;
				tempColor.green += (light->diffuseLight->green*diffuseReflection->green)*cosAlph;
				tempColor.blue += (light->diffuseLight->blue*diffuseReflection->blue)*cosAlph;
			}
			tempColor.red += light->ambientLight->red*ambientReflection->red;
			tempColor.green += light->ambientLight->green*ambientReflection->green;
			tempColor.blue += light->ambientLight->blue*ambientReflection->blue;
		} else {
			nonInfiniteLights.push_back(light);
		}

	}




	for (int yValue = yMin; yValue <= yMax; ++yValue) {
		double xLAB = std::numeric_limits<double>::infinity();
		double xLAC = xLAB;
		double xLBC = xLAB;
		double xRAB = -xLAB;
		double xRAC = xRAB;
		double xRBC = xRAB;

		if (((yValue-projectedA->y)*(yValue-projectedB->y)<=0) and (projectedA->y != projectedB->y)) {
			double xi = projectedB->x+(projectedA->x-projectedB->x)*((yValue-projectedB->y)/(projectedA->y-projectedB->y));
			xLAB = xi;
			xRAB = xi;
		}
		if (((yValue-projectedA->y)*(yValue-projectedC->y)<=0) and (projectedA->y != projectedC->y)) {
			double xi = projectedC->x+(projectedA->x-projectedC->x)*((yValue-projectedC->y)/(projectedA->y-projectedC->y));
			xLAC = xi;
			xRAC = xi;
		}
		if (((yValue-projectedB->y)*(yValue-projectedC->y)<=0) and (projectedB->y != projectedC->y)) {
			double xi = projectedC->x+(projectedB->x-projectedC->x)*((yValue-projectedC->y)/(projectedB->y-projectedC->y));
			xLBC = xi;
			xRBC = xi;
		}
		int xL = roundToInt(min(min(xLAB,xLAC),xLBC)+0.5);
		int xR = roundToInt(max(max(xRAB,xRAC),xRBC)-0.5);

		//special case for y0 == y1
		for (unsigned int i = std::min(xL, xR); i <= std::max(xL, xR); i++)
		{
			double zValue = 1.0001*inverseZG+(i-xG)*dzdx+(yValue-yG)*dzdy;
			for (auto& light:nonInfiniteLights) {
				double originalZ = 1/zValue;
				double originalX = -(i/d)*originalZ;
				double originalY = -(yValue/d)*originalZ;
				Vector3D originalPoint = Vector3D::point(originalX,originalY,originalZ);
				Vector3D l = -light->ldVector-originalPoint;
				l.normalise();
				double cosAlph = Vector3D::dot(l,normal);
				if (cosAlph >= 0) {
					tempColor.red += (light->diffuseLight->red*diffuseReflection->red)*cosAlph;
					tempColor.green += (light->diffuseLight->green*diffuseReflection->green)*cosAlph;
					tempColor.blue += (light->diffuseLight->blue*diffuseReflection->blue)*cosAlph;
				}
				tempColor.red += light->ambientLight->red*ambientReflection->red;
				tempColor.green += light->ambientLight->green*ambientReflection->green;
				tempColor.blue += light->ambientLight->blue*ambientReflection->blue;
			}
			img::Color lijnkleur = img::Color(clamp(tempColor.red)*255,clamp(tempColor.green)*255,clamp(tempColor.blue)*255);
			if (buffer->getBufferValue(i,yValue) >= zValue) {
				(img)(i, yValue) = lijnkleur;
				buffer->setBufferValue(i, yValue,zValue);
			}
		}
	}
}

img::EasyImage draw2DLines(const int size, Lines2D& lines, img::Color backgroundColor, bool ZBufferOn) {
	vector<double> maxmin = lines.getMinMax();
	double maxX = maxmin[0];
	double maxY = maxmin[1];
	double minX = maxmin[2];
	double minY = maxmin[3];
	if (maxX == 0) {
		maxX+=0.001;
	}
	if (maxY == 0) {
		maxY+=0.001;
	}
	double xrange = maxX-minX;
	double yrange = maxY-minY;
	double imageX = size*(xrange/max(xrange,yrange));
	double imageY = size*(yrange/max(xrange,yrange));

	img::EasyImage img = img::EasyImage(ceil(imageX), ceil(imageY), backgroundColor);

	double schaalfactor = 0.95*(imageX/xrange);

	double DCx = schaalfactor*(minX+maxX)/2;
	double DCy = schaalfactor*(minY+maxY)/2;
	double dx = (imageX/2)-DCx;
	double dy = (imageY/2)-DCy;


	ZBuffer* buffer = new ZBuffer(img.get_width(),img.get_height());

	for(auto line:lines.getLines()) {
		double newP1X = line->p1->x*schaalfactor;
		double newP2X = line->p2->x*schaalfactor;
		double newP1Y = line->p1->y*schaalfactor;
		double newP2Y = line->p2->y*schaalfactor;

		newP1X += dx;
		newP2X += dx;
		newP1Y += dy;
		newP2Y += dy;



		img::Color lijnkleur = img::Color(line->color->red,line->color->green, line->color->blue);
		if (ZBufferOn) {
			draw_zbuff_line(buffer,img,roundToInt(newP1X),roundToInt(newP1Y),line->z1,roundToInt(newP2X),roundToInt(newP2Y),line->z2,lijnkleur);
		} else {
			img.draw_line(roundToInt(newP1X),roundToInt(newP1Y),
							roundToInt(newP2X), roundToInt(newP2Y),
							lijnkleur);
		}
	}
	return img;
}



string replaceRule(string currentRule, unsigned int currentIteration, LParser::LSystem2D& lSystem) {
	string replacedRule = "";
	currentIteration += 1;
	for (char current:currentRule) {
		if ((current == '+') or (current == '-') or (current == '(') or (current == ')')) {
			replacedRule += current;
		} else if (find(lSystem.get_alphabet().begin(),lSystem.get_alphabet().end(),current) != lSystem.get_alphabet().end()) {
			replacedRule += lSystem.get_replacement(current);
		}
	}
	if (currentIteration != lSystem.get_nr_iterations()) {
		replacedRule = replaceRule(replacedRule,currentIteration,lSystem);
	}
	return replacedRule;
}

vector<Line2D*> parse_rule(LParser::LSystem2D& lSystem, vector<double> currentPoint, Color* lijnkleur, double current_angle, double angle_change) {
	vector<Line2D*> lines;
	string startRule = lSystem.get_initiator();
	string rule = "";
	if (lSystem.get_nr_iterations() == 0) {
		rule = startRule;
	} else {
		rule = replaceRule(startRule,0,lSystem);
	}
	double currentX = currentPoint[0];
	double currentY = currentPoint[1];
	double angleRightNow = current_angle;
	stack<vector<double>, vector<vector<double>>> bracketPositions;
	for (auto current:rule) {
		if (current == '+') {
			angleRightNow += angle_change;
		} else if (current == '-') {
			angleRightNow -= angle_change;
		} else if (current == '(') {
			vector<double> currentPos {currentX, currentY, angleRightNow};
			bracketPositions.push(currentPos);
		} else if (current == ')') {
			if (!bracketPositions.empty()) {
				vector<double> currentPos = bracketPositions.top();
				bracketPositions.pop();
				currentX = currentPos[0];
				currentY = currentPos[1];
				angleRightNow = currentPos[2];
			}
		} else if (lSystem.draw(current)) {
			Point2D* p1 = new Point2D(currentX, currentY);
			currentX += cos(angleRightNow);
			currentY += sin(angleRightNow);
			Point2D* p2 = new Point2D(currentX, currentY);
			lines.push_back(new Line2D(p1,p2,lijnkleur));
		}
		else {
			currentX += cos(angleRightNow);
			currentY += sin(angleRightNow);
		}
	}
	return lines;
}

img::EasyImage generate2DLSys(const ini::Configuration &configuration) {
	int size = configuration["General"]["size"].as_int_or_die();
	vector<double> background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
	string inputfile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
	vector<double> lijnkleur = configuration["2DLSystem"]["color"].as_double_tuple_or_die();

	img::Color backgroundColor = img::Color(background[0]*255,background[1]*255,background[2]*255);

	LParser::LSystem2D lSystem;
	ifstream input_stream(inputfile);
	input_stream >> lSystem;
	input_stream.close();

	double starting_angle = lSystem.get_starting_angle()*(M_PI/180);
	double angle = lSystem.get_angle()*(M_PI/180);

	Color* lijnkleurColor = new Color(lijnkleur[0]*255, lijnkleur[1]*255, lijnkleur[2]*255);

	vector<Line2D*> lines = parse_rule(lSystem,vector<double>(2,0),lijnkleurColor,starting_angle,angle);

	Lines2D lines2D = Lines2D(lines);

	img::EasyImage img = draw2DLines(size,lines2D,backgroundColor,false);

	return img;
}

img::EasyImage generateZBuffTriangulatedImage(Lines2D& projectie,int size, img::Color backgroundColor, Figures3D& figuren, Lights3D& lights, Matrix eyePointMatrix) {
	vector<double> maxmin = projectie.getMinMax();
	double maxX = maxmin[0];
	double maxY = maxmin[1];
	double minX = maxmin[2];
	double minY = maxmin[3];
	double xrange = maxX-minX;
	double yrange = maxY-minY;
	double imageX = size*(xrange/max(xrange,yrange));
	double imageY = size*(yrange/max(xrange,yrange));

	img::EasyImage img = img::EasyImage(ceil(imageX), ceil(imageY), backgroundColor);

	double schaalfactor = 0.95*(imageX/xrange);

	double DCx = schaalfactor*(minX+maxX)/2;
	double DCy = schaalfactor*(minY+maxY)/2;
	double dx = (imageX/2)-DCx;
	double dy = (imageY/2)-DCy;

	ZBuffer* buffer = new ZBuffer(img.get_width(),img.get_height());
	for (auto figuur:figuren.getFigures()) {
		vector<Face*> newFaces;
		for (auto face:figuur->faces) {
			auto faceTriangles = figuur->triangulateFace(face);
			for (auto faceInOriginalFace:faceTriangles) {
				newFaces.push_back(faceInOriginalFace);
			}
		}
		for (auto newFace:newFaces) {
			drawZBuffTriangle(buffer,img,figuur->points[newFace->getPointIndexes()[0]],figuur->points[newFace->getPointIndexes()[1]],
					figuur->points[newFace->getPointIndexes()[2]],schaalfactor,dx,dy,figuur->ambientReflection, figuur->diffusereflection, lights, eyePointMatrix);
		}
	}
	return img;
}

void generateFractal(Figure3D figuur, Figures3D& fractal,const int nr_iterations, const double scale) {
	for (int times = 1; times <= nr_iterations; ++times) {
		vector<Figure3D*> newFigures;
		vector<Figure3D*> oldFigures;
		for (auto& figureInFractal:fractal.getFigures()) {
			for (int pointIndex = 0; pointIndex < figureInFractal->points.size(); ++pointIndex) {
				Figure3D* copyFiguur = new Figure3D(figuur);
				double scalefactor = 1/pow(scale,times);
				copyFiguur->scaleFigure(scalefactor);
				copyFiguur->applyTransformations();
				Vector3D translateVector = figureInFractal->points[pointIndex]-copyFiguur->points[pointIndex];
				copyFiguur->translate(translateVector);
				copyFiguur->applyTransformations();
				newFigures.push_back(copyFiguur);
			}
			oldFigures.push_back(figureInFractal);
		}
		for (auto figure:newFigures) {
			fractal.addFigure(figure);
		}
		for (auto oldFigure:oldFigures) {
			fractal.deleteFigure(oldFigure);
		}
	}

}

img::EasyImage generate3DLines(const ini::Configuration &configuration) {
	int size = configuration["General"]["size"].as_int_or_die();

	vector<double> background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
	img::Color backgroundColor = img::Color(background[0]*255,background[1]*255,background[2]*255);

	unsigned int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
	vector<double> eye = configuration["General"]["eye"].as_double_tuple_or_die();
	Vector3D eyePoint = Vector3D::point(eye[0],eye[1],eye[2]);

	Figures3D figuren;
	unsigned int currentFigure = 0;
	while (currentFigure < nrFigures) {
		string currentFigureString = "Figure" + to_string(currentFigure);
		Color* lijnkleur;
		if (configuration["General"]["type"].as_string_or_die() == "LightedZBuffering") {
			vector<double> ambientReflVector = configuration[currentFigureString]["ambientReflection"].as_double_tuple_or_die();
			lijnkleur = new Color(ambientReflVector[0], ambientReflVector[1], ambientReflVector[2]);
		} else {
			vector<double> lijnkleurVector = configuration[currentFigureString]["color"].as_double_tuple_or_die();
			lijnkleur = new Color(lijnkleurVector[0]*255, lijnkleurVector[1]*255, lijnkleurVector[2]*255);
		}
		Figure3D* figuur;
		string type = configuration[currentFigureString]["type"].as_string_or_die();
		if (type == "LineDrawing") {
			vector<Vector3D> punten;
			unsigned int nrPunten = configuration[currentFigureString]["nrPoints"].as_int_or_die();
			unsigned int huidigPunt = 0;
			while (huidigPunt < nrPunten) {
				vector<double> puntVector = configuration[currentFigureString]["point"+to_string(huidigPunt)].as_double_tuple_or_die();
				Vector3D punt =  Vector3D::point(puntVector[0], puntVector[1], puntVector[2]);
				punten.push_back(punt);
				huidigPunt += 1;
			}
			vector<Face*> faces;
			unsigned int nrLijnen = configuration[currentFigureString]["nrLines"].as_int_or_die();
			unsigned int huidigeLijn = 0;
			while (huidigeLijn < nrLijnen) {
				vector<int> faceVector = configuration[currentFigureString]["line"+to_string(huidigeLijn)].as_int_tuple_or_die();
				Face* face = new Face(faceVector);
				faces.push_back(face);
				huidigeLijn += 1;
			}
			figuur = new Figure3D(punten,faces,lijnkleur);
		} else if (type == "Tetrahedron" or type == "FractalTetrahedron") {
			figuur = figuren.createTetrahedron(lijnkleur);
		} else if (type == "Cube" or type == "FractalCube") {
			figuur = figuren.createCube(lijnkleur);
		} else if (type == "Octahedron" or type == "FractalOctahedron") {
			figuur = figuren.createOctahedron(lijnkleur);
		} else if (type == "Icosahedron" or type == "FractalIcosahedron") {
			figuur = figuren.createIcosahedron(lijnkleur);
		} else if (type == "Dodecahedron" or type == "FractalDodecahedron") {
			figuur = figuren.createDodecahedron(lijnkleur);
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Sphere") {
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			figuur = figuren.createSphere(n,lijnkleur);
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Cone") {
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			double height = configuration[currentFigureString]["height"].as_double_or_die();
			figuur = figuren.createCone(n,height,lijnkleur);
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Cylinder") {
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			double height = configuration[currentFigureString]["height"].as_double_or_die();
			figuur = figuren.createCylinder(n,height,lijnkleur);
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Torus") {
			double R = configuration[currentFigureString]["R"].as_double_or_die();
			double r = configuration[currentFigureString]["r"].as_double_or_die();
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			int m = configuration[currentFigureString]["m"].as_int_or_die();
			figuur = figuren.createTorus(R,r,n,m,lijnkleur);
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "3DLSystem") {
			string inputFile = configuration[currentFigureString]["inputfile"].as_string_or_die();
			figuur = figuren.create3DLSystem(inputFile,lijnkleur);
		}
		if (configuration["General"]["type"].as_string_or_die() == "LightedZBuffering") {
			vector<double> diffRefl;
			diffRefl = configuration[currentFigureString]["diffuseReflection"].as_double_tuple_or_default(vector<double>(3,0));
			figuur->diffusereflection = new Color(diffRefl[0],diffRefl[1],diffRefl[2]);
		}
		currentFigure += 1;
		if (type == "FractalCube" or type == "FractalDodecahedron" or type == "FractalIcosahedron" or type == "FractalOctahedron" or type == "FractalTetrahedron") {
			Figures3D fractal;
			fractal.addFigure(figuur);
			double scale = configuration[currentFigureString]["fractalScale"].as_double_or_die();
			int nr_iterations = configuration[currentFigureString]["nrIterations"].as_int_or_die();
			generateFractal(*figuur,fractal,nr_iterations,scale);
			for (auto figure:fractal.getFigures()) {
				figure->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
				figure->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
				figure->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
				figure->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
				vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
				Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
				figure->translate(center);
				figuren.addFigure(figure);
			}
		} else {
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
		}
	}
	Matrix eyePointMatrix = figuren.eyepointTrans(eyePoint);
	figuren.applyTransformations(eyePointMatrix);
	Lines2D projectie = figuren.doProjection();
	Lights3D lights;
	string typeString = configuration["General"]["type"].as_string_or_die();
	if (typeString == "LightedZBuffering") {
		unsigned int nrLights = configuration["General"]["nrLights"].as_int_or_die();
		unsigned int currentLight = 0;
		while (currentLight < nrLights) {
			string currentLightString = "Light" + to_string(currentLight);
			vector<double> ambientLightVector = configuration[currentLightString]["ambientLight"].as_double_tuple_or_die();
			Color* ambientLight = new Color(ambientLightVector[0],ambientLightVector[1],ambientLightVector[2]);
			bool inf = configuration[currentLightString]["infinity"].as_bool_or_default(true);
			vector<double> directionVector = configuration[currentLightString]["direction"].as_double_tuple_or_default(vector<double>(3,0));
			vector<double> diffuseLightVector = configuration[currentLightString]["diffuseLight"].as_double_tuple_or_default(vector<double>(3,0));
			Color* diffuseLight = new Color(diffuseLightVector[0],diffuseLightVector[1],diffuseLightVector[2]);
			Vector3D direction = Vector3D::vector(directionVector[0],directionVector[1],directionVector[2]);
			direction *= eyePointMatrix;
			vector<double> positionVector = configuration[currentLightString]["position"].as_double_tuple_or_default(vector<double>(3,0));
			Vector3D position = Vector3D::point(positionVector[0],positionVector[1],positionVector[2]);
			position *= eyePointMatrix;
			Light* newLight;
			if (inf) {
				newLight = new Light(ambientLight,diffuseLight,direction,inf);
			} else {
				newLight = new Light(ambientLight,diffuseLight,position,inf);
			}
			lights.addLight(newLight);
			currentLight++;
		}
	}
	if (typeString == "Wireframe") {
		return draw2DLines(size,projectie,backgroundColor,false);
	} else if (typeString == "ZBufferedWireframe") {
		return draw2DLines(size,projectie,backgroundColor,true);
	} else if (typeString == "ZBuffering" or typeString == "LightedZBuffering") {
		return generateZBuffTriangulatedImage(projectie,size,backgroundColor,figuren,lights,eyePointMatrix);
	}
	img::EasyImage img = img::EasyImage(10,10);
	return img;
}

img::EasyImage generate_image(const ini::Configuration &configuration)
{
	const std::string typeString = configuration["General"]["type"].as_string_or_die();
	if (typeString == "IntroColorRectangle") {
		return generate_ColorRect(configuration);
	} else if (typeString == "IntroBlocks") {
		return generateBlocks(configuration);
	} else if (typeString == "IntroLines") {
		const std::string figure = configuration["LineProperties"]["figure"].as_string_or_die();
		if (figure=="QuarterCircle") {
			return QuarterCircle(configuration);
		}
	} else if (typeString == "2DLSystem") {
		return generate2DLSys(configuration);
	} else if (typeString == "Wireframe" or typeString == "ZBufferedWireframe" or typeString == "ZBuffering" or typeString == "LightedZBuffering") {
		return generate3DLines(configuration);
	}
	img::EasyImage img = img::EasyImage(10,10);
	return img;
}

int main(int argc, char const* argv[])
{
		int retVal = 0;
        try
        {
                for(int i = 1; i < argc; ++i)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(argv[i]);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string fileName(argv[i]);
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << argv[i] << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}

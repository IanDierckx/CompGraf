#include "easy_image.h"
#include "ini_configuration.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stack>
#include "figure3D.h"
#include "l_parser.h"
#include "lines2D.h"
#include "vector3d.h"

using namespace std;

inline int roundToInt(double d) {
	return static_cast<int>(round(d));
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

img::EasyImage draw2DLines(const int size, Lines2D& lines, img::Color backgroundColor) {

	vector<double> maxmin = lines.getMinMax();
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

		img.draw_line(roundToInt(newP1X),roundToInt(newP1Y),
				roundToInt(newP2X), roundToInt(newP2Y),
				lijnkleur);
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
	string rule = replaceRule(startRule,0,lSystem);
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

	img::EasyImage img = draw2DLines(size,lines2D,backgroundColor);

	return img;
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
		vector<double> lijnkleurVector = configuration[currentFigureString]["color"].as_double_tuple_or_die();
		Color* lijnkleur = new Color(lijnkleurVector[0]*255, lijnkleurVector[1]*255, lijnkleurVector[2]*255);
		if (configuration[currentFigureString]["type"].as_string_or_die() == "LineDrawing") {
			vector<Vector3D> punten;
			unsigned int nrPunten = configuration[currentFigureString]["nrPoints"].as_int_or_die();
			unsigned int huidigPunt = 0;
			while (huidigPunt < nrPunten) {
				vector<double> puntVector = configuration[currentFigureString]["point"+(huidigPunt)].as_double_tuple_or_die();
				Vector3D punt =  Vector3D::point(puntVector[0], puntVector[1], puntVector[2]);
				punten.push_back(punt);
				huidigPunt += 1;
			}
			vector<Face*> faces;
			unsigned int nrLijnen = configuration[currentFigureString]["nrLines"].as_int_or_die();
			unsigned int huidigeLijn = 0;
			while (huidigeLijn < nrLijnen) {
				vector<int> faceVector = configuration[currentFigureString]["line"+(huidigeLijn)].as_int_tuple_or_die();
				Face* face = new Face(faceVector);
				faces.push_back(face);
				huidigeLijn += 1;
			}
			Figure3D* figuur = new Figure3D(punten,faces,lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Tetrahedron") {
			Figure3D* figuur = figuren.createTetrahedron(lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Cube") {
			Figure3D* figuur = figuren.createCube(lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Octahedron") {
			Figure3D* figuur = figuren.createOctahedron(lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Icosahedron") {
			Figure3D* figuur = figuren.createIcosahedron(lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Dodecahedron") {
			Figure3D* figuur = figuren.createDodecahedron(lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Sphere") {
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			Figure3D* figuur = figuren.createSphere(n,lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Cone") {
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			double height = configuration[currentFigureString]["height"].as_double_or_die();
			Figure3D* figuur = figuren.createCone(n,height,lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Cylinder") {
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			double height = configuration[currentFigureString]["height"].as_double_or_die();
			Figure3D* figuur = figuren.createCylinder(n,height,lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "Torus") {
			double R = configuration[currentFigureString]["R"].as_double_or_die();
			double r = configuration[currentFigureString]["r"].as_double_or_die();
			int n = configuration[currentFigureString]["n"].as_int_or_die();
			int m = configuration[currentFigureString]["m"].as_int_or_die();
			Figure3D* figuur = figuren.createTorus(R,r,n,m,lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		} else if (configuration[currentFigureString]["type"].as_string_or_die() == "3DLSystem") {
			string inputFile = configuration[currentFigureString]["inputfile"].as_string_or_die();
			Figure3D* figuur = figuren.create3DLSystem(inputFile,lijnkleur);
			figuur->scaleFigure(configuration[currentFigureString]["scale"].as_double_or_die());
			figuur->rotateX(configuration[currentFigureString]["rotateX"].as_double_or_die());
			figuur->rotateY(configuration[currentFigureString]["rotateY"].as_double_or_die());
			figuur->rotateZ(configuration[currentFigureString]["rotateZ"].as_double_or_die());
			vector<double> centerVector = configuration[currentFigureString]["center"].as_double_tuple_or_die();
			Vector3D center = Vector3D::vector(centerVector[0], centerVector[1], centerVector[2]);
			figuur->translate(center);
			figuren.addFigure(figuur);
			currentFigure += 1;
		}
	}

	Matrix eyePointMatrix = figuren.eyepointTrans(eyePoint);
	figuren.applyTransformations(eyePointMatrix);
	Lines2D projectie = figuren.doProjection();

	return draw2DLines(size,projectie,backgroundColor);
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
	} else if (typeString == "Wireframe") {
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

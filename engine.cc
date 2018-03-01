#include "easy_image.h"
#include "ini_configuration.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include "Lines2D.cpp"

inline int roundToInt(double d) {
	return static_cast<int>(std::round(d));
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
img::EasyImage EyeDraw(const ini::Configuration &configuration) {
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
		} else if (figure=="Eye") {
			return EyeDraw(configuration);
		}
	} else if (typeString == "test") {
		Color* col = new Color(0.5,0.1,0);
		Point2D* p1 = new Point2D(1.5,1);
		Point2D* p2 = new Point2D(5.8,1);
		Point2D* p3 = new Point2D(5.3,1);
		Point2D* p4 = new Point2D(5.3,5);
		Line2D* line1 = new Line2D(p1,p2,col);
		Line2D* line2 = new Line2D(p3,p4,col);
		vector<Line2D*> linesvector = {line1, line2};
		Lines2D lines = Lines2D(linesvector);
		return lines.drawLines(500);
	}
	img::EasyImage img;
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

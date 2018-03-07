<<<<<<< HEAD
engine.o: engine.cc easy_image.h ini_configuration.h Lines2D.cpp \
 Lines2D.h l_parser.h vector3d.h Figure3D.h
=======
engine.o: engine.cc easy_image.h ini_configuration.h Lines2D.h l_parser.h
>>>>>>> ca9edabcd5edfedcef798d9221314bd74d973b3b
	$(CC) $(CXXFLAGS) -c $< -o $@

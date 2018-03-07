engine.o: engine.cc easy_image.h ini_configuration.h Lines2D.cpp \
 Lines2D.h l_parser.h vector3d.h Figure3D.h
	$(CC) $(CXXFLAGS) -c $< -o $@

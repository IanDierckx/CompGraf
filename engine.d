engine.o: engine.cc easy_image.h ini_configuration.h figure3D.h \
 vector3d.h lines2D.h l_parser.h
	$(CC) $(CXXFLAGS) -c $< -o $@

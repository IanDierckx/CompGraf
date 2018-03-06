engine.o: engine.cc easy_image.h ini_configuration.h Lines2D.h l_parser.h
	$(CC) $(CXXFLAGS) -c $< -o $@

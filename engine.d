engine.o: engine.cc easy_image.h ini_configuration.h Lines2D.cpp \
 Lines2D.h
	$(CC) $(CXXFLAGS) -c $< -o $@

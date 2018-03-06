Lines2D.o: Lines2D.cc Lines2D.h easy_image.h
	$(CC) $(CXXFLAGS) -c $< -o $@

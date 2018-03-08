lines2D.o: lines2D.cc lines2D.h easy_image.h
	$(CC) $(CXXFLAGS) -c $< -o $@

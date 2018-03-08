figure3D.o: figure3D.cc figure3D.h vector3d.h lines2D.h easy_image.h
	$(CC) $(CXXFLAGS) -c $< -o $@

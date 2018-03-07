vector3d.o: vector3d.cc vector3d.h
	$(CC) $(CXXFLAGS) -c $< -o $@

ini_configuration.o: ini_configuration.cc ini_configuration.h
	$(CC) $(CXXFLAGS) -c $< -o $@

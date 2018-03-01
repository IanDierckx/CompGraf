CXXFLAGS   = -Wall -Wextra -g -fstack-protector-all -std=c++11
LDFLAGS    =
EXECUTABLE = engine
EXTENSION  = cc
SOURCES    = $(basename $(shell find . -name '*.$(EXTENSION)'))

.PHONY: all
all: $(EXECUTABLE)

$(EXECUTABLE): $(addsuffix .o,$(SOURCES))
	$(CXX) $(LDFLAGS) $^ -o $@

%.d: %.$(EXTENSION)
	$(CXX) $(CXXFLAGS) -MM $< -o $@
	echo sed -i 's/:/ $@:/' $@
	printf '\t$$(CC) $$(CXXFLAGS) -c $$< -o $$@\n' >>$@

-include $(addsuffix .d,$(SOURCES))

.PHONY: clean
clean:
	find . -name '*.o' -delete
	find . -name '*.d' -delete
	find . -name '*~'  -delete

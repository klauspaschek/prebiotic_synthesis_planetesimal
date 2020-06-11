CXXFLAGS = -std=c++11 -Wall -shared -lgfortran -O3

PYBIND11FLAGS = -fPIC `python3 -m pybind11 --includes`
PYBIND11SUFFIX = `python3-config --extension-suffix`

LIBSDIR = ./chemapp-v740/
LIBS = $(LIBSDIR)libChemAppCS.so $(LIBSDIR)libEChemAppS.so


TARGET	= ChemApp

$(TARGET)$(PYBIND11SUFFIX): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(PYBIND11FLAGS) $(TARGET).cpp -o $(TARGET)$(PYBIND11SUFFIX) $(LIBS)

# In the rule below "|more" is a workaround to overcome a
# problem noticed on some Linux systems. Without it,
# output from the ChemApp library (tqshow, tqcel, etc.)
# and output via the C runtime library (printf etc.)
# are not in the proper sequence when the output is redirected
# to a file using ">", despite the fact that fflush(NULL)
# is used in the code.
#test: ./$(TARGET)
#	./$(TARGET) |more > $(TARGET).rst

clean:
	-rm $(TARGET)$(PYBIND11SUFFIX)

CXXFLAGS= -Wall -std=c++11

LD= gfortran
EXTRACFLAGS= -DUNIX
EXTRACLIBS= -lf2c -lm -lgfortran
EXTRALDFLAGS =
EXTRALDLIBS = -lstdc++

LIBS = ../calight-v740/libChemAppC.a ../calight-v740/libLChemApp.a

# This examples Makefile uses the STATIC version of ChemApp. If your
# distribution also contains the SHARED version of the ChemApp
# libraries, you can also link them dynamically. In this case,
# uncomment the "LIBS"-line below instead.
# Make sure though that the environment variable LD_LIBRARY_PATH
# contains the directory where you choose to put the ChemApp shared
# libs, otherwise your application program won't be able to find them!
# For instance, for you to be able to run the example program cademo1
# from this directory, LD_LIBRARY_PATH needs to contain the path ".."!

#LIBS = ../libChemAppCS.so ../libLChemAppS.so


TARGET	= nucleobaseSynthesis

default: $(TARGET)

all: $(TARGET) test

$(TARGET).o: $(TARGET).cpp
	$(CXX) $(CXXFLAGS) -I. $(CFLAGS) $(EXTRACFLAGS) -c $(TARGET).cpp

$(TARGET): $(TARGET).o
	$(LD) $(EXTRALDFLAGS) -o $(TARGET) $(TARGET).o $(LIBS) $(EXTRALDLIBS)

# In the rule below "|more" is a workaround to overcome a
# problem noticed on some Linux systems. Without it,
# output from the ChemApp library (tqshow, tqcel, etc.)
# and output via the C runtime library (printf etc.)
# are not in the proper sequence when the output is redirected
# to a file using ">", despite the fact that fflush(NULL)
# is used in the code.
test: ./$(TARGET)
	./$(TARGET) |more > $(TARGET).rst

clean:
	-rm $(TARGET).o




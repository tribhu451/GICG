GSLLIBS      := $(shell gsl-config --libs)
EXTRA_FLAGS  = -D TABLE 
CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     +=  $(EXTRA_FLAGS)
LIBS          =  $(SYSLIBS) $(GSLLIBS)

vpath %.cpp src
objdir     = obj

SRC        = main.cpp job.cpp grid.cpp mc_glau.cpp cell.cpp random.cpp mc_glau_smear.cpp
             
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	 = mc_glauber
#-------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

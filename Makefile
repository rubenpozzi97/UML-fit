#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore -lMinuit
ROOTCINT  := $(shell which rootcint)

#directories
SOURCEDIR   := ./src
INCLUDEDIR  := ./interface

#exe_files
EXECUTABLE0 := simfit_recoMC_singleComponent
EXECUTABLE1 := simfit4d_recoMC_singleComponent
EXECUTABLE2 := simfit_recoMC_fullAngular

EXTRACLASS := RooDataHist.cxx
CLASS0     := PdfRT
CLASS1     := PdfWT
CLASS2     := DecayRate
CLASS3     := PdfSigAng
CLASS4     := RooDoubleCBFast
CLASS5     := ParBound
CLASSDICT  := AngDict
CLASSDICT2 := RooDoubleCBDict

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(SOURCEDIR)/$(CLASS0).cc $(SOURCEDIR)/$(CLASS1).cc $(SOURCEDIR)/$(CLASS2).cc $(SOURCEDIR)/$(CLASS3).cc $(SOURCEDIR)/$(CLASS5).cc $(CLASSDICT).cc $(SOURCEDIR)/$(EXTRACLASS)
# LIBS := $(SOURCEDIR)/$(CLASS0).cc $(SOURCEDIR)/$(CLASS1).cc $(SOURCEDIR)/$(CLASS2).cc $(SOURCEDIR)/$(CLASS3).cc $(SOURCEDIR)/$(CLASS4).cc $(SOURCEDIR)/$(CLASS5).cc $(CLASSDICT).cc $(SOURCEDIR)/$(EXTRACLASS)

$(CLASSDICT): $(INCLUDEDIR)/$(CLASS0).h $(INCLUDEDIR)/$(CLASS1).h $(INCLUDEDIR)/$(CLASS2).h $(INCLUDEDIR)/$(CLASS3).h $(INCLUDEDIR)/$(CLASS5).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^

$(CLASSDICT2): $(INCLUDEDIR)/$(CLASS4).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^ -I./vdt	

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE1): $(EXECUTABLE1).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT2).cc $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE2): $(EXECUTABLE2).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)



#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE0) $(EXECUTABLE1) $(EXECUTABLE2)

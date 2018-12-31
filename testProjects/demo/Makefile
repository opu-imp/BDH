CXX = /usr/local/gcc5.1/bin/g++ -lstdc++ -std=c++14 -Wall -O3 
#CXX = g++ -Wall -g -std=c++11 -I./../../build/include -I./../Utility/include -L./../../build/lib -BDH_x64_Debug.lib

PACKAGE = NN_Search

SRCS=	demoMain.cpp\
	./../Utility/src/FileHandle.cpp\
	./../Utility/src/measure.cpp\
	./../Utility/src/superPCA.cpp\
	./../../BDH/src/BDH.cpp\
	./../../BDH/src/BDHtraining.cpp\
	./../../BDH/src/HashTable.cpp\
	./../../BDH/src/Subspace.cpp

HEDS =	

CXXINCLUDES =	-I./../Utility/include\
		-I./../include\
		-I./../../build/include \
		  -I /usr/include/octave-3.8.2

LIBS = -loctave

OBJS = $(subst .c,.o, $(subst .cpp,.o, $(SRCS) ) )
#OBJS = $(SRCS:.c=.o)
FILES = Makefile $(HEADS) $(SRCS)

all:	$(PACKAGE)

$(PACKAGE):	$(OBJS)
	$(CXX) -o $(PACKAGE) $(OBJS) $(LIBS)
.cpp.o:
	(cd `dirname $<` && \
	$(CXX) $(CXXINCLUDES) -c `basename $<` )
.cc.o:
	(cd `dirname $<` && \
	$(CXX) $(CXXINCLUDES) -c `basename $<` )
.c.o:
	(cd `dirname $<` && \
	$(CXX) $(CXXINCLUDES) -c `basename $<` )

#$(OBJS):	$(SRCS)
#	$(CXX) -c $(SRCS)
$(OBJS): $(HEADS) Makefile

clean:	
	rm -f $(OBJS)

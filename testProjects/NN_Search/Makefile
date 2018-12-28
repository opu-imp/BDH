CXX = /usr/local/gcc5.1/bin/gcc -lstdc++ -std=c++14 -Wall -O3 
#CXX = g++ -Wall -g -std=c++11 -I./../../build/include -I./../Utility/include -L./../../build/lib -BDH_x64_Debug.lib

PACKAGE = NN_Search

SRCS=	main_NN_Search.cpp\
	./../Utility/src/FileHandle.cpp\
	./../Utility/src/measure.cpp\
	./../../BDH/src/BDH.cpp\
	./../../BDH/src/BDHtraining.cpp\
	./../../BDH/src/HashTable.cpp\
	./../../BDH/src/Subspace.cpp

HEDS =	

CXXINCLUDES =	-I./../Utility/include\
		-I./../include\
		-I./../../build/include
              

OBJS = $(subst .c,.o, $(subst .cpp,.o, $(SRCS) ) )
#OBJS = $(SRCS:.c=.o)
FILES = Makefile $(HEADS) $(SRCS)

all:	$(PACKAGE)

$(PACKAGE):	$(OBJS)
	$(CXX) -o $(PACKAGE) $(OBJS)
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

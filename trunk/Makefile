OBJS = datastructures.o random_graph.o centrality.o community_tools.o graph_operations.o paths_and_components.o graphio.o graph.o graph_properties.o statistics.o graph_summary.o community.o
OBJS_D = datastructures_d.o random_graph_d.o centrality_d.o community_tools_d.o graph_operations_d.o paths_and_components_d.o graphio_d.o graph_d.o graph_properties_d.o statistics_d.o graph_summary_d.o community_d.o

CC = g++
CFLAGS = -O3 -fPIC -fopenmp -std=c++0x -DNDEBUG -DENABLE_MULTITHREADING
CFLAGS_D = -g -fPIC -fopenmp -DENABLE_MULTITHREADING -std=c++0x -Wall
LIBS = -lm -lgomp -lrt -lpthread

libcdlib.so : $(OBJS)
	$(CC) -std=c++0x -shared -o libcdlib.so $(OBJS) $(LIBS)
	$(CC) -std=c++0x -shared -o libcdlib_d.so $(OBJS_D) $(LIBS)
	

community_tools.o : graph.o graph_operations.o statistics.o paths_and_components.o
	$(CC) $(CFLAGS) -o community_tools.o  -c community_tools.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o community_tools_d.o  -c community_tools.cpp $(LIBS)
	
graph_properties.o : graph.o paths_and_components.o statistics.o
	$(CC) $(CFLAGS) -o graph_properties.o  -c graph_properties.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_properties_d.o  -c graph_properties.cpp $(LIBS)

graph_summary.o : graph.o paths_and_components.o graph_properties.o statistics.o
	$(CC) $(CFLAGS) -o graph_summary.o  -c graph_summary.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_summary_d.o  -c graph_summary.cpp $(LIBS)

community.o : graph.o community_tools.o graphio.o graph_operations.o centrality.o paths_and_components.o 
	$(CC) $(CFLAGS) -o community.o  -c community.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o community_d.o  -c community.cpp $(LIBS)

random_graph.o : graph.o community_tools.o paths_and_components.o
	$(CC) $(CFLAGS) -o random_graph.o  -c random_graph.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o random_graph_d.o  -c random_graph.cpp $(LIBS)	

centrality.o : graph.o datastructures.o
	$(CC) $(CFLAGS) -o centrality.o  -c centrality.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o centrality_d.o  -c centrality.cpp $(LIBS)

graph_operations.o  : graph.o
	$(CC) $(CFLAGS) -o graph_operations.o  -c graph_operations.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_operations_d.o  -c graph_operations.cpp $(LIBS)

paths_and_components.o  : graph.o datastructures.o
	$(CC) $(CFLAGS) -o paths_and_components.o  -c paths_and_components.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o paths_and_components_d.o  -c paths_and_components.cpp $(LIBS)

graphio.o : graph.o paths_and_components.o graph_operations.o
	$(CC) $(CFLAGS) -o graphio.o -c graphio.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graphio_d.o -c graphio.cpp $(LIBS)

graph.o : datastructures.o
	$(CC) $(CFLAGS) -o graph.o -c graph.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_d.o -c graph.cpp $(LIBS)

datastructures.o : 
	$(CC) $(CFLAGS) -o datastructures.o -c datastructures.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o datastructures_d.o -c datastructures.cpp $(LIBS)

statistics.o : 
	$(CC) $(CFLAGS) -o statistics.o -c statistics.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o statistics_d.o -c statistics.cpp $(LIBS)

clean:
	rm *.o *.so


#CC = g++
#CFLAGS = -O3 -fPIC -fopenmp -std=c++0x -DNDEBUG
#CFLAGS_D = -g -fPIC -fopenmp -std=c++0x -Wall
#LIBS = -lm -lgomp -lrt -lpthread
#
#all: libcdlib.so libcdlib_d.so
#	
#
#libcdlib.so:
#	$(CC) -std=c++0x -fPIC -O3 -c *.cpp $(LIBS)
#	$(CC) -std=c++0x -O3 -shared -o libcdlib.so $(LIBS)
#
#libcdlib_d.so:
#	$(CC) -std=c++0x -fPIC -g -c *.cpp $(LIBS)
#	$(CC) -std=c++0x -g -shared -o libcdlib_d.so $(LIBS)	
#
#clean:
#	rm *.o *.so
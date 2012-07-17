OBJS = label_propagation.o random_graph.o divisive_algorithms.o disjoint_set.o centrality.o community_tools.o graph_operations.o paths_and_components.o graphio.o graph.o binary_heap.o double_adjacency_map.o bidirectional_label_map.o graph_properties.o local_community.o
OBJS_D = label_propagation_d.o random_graph_d.o divisive_algorithms_d.o disjoint_set_d.o centrality_d.o community_tools_d.o graph_operations_d.o paths_and_components_d.o graphio_d.o graph_d.o binary_heap_d.o double_adjacency_map_d.o bidirectional_label_map_d.o graph_properties_d.o local_community_d.o

CC = g++
CFLAGS = -O3 -fPIC -fopenmp -std=c++0x -DNDEBUG
CFLAGS_D = -g -fPIC -fopenmp -std=c++0x -Wall
LIBS = -lm -lgomp -lrt -lpthread

libcdlib.so : $(OBJS)
	$(CC) -std=c++0x -shared -o libcdlib.so $(OBJS) $(LIBS)
	$(CC) -std=c++0x -shared -o libcdlib_d.so $(OBJS_D) $(LIBS)
	
graph_properties.o : graph.o
	$(CC) $(CFLAGS) -o graph_properties.o  -c graph_properties.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_properties_d.o  -c graph_properties.cpp $(LIBS)

dynamic_communities.o : graph.o graphio.o label_propagation.o
	$(CC) $(CFLAGS) -o dynamic_communities.o  -c dynamic_communities.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o dynamic_communities_d.o  -c dynamic_communities.cpp $(LIBS)

label_propagation.o : graph.o graph_operations.o community_tools.o
	$(CC) $(CFLAGS) -o label_propagation.o  -c label_propagation.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o label_propagation_d.o  -c label_propagation.cpp $(LIBS)

random_graph.o : graph.o
	$(CC) $(CFLAGS) -o random_graph.o  -c random_graph.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o random_graph_d.o  -c random_graph.cpp $(LIBS)
	
divisive_algorithms.o : graph.o centrality.o paths_and_components.o 
	$(CC) $(CFLAGS) -o divisive_algorithms.o  -c divisive_algorithms.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o divisive_algorithms_d.o  -c divisive_algorithms.cpp $(LIBS)

disjoint_set.o : graph.o
	$(CC) $(CFLAGS) -o disjoint_set.o  -c disjoint_set.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o disjoint_set_d.o  -c disjoint_set.cpp $(LIBS)

centrality.o : graph.o binary_heap.o
	$(CC) $(CFLAGS) -o centrality.o  -c centrality.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o centrality_d.o  -c centrality.cpp $(LIBS)

community_tools.o : graph.o graph_operations.o
	$(CC) $(CFLAGS) -o community_tools.o  -c community_tools.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o community_tools_d.o  -c community_tools.cpp $(LIBS)

graph_operations.o  : graph.o
	$(CC) $(CFLAGS) -o graph_operations.o  -c graph_operations.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_operations_d.o  -c graph_operations.cpp $(LIBS)

paths_and_components.o  : graph.o binary_heap.o random_graph.o
	$(CC) $(CFLAGS) -o paths_and_components.o  -c paths_and_components.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o paths_and_components_d.o  -c paths_and_components.cpp $(LIBS)

graphio.o : graph.o
	$(CC) $(CFLAGS) -o graphio.o -c graphio.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graphio_d.o -c graphio.cpp $(LIBS)


local_community.o : graph.o
	$(CC) $(CFLAGS) -o local_community.o -c local_community.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o local_community_d.o -c local_community.cpp $(LIBS)

graph.o : double_adjacency_map.o bidirectional_label_map.o
	$(CC) $(CFLAGS) -o graph.o -c graph.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o graph_d.o -c graph.cpp $(LIBS)

binary_heap.o : 
	$(CC) $(CFLAGS) -o binary_heap.o -c binary_heap.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o binary_heap_d.o -c binary_heap.cpp $(LIBS)

double_adjacency_map.o : 
	$(CC) $(CFLAGS) -o double_adjacency_map.o -c double_adjacency_map.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o double_adjacency_map_d.o -c double_adjacency_map.cpp $(LIBS)

bidirectional_label_map.o : 
	$(CC) $(CFLAGS) -o bidirectional_label_map.o -c bidirectional_label_map.cpp $(LIBS)
	$(CC) $(CFLAGS_D) -o bidirectional_label_map_d.o -c bidirectional_label_map.cpp $(LIBS)



clean:
	rm *.o *.so

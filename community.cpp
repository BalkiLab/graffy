#include "community.h"

using namespace std;
using namespace CDLib;



/*
 * Implements Aaron Clauset, "Finding Local Community Structure in Networks" 2005
 * 
 * Input: The graph g
 *        The seed/source vertex src
 *        The size of the community desired k
 * 
 * Pseudocode:
 * -------------------------------
 * add src to C
 * add all neighbors of src to U
 * set B = src
 * set R=0, T = degree(src)
 * 
 * while |C| < k do
 *   for each vj belonging to U do
 *     compute delta_Rj (R,T,x,y,z used here)
 *   end for
 *   find vj such that delta_Rj is maximum (breaking ties randomly if the need arises)
 *   add that vj to C
 *   add all new neighbors of that vj to U and remove vj from U
 *   update R, T and B
 * end while
 * -------------------------------
 */

//Other functions go here

/*
 * Main Function
 * 
 * Input: The graph g
 *        The seed vertex src
 *        The size of the community desired k
 * 
 * Output: 
 */

double get_neighs_in_set(id_type v, const graph& g, node_set& B, node_set& X) {
    X.clear();
    double retval = 0;
    for (adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++) {
        if (B.find(aeit->first) != B.end()) {
            X.insert(aeit->first);
            retval += 1;
        }
    }
    return retval;
}

double get_neighs_not_in_sets(id_type v, const graph& g, node_set& C, node_set& U, node_set& X) {
    X.clear();
    double retval = 0;
    for (adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++) {
        if (C.find(aeit->first) == C.end() && U.find(aeit->first) == U.end()) {
            X.insert(aeit->first);
            retval += 1;
        }
    }
    return retval;
}

bool membership(id_type vk, id_type v, const graph& g, node_set& U) {
    for (adjacent_edges_iterator aeit = g.out_edges_begin(vk); aeit != g.out_edges_end(vk); aeit++)
        if ((aeit->first) != v && U.find(aeit->first) != U.end())
            return 1;
    return 0;
}

double no_BI_edges(id_type v, const graph& g, node_set& C, node_set& B) {
    double count = 0;
    for (adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
        if (C.find(aeit->first) != C.end() && B.find(aeit->first) == B.end())
            count += 1;
    return count;
}

double no_BY_edges(id_type v, const graph& g, node_set& Y) {
    double count = 0;
    for (adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
        if (Y.find(aeit->first) != Y.end())
            count += 1;
    return count;
}



//-----------------------------------------------------------------------------------------------------------------------------------------

bool CDLib::local_community_clauset(const graph& g, id_type src, size_t k, node_set& output) {
    vector< pair<id_type, double> > output2;
    node_set C, U, B;
    C.insert(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first);
    B.insert(src);

    double R = 0;
    double T = g.get_node_out_degree(src);

    //Update output2
    output2.push_back(make_pair(src, R));
    output.insert(src);

    while (C.size() < k && !U.empty()) {
        double max = -numeric_limits<double>::infinity();
        long cntr = -1;
        long no_max = 1;
        vector<id_type> V; //DOUBT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        unordered_map<id_type, node_set> temp_Y; //To store the necessary v and Y pairs 
        unordered_map<id_type, double> temp_T; //To store the necessary v and T pairs        

        node_set X, Y;
        double x, y, z;

        for (node_set::iterator iter = U.begin(); iter != U.end(); iter++) {
            //compute delta_Rj
            double delta;

            id_type v = *iter; //vertex vj

            x = get_neighs_in_set(v, g, B, X); //X = set of all neighbors of vj in B   //x = No. of edges in T that terminate at vj

            y = g.get_node_out_degree(v) - x; //No. of edges that will be added to T by the agglomeration of vj

            z = 0; //No. of edges that will be removed from T by the agglomeration of vj
            Y.clear(); // DEBUGGED!!!!
            for (node_set::iterator iter1 = X.begin(); iter1 != X.end(); iter1++) {
                id_type vk = *iter1;
                if (membership(vk, v, g, U) == 0) //vk has no neighbor in U other than vj, i.e., vk no longer belongs to B 
                {
                    Y.insert(vk);
                    z += no_BI_edges(vk, g, C, B) + no_BY_edges(vk, g, Y);
                    if (y == 0)
                        z = z + 1;
                }
            }

            delta = (x - (R * y) - z * (1 - R)) / (T - z + y);

            if (delta > max) {
                V.push_back(v);
                cntr += no_max;
                max = delta;
                no_max = 1;
                //store v and corresponding Y
                temp_Y.insert(make_pair(v, Y));
                //store v and corresponding new value of T
                temp_T.insert(make_pair(v, T - z + y));
            } else if (delta == max) {
                V.push_back(v);
                no_max += 1;
                //store v and corresponding Y
                temp_Y.insert(make_pair(v, Y));
                //store v and corresponding new value of T
                temp_T.insert(make_pair(v, T - z + y));
            }
        }

        //Selecting the node to be agglomerated v_aggl (breaking ties randomly if needed)
        //Choose a random number rand in [cntr,V.size()-1] such that v_aggl=V[rand] 
        RandomGenerator<long> rnd_gen(cntr, V.size() - 1);
        id_type v_aggl = V[rnd_gen.next()];

        //Update B
        node_set Y_aggl; //corresponding Y of v_aggl
        /*vector<pair<id_type,node_set>>::iterator iter2; //Doubt
        for(iter2=temp_Y.begin();iter2<temp_Y.end();iter2++) //Doubt
            if(iter2->first==v_aggl)
            {
                Y_aggl=iter2->second;
                break;
            }*/
        Y_aggl = (temp_Y.find(v_aggl))->second;

        x = get_neighs_in_set(v_aggl, g, B, X);
        y = g.get_node_out_degree(v_aggl) - x;

        for (node_set::iterator iter3 = Y_aggl.begin(); iter3 != Y_aggl.end(); iter3++)
            B.erase(*iter3);

        if (y > 0)
            B.insert(v_aggl);

        //Update U
        U.erase(v_aggl);
        get_neighs_not_in_sets(v_aggl, g, C, U, X); //X = all neighbors of v_aggl not in C
        for (node_set::iterator iter4 = X.begin(); iter4 != X.end(); iter4++)
            U.insert(*iter4);

        //Update C
        C.insert(v_aggl);

        //Update R
        R = R + max;

        //Update T
        /*vector<pair<id_type,size_t>>::iterator iter5; //Doubt
        for(iter5=temp_T.begin();iter5<temp_T.end();iter5++) //Doubt
            if(iter5->first==v_aggl)
                T=iter5->second;*/
        T = (temp_T.find(v_aggl))->second;

        //Update output2
        output2.push_back(make_pair(v_aggl, R));
        output.insert(v_aggl);
    }

    if (output.size() == k)
        return 1; //comm found/exists
    else
        return 0; //comm not found/doesn't exist

}
//-----------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------

bool CDLib::local_community_clauset_modified(const graph& g, id_type src, size_t k, node_set& output) {
    vector< pair<id_type, double> > output2;
    node_set C, U;
    C.insert(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first);

    double Z = 0;
    double I = 0;
    double T = g.get_node_out_degree(src);

    //Update output2
    output2.push_back(make_pair(src, Z));
    output.insert(src);

    while (C.size() < k && !U.empty()) {
        double max = -numeric_limits<double>::infinity();
        long cntr = -1;
        long no_max = 1;
        vector<id_type> V;

        node_set X;
        double x, y;

        for (node_set::iterator iter = U.begin(); iter != U.end(); iter++) {
            //compute delta_Zj
            double delta;

            id_type v = *iter; //vertex vj

            x = get_neighs_in_set(v, g, C, X); //X = set of all neighbors of vj in C  //x = No. of edges in T that terminate at vj

            y = g.get_node_out_degree(v) - x; //No. of edges that will be added to T by the agglomeration of vj

            delta = ((I + x) / (T + y)) - Z;

            if (delta > max) {
                V.push_back(v);
                cntr += no_max;
                max = delta;
                no_max = 1;
                //store v and corresponding new value of T
                //temp_T.insert(make_pair(v,T-z+y));
            } else if (delta == max) {
                V.push_back(v);
                no_max += 1;
                //store v and corresponding new value of T
                //temp_T.insert(make_pair(v,T-z+y));
            }
        }

        //Selecting the node to be agglomerated v_aggl (breaking ties randomly if needed)
        //Choose a random number rand in [cntr,V.size()-1] such that v_aggl=V[rand] 
        RandomGenerator<long> rnd_gen(cntr, V.size() - 1); //Doubt
        id_type v_aggl = V[rnd_gen.next()]; //Doubt

        //Update I and T
        x = get_neighs_in_set(v_aggl, g, C, X);

        y = g.get_node_out_degree(v_aggl) - x;

        I = I + x;
        T = T + y;

        //Update U
        U.erase(v_aggl);
        get_neighs_not_in_sets(v_aggl, g, C, U, X); //X = all neighbors of v_aggl not in C
        for (node_set::iterator iter4 = X.begin(); iter4 != X.end(); iter4++)
            U.insert(*iter4);

        //Update C
        C.insert(v_aggl);

        //Update Z
        Z = Z + max;

        //Update output2
        output2.push_back(make_pair(v_aggl, Z));
        output.insert(v_aggl);
    }

    if (output.size() == k)
        return 1; //comm found/exists
    else
        return 0; //comm not found/doesn't exist
}
//-----------------------------------------------------------------------------------------------------------------------------------------

struct deg_comp {
    const graph* graph_ref;

    deg_comp(const graph * g) {
        graph_ref = g;
    }

    bool operator() (id_type i, id_type j) {
        return (graph_ref->get_node_out_degree(i) <= graph_ref->get_node_out_degree(j));
    }
};

bool degpr_comp(const pair<id_type, id_type>& lhs, const pair<id_type, id_type>& rhs) {
    return lhs.first < rhs.first;
}

void sort_acc_degrees(vector<id_type>& N, const graph& g) {
    //deg_comp this_object_compares_degrees(&g);
    //sort(N.begin(),N.end(),degpr_comp);
    vector<pair<id_type, id_type> > degpr(g.get_num_nodes(), make_pair(0, 0));
    for (id_type i = 0; i < N.size(); i++) {
        degpr[i].first = g.get_node_out_degree(N[i]);
        degpr[i].second = N[i];
    }
    for (id_type i = 0; i < N.size(); i++)N[i] = degpr[i].second;

}

void bubble_sort_acc_degrees(vector<id_type>& N, const graph& g) {
    id_type temp;
    for (size_t i = 0; i < (N.size() - 1); i++)
        for (size_t j = N.size() - 1; j > i; j--)
            if (g.get_node_out_degree(N[j]) < g.get_node_out_degree(N[j - 1])) {
                //swap
                temp = N[j];
                N[j] = N[j - 1];
                N[j - 1] = temp;
            }
}

bool found_in_vector(id_type v, const vector<id_type>& S, size_t& index) {
    for (size_t i = 0; i < S.size(); i++) //Linear Search
        if (v == S[i]) {
            index = i;
            return 1;
        }
    return 0;
}

double get_no_neighs_in_vector(id_type v, const graph& g, const vector<id_type>& subgraph) {
    //create an unordered set copy of the subgraph vector for time complexity issues
    //    node_set S;
    //    for(size_t i=0; i<subgraph.size(); i++)
    //        S.insert(subgraph[i]);
    //    
    //    double retval = 0;
    //    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
    //    {
    //        if(S.find(aeit->first) != S.end())
    //            retval += 1;
    //    }
    //    return retval;
    double retval = 0;
    for (size_t i = 0; i < subgraph.size(); i++)
        retval += g.get_edge_weight(v, subgraph[i]);
    return retval;
}

void delete_from_vector_at_position(id_type v, vector<id_type>& N, size_t i) {
    if (i >= 0 && N.size() && i < N.size()) {
        //for(size_t j=i; j<N.size()-1; j++)
        //    N[j]=N[j+1];
        // N.pop_back();
        N.erase(N.begin() + i);
    }

}

void bfs_visitor_in_subgraph(const graph& g, const vector<id_type>& subgraph, id_type source, node_set& visited) {
    //create an unordered set copy of the subgraph vector for time complexity issues
    node_set S;
    for (size_t i = 0; i < subgraph.size(); i++)
        S.insert(subgraph[i]);

    queue<id_type> q_bfs;
    q_bfs.push(source); //push the source node to be visited first
    while (!q_bfs.empty()) //you have nodes to visit
    {
        id_type current = q_bfs.front();
        q_bfs.pop();
        if (visited.find(current) == visited.end()) //current node has not yet been visited
        {
            //visit node
            visited.insert(current);
            //expand node
            for (adjacent_edges_iterator aeit = g.out_edges_begin(current); aeit != g.out_edges_end(current); aeit++) {
                if (visited.find(aeit->first) == visited.end() && S.find(aeit->first) != S.end()) {
                    //push neighbors to be visited
                    q_bfs.push(aeit->first); //insert only those neighbors which have not been "visited" yet
                }
            }
        }
    }
}

size_t get_component_around_node_in_subgraph(const graph& g, const vector<id_type>& S, id_type source, node_set& visited) {
    bfs_visitor_in_subgraph(g, S, source, visited);
    return visited.size();
}

bool connected_on_removal_at_position(const graph& g, const vector<id_type>& S, id_type u, size_t i) {
    //Note: S has been passed by value to avoid reflection of changes to the original subgraph S after deletion in this function
    //removing node u at position i in subgraph vector S in graph G
    vector<id_type> Sc(S.size(), 0);
    for (id_type k = 0; k < S.size(); k++)Sc[k] = S[k];
    delete_from_vector_at_position(u, Sc, i);
    node_set visited;
    return (get_component_around_node_in_subgraph(g, Sc, Sc[0], visited) == Sc.size());
}



//-----------------------------------------------------------------------------------------------------------------------------------------

bool CDLib::LWP_2006(const graph& g, id_type src, node_set& output) {
    output.clear();
    vector<id_type> S, N, Q, deleteQ;
    size_t index;
    //initialization
    S.push_back(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        N.push_back(aeit->first);

    double I = 0, E = g.get_node_out_degree(src), M = 0;
    size_t iter = 0;

    do {
        //vector<id_type> Q;
        iter = iter + 1;
        //cout<<"\nITERATION "<<iter<<":"<<endl;

        Q.clear();

        //ADDITION STEP
        sort_acc_degrees(N, g); //sort N in order of increasing degrees
        //sort_acc_degrees(N,g); //sort N in order of increasing degrees

        for (long i = 0; i < (long) N.size(); i++) //Doubt (i's datatype to allow for value -1)
        {
            id_type v = N[i]; //vertex vj
            //compute delta_M
            double x = get_no_neighs_in_vector(v, g, S);
            double y = g.get_node_out_degree(v) - x;
            double delta = (I + x) / (E - x + y) - M;

            if (delta > 0) {
                //cout<<"\ndelta is (insertion):"<<delta<<endl;
                //cout<<"Node: "<<g.get_node_label(v)<<endl;
                S.push_back(v); //push into S
                delete_from_vector_at_position(v, N, i); //pop from N
                Q.push_back(v);
                I = I + x;
                E = E - x + y;
                M += delta;
                i = i - 1; //shift needed to prevent skipping next element in N
            }
        }

        //S, N and Q have been updated at this point

        //DELETION STEP
        size_t counter = 0;
        do {
            counter += 1;
            deleteQ.clear();
            for (long i = 0; i < (long) S.size(); i++) //Doubt (i's datatype to allow for value -1)
            {
                //if(S[i]==src) //Modification to the original algorithm: Deletion of the source vertex is never allowed
                //  continue; 

                id_type u = S[i]; //vertex ui
                //compute delta_M
                double x = get_no_neighs_in_vector(u, g, S);
                double y = g.get_node_out_degree(u) - x;
                double delta = (I - x) / (E - y + x) - M;

                if (delta > 0 && connected_on_removal_at_position(g, S, u, i)) {
                    //cout<<"\ndelta is (deletion):"<<delta<<endl;
                    //cout<<"Node: "<<g.get_node_label(u)<<endl;
                    delete_from_vector_at_position(u, S, i);

                    //N.push_back(u); //Not a part of this algorithm as per Bagrow's paper

                    deleteQ.push_back(u);
                    I = I - x;
                    E = E - y + x; //debugged!!
                    M += delta;
                    if (found_in_vector(u, Q, index))
                        delete_from_vector_at_position(u, Q, index);
                    i = i - 1; //shift needed to prevent skipping next element in S
                }
            }
        } while (!deleteQ.empty());

        //cout<<"\nNo. of deletion loops in iteration "<<iter<<" are "<<counter<<endl;

        //ADD VERTICES TO N
        for (size_t i = 0; i < Q.size(); i++)
            for (adjacent_edges_iterator aeit = g.out_edges_begin(Q[i]); aeit != g.out_edges_end(Q[i]); aeit++)
                if (found_in_vector(aeit->first, S, index) == 0 && found_in_vector(aeit->first, N, index) == 0)
                    N.push_back(aeit->first);

    } while (!Q.empty());

    //cout<<"\nNo. of iterations: "<<iter<<endl;

    if (M > 1 && found_in_vector(src, S, index)) {
        output.clear();
        //cout<<endl<<I<<" "<<E<<" "<<M<<endl;
        for (size_t i = 0; i < S.size(); i++)
            output.insert(S[i]);
        deleteQ.clear();
        Q.clear();
        N.clear();
        S.clear();
        return 1; //comm found/exists
    } else {
        //cout<<"\nNO MODULE FOUND BY LWP FOR SOURCE VERTEX "<<g.get_node_label(src)<<endl;
        //cout<<"\nNO MODULE FOUND BY LWP FOR SOURCE VERTEX "<<src<<endl;        
        //output2 = S;
        output.clear();
        output.insert(src);
        deleteQ.clear();
        Q.clear();
        N.clear();
        S.clear();
        return 0; //comm not found/doesn't exist
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------

double get_no_neighs_in_set(id_type v, const graph& g, node_set& C) {
    double retval = 0;
    for (adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
        if (C.find(aeit->first) != C.end())
            retval += aeit->second;
    return retval;
}

bool not_found_in_deque(id_type w, deque<id_type>& Q) {
    for (size_t i = 0; i < Q.size(); i++)
        if (Q[i] == w)
            return 0;
    return 1;
}

void delete_from_deque(id_type u, deque<id_type>& Q) {
    for (size_t i = 0; i < Q.size(); i++) {
        if (Q[i] == u) {
            //delete from deque
            for (size_t j = i; j < Q.size() - 1; j++)
                Q[j] = Q[j + 1];
            Q.pop_back();
            i = i - 1; //shift needed to prevent skipping next element in Q
        }
    }
}



//-----------------------------------------------------------------------------------------------------------------------------------------

void CDLib::VD_2011(const graph& g, id_type src, node_set& output) {
    deque<id_type> Q;
    node_set C;
    node_set visited;
    id_type u = src;
    C.insert(u);
    double I = 0, E = g.get_node_out_degree(u), N = 1, f, x, y, I2, E2, N2, f2;
    size_t iter_count = 0, N_at_start;

    //Initialize C with the source vertex and its neighbors AND update I, E and N
    for (adjacent_edges_iterator aeit = g.out_edges_begin(u); aeit != g.out_edges_end(u); aeit++) {
        x = get_no_neighs_in_set(aeit->first, g, C);
        y = g.get_node_out_degree(aeit->first) - x;
        C.insert(aeit->first);
        I = I + x; //No. of internal edges of C
        E = E - x + y; //No. of external edges of C
        N = N + 1; //No. of nodes in C
    }

    f = (2 * I * I) / (N * (N - 1)*(I + E)); //the objective function

    visited.insert(u); //the source vertex has now been visited

    //push all neighbors of src in Q for visiting them
    for (adjacent_edges_iterator aeit = g.out_edges_begin(u); aeit != g.out_edges_end(u); aeit++) {
        Q.push_back(aeit->first);
    }

    //cout<<"\n\nOutput C before iterations begin:\n\n";

    //for(node_set::iterator iter=C.begin(); iter!=C.end(); iter++)
    //cout<<g.get_node_label(*iter)<<endl<<endl;

    do {
        iter_count += 1;

        //cout<<"\n\nITERATION NO. "<<iter_count<<endl<<endl;

        N_at_start = C.size(); //size of C at the start of the iteration

        //ADDITION PHASE
        if (iter_count != 1 && !Q.empty()) {
            u = Q.front();
            Q.pop_front();
            //}

            visited.insert(u); //visit the node

            //push the required neighbors of u in Q and try agglomerating them into C
            for (adjacent_edges_iterator aeit = g.out_edges_begin(u); aeit != g.out_edges_end(u); aeit++) {
                id_type w = aeit->first;
                if (visited.find(w) == visited.end() && not_found_in_deque(w, Q)) //Doubt
                    Q.push_back(w);
                if (C.find(w) == C.end()) {
                    x = get_no_neighs_in_set(w, g, C);
                    y = g.get_node_out_degree(w) - x;
                    //C.insert(w);

                    I2 = I + x;
                    E2 = E - x + y;
                    N2 = N + 1;
                    f2 = (2 * I2 * I2) / (N2 * (N2 - 1)*(I2 + E2));

                    if (f2 < f) {
                        //C.erase(w);
                    } else {
                        C.insert(w);
                        I = I2;
                        E = E2;
                        N = N2;
                        f = f2;
                        //cout<<"\nNode inserted: "<<g.get_node_label(w)<<endl;
                    }
                }
            }
        }

        //DELETION PHASE
        for (node_set::iterator iter = C.begin(); iter != C.end();) {
            u = *iter;
            x = get_no_neighs_in_set(u, g, C);
            y = g.get_node_out_degree(u) - x;

            //cout<<g.get_node_label(u)<<"  ";
            //cout<<"\n\nInitial erase from C successful for node \n\n"<<g.get_node_label(u)<<endl<<endl;

            I2 = I - x;
            E2 = E + x - y;
            N2 = N - 1;
            f2 = (2 * I2 * I2) / (N2 * (N2 - 1)*(I2 + E2));

            if (f2 < f) {
                iter++;
            } else {
                iter++;
                C.erase(u);
                I = I2;
                E = E2;
                N = N2;
                f = f2;
                //cout<<"\nNode deleted: "<<g.get_node_label(u)<<endl;
                delete_from_deque(u, Q);
            }
        }

        //cout<<"\nN is "<<N<<endl;
        //cout<<"\nSize of C is "<<C.size()<<endl;

        //cout<<"\nQueue content at the end of iteration "<<iter_count<<endl;
        //for(size_t i=0; i<Q.size(); i++)
        //cout<<g.get_node_label(Q[i])<<"  ";

    } while (N_at_start != C.size()); //go to next iteration if community size changes, else exit

    //cout<<"\n\nNo. of iterations: "<<iter_count<<endl<<endl;

    //return C;
    output = C;
}
//-----------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------

/*pair<double,node_set> CDLib::Bagrow_2007(const graph& g, id_type src)
{
    //const double max_threshold = 0.1;
    node_set C, B;
    C.insert(src);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        B.insert(aeit->first);
    
    //map<double,id_type> omegas;
    id_dbl_min_heap omegas;
    double omega_val, p_eff = 0, x, y, xi, yi, xf, yf, p_eff_numerator = 0;
    double M_out = g.get_node_out_degree(src);
    id_type v, node_aggl, u, z;
    double threshold = 0.10;
    map<double,node_set> choices; //no duplicates allowed--flaw
    bool flag=1;
    
    //Initializing omegas
    for(node_set::iterator iter=B.begin();iter!=B.end();iter++)
    {
        v = *iter;
        omega_val = 1 - 2/g.get_node_out_degree(v);
        omegas.insert(make_pair(v,omega_val));   
    }
    
    //cout<<"initializations complete";
    //cout.flush();
    
    do
    {
        //extract the node with the minimum omega_val
        id_dbl_min_heap::iterator iter2 = omegas.begin();
        node_aggl = iter2->first; //no breaking ties randomly??????!!!!!!
        omegas.pop();
        
        x = get_no_neighs_in_set(node_aggl, g, C);
        y = g.get_node_out_degree(node_aggl) - x;
        
        for(adjacent_edges_iterator aeit = g.out_edges_begin(node_aggl); aeit != g.out_edges_end(node_aggl); aeit++)
            if(C.find(aeit->first) != C.end())
            {
                z = aeit->first;
                xi = get_no_neighs_in_set(z, g, C);
                yi = g.get_node_out_degree(z) - xi;
                xf = xi+1;
                yf = yi-1;
                if(xi<yi && xf>yf)
                    p_eff_numerator+=1;
            }
                
        
        C.insert(node_aggl);
        
        M_out = M_out-x+y;
        
        //cout<<"\n\nNode aggl: "<<g.get_node_label(node_aggl);
        //cout<<"\nnew M_out: "<<M_out;
        
        if(x>y)
            p_eff_numerator+=1;
        
        p_eff = p_eff_numerator/C.size();
        cout<<"\nnew p_eff: "<<p_eff;
        
        if(p_eff>=threshold)
        {
            //if(p_eff!=1.0)    //equality!!!!
            if(M_out!=0)   
                choices.insert(make_pair(M_out,C));
            cout<<"\nHIT P_EFF: C SAVED";
            cout<<"\ncurr threshold: "<<threshold <<endl;
            if(threshold==1) //equality!!!!!
            {
                flag=0;
                cout<<"flag is now 0";
                cout.flush();
            }
            else
            { 
                threshold+=0.01*((floor((p_eff-threshold)/0.01))+1); //check!!!!!!!!!
                cout<<"\nNEW threshold: "<<threshold <<endl;
                //threshold+=0.01;
            }
        }
        
        B.erase(node_aggl);
        //cout<<"\nnew content of B: "<<endl;
        //for(node_set::iterator iter5=B.begin(); iter5!=B.end(); iter5++)
        //    cout<<g.get_node_label(*iter5)<<endl<<endl;
        //cout<<"\nnew flag: "<<flag;
        
        //Update omegas
        for(adjacent_edges_iterator aeit = g.out_edges_begin(node_aggl); aeit != g.out_edges_end(node_aggl); aeit++)
        {
            u = aeit->first;
            if(B.find(u) == B.end() && C.find(u) == C.end())
            {
                omega_val = 1 - 2/g.get_node_out_degree(u);
                omegas.insert(make_pair(u,omega_val));
                B.insert(u);
            }
            else if (B.find(u) != B.end())
            {
                id_dbl_min_heap::iterator iter3 = omegas.find(u);
                omega_val = iter3->second;
                omega_val-=2/g.get_node_out_degree(u);
                omegas.update_key(make_pair(u,omega_val)); //redundant search (but constant time)
            }
        }
        
    } while(flag==1 && !B.empty()) ;
    
    //return C;
    
    map<double,node_set>::iterator iter4 = choices.begin();
    
    return *iter4;
    //output = make_pair(iter4->first,iter4->second); //in case of duplicates????????
}
 */

bool CDLib::Bagrow_2007(const graph& g, id_type src, node_set& output) {
    //const double max_threshold = 0.1;
    node_set C, B, temp;
    C.insert(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        B.insert(aeit->first);

    //map<double,id_type> omegas;
    id_dbl_min_heap omegas;
    double omega_val, p_eff = 0, x, y, xi, yi, xf, yf, p_eff_numerator = 0;
    double M_out = g.get_node_out_degree(src);
    id_type v, node_aggl, u, z;
    double start_val = 0.75;
    double threshold = start_val; //dependency on arbitrary parameters
    multimap<double, node_set> choices; //no duplicates allowed--flaw------------RESOLVED!!
    bool flag = 1;

    //Initializing omegas
    for (node_set::iterator iter = B.begin(); iter != B.end(); iter++) {
        v = *iter;
        omega_val = 1 - 2 / g.get_node_out_degree(v);
        omegas.insert(make_pair(v, omega_val));
    }

    //cout<<"initializations complete";
    //cout.flush();

    do {
        //extract the node with the minimum omega_val
        id_dbl_min_heap::iterator iter2 = omegas.begin();
        node_aggl = iter2->first; //no breaking ties randomly??????!!!!!!-------NOT NEEDED I GUESS
        omegas.pop();

        x = get_no_neighs_in_set(node_aggl, g, C);
        y = g.get_node_out_degree(node_aggl) - x;

        for (adjacent_edges_iterator aeit = g.out_edges_begin(node_aggl); aeit != g.out_edges_end(node_aggl); aeit++)
            if (C.find(aeit->first) != C.end()) {
                z = aeit->first;
                xi = get_no_neighs_in_set(z, g, C);
                yi = g.get_node_out_degree(z) - xi;
                xf = xi + 1;
                yf = yi - 1;
                if (xi < yi && xf > yf)
                    p_eff_numerator += 1;
            }


        C.insert(node_aggl);

        M_out = M_out - x + y;

        //cout<<"\n\nNode aggl: "<<g.get_node_label(node_aggl);
        //cout<<"\nnew M_out: "<<M_out;

        if (x > y)
            p_eff_numerator += 1;

        p_eff = p_eff_numerator / C.size();
        //cout<<"\nnew p_eff: "<<p_eff;

        if (p_eff >= threshold) {
            //if(p_eff!=1.0)    //equality!!!!
            if (M_out != 0)
                choices.insert(make_pair(M_out, C)); //state of C saved
            //cout<<"\nHIT P_EFF: C SAVED";
            //cout<<"\ncurr threshold: "<<threshold <<endl;
            if (threshold == 1) //equality!!!!!
            {
                flag = 0;
                //cout<<"flag is now 0";
                //cout.flush();
            } else {
                threshold += 0.01 * ((floor((p_eff - threshold) / 0.01)) + 1); //checked
                //cout<<"\nNEW threshold: "<<threshold <<endl;
                //threshold+=0.01;
            }
        }

        B.erase(node_aggl);
        //cout<<"\nnew content of B: "<<endl;
        //for(node_set::iterator iter5=B.begin(); iter5!=B.end(); iter5++)
        //    cout<<g.get_node_label(*iter5)<<endl<<endl;
        //cout<<"\nnew flag: "<<flag;

        //Update omegas and B
        for (adjacent_edges_iterator aeit = g.out_edges_begin(node_aggl); aeit != g.out_edges_end(node_aggl); aeit++) {
            u = aeit->first;
            if (B.find(u) == B.end() && C.find(u) == C.end()) {
                omega_val = 1 - 2 / g.get_node_out_degree(u);
                omegas.insert(make_pair(u, omega_val));
                B.insert(u);
            } else if (B.find(u) != B.end()) {
                id_dbl_min_heap::iterator iter3 = omegas.find(u);
                omega_val = iter3->second;
                omega_val -= 2 / g.get_node_out_degree(u);
                omegas.update_key(make_pair(u, omega_val)); //redundant search (but constant time)
            }
        }

    } while (flag == 1 && !B.empty());

    //return C;

    //map<double,node_set>::iterator iter4 = choices.begin();
    if (threshold > start_val) {
        multimap<double, node_set>::iterator iter4 = choices.begin();
        double min_val = iter4->first;
        int min_count = choices.count(min_val);
        if (min_count == 1) {
            output = iter4->second;
        } else {
            double max = -numeric_limits<double>::infinity();
            for (int i = 1; i <= min_count; i++) {
                //calculate p_eff for iter4->second
                temp.clear();
                temp = iter4->second;

                p_eff_numerator = 0;
                for (node_set::iterator iter5 = temp.begin(); iter5 != temp.end(); iter5++) {
                    v = *iter5;
                    x = get_no_neighs_in_set(v, g, temp);
                    y = g.get_node_out_degree(v) - x;
                    if (x > y)
                        p_eff_numerator += 1;
                }
                p_eff = p_eff_numerator / temp.size();

                if (p_eff > max) {
                    max = p_eff;
                    output.clear();
                    output = temp;
                }

                iter4++;
            }
        }
        return 1; //comm found/exists
    } else {
        output.insert(src);
        return 0; //comm not found/doesn't exist
    }

    //output = make_pair(iter4->first,iter4->second); //in case of duplicates????????
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void bfs_visitor_in_subgraph_set(const graph& g, node_set& C, id_type source, node_set& visited) {
    //create an unordered set copy of the subgraph vector for time complexity issues
    //node_set S;
    //for(size_t i=0; i<subgraph.size(); i++)
    //    S.insert(subgraph[i]);

    queue<id_type> q_bfs;
    q_bfs.push(source); //push the source node to be visited first
    while (!q_bfs.empty()) //you have nodes to visit
    {
        id_type current = q_bfs.front();
        q_bfs.pop();
        if (visited.find(current) == visited.end()) //current node has not yet been visited
        {
            //visit node
            visited.insert(current);
            //expand node
            for (adjacent_edges_iterator aeit = g.out_edges_begin(current); aeit != g.out_edges_end(current); aeit++) {
                if (visited.find(aeit->first) == visited.end() && C.find(aeit->first) != C.end()) {
                    //push neighbors to be visited
                    q_bfs.push(aeit->first); //insert only those neighbors which have not been "visited" yet
                }
            }
        }
    }
}

size_t get_component_around_node(const graph& g, node_set& C, id_type source, node_set& visited) {
    bfs_visitor_in_subgraph_set(g, C, source, visited);
    return visited.size();
}

bool connected_on_removal(const graph& g, node_set C, id_type u) {
    //Note: C has been passed by value to avoid reflection of changes to the original subgraph C after deletion in this function
    //deleting
    C.erase(u);
    node_set visited;
    id_type start_node = *(C.begin());
    return (get_component_around_node(g, C, start_node, visited) == C.size());
}
//-----------------------------------------------------------------------------------------------------------------------------------------

bool CDLib::My_Algorithm(const graph& g, id_type src, node_set& output) {
    node_set C, U;
    vector<id_type> U_vector;
    C.insert(src);
    //C_vector.push_back(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first);
    //for(node_set::iterator iter=U.begin(); iter!=U.end();iter++)
    //    U_vector.push_back(*iter);

    double I = 0, E = g.get_node_out_degree(src), n = 1, nb = 1, nu = g.get_node_out_degree(src);
    double I_new, E_new, n_new, nb_new, nu_new;
    double rho = 0, delta_rho, x, y;

    vector<id_type> Q, deleteQ;
    bool flag;
    id_type vj, ui, v;
    size_t index;

    do {
        //cout<<"\n\nNEW ITERATION STARTS:\n\n";
        //ADDITION PHASE
        //cout<<"\nentered addition phase\n";
        //cout.flush();
        Q.clear();

        U_vector.clear();
        for (node_set::iterator iter = U.begin(); iter != U.end(); iter++)
            U_vector.push_back(*iter);

        //ADDITION PHASE
        do {
            flag = 0;

            for (long i = 0; i < (long) U.size(); i++)
                //for(node_set::iterator iter=U.begin(); iter!=U.end();) //change way of iterating over U!!!! 
            {
                vj = U_vector[i];
                if (U.find(vj) != U.end()) {
                    //vj = *iter;
                    //compute delta_rho
                    x = get_no_neighs_in_set(vj, g, C);
                    y = g.get_node_out_degree(vj) - x;

                    I_new = I + x;
                    E_new = E - x + y;
                    n_new = n + 1;


                    nu_new = nu - 1;
                    for (adjacent_edges_iterator aeit = g.out_edges_begin(vj); aeit != g.out_edges_end(vj); aeit++) {
                        v = aeit->first;
                        if (C.find(v) == C.end() && U.find(v) == U.end())
                            nu_new = nu_new + 1;
                    }

                    if (y > 0)
                        nb_new = nb + 1;
                    for (adjacent_edges_iterator aeit = g.out_edges_begin(vj); aeit != g.out_edges_end(vj); aeit++) {
                        v = aeit->first;
                        if (C.find(v) != C.end() && get_no_neighs_in_set(v, g, U) == 1)
                            nb_new = nb_new - 1;
                    }

                    //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*(n_new-1)*E_new) - rho;
                    //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*E_new) - rho;
                    //delta_rho = ((pow(I_new,0.5))/(n_new*(n_new-1)/2))/(E_new/(nu_new*nb_new)) -rho;
                    //delta_rho = (I_new)/(E_new) - rho;
                    //delta_rho = E_new/(nu_new*nb_new);
                    delta_rho = (2 * I_new * nb_new) / (n_new * E_new) - rho;

                    if (delta_rho > 0) {
                        //cout<<"\nnode "<<g.get_node_label(vj)<<" agglomerated"<<endl;
                        C.insert(vj);
                        //iter++;
                        U.erase(vj);
                        for (adjacent_edges_iterator aeit = g.out_edges_begin(vj); aeit != g.out_edges_end(vj); aeit++) //REQUIRED!!!!
                        {
                            v = aeit->first;
                            if (C.find(v) == C.end() && U.find(v) == U.end())
                                U.insert(v);
                        }
                        Q.push_back(vj);
                        flag = 1;
                        rho += delta_rho;
                        I = I_new;
                        E = E_new;
                        n = n_new;
                        nb = nb_new;
                        nu = nu_new;
                    }
                }
                //else
                //    iter++;
            }

        } while (flag);

        //addition over

        //C_vector.clear();
        //for(node_set::iterator iter=C.begin(); iter!=C.end();iter++) //can be done incrementally inside the addition phase!
        //    C_vector.push_back(*iter);

        //DELETION PHASE
        //cout<<"entered deletion phase";
        //cout.flush();
        do {
            deleteQ.clear();

            for (node_set::iterator iter = C.begin(); iter != C.end();) {
                ui = *iter;
                //if(ui!=src)
                //{
                //compute delta_rho
                x = get_no_neighs_in_set(ui, g, C);
                y = g.get_node_out_degree(ui) - x;

                I_new = I - x;
                E_new = E + x - y;
                n_new = n - 1;

                nu_new = nu + 1;

                for (adjacent_edges_iterator aeit = g.out_edges_begin(ui); aeit != g.out_edges_end(ui); aeit++) {
                    v = aeit->first;
                    if (C.find(v) == C.end() && get_no_neighs_in_set(v, g, C) == 1)
                        nu_new = nu_new - 1;
                }

                if (y > 0)
                    nb_new = nb - 1;

                for (adjacent_edges_iterator aeit = g.out_edges_begin(ui); aeit != g.out_edges_end(ui); aeit++) {
                    v = aeit->first;
                    if (C.find(v) != C.end() && get_no_neighs_in_set(v, g, U) == 0)
                        nb_new = nb_new + 1;
                }

                //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*(n_new-1)*E_new) - rho;
                //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*E_new) - rho;
                //delta_rho = ((pow(I_new,0.5))/(n_new*(n_new-1)/2))/(E_new/(nu_new*nb_new)) -rho;
                //delta_rho = (I_new)/(E_new) - rho;
                delta_rho = (2 * I_new * nb_new) / (n_new * E_new) - rho;

                if (delta_rho > 0 && connected_on_removal(g, C, ui)) {
                    //cout<<"\nnode "<<g.get_node_label(ui)<<" deleted"<<endl;
                    iter++;
                    C.erase(ui);
                    U.insert(ui);
                    for (adjacent_edges_iterator aeit = g.out_edges_begin(ui); aeit != g.out_edges_end(ui); aeit++) {
                        v = aeit->first;
                        if (U.find(v) != U.end() && get_no_neighs_in_set(v, g, C) == 0)
                            U.erase(v);
                    }

                    deleteQ.push_back(ui);

                    if (found_in_vector(ui, Q, index))
                        delete_from_vector_at_position(ui, Q, index);

                    rho += delta_rho;
                    I = I_new;
                    E = E_new;
                    n = n_new;
                    nb = nb_new;
                    nu = nu_new;
                } else
                    iter++;
                //}
                //else
                //    iter++;
            }
        } while (!deleteQ.empty());

        //ADD REQUIRED VERTICES TO U
        //for(size_t k=0; k<Q.size(); k++)
        //    for(adjacent_edges_iterator aeit = g.out_edges_begin(Q[k]); aeit != g.out_edges_end(Q[k]); aeit++)
        //        if(C.find(aeit->first)==C.end() && U.find(aeit->first)==U.end())
        //            U.insert(aeit->first);

    } while (!Q.empty());

    if (rho > 1 && C.find(src) != C.end()) {
        //cout<<endl<<I<<" "<<E<<" "<<M<<endl;
        output = C;
        return 1; //comm found/exists
    } else {
        //cout<<"\nNO MODULE FOUND BY MY_ALGO FOR SOURCE VERTEX "<<g.get_node_label(src)<<endl;
        //output = C;
        output.insert(src);
        return 0; //comm not found/doesn't exist
    }

}
//-----------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------

bool CDLib::CZR(const graph& g, id_type src, node_set& output) {
    //DISCOVERY PHASE
    node_set C, U;
    C.insert(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first);

    double I = 0, E = g.get_node_out_degree(src), n = 1, nb = 1, L_in = 0, L_ex = E / nb, L = 0;
    double I_new, E_new, n_new, nb_new, L_in_new, L_ex_new, L_new;
    double L_at_start, L_at_end;
    double x, y;

    id_type ni, v, v_aggl;

    do {
        L_at_start = L;
        //double max = L_at_start;
        double max = -numeric_limits<double>::infinity(); //CZR being really unclear!!!!
        long no_max = 1;
        long cntr = -1;
        vector<id_type> V;

        for (node_set::iterator iter = U.begin(); iter != U.end(); iter++) {
            ni = *iter;
            //compute L_new
            x = get_no_neighs_in_set(ni, g, C);
            y = g.get_node_out_degree(ni) - x;

            I_new = I + x;
            E_new = E - x + y;
            n_new = n + 1;

            if (y > 0)
                nb_new = nb + 1;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(ni); aeit != g.out_edges_end(ni); aeit++) {
                v = aeit->first;
                //if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==1)
                if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 1)
                    nb_new = nb_new - 1;
            }

            L_new = (2 * I_new * nb_new) / (E_new * n_new);

            if (L_new > max) {
                V.push_back(ni);
                cntr += no_max;
                max = L_new;
                no_max = 1;
            } else if (L_new == max) {
                V.push_back(ni);
                no_max += 1;
            }
        }

        //cout<<V.size();

        if (!V.empty()) {
            //Selecting the node to be agglomerated v_aggl (breaking ties randomly if needed)
            //Choose a random number rand in [cntr,V.size()-1] such that v_aggl=V[rand] 
            RandomGenerator<long> rnd_gen(cntr, V.size() - 1);
            v_aggl = V[rnd_gen.next()];

            //compute and recompute
            x = get_no_neighs_in_set(v_aggl, g, C);
            y = g.get_node_out_degree(v_aggl) - x;

            I_new = I + x;
            E_new = E - x + y;
            n_new = n + 1;

            if (y > 0)
                nb_new = nb + 1;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(v_aggl); aeit != g.out_edges_end(v_aggl); aeit++) {
                v = aeit->first;
                //if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==1)
                if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 1)
                    nb_new = nb_new - 1;
            }

            L_new = (2 * I_new * nb_new) / (E_new * n_new);
            L_in_new = (2 * I_new) / n_new;
            L_ex_new = E_new / nb_new;

            //check extra conditions
            if ((L_in_new > L_in && L_ex_new < L_ex) || (L_in_new > L_in && L_ex_new > L_ex)) {
                //cout<<"node agglomerated";
                C.insert(v_aggl);
                U.erase(v_aggl);
                for (adjacent_edges_iterator aeit = g.out_edges_begin(v_aggl); aeit != g.out_edges_end(v_aggl); aeit++) {
                    v = aeit->first;
                    if (C.find(v) == C.end() && U.find(v) == U.end())
                        U.insert(v);
                }

                I = I_new, E = E_new, n = n_new, nb = nb_new;
                L = L_new, L_in = L_in_new, L_ex = L_ex_new;
            } else {
                //cout<<"node deleted from U";
                U.erase(v_aggl); //(DONE INTENTIONALLY----CAUTION!!!!)
            }

        } else //artificially added else condition (prevents infinite loop when U becomes empty) because CZR write pathetic pseudocode
            break;

        L_at_end = L;

    } while (L_at_end >= L_at_start); //again!! an example of CZR's irritating pseudocode! modification: >= instead of > to make the best sense

    //EXAMINATION PHASE
    for (node_set::iterator iter = C.begin(); iter != C.end();) {
        ni = *iter;
        //compute
        x = get_no_neighs_in_set(ni, g, C);
        y = g.get_node_out_degree(ni) - x;

        I_new = I - x;
        E_new = E + x - y;
        n_new = n - 1;

        if (y > 0)
            nb_new = nb - 1;

        for (adjacent_edges_iterator aeit = g.out_edges_begin(ni); aeit != g.out_edges_end(ni); aeit++) {
            v = aeit->first;
            //if (C.find(v) != C.end() && get_no_neighs_in_set(v, g, U) == 0)
            if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 0)
                nb_new = nb_new + 1;
        }

        L_new = (2 * I_new * nb_new) / (E_new * n_new);
        L_in_new = (2 * I_new) / n_new;
        L_ex_new = E_new / nb_new;

        if (!(L_in > L_in_new && L_ex < L_ex_new)) {
            iter++;
            //cout<<endl<<"node "<<g.get_node_label(ni)<<" deleted";
            C.erase(ni);
            U.insert(ni);
            for (adjacent_edges_iterator aeit = g.out_edges_begin(ni); aeit != g.out_edges_end(ni); aeit++) {
                v = aeit->first;
                //if (U.find(v) != U.end() && get_no_neighs_in_set(v, g, C) == 0) //DEPENDENT ON U----CAUTION!!!!
                if (C.find(v) == C.end() && get_no_neighs_in_set(v, g, C) == 0)
                    U.erase(v);
            }

            I = I_new, E = E_new, n = n_new, nb = nb_new;
            L = L_new, L_in = L_in_new, L_ex = L_ex_new;
        } else
            iter++;
    }

    //LAST PHASE
    if (C.find(src) != C.end()) {
        //cout<<endl<<I<<" "<<E<<" "<<M<<endl;
        output = C;
    } else {
        //cout<<"\nNO MODULE FOUND BY CZR FOR SOURCE VERTEX "<<g.get_node_label(src)<<endl;
        //output = C;
        output.insert(src);
    }

    if (output.size() > 1)
        return 1; //comm found/exists
    else
        return 0; //comm not found/doesn't exist

}
//-----------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------

bool CDLib::CZR_Beta(const graph& g, id_type src, node_set& output) {
    //DISCOVERY PHASE
    node_set C, U;
    C.insert(src);
    for (adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first);

    double I = 0, E = g.get_node_out_degree(src), n = 1, nb = 1, L_in = 0, L_ex = E / nb, L = 0;
    double I_new, E_new, n_new, nb_new, L_in_new, L_ex_new, L_new;
    double L_at_start, L_at_end;
    double x, y;

    id_type ni, v, v_aggl;

    do {
        L_at_start = L;
        //double max = L_at_start;
        double max = -numeric_limits<double>::infinity(); //CZR being really unclear!!!!
        long no_max = 1;
        long cntr = -1;
        vector<id_type> V;

        for (node_set::iterator iter = U.begin(); iter != U.end(); iter++) {
            ni = *iter;
            //compute L_new
            x = get_no_neighs_in_set(ni, g, C);
            y = g.get_node_out_degree(ni) - x;

            I_new = I + x;
            E_new = E - x + y;
            n_new = n + 1;

            if (y > 0)
                nb_new = nb + 1;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(ni); aeit != g.out_edges_end(ni); aeit++) {
                v = aeit->first;
                //if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==1)
                if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 1)
                    nb_new = nb_new - 1;
            }

            L_new = (2 * I_new * nb_new) / (E_new * n_new);

            if (L_new > max) {
                V.push_back(ni);
                cntr += no_max;
                max = L_new;
                no_max = 1;
            } else if (L_new == max) {
                V.push_back(ni);
                no_max += 1;
            }
        }

        //cout<<V.size();

        if (!V.empty()) {
            //Selecting the node to be agglomerated v_aggl (breaking ties randomly if needed)
            //Choose a random number rand in [cntr,V.size()-1] such that v_aggl=V[rand] 
            RandomGenerator<long> rnd_gen(cntr, V.size() - 1);
            v_aggl = V[rnd_gen.next()];

            //compute and recompute
            x = get_no_neighs_in_set(v_aggl, g, C);
            y = g.get_node_out_degree(v_aggl) - x;

            I_new = I + x;
            E_new = E - x + y;
            n_new = n + 1;

            if (y > 0)
                nb_new = nb + 1;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(v_aggl); aeit != g.out_edges_end(v_aggl); aeit++) {
                v = aeit->first;
                //if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==1)
                if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 1)
                    nb_new = nb_new - 1;
            }

            L_new = (2 * I_new * nb_new) / (E_new * n_new);
            L_in_new = (2 * I_new) / n_new;
            L_ex_new = E_new / nb_new;

            //check extra conditions
            if ((L_in_new > L_in && L_ex_new < L_ex) || (L_in_new > L_in && L_ex_new > L_ex)) {
                //cout<<"node agglomerated";
                C.insert(v_aggl);
                U.erase(v_aggl);
                for (adjacent_edges_iterator aeit = g.out_edges_begin(v_aggl); aeit != g.out_edges_end(v_aggl); aeit++) {
                    v = aeit->first;
                    if (C.find(v) == C.end() && U.find(v) == U.end())
                        U.insert(v);
                }

                I = I_new, E = E_new, n = n_new, nb = nb_new;
                L = L_new, L_in = L_in_new, L_ex = L_ex_new;
            } else {
                //cout<<"node deleted from U";
                U.erase(v_aggl); //(DONE INTENTIONALLY----CAUTION!!!!)
                //update E, nb, L_ex, L
                E = E - x;
                for (adjacent_edges_iterator aeit = g.out_edges_begin(v_aggl); aeit != g.out_edges_end(v_aggl); aeit++) {
                    v = aeit->first;
                    //if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==1)
                    if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 1)
                        nb = nb - 1;
                }
                L_ex = E / nb;
                L = L_in / L_ex;
            }

        } else //artificially added else condition (prevents infinite loop when U becomes empty) because CZR write pathetic pseudocode
            break;

        L_at_end = L;

    } while (L_at_end > L_at_start); //again!! an example of CZR's irritating pseudocode! modification: >= instead of > to make the best sense

    //EXAMINATION PHASE
    for (node_set::iterator iter = C.begin(); iter != C.end();) {
        ni = *iter;
        //compute
        x = get_no_neighs_in_set(ni, g, C);
        y = g.get_node_out_degree(ni) - x;

        I_new = I - x;
        E_new = E + x - y;
        n_new = n - 1;

        if (y > 0)
            nb_new = nb - 1;

        for (adjacent_edges_iterator aeit = g.out_edges_begin(ni); aeit != g.out_edges_end(ni); aeit++) {
            v = aeit->first;
            //if (C.find(v) != C.end() && get_no_neighs_in_set(v, g, U) == 0)
            if (C.find(v) != C.end() && (g.get_node_out_degree(v) - get_no_neighs_in_set(v, g, C)) == 0)
                nb_new = nb_new + 1;
        }

        L_new = (2 * I_new * nb_new) / (E_new * n_new);
        L_in_new = (2 * I_new) / n_new;
        L_ex_new = E_new / nb_new;

        if (!(L_in > L_in_new && L_ex < L_ex_new)) {
            iter++;
            //cout<<endl<<"node "<<g.get_node_label(ni)<<" deleted";
            C.erase(ni);
            U.insert(ni);
            for (adjacent_edges_iterator aeit = g.out_edges_begin(ni); aeit != g.out_edges_end(ni); aeit++) {
                v = aeit->first;
                //if (U.find(v) != U.end() && get_no_neighs_in_set(v, g, C) == 0) //DEPENDENT ON U----CAUTION!!!!
                if (C.find(v) == C.end() && get_no_neighs_in_set(v, g, C) == 0)
                    U.erase(v);
            }

            I = I_new, E = E_new, n = n_new, nb = nb_new;
            L = L_new, L_in = L_in_new, L_ex = L_ex_new;
        } else
            iter++;
    }

    //LAST PHASE
    if (C.find(src) != C.end()) {
        //cout<<endl<<I<<" "<<E<<" "<<M<<endl;
        output = C;
    } else {
        //cout<<"\nNO MODULE FOUND BY CZR FOR SOURCE VERTEX "<<g.get_node_label(src)<<endl;
        //output = C;
        output.insert(src);
    }

    if (output.size() > 1)
        return 1; //comm found/exists
    else
        return 0; //comm not found/doesn't exist

}

struct edge_radicchi {
    id_type from_id;
    id_type to_id;
    double ecc;

    edge_radicchi() : from_id(0), to_id(0), ecc(0) {
    }

    edge_radicchi(id_type from, id_type to, double ec) : from_id(from), to_id(to), ecc(ec) {
    }
};

void CDLib::girvan_newman_2002(const graph& g, dendrogram& dendro) {
    if (is_connected_weakly(g)) {
        graph gc(g);
        id_type component_counter = 1;
        while (gc.get_num_edges()) {
            vector<double> bc;
            betweeness_centralities(gc, bc);
            double max_edge_bet = 0;
            id_type from_id = 0, to_id = 0;
            for (id_type i = 0; i < gc.get_num_nodes(); i++) {
                for (adjacent_edges_iterator aeit = gc.out_edges_begin(i); aeit != gc.out_edges_end(i); aeit++) {
                    double edge_bet = bc[i] + bc[aeit->first];
                    if (edge_bet >= max_edge_bet) {
                        max_edge_bet = edge_bet;
                        from_id = i;
                        to_id = aeit->first;
                    }
                }
            }
            gc.remove_edge(from_id, to_id);
            vector<node_set> components;
            id_type num_comps = get_weakly_connected_components(gc, components);
            if (num_comps - component_counter) {
                dendro.push_back(components);
                component_counter = num_comps;
            }
        }
    }
}

bool ec_comp(const edge_radicchi& lhs, const edge_radicchi& rhs) {
    return lhs.ecc < rhs.ecc;
}

void CDLib::radicchi_et_al_2004(const graph& g, dendrogram& dendro) {
    vector<node_set> components;
    if (get_weakly_connected_components(g, components) == 1) {
        graph gc(g);
        vector<edge_radicchi> ec(g.get_num_edges(), edge_radicchi());
        id_type count = 0;
        for (id_type i = 0; i < gc.get_num_nodes(); i++)
            for (adjacent_edges_iterator aeit = gc.out_edges_begin(i); aeit != gc.out_edges_end(i); aeit++)
                ec[count++] = edge_radicchi(i, aeit->first, edge_clustering_coefficient(gc, i, aeit->first));
        sort(ec.begin(), ec.end(), ec_comp);
        for (id_type i = 0; i < ec.size(); i++) {
            gc.remove_edge(ec[i].from_id, ec[i].to_id);
            vector<node_set> components;
            id_type num_comps = get_weakly_connected_components(gc, components);
            if (num_comps > 2)dendro.push_back(components);
        }
    }
}

id_type CDLib::label_propagation_run(const graph&g, vector<id_type>& communities, id_type max_iters, lp_raghavan_2007* algoman) {
    id_type num_iters = 0;
    for (num_iters = 0; num_iters < max_iters; num_iters++)
        if (algoman->do_iteration(g, num_iters)) break;
    algoman->finalize(g, num_iters, communities);
    return num_iters;
}

bool lp_raghavan_2007::is_node_lplabel_fixed(id_type node_id) const {
    if (node_id > ids.size()) return false;
    else return (nodes_with_fixed_lplabels.find(node_id) != nodes_with_fixed_lplabels.end());
}

bool lp_raghavan_2007::fix_node_lplabel(id_type node_id) {
    if (node_id > ids.size()) return false;
    nodes_with_fixed_lplabels.insert(node_id);
    return true;
}

bool lp_raghavan_2007::free_node_lplabel(id_type node_id) {
    if (node_id > ids.size()) return false;
    nodes_with_fixed_lplabels.erase(node_id);
    return true;
}

bool lp_raghavan_2007::set_node_lplabel(id_type node_id, id_type lplabel) {
    if (node_id > lplabels.size() && lplabel < lplabels.size() && is_node_lplabel_fixed(node_id)) return false;
    lplabels[node_id] = lplabel;
    return true;
}

id_type lp_raghavan_2007::get_node_lplabel(id_type node_id) const {
    if (node_id > lplabels.size()) return lplabels.size();
    return lplabels[node_id];
}

void lp_raghavan_2007::post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels) {
}

void lp_raghavan_2007::post_iteration(const graph& g, id_type num_iters) {
}

lp_raghavan_2007::lp_raghavan_2007(const graph& g, bool synchronous_val) {
    synchronous = synchronous_val;
    ids.assign(g.get_num_nodes(), 0);
    positions.assign(g.get_num_nodes(), 0);
    lplabels.assign(g.get_num_nodes(), 0);
    if (synchronous)lpnextlabels.assign(g.get_num_nodes(), 0);
    for (id_type i = 0; i < g.get_num_nodes(); i++) {
        ids[i] = i;
        positions[i] = i;
        lplabels[i] = i;
        if (synchronous)lpnextlabels[i] = i;
    }
}

void lp_raghavan_2007::finalize(const graph& g, id_type num_iters, vector<id_type>& communities) {
    communities.assign(g.get_num_nodes(), 0);
    copy(lplabels.begin(), lplabels.end(), communities.begin());
}

bool lf_comp(const pair<id_type, double>& lhs, const pair<id_type, double>& rhs) {
    return lhs.second < rhs.second;
}

void lp_raghavan_2007::get_max_lplabels(const graph& g, id_type current_node, max_lplabel_container& max_labels) {
    lplabel_fitness_container lplabel_fitness_values;
    for (adjacent_edges_iterator aeit = g.in_edges_begin(current_node); aeit != g.in_edges_end(current_node); aeit++) {
        double fitness_value = get_node_score(g, aeit->first) * get_label_fitness_for_edge(g, current_node, aeit->first, aeit->second);
        pair < lplabel_fitness_container::iterator, bool> ret = lplabel_fitness_values.insert(make_pair(lplabels[aeit->first], fitness_value));
        if (!ret.second) ret.first->second += fitness_value;
    }
    lplabel_fitness_container::iterator mit = max_element(lplabel_fitness_values.begin(), lplabel_fitness_values.end(), lf_comp);
    double max_val = mit->second;
    while (mit != lplabel_fitness_values.end() && mit->second == max_val) {
        max_labels.push_back(mit->first);
        lplabel_fitness_values.erase(mit->first);
        mit = max_element(lplabel_fitness_values.begin(), lplabel_fitness_values.end(), lf_comp);
    }
}

bool lp_raghavan_2007::check_lplabels(const graph& g) {
    for (id_type i = 0; i < g.get_num_nodes(); i++) {
        vector<id_type> max_labels;
        get_max_lplabels(g, i, max_labels);
        if (find(max_labels.begin(), max_labels.end(), lplabels[i]) == max_labels.end())
            return false;
    }
    return true;
}

id_type lp_raghavan_2007::new_label_break_ties_randomly(const graph& g, id_type current_node, max_lplabel_container& max_labels) {
    get_max_lplabels(g, current_node, max_labels);
    if (max_labels.empty()) return lplabels[current_node];
    CDLib::RandomGenerator<id_type> rnd_gen(0, max_labels.size() - 1);
    return max_labels[rnd_gen.next()];
}

bool lp_raghavan_2007::do_iteration(const graph& g, id_type num_iters) {
    for (vector<id_type>::iterator it = ids.begin(); it != ids.end(); it++) {
        reorder(g, num_iters);
        max_lplabel_container max_labels;
        id_type new_label = new_label_break_ties_randomly(g, *it, max_labels);
        post_node_assign(g, *it, new_label, num_iters, max_labels);
        if (synchronous) lpnextlabels[*it] = new_label;
        else lplabels[*it] = new_label;
    }
    if (synchronous)
        for (id_type i = 0; i < g.get_num_nodes(); i++)
            lplabels[i] = lpnextlabels[i];
    post_iteration(g, num_iters);
    return check_lplabels(g);
}

void lp_raghavan_2007::reorder(const graph& g, id_type num_iters) {
    random_shuffle(ids.begin(), ids.end());
    for (vector<id_type>::iterator it = ids.begin(); it != ids.end(); it++)
        positions[*it] = it - ids.begin();
}

void lp_offensive_lpa::post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels) {
    lp_dyn_hop_att::post_node_assign(g, current, new_label, num_iters, max_labels);
    if (lplabels[current] != new_label)
        probabilities[current] = compute_probability_for_id(g, current, new_label, max_labels);

}

double lp_offensive_lpa::compute_probability_for_id(const graph& g, id_type current, id_type new_label, max_lplabel_container& max_labels) {
    double prod_num = 0;
    for (adjacent_edges_iterator aeit = g.in_edges_begin(current); aeit != g.in_edges_end(current); aeit++) {
        if (lplabels[aeit->first] == new_label)prod_num += probabilities[aeit->first];
    }
    return prod_num / g.get_node_in_weight(current);
}

double lp_offensive_lpa::get_edge_weight_function(const graph&g, id_type from_id, id_type to_id, double edge_weight) {
    return edge_weight;
}

double lp_offensive_lpa::get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
    return get_node_score(g, to_id)*(1 - probabilities[to_id]) * get_edge_weight_function(g, from_id, to_id, edge_weight);
}

double lp_offensive_lpa::get_node_probability(id_type id) const {
    if (id > probabilities.size()) return 0.0;
    else return probabilities[id];
}

lp_offensive_lpa::lp_offensive_lpa(const graph& g, bool synchronous_val, double hop_att_max_val) : lp_dyn_hop_att(g, synchronous_val, hop_att_max_val) {
    probabilities.assign(g.get_num_nodes(), (double) 1 / (double) g.get_num_nodes());
}

double lp_defensive_lpa::get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
    return get_node_score(g, to_id) * probabilities[to_id] * get_edge_weight_function(g, from_id, to_id, edge_weight);
}

double lp_defensive_lpa::compute_probability_for_id(const graph& g, id_type current, id_type new_label, max_lplabel_container& max_labels) {
    double prob_denom = 0, prod_num = 0;
    for (adjacent_edges_iterator aeit = g.in_edges_begin(current); aeit != g.in_edges_end(current); aeit++) {
        if (lplabels[aeit->first] == new_label) {
            prod_num += probabilities[aeit->first];
            prob_denom += aeit->second;
        }
    }
    return prod_num / prob_denom;
}

lp_defensive_lpa::lp_defensive_lpa(const graph& g, bool synchronous_val, double hop_att_max_val) : lp_offensive_lpa(g, synchronous_val, hop_att_max_val) {
}

id_type CDLib::label_propagation_raghavan_2007(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous) {
    lp_raghavan_2007 algoman(g, synchronous);
    return label_propagation_run(g, labels, max_iters, &algoman);
}

id_type CDLib::label_propagation_leung_2009(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att) {
    lp_leung_2009 algoman(g, synchronous, hop_att);
    return label_propagation_run(g, labels, max_iters, &algoman);
}

id_type CDLib::label_propagation_dyn_hop_2010(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att_max) {
    lp_dyn_hop_att algoman(g, synchronous, hop_att_max);
    return label_propagation_run(g, labels, max_iters, &algoman);
}

id_type CDLib::label_propagation_olpa_2010(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att_max) {
    lp_offensive_lpa algoman(g, synchronous, hop_att_max);
    return label_propagation_run(g, labels, max_iters, &algoman);
}

id_type CDLib::label_propagation_dlpa_2010(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att_max) {
    lp_defensive_lpa algoman(g, synchronous, hop_att_max);
    return label_propagation_run(g, labels, max_iters, &algoman);
}

id_type CDLib::label_propagation_track_changes_2012(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous) {
    lp_track_changes algoman(g, synchronous);
    return label_propagation_run(g, labels, max_iters, &algoman);
}

struct dynamic_lp_temp {
    vector<id_type> labels_to_pass;
    vector<id_type> label_changes_to_pass;
};

void first_snapshot_init(string filepath, bool directed, bool weighted, dynamic_lp_output & output) {
    output.lp_iters.push_back(0);
    output.lp_times.push_back(0.0);
    output.book_times.push_back(0.0);
    output.lp_labels.push_back(vector<id_type > ());
    output.graphs.push_back(graph(directed, weighted));
    read_edgelist(output.graphs[output.graphs.size() - 1], filepath, directed, weighted);
    output.lp_num_label_changes.push_back(vector<id_type > (output.graphs[output.graphs.size() - 1].get_num_nodes(), 0));
}

void subsequent_snapshot_init(string filepath, bool directed, bool weighted, dynamic_lp_output& output, dynamic_lp_temp& temp, bool copy_labels, timer_rt& op_timer) {
    first_snapshot_init(filepath, directed, weighted, output);
    temp.labels_to_pass.assign(output.graphs[output.graphs.size() - 1].get_num_nodes(), 0);
    temp.label_changes_to_pass.assign(output.graphs[output.graphs.size() - 1].get_num_nodes(), 0);
    op_timer.start_clock();
    node_set seen_nodes;
    id_type max_label_value = 0;
    for (id_type i = 0; i < output.graphs[output.graphs.size() - 2].get_num_nodes(); i++) {
        id_type current_node_id = output.graphs[output.graphs.size() - 1].get_node_id(output.graphs[output.graphs.size() - 2].get_node_label(i));
        if (current_node_id < output.graphs[output.graphs.size() - 1].get_num_nodes()) //Node has been found
        {
            temp.labels_to_pass[current_node_id] = output.lp_labels[output.lp_labels.size() - 2][i];
            max_label_value = (output.lp_labels[output.lp_labels.size() - 2][i] > max_label_value) ? output.lp_labels[output.lp_labels.size() - 2][i] : max_label_value;
            temp.label_changes_to_pass[current_node_id] = output.lp_num_label_changes[output.lp_num_label_changes.size() - 2][i];
            seen_nodes.insert(current_node_id);
        }
    }
    max_label_value++;
    for (id_type i = 0; i < output.graphs[output.graphs.size() - 1].get_num_nodes(); i++) {
        if (seen_nodes.find(i) == seen_nodes.end()) { //Unseen Node
            temp.labels_to_pass[i] = max_label_value + i;
            temp.label_changes_to_pass[i] = 0;
        }
    }
    op_timer.stop_clock();
}

void first_run(timer_rt& op_timer, const string& filepath, bool directed, bool weighted, dynamic_lp_output& output) {
    first_snapshot_init(filepath, directed, weighted, output);
    op_timer.start_clock();
    lp_dyn_hop_att algoman(output.graphs[output.graphs.size() - 1], 0, 0.0);
    output.lp_iters[output.lp_iters.size() - 1] = label_propagation_run(output.graphs[output.graphs.size() - 1], output.lp_labels[output.lp_labels.size() - 1], 10000, &algoman);
    op_timer.stop_clock();
    output.lp_times[output.lp_times.size() - 1] = op_timer.run_time();
    for (id_type i = 0; i < output.graphs[output.graphs.size() - 1].get_num_nodes(); i++)output.lp_num_label_changes[output.lp_num_label_changes.size() - 1][i] = algoman.get_node_num_lplabel_changes(i);
}

void subsequent_run(timer_rt& op_timer, const string& filepath, bool directed, bool weighted, dynamic_lp_input& input, dynamic_lp_output& output) {
    dynamic_lp_temp temp;
    subsequent_snapshot_init(filepath, directed, weighted, output, temp, input.copy_labels, op_timer);
    op_timer.start_clock();
    evolutionary_label_propagation algoman(output.graphs[output.graphs.size() - 1], 0, input.activate, input.copy_labels, input.alpha, output.lp_iters[output.lp_iters.size() - 2], temp.label_changes_to_pass, temp.labels_to_pass);
    output.lp_iters[output.lp_iters.size() - 1] = label_propagation_run(output.graphs[output.graphs.size() - 1], output.lp_labels[output.lp_labels.size() - 1], 10000, &algoman);
    op_timer.stop_clock();
    output.lp_times[output.lp_times.size() - 1] = op_timer.run_time();
    for (id_type i = 0; i < output.graphs[output.graphs.size() - 1].get_num_nodes(); i++)output.lp_num_label_changes[output.lp_num_label_changes.size() - 1][i] = algoman.get_node_num_lplabel_changes(i);
}

double CDLib::evolutionary_label_propagation_edgelists(const string& snapshot_filepath, bool directed, bool weighted, dynamic_lp_input& input, dynamic_lp_output& output) {
    ifstream ifs;
    ifs.open(snapshot_filepath.c_str());
    timer_rt op_timer;
    if (ifs.is_open()) {
        string filepath;
        ifs >> filepath;
        first_run(op_timer, filepath, directed, weighted, output);
        while (!ifs.eof()) {
            ifs >> filepath;
            subsequent_run(op_timer, filepath, directed, weighted, input, output);
        }
    }
    return op_timer.total_time();
}

void pre_run(ifstream& ifs, bool directed, bool weighted, dynamic_lp_output& output) {
    output.lp_iters.push_back(0);
    output.lp_times.push_back(0.0);
    output.book_times.push_back(0.0);
    output.lp_labels.push_back(vector<id_type > ());
    output.graphs.push_back(graph(directed, weighted));
    string filepath, outfilepath;
    ifs >> filepath >> outfilepath;
    read_edgelist(output.graphs[output.graphs.size() - 1], filepath, directed, weighted);
}

double prepare_labels_to_pass(dynamic_lp_output& output, vector<id_type>::const_iterator begin_it, vector<id_type>& temp_changes) {
    timer_rt lp_timer;
    lp_timer.start_clock();
    temp_changes.assign(output.graphs[output.graphs.size() - 1].get_num_nodes(), 0);
    for (id_type i = 0; i < output.graphs[output.graphs.size() - 2].get_num_nodes(); i++) {
        id_type current_node_id = output.graphs[output.graphs.size() - 1].get_node_id(output.graphs[output.graphs.size() - 2].get_node_label(i));
        if (current_node_id < output.graphs[output.graphs.size() - 1].get_num_nodes())
            temp_changes[current_node_id] = *(begin_it + i);
    }
    lp_timer.stop_clock();
    return lp_timer.run_time();
}

void new_subsequent_run(ifstream& ifs, bool directed, bool weighted, double alpha, dynamic_lp_output& output, vector<id_type>& temp_changes) {
    timer_rt lp_timer;
    lp_timer.start_clock();
    evol_label_prop_new algoman(output.graphs[output.graphs.size() - 1], 0, alpha, output.lp_iters[output.lp_iters.size() - 2], temp_changes.begin());
    output.lp_iters[output.lp_iters.size() - 1] = label_propagation_run(output.graphs[output.graphs.size() - 1], output.lp_labels[output.lp_labels.size() - 1], 10000, &algoman);
    lp_timer.stop_clock();
    output.lp_times[output.lp_times.size() - 1] = lp_timer.run_time();
    pre_run(ifs, directed, weighted, output);
    output.book_times[output.book_times.size() - 1] = prepare_labels_to_pass(output, algoman.label_changes_begin(), temp_changes);
}

void CDLib::new_evolutionary_label_propagation_edgelists(const string& snapshot_filepath, bool directed, bool weighted, double alpha, bool basic, dynamic_lp_output& output) {
    ifstream ifs;
    ifs.open(snapshot_filepath.c_str());

    if (ifs.is_open()) {
        id_type num_snapshots = 0;
        vector<id_type> temp_changes;
        while (!ifs.eof()) {
            pre_run(ifs, directed, weighted, output);
            timer_rt lp_timer;
            if (!basic) {
                if (!num_snapshots) {
                    lp_timer.start_clock();
                    lp_track_changes algoman(output.graphs[output.graphs.size() - 1], 0);
                    output.lp_iters[output.lp_iters.size() - 1] = label_propagation_run(output.graphs[output.graphs.size() - 1], output.lp_labels[output.lp_labels.size() - 1], 10000, &algoman);
                    lp_timer.stop_clock();
                    output.lp_times[output.lp_times.size() - 1] = lp_timer.run_time();
                    pre_run(ifs, directed, weighted, output);
                    temp_changes.assign(output.graphs[output.graphs.size() - 1].get_num_nodes(), 0);
                    output.book_times[output.book_times.size() - 1] = prepare_labels_to_pass(output, algoman.label_changes_begin(), temp_changes);
                } else new_subsequent_run(ifs, directed, weighted, alpha, output, temp_changes);
            } else {
                lp_timer.start_clock();
                lp_raghavan_2007 algoman(output.graphs[output.graphs.size() - 1], 0);
                output.lp_iters[output.lp_iters.size() - 1] = label_propagation_run(output.graphs[output.graphs.size() - 1], output.lp_labels[output.lp_labels.size() - 1], 10000, &algoman);
                lp_timer.stop_clock();
                output.lp_times[output.lp_times.size() - 1] = lp_timer.run_time();
            }
        }
    }
}

id_type bgll_find_best_community_for_node(const graph& curr_graph, const vector<id_type>& labels, const bgll_objective& book, id_type vertex) {
    double max_node_gain = -numeric_limits<double>::infinity();
    id_type node_max_gain_comm_index = labels[vertex];
    for (adjacent_edges_iterator aeit = curr_graph.out_edges_begin(vertex); aeit != curr_graph.out_edges_end(vertex); aeit++) {
        if (labels[vertex] != labels[aeit->first]) {
            double curr_node_gain = book.compute_gain(curr_graph, labels, vertex, labels[aeit->first]);
            if (curr_node_gain > max_node_gain) {
                max_node_gain = curr_node_gain;
                node_max_gain_comm_index = labels[aeit->first];
            }
        }
    }
    return node_max_gain_comm_index;
}

id_type bgll_vertex_mover_single_pass(const graph& curr_graph, vector<id_type>& labels, bgll_objective& book) {
    id_type num_nodes_moved = 0;
    for (id_type i = 0; i < curr_graph.get_num_nodes(); i++) {
        book.detach_node(curr_graph, labels, i); //Detach node from its community
        id_type new_node_community = bgll_find_best_community_for_node(curr_graph, labels, book, i);
        if (new_node_community != labels[i]) {
            num_nodes_moved++;
            //Insert node into new community
            book.attach_node(curr_graph, labels, i, new_node_community);
            labels[i] = new_node_community;
        }
    }
    return num_nodes_moved;
}

double CDLib::bgll_vertex_mover_optimizer(const graph& g, vector<id_type>& labels, bgll_objective& book) {
    double start_obj_val = book.objval(g, labels), curr_objval = -numeric_limits<double>::infinity();
    id_type num_nodes_moved = 0;
    do {
        vector<id_type> temp_comms(labels);
        num_nodes_moved += bgll_vertex_mover_single_pass(g, temp_comms, book);
        curr_objval = book.objval(g, temp_comms);
        //Rollback communities and exit if modularity decreases
        if (curr_objval < start_obj_val) return start_obj_val;
        //Commit changes, book should have already been updated
        copy(temp_comms.begin(), temp_comms.end(), labels.begin());
    } while (num_nodes_moved && curr_objval < start_obj_val);
    return curr_objval;
}

void bgll_recover_communities(vector<id_type>& labels, vector< vector<id_type> >& hier_comms) {
    vector<node_set> last_comms;
    convert_labels_to_communities(hier_comms[hier_comms.size() - 1], last_comms);
    reindex_communities(labels);
    hier_comms.push_back(vector<id_type > (hier_comms[hier_comms.size() - 1].size(), 0));
    for (id_type i = 0; i < labels.size(); i++) {
        for (node_set::iterator nit = last_comms[i].begin(); nit != last_comms[i].end(); nit++) {
            hier_comms[hier_comms.size() - 1][*nit] = labels[i];
        }
    }
}

void bgll_collapse_nodes(graph& orig_graph, vector<id_type>& labels) {
    vector<unordered_map<id_type, double> > edges;
    for (id_type i = 0; i < labels.size(); i++)
        while (edges.size() <= labels[i]) edges.push_back(unordered_map<id_type, double>());
    for (id_type i = 0; i < orig_graph.get_num_nodes(); i++) {
        id_type node_comm = labels[i];
        for (adjacent_edges_iterator aeit = orig_graph.out_edges_begin(i); aeit != orig_graph.out_edges_end(i); aeit++) {
            id_type neigh_comm = labels[aeit->first];
            pair < unordered_map<id_type, double>::iterator, bool> ret = edges[node_comm].insert(make_pair(neigh_comm, 0));
            ret.first->second += aeit->second;
        }
    }
    orig_graph.clear();
    for (id_type i = 0; i < edges.size(); i++) orig_graph.add_node();
    for (id_type i = 0; i < edges.size(); i++)
        for (unordered_map<id_type, double>::iterator it = edges[i].begin(); it != edges[i].end(); it++)
            orig_graph.add_edge(i, it->first, it->second);
    labels.assign(orig_graph.get_num_nodes(), 0);
    for (id_type i = 0; i < orig_graph.get_num_nodes(); i++)
        labels[i] = i;
}

void CDLib::cda_bgll_generic(const graph&g, const vector<id_type>& init_comms, vector< vector<id_type> >& hier_comms, bgll_objective& book) {
    hier_comms.clear();
    //Initialize the initial_partition
    vector<id_type> curr_comms(g.get_num_nodes(), 0);
    if (init_comms.size() == g.get_num_nodes()) copy(init_comms.begin(), init_comms.end(), curr_comms.begin());
    else for (id_type i = 0; i < g.get_num_nodes(); i++) curr_comms[i] = i;
    hier_comms.push_back(curr_comms);
    //Create a copy of the graph
    graph curr_graph(g);
    curr_graph.convert_to_weighted();
    double orig_obj_val = book.objval(curr_graph, curr_comms);
    //The Main Loop
    while (1) {
        id_type curr_graph_size = curr_graph.get_num_nodes();
        //Initialize the book for the new level
        book.init(g, hier_comms, curr_graph, curr_comms);
        //Run the vertex mover optimization
        double new_obj_val = bgll_vertex_mover_optimizer(curr_graph, curr_comms, book);
        //Reindex the community memberships
        bgll_recover_communities(curr_comms, hier_comms);
        //Shrink the graph
        bgll_collapse_nodes(curr_graph, curr_comms);
        if (curr_graph.get_num_nodes() == curr_graph_size || new_obj_val < orig_obj_val) break;
    }
}

class bgll_modularity : public bgll_objective_w_internal {
private:
    unordered_map<id_type, double> total_weight;
public:

    bgll_modularity() : bgll_objective_w_internal() {
    }

    bgll_modularity(double resolution_p) : bgll_objective_w_internal(resolution_p) {
    }

    void init(const graph& orig_graph, const vector<vector<id_type> >& hiercomms, const graph& curr_graph, const vector<id_type>& curr_comms) {
        init_internal_params();
        total_weight.clear();
        for (id_type i = 0; i < curr_graph.get_num_nodes(); i++) {
            fill_internal_edges_for_vertex(curr_graph, curr_comms, i);
            map_insert_and_increment<id_type, double>(total_weight, curr_comms[i], curr_graph.get_node_out_weight(i));
        }
    }

    void detach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) {
        update_internal_edges_on_vertex_detach(curr_graph, labels, vertex);
        map_find_and_modify<id_type, double>(total_weight, labels[vertex], -curr_graph.get_node_out_weight(vertex));
    }

    void attach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) {
        if (labels[vertex] != dst_comm) {
            update_internal_edges_on_vertex_attach(dst_comm);
            map_find_and_modify<id_type, double>(total_weight, dst_comm, curr_graph.get_node_out_weight(vertex));
        }
    }

    double objval(const graph& curr_graph, const vector<id_type>& curr_comms) const {
        double modularity = 0;
        for (unordered_map<id_type, double>::const_iterator it = total_weight.begin(); it != total_weight.end(); it++) {
            modularity += objval_internal_term_running(curr_graph, it->first);
            modularity -= (pow(it->second / (2 * curr_graph.get_total_weight()), 2));
        }
        return modularity;
    }

    double compute_gain(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) const {
        double gain = 0;
        gain += gain_internal_term(dst_comm);
        unordered_map<id_type, double>::const_iterator it = total_weight.find(dst_comm);
        if (it != total_weight.end()) gain -= ((curr_graph.get_node_out_weight(vertex) * it->second) / curr_graph.get_total_weight());
        return gain;
    }

};

void CDLib::cda_bgll_modularity(const graph& g, const vector<id_type>& init_comms, vector< vector<id_type> >& hier_comms, double resolution_param) {
    bgll_modularity book(resolution_param);
    cda_bgll_generic(g, init_comms, hier_comms, static_cast<bgll_objective&> (book));
}


/* 
 * File:   local_community.cpp
 * Author: prashant
 * 
 * Created on June 16, 2012, 12:09 PM
 */

#include "local_community.h"
#include "utility.h"
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

double get_neighs_in_set(id_type v, const graph& g, node_set& B, node_set& X)
{
    X.clear();
    double retval = 0;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
    {
        if(B.find(aeit->first) != B.end())
        {
            X.insert(aeit->first);
            retval += aeit->second;
        }
    }
    return retval;
}

double get_neighs_not_in_set(id_type v, const graph& g, node_set& C, node_set& X)
{
    X.clear();
    double retval = 0;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
    {
        if(C.find(aeit->first) == C.end())   
        {
            X.insert(aeit->first);
            retval += aeit->second;
        }
    }
    return retval;
}


bool membership(id_type vk, id_type v, const graph& g, node_set& U)
{
    for(adjacent_edges_iterator aeit = g.out_edges_begin(vk); aeit != g.out_edges_end(vk); aeit++)
        if((aeit->first)!=v && U.find(aeit->first)!=U.end())
            return 1;
    return 0;
}

double no_BI_edges(id_type v, const graph& g, node_set& C, node_set& B)
{
    double count=0;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
        if(C.find(aeit->first)!=C.end() && B.find(aeit->first)==B.end())
            count+=aeit->second;
    return count;
}

double no_BY_edges(id_type v, const graph& g, node_set& Y)
{
    double count=0;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
        if(Y.find(aeit->first)!=Y.end())
            count+=aeit->second;
    return count;
}



//-----------------------------------------------------------------------------------------------------------------------------------------
void CDLib::local_community_clauset(const graph& g, id_type src, size_t k, vector< pair<id_type,double> >& output)
{
    node_set C,U,B;
    C.insert(src);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first); 
    B.insert(src);
    
    double R=0;
    double T=g.get_node_out_weight(src);
    
    //Update output
    output.push_back(make_pair(src,R));
    
    while (C.size() < k)
    {
        double max=-numeric_limits<double>::infinity();
        long cntr=-1; //Doubt
        long no_max=1;
        vector<id_type> V;
        
        unordered_map<id_type,node_set> temp_Y; //To store the necessary v and Y pairs 
        unordered_map<id_type,double> temp_T; //To store the necessary v and T pairs        
        
        node_set X,Y;
        double x,y,z;
        
        for(node_set::iterator iter=U.begin();iter!=U.end();iter++)
        {
            //compute delta_Rj
            double delta;

            id_type v=*iter; //vertex vj
            
            x = get_neighs_in_set(v,g,B,X); //X = set of all neighbors of vj in B   //x = No. of edges in T that terminate at vj
            
            y = g.get_node_out_weight(v)-x; //No. of edges that will be added to T by the agglomeration of vj
            
            z=0; //No. of edges that will be removed from T by the agglomeration of vj
            for(node_set::iterator iter1=X.begin();iter1!=X.end();iter1++)
            {
                id_type vk=*iter1;
                if(membership(vk,v,g,U)==0) //vk has no neighbor in U other than vj, i.e., vk no longer belongs to B 
                {
                    Y.insert(vk);
                    z+=no_BI_edges(vk,g,C,B)+no_BY_edges(vk,g,Y);
                    if(y==0)
                        z=z+g.get_edge_weight(vk,v);
                }     
            }
            
            delta=(x-(R*y)-z*(1-R))/(T-z+y);
            
            if(delta>max)
            {
                V.push_back(v);
                cntr+=no_max;
                max=delta;
                no_max=1;
                //store v and corresponding Y
                temp_Y.insert(make_pair(v,Y));
                //store v and corresponding new value of T
                temp_T.insert(make_pair(v,T-z+y));
            }
            else if(delta==max)
            {
                V.push_back(v);
                no_max+=1;
                //store v and corresponding Y
                temp_Y.insert(make_pair(v,Y));
                //store v and corresponding new value of T
                temp_T.insert(make_pair(v,T-z+y));
            }
        }
        
        //Selecting the node to be agglomerated v_aggl (breaking ties randomly if needed)
        //Choose a random number rand in [cntr,V.size()-1] such that v_aggl=V[rand] 
        RandomGenerator<long> rnd_gen(cntr,V.size()-1); //Doubt
        id_type v_aggl=V[rnd_gen.next()]; //Doubt
        
        //Update B
        node_set Y_aggl; //corresponding Y of v_aggl
        /*vector<pair<id_type,node_set>>::iterator iter2; //Doubt
        for(iter2=temp_Y.begin();iter2<temp_Y.end();iter2++) //Doubt
            if(iter2->first==v_aggl)
            {
                Y_aggl=iter2->second;
                break;
            }*/
        Y_aggl=(temp_Y.find(v_aggl))->second;
        
        x = get_neighs_in_set(v_aggl,g,B,X);
        y =g.get_node_out_weight(v_aggl)-x;
        
        for(node_set::iterator iter3=Y_aggl.begin();iter3!=Y_aggl.end();iter3++)
            B.erase(*iter3);
        
        if(y>0) //Doubt (complexity)
            B.insert(v_aggl);
        
        //Update U
        U.erase(v_aggl);
        get_neighs_not_in_set(v_aggl,g,C,X); //X = all neighbors of v_aggl not in C
        for(node_set::iterator iter4=X.begin();iter4!=X.end();iter4++)
            U.insert(*iter4);
        
        //Update C
        C.insert(v_aggl);
        
        //Update R
        R=R+max;
        
        //Update T
        /*vector<pair<id_type,size_t>>::iterator iter5; //Doubt
        for(iter5=temp_T.begin();iter5<temp_T.end();iter5++) //Doubt
            if(iter5->first==v_aggl)
                T=iter5->second;*/
        T=(temp_T.find(v_aggl))->second;
        
        //Update output
        output.push_back(make_pair(v_aggl,R));
    }  
}
//-----------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------
void CDLib::local_community_clauset_modified(const graph& g, id_type src, size_t k, vector< pair<id_type,double> >& output)
{
    node_set C,U;
    C.insert(src);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first); 
    
    double Z=0;
    double I=0;
    double T=g.get_node_out_weight(src);
    
    //Update output
    output.push_back(make_pair(src,Z));
    
    while (C.size() < k)
    {
        double max=-numeric_limits<double>::infinity();
        long cntr=-1; //Doubt
        long no_max=1;
        vector<id_type> V;
        
        node_set X;
        double x,y;
        
        for(node_set::iterator iter=U.begin();iter!=U.end();iter++)
        {
            //compute delta_Zj
            double delta;

            id_type v=*iter; //vertex vj
            
            x = get_neighs_in_set(v,g,C,X); //X = set of all neighbors of vj in C  //x = No. of edges in T that terminate at vj
            
            y = g.get_node_out_weight(v)-x; //No. of edges that will be added to T by the agglomeration of vj
            
            delta=((I+x)/(T+y)) - Z;
            
            if(delta>max)
            {
                V.push_back(v);
                cntr+=no_max;
                max=delta;
                no_max=1;
                //store v and corresponding new value of T
                //temp_T.insert(make_pair(v,T-z+y));
            }
            else if(delta==max)
            {
                V.push_back(v);
                no_max+=1;
                //store v and corresponding new value of T
                //temp_T.insert(make_pair(v,T-z+y));
            } 
        }
        
        //Selecting the node to be agglomerated v_aggl (breaking ties randomly if needed)
        //Choose a random number rand in [cntr,V.size()-1] such that v_aggl=V[rand] 
        RandomGenerator<long> rnd_gen(cntr,V.size()-1); //Doubt
        id_type v_aggl=V[rnd_gen.next()]; //Doubt
        
        //Update I and T
        x = get_neighs_in_set(v_aggl,g,C,X); 
            
        y = g.get_node_out_weight(v_aggl)-x; 
        
        I=I+x;
        T=T+y;
                
        //Update U
        U.erase(v_aggl);
        get_neighs_not_in_set(v_aggl,g,C,X); //X = all neighbors of v_aggl not in C
        for(node_set::iterator iter4=X.begin();iter4!=X.end();iter4++)
            U.insert(*iter4);
        
        //Update C
        C.insert(v_aggl);
        
        //Update Z
        Z=Z+max;
        
        //Update output
        output.push_back(make_pair(v_aggl,Z));        
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------



void bubble_sort_acc_degrees(vector<id_type>& N, const graph& g)
{
    id_type temp;
    for(size_t i=0; i<(N.size()-1); i++)
        for(size_t j=N.size()-1; j>i; j--)
            if(g.get_node_out_degree(N[j]) < g.get_node_out_degree(N[j-1]))
            {
                //swap
                temp=N[j];
                N[j]=N[j-1];
                N[j-1]=temp;
            }
}

bool found_in_vector(id_type v, vector<id_type>& S, size_t& index)
{
    for(size_t i=0; i<S.size(); i++) //Linear Search
        if(v==S[i])
        {
            index=i;
            return 1;
        }
    return 0;
}

double get_no_neighs_in_vector(id_type v, const graph& g, vector<id_type>& subgraph)
{
    //create an unordered set copy of the subgraph vector for time complexity issues
    node_set S;
    for(size_t i=0; i<subgraph.size(); i++)
        S.insert(subgraph[i]);
    
    double retval = 0;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
    {
        if(S.find(aeit->first) != S.end())
            retval += 1;
    }
    return retval;
}

void delete_from_vector_at_position(id_type v, vector<id_type>& N, size_t i)
{
    for(size_t j=i; j<N.size()-1; j++)
        N[j]=N[j+1];
    N.pop_back();
}

void bfs_visitor_in_subgraph(const graph& g, vector<id_type>& subgraph, id_type source, node_set& visited)
{
    //create an unordered set copy of the subgraph vector for time complexity issues
    node_set S;
    for(size_t i=0; i<subgraph.size(); i++)
        S.insert(subgraph[i]);

    queue<id_type> q_bfs;
    q_bfs.push(source); //push the source node to be visited first
    while(!q_bfs.empty()) //you have nodes to visit
    {
        id_type current  = q_bfs.front();
        q_bfs.pop();
        if(visited.find(current) == visited.end()) //current node has not yet been visited
        {    
           //visit node
           visited.insert(current);
           //expand node
           for(adjacent_edges_iterator aeit = g.out_edges_begin(current);aeit != g.out_edges_end(current);aeit++)
           {
               if(visited.find(aeit->first) == visited.end() && S.find(aeit->first) != S.end())
               {
                   //push neighbors to be visited
                   q_bfs.push(aeit->first); //insert only those neighbors which have not been "visited" yet
               }
           }
        }
    }
}

size_t get_component_around_node_in_subgraph(const graph& g, vector<id_type>& S, id_type source, node_set& visited)
{
    bfs_visitor_in_subgraph(g,S,source,visited);
    return visited.size();
}

bool connected_on_removal_at_position(const graph& g, vector<id_type> S, id_type u, size_t i)
{
    //Note: S has been passed by value to avoid reflection of changes to the original subgraph S after deletion in this function
    //removing node u at position i in subgraph vector S in graph G
    delete_from_vector_at_position(u,S,i);
    node_set visited;
    return (get_component_around_node_in_subgraph(g,S,S[0],visited) == S.size());
}



//-----------------------------------------------------------------------------------------------------------------------------------------
void CDLib::LWP_2006(const graph& g, id_type src, vector<id_type>& output)
{
    vector<id_type> S, N, Q, deleteQ;
    size_t index;
    //initialization
    S.push_back(src);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        N.push_back(aeit->first);
    
    double I=0, E=g.get_node_out_degree(src), M=0;
    size_t iter=0;
    
    do
    {
        //vector<id_type> Q;
        iter=iter+1;
        //cout<<"\nITERATION "<<iter<<":"<<endl;
        
        Q.clear();
        
        //ADDITION STEP
        bubble_sort_acc_degrees(N,g); //sort N in order of increasing degrees
        
        for(long i=0; i<(long)N.size(); i++) //Doubt (i's datatype to allow for value -1)
        {
            id_type v=N[i]; //vertex vj
            //compute delta_M
            double x = get_no_neighs_in_vector(v,g,S);
            double y = g.get_node_out_degree(v)-x;
            double delta = (I+x)/(E-x+y) - M;
            
            if(delta>0)
            {
                //cout<<"\ndelta is (insertion):"<<delta<<endl;
                //cout<<"Node: "<<g.get_node_label(v)<<endl;
                S.push_back(v); //push into S
                delete_from_vector_at_position(v,N,i); //pop from N
                Q.push_back(v);
                I=I+x;
                E=E-x+y;
                M+=delta;
                i=i-1; //shift needed to prevent skipping next element in N
            }
        }
        
        //S, N and Q have been updated at this point
        
        //DELETION STEP
        size_t counter=0;
        do
        {
            counter+=1;
            deleteQ.clear(); 
            for(long i=0; i<(long)S.size(); i++) //Doubt (i's datatype to allow for value -1)
            {
                //if(S[i]==src) //Modification to the original algorithm: Deletion of the source vertex is never allowed
                //  continue; 
            
                id_type u=S[i]; //vertex ui
                //compute delta_M
                double x = get_no_neighs_in_vector(u,g,S);
                double y = g.get_node_out_degree(u)-x;
                double delta = (I-x)/(E-y+x) - M;
            
                if(delta>0 && connected_on_removal_at_position(g,S,u,i))
                {
                    //cout<<"\ndelta is (deletion):"<<delta<<endl;
                    //cout<<"Node: "<<g.get_node_label(u)<<endl;
                    delete_from_vector_at_position(u,S,i);
                    
                    //N.push_back(u); //Not a part of this algorithm as per Bagrow's paper
                    
                    deleteQ.push_back(u);
                    I=I-x;
                    E=E-y;
                    M+=delta;
                    if(found_in_vector(u,Q,index))
                        delete_from_vector_at_position(u,Q,index);
                    i=i-1; //shift needed to prevent skipping next element in S
                } 
            }
        } while(!deleteQ.empty());
        
        //cout<<"\nNo. of deletion loops in iteration "<<iter<<" are "<<counter<<endl;
        
        //ADD VERTICES TO N
        for(size_t i=0; i<Q.size(); i++)
            for(adjacent_edges_iterator aeit = g.out_edges_begin(Q[i]); aeit != g.out_edges_end(Q[i]); aeit++)
                if(found_in_vector(aeit->first, S, index)==0 && found_in_vector(aeit->first, N, index)==0)
                    N.push_back(aeit->first);
             
    } while (!Q.empty());
        
    //cout<<"\nNo. of iterations: "<<iter<<endl;
    output = S;
    if(M>1 && found_in_vector(src,S,index))
    {
        //cout<<endl<<I<<" "<<E<<" "<<M<<endl;
        //output = S;
    }
    else
    {
        //cout<<"\nNO MODULE FOUND BY LWP FOR SOURCE VERTEX "<<g.get_node_label(src)<<endl;
        //output = S;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------



double get_no_neighs_in_set(id_type v, const graph& g, node_set& C)
{
    double retval = 0;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(v); aeit != g.out_edges_end(v); aeit++)
        if(C.find(aeit->first) != C.end())
            retval += aeit->second;
    return retval;
}

bool not_found_in_deque(id_type w, deque<id_type>& Q)
{
    for(size_t i=0; i<Q.size(); i++)
        if(Q[i] == w)
            return 0;
    return 1;
}

void delete_from_deque(id_type u, deque<id_type>& Q)
{
    for(size_t i=0; i<Q.size(); i++)
    {
        if(Q[i] == u)
        {
            //delete from deque
            for(size_t j=i; j<Q.size()-1; j++)
                Q[j]=Q[j+1];
            Q.pop_back();
            i=i-1; //shift needed to prevent skipping next element in Q
        }
    }
}



//-----------------------------------------------------------------------------------------------------------------------------------------
void CDLib::VD_2011(const graph& g, id_type src, node_set& output,id_type k)
{
    deque<id_type> Q;
    node_set C;
    node_set visited;
    id_type u = src;
    C.insert(u);
    double I=0, E=g.get_node_out_degree(u), N=1, f, x, y, I2, E2, N2, f2;
    size_t iter_count = 0;// N_at_start;
    
    //Initialize C with the source vertex and its neighbors AND update I, E and N
    for(adjacent_edges_iterator aeit = g.out_edges_begin(u); aeit != g.out_edges_end(u); aeit++)
    {
        x = get_no_neighs_in_set(aeit->first,g,C);
        y = g.get_node_out_degree(aeit->first) - x;
        C.insert(aeit->first);
        I=I+x; //No. of internal edges of C
        E=E-x+y; //No. of external edges of C
        N=N+1; //No. of nodes in C
    }
    
    f = (2*I*I)/(N*(N-1)*(I+E)); //the objective function
    
    visited.insert(u); //the source vertex has now been visited
    
    //push all neighbors of src in Q for visiting them
    for(adjacent_edges_iterator aeit = g.out_edges_begin(u); aeit != g.out_edges_end(u); aeit++)
    {
            Q.push_back(aeit->first);
    }
    
    //cout<<"\n\nOutput C before iterations begin:\n\n";
    
    //for(node_set::iterator iter=C.begin(); iter!=C.end(); iter++)
        //cout<<g.get_node_label(*iter)<<endl<<endl;
    
    do
    {
        iter_count+=1;
        
        //cout<<"\n\nITERATION NO. "<<iter_count<<endl<<endl;
        
        //N_at_start = C.size(); //size of C at the start of the iteration
        
        //ADDITION PHASE
        if(iter_count!=1 && !Q.empty())
        {
            u = Q.front(); 
            Q.pop_front();
        //}
        
        visited.insert(u); //visit the node
        
        //push the required neighbors of u in Q and try agglomerating them into C
        for(adjacent_edges_iterator aeit = g.out_edges_begin(u); aeit != g.out_edges_end(u); aeit++)
        {
            id_type w = aeit->first;
            if(visited.find(w) == visited.end() && not_found_in_deque(w,Q))  //Doubt
                Q.push_back(w);
            if(C.find(w) == C.end())
            {
              x = get_no_neighs_in_set(w,g,C);
              y = g.get_node_out_degree(w) - x;
              //C.insert(w);
            
              I2=I+x;
              E2=E-x+y;
              N2=N+1;
              f2 = (2*I2*I2)/(N2*(N2-1)*(I2+E2));
            
              if(f2<f)
              {
                  //C.erase(w);
              }
              else
              {
                  C.insert(w);
                  I=I2; E=E2; N=N2;
                  f=f2;
                  //cout<<"\nNode inserted: "<<g.get_node_label(w)<<endl;
              }
            }
        }
        }
        
        //DELETION PHASE
        for(node_set::iterator iter=C.begin(); iter!=C.end();)
        {
            u = *iter;
            x = get_no_neighs_in_set(u,g,C);
            y = g.get_node_out_degree(u) - x;
            
            //cout<<g.get_node_label(u)<<"  ";
            //cout<<"\n\nInitial erase from C successful for node \n\n"<<g.get_node_label(u)<<endl<<endl;
            
            I2=I-x;
            E2=E+x-y;
            N2=N-1;
            f2 = (2*I2*I2)/(N2*(N2-1)*(I2+E2));
            
            if(f2<f)
            {
                iter++;
            }
            else
            {
                iter++;
                C.erase(u);
                I=I2; E=E2; N=N2;
                f=f2;
                //cout<<"\nNode deleted: "<<g.get_node_label(u)<<endl;
                delete_from_deque(u,Q);
            }
        }
        
        //cout<<"\nN is "<<N<<endl;
        //cout<<"\nSize of C is "<<C.size()<<endl;
        
        //cout<<"\nQueue content at the end of iteration "<<iter_count<<endl;
        //for(size_t i=0; i<Q.size(); i++)
            //cout<<g.get_node_label(Q[i])<<"  ";
    
    } while(C.size() < k);//while (N_at_start != C.size()); //go to next iteration if community size changes, else exit
    
    //cout<<"\n\nNo. of iterations: "<<iter_count<<endl<<endl;
    
    //return C;
    output = C;    
}
//-----------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------
pair<double,node_set> CDLib::Bagrow_2007(const graph& g, id_type src)
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
//-----------------------------------------------------------------------------------------------------------------------------------------

void bfs_visitor_in_subgraph_set(const graph& g, node_set& C, id_type source, node_set& visited)
{
    //create an unordered set copy of the subgraph vector for time complexity issues
    //node_set S;
    //for(size_t i=0; i<subgraph.size(); i++)
    //    S.insert(subgraph[i]);

    queue<id_type> q_bfs;
    q_bfs.push(source); //push the source node to be visited first
    while(!q_bfs.empty()) //you have nodes to visit
    {
        id_type current  = q_bfs.front();
        q_bfs.pop();
        if(visited.find(current) == visited.end()) //current node has not yet been visited
        {    
           //visit node
           visited.insert(current);
           //expand node
           for(adjacent_edges_iterator aeit = g.out_edges_begin(current);aeit != g.out_edges_end(current);aeit++)
           {
               if(visited.find(aeit->first) == visited.end() && C.find(aeit->first) != C.end())
               {
                   //push neighbors to be visited
                   q_bfs.push(aeit->first); //insert only those neighbors which have not been "visited" yet
               }
           }
        }
    }
}

size_t get_component_around_node(const graph& g, node_set& C, id_type source, node_set& visited)
{
    bfs_visitor_in_subgraph_set(g,C,source,visited);
    return visited.size();
}

bool connected_on_removal(const graph& g, node_set C, id_type u)
{
    //Note: C has been passed by value to avoid reflection of changes to the original subgraph C after deletion in this function
    //deleting
    C.erase(u);
    node_set visited;
    id_type start_node = *(C.begin());
    return (get_component_around_node(g,C,start_node,visited) == C.size());
}
//-----------------------------------------------------------------------------------------------------------------------------------------
void CDLib::My_Algorithm(const graph& g, id_type src, node_set& output)
{
    node_set C,U;
    vector<id_type> U_vector;
    C.insert(src);
    //C_vector.push_back(src);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(src); aeit != g.out_edges_end(src); aeit++)
        U.insert(aeit->first);
    //for(node_set::iterator iter=U.begin(); iter!=U.end();iter++)
    //    U_vector.push_back(*iter);
        
    double I=0, E=g.get_node_out_degree(src), n=1, nb=1, nu=g.get_node_out_degree(src);
    double I_new, E_new, n_new, nb_new, nu_new;
    double rho = 0, delta_rho, x, y;
    
    vector<id_type> Q, deleteQ;
    bool flag;
    id_type vj, ui, v;
    size_t index;
    
    do
    {
        //cout<<"\n\nNEW ITERATION STARTS:\n\n";
        //ADDITION PHASE
        //cout<<"\nentered addition phase\n";
        //cout.flush();
        Q.clear();
        
        U_vector.clear();
        for(node_set::iterator iter=U.begin(); iter!=U.end();iter++)
            U_vector.push_back(*iter);
        
        do
        {
            flag = 0;
            
            for(long i=0; i<(long)U.size(); i++)
            //for(node_set::iterator iter=U.begin(); iter!=U.end();) //change way of iterating over U!!!! 
            {
                vj = U_vector[i];
                if(U.find(vj)!=U.end())
                {
                //vj = *iter;
                //compute delta_rho
                x = get_no_neighs_in_set(vj,g,C);
                y = g.get_node_out_degree(vj) - x;
                
                I_new = I+x; E_new = E-x+y; n_new = n+1;
                
                
                nu_new = nu-1;
                for(adjacent_edges_iterator aeit = g.out_edges_begin(vj); aeit != g.out_edges_end(vj); aeit++)
                {
                    v = aeit->first;
                    if(C.find(v)==C.end() && U.find(v)==U.end())
                        nu_new = nu_new+1;
                }
                
                if(y>0)
                    nb_new = nb+1;
                for(adjacent_edges_iterator aeit = g.out_edges_begin(vj); aeit != g.out_edges_end(vj); aeit++)
                {
                    v = aeit->first;
                    if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==1)
                        nb_new = nb_new-1;
                }
                
                //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*(n_new-1)*E_new) - rho;
                //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*E_new) - rho;
                //delta_rho = ((pow(I_new,0.5))/(n_new*(n_new-1)/2))/(E_new/(nu_new*nb_new)) -rho;
                //delta_rho = (I_new)/(E_new) - rho;
                //delta_rho = E_new/(nu_new*nb_new);
                delta_rho = (2*I_new*nb_new)/(n_new*E_new) - rho;
                
                if(delta_rho > 0)
                {
                    //cout<<"\nnode "<<g.get_node_label(vj)<<" agglomerated"<<endl;
                    C.insert(vj);
                    //iter++;
                    U.erase(vj);
                    for(adjacent_edges_iterator aeit = g.out_edges_begin(vj); aeit != g.out_edges_end(vj); aeit++) //REQUIRED!!!!
                    {
                        v = aeit->first;
                        if(C.find(v)==C.end() && U.find(v)==U.end())
                            U.insert(v);
                    }
                    Q.push_back(vj);
                    flag = 1;
                    rho+=delta_rho;
                    I = I_new; E = E_new; n = n_new;
                    nb = nb_new; nu = nu_new;
                }
                }
                //else
                //    iter++;
            }
            
        } while(flag);
        
        //addition over
        
        //C_vector.clear();
        //for(node_set::iterator iter=C.begin(); iter!=C.end();iter++) //can be done incrementally inside the addition phase!
        //    C_vector.push_back(*iter);
        
        //DELETION PHASE
        //cout<<"entered deletion phase";
        //cout.flush();
        do
        {
            deleteQ.clear();
            
            for(node_set::iterator iter=C.begin(); iter!=C.end();)
            {
                ui = *iter;
                //if(ui!=src)
                //{
                //compute delta_rho
                x = get_no_neighs_in_set(ui,g,C);
                y = g.get_node_out_degree(ui) - x;
                
                I_new = I-x; E_new = E+x-y; n_new = n-1;
                
                nu_new = nu+1;
                
                for(adjacent_edges_iterator aeit = g.out_edges_begin(ui); aeit != g.out_edges_end(ui); aeit++)
                {
                    v = aeit->first;
                    if(C.find(v)==C.end() && get_no_neighs_in_set(v,g,C)==1)
                        nu_new = nu_new-1;
                }
                
                if(y>0)
                    nb_new = nb-1;
                
                for(adjacent_edges_iterator aeit = g.out_edges_begin(ui); aeit != g.out_edges_end(ui); aeit++)
                {
                    v = aeit->first;
                    if(C.find(v)!=C.end() && get_no_neighs_in_set(v,g,U)==0)
                        nb_new = nb_new+1;
                }
                
                //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*(n_new-1)*E_new) - rho;
                //delta_rho = (2*I_new*nb_new*nu_new)/(n_new*E_new) - rho;
                //delta_rho = ((pow(I_new,0.5))/(n_new*(n_new-1)/2))/(E_new/(nu_new*nb_new)) -rho;
                //delta_rho = (I_new)/(E_new) - rho;
                delta_rho = (2*I_new*nb_new)/(n_new*E_new) - rho;
                
                if(delta_rho>0 && connected_on_removal(g,C,ui))
                {
                    //cout<<"\nnode "<<g.get_node_label(ui)<<" deleted"<<endl;
                    iter++;
                    C.erase(ui);
                    U.insert(ui);
                    for(adjacent_edges_iterator aeit = g.out_edges_begin(ui); aeit != g.out_edges_end(ui); aeit++)
                    {
                        v = aeit->first;
                        if(U.find(v)!=U.end() && get_no_neighs_in_set(v,g,C)==0)
                            U.erase(v);
                    }               
                
                    deleteQ.push_back(ui);
                
                    if(found_in_vector(ui,Q,index))
                        delete_from_vector_at_position(ui,Q,index);
                
                    rho+=delta_rho;
                    I = I_new; E = E_new; n = n_new;
                    nb = nb_new; nu = nu_new;
                }
                else
                    iter++;
                //}
                //else
                //    iter++;
            }
        } while(!deleteQ.empty());
        
        //ADD REQUIRED VERTICES TO U
        //for(size_t k=0; k<Q.size(); k++)
        //    for(adjacent_edges_iterator aeit = g.out_edges_begin(Q[k]); aeit != g.out_edges_end(Q[k]); aeit++)
        //        if(C.find(aeit->first)==C.end() && U.find(aeit->first)==U.end())
        //            U.insert(aeit->first);
    
    } while(!Q.empty());
    
    if(rho>1 && C.find(src)!=C.end())
    {
        //cout<<endl<<I<<" "<<E<<" "<<M<<endl;
        output = C;
    }
    else
    {
        //cout<<"\nNO MODULE FOUND BY MY_ALGO FOR SOURCE VERTEX "<<g.get_node_label(src)<<endl;
        //output = C;
    }
}

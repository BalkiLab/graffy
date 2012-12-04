/* 
 * File:   random_graph.cpp
 * Author: bharath
 * 
 * Created on 16 April, 2012, 4:25 PM
 */

#include "random_graph.h"


using namespace CDLib;

void CDLib::init_empty_graph(graph& g,size_t size)
{
    g.clear();
    if(g.is_directed()) g.convert_to_undirected();
    if(g.is_weighted()) g.convert_to_unweighted(0);
    for(id_type i=0; i<size;i++)  g.add_node();
    
}

void CDLib::generate_erdos_renyi_graph(graph& g, id_type num_nodes,double p)
{
    if(p >= 0 && p<=1)
    {
        init_empty_graph(g,num_nodes);
        RandomGenerator<double> p_gen(0,1,1);
        for(id_type i=0; i<num_nodes;i++)
            for(id_type j=0; j<num_nodes;j++)
                if(p_gen.next() <= p) g.add_edge(i,j,1); 
    }
}

void CDLib::generate_erdos_renyi_graph(graph& g, id_type num_nodes,id_type num_edges)
{
    double p = (double)num_edges/(num_nodes * (num_nodes - 1));
    generate_erdos_renyi_graph(g,num_nodes,p);
}

void CDLib::generate_scale_free_graph(graph& g, id_type num_nodes,id_type num_edges,double alpha, double beta)
{
    init_empty_graph(g,num_nodes);
    RandomGenerator<id_type> x_gen(0,num_nodes),from_gen(0,num_nodes),to_gen(0,num_nodes);
    vector<id_type> outdegrees(num_nodes,0.0);
    for(id_type i=0; i<num_nodes;i++)
        outdegrees[i] = (id_type)x_gen.exp_next(alpha,beta);
    for(id_type i=0 ;i<num_edges;i++)
    {
        id_type from_id = from_gen.next(),to_id = to_gen.next();
        if(outdegrees[from_id] && outdegrees[to_id] && from_id != to_id && !g.get_edge_weight(from_id,to_id)) 
            g.add_edge(from_id,to_id,1);
    }
}


void CDLib::generate_planted_partition_graph(graph& g, id_type num_comms, id_type comm_size, double pin, double pout,vector< node_set>& communities)
{
    if(pin >=0 && pout >=0 && pin<= 1 && pout <=1)
    {
        id_type num_nodes = comm_size*num_comms;
        init_empty_graph(g,num_nodes);
        communities.assign(num_comms,node_set());
        Uniform01RandomGeneratorMT p_gen;
        for(id_type i=0; i<num_nodes; i++)
            if(i % comm_size) communities[i/comm_size].insert(i);
        for(id_type i=0;i<g.get_num_nodes();i++)
        {
            id_type comm_id_i = i/comm_size;
            for(id_type j=0;j<g.get_num_nodes();j++)
            {
                id_type comm_id_j = j/comm_size;
                double p = p_gen.next();
                if(comm_id_i == comm_id_j && p > pin) g.add_edge(i,j,1);
                else if(p > pout) g.add_edge(i,j,1);
            }
        }
    }
}

void CDLib::generate_ring_graph(graph& g, id_type size)
{
    init_empty_graph(g,size);
    for(id_type i=0;i<g.get_num_nodes();i++)g.add_edge(i,(i+1)%size,1);
}

void CDLib::generate_star_graph(graph& g, id_type size)
{
    init_empty_graph(g,size);
    for(id_type i=1;i<g.get_num_nodes();i++)g.add_edge(0,i,1);
}

void CDLib::generate_clique_graph(graph& g, id_type size)
{
    init_empty_graph(g,size);
    for(id_type i=0;i<g.get_num_nodes();i++)
        for(id_type j=0;j<i;j++)
            if(i!=j) g.add_edge(i,j,1);
}

void CDLib::generate_spoke_graph(graph& g,id_type size)
{
    init_empty_graph(g,size);
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        if(i)g.add_edge(i,(i+1)%size,1);
        g.add_edge(i,(i+1)%size,1);
    }
}

void CDLib::generate_de_bruijn_graph(graph& g,id_type num_symbols,id_type sequence_length)
{
    if(num_symbols && sequence_length)
    {
        id_type size = (unsigned long)pow((double)num_symbols,(double)sequence_length);
        init_empty_graph(g,size);
        for (id_type i=0; i<g.get_num_nodes(); i++)
        {
            id_type basis=(i*num_symbols) % g.get_num_nodes();
            for (id_type j=0; j<num_symbols; j++) g.add_edge(i,basis+j,1);    
        }
    }
}

void CDLib::generate_chord_graph(graph& g,id_type num_nodes)
{
    init_empty_graph(g,num_nodes);
    for (id_type i=0; i<g.get_num_nodes(); i++)
    {
        g.add_edge(i,(i-1) %  g.get_num_nodes(),1);
        for(id_type j=1; j<= g.get_num_nodes()/2; j*=2 )
            g.add_edge(i,(i+j) % g.get_num_nodes(),1);
    }
}

void CDLib::generate_LEET_chord_graph(graph& g, id_type num_nodes)
{
    init_empty_graph(g,num_nodes);
    id_type cluster_size = log2(g.get_num_nodes());
    for (id_type i=0; i<g.get_num_nodes(); i++) 
    {       
        g.add_edge(i,(i-1) %  g.get_num_nodes(),1);
        g.add_edge(i,(i+1) %  g.get_num_nodes(),1);
        id_type id_in_cluster = g.get_num_nodes() % static_cast<id_type>(cluster_size);
        g.add_edge(i,id_in_cluster*(1+static_cast<id_type>(cluster_size)),1);
    }
}

bool cmp(const vector<id_type>& lhs, const vector<id_type>& rhs) { return lhs.size() < rhs.size();}

void CDLib::generate_kademlia_graph(graph& g,id_type num_nodes)
{
    init_empty_graph(g,num_nodes);
    unsigned int num_bits = log2(num_nodes);
    for(unsigned int i=0;i<num_nodes;i++)
        for(unsigned int j=0; j<num_bits;j++)
                g.add_edge(i,(i ^ ((1 << j)))%num_nodes,1) ;
}

//void CDLib::generate_kademlia_graph(graph& g,id_type num_nodes, id_type bucket_length)
//{
//    init_empty_graph(g,num_nodes);
//    id_type bit_lmt = 0, edg;
//    for (id_type i=0; i<g.get_num_nodes();i++)
//    {
//        for(id_type j= 1, l = 0; j<= g.get_num_nodes()/2; j*=2,l++ )
//        {
//            bit_lmt = i >> l;
//            bit_lmt = (bit_lmt % 2 == 0) ? bit_lmt + 1 : bit_lmt -1;
//            bit_lmt <<=  l;
//            for(id_type k = 0; k < j && k < bucket_length && k < (num_nodes - (1 << l)); k++)
//            {
//                edg = bit_lmt + (rand() % (1 << l));
//                while (g.get_edge_weight(i,edg))
//                    edg = bit_lmt + (rand() % (1 << l));
//                g.add_edge(i, edg,1);
//            }
//        }
//    }
//}
/*void CDLib::generate_kademlia_graph(graph& g,id_type num_bits, id_type bucket_length)
{
    id_type num_nodes = static_cast<id_type>(pow(2,num_bits));
    init_empty_graph(g,num_nodes);
    vector< vector<node_set> > k_buckets(num_nodes,vector<node_set>(num_bits,node_set()));
    for(id_type i=0;i<num_nodes;i++)
    {
        vector< vector<id_type> > tree_part(num_bits,vector<id_type>());
        for(id_type j=0;j<num_nodes;j++)
        {
            for(id_type k=1;k<=num_bits;k++)
            {
                id_type iand = i & ((num_nodes - 1) << k);
                id_type jand = j & ((num_nodes - 1) << k);
                if(iand == jand)
                {
                    tree_part[k-1].push_back(j);
                    if(k<= bucket_length)
                    {
                        k_buckets[i][k-1].insert(j);
                        k_buckets[j][k-1].insert(i);
                    }
                }
            }
        }
        for(id_type k=bucket_length+1;k<num_bits;k++)
        {
            random_shuffle(tree_part[k-1].begin(),tree_part[k-1].end());
            id_type index = 0;
            while(k_buckets[i][k-1].size()<bucket_length && index <tree_part[k-1].size() )
            {
                id_type j = tree_part[k-1][index];
                if(k_buckets[j][k-1].size() < bucket_length)
                {
                    k_buckets[i][k-1].insert(j);
                    k_buckets[j][k-1].insert(i);
                }
            }
        }
    }
    for(id_type i=0;i<num_nodes;i++)
        for(id_type k=0;k<num_bits;k++)
            for(node_set::iterator nit = k_buckets[i][k].begin(); nit != k_buckets[i][k].end();nit++)
                g.add_edge(i,*nit,1);
}*/



//
//long CDLib::generate_lfr_graph(graph& g,id_type num_nodes,id_type num_edges, id_type max_degree,double tau, double tau2,double mixing_parameter,vector<node_set>& comms)
//{
//    init_empty_graph(g,num_nodes);
//    double average_k = static_cast<double>(num_edges)/static_cast<double>(num_nodes);
//    long overlapping_nodes = 0;
//    long overlap_membership = 0;
//    long nmin = unlikely;
//    long nmax = unlikely;
//    bool fixed_range = false;
//    bool excess = 0;
//    bool defect = 0;
//    deque<set<long> > E; 
//    deque<deque<long> > member_list; 
//    deque<deque<long> > link_list;
//    init_empty_graph(g,num_nodes);
//    long retval = lfr_binary_networks::benchmark(excess,defect,num_nodes,average_k,max_degree,tau,tau2,mixing_parameter,overlapping_nodes,overlap_membership,nmin,nmax,fixed_range,E,member_list,link_list);
//    for(id_type u = 0; u < E.size(); u++) 
//    {
//        set<long>::iterator itb = E[u].begin();
//        while (itb != E[u].end()) g.add_edge(u,*itb++,1);
//    }
//    vector<id_type> labels;
//    for (id_type i = 0; i < member_list.size(); i++) 
//        for (id_type j = 0; j < member_list[i].size(); j++)
//            labels[i] = member_list[i][j];
//    convert_labels_to_communities(labels,comms);
//    return retval;
//}


void CDLib::generate_configuration_model(graph& g, vector<id_type>& degree_sequence)
{
    if(degree_sequence.size()>0)
    {
        init_empty_graph(g,degree_sequence.size());
        vector<id_type> nodes_with_non_0_degree;
               
        for(id_type i=0;i<degree_sequence.size();i++)
        {
            if(degree_sequence[i]!=0)
            {
                nodes_with_non_0_degree.push_back(i);
            }
        }
        
        UniformRandomGeneratorAkash<id_type> rand;
        for(id_type i=0;i<degree_sequence.size();i++)
        {
            id_type remaining_degree=degree_sequence[i];
            for(id_type j=0;j<remaining_degree;j++)
            {
                back:
                id_type R=rand.next(nodes_with_non_0_degree.size());
                if(degree_sequence[i]==1&&i==nodes_with_non_0_degree[R])
                    goto back;
                
                g.add_edge(i,nodes_with_non_0_degree[R],1);
                
                degree_sequence[i]--;
                degree_sequence[nodes_with_non_0_degree[R]]--;
                
                if(degree_sequence[nodes_with_non_0_degree[R]]==0)
                {
                    if(nodes_with_non_0_degree[R]==i)
                        break;
                    nodes_with_non_0_degree.erase(nodes_with_non_0_degree.begin()+R);                       
                }                
                if(degree_sequence[i]==0)
                    break;
            }
            if(remaining_degree!=0)
            {
                nodes_with_non_0_degree.erase(nodes_with_non_0_degree.begin());
            }
        } 
    }
}

/*this is preferential attachment model. it follows power law degree distribution.
  this is directed graph model. */

void CDLib::generate_prices_model(graph& g,size_t num_nodes, size_t num_of_out_degree, size_t in_degree_constant)
{
    init_empty_graph(g,num_nodes);
    vector<id_type> vertices_pointed_by_edges;
    UniformRandomGeneratorAkash<id_type> randint;
    UniformRandomGeneratorAkash<double> randdouble;
    
    double probability=num_of_out_degree/(num_of_out_degree + in_degree_constant);
    
    for(id_type i=0;i<num_nodes;i++)
    {
        while(g.get_node_out_degree(i)<num_of_out_degree)
        {
            double R1=randdouble.next(1);
            if(R1<probability)
            {
                id_type R2=randint.next(vertices_pointed_by_edges.size());
                g.add_edge(i,vertices_pointed_by_edges[R2],1);
                vertices_pointed_by_edges.push_back(vertices_pointed_by_edges[R2]);                      
            }
            else
            {
                id_type R3=randint.next(num_nodes);
                g.add_edge(i,R3,1);
                vertices_pointed_by_edges.push_back(R3);
            }
        }
    }
}

/*this is preferential attachment model. it follows power law degree distribution.
  it is undirected graph model. */

void CDLib::generate_barabasi_albert_model(graph& g, size_t num_nodes, size_t min_degree_of_node)
{
    init_empty_graph(g,num_nodes);
    vector<id_type> vertices_with_edges;
    UniformRandomGeneratorAkash<id_type> randint;
    UniformRandomGeneratorAkash<double> randdouble;
    
    for(id_type k=1;k<=min_degree_of_node;k++)
    {
        g.add_edge(0,k,1);
        vertices_with_edges.push_back(k);
    }
        
    for(id_type i=1;i<num_nodes;i++)
    {
        while(g.get_node_in_degree(i)<min_degree_of_node)
        {
            double R1=randdouble.next(1);
            if(R1<0.5)
            {
                id_type R2=randint.next(vertices_with_edges.size());
                g.add_edge(i,vertices_with_edges[R2],1);
                vertices_with_edges.push_back(vertices_with_edges[R2]);  
            }
            else
            {
                id_type R3=randint.next(i);
                g.add_edge(i,R3,1);
                vertices_with_edges.push_back(R3);
            }
        }
    }
}

/*degree distribution of this model follows power law distribution. 
 * this is a most suitable synthetic model to real-world network like peer to peer networks and citation networks.*/
void CDLib::generate_vertex_copying_model(graph& g,size_t num_nodes,size_t num_of_out_degree,size_t num_of_vertices_at_initial,double probability_to_copy_from_existing_vertex)
{
    if (num_of_vertices_at_initial > num_of_out_degree && probability_to_copy_from_existing_vertex>=0 && probability_to_copy_from_existing_vertex<=1) 
    {
        UniformRandomGeneratorAkash<wt_t> randdouble;
        UniformRandomGeneratorAkash<id_type> randint;

        init_empty_graph(g, num_nodes);

        for (id_type i = 0; i < num_of_vertices_at_initial; i++) 
        {
            while(g.get_node_out_degree(i)<num_of_out_degree)
            {
                back:
                id_type R1=randint.next(num_of_vertices_at_initial);
                if(R1==i)
                    goto back;
                else
                {
                    g.add_edge(i,R1,1);
                }
            }
        }
        
        for(id_type i=num_of_vertices_at_initial;i<num_nodes;i++)
        {
            id_type R2=randint.next(i);
            vector<id_type> vertices_pointed_by_R2;
            for(adjacent_edges_iterator aeit=g.out_edges_begin(R2);aeit!=g.out_edges_end(R2);aeit++)
            {
                vertices_pointed_by_R2.push_back(aeit->first);
            }
            while(g.get_node_out_degree(i)<num_of_out_degree)
            {
                wt_t R3=randdouble.next(1);
                if(R3<probability_to_copy_from_existing_vertex)
                {
                    g.add_edge(i,vertices_pointed_by_R2[g.get_node_out_degree(i)],1);
                    
                }
                else
                {
                    A:
                    id_type R4=randint.next(num_nodes);
                    if(R4!=i)
                    {
                        g.add_edge(i,R4,1);
                    }    
                    else 
                        goto A;
                }
            }
        }
    }
    else
        cout<<"\nnum_of_vertices_at_initial should greater than num_of_out_degree\n";
}

/*this model is used to generate graph with high clustering coefficient.
 this model matches sports network like american football league.*/
void CDLib::generate_small_world_model(graph& g,size_t num_nodes,size_t degree_of_each_vertex,double probability_to_replace_edge)
{
    if(probability_to_replace_edge>=0 && probability_to_replace_edge<=1)
    {
        UniformRandomGeneratorAkash<id_type> randint;
        UniformRandomGeneratorAkash<double> randdouble;
        
        init_empty_graph(g,num_nodes);
        
        for(id_type i=0;i<num_nodes;i++)
        {
            for(id_type j=0;j<degree_of_each_vertex/2;j++)
            {
                g.add_edge(i,(i+j+1)%num_nodes,1);
            }
            if(degree_of_each_vertex%2==1)
            {
                g.add_edge(i,(i+(degree_of_each_vertex/2)+1)%num_nodes,1);
            }
        }
        for(id_type i=0;i<num_nodes;i++)
        {
            for(id_type j=0;j<degree_of_each_vertex;j++)
            {
                double R1=randdouble.next(1);
                if(R1<probability_to_replace_edge)
                {
                    back:
                    id_type R2=randint.next(num_nodes);
                    if(R2==i)
                        goto back;
                    else
                    {
                        g.add_edge(i,R2,1);
                    }
                }
            }
        }
        
    }
    else
        cout<<"\nprobability value might be wrong";
}


/*this model generates graph with 128 nodes distributed among 4 communities at initial time step.
 at each time step each community exchange three nodes randomly among them.*/
void CDLib::generate_evolutionary_model_128_nodes_4_communities(vector<graph>& g,size_t num_inter_community_edges,size_t num_timesteps)
{
    if(num_inter_community_edges<16)
    {
        vector<id_type> comm[4];
        graph g_t(0,0);
        init_empty_graph(g_t,128);
        
        double probability_inter_comm_edge = num_inter_community_edges/16.0;
        double probability_intra_comm_edge = 1-probability_inter_comm_edge;
        
        UniformRandomGeneratorAkash<double> randdouble;
        UniformRandomGeneratorAkash<id_type> randint;
        
        for(id_type t=0;t<num_timesteps;t++)
        {
            if(t==0)
            {
                for(id_type k=0;k<32;k++)
                {
                    comm[0].push_back(k);
                    comm[1].push_back(k+32);
                    comm[2].push_back(k+64);
                    comm[3].push_back(k+96);
                }
                for(id_type i=0;i<4;i++)
                {
                    for(id_type j=0;j<32;j++)
                    {
                        id_type degree_of_node = g_t.get_node_in_degree(comm[i][j]);
                        while (g_t.get_node_in_degree(comm[i][j]) < 16 - degree_of_node) 
                        {
                            back:
                            id_type R1 = randint.next(32);
                            if (R1 == comm[i][j])
                                goto back;
                            double R2 = randdouble.next(1);

                            if (R2 < probability_intra_comm_edge) 
                            {
                                g_t.add_edge(comm[i][j],comm[i][R1],1);
                            }
                            else
                            {
                                id_type R3=randint.next(4);
                                g_t.add_edge(comm[i][j],comm[R3][R1],1);
                            }
                        }
                    }
                }
            }
            else
            {
                for(id_type i=0;i<4;i++)
                {
                    id_type R5=randint.next(4);
                    for(id_type j=0;j<3;j++)
                    {
                        st:
                        id_type R4=randint.next(comm[i].size());
                        if(R4>31)
                            goto st;
                        if(R5==i)
                            R5=(R5+1)%4;
                        vector<id_type> node_list;
                        for(adjacent_edges_iterator aeit=g_t.in_edges_begin(comm[i][R4]);aeit!=g_t.in_edges_end(comm[i][R4]);aeit++)
                        {
                            node_list.push_back(aeit->first);
                        }
                        for(id_type k=0;k<node_list.size();k++)
                        {
                            g_t.remove_edge(comm[i][R4],node_list[k]);
                        }
                        comm[R5].push_back(comm[i][R4]);
                        comm[i].erase(comm[i].begin()+R4);
                        R5=(R5+1)%4;
                    }
                }
                for(id_type i=0;i<4;i++)
                {
                    for(id_type j=0;j<3;j++)
                    {
                        while (g_t.get_node_in_degree(comm[i][31-j]) < 16) 
                        {
                            back1:
                            id_type R1 = randint.next(32);
                            if (R1 == comm[i][31-j])
                                goto back1;
                            double R2 = randdouble.next(1);

                            if (R2 < probability_intra_comm_edge) 
                            {
                                g_t.add_edge(comm[i][31-j],comm[i][R1],1);
                            }
                            else
                            {
                                id_type R3=randint.next(4);
                                g_t.add_edge(comm[i][31-j],comm[R3][R1],1);
                            }
                        }
                    }
                }
            }
            g.push_back(g_t);
        }
    }
    else
        cout<<"\n num_inter_community_edges should less then 16.";
}

void CDLib::generate_evolutionary_model_800_nodes_4_communities(vector < vector < vector < double > > >& points,vector<double>& x_coordinates,vector<double>& y_coordinates,double variance,size_t num_timesteps)
{
    if(variance>0)
    {
        vector< vector <double> > points_t;
        
        GaussianRandomGenerator<double> gaussrand;
        for(id_type i=0;i<800;i++)
        {
            double x,y;
            vector<double> point;
            x = gaussrand.next(x_coordinates[(int)(i/200)],variance);
            y = gaussrand.next(y_coordinates[(int)(i/200)],variance);
            point.push_back(x);
            point.push_back(y);
            points_t.push_back(point);
        }
        points.push_back(points_t);
        int c;
        back:
        cout<<"\n enter your choice : \n 1:: adding just stationary noise. \n 2:: adding non stationary noise. \n";
        cin>>c;
        if(c==1)
        {
           for (id_type t = 1; t < num_timesteps; t++) 
           {
               vector< vector <double> > points_t;
               for(id_type i=0;i<800;i++)
               {
                   double Rx1,Ry1,xi,yi;
                   vector<double> point;
                   Rx1 = gaussrand.next(0,0.5);
                   Ry1 = gaussrand.next(0,0.5);
                   
                   xi = points[0][i][0] + Rx1;
                   yi = points[0][i][1] + Ry1;
                   
                   point.push_back(xi);
                   point.push_back(yi);
                   
                   points_t.push_back(point);
               }
               points.push_back(points_t);
           }
           
        }
        else if(c==2)
        {
            for (id_type t = 1; t < num_timesteps; t++) 
           {
               vector< vector <double> > points_t;
               for(id_type i=0;i<800;i++)
               {
                   double Rx1,Ry1,xi,yi,angle;
                   vector<double> point;
                   
                   angle = gaussrand.next(0,(3.14159265/4.0));
                   
                   xi = points[0][i][0] * cos(angle) - points[0][i][1] * sin(angle);
                   yi = points[0][i][0] * sin(angle) + points[0][i][1] * cos(angle);
                   
                   Rx1 = gaussrand.next(0,0.5);
                   Ry1 = gaussrand.next(0,0.5);
                   
                                      
                   point.push_back(xi+Rx1);
                   point.push_back(yi+Ry1);
                   
                   points_t.push_back(point);
               }
               points.push_back(points_t);
           }
        }
        else
        {
            cout<<" \n enter valid choice:\n";
            goto back;
        }
    }
    else
        cout<<"\n variance should be greater than 0";
}

/*degree_dist_controlling_parameter controls degree distribution between power law and exponential degree distribution as it varies between 0 to 1.
  it is undirected graph model.*/

double energy(graph& g,double degree_dist_controlling_parameter,id_type nC2,double Dlinear);



double energy(graph& g,double degree_dist_controlling_parameter,id_type nC2,double Dlinear)
{
    id_type num_edges = g.get_num_edges();
    id_type num_nodes = g.get_num_nodes();
    
    vector<double> distance;
    vector< vector <id_type> > seq;
    
    id_type total_of_min_distance=0;
    
    for(id_type i=0;i<num_nodes;i++)
    {
        single_source_shortest_paths_bfs(g,i,distance,seq);
        for(id_type j=0;j<num_nodes;j++)
        {
            total_of_min_distance = total_of_min_distance + distance[j];
        }
    }
    double normalized_num_of_link = (double)num_edges/nC2;
    double avg_min_distance = (double)total_of_min_distance/(nC2*2);
    
    double normalized_distance = avg_min_distance/Dlinear;
    double Energy = degree_dist_controlling_parameter*normalized_distance + (1-degree_dist_controlling_parameter)*normalized_num_of_link;
    
    return Energy;
}

void CDLib::generate_ferrer_i_cancho_model(graph& g,size_t num_nodes, size_t max_failure_allowed,double degree_dist_controlling_parameter,double probability_to_alter_edge,double initial_probability_of_edge)
{
    if(degree_dist_controlling_parameter<=1 && degree_dist_controlling_parameter>=0 && initial_probability_of_edge<=1 && initial_probability_of_edge>=0 && probability_to_alter_edge<=1 && probability_to_alter_edge>=0)
    {
        double Dlinear = ((num_nodes + 1) / 3.0);
        id_type nC2 = ((num_nodes * (num_nodes - 1)) / 2);
        
        UniformRandomGenerator<double> rand;

        graph g_t(0, 0);
        init_empty_graph(g,num_nodes);
        init_empty_graph(g_t,num_nodes);
      
        for(id_type i=0;i<num_nodes;i++)
        {
            for(id_type j=i+1;j<num_nodes;j++)
            {
                double R2=rand.next(1);
                if(R2<initial_probability_of_edge && i!=j)
                {
                   g.add_edge(i,j,1);
                } 
            }
        }

        id_type failure = 0;

        double Energy_old = energy(g, degree_dist_controlling_parameter, nC2, Dlinear);
        
        while (failure < max_failure_allowed) 
        {
            g_t=g;

            for (id_type i = 0; i < num_nodes; i++) 
            {
                for (id_type j = i+1 ; j < num_nodes; j++) 
                {
                    double R1=rand.next(1);
                    if(R1<probability_to_alter_edge)
                    {
                        if(g_t.get_edge_weight(i,j)==0)
                            g_t.add_edge(i,j,1);
                        else
                            g_t.remove_edge(i,j);
                    }
                }
            }

            double Energy_new = energy(g_t, degree_dist_controlling_parameter, nC2, Dlinear);

            if (Energy_new < Energy_old) 
            {
                g = g_t;
                Energy_old = Energy_new;
                failure = 0;
            } else
                failure++;
        }
    }
    else
        cout<<"\n last three parameter value should be between 0 and 1";
    
}

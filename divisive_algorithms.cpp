/* 
 * File:   divisive_algorithms.cpp
 * Author: bharath
 * 
 * Created on 5 April, 2012, 12:42 PM
 */

#include "divisive_algorithms.h"
using namespace CDLib;

CDLib::edge_radicchi::edge_radicchi(id_type from,id_type to,double ec): from_id(from),to_id(to),ecc(ec) {}
CDLib::edge_radicchi::edge_radicchi(): from_id(0),to_id(0),ecc(0) {}

void CDLib::girvan_newman_2002(const graph& g, dendrogram& dendro)
{
    if(is_connected_weakly(g))
    {
        graph gc(g);
        id_type component_counter = 1;
        while(gc.get_num_edges())
        {
                vector<double> bc;
                betweeness_centralities(gc,bc);
                double max_edge_bet = 0;
                id_type from_id=0,to_id=0;
                for(id_type i=0;i<gc.get_num_nodes(); i++)
                {
                        for(adjacent_edges_iterator aeit = gc.out_edges_begin(i); aeit != gc.out_edges_end(i); aeit++)
                        {
                            double edge_bet = bc[i] + bc[aeit->first];
                            if(edge_bet >= max_edge_bet)
                            {
                                max_edge_bet = edge_bet;
                                from_id = i;
                                to_id = aeit->first;
                            }
                        }
                }
                gc.remove_edge(from_id,to_id);
                vector<node_set> components;
                id_type num_comps = get_weakly_connected_components(gc,components);
                if(num_comps - component_counter)
                {
                    dendro.push_back(components);
                    component_counter = num_comps;
                }
        }
    }
}

bool ec_comp (const edge_radicchi& lhs,const edge_radicchi& rhs) { return lhs.ecc < rhs.ecc; }

void CDLib::radicchi_et_al_2004(const graph& g, dendrogram& dendro)
{
    vector<node_set> components;
    if(get_weakly_connected_components(g,components) == 1)
    {
        graph gc(g);
        vector<edge_radicchi> ec(g.get_num_edges(),edge_radicchi());
        id_type count = 0;
        for(id_type i=0;i<gc.get_num_nodes(); i++)
                for(adjacent_edges_iterator aeit = gc.out_edges_begin(i); aeit != gc.out_edges_end(i); aeit++)
                    ec[count++] = edge_radicchi(i,aeit->first,edge_clustering_coefficient(gc,i,aeit->first));
        sort(ec.begin(),ec.end(),ec_comp);
        for(id_type i=0;i<ec.size();i++)
        {
                gc.remove_edge(ec[i].from_id,ec[i].to_id);
                vector<node_set> components;
                id_type num_comps = get_weakly_connected_components(gc,components);
                if(num_comps > 2)dendro.push_back(components);
        }
    }
}
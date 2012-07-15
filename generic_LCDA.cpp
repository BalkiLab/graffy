/* 
 * File:   generic_LCDA.cpp
 * Author: prashant
 * 
 * Created on July 11, 2012, 5:51 PM
 */

#include "generic_LCDA.h"
using namespace CDLib;

namespace temp_funcs
{

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
};

using namespace temp_funcs;


/*
Implement

lcd_obj_clauset()

<New obj functions>
 */

//-----------------------------------------------------------------------------------------------------------------------------------------
/*double lcd_obj_clauset(const graph& g, lcd_new& new_vars, id_type node)
{
    
}*/


double lcd_obj_clauset_mod(const graph& g, lcd_new& new_vars, id_type node)
{
    double retval;
    retval = new_vars.I_new/(new_vars.I_new+new_vars.E_new);
    return retval;
}

double lcd_obj_lwp(const graph& g, lcd_new& new_vars, id_type node)
{
    double retval;
    retval = new_vars.I_new/new_vars.E_new;
    return retval;
}

double lcd_obj_vd(const graph& g, lcd_new& new_vars, id_type node)
{
    double retval;
    retval = (2*new_vars.I_new*new_vars.I_new)/(new_vars.n_new*(new_vars.n_new-1)*(new_vars.I_new+new_vars.E_new));
    return retval;
}

double lcd_obj_bagrow(const graph& g, lcd_new& new_vars, id_type node)
{
    double retval;
    retval = 1 - 2/g.get_node_out_degree(node);
    return retval;
}

double lcd_obj_czr(const graph& g, lcd_new& new_vars, id_type node)
{
    double retval;
    retval = (2*new_vars.I_new*new_vars.n_b_new)/(new_vars.E_new*new_vars.n_new);
    return retval;
}
//-----------------------------------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------------------------------

double lcd_delta_obj_clauset(const graph& g, lcd_core& coredata, lcd_new& new_vars, id_type node, double curr_val, node_set& Y)
{
    //compute this in the "different" way
    double delta, x, y, z;
    double T;
    node_set X;

    x = get_neighs_in_set(node, g, coredata.B, X); //X = set of all neighbors of vj in B   //x = No. of edges in T that terminate at vj

    y = g.get_node_out_weight(node) - x; //No. of edges that will be added to T by the agglomeration of vj

    z = 0; //No. of edges that will be removed from T by the agglomeration of vj
    for (node_set::iterator iter1 = X.begin(); iter1 != X.end(); iter1++)
    {
        id_type vk = *iter1;
        if (membership(vk, node, g, coredata.U) == 0) //vk has no neighbor in U other than vj, i.e., vk no longer belongs to B 
        {
            Y.insert(vk); //doubt
            z += no_BI_edges(vk, g, coredata.C, coredata.B) + no_BY_edges(vk, g, Y);
            if (y == 0)
                z = z + g.get_edge_weight(vk, node);
        }
    }

    delta = (x - (curr_val * y) - z * (1 - curr_val)) / (T - z + y); //doubt
    
    return delta;
}


double lcd_delta_obj_clauset_mod(const graph& g, lcd_new& new_vars, id_type node, double curr_val)
{
    double retval;
    retval = lcd_obj_clauset_mod(g,new_vars,node) - curr_val;
    return retval;
}

double lcd_delta_obj_lwp(const graph& g, lcd_new& new_vars, id_type node, double curr_val)
{
    double retval;
    retval = lcd_obj_lwp(g,new_vars,node) - curr_val;
    return retval;
}

double lcd_delta_obj_vd(const graph& g, lcd_new& new_vars, id_type node, double curr_val)
{
    double retval;
    retval = lcd_obj_vd(g,new_vars,node) - curr_val;
    return retval;
}

/*double lcd_delta_obj_bagrow(const graph& g, lcd_new& new_vars, id_type node, double curr_val)
{
    double retval;
    retval = ;
    return retval;
}*/

double lcd_delta_obj_czr(const graph& g, lcd_new& new_vars, id_type node, double curr_val)
{
    double retval;
    retval = lcd_obj_czr(g,new_vars,node) - curr_val;
    return retval;
}
//-----------------------------------------------------------------------------------------------------------------------------------------


bool should_delete_node(const graph& g, lcd_core& coredata, id_type i, lcd_obj_fn deletion_obj_func) {
    ////test deletion here -> instatiates lcd_new() and passes tp obj_func
}

void lcd_deletion(const graph& g, lcd_core& coredata, lcd_obj_fn deletion_obj_func, node_set& nodesDeleted) {
    //For each node in C call should_delete_node()
}

//Return true if you have to stop, the 2 node sets taken as inputs, hold the inserted nodes and the deleted nodes , the vector<double> takes in a bunch of parameters for the stopping function

/*
Implement


lcd_stop_best_p_strong()
lcd_stop_trailing_least_squares()
 * 
lcd_stop_<prashant’s>()
 */

bool lcd_stop_fixed_size_k(const graph& g, lcd_core& coredata, node_set& Q, node_set& deleteQ, vector<double>& params)
{
    if(coredata.C.size()==params[0])
        return 1;
    else
        return 0;
}

bool lcd_stop_no_net_insertions(const graph& g, lcd_core& coredata, node_set& Q, node_set& deleteQ, vector<double>& params)
{
    if(Q.empty())
        return 1;
    else
        return 0;
}

bool lcd_stop_no_change_in_size(const graph& g, lcd_core& coredata, node_set& Q, node_set& deleteQ, vector<double>& params)
{
    if(coredata.C.size()==params[0])
        return 1;
    else
        return 0;
}

bool lcd_stop_no_increase_in_objfunc(const graph& g, lcd_core& coredata, node_set& Q, node_set& deleteQ, vector<double>& params)
{
    if(params[1]>params[0]) //new value of objfunc > old value of objfunc
        return 0;
    else
        return 1;
}



bool should_agglomerate_node(const graph& g, lcd_core& coredata, id_type i) {
    //// ONLY FOR LWPish algorithms - test agglomeration here -> instatiates lcd_new() and passes tp obj_func

}


// last node_set for holding inserted elements

/*
Implement
        lcd_traversal_clauset_lwp_sorted()
        lcd_traversal_clauset_lwp_multipass()
        lcd_traversal_clauset_lwp_backtrack()
        lcd_traversal_clauset_aggl_minmax()
        lcd_traversal_bfs()
        //More traversal functions
 */
void local_community_detection_std(
        const graph& g,
        id_type seed,
        vector<double>& stop_fn_params,
        lcd_traversal_fn traversal,
        bool include_deletion_phase,
        lcd_obj_fn agg_obj_fn,
        lcd_obj_fn del_obj_fn,
        lcd_stop_fn stop_func,
        node_set& comm) {
    lcd_core coredata(g, seed);
    node_set inserted_nodes, deleted_nodes;
    do {
        traversal(g, coredata, agg_obj_fn, inserted_nodes);
        if (include_deletion_phase) lcd_deletion(g, coredata, del_obj_fn, deleted_nodes);
    } while (stop_func(g, coredata, inserted_nodes, deleted_nodes, stop_fn_params));
}

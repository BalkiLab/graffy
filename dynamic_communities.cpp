/* 
 * File:   dynamic_communities.cpp
 * Author: bharath
 * 
 * Created on 10 May, 2012, 11:25 AM
 */

#include "dynamic_communities.h"
using namespace CDLib;

    struct dynamic_lp_temp
    {
        vector<id_type> labels_to_pass;
        vector<id_type> label_changes_to_pass;
    };

void first_snapshot_init(string filepath,bool directed,bool weighted,dynamic_lp_output & output)
{
    output.lp_iters.push_back(0);
    output.lp_times.push_back(0.0);
    output.book_times.push_back(0.0);
    output.lp_labels.push_back(vector<id_type>());
    output.graphs.push_back(graph(directed,weighted));
    read_edgelist(output.graphs[output.graphs.size()-1],filepath,directed,weighted);
    output.lp_num_label_changes.push_back(vector<id_type>(output.graphs[output.graphs.size()-1].get_num_nodes(),0));
}

void subsequent_snapshot_init(string filepath,bool directed,bool weighted,dynamic_lp_output& output, dynamic_lp_temp& temp,bool copy_labels,timer_rt& op_timer)
{
    first_snapshot_init(filepath,directed,weighted,output);
    temp.labels_to_pass.assign(output.graphs[output.graphs.size()-1].get_num_nodes(),0);
    temp.label_changes_to_pass.assign(output.graphs[output.graphs.size()-1].get_num_nodes(),0);
    op_timer.start_clock();
    node_set seen_nodes;
    id_type max_label_value = 0;
    for(id_type i=0;i<output.graphs[output.graphs.size()-2].get_num_nodes();i++)
    {
        id_type current_node_id = output.graphs[output.graphs.size()-1].get_node_id(output.graphs[output.graphs.size()-2].get_node_label(i));
        if(current_node_id < output.graphs[output.graphs.size()-1].get_num_nodes()) //Node has been found
        {
            temp.labels_to_pass[current_node_id] = output.lp_labels[output.lp_labels.size()-2][i];
            max_label_value = (output.lp_labels[output.lp_labels.size()-2][i] > max_label_value) ? output.lp_labels[output.lp_labels.size()-2][i] : max_label_value;
            temp.label_changes_to_pass[current_node_id] = output.lp_num_label_changes[output.lp_num_label_changes.size()-2][i];
            seen_nodes.insert(current_node_id);
        }
    }
    max_label_value++;
    for(id_type i=0;i<output.graphs[output.graphs.size()-1].get_num_nodes();i++)
    {
        if(seen_nodes.find(i) == seen_nodes.end())
        {                                              //Unseen Node
            temp.labels_to_pass[i] =max_label_value + i;
            temp.label_changes_to_pass[i] = 0;
        }
    }
    op_timer.stop_clock();
}

void first_run(timer_rt& op_timer,const string& filepath,bool directed, bool weighted,dynamic_lp_output& output)
{
    first_snapshot_init(filepath,directed,weighted,output);
    op_timer.start_clock();
    lp_dyn_hop_att algoman(output.graphs[output.graphs.size()-1],0,0.0);
    output.lp_iters[output.lp_iters.size()-1] = label_propagation_run(output.graphs[output.graphs.size()-1],output.lp_labels[output.lp_labels.size()-1],10000,&algoman);
    op_timer.stop_clock();
    output.lp_times[output.lp_times.size()-1] = op_timer.run_time();
    for(id_type i=0;i<output.graphs[output.graphs.size()-1].get_num_nodes();i++)output.lp_num_label_changes[output.lp_num_label_changes.size()-1][i] = algoman.get_node_num_lplabel_changes(i);
}

void subsequent_run(timer_rt& op_timer,const string& filepath,bool directed, bool weighted,dynamic_lp_input& input, dynamic_lp_output& output)
{
    dynamic_lp_temp temp;
    subsequent_snapshot_init(filepath,directed,weighted,output,temp,input.copy_labels,op_timer);
    op_timer.start_clock();
    evolutionary_label_propagation algoman(output.graphs[output.graphs.size()-1],0,input.activate,input.copy_labels,input.alpha,output.lp_iters[output.lp_iters.size()-2],temp.label_changes_to_pass,temp.labels_to_pass);
    output.lp_iters[output.lp_iters.size()-1] = label_propagation_run(output.graphs[output.graphs.size()-1],output.lp_labels[output.lp_labels.size()-1],10000,&algoman);
    op_timer.stop_clock();
    output.lp_times[output.lp_times.size()-1] = op_timer.run_time();
    for(id_type i=0;i<output.graphs[output.graphs.size()-1].get_num_nodes();i++)output.lp_num_label_changes[output.lp_num_label_changes.size()-1][i] = algoman.get_node_num_lplabel_changes(i);
}

double CDLib::evolutionary_label_propagation_edgelists(const string& snapshot_filepath,bool directed, bool weighted,dynamic_lp_input& input, dynamic_lp_output& output)
{
    ifstream ifs;
    ifs.open(snapshot_filepath.c_str());
    timer_rt op_timer;
    if(ifs.is_open())
    {
        string filepath;
        ifs >> filepath;
        first_run(op_timer,filepath,directed,weighted,output);
        while(!ifs.eof()) 
        {
            ifs >> filepath;
            subsequent_run(op_timer,filepath,directed,weighted,input,output);
        }
    }
    return op_timer.total_time();
}

void pre_run(ifstream& ifs, bool directed, bool weighted,dynamic_lp_output& output)
{
        output.lp_iters.push_back(0);
        output.lp_times.push_back(0.0);
        output.book_times.push_back(0.0);
        output.lp_labels.push_back(vector<id_type>());
        output.graphs.push_back(graph(directed,weighted));
        string filepath;
        ifs >> filepath;
        read_edgelist(output.graphs[output.graphs.size()-1],filepath,directed,weighted);
}

double prepare_labels_to_pass(dynamic_lp_output& output,vector<id_type>::const_iterator begin_it,vector<id_type>& temp_changes)
{
    timer_rt lp_timer;
    lp_timer.start_clock();
    temp_changes.assign(output.graphs[output.graphs.size()-1].get_num_nodes(),0);
    for(id_type i=0;i<output.graphs[output.graphs.size()-2].get_num_nodes();i++)
    {
        id_type current_node_id = output.graphs[output.graphs.size()-1].get_node_id(output.graphs[output.graphs.size()-2].get_node_label(i));
        if(current_node_id < output.graphs[output.graphs.size()-1].get_num_nodes())
            temp_changes[current_node_id] = *(begin_it + i);
    }
    lp_timer.stop_clock();
    return lp_timer.run_time();
}

void new_subsequent_run(ifstream& ifs,bool directed,bool weighted,double alpha,dynamic_lp_output& output,vector<id_type>& temp_changes)
{
    timer_rt lp_timer;
    lp_timer.start_clock();
    evol_label_prop_new algoman(output.graphs[output.graphs.size()-1],0,alpha,output.lp_iters[output.lp_iters.size()-2], temp_changes.begin());
    output.lp_iters[output.lp_iters.size()-1] = label_propagation_run(output.graphs[output.graphs.size()-1],output.lp_labels[output.lp_labels.size()-1],10000,&algoman);
    lp_timer.stop_clock();
    output.lp_times[output.lp_times.size()-1] = lp_timer.run_time();
    pre_run(ifs,directed,weighted,output);
    output.book_times[output.book_times.size()-1] = prepare_labels_to_pass(output,algoman.label_changes_begin(),temp_changes);
}

void CDLib::new_evolutionary_label_propagation_edgelists(const string& snapshot_filepath,bool directed, bool weighted,double alpha,bool basic, dynamic_lp_output& output)
{
    ifstream ifs;
    ifs.open(snapshot_filepath.c_str());

    if(ifs.is_open())
    {
        id_type num_snapshots = 0;
        vector<id_type> temp_changes;
        while(!ifs.eof()) 
        {
            pre_run(ifs,directed,weighted,output);
            timer_rt lp_timer;
            if(!basic)
            {
                if(!num_snapshots)
                {
                        lp_timer.start_clock();
                        lp_track_changes algoman(output.graphs[output.graphs.size()-1],0);
                        output.lp_iters[output.lp_iters.size()-1] = label_propagation_run(output.graphs[output.graphs.size()-1],output.lp_labels[output.lp_labels.size()-1],10000,&algoman);
                        lp_timer.stop_clock();
                        output.lp_times[output.lp_times.size()-1] = lp_timer.run_time();
                        pre_run(ifs,directed,weighted,output);
                        temp_changes.assign(output.graphs[output.graphs.size()-1].get_num_nodes(),0);
                        output.book_times[output.book_times.size()-1] = prepare_labels_to_pass(output,algoman.label_changes_begin(),temp_changes);
                }
                else new_subsequent_run(ifs,directed,weighted,alpha,output,temp_changes);
            }
            else
            {
                lp_timer.start_clock();
                lp_raghavan_2007 algoman(output.graphs[output.graphs.size()-1],0);
                output.lp_iters[output.lp_iters.size()-1] = label_propagation_run(output.graphs[output.graphs.size()-1],output.lp_labels[output.lp_labels.size()-1],10000,&algoman);
                lp_timer.stop_clock();
                output.lp_times[output.lp_times.size()-1] = lp_timer.run_time();
            }
        }
    }
}
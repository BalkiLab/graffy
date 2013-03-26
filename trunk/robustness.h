/*
 * File:   robustness.h
 * Author: sudip
 *
 * Created on 26 March, 2013, 12:58 PM
 */

#ifndef ROBUSTNESS_H
#define	ROBUSTNESS_H

#include "typedefs.h"
#include "statistics.h"

namespace CDLib {

    struct stability {
        double efficiency, fraction_in_lcc;
        double degree_entropy, path_entropy, connectivity_entropy;
        double assortativity, reachability, regularity;
        double avg_clustering_coefficient, avg_path_length_lcc;
        id_type size; // This represents the number of parameters.

        stability() {
            // Default Indexing: 0. efficiency, 1. fraction_in_lcc, 2. assortativity, 3. degree_entropy, 4. path_entropy,
            //                   5. connectivity_entropy 6. reachability, 7. regularity, 8. avg_clustering_coefficient, 9. avg_path_length_lcc
            init();
        }

        stability(double value) {
            init();
            assign(value);
        }

        stability(const graph & g) {
            init();
            evaluations(g);
        }

        stability(const vector<double>& params) {
            init();
            load_from_vector(params);
        }

        void clear() {
            efficiency = 0;
            fraction_in_lcc = 0;
            assortativity = 0;
            degree_entropy = 0;
            path_entropy = 0;
            connectivity_entropy = 0;
            reachability = 0;
            regularity = 0;
            avg_clustering_coefficient = 0;
            avg_path_length_lcc = 0;
        }

        bool assign(const vector<double>& evaluated_values) {
            init();
            return load_from_vector(evaluated_values);
        }

        void assign(const double value) {
            efficiency = value;
            fraction_in_lcc = value;
            assortativity = value;
            degree_entropy = value;
            path_entropy = value;
            connectivity_entropy = value;
            reachability = value;
            regularity = value;
            avg_clustering_coefficient = value;
            avg_path_length_lcc = value;
        }

        inline void evaluations(const graph & g) {
            /* The code works properly only for undirected graph. */
            vector<double> reachability_distt(g.get_num_nodes(), 0);
            vector<double> total_distance(g.get_num_nodes(), 0.0);
            vector<double> degree_distt(g.get_num_nodes(), 0);
            id_type all_sources = 0, number_of_isolates = 0;
            double loc_eff = 0, loc_cc = 0, loc_reg = 0;
            double excess_mean = 0, excess_variance = 0, loc_ass = 0, max_degree = 0, deg_mean = 0;
            vector<bool> has_node_visited(g.get_num_nodes(), 0);
            vector<node_set> components;
            vector<double>avg_nei_deg(g.get_num_nodes(), 0);
            double normalization = g.get_num_nodes() * (g.get_num_nodes() - 1);
            double norm_edges = 2 * g.get_num_edges();
#pragma omp parallel for shared(g,total_distance,reachability_distt,has_node_visited,components,avg_nei_deg,degree_distt) reduction(+:all_sources,loc_eff,loc_cc, loc_reg, excess_mean, excess_variance,number_of_isolates,deg_mean)
            for (id_type i = 0; i < g.get_num_nodes(); i++) {
                deg_mean += g.get_node_out_degree(i);
                if (!g.get_node_out_degree(i))
                    number_of_isolates++;
                /* This is the BFS part from each node. */
                node_set current_comm;
                current_comm.insert(i);
                vector<double> distances(g.get_num_nodes(), numeric_limits<double>::infinity());
                queue<id_type> q_bfs;
                q_bfs.push(i);
                distances[i] = 0;
                while (!q_bfs.empty()) {
                    id_type current = q_bfs.front();
                    q_bfs.pop();
                    for (adjacent_edges_iterator aeit = g.out_edges_begin(current); aeit != g.out_edges_end(current); aeit++) {
                        if (distances[aeit->first] == numeric_limits<double>::infinity() && current != aeit->first) {
                            distances[aeit->first] = distances[current] + 1;
                            q_bfs.push(aeit->first);
                            loc_eff += (1 / distances[aeit->first]);
                            total_distance[i] += distances[aeit->first];
                            reachability_distt[i]++;
                            all_sources++;
                            current_comm.insert(aeit->first);
                            /* This will iterate over all edges. */
                            if (distances[aeit->first] == 1) {
                                loc_reg += 1 / (1 + abs(g.get_node_out_degree(i) - g.get_node_in_degree(aeit->first)));
                                avg_nei_deg[i] += g.get_node_in_degree(aeit->first);
                            }
                        }
                    }
                }
                /* The BFS part from each node ends here. */
                loc_cc += node_clustering_coefficient(g, i);
                if (g.get_node_out_degree(i) > 0)
                    avg_nei_deg[i] = (avg_nei_deg[i] / g.get_node_out_degree(i)) - 1;
#pragma omp critical
                {
                    max_degree = (max_degree < g.get_node_out_degree(i)) ? g.get_node_out_degree(i) : max_degree;
                    degree_distt[g.get_node_out_degree(i)]++;
                    if (!has_node_visited[i]) {
                        components.push_back(current_comm);
                        has_node_visited[i] = true;
                        for (node_set::iterator it = current_comm.begin(); it != current_comm.end(); it++) {
                            has_node_visited[*it] = true;
                        }
                    }
                }

            }
            /* Getting degree distribution and excess degree details */
            degree_distt.resize(max_degree + 1);
            if (degree_distt.size())
                degree_distt[0] = 0;
            double denom = g.get_num_nodes() - number_of_isolates;
            for (id_type i = 0; i < degree_distt.size(); i++)
                degree_distt[i] /= denom;
            deg_mean /= denom;
            excess_degree(degree_distt, deg_mean, excess_mean, excess_variance);
            /* This consider only when there is atleast one edge. */
            id_type lcc = 0;
            if (all_sources > 0) {
                vector<id_type> index_max;
                for (id_type i = 0; i < g.get_num_nodes(); i++) {
                    reachability += reachability_distt[i];
                    double reach_entropy = reachability_distt[i] / all_sources;
                    if (reach_entropy != 0) {
                        path_entropy -= reach_entropy * (log(reach_entropy) / log(2));
                    }
                    if (i < components.size()) {
                        if (components[i].size() >= lcc) {
                            if (components[i].size() > lcc) {
                                index_max.clear();
                                lcc = components[i].size();
                                index_max.push_back(i);
                            } else {
                                index_max.push_back(i);
                            }
                        }
                    }
                    /* Getting degree histogram for degree distribution. */
                    double degree = (double) g.get_node_out_degree(i);
                    if (degree > 1) {
                        loc_ass += degree * (degree - 1) * (avg_nei_deg[i] - excess_mean);
                        double prob = degree / norm_edges;
                        connectivity_entropy += prob * log2(prob);
                    }
                }
                connectivity_entropy *= -1;
                /* Any measures for LCC can be calculated in the below loop */
                for (id_type i = 0; i < index_max.size(); i++) {
                    for (node_set::iterator it = components[index_max[i]].begin(); it != components[index_max[i]].end(); it++) {
                        avg_path_length_lcc += total_distance[*it];
                    }
                }
                if (index_max.size() > 0) {
                    avg_path_length_lcc /= (lcc * (lcc - 1) * index_max.size());
                }
                /* 2 x |E| normalizations. Normalizing degree histogram to get degree distribution. */
                regularity = loc_reg / norm_edges;
                if (excess_variance != 0)
                    assortativity = loc_ass / (norm_edges * excess_variance);
                else
                    assortativity = 1;
            }
            /* This consider only when there is atleast one node. */
            if (g.get_num_nodes()) {
                avg_clustering_coefficient = loc_cc / g.get_num_nodes();
                efficiency = loc_eff / normalization;
                reachability /= normalization;
                fraction_in_lcc = (double) lcc / g.get_num_nodes();
                degree_entropy = distribution_entropy(degree_distt);
                connectivity_entropy /= log2(g.get_num_nodes());
            }
        }

        template <typename T>
                void copy_to_vector(vector<T>& parameters) {
            parameters.clear();
            parameters.push_back(efficiency);
            parameters.push_back(fraction_in_lcc);
            parameters.push_back(assortativity);
            parameters.push_back(degree_entropy);
            parameters.push_back(path_entropy);
            parameters.push_back(connectivity_entropy);
            parameters.push_back(reachability);
            parameters.push_back(regularity);
            parameters.push_back(avg_clustering_coefficient);
            parameters.push_back(avg_path_length_lcc);
        }

        void copy(struct stability & get) {
            size = get.size;
            efficiency = get.efficiency;
            fraction_in_lcc = get.fraction_in_lcc;
            assortativity = get.assortativity;
            degree_entropy = get.degree_entropy;
            path_entropy = get.path_entropy;
            connectivity_entropy = get.connectivity_entropy;
            reachability = get.reachability;
            regularity = get.regularity;
            avg_clustering_coefficient = get.avg_clustering_coefficient;
            avg_path_length_lcc = get.avg_path_length_lcc;
        }

        template <typename T>
                bool add_with_current(vector<T>& parameters) {
            if (parameters.size() != size)
                return false;
            parameters[0] += efficiency;
            parameters[1] += fraction_in_lcc;
            parameters[2] += assortativity;
            parameters[3] += degree_entropy;
            parameters[4] += path_entropy;
            parameters[5] += connectivity_entropy;
            parameters[6] += reachability;
            parameters[7] += regularity;
            parameters[8] += avg_clustering_coefficient;
            parameters[9] += avg_path_length_lcc;
            return true;
        }

        inline bool load_from_vector(const vector<double>& params) {
            if (params.size() != size)
                return false;
            efficiency = params[0];
            fraction_in_lcc = params[1];
            assortativity = params[2];
            degree_entropy = params[3];
            path_entropy = params[4];
            connectivity_entropy = params[5];
            reachability = params[6];
            regularity = params[7];
            avg_clustering_coefficient = params[8];
            avg_path_length_lcc = params[9];
            return true;
        }

        inline bool load_from_string(const string & line) {
            vector<double> params;
            split(line, params);
            return load_from_vector(params);
        }

        bool load_from_file(const string & infile) {
            ifstream ifs(infile);
            string line;
            getline(ifs, line);
            return load_from_string(line);
        }

        string save_to_string() {
            ostringstream oss;
            oss.precision(6);
            oss << fixed << efficiency << "\t" << fixed << fraction_in_lcc << "\t" << fixed << assortativity << "\t" << fixed << degree_entropy << "\t" << fixed << path_entropy << "\t";
            oss << fixed << connectivity_entropy << "\t" << fixed << reachability << "\t" << fixed << regularity << "\t" << fixed << avg_clustering_coefficient << "\t" << fixed << avg_path_length_lcc;
            return oss.str();
        }

        bool save_to_file(const string & outfile) {
            ofstream ofs;
            ofs.open(outfile.c_str());
            if (ofs.is_open()) {
                ofs << save_to_string();
                return true;
            }
            return false;
        }
    private:

        void init() {
            size = 10;
            clear();
        }

        void excess_degree(vector<double>& degree_distt, double dmean, double& ex_mean, double& ex_variance) {
            vector<double> excess_degree_distt(degree_distt.size() - 1, 0);
            ex_mean = 0;
            ex_variance = 0;
            for (id_type i = 0; i < excess_degree_distt.size(); i++) {
                excess_degree_distt[i] = ((i + 1) * degree_distt[i + 1]) / dmean;
                double imean = i * excess_degree_distt[i];
                ex_mean += imean;
                ex_variance += i*imean;
            }
            ex_variance -= (ex_mean * ex_mean);
        }
    };

    struct stable {
        vector<string> label;
        vector<double> value;
        struct stability init;
        vector<struct stability> removed;
        id_type size; // This represents the number of times it is re-evaluated after graph change operation.

        stable() {
            size = 0;
        }

        inline bool is_empty() {
            return (size > 0) ? false : true;
        }

        inline void clear_init() {
            init.clear();
        }

        inline void clear_evaluations() {
            removed.clear();
            size = 0;
        }

        inline void clear_non_init() {
            label.clear();
            value.clear();
            clear_evaluations();
        }

        void clear() {
            clear_init();
            clear_non_init();
        }

        bool copy_all(const struct stable& get, id_type limit) {
            if (get.size < limit)
                return false;
            clear_non_init();
            copy(get.label.begin(), get.label.begin() + limit, back_inserter(label));
            copy(get.value.begin(), get.value.begin() + limit, back_inserter(value));
            copy(get.removed.begin(), get.removed.begin() + limit, back_inserter(removed));
            size = limit;
            return true;
        }

        bool copy_evaluations(const struct stable& get, id_type limit) {
            if (get.size < limit)
                return false;
            clear_evaluations();
            copy(get.removed.begin(), get.removed.begin() + limit, back_inserter(removed));
            size = limit;
            return true;
        }

        inline bool set(id_type index, vector<double>& evaluated_values) {
            if (index >= size)
                return false;
            removed[index].assign(evaluated_values);
            return true;
        }

        inline void assign(id_type i, double value) {
            clear_evaluations();
            size = i;
            removed.assign(size, stability(value));
        }

        void assign_evaluations(vector< vector<double> >& evaluated_values) {
            clear_evaluations();
            size = evaluated_values.size();
            removed.assign(size, stability());
            for (id_type i = 0; i < evaluated_values.size(); i++) {
                set(i, evaluated_values[i]);
            }
        }

        inline void assign(id_type i) {
            assign(i, 0);
            label.clear();
            value.clear();
            label.assign(i, "");
            value.assign(i, 0);
        }

        inline bool insert_from_previous(const id_type previous, const id_type current) {
            if ((previous >= size) || (current >= size))
                return false;
            removed[current].copy(removed[previous]);
            return true;
        }

        inline bool evaluate_and_insert(const id_type i, const string& ilabel, const double& ivalue, const graph & g) {
            if (i >= size)
                return false;
            label[i] = ilabel;
            value[i] = ivalue;
            evaluate_and_insert(i, g);
            return true;
        }

        inline bool evaluate_and_insert(const id_type i, const pair<string, double>& element, const graph & g) {
            return evaluate_and_insert(i, element.first, element.second, g);
        }

        inline bool evaluate_and_insert(id_type index, const graph & g) {
            if (index >= size)
                return false;
            removed[index].evaluations(g);
            return true;
        }

        inline bool evaluate_and_insert(const pair<string, double>& element, const graph & g) {
            return evaluate_and_insert(element.first, element.second, g);
        }

        bool evaluate_and_insert(const string& ilabel, const double& ivalue, const graph & g) {
            if (ivalue > 0) {
                value.push_back(ivalue);
                if (!(ilabel == ""))
                    label.push_back(string(ilabel));
            }
            removed.push_back(stability(g));
            size++;
            return size;
        }

        inline void graph_init(const graph & g) {
            init.evaluations(g);
        }

        inline void copy_graph_init(vector<double>& vec_init) {
            init.copy_to_vector(vec_init);
        }

        inline bool insert_graph_init(const vector<double>& vec_init) {
            return init.assign(vec_init);
        }

        id_type load_from_vector(const vector< vector<double> >& parameters) {
            if (parameters.size() == 0)
                return 0;
            else {
                clear(); // Need to put clear in rest of the load functions.
                init.load_from_vector(parameters[0]);
                for (id_type i = 1; i < parameters.size(); i++) {
                    removed.push_back(stability(parameters[i]));
                    size++;
                }
            }
            return size;
        }

        id_type load_from_vector(const vector< vector<double> >& parameters, const vector<double>& vals, const vector<string>& labs) {
            if (!(parameters.size() && vals.size() && labs.size()))
                return 0;
            else {
                if (vals.size() != labs.size())
                    return 0;
                id_type size = load_from_vector(parameters); // This also clears the old values.
                copy(vals.begin(), vals.end(), back_inserter(value));
                copy(labs.begin(), labs.end(), back_inserter(label));
                return size;
            }
        }

        id_type load_from_vector(const vector< vector<double> >& parameters, const vector<double>& vals) {
            if (!(parameters.size() && vals.size()))
                return 0;
            else {
                id_type size = load_from_vector(parameters); // This also clears the old values.
                value.assign(vals.size(), 0);
                for (id_type i = 0; i < vals.size(); i++) {
                    value[i] = vals[i];
                    label.push_back(T2str(i)); // The value ID's are assigned as labels.
                }
                return size;
            }
        }

        id_type load_from_string(const string & iline) {
            /* LINE 1: L label1 label2 label3 ....*/
            /* LINE 2: V value1 value2 value3 ....*/
            /* PARAM 0: init */
            clear(); // Clears all previous values
            bool flag_label = false;
            istringstream iss(iline);
            string line;
            getline(iss, line);
            while (line[0] == 'L' || line[0] == 'V') {
                if (line[0] == 'L') {
                    line.erase(0, 2);
                    split(line, label);
                    flag_label = true;
                } else {
                    line.erase(0, 2);
                    split(line, value);
                    if (!flag_label) {
                        for (id_type i = 0; i < value.size(); i++)
                            label.push_back(T2str(i));
                    }
                }
                line.clear();
                getline(iss, line);
            }
            init.load_from_string(line);
            while (!iss.eof()) {
                line.clear();
                getline(iss, line);
                vector<double> params;
                split(line, params);
                if (params.size() == init.size) {
                    removed.push_back(stability(params));
                    size++;
                }
            }
            return size;
        }

        id_type load_from_file(const string & infile) {
            ifstream ifs;
            ifs.open(infile.c_str());
            if (ifs.is_open()) {
                string params_to_fill;
                while (!ifs.eof()) {
                    string line;
                    getline(ifs, line);
                    params_to_fill.append(line);
                }
                return load_from_string(params_to_fill);
            }
            return 0;
        }

        string save_to_string() {
            ostringstream oss;
            if (label.size()) {
                oss << "L " << label[0]; // Starting with 'L ' specifies that the line contains label
                for (id_type i = 1; i < label.size(); i++)
                    oss << "\t" << label[i];
                oss << endl;
            }
            if (value.size()) {
                oss << "V " << value[0]; // Starting with 'V ' specifies that the line contains value
                for (id_type i = 1; i < value.size(); i++)
                    oss << "\t" << value[i];
                oss << endl;
            }
            oss << init.save_to_string() << endl;
            for (id_type i = 0; i < size; i++)
                oss << removed[i].save_to_string() << endl;
            return oss.str();
        }

        bool save_to_file(const string & outfile) {
            ofstream ofs;
            ofs.open(outfile.c_str());
            if (ofs.is_open()) {
                ofs << save_to_string();
                return true;
            }
            return false;
        }
    };

    struct node_attack_strategy {
        string graph_name;
        struct stable degree, current_degree, betweenness, current_betweenness;
        struct stable random, efficiency_score, current_efficiency_score;
        struct stable reduced_degree, reduced_betweenness, reduced_efficiency_score, reduced_random;
        struct stable reduced_current_degree, reduced_current_betweenness, reduced_current_efficiency_score;

        node_attack_strategy() {

            initialized = false;
        }

        node_attack_strategy(graph & g) {

            graph_init(g);
        }

        void perform_all_attacks(graph & g) {

            degree_attack(g);
            current_degree_attack(g);
            betweenness_attack(g);
            current_betweenness_attack(g);
            random_attack(g);
            efficiency_attack(g);
            current_efficiency_attack(g);
        }

        void clear() {
            initialized = false;
            graph_name.clear();
            degree.clear();
            current_degree.clear();
            betweenness.clear();
            current_betweenness.clear();
            random.clear();
            efficiency_score.clear();
            current_efficiency_score.clear();
            reduced_degree.clear();
            reduced_current_degree.clear();
            reduced_betweenness.clear();
            reduced_current_betweenness.clear();
            reduced_random.clear();
            reduced_efficiency_score.clear();
            reduced_current_efficiency_score.clear();
        }

        void graph_init(const graph & g) {
            graph_name = g.get_graph_name();
            degree.graph_init(g);
            vector<double> init;
            degree.copy_graph_init(init);
            current_degree.insert_graph_init(init);
            betweenness.insert_graph_init(init);
            current_betweenness.insert_graph_init(init);
            random.insert_graph_init(init);
            efficiency_score.insert_graph_init(init);
            current_efficiency_score.insert_graph_init(init);
            reduced_degree.insert_graph_init(init);
            reduced_current_degree.insert_graph_init(init);
            reduced_betweenness.insert_graph_init(init);
            reduced_current_betweenness.insert_graph_init(init);
            reduced_random.insert_graph_init(init);
            reduced_efficiency_score.insert_graph_init(init);
            reduced_current_efficiency_score.insert_graph_init(init);
            initialized = true;
        }

        void perform_attack(const graph& g, void (*get_params)(const graph& g, vector<double>& params), struct stable & attack_type) {
            if (!initialized)
                graph_init(g);
            attack_type.assign(g.get_num_nodes());
            graph tmpg(g);
            vector< pair<id_type, double> > rearrange;
            vector<double> pre_calc;
            get_params(g, pre_calc);
            for (id_type i = 0; i < g.get_num_nodes(); i++)
                rearrange.push_back(make_pair(i, pre_calc[i]));
            sort(rearrange, 0);
            for (id_type i = 0; i < rearrange.size(); i++) {
                string label = g.get_node_label(rearrange[i].first);
                tmpg.isolate_node(label);
                attack_type.evaluate_and_insert(i, label, (double) rearrange[i].second, tmpg);
            }
        }

        void perform_attack(const graph& g, void (*get_params)(const graph& g, vector<double>& params), id_type attack_size, struct stable& pre_calc, struct stable & attack_type) {
            if (!initialized)
                graph_init(g);
            attack_size = (attack_size > g.get_num_nodes()) ? g.get_num_nodes() : attack_size;
            if (pre_calc.size >= attack_size) {
                attack_type.copy_all(pre_calc, attack_size);
            } else {
                attack_type.assign(attack_size);
                graph tmpg(g);
                vector< pair<id_type, double> > rearrange;
                vector<double> params;
                get_params(g, params);
                for (id_type i = 0; i < g.get_num_nodes(); i++)
                    rearrange.push_back(make_pair(i, params[i]));
                sort(rearrange, 0);
                for (id_type i = 0; i < attack_size; i++) {
                    string label = g.get_node_label(rearrange[i].first);
                    tmpg.isolate_node(label);
                    attack_type.evaluate_and_insert(i, label, (double) rearrange[i].second, tmpg);
                }
            }
        }

        void perform_current_attack(const graph& g, pair<string, double> (*get_max_node)(const graph & g), id_type attack_size, struct stable & attack_type) {
            if (!initialized)
                graph_init(g);
            attack_type.clear_non_init();
            graph tmpg(g);
            pair<string, double> element;
            id_type i = 0;
            attack_size = (attack_size > g.get_num_nodes()) ? g.get_num_nodes() : attack_size;
            for (i = 0; i < attack_size; i++) {
                element = get_max_node(tmpg);
                tmpg.isolate_node(element.first);
                attack_type.evaluate_and_insert(element, tmpg);
                if (element.second <= 0)
                    break;
            }
            for (id_type k = i + 1; k < attack_size; k++) {
                attack_type.insert_from_previous(i, k);
            }
        }

        void perform_reduced_current_attack(const graph& g, pair<string, double> (*get_max_node)(const graph & g, const node_set_string & members), id_type attack_size, node_set_string& members, struct stable& pre_calc, struct stable & attack_type) {
            if (!initialized)
                graph_init(g);
            attack_type.clear_non_init();
            node_set_string removal_set(members);
            graph tmpg(g);
            pair<string, double> element;
            id_type i = 0;
            assert(removal_set.size() <= g.get_num_nodes());
            attack_size = (attack_size > removal_set.size()) ? removal_set.size() : attack_size;
            if (pre_calc.size) {
                node_set_string::iterator it = removal_set.find(pre_calc.label[i]);
                while ((it != removal_set.end()) && (i < attack_size)) {
                    removal_set.erase(pre_calc.label[i]);
                    tmpg.isolate_node(pre_calc.label[i]);
                    i++;
                    it = removal_set.find(pre_calc.label[i]);
                }
                attack_type.copy_all(pre_calc, i);
            }
            for (; i < attack_size; i++) {
                element = get_max_node(tmpg, removal_set);
                tmpg.isolate_node(element.first);
                attack_type.evaluate_and_insert(element, tmpg);
                if (element.second <= 0)
                    break;
                removal_set.erase(element.first);
            }
            for (id_type k = i + 1; k < attack_size; k++) {
                attack_type.insert_from_previous(i, k);
            }
        }

        void perform_random_attack(graph& g, id_type monte_carlo, id_type attack_size, struct stable& pre_calc, struct stable & attack_type) {
            if (!initialized)
                graph_init(g);
            attack_size = (attack_size > g.get_num_nodes()) ? g.get_num_nodes() : attack_size;
            if (attack_size <= pre_calc.size) {
                attack_type.copy_evaluations(pre_calc, attack_size);
            } else {
                struct stability g_stab;
                vector< vector<double> > evals(attack_size, vector<double>(g_stab.size, 0));
                vector<string> nodes;
                for (id_type i = 0; i < g.get_num_nodes(); i++)
                    nodes.push_back(g.get_node_label(i));
                for (id_type mc = 0; mc < monte_carlo; mc++) {
                    graph tmpg(g);
                    random_shuffle(nodes.begin(), nodes.end());
                    for (id_type j = 0; j < attack_size; j++) {
                        tmpg.isolate_node(nodes[j]);
                        g_stab.evaluations(tmpg);
                        g_stab.add_with_current(evals[j]);
                    }
                }
                for (id_type i = 0; i < attack_size; i++) {
                    for (id_type j = 0; j < g_stab.size; j++) {
                        evals[i][j] /= monte_carlo;
                    }
                }
                attack_type.assign_evaluations(evals);
            }
        }

        void degree_attack(const graph & g) {
            perform_attack(g, degree_centralities, degree);
        }

        void reduced_degree_attack(const graph& g, id_type attack_size) {
            perform_attack(g, degree_centralities, attack_size, degree, reduced_degree);
        }

        void current_degree_attack(const graph & g) {
            perform_current_attack(g, get_max_degree_node, g.get_num_nodes(), current_degree);
        }

        void reduced_current_degree_attack(const graph & g, node_set_string & members) {
            perform_reduced_current_attack(g, get_max_degree_node, members.size(), members, current_degree, reduced_current_degree);
        }

        void betweenness_attack(graph & g) {
            perform_attack(g, betweeness_centralities, betweenness);
        }

        void reduced_betweenness_attack(graph & g, id_type attack_size) {
            perform_attack(g, betweeness_centralities, attack_size, betweenness, reduced_betweenness);
        }

        void current_betweenness_attack(const graph & g) {
            perform_current_attack(g, get_max_betweenness_node, g.get_num_nodes(), current_betweenness);
        }

        void reduced_current_betweenness_attack(const graph & g, node_set_string & members) {
            perform_reduced_current_attack(g, get_max_betweenness_node, members.size(), members, current_betweenness, reduced_current_betweenness);
        }

        void efficiency_attack(graph & g) {
            perform_attack(g, efficiency_centralities, efficiency_score);
        }

        void reduced_efficiency_attack(graph & g, id_type attack_size) {
            perform_attack(g, efficiency_centralities, attack_size, efficiency_score, reduced_efficiency_score);
        }

        void current_efficiency_attack(graph & g) {
            perform_current_attack(g, get_max_efficiency_centrality_node, g.get_num_nodes(), current_efficiency_score);
        }

        void reduced_current_efficiency_attack(graph & g, node_set_string & members) {
            perform_reduced_current_attack(g, get_max_efficiency_centrality_node, members.size(), members, current_efficiency_score, reduced_current_efficiency_score);
        }

        inline void reduced_random_attack(graph& g, id_type monte_carlo, id_type attack_size) {
            perform_random_attack(g, monte_carlo, attack_size, random, reduced_random);
        }

        inline void random_attack(graph& g, id_type monte_carlo) {
            random.clear_evaluations();
            perform_random_attack(g, monte_carlo, g.get_num_nodes(), random, random);
        }

        inline void random_attack(graph & g) {

            id_type monte_carlo = 100; // This is the default value of monte carlo simulation set.
            random_attack(g, monte_carlo);
        }

        string save_to_string() {
            ostringstream oss;
            oss << "G " << graph_name << endl;
            if (degree.size) {
                oss << "A DEGREE" << endl;
                oss << degree.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (current_degree.size) {
                oss << "A CURRENT_DEGREE" << endl;
                oss << current_degree.save_to_string();
                oss << "E" << endl;
            }
            if (betweenness.size) {
                oss << "A BETWEENNESS" << endl;
                oss << betweenness.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (current_betweenness.size) {
                oss << "A CURRENT_BETWEENNESS" << endl;
                oss << current_betweenness.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (random.size) {
                oss << "A RANDOM" << endl;
                oss << random.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (efficiency_score.size) {
                oss << "A EFFICIENCY_SCORE" << endl;
                oss << efficiency_score.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (current_efficiency_score.size) {
                oss << "A CURRENT_EFFICIENCY_SCORE" << endl;
                oss << current_efficiency_score.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (reduced_degree.size) {
                oss << "A REDUCED_DEGREE" << endl;
                oss << reduced_degree.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (reduced_current_degree.size) {
                oss << "A REDUCED_CURRENT_DEGREE" << endl;
                oss << reduced_current_degree.save_to_string();
                oss << "E" << endl;
            }
            if (reduced_betweenness.size) {
                oss << "A REDUCED_BETWEENNESS" << endl;
                oss << reduced_betweenness.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (reduced_current_betweenness.size) {
                oss << "A REDUCED_CURRENT_BETWEENNESS" << endl;
                oss << reduced_current_betweenness.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (reduced_random.size) {
                oss << "A REDUCED_RANDOM" << endl;
                oss << reduced_random.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (reduced_efficiency_score.size) {
                oss << "A REDUCED_EFFICIENCY_SCORE" << endl;
                oss << reduced_efficiency_score.save_to_string() << endl;
                oss << "E" << endl;
            }
            if (reduced_current_efficiency_score.size) {
                oss << "A REDUCED_CURRENT_EFFICIENCY_SCORE" << endl;
                oss << reduced_current_efficiency_score.save_to_string() << endl;
                oss << "E" << endl;
            }
            return oss.str();
        }

        bool save_to_file(const string & outfile) {
            ofstream ofs;
            ofs.open(outfile.c_str());
            if (ofs.is_open()) {
                ofs << save_to_string();
                return true;
            }
            return false;
        }

        bool load_from_file(const string & filepath) {
            ifstream ifs;
            ifs.open(filepath.c_str());
            if (ifs.is_open()) {
                clear();
                ostringstream block;
                string line;
                getline(ifs, line);
                copy(line.begin() + 2, line.end(), back_inserter(graph_name));
                while (!ifs.eof()) {
                    line.clear();
                    getline(ifs, line);
                    if (line[0] == 'A') {
                        string attack(line.begin() + 2, line.end());
                        line.clear();
                        getline(ifs, line);
                        block.str(string());
                        while (line[0] != 'E') {
                            block << line << endl;
                            line.clear();
                            getline(ifs, line);
                        }
                        if (attack == "DEGREE")
                            degree.load_from_string(block.str());
                        else if (attack == "CURRENT_DEGREE")
                            current_degree.load_from_string(block.str());
                        else if (attack == "BETWEENNESS")
                            betweenness.load_from_string(block.str());
                        else if (attack == "CURRENT_BETWEENNESS")
                            current_betweenness.load_from_string(block.str());
                        else if (attack == "RANDOM")
                            random.load_from_string(block.str());
                        else if (attack == "EFFICIENCY_SCORE")
                            efficiency_score.load_from_string(block.str());
                        else if (attack == "CURRENT_EFFICIENCY_SCORE")
                            current_efficiency_score.load_from_string(block.str());
                        else if (attack == "REDUCED_DEGREE")
                            reduced_degree.load_from_string(block.str());
                        else if (attack == "REDUCED_CURRENT_DEGREE")
                            reduced_current_degree.load_from_string(block.str());
                        else if (attack == "REDUCED_BETWEENNESS")
                            reduced_betweenness.load_from_string(block.str());
                        else if (attack == "REDUCED_CURRENT_BETWEENNESS")
                            reduced_current_betweenness.load_from_string(block.str());
                        else if (attack == "REDUCED_RANDOM")
                            reduced_random.load_from_string(block.str());
                        else if (attack == "REDUCED_EFFICIENCY_SCORE")
                            reduced_efficiency_score.load_from_string(block.str());
                        else if (attack == "REDUCED_CURRENT_EFFICIENCY_SCORE")
                            reduced_current_efficiency_score.load_from_string(block.str());
                    }
                    initialized = true;
                }
                return true;
            }
            return false;
        }

    private:
        bool initialized;

    };

    struct edge_attack_strategy {
        struct stable betweenness, current_betweenness;
        struct stable random, eff_score, current_eff_score;
    };
};

#endif	/* ROBUSTNESS_H */


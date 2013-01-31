#include "datastructures.h"
using namespace std;
using namespace CDLib;

id_type bidirectional_label_map::size() const {
    return fm_labels.size();
}

node_label_iterator bidirectional_label_map::begin() const {
    return fm_labels.begin();
}

node_label_iterator bidirectional_label_map::end() const {
    return fm_labels.end();
}

string bidirectional_label_map::get_label(id_type id) const {
    unordered_map<id_type, string>::const_iterator it = rm_ids.find(id);
    if (it != rm_ids.end()) return it->second;
    return string("");
}

id_type bidirectional_label_map::get_id(const string& label) const {
    node_label_iterator it = fm_labels.find(label);
    if (it != end()) return it->second;
    return size();
}

bool bidirectional_label_map::insert(const string& label) {
    id_type curr_size = size();
    if (get_id(label) == curr_size) {
        fm_labels.insert(make_pair(label, curr_size));
        rm_ids.insert(make_pair(curr_size, label));
        return true;
    }
    return false;
}

bool bidirectional_label_map::erase(const string& label) {
    id_type del_id = get_id(label);
    if (del_id != size()) {
        fm_labels.erase(label);
        rm_ids.erase(del_id);
        return true;
    }
    return false;
}

bool bidirectional_label_map::erase(id_type id) {
    id_type last_id = size() - 1;
    string del_label = get_label(id), rep_label = get_label(last_id);
    if (del_label != "") {
        fm_labels[rep_label] = id;
        rm_ids[id] = rep_label;
        fm_labels.erase(del_label);
        rm_ids.erase(last_id);
        return true;
    }
    return false;
}

bool bidirectional_label_map::clear() {
    if (!fm_labels.size()) return false;
    fm_labels.clear();
    rm_ids.clear();
    return true;
}

double_adjacency_map::double_adjacency_map() : st_num_edges(0), st_num_self_edges(0), wt_total_wt(0), wt_self_edge_wt(0), am_in_edges(), am_out_edges(), vw_in_degree(), vw_out_degree() {
}

id_type double_adjacency_map::num_nodes() const {
    return am_in_edges.size();
}

id_type double_adjacency_map::num_edges() const {
    return st_num_edges;
}

id_type double_adjacency_map::num_self_edges() const {
    return st_num_self_edges;
}

wt_t double_adjacency_map::total_weight() const {
    return wt_total_wt;
}

wt_t double_adjacency_map::self_edges_weight() const {
    return wt_self_edge_wt;
}

bool double_adjacency_map::is_valid_node(id_type id) const {
    return (id >= 0 && id < am_in_edges.size());
}

id_type double_adjacency_map::in_degree(id_type id) const {
    if (is_valid_node(id)) return am_in_edges[id].size();
    return 0;
}

id_type double_adjacency_map::out_degree(id_type id) const {
    if (is_valid_node(id)) return am_out_edges[id].size();
    return 0;
}

wt_t double_adjacency_map::in_degree_wt(id_type id) const {
    if (is_valid_node(id)) return vw_in_degree[id];
    return 0;
}

wt_t double_adjacency_map::out_degree_wt(id_type id) const {
    if (is_valid_node(id)) return vw_out_degree[id];
    return 0;
}

adjacent_edges_iterator double_adjacency_map::in_edges_begin(id_type id) const {
    if (is_valid_node(id)) return am_in_edges[id].begin();
    return am_in_edges[am_in_edges.size() - 1].end();
}

adjacent_edges_iterator double_adjacency_map::in_edges_end(id_type id) const {
    if (is_valid_node(id)) return am_in_edges[id].end();
    return am_in_edges[am_in_edges.size() - 1].end();
}

adjacent_edges_iterator double_adjacency_map::out_edges_begin(id_type id) const {
    if (is_valid_node(id)) return am_out_edges[id].begin();
    return am_out_edges[am_out_edges.size() - 1].end();
}

adjacent_edges_iterator double_adjacency_map::out_edges_end(id_type id) const {
    if (is_valid_node(id)) return am_out_edges[id].end();
    return am_out_edges[am_out_edges.size() - 1].end();
}

wt_t double_adjacency_map::edge_weight(id_type from_id, id_type to_id) const {
    if (is_valid_node(from_id) && is_valid_node(to_id)) {
        adjacent_edges_iterator aeit = am_out_edges[from_id].find(to_id);
        if (aeit != am_out_edges[from_id].end()) return aeit->second;
    }
    return 0;
}

id_type double_adjacency_map::insert_node() {
    am_in_edges.push_back(adjacent_edge_sequence());
    am_out_edges.push_back(adjacent_edge_sequence());
    vw_in_degree.push_back(0);
    vw_out_degree.push_back(0);
    return am_in_edges.size();
}

bool double_adjacency_map::insert_edge(id_type from_id, id_type to_id, wt_t weight) {
    if (!weight) return false;
    if (edge_weight(from_id, to_id)) return false;
    am_in_edges[to_id].insert(make_pair(from_id, weight));
    am_out_edges[from_id].insert(make_pair(to_id, weight));
    st_num_edges++;
    wt_total_wt += weight;
    vw_in_degree[to_id] += weight;
    vw_out_degree[from_id] += weight;
    if (from_id == to_id) {
        st_num_self_edges++;
        wt_self_edge_wt += weight;
    }
    return true;
}

bool double_adjacency_map::delete_edge(id_type from_id, id_type to_id) {
    wt_t weight = edge_weight(from_id, to_id);
    if (!weight) return false;
    am_in_edges[to_id].erase(from_id);
    am_out_edges[from_id].erase(to_id);
    st_num_edges--;
    wt_total_wt -= weight;
    vw_in_degree[to_id] -= weight;
    vw_out_degree[from_id] -= weight;
    if (from_id == to_id) {
        st_num_self_edges--;
        wt_self_edge_wt -= weight;
    }
    return true;
}

wt_t double_adjacency_map::set_edge_wt(id_type from_id, id_type to_id, wt_t weight) {
    wt_t old_weight = edge_weight(from_id, to_id);
    if (old_weight && weight) {
        wt_total_wt += (weight - old_weight);
        vw_in_degree[to_id] += (weight - old_weight);
        vw_out_degree[from_id] += (weight - old_weight);
        if (from_id == to_id) wt_self_edge_wt += (weight - old_weight);
        am_in_edges[to_id][from_id] += (weight - old_weight);
        am_out_edges[from_id][to_id] += (weight - old_weight);
    } else {
        if (!old_weight && weight) insert_edge(from_id, to_id, weight);
        else if (old_weight && !weight) delete_edge(from_id, to_id);
    }
    return old_weight;
}

bool double_adjacency_map::delete_node(id_type id) {
    if (!is_valid_node(id)) return false;
    id_type last_id = am_in_edges.size() - 1;
    deque< pair<id_type, id_type> > edges_to_delete;
    deque< pair< pair<id_type, id_type>, wt_type> >edges_to_add;
    if (id != last_id) {
        for (adjacent_edges_iterator aeit = am_in_edges[id].begin(); aeit != am_in_edges[id].end(); aeit++)
            edges_to_delete.push_back(make_pair(aeit->first, id));
        for (adjacent_edges_iterator aeit = am_out_edges[id].begin(); aeit != am_out_edges[id].end(); aeit++)
            edges_to_delete.push_back(make_pair(id, aeit->first));
        for (adjacent_edges_iterator aeit = am_in_edges[last_id].begin(); aeit != am_in_edges[last_id].end(); aeit++) {
            if (aeit->first != id)edges_to_add.push_back(make_pair(make_pair(aeit->first, id), aeit->second));
            edges_to_delete.push_back(make_pair(aeit->first, last_id));
        }
        for (adjacent_edges_iterator aeit = am_out_edges[last_id].begin(); aeit != am_out_edges[last_id].end(); aeit++) {
            if (aeit->first != id)edges_to_add.push_back(make_pair(make_pair(id, aeit->first), aeit->second));
            edges_to_delete.push_back(make_pair(last_id, aeit->first));
        }
    } else {
        for (adjacent_edges_iterator aeit = am_in_edges[last_id].begin(); aeit != am_in_edges[last_id].end(); aeit++)
            edges_to_delete.push_back(make_pair(aeit->first, last_id));
        for (adjacent_edges_iterator aeit = am_out_edges[last_id].begin(); aeit != am_out_edges[last_id].end(); aeit++)
            edges_to_delete.push_back(make_pair(last_id, aeit->first));
    }
    for (id_type i = 0; i < edges_to_delete.size(); i++) delete_edge(edges_to_delete[i].first, edges_to_delete[i].second);
    if (id != last_id)for (id_type i = 0; i < edges_to_add.size(); i++) insert_edge(edges_to_add[i].first.first, edges_to_add[i].first.second, edges_to_add[i].second);
    am_in_edges.pop_back();
    am_out_edges.pop_back();
    vw_in_degree.pop_back();
    vw_out_degree.pop_back();
    return true;
}

bool double_adjacency_map::delete_all_edges() {

    id_type num_nodes = am_in_edges.size();
    if (!num_nodes || !st_num_edges) return false;
    st_num_edges = 0;
    st_num_self_edges = 0;
    wt_total_wt = 0;
    wt_self_edge_wt = 0;
    am_in_edges.assign(num_nodes, adjacent_edge_sequence());
    am_out_edges.assign(num_nodes, adjacent_edge_sequence());
    vw_in_degree.assign(num_nodes, wt_t());
    vw_out_degree.assign(num_nodes, wt_t());
    return true;
}

bool double_adjacency_map::clear() {
    delete_all_edges();
    id_type num_nodes = am_in_edges.size();
    if (!num_nodes) return false;
    am_in_edges.clear();
    am_out_edges.clear();
    vw_in_degree.clear();
    vw_out_degree.clear();
    return true;
}

bool binary_heap::compare(const pair<id_type, wt_type>& left, const pair<id_type, wt_type>& right) const {
    return (b_max) ? (left.second > right.second) : (left.second < right.second);
}

id_type binary_heap::size() const {
    return vpiw_heap.size();
}

bool binary_heap::empty() const {
    return !vpiw_heap.size();
}

pair<id_type, wt_type> binary_heap::top() const {
    return vpiw_heap[0];
}

void binary_heap::heapify_up(id_type pos) {
    while (pos > 0 && compare(vpiw_heap[pos], vpiw_heap[(pos - 1) / 2])) {
        swap(umis_pos[vpiw_heap[pos].first], umis_pos[vpiw_heap[(pos - 1) / 2].first]);
        //cout << "Heapify_Up("<< pos <<") Swapping (" << vpiw_heap[pos].first << "," << vpiw_heap[pos].second << ") and (" << vpiw_heap[(pos-1)/2].first << "," << vpiw_heap[(pos-1)/2].second << ")" <<endl;
        swap(vpiw_heap[pos], vpiw_heap[(pos - 1) / 2]);
        pos = (pos - 1) / 2;
    }
}

void binary_heap::update_key(const pair<id_type, wt_type>& piw_in) {
    id_type pos = umis_pos[piw_in.first];
    double old_val = vpiw_heap[pos].second;
    vpiw_heap[pos].second = piw_in.second;
    bool logic = b_max ? (old_val < piw_in.second) : (old_val > piw_in.second);
    if (logic)heapify_up(pos);
    else heapify_down(pos);
}

void binary_heap::insert(const pair<id_type, wt_type>& piw_in) {
    vpiw_heap.push_back(piw_in);
    id_type pos = vpiw_heap.size() - 1;
    umis_pos.insert(make_pair(piw_in.first, pos));
    update_key(piw_in);
}

void binary_heap::heapify_down(id_type root) {
    if (root < vpiw_heap.size()) {
        id_type left = 2 * root + 1, right = 2 * root + 2, top = root;
        if (left < vpiw_heap.size() && compare(vpiw_heap[left], vpiw_heap[top])) top = left;
        if (right < vpiw_heap.size() && compare(vpiw_heap[right], vpiw_heap[top])) top = right;
        if (top != root) {
            swap(umis_pos[vpiw_heap[top].first], umis_pos[vpiw_heap[root].first]);
            swap(vpiw_heap[top], vpiw_heap[root]);
            heapify_down(top);
        }
    }
}

void binary_heap::pop() {
    umis_pos.erase(vpiw_heap[0].first);
    umis_pos[vpiw_heap[vpiw_heap.size() - 1].first] = 0;
    vpiw_heap[0] = vpiw_heap[vpiw_heap.size() - 1];
    vpiw_heap.pop_back();
    heapify_down(0);
}

binary_heap::binary_heap(const vector<wt_type>& v, bool max) : vpiw_heap(), umis_pos(), b_max(max) {
    for (id_type i = 0; i < v.size(); i++) {
        vpiw_heap.push_back(make_pair(i, v[i]));
        umis_pos.insert(make_pair(i, i));
    }
    for (id_type i = (vpiw_heap.size() + 1) / 2; i > 0; i--) heapify_down(i);
    heapify_down(0);

}

binary_heap::binary_heap(bool max) : vpiw_heap(), umis_pos(), b_max(max) {
}

disjoint_set::disjoint_set() : st_num_sets(0), inm_elems() {
}

id_type disjoint_set::size() const {
    return inm_elems.size();
}

id_type disjoint_set::num_sets() const {
    return st_num_sets;
}

ds_iterator disjoint_set::begin() const {
    return inm_elems.begin();
}

ds_iterator disjoint_set::end() const {
    return inm_elems.end();
}

pair<ds_iterator, bool> disjoint_set::make_set(id_type x) {
    return inm_elems.insert(make_pair(x, make_pair(x, 0)));
}

pair<ds_iterator, bool> disjoint_set::find(id_type x) {
    ds_iterator dit = inm_elems.find(x);
    if (dit == inm_elems.end()) return make_pair(dit, false);
    if (inm_elems[x].first != x)
        inm_elems[x].first = find(inm_elems[x].first).first->second.first;
    return make_pair(inm_elems.find(dit->second.first), true);
}

pair<ds_iterator, bool> disjoint_set::join(id_type x, id_type y) {
    pair<ds_iterator, bool> ret_x = find(x);
    pair<ds_iterator, bool> ret_y = find(y);
    if (!ret_x.second || !ret_y.second) return make_pair(inm_elems.end(), false);
    if (ret_x.first == ret_y.first) return ret_x;
    if (ret_x.first->second.second < ret_y.first->second.second) {
        inm_elems[ret_x.first->first].first = ret_y.first->first;
        return ret_y;
    } else {
        inm_elems[ret_y.first->first].second++;
        inm_elems[ret_y.first->first].first = ret_x.first->first;
        return ret_x;
    }
}

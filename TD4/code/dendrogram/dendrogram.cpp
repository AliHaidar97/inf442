#include "dendrogram.hpp" // This is the header for the class implemented here

#include "cloud.hpp"      // Variable c of type cloud& is used in the constructor and several other methods
#include "edge.hpp"       // Used in several methods, e.g. dendrogram::merge(edge *e)
#include "graph.hpp"      // Variable g of type graph* is used in the constructor and several other methods

#include <iostream> // This provides the input/output functionality
#include <cassert>  // This provides the assert() method

#include<map>
using namespace std;

dendrogram::dendrogram(graph& _g) {
    g = &_g;

    int n = g->get_num_nodes();

    parent = new int[n];
    rank = new int[n];

    left = new int[n];
    down = new int[n];

    height = new double[n];

    clusters = new int[n];

    for (int i = 0; i < n; i++) {
        parent[i] = -1;
        rank[i] = 0;

        left[i] = -1;
        down[i] = -1;

        height[i] = -1;

        clusters[i] = -1;
    }

    cut_height = -1;
    total_clusters = 0;
    nonsingle_clusters = 0;
}

dendrogram::~dendrogram() {
    delete g;

    delete[] parent;
    delete[] left;
    delete[] down;
    delete[] rank;
    delete[] height;
    delete[] clusters;

    if (sign_heights != nullptr)
        delete[] sign_heights;
}

int dendrogram::find(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    int result = i;
    while (parent[result] != -1) {
        result = parent[result];
    }
  
    return result;
}

void dendrogram::merge(edge *e) {
    // TODO: Exercise 3.2
    // Plan:
    // 1. find the representatives
    // 2. choose the highest
    // 3. adjust parent, left, and down
    // 4. update ranks
    // 5. update heights
    int u = e->get_p1();
    int v = e->get_p2();
    
    u = find(u);
    v = find(v);

    if (u == v) {
        return;
    }
    
    if (rank[u] < rank[v]) {
        swap(u, v);
    }
    
    if (rank[u] > rank[v]) {
        parent[v] = u;
        left[v] = down[u];
        down[u] = v;
        height[v] = e->get_length() / 2;
    }

    else {
        if (u < v) {
            swap(u, v);
        }
        rank[u]++;
        parent[v] = u;
        left[v] = down[u];
        down[u] = v;
        height[v] = e->get_length() / 2;
    }
   
}

void dendrogram::build() {
    g->start_iteration();

    while (1) {
        edge* e = g->get_next();
        if (e == NULL)break;
        merge(e);
    }

}

double dendrogram::get_dendro_height() {
    return height[down[find(0)]];
}

void dendrogram::set_clusters(int i, double h) {
    assert(i >= -1 && i < g->get_num_nodes());

    if (i == -1) {
        return;
    }

    if (parent[i] == -1) {
        clusters[i] = i;
        total_clusters++;
    }
    else {
        if (height[i] <= h) {
            clusters[i] = clusters[parent[i]];
        }
        else {
            clusters[i] = i;
            total_clusters++;
        }
    }

    set_clusters(get_down(i), h);
    set_clusters(get_left(i), h);

}

void dendrogram::set_clusters(double h) {
    if (cut_height != h) {
        cut_height = h;
        set_clusters(find(0), h);
    }
}

int dendrogram::_count_ns_clusters() {
    int count = 0;
    map<int, int>freq;
    for (int i = 0; i < g->get_num_nodes(); i++) {
        freq[clusters[i]]++;
        if (freq[clusters[i]] == 2) {
            count++;
        }
    }

    return count;
}

int dendrogram::count_ns_clusters() {
    if (nonsingle_clusters == 0)
        nonsingle_clusters = _count_ns_clusters();
    return nonsingle_clusters;
}

void dendrogram::clear_clusters() {
    for (int i = 0; i < g->get_num_nodes(); i++) {
        clusters[i] = -1;
    }
    nonsingle_clusters = 0;
    total_clusters = 0;
    cut_height = -1;
}

double dendrogram::get_cluster_height(int cluster) {
    assert(0 <= cluster && cluster < g->get_num_nodes());


    double best = 0;
   
    for (int i = 0; i < g->get_num_nodes(); i++) {
        if (i == cluster) {  continue; }
        if (clusters[i] == cluster) {
            best = max(best, max(0.0, height[i]));
        }
    }

    return best; 
}

/******** Significant heights ********/

int count_non_zero(double *unfiltered, int countu) {
    int countf = 0;
    for (int i = 0; i < countu; i++) {
        if (unfiltered[i] > 0)
            countf++;
    }

    return countf;
}

double *filter_double_array(double *unfiltered, int countu, int countf) {
    double *filtered = new double[countf];
    int f = 0;
    for (int i = 0; i < countu; i++) {
        if (unfiltered[i] > 0)
            filtered[f++] = unfiltered[i];
    }
    return filtered;
}

void dendrogram::find_heights(double eps) {
    double h = get_dendro_height();
    int slots = 1 / eps + 1;
    double *buckets = new double[slots];
    for (int i = 0; i < slots; i++)
        buckets[i] = 0;

    for (int i = 0; i < g->get_num_nodes(); i++) {
        if (height[i] != -1 && height[i] > buckets[(int)(height[i] / eps / h)])
            buckets[(int)(height[i] / eps / h)] = height[i];
    }

    count_sign_heights = count_non_zero(buckets, slots);

    if (sign_heights != nullptr)
        delete[] sign_heights;

    sign_heights = filter_double_array(buckets, slots, count_sign_heights);
    delete[] buckets;
}

/***** GETTERS *****/

long dendrogram::get_n() {
    return g->get_num_nodes();
}

int dendrogram::get_parent(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return parent[i];
}

int dendrogram::get_down(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return down[i];
}

int dendrogram::get_left(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return left[i];
}

int dendrogram::get_rank(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return rank[i];
}

int dendrogram::get_cluster(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return clusters[i];
}

std::string& dendrogram::get_name(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return g->get_name(i);
}

double dendrogram::get_height(int i) {
    assert(0 <= i && i < g->get_num_nodes());
    return height[i];
}

int dendrogram::get_count_sign_heights() {
    return count_sign_heights;
}

double dendrogram::get_sign_height(int i) {
    assert(0 <= i && i < count_sign_heights);
    return sign_heights[i];
}

int dendrogram::get_total_clusters() {
    return total_clusters;
}

double dendrogram::get_cut_height() {
    return cut_height;
}

/****** Methods for testing ******/

void dendrogram::print_node(int i) {
    cout << "Node " << i
         << "(parent =  " << parent[i]
         << ", down = " << down[i]
         << ", left = " << left[i]
         << ", rank = " << rank[i]
         << ", height = " << height[i]
         << ", cluster = " << clusters[i]
         << ")";
}

void dendrogram::print_dendrogram() {
    cout << "node\tparent\trank\tleft\tdown\theight\tcluster" << endl;
    for (int i = 0; i < get_n(); i++)
        cout << i << "\t"
             << parent[i] << "\t"
             << rank[i] << "\t"
             << left[i] << "\t"
             << down[i] << "\t"
             << height[i] << "\t"
             << clusters[i] << "\t"
             << endl;
}

void dendrogram::print_clusters() {
    int print_after[g->get_num_nodes()];
    for (size_t i = 0; i < g->get_num_nodes(); i++) {
        print_after[i] = -1;
    }
    // Go through all nodes
    for (size_t i = 0; i < g->get_num_nodes(); i++) {
        if (clusters[i] != i) {
            print_after[i] = print_after[clusters[i]];
            print_after[clusters[i]] = i;
        }
    }

    // Print the clusters
    for (int i = 0; i < g->get_num_nodes(); i++) {
        // For each encountered cluster representative
        if (clusters[i] == i) {
            // Print the corresponding cluster
            cout << "Cluster \"" << clusters[i]
                 << "\" (node: " << i << ";"
                 << " height: " << get_cluster_height(i)
                 << ")" << endl;

            int next = print_after[i];
            // Print all points in the cluster
            while (next != -1) {
                std::cout << get_name(next) << std::endl;
                next = print_after[next];
            }
        }
    }
}

void dendrogram::iterate_sign_heights() {
    static int i = 0;
    double cut = get_sign_height(i);
    cout << "Setting clusters at height " << cut
         << "..." << endl;
    clear_clusters();
    set_clusters(cut);
    cout << "\t" << count_ns_clusters() << " non-singleton cluster"
         << (count_ns_clusters() > 1 ? "s" : "")
         << " below " << cut << endl;
    i = (i + 1) % get_count_sign_heights();
}

void dendrogram::set_parent(int *_parent) {
    for (int i = 0; i < g->get_num_nodes(); i++)
        parent[i] = _parent[i];
}

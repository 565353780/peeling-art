#include <math.h>

#include <Mesh/surface_manager.h>

#include <Windows.h>
#include <GL/GL.h>
#include <GL/glu.h>

#define pi 3.14159265358979323846

#define tm "-------- test_msg : "
#define dou " , "

using namespace std;

struct easy_point
{
    double x;
    double y;
    double z;
};

struct graph_point
{
    double x;
    double y;
    double z;

    int deg = 1;
};

class EdgeToGraph
{
public:
    EdgeToGraph(Mesh& mesh_in);

    ~EdgeToGraph();

    bool is_the_same_point(easy_point point_1, easy_point point_2);
    bool is_the_same_point(easy_point point_1, graph_point point_2);
    bool is_the_same_point(graph_point point_1, easy_point point_2);
    bool is_the_same_point(graph_point point_1, graph_point point_2);

    double get_dist(easy_point point_1, easy_point point_2);
    double get_dist(easy_point point_1, graph_point point_2);
    double get_dist(graph_point point_1, easy_point point_2);
    double get_dist(graph_point point_1, graph_point point_2);

    double get_dist_2(easy_point point_1, easy_point point_2);
    double get_dist_2(easy_point point_1, graph_point point_2);
    double get_dist_2(graph_point point_1, easy_point point_2);
    double get_dist_2(graph_point point_1, graph_point point_2);

    double get_dist_of_line(vector<graph_point> points_of_line, int start_idx, int end_idx);

    double get_cross_dist(easy_point point_11, easy_point point_12, easy_point point_21, easy_point point_22);

    void get_halfedge_handle_on_bound();

    Mesh::HalfedgeHandle get_last_halfedge_handle_on_bound(Mesh::HalfedgeHandle halfedge_handle);

    Mesh::HalfedgeHandle get_next_halfedge_handle_on_bound(Mesh::HalfedgeHandle halfedge_handle);

    Mesh::HalfedgeHandle get_nearest_halfedge_handle(Mesh::HalfedgeHandle halfedge_handle);

    Mesh::HalfedgeHandle get_nearest_drawn_halfedge_handle(Mesh::HalfedgeHandle halfedge_handle);

    void get_halfedge_tag_to_draw_on_bound();

    void remove_discontinuous_lines_with_simple_dist(int simple_dist);

    void repair_discontinuous_with_simple_dist(int simple_dist);

    void connect_lines_and_get_new_lines();

    void merge_new_lines();

    void get_graph_points_and_degree();

    void remove_large_angle_on_edges(double theta);

    void get_graph();

    void get_connected_nodes(int node);

    void connect_graph();

    vector<int> get_single_farthest_path(int node, int node_old);

    void get_farthest_path();

    void show_points();

    void show_idx_of_points();

    void set_color_mode();

    void draw_selected_edge_with_color(vector<int> nodes, float red, float green, float blue);

    void draw_child_edges(vector<int> nodes, int color_idx);

    void draw_color_on_total_graph();

    void draw_edges_on_bound(float red, float green, float blue);

    void draw_all_edges(float red, float green, float blue);

    void draw_all_edges_on_bound(float red, float green, float blue);

    easy_point get_surface_point(double theta, double phi);

    void get_surface_points(double warp_angle, double weft_angle, int density);

    void draw_surface_map(float red, float green, float blue);

    void show_idx_of_surface_map();

public:
    Mesh mesh;

    vector<Mesh::HalfedgeHandle> halfedge_handle_on_bound;

    vector<vector<easy_point>> points_of_new_lines;

    vector<vector<graph_point>> points_of_lines_on_graph;

    vector<vector<easy_point>> surface_points;

    vector<int> tag_arr;

    vector<vector<float>> color_mode;

    vector<vector<double>> graph;

    vector<int> node_deg;

    vector<vector<vector<int>>> edge_tag;

    vector<int> farthest_path;

    vector<int> go_through_tag;

    vector<int> nodes_connected_tag;

    vector<vector<int>> connected_nodes;

    double dist_of_edges = 0;

    int surface_lines_density;

    GLuint TextFont;
};

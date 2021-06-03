#include "EdgeToGraph.h"

extern OpenMesh::EPropHandleT<int> edge_type;

EdgeToGraph::EdgeToGraph(Mesh& mesh_in)
{
    mesh = mesh_in;

    get_halfedge_handle_on_bound();

    get_halfedge_tag_to_draw_on_bound();

    remove_discontinuous_lines_with_simple_dist(5);

    repair_discontinuous_with_simple_dist(3);

    connect_lines_and_get_new_lines();

    merge_new_lines();

    get_graph_points_and_degree();

    remove_large_angle_on_edges(150);

    get_graph();

    get_farthest_path();

    connect_graph();

    set_color_mode();

    get_surface_points(10, 10, 10);

    int MAX_CHAR = 5120;

    TextFont = glGenLists(MAX_CHAR);

    wglUseFontBitmaps(wglGetCurrentDC(), 0, MAX_CHAR, TextFont);
}

EdgeToGraph::~EdgeToGraph()
{

}

bool EdgeToGraph::is_the_same_point(easy_point point_1, easy_point point_2)
{
    if(point_1.x == point_2.x && point_1.y == point_2.y && point_1.z == point_2.z)
    {
        return true;
    }

    return false;
}
bool EdgeToGraph::is_the_same_point(easy_point point_1, graph_point point_2)
{
    if(point_1.x == point_2.x && point_1.y == point_2.y && point_1.z == point_2.z)
    {
        return true;
    }

    return false;
}
bool EdgeToGraph::is_the_same_point(graph_point point_1, easy_point point_2)
{
    if(point_1.x == point_2.x && point_1.y == point_2.y && point_1.z == point_2.z)
    {
        return true;
    }

    return false;
}
bool EdgeToGraph::is_the_same_point(graph_point point_1, graph_point point_2)
{
    if(point_1.x == point_2.x && point_1.y == point_2.y && point_1.z == point_2.z)
    {
        return true;
    }

    return false;
}


double EdgeToGraph::get_dist(easy_point point_1, easy_point point_2)
{
    return sqrt(get_dist_2(point_1, point_2));
}
double EdgeToGraph::get_dist(easy_point point_1, graph_point point_2)
{
    return sqrt(get_dist_2(point_1, point_2));
}
double EdgeToGraph::get_dist(graph_point point_1, easy_point point_2)
{
    return sqrt(get_dist_2(point_1, point_2));
}
double EdgeToGraph::get_dist(graph_point point_1, graph_point point_2)
{
    return sqrt(get_dist_2(point_1, point_2));
}

double EdgeToGraph::get_dist_2(easy_point point_1, easy_point point_2)
{
    return (point_1.x - point_2.x)*(point_1.x - point_2.x) + (point_1.y - point_2.y)*(point_1.y - point_2.y) + (point_1.z - point_2.z)*(point_1.z - point_2.z);
}
double EdgeToGraph::get_dist_2(easy_point point_1, graph_point point_2)
{
    return (point_1.x - point_2.x)*(point_1.x - point_2.x) + (point_1.y - point_2.y)*(point_1.y - point_2.y) + (point_1.z - point_2.z)*(point_1.z - point_2.z);
}
double EdgeToGraph::get_dist_2(graph_point point_1, easy_point point_2)
{
    return (point_1.x - point_2.x)*(point_1.x - point_2.x) + (point_1.y - point_2.y)*(point_1.y - point_2.y) + (point_1.z - point_2.z)*(point_1.z - point_2.z);
}
double EdgeToGraph::get_dist_2(graph_point point_1, graph_point point_2)
{
    return (point_1.x - point_2.x)*(point_1.x - point_2.x) + (point_1.y - point_2.y)*(point_1.y - point_2.y) + (point_1.z - point_2.z)*(point_1.z - point_2.z);
}

double EdgeToGraph::get_dist_of_line(vector<graph_point> points_of_line, int start_idx, int end_idx)
{
    if(start_idx >= end_idx)
    {
        return 0;
    }

    double dist_of_line = 0;

    for(int i = start_idx; i < end_idx; ++i)
    {
        dist_of_line += get_dist(points_of_line[i], points_of_line[i + 1]);
    }

    return dist_of_line;
}

double EdgeToGraph::get_cross_dist(easy_point point_11, easy_point point_12, easy_point point_21, easy_point point_22)
{
    return get_dist(point_11, point_21) + get_dist(point_11, point_22) + get_dist(point_12, point_21) + get_dist(point_12, point_22);
}

void EdgeToGraph::get_halfedge_handle_on_bound()
{
    if(halfedge_handle_on_bound.size())
    {
        halfedge_handle_on_bound.clear();
    }
    for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
    {
        if (mesh.property(edge_type, e_it) == 0)
        {
            halfedge_handle_on_bound.push_back(mesh.halfedge_handle(e_it.handle(), 0));
        }
    }
}

Mesh::HalfedgeHandle EdgeToGraph::get_last_halfedge_handle_on_bound(Mesh::HalfedgeHandle halfedge_handle)
{
    Mesh::HalfedgeHandle old_handle = halfedge_handle;
    Mesh::HalfedgeHandle temp_handle = mesh.next_halfedge_handle(mesh.next_halfedge_handle(old_handle));

    while(mesh.next_halfedge_handle(temp_handle) != old_handle)
    {
        temp_handle = mesh.next_halfedge_handle(temp_handle);
    }

    while(mesh.face_handle(mesh.opposite_halfedge_handle(temp_handle)).idx() != -1)
    {
        old_handle = mesh.opposite_halfedge_handle(temp_handle);
        temp_handle = mesh.next_halfedge_handle(mesh.next_halfedge_handle(old_handle));

        while(mesh.next_halfedge_handle(temp_handle) != old_handle)
        {
            temp_handle = mesh.next_halfedge_handle(temp_handle);
        }
    }

    return temp_handle;
}

Mesh::HalfedgeHandle EdgeToGraph::get_next_halfedge_handle_on_bound(Mesh::HalfedgeHandle halfedge_handle)
{
    Mesh::HalfedgeHandle temp_handle = mesh.next_halfedge_handle(halfedge_handle);

    while(mesh.face_handle(mesh.opposite_halfedge_handle(temp_handle)).idx() != -1)
    {
        temp_handle = mesh.opposite_halfedge_handle(temp_handle);

        temp_handle = mesh.next_halfedge_handle(temp_handle);
    }

    return temp_handle;
}

Mesh::HalfedgeHandle EdgeToGraph::get_nearest_halfedge_handle(Mesh::HalfedgeHandle halfedge_handle)
{
    easy_point p1, p2, p3, p4;

    Mesh::HalfedgeHandle min_dist_handle;

    double min_dist = 10000;
    double current_dist = 0;

    p1.x = mesh.point(mesh.from_vertex_handle(halfedge_handle))[0];
    p1.y = mesh.point(mesh.from_vertex_handle(halfedge_handle))[1];
    p1.z = mesh.point(mesh.from_vertex_handle(halfedge_handle))[2];

    p2.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
    p2.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
    p2.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

    for(int i = 0; i < halfedge_handle_on_bound.size(); ++i)
    {
        if(halfedge_handle_on_bound[i] == halfedge_handle
                || halfedge_handle_on_bound[i] == get_last_halfedge_handle_on_bound(halfedge_handle)
                || halfedge_handle_on_bound[i] == get_next_halfedge_handle_on_bound(halfedge_handle)
                || halfedge_handle_on_bound[i] == get_last_halfedge_handle_on_bound(get_last_halfedge_handle_on_bound(halfedge_handle))
                || halfedge_handle_on_bound[i] == get_next_halfedge_handle_on_bound(get_next_halfedge_handle_on_bound(halfedge_handle)))
        {
            continue;
        }

        p3.x = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[0];
        p3.y = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[1];
        p3.z = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[2];

        p4.x = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[0];
        p4.y = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[1];
        p4.z = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[2];

        current_dist = get_cross_dist(p1, p2, p3, p4);

        if(current_dist < min_dist)
        {
            min_dist_handle = halfedge_handle_on_bound[i];
            min_dist = current_dist;
        }
    }

    dist_of_edges = min_dist;

    return min_dist_handle;
}

Mesh::HalfedgeHandle EdgeToGraph::get_nearest_drawn_halfedge_handle(Mesh::HalfedgeHandle halfedge_handle)
{
    easy_point p1, p2, p3, p4;

    Mesh::HalfedgeHandle min_dist_drawn_handle;

    double min_dist = 10000;
    double current_dist = 0;

    p1.x = mesh.point(mesh.from_vertex_handle(halfedge_handle))[0];
    p1.y = mesh.point(mesh.from_vertex_handle(halfedge_handle))[1];
    p1.z = mesh.point(mesh.from_vertex_handle(halfedge_handle))[2];

    p2.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
    p2.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
    p2.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

    for(int i = 0; i < halfedge_handle_on_bound.size(); ++i)
    {
        if(halfedge_handle_on_bound[i] == halfedge_handle
                || halfedge_handle_on_bound[i] == get_last_halfedge_handle_on_bound(halfedge_handle)
                || halfedge_handle_on_bound[i] == get_next_halfedge_handle_on_bound(halfedge_handle)
                || halfedge_handle_on_bound[i] == get_last_halfedge_handle_on_bound(get_last_halfedge_handle_on_bound(halfedge_handle))
                || halfedge_handle_on_bound[i] == get_next_halfedge_handle_on_bound(get_next_halfedge_handle_on_bound(halfedge_handle))
                || tag_arr[halfedge_handle_on_bound[i].idx()] != 1)
        {
            continue;
        }

        p3.x = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[0];
        p3.y = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[1];
        p3.z = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[2];

        p4.x = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[0];
        p4.y = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[1];
        p4.z = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[2];

        current_dist = get_cross_dist(p1, p2, p3, p4);

        if(current_dist < min_dist)
        {
            min_dist_drawn_handle = halfedge_handle_on_bound[i];
            min_dist = current_dist;
        }
    }

    return min_dist_drawn_handle;
}

void EdgeToGraph::get_halfedge_tag_to_draw_on_bound()
{
    double dist_cut = 1.5;

    Mesh::HalfedgeHandle halfedge_handle = halfedge_handle_on_bound[0];

    tag_arr.resize(mesh.n_halfedges());

    vector<int> tag_arr_positive;
    vector<int> tag_arr_reverse;
    tag_arr_positive.resize(mesh.n_halfedges());
    tag_arr_reverse.resize(mesh.n_halfedges());
    int positive_edge_num = 0;
    int reverse_edge_num = 0;

    //Positive
    if(tag_arr_positive[halfedge_handle.idx()] == 0)
    {
        if(tag_arr_positive[get_nearest_halfedge_handle(halfedge_handle).idx()] != 1)
        {

            tag_arr_positive[halfedge_handle.idx()] = 1;
            ++positive_edge_num;
        }
        else
        {
            if(dist_of_edges > dist_cut)
            {
                tag_arr_positive[halfedge_handle.idx()] = 1;
                ++positive_edge_num;
            }
        }
    }

    halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

    while(halfedge_handle != halfedge_handle_on_bound[0])
    {
        if(tag_arr_positive[halfedge_handle.idx()] == 0)
        {
            if(tag_arr_positive[get_nearest_halfedge_handle(halfedge_handle).idx()] != 1)
            {

                tag_arr_positive[halfedge_handle.idx()] = 1;
                ++positive_edge_num;
            }
            else
            {
                if(dist_of_edges > dist_cut)
                {
                    tag_arr_positive[halfedge_handle.idx()] = 1;
                    ++positive_edge_num;
                }
            }
        }

        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
    }

    //Reverse
    if(tag_arr_reverse[halfedge_handle.idx()] == 0)
    {
        if(tag_arr_reverse[get_nearest_halfedge_handle(halfedge_handle).idx()] != 1)
        {

            tag_arr_reverse[halfedge_handle.idx()] = 1;
            ++reverse_edge_num;
        }
        else
        {
            if(dist_of_edges > dist_cut)
            {
                tag_arr_reverse[halfedge_handle.idx()] = 1;
                ++reverse_edge_num;
            }
        }
    }

    halfedge_handle = get_last_halfedge_handle_on_bound(halfedge_handle);

    while(halfedge_handle != halfedge_handle_on_bound[0])
    {
        if(tag_arr_reverse[halfedge_handle.idx()] == 0)
        {
            if(tag_arr_reverse[get_nearest_halfedge_handle(halfedge_handle).idx()] != 1)
            {

                tag_arr_reverse[halfedge_handle.idx()] = 1;
                ++reverse_edge_num;
            }
            else
            {
                if(dist_of_edges > dist_cut)
                {
                    tag_arr_reverse[halfedge_handle.idx()] = 1;
                    ++reverse_edge_num;
                }
            }
        }

        halfedge_handle = get_last_halfedge_handle_on_bound(halfedge_handle);
    }

    if(positive_edge_num > reverse_edge_num)
    {
        for(int i = 0; i < tag_arr.size(); ++i)
        {
            tag_arr[i] = tag_arr_positive[i];
        }
    }
    else
    {
        for(int i = 0; i < tag_arr.size(); ++i)
        {
            tag_arr[i] = tag_arr_reverse[i];
        }
    }
}

void EdgeToGraph::remove_discontinuous_lines_with_simple_dist(int simple_dist)
{
    if(simple_dist < 1)
    {
        return;
    }

    Mesh::HalfedgeHandle halfedge_handle = halfedge_handle_on_bound[0];
    Mesh::HalfedgeHandle find_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

    if(tag_arr[halfedge_handle.idx()] != 1 && tag_arr[find_handle.idx()] == 1)
    {
        halfedge_handle = find_handle;

        for(int i = 0; i < simple_dist; ++i)
        {
            find_handle = get_next_halfedge_handle_on_bound(find_handle);

            if(tag_arr[find_handle.idx()] != 1)
            {
                break;
            }
        }

        if(tag_arr[find_handle.idx()] != 1)
        {
            while(halfedge_handle != find_handle && halfedge_handle != halfedge_handle_on_bound[0])
            {
                tag_arr[halfedge_handle.idx()] = 0;

                halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
            }
        }
        else
        {
            while(halfedge_handle != find_handle && halfedge_handle != halfedge_handle_on_bound[0])
            {
                halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
            }
        }
    }
    else
    {
        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
    }

    if(halfedge_handle == halfedge_handle_on_bound[0])
    {
        return;
    }

    while(true)
    {
        find_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

        if(tag_arr[halfedge_handle.idx()] != 1 && tag_arr[find_handle.idx()] == 1)
        {
            halfedge_handle = find_handle;

            for(int i = 0; i < simple_dist; ++i)
            {
                find_handle = get_next_halfedge_handle_on_bound(find_handle);

                if(tag_arr[find_handle.idx()] != 1)
                {
                    break;
                }
            }

            if(tag_arr[find_handle.idx()] != 1)
            {
                while(halfedge_handle != find_handle && halfedge_handle != halfedge_handle_on_bound[0])
                {
                    tag_arr[halfedge_handle.idx()] = 0;

                    halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
                }
            }
            else
            {
                while(halfedge_handle != find_handle && halfedge_handle != halfedge_handle_on_bound[0])
                {
                    halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
                }
            }
        }
        else
        {
            halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
        }

        if(halfedge_handle == halfedge_handle_on_bound[0])
        {
            break;
        }
    }
}

void EdgeToGraph::repair_discontinuous_with_simple_dist(int simple_dist)
{
    if(simple_dist < 1)
    {
        return;
    }

    bool have_last_edge = false;
    bool have_next_edge = false;

    bool edge_tag_changed = true;

    Mesh::HalfedgeHandle halfedge_handle;
    Mesh::HalfedgeHandle find_handle;

    while(edge_tag_changed)
    {
        edge_tag_changed = false;

        halfedge_handle = halfedge_handle_on_bound[0];
        find_handle = halfedge_handle;

        if(tag_arr[halfedge_handle.idx()] != 1)
        {
            for(int i = 0; i < simple_dist; ++i)
            {
                find_handle = get_last_halfedge_handle_on_bound(find_handle);
                if(tag_arr[find_handle.idx()] == 1)
                {
                    have_last_edge = true;

                    break;
                }
            }

            find_handle = halfedge_handle;

            for(int i = 0; i < simple_dist; ++i)
            {
                find_handle = get_next_halfedge_handle_on_bound(find_handle);
                if(tag_arr[find_handle.idx()] == 1)
                {
                    have_next_edge = true;

                    break;
                }
            }

            if(have_last_edge && have_next_edge)
            {
                tag_arr[halfedge_handle.idx()] = 1;

                edge_tag_changed = true;

                continue;
            }
        }

        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
        find_handle = halfedge_handle;
        have_last_edge = false;
        have_next_edge = false;

        while(halfedge_handle != halfedge_handle_on_bound[0])
        {
            if(tag_arr[halfedge_handle.idx()] != 1)
            {
                for(int i = 0; i < simple_dist; ++i)
                {
                    find_handle = get_last_halfedge_handle_on_bound(find_handle);
                    if(tag_arr[find_handle.idx()] == 1)
                    {
                        have_last_edge = true;

                        break;
                    }
                }

                find_handle = halfedge_handle;

                for(int i = 0; i < simple_dist; ++i)
                {
                    find_handle = get_next_halfedge_handle_on_bound(find_handle);
                    if(tag_arr[find_handle.idx()] == 1)
                    {
                        have_next_edge = true;

                        break;
                    }
                }

                if(have_last_edge && have_next_edge)
                {
                    tag_arr[halfedge_handle.idx()] = 1;

                    edge_tag_changed = true;

                    continue;
                }
            }

            halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
            find_handle = halfedge_handle;
            have_last_edge = false;
            have_next_edge = false;
        }
    }
}

void EdgeToGraph::connect_lines_and_get_new_lines()
{
    double max_connect_dist = 0.05;

    Mesh::HalfedgeHandle halfedge_handle = halfedge_handle_on_bound[0];

    Mesh::HalfedgeHandle find_handle;

    Mesh::HalfedgeHandle nearest_drawn_handle;

    Mesh::HalfedgeHandle start_handle;

    bool self_cross = false;

    vector<easy_point> points_of_line;

    easy_point point, p1, p2;

    if(tag_arr[halfedge_handle.idx()] == 1)
    {
        while(tag_arr[halfedge_handle.idx()] == 1)
        {
            halfedge_handle = get_last_halfedge_handle_on_bound(halfedge_handle);
        }

        start_handle = halfedge_handle;
    }
    else
    {
        while(tag_arr[halfedge_handle.idx()] != 1)
        {
            halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
        }

        halfedge_handle = get_last_halfedge_handle_on_bound(halfedge_handle);

        start_handle = halfedge_handle;
    }

    find_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

    nearest_drawn_handle = get_nearest_drawn_halfedge_handle(find_handle);
    //nearest_drawn_handle = get_nearest_drawn_halfedge_handle(halfedge_handle);

    for(int i = 0; i < 50; ++i)
    {
        find_handle = get_next_halfedge_handle_on_bound(find_handle);

        if(find_handle == nearest_drawn_handle)
        {
            self_cross = true;

            break;
        }
    }

    point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
    point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
    point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

    if(!self_cross)
    {
        p1.x = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[0];
        p1.y = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[1];
        p1.z = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[2];

        p2.x = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[0];
        p2.y = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[1];
        p2.z = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[2];

        if(get_dist_2(point, p1) < get_dist_2(point, p2) && get_dist(point, p1) < max_connect_dist)
        {
            points_of_line.push_back(p1);
            points_of_line.push_back(point);
        }
        else if(get_dist(point, p2) < max_connect_dist)
        {
            points_of_line.push_back(p2);
            points_of_line.push_back(point);
        }
    }
    else
    {
        points_of_line.push_back(point);
    }

    self_cross = false;

    halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

    point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
    point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
    point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

    points_of_line.push_back(point);

    while(tag_arr[halfedge_handle.idx()] + tag_arr[get_next_halfedge_handle_on_bound(halfedge_handle).idx()] != 1)
    {
        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

        point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
        point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
        point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

        points_of_line.push_back(point);
    }

    find_handle = halfedge_handle;

    nearest_drawn_handle = get_nearest_drawn_halfedge_handle(find_handle);
    //nearest_drawn_handle = get_next_halfedge_handle_on_bound(get_nearest_drawn_halfedge_handle(halfedge_handle));

    for(int i = 0; i < 50; ++i)
    {
        find_handle = get_last_halfedge_handle_on_bound(find_handle);

        if(find_handle == nearest_drawn_handle)
        {
            self_cross = true;

            break;
        }
    }

    point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
    point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
    point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

    if(!self_cross)
    {
        p1.x = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[0];
        p1.y = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[1];
        p1.z = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[2];

        p2.x = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[0];
        p2.y = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[1];
        p2.z = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[2];

        if(get_dist_2(point, p1) < get_dist_2(point, p2) && get_dist(point, p1) < max_connect_dist)
        {
            points_of_line.push_back(p1);
        }
        else if(get_dist(point, p2) < max_connect_dist)
        {
            points_of_line.push_back(p2);
        }
    }

    self_cross = false;

    halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

    while(tag_arr[halfedge_handle.idx()] + tag_arr[get_next_halfedge_handle_on_bound(halfedge_handle).idx()] != 1)
    {
        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
    }

    points_of_new_lines.push_back(points_of_line);

    points_of_line.clear();

    while(halfedge_handle != start_handle)
    {
        find_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

        nearest_drawn_handle = get_nearest_drawn_halfedge_handle(find_handle);
        //nearest_drawn_handle = get_nearest_drawn_halfedge_handle(halfedge_handle);

        for(int i = 0; i < 10; ++i)
        {
            find_handle = get_next_halfedge_handle_on_bound(find_handle);

            if(find_handle == nearest_drawn_handle)
            {
                self_cross = true;

                break;
            }
        }

        point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
        point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
        point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

        if(!self_cross)
        {
            p1.x = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[0];
            p1.y = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[1];
            p1.z = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[2];

            p2.x = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[0];
            p2.y = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[1];
            p2.z = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[2];

            if(get_dist_2(point, p1) < get_dist_2(point, p2) && get_dist(point, p1) < max_connect_dist)
            {
                points_of_line.push_back(p1);
                points_of_line.push_back(point);
            }
            else if(get_dist(point, p2) < max_connect_dist)
            {
                points_of_line.push_back(p2);
                points_of_line.push_back(point);
            }
        }
        else
        {
            points_of_line.push_back(point);
        }

        self_cross = false;

        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

        point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
        point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
        point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

        points_of_line.push_back(point);

        while(tag_arr[halfedge_handle.idx()] + tag_arr[get_next_halfedge_handle_on_bound(halfedge_handle).idx()] != 1)
        {
            halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

            point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
            point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
            point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

            points_of_line.push_back(point);
        }

        find_handle = halfedge_handle;

        nearest_drawn_handle = get_nearest_drawn_halfedge_handle(find_handle);
        //nearest_drawn_handle = get_next_halfedge_handle_on_bound(get_nearest_drawn_halfedge_handle(halfedge_handle));

        for(int i = 0; i < 10; ++i)
        {
            find_handle = get_last_halfedge_handle_on_bound(find_handle);

            if(find_handle == nearest_drawn_handle)
            {
                self_cross = true;

                break;
            }
        }

        point.x = mesh.point(mesh.to_vertex_handle(halfedge_handle))[0];
        point.y = mesh.point(mesh.to_vertex_handle(halfedge_handle))[1];
        point.z = mesh.point(mesh.to_vertex_handle(halfedge_handle))[2];

        if(!self_cross)
        {
            p1.x = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[0];
            p1.y = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[1];
            p1.z = mesh.point(mesh.from_vertex_handle(nearest_drawn_handle))[2];

            p2.x = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[0];
            p2.y = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[1];
            p2.z = mesh.point(mesh.to_vertex_handle(nearest_drawn_handle))[2];

            if(get_dist_2(point, p1) < get_dist_2(point, p2) && get_dist(point, p1) < max_connect_dist)
            {
                points_of_line.push_back(p1);
            }
            else if(get_dist(point, p2) < max_connect_dist)
            {
                points_of_line.push_back(p2);
            }
        }

        self_cross = false;

        halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);

        while(tag_arr[halfedge_handle.idx()] + tag_arr[get_next_halfedge_handle_on_bound(halfedge_handle).idx()] != 1)
        {
            halfedge_handle = get_next_halfedge_handle_on_bound(halfedge_handle);
        }

        points_of_new_lines.push_back(points_of_line);

        points_of_line.clear();
    }

    bool have_repeated_edge = true;

    while(have_repeated_edge)
    {
        have_repeated_edge = false;

        for(int i = 0; i < points_of_new_lines.size() - 1; ++i)
        {
            if(have_repeated_edge)
            {
                break;
            }

            for(int j = i + 1; j < points_of_new_lines.size(); ++j)
            {
                if(is_the_same_point(points_of_new_lines[i][0], points_of_new_lines[j][1]))
                {
                    points_of_new_lines[j].erase(points_of_new_lines[j].begin());

                    have_repeated_edge = true;

                    break;
                }
                if(is_the_same_point(points_of_new_lines[i][0], points_of_new_lines[j][points_of_new_lines[j].size() - 2]))
                {
                    points_of_new_lines[j].pop_back();

                    have_repeated_edge = true;

                    break;
                }
                if(is_the_same_point(points_of_new_lines[i][points_of_new_lines[i].size() - 1], points_of_new_lines[j][1]))
                {
                    points_of_new_lines[i].pop_back();

                    have_repeated_edge = true;

                    break;
                }
                if(is_the_same_point(points_of_new_lines[i][points_of_new_lines[i].size() - 1], points_of_new_lines[j][points_of_new_lines[j].size() - 2]))
                {
                    points_of_new_lines[i].pop_back();

                    have_repeated_edge = true;

                    break;
                }
            }
        }
    }
}

void EdgeToGraph::merge_new_lines()
{
    bool lines_changed = true;

    bool cross_error = false;

    vector<easy_point> new_line;

    while(lines_changed)
    {
        lines_changed = false;

        for(int i = 0; i < points_of_new_lines.size() - 1; ++i)
        {
            if(lines_changed)
            {
                break;
            }

            for(int j = i + 1; j < points_of_new_lines.size(); ++j)
            {
                if(is_the_same_point(points_of_new_lines[i][0], points_of_new_lines[j][0]))
                {
                    for(int ii = 0; ii < points_of_new_lines.size(); ++ii)
                    {
                        if(cross_error)
                        {
                            break;
                        }

                        for(int jj = 1; jj < points_of_new_lines[ii].size() - 1; ++jj)
                        {
                            if(is_the_same_point(points_of_new_lines[ii][jj], points_of_new_lines[i][0]))
                            {
                                cross_error = true;

                                break;
                            }
                        }
                    }

                    if(cross_error)
                    {
                        cross_error = false;
                    }
                    else
                    {
                        new_line.clear();

                        for(int ii = points_of_new_lines[i].size() - 1; ii > -1; --ii)
                        {
                            new_line.push_back(points_of_new_lines[i][ii]);
                        }
                        for(int jj = 1; jj < points_of_new_lines[j].size(); ++jj)
                        {
                            new_line.push_back(points_of_new_lines[j][jj]);
                        }

                        points_of_new_lines.push_back(new_line);

                        points_of_new_lines.erase(points_of_new_lines.begin()+j);
                        points_of_new_lines.erase(points_of_new_lines.begin()+i);

                        lines_changed = true;

                        break;
                    }
                }
                else if(is_the_same_point(points_of_new_lines[i][0], points_of_new_lines[j][points_of_new_lines[j].size() - 1]))
                {
                    for(int ii = 0; ii < points_of_new_lines.size(); ++ii)
                    {
                        if(cross_error)
                        {
                            break;
                        }

                        for(int jj = 1; jj < points_of_new_lines[ii].size() - 1; ++jj)
                        {
                            if(is_the_same_point(points_of_new_lines[ii][jj], points_of_new_lines[i][0]))
                            {
                                cross_error = true;

                                break;
                            }
                        }
                    }

                    if(cross_error)
                    {
                        cross_error = false;
                    }
                    else
                    {
                        new_line.clear();

                        for(int jj = 0; jj < points_of_new_lines[j].size(); ++jj)
                        {
                            new_line.push_back(points_of_new_lines[j][jj]);
                        }
                        for(int ii = 1; ii < points_of_new_lines[i].size(); ++ii)
                        {
                            new_line.push_back(points_of_new_lines[i][ii]);
                        }

                        points_of_new_lines.push_back(new_line);

                        points_of_new_lines.erase(points_of_new_lines.begin()+j);
                        points_of_new_lines.erase(points_of_new_lines.begin()+i);

                        lines_changed = true;

                        break;
                    }
                }
                else if(is_the_same_point(points_of_new_lines[i][points_of_new_lines[i].size() - 1], points_of_new_lines[j][0]))
                {
                    for(int ii = 0; ii < points_of_new_lines.size(); ++ii)
                    {
                        if(cross_error)
                        {
                            break;
                        }

                        for(int jj = 1; jj < points_of_new_lines[ii].size() - 1; ++jj)
                        {
                            if(is_the_same_point(points_of_new_lines[ii][jj], points_of_new_lines[i][points_of_new_lines[i].size() - 1]))
                            {
                                cross_error = true;

                                break;
                            }
                        }
                    }

                    if(cross_error)
                    {
                        cross_error = false;
                    }
                    else
                    {
                        new_line.clear();

                        for(int ii = 0; ii < points_of_new_lines[i].size(); ++ii)
                        {
                            new_line.push_back(points_of_new_lines[i][ii]);
                        }
                        for(int jj = 1; jj < points_of_new_lines[j].size(); ++jj)
                        {
                            new_line.push_back(points_of_new_lines[j][jj]);
                        }

                        points_of_new_lines.push_back(new_line);

                        points_of_new_lines.erase(points_of_new_lines.begin()+j);
                        points_of_new_lines.erase(points_of_new_lines.begin()+i);

                        lines_changed = true;

                        break;
                    }
                }
                else if(is_the_same_point(points_of_new_lines[i][points_of_new_lines[i].size() - 1], points_of_new_lines[j][points_of_new_lines[j].size() - 1]))
                {
                    for(int ii = 0; ii < points_of_new_lines.size(); ++ii)
                    {
                        if(cross_error)
                        {
                            break;
                        }

                        for(int jj = 1; jj < points_of_new_lines[ii].size() - 1; ++jj)
                        {
                            if(is_the_same_point(points_of_new_lines[ii][jj], points_of_new_lines[i][points_of_new_lines[i].size() - 1]))
                            {
                                cross_error = true;

                                break;
                            }
                        }
                    }

                    if(cross_error)
                    {
                        cross_error = false;
                    }
                    else
                    {
                        new_line.clear();

                        for(int ii = 0; ii < points_of_new_lines[i].size(); ++ii)
                        {
                            new_line.push_back(points_of_new_lines[i][ii]);
                        }
                        for(int jj = points_of_new_lines[j].size() - 2; jj > -1; --jj)
                        {
                            new_line.push_back(points_of_new_lines[j][jj]);
                        }

                        points_of_new_lines.push_back(new_line);

                        points_of_new_lines.erase(points_of_new_lines.begin()+j);
                        points_of_new_lines.erase(points_of_new_lines.begin()+i);

                        lines_changed = true;

                        break;
                    }
                }
            }
        }
    }
}

void EdgeToGraph::get_graph_points_and_degree()
{
    points_of_lines_on_graph.clear();

    vector<graph_point> points_of_line;
    graph_point point;

    for(int i = 0; i < points_of_new_lines.size(); ++i)
    {
        points_of_line.clear();

        for(int j = 0; j < points_of_new_lines[i].size(); ++j)
        {
            point.x = points_of_new_lines[i][j].x;
            point.y = points_of_new_lines[i][j].y;
            point.z = points_of_new_lines[i][j].z;

            points_of_line.push_back(point);
        }

        points_of_lines_on_graph.push_back(points_of_line);
    }

    for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
    {
        for(int j = 0; j < points_of_lines_on_graph[i].size(); ++j)
        {
            for(int ii = 0; ii < points_of_new_lines.size(); ++ii)
            {
                for(int jj = 0; jj < points_of_new_lines[ii].size(); ++jj)
                {
                    if(is_the_same_point(points_of_lines_on_graph[i][j], points_of_new_lines[ii][jj]))
                    {
                        ++points_of_lines_on_graph[i][j].deg;
                    }
                }
            }
        }
    }

    node_deg.resize(2 * points_of_lines_on_graph.size());

    for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
    {
        if(points_of_lines_on_graph[i][0].deg == 2)
        {
            --points_of_lines_on_graph[i][0].deg;
        }

        if(points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].deg == 2)
        {
            --points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].deg;
        }

        node_deg[2 * i] = points_of_lines_on_graph[i][0].deg;
        node_deg[2 * i + 1] = points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].deg;
    }
}

void EdgeToGraph::remove_large_angle_on_edges(double theta)
{
    bool edges_changed = true;

    while(edges_changed)
    {
        edges_changed = false;

        for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
        {
            if(edges_changed)
            {
                break;
            }

            for(int j = 1; j < points_of_lines_on_graph[i].size() - 1; ++j)
            {
                if(points_of_lines_on_graph[i][j].deg == 2 &&
                        ((points_of_lines_on_graph[i][j].x - points_of_lines_on_graph[i][j-1].x) * (points_of_lines_on_graph[i][j+1].x - points_of_lines_on_graph[i][j].x) +
                        (points_of_lines_on_graph[i][j].y - points_of_lines_on_graph[i][j-1].y) * (points_of_lines_on_graph[i][j+1].y - points_of_lines_on_graph[i][j].y) +
                        (points_of_lines_on_graph[i][j].z - points_of_lines_on_graph[i][j-1].z) * (points_of_lines_on_graph[i][j+1].z - points_of_lines_on_graph[i][j].z)) /
                        (get_dist(points_of_lines_on_graph[i][j], points_of_lines_on_graph[i][j-1]) * get_dist(points_of_lines_on_graph[i][j], points_of_lines_on_graph[i][j+1])) <
                        cos(theta / pi))
                {
                    points_of_lines_on_graph[i].erase(points_of_lines_on_graph[i].begin()+j);

                    points_of_new_lines[i].erase(points_of_new_lines[i].begin()+j);

                    get_graph_points_and_degree();

                    edges_changed = true;

                    break;
                }
            }
        }
    }
}

void EdgeToGraph::get_graph()
{
    graph.resize(2 * points_of_lines_on_graph.size());

    for(int i = 0; i < graph.size(); ++i)
    {
        graph[i].resize(graph.size());

        for(int j = 0; j < graph.size(); ++j)
        {
            graph[i][j] = -1;
        }
    }

    edge_tag.resize(graph.size());

    for(int i = 0; i < graph.size(); ++i)
    {
        edge_tag[i].resize(graph.size());

        for(int j = 0; j < graph.size(); ++j)
        {
            edge_tag[i][j].resize(3);

            edge_tag[i][j][0] = -1;
            edge_tag[i][j][1] = -1;
            edge_tag[i][j][2] = -1;
        }
    }

    int deg_3_idx = 1;
    int deg_3_idx_old = 0;
    int connect_node_idx = -1;
    int connect_node_idx_old = -1;

    for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
    {
        graph[2 * i][2 * i + 1] = get_dist_of_line(points_of_lines_on_graph[i], 0, points_of_lines_on_graph[i].size() - 1);
        graph[2 * i + 1][2 * i] = graph[2 * i][2 * i + 1];

        edge_tag[2 * i][2 * i + 1][0] = i;
        edge_tag[2 * i][2 * i + 1][1] = 0;
        edge_tag[2 * i][2 * i + 1][2] = points_of_lines_on_graph[i].size() - 1;
        edge_tag[2 * i + 1][2 * i][0] = edge_tag[2 * i][2 * i + 1][0];
        edge_tag[2 * i + 1][2 * i][1] = edge_tag[2 * i][2 * i + 1][1];
        edge_tag[2 * i + 1][2 * i][2] = edge_tag[2 * i][2 * i + 1][2];

        deg_3_idx = 1;
        deg_3_idx_old = 0;
        connect_node_idx = -1;
        connect_node_idx_old = -1;

        while(deg_3_idx != points_of_lines_on_graph[i].size() - 1)
        {
            if(points_of_lines_on_graph[i][deg_3_idx].deg >= 3)
            {
                graph[2 * i][2 * i + 1] = -1;
                graph[2 * i + 1][2 * i] = -1;

                edge_tag[2 * i][2 * i + 1][0] = -1;
                edge_tag[2 * i][2 * i + 1][1] = -1;
                edge_tag[2 * i][2 * i + 1][2] = -1;
                edge_tag[2 * i + 1][2 * i][0] = -1;
                edge_tag[2 * i + 1][2 * i][1] = -1;
                edge_tag[2 * i + 1][2 * i][2] = -1;

                for(int j = 0; j < points_of_lines_on_graph.size(); ++j)
                {
                    if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx], points_of_lines_on_graph[j][0]))
                    {
                        connect_node_idx = 2 * j;
                    }
                    if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx], points_of_lines_on_graph[j][points_of_lines_on_graph[j].size() - 1]))
                    {
                        connect_node_idx = 2 * j + 1;
                    }

                    if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx_old], points_of_lines_on_graph[j][0]))
                    {
                        connect_node_idx_old = 2 * j;
                    }
                    if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx_old], points_of_lines_on_graph[j][points_of_lines_on_graph[j].size() - 1]))
                    {
                        connect_node_idx_old = 2 * j + 1;
                    }
                }

                graph[connect_node_idx][connect_node_idx_old] = get_dist_of_line(points_of_lines_on_graph[i], deg_3_idx_old, deg_3_idx);
                graph[connect_node_idx_old][connect_node_idx] = graph[connect_node_idx][connect_node_idx_old];

                edge_tag[connect_node_idx][connect_node_idx_old][0] = i;
                edge_tag[connect_node_idx][connect_node_idx_old][1] = deg_3_idx_old;
                edge_tag[connect_node_idx][connect_node_idx_old][2] = deg_3_idx;
                edge_tag[connect_node_idx_old][connect_node_idx][0] = edge_tag[connect_node_idx][connect_node_idx_old][0];
                edge_tag[connect_node_idx_old][connect_node_idx][1] = edge_tag[connect_node_idx][connect_node_idx_old][1];
                edge_tag[connect_node_idx_old][connect_node_idx][2] = edge_tag[connect_node_idx][connect_node_idx_old][2];

                deg_3_idx_old = deg_3_idx;
            }

            ++deg_3_idx;
        }

        for(int j = 0; j < points_of_lines_on_graph.size(); ++j)
        {
            if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx], points_of_lines_on_graph[j][0]))
            {
                connect_node_idx = 2 * j;
            }
            if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx], points_of_lines_on_graph[j][points_of_lines_on_graph[j].size() - 1]))
            {
                connect_node_idx = 2 * j + 1;
            }

            if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx_old], points_of_lines_on_graph[j][0]))
            {
                connect_node_idx_old = 2 * j;
            }
            if(is_the_same_point(points_of_lines_on_graph[i][deg_3_idx_old], points_of_lines_on_graph[j][points_of_lines_on_graph[j].size() - 1]))
            {
                connect_node_idx_old = 2 * j + 1;
            }
        }

        graph[connect_node_idx][connect_node_idx_old] = get_dist_of_line(points_of_lines_on_graph[i], deg_3_idx_old, deg_3_idx);
        graph[connect_node_idx_old][connect_node_idx] = graph[connect_node_idx][connect_node_idx_old];

        edge_tag[connect_node_idx][connect_node_idx_old][0] = i;
        edge_tag[connect_node_idx][connect_node_idx_old][1] = deg_3_idx_old;
        edge_tag[connect_node_idx][connect_node_idx_old][2] = deg_3_idx;
        edge_tag[connect_node_idx_old][connect_node_idx][0] = edge_tag[connect_node_idx][connect_node_idx_old][0];
        edge_tag[connect_node_idx_old][connect_node_idx][1] = edge_tag[connect_node_idx][connect_node_idx_old][1];
        edge_tag[connect_node_idx_old][connect_node_idx][2] = edge_tag[connect_node_idx][connect_node_idx_old][2];
    }
}

vector<int> EdgeToGraph::get_single_farthest_path(int node, int node_old)
{
    vector<int> single_farthest_path;
    vector<int> single_farthest_path_return;

    double length_of_path = 0;
    double length_of_path_return;

    bool can_go = false;

    for(int i = 0; i < graph.size(); ++i)
    {
        if(graph[node][i] != -1 && i != node_old)
        {
            can_go = true;

            single_farthest_path_return = get_single_farthest_path(i, node);

            length_of_path_return = 0;

            for(int j = 0; j < single_farthest_path_return.size() - 1; ++j)
            {
                length_of_path_return += graph[single_farthest_path_return[j]][single_farthest_path_return[j + 1]];
            }

            if(length_of_path_return > length_of_path)
            {
                length_of_path = length_of_path_return;

                single_farthest_path = single_farthest_path_return;
            }
        }
    }

    if(!can_go)
    {
        single_farthest_path.push_back(node);
    }

    if(node != node_old)
    {
        single_farthest_path.push_back(node_old);
    }

    return single_farthest_path;
}

void EdgeToGraph::get_farthest_path()
{
    vector<int> farthest_path_return;

    double length_of_path = 0;
    double length_of_path_return;

    for(int i = 0; i < edge_tag.size(); ++i)
    {
        if(node_deg[i] == 1)
        {
            farthest_path_return = get_single_farthest_path(i, i);

            length_of_path_return = 0;

            for(int j = 0; j < farthest_path_return.size() - 1; ++j)
            {
                length_of_path_return += graph[farthest_path_return[j]][farthest_path_return[j + 1]];
            }

            if(length_of_path_return > length_of_path)
            {
                length_of_path = length_of_path_return;

                farthest_path = farthest_path_return;
            }
        }
    }
}

void EdgeToGraph::get_connected_nodes(int node)
{
    nodes_connected_tag[node] = 1;

    vector<int> connected_node;

    connected_node.resize(2);

    for(int i = 0; i < graph.size(); ++i)
    {
        if(graph[node][i] != -1 && nodes_connected_tag[i] == 0)
        {            
            connected_node[0] = node;
            connected_node[1] = i;
            connected_nodes.push_back(connected_node);

            get_connected_nodes(i);
        }
    }
}

void EdgeToGraph::connect_graph()
{
    bool graph_connected = false;

    double min_dist = 1000;
    double current_dist = 0;

    int node_on_edge_1_idx_1 = -1;
    int node_on_edge_1_idx_2 = -1;
    int node_on_edge_2_idx_1 = -1;
    int node_on_edge_2_idx_2 = -1;

    int line_idx_1 = -1;
    int line_idx_2 = -1;

    int point_idx_1 = -1;
    int point_idx_2 = -1;

    vector<vector<vector<int>>> connected_nodes_list;

    while(!graph_connected)
    {
        graph_connected = true;

        min_dist = 1000;
        current_dist = 0;

        node_on_edge_1_idx_1 = -1;
        node_on_edge_1_idx_2 = -1;
        node_on_edge_2_idx_1 = -1;
        node_on_edge_2_idx_2 = -1;

        line_idx_1 = -1;
        line_idx_2 = -1;

        point_idx_1 = -1;
        point_idx_2 = -1;

        connected_nodes_list.clear();

        nodes_connected_tag.resize(graph.size());

        for(int i = 0; i < nodes_connected_tag.size(); ++i)
        {
            nodes_connected_tag[i] = 0;
        }

        connected_nodes.clear();

        get_connected_nodes(farthest_path[0]);

        int disconnected_node = -1;

        for(int i = 0; i < graph.size(); ++i)
        {
            if(nodes_connected_tag[i] == 0)
            {
                graph_connected = false;

                disconnected_node = i;

                break;
            }
        }

        if(graph_connected)
        {
            break;
        }

        if(connected_nodes.size())
        {
            connected_nodes_list.push_back(connected_nodes);
        }

        while(!graph_connected)
        {
            graph_connected = true;

            connected_nodes.clear();

            get_connected_nodes(disconnected_node);

            if(connected_nodes.size())
            {
                connected_nodes_list.push_back(connected_nodes);
            }

            for(int i = 0; i < graph.size(); ++i)
            {
                if(nodes_connected_tag[i] == 0)
                {
                    graph_connected = false;

                    disconnected_node = i;

                    break;
                }
            }
        }

        if(connected_nodes_list.size() == 1)
        {
            break;
        }

        graph_connected = false;

        for(int i = 0; i < connected_nodes_list.size() - 1; ++i)
        {
            for(int j = i + 1; j < connected_nodes_list.size(); ++j)
            {
                for(int i_1 = 0; i_1 < connected_nodes_list[i].size(); ++i_1)
                {
                    for(int j_1 = edge_tag[connected_nodes_list[i][i_1][0]][connected_nodes_list[i][i_1][1]][1];
                        j_1 <= edge_tag[connected_nodes_list[i][i_1][0]][connected_nodes_list[i][i_1][1]][2];
                        ++j_1)
                    {
                        for(int i_2 = 0; i_2 < connected_nodes_list[j].size(); ++i_2)
                        {
                            for(int j_2 = edge_tag[connected_nodes_list[j][i_2][0]][connected_nodes_list[j][i_2][1]][1];
                                j_2 <= edge_tag[connected_nodes_list[j][i_2][0]][connected_nodes_list[j][i_2][1]][2];
                                ++j_2)
                            {
                                if(points_of_lines_on_graph[edge_tag[connected_nodes_list[i][i_1][0]][connected_nodes_list[i][i_1][1]][0]][j_1].deg == 3 ||
                                        points_of_lines_on_graph[edge_tag[connected_nodes_list[j][i_2][0]][connected_nodes_list[j][i_2][1]][0]][j_2].deg == 3)
                                {
                                    continue;
                                }

                                current_dist = get_dist(points_of_lines_on_graph[edge_tag[connected_nodes_list[i][i_1][0]][connected_nodes_list[i][i_1][1]][0]][j_1],
                                        points_of_lines_on_graph[edge_tag[connected_nodes_list[j][i_2][0]][connected_nodes_list[j][i_2][1]][0]][j_2]);

                                if(current_dist < min_dist)
                                {
                                    min_dist = current_dist;

                                    node_on_edge_1_idx_1 = connected_nodes_list[i][i_1][0];
                                    node_on_edge_1_idx_2 = connected_nodes_list[i][i_1][1];
                                    node_on_edge_2_idx_1 = connected_nodes_list[j][i_2][0];
                                    node_on_edge_2_idx_2 = connected_nodes_list[j][i_2][1];

                                    line_idx_1 = edge_tag[node_on_edge_1_idx_1][node_on_edge_1_idx_2][0];
                                    line_idx_2 = edge_tag[node_on_edge_2_idx_1][node_on_edge_2_idx_2][0];

                                    point_idx_1 = j_1;
                                    point_idx_2 = j_2;
                                }
                            }
                        }
                    }
                }
            }
        }

        node_deg.push_back(1);
        node_deg.push_back(1);

        vector<graph_point> points_of_line;

        points_of_line.push_back(points_of_lines_on_graph[line_idx_1][point_idx_1]);
        points_of_line.push_back(points_of_lines_on_graph[line_idx_2][point_idx_2]);

        points_of_lines_on_graph.push_back(points_of_line);

        points_of_line.clear();

        points_of_new_lines.clear();

        points_of_new_lines.resize(points_of_lines_on_graph.size());

        for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
        {
            points_of_new_lines[i].resize(points_of_lines_on_graph[i].size());

            for(int j = 0; j < points_of_lines_on_graph[i].size(); ++j)
            {
                points_of_new_lines[i][j].x = points_of_lines_on_graph[i][j].x;
                points_of_new_lines[i][j].y = points_of_lines_on_graph[i][j].y;
                points_of_new_lines[i][j].z = points_of_lines_on_graph[i][j].z;
            }
        }

        merge_new_lines();

        get_graph_points_and_degree();

        get_graph();

        get_farthest_path();
    }
}

void EdgeToGraph::show_points()
{
    double delta = 1.005;

    glPointSize(6);

    glColor3f(1.0f, 0.0f, 0.0f);

    glBegin(GL_POINTS);

    for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
    {
        glVertex3f(points_of_lines_on_graph[i][0].x * delta,
                points_of_lines_on_graph[i][0].y * delta,
                points_of_lines_on_graph[i][0].z * delta);

        glVertex3f(points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].x * delta,
                points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].y * delta,
                points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].z * delta);
    }

    glEnd();

/*    glColor3f(0.0f, 0.0f, 0.0f);

    glBegin(GL_POINTS);

    int i = 27;
    int j = 12;

    glVertex3f(points_of_lines_on_graph[i][j].x * delta,
            points_of_lines_on_graph[i][j].y * delta,
            points_of_lines_on_graph[i][j].z * delta);

    i = 12;
    j = 21;

    glVertex3f(points_of_lines_on_graph[i][j].x * delta,
            points_of_lines_on_graph[i][j].y * delta,
            points_of_lines_on_graph[i][j].z * delta);

    glEnd();*/

    glFlush();
}

void EdgeToGraph::show_idx_of_points()
{
    double delta = 1.05;

    for(int i = 0; i < points_of_lines_on_graph.size(); ++i)
    {
        char s1[10];
        _itoa(2 * i, s1, 10);

        glRasterPos3f(points_of_lines_on_graph[i][0].x * delta,
                points_of_lines_on_graph[i][0].y * delta,
                points_of_lines_on_graph[i][0].z * delta);

        glListBase(TextFont);
        glCallLists(strlen(s1), GL_UNSIGNED_BYTE, s1);

        char s2[10];
        _itoa(2 * i + 1, s2, 10);

        glRasterPos3f(points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].x * delta,
                points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].y * delta,
                points_of_lines_on_graph[i][points_of_lines_on_graph[i].size() - 1].z * delta);

        glListBase(TextFont);
        glCallLists(strlen(s2), GL_UNSIGNED_BYTE, s2);
    }
}

void EdgeToGraph::set_color_mode()
{
    color_mode.resize(4);

    for(int i = 0; i < color_mode.size(); ++i)
    {
        color_mode[i].resize(3);
    }

    color_mode[0][0] = 255;
    color_mode[0][1] = 0;
    color_mode[0][2] = 0;

    color_mode[1][0] = 0;
    color_mode[1][1] = 255;
    color_mode[1][2] = 0;

    color_mode[2][0] = 255;
    color_mode[2][1] = 255;
    color_mode[2][2] = 0;

    color_mode[3][0] = 0;
    color_mode[3][1] = 255;
    color_mode[3][2] = 255;
}

void EdgeToGraph::draw_selected_edge_with_color(vector<int> nodes, float red, float green, float blue)
{
    double delta = 1.002;

    glLineWidth(1.5);

    glEnable(GL_LINE_SMOOTH);

    double* p = new double[3];

    glColor3f(red/255.0, green/255.0, blue/255.0);

    glBegin(GL_LINES);

    for(int i = 0; i < nodes.size() - 1; ++i)
    {
        for(int j = edge_tag[nodes[i]][nodes[i + 1]][1]; j < edge_tag[nodes[i]][nodes[i + 1]][2]; ++j)
        {
            p[0] = points_of_lines_on_graph[edge_tag[nodes[i]][nodes[i + 1]][0]][j].x * delta;
            p[1] = points_of_lines_on_graph[edge_tag[nodes[i]][nodes[i + 1]][0]][j].y * delta;
            p[2] = points_of_lines_on_graph[edge_tag[nodes[i]][nodes[i + 1]][0]][j].z * delta;
            glVertex3dv(&p[0]);
            p[0] = points_of_lines_on_graph[edge_tag[nodes[i]][nodes[i + 1]][0]][j + 1].x * delta;
            p[1] = points_of_lines_on_graph[edge_tag[nodes[i]][nodes[i + 1]][0]][j + 1].y * delta;
            p[2] = points_of_lines_on_graph[edge_tag[nodes[i]][nodes[i + 1]][0]][j + 1].z * delta;
            glVertex3dv(&p[0]);
        }
    }

    glEnd();

    delete[] p;
}

void EdgeToGraph::draw_child_edges(vector<int> nodes, int color_idx)
{
    for(int i = 0; i < nodes.size(); ++i)
    {
        go_through_tag[nodes[i]] = 1;
    }

    draw_selected_edge_with_color(nodes, color_mode[color_idx][0], color_mode[color_idx][1], color_mode[color_idx][2]);

    if(nodes.size() < 3)
    {
        return;
    }

    for(int i = 1; i < nodes.size() - 1; ++i)
    {
        if(node_deg[nodes[i]] == 3)
        {
            for(int j = 0; j < graph.size(); ++j)
            {
                if(graph[nodes[i]][j] != -1 && go_through_tag[j] == 0)
                {
                    vector<int> child_farthest_path;

                    child_farthest_path = get_single_farthest_path(j, nodes[i]);

                    draw_child_edges(child_farthest_path, color_idx + 1);
                }
            }
        }
    }
}

void EdgeToGraph::draw_color_on_total_graph()
{
    go_through_tag.resize(node_deg.size());

    draw_child_edges(farthest_path, 0);

    for(int i = 0; i < go_through_tag.size(); ++i)
    {
        go_through_tag[i] = 0;
    }
}

void EdgeToGraph::draw_edges_on_bound(float red, float green, float blue)
{
    double delta = 1.001;

    glLineWidth(1.5);

    glEnable(GL_LINE_SMOOTH);

    double* p = new double[3];

    /*
    for(int i = 0; i < halfedge_handle_on_bound.size(); ++i)
    {
        if(tag_arr[halfedge_handle_on_bound[i].idx()] == 1)
        {
            p[0] = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[0] * delta;
            p[1] = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[1] * delta;
            p[2] = mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[2] * delta;
            glVertex3dv(&p[0]);
            p[0] = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[0] * delta;
            p[1] = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[1] * delta;
            p[2] = mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[2] * delta;
            glVertex3dv(&p[0]);
        }
    }
    */

    for(int i = 0; i < points_of_new_lines.size(); ++i)
    {
        glColor3f(red/255.0, green/255.0, blue/255.0);

        glBegin(GL_LINES);

        for(int j = 0; j < points_of_new_lines[i].size() - 1; ++j)
        {
            p[0] = points_of_new_lines[i][j].x * delta;
            p[1] = points_of_new_lines[i][j].y * delta;
            p[2] = points_of_new_lines[i][j].z * delta;
            glVertex3dv(&p[0]);
            p[0] = points_of_new_lines[i][j + 1].x * delta;
            p[1] = points_of_new_lines[i][j + 1].y * delta;
            p[2] = points_of_new_lines[i][j + 1].z * delta;
            glVertex3dv(&p[0]);
        }

        glEnd();
    }

    delete[] p;
}

void EdgeToGraph::draw_all_edges(float red, float green, float blue)
{
    for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
    {
        if (mesh.property(edge_type, e_it) <= 0)
        {
            glColor3f(red/255.0, green/255.0, blue/255.0);
            Mesh::EdgeHandle e_handle;
            Mesh::HalfedgeHandle he_handle;

            e_handle = e_it.handle();
            he_handle = mesh.halfedge_handle(e_handle, 0);
            glBegin(GL_LINES);
            glVertex3dv(&mesh.point(mesh.from_vertex_handle(he_handle))[0]);
            glVertex3dv(&mesh.point(mesh.to_vertex_handle(he_handle))[0]);
            glEnd();
        }
    }
}

void EdgeToGraph::draw_all_edges_on_bound(float red, float green, float blue)
{
    glColor3f(red/255.0, green/255.0, blue/255.0);

    for(int i = 0; i < halfedge_handle_on_bound.size(); ++i)
    {
        glBegin(GL_LINES);
        glVertex3dv(&mesh.point(mesh.from_vertex_handle(halfedge_handle_on_bound[i]))[0]);
        glVertex3dv(&mesh.point(mesh.to_vertex_handle(halfedge_handle_on_bound[i]))[0]);
        glEnd();
    }
}

easy_point EdgeToGraph::get_surface_point(double theta, double phi)
{
    easy_point point;

    double lambda = 0.76;
    double a = 1.618;

    double r = (1 - lambda)* pow(cos(pi * phi / 180), a) + lambda;
    point.x = r*cos(pi * theta / 180)*cos(pi * phi / 180);
    point.y = r*sin(pi * theta / 180)*cos(pi * phi / 180);
    point.z = r*sin(pi * phi / 180);

    return point;
}

void EdgeToGraph::get_surface_points(double warp_angle, double weft_angle, int density)
{
    surface_lines_density = density;

    int warp_num = 360 / warp_angle;
    int weft_num = 180 / weft_angle;

    warp_angle = 360.0 / warp_num;
    weft_angle = 180.0 / weft_num;

    double delta_warp = warp_angle / surface_lines_density;
    double delta_weft = weft_angle / surface_lines_density;

    surface_points.resize(warp_num * surface_lines_density + 1);

    for(int i = 0; i < surface_points.size(); ++i)
    {
        surface_points[i].resize(weft_num * surface_lines_density + 1);

        for(int j = 0; j < surface_points[i].size(); ++j)
        {
            surface_points[i][j] = get_surface_point(i*delta_warp, j*delta_weft - 90);
        }
    }
}

void EdgeToGraph::draw_surface_map(float red, float green, float blue)
{
    double delta = 1.001;

    glLineWidth(1.5);

    glEnable(GL_LINE_SMOOTH);

    double* p = new double[3];

    glColor3f(red/255.0, green/255.0, blue/255.0);

    glBegin(GL_LINES);

    for(int i = 0; i < surface_points.size(); ++i)
    {
        for(int j = 0; j < surface_points[i].size(); ++j)
        {
            if(j%surface_lines_density == 0 && i < surface_points.size() - 1)
            {
                p[0] = surface_points[i][j].x * delta;
                p[1] = surface_points[i][j].y * delta;
                p[2] = surface_points[i][j].z * delta;
                glVertex3dv(&p[0]);

                p[0] = surface_points[i+1][j].x * delta;
                p[1] = surface_points[i+1][j].y * delta;
                p[2] = surface_points[i+1][j].z * delta;
                glVertex3dv(&p[0]);
            }

            if(i%surface_lines_density == 0 && j < surface_points[i].size() - 1)
            {
                p[0] = surface_points[i][j].x * delta;
                p[1] = surface_points[i][j].y * delta;
                p[2] = surface_points[i][j].z * delta;
                glVertex3dv(&p[0]);

                p[0] = surface_points[i][j+1].x * delta;
                p[1] = surface_points[i][j+1].y * delta;
                p[2] = surface_points[i][j+1].z * delta;
                glVertex3dv(&p[0]);
            }
        }
    }

    glEnd();

    delete[] p;
}

void EdgeToGraph::show_idx_of_surface_map()
{
    double delta = 1.05;

    double warp_angle = 10.0;
    double weft_angle = 10.0;

    for(int i = 0; i < 360.0 / warp_angle; ++i)
    {
        double warp_pose = -2;

        easy_point point;

        point = get_surface_point(1.0*i*warp_angle, 1.0*warp_pose);

        char s1[10];
        _itoa(i * warp_angle, s1, 10);

        glRasterPos3f(point.x * delta, point.y * delta, point.z * delta);

        glListBase(TextFont);
        glCallLists(strlen(s1), GL_UNSIGNED_BYTE, s1);
    }

    for(int i = 1; i < 180.0 / weft_angle; ++i)
    {
        double weft_pose = 0;

        easy_point point;

        point = get_surface_point(1.0*weft_pose, 1.0*i*weft_angle - 90);

        char s2[10];
        _itoa(i * weft_angle - 90, s2, 10);

        glRasterPos3f(point.x * delta, point.y * delta, point.z * delta);

        glListBase(TextFont);
        glCallLists(strlen(s2), GL_UNSIGNED_BYTE, s2);
    }
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Object.h>
#include "arr_print.h"
#include "vertical_decomposition.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Quotient<CGAL::MP_Float> Number_type;
typedef CGAL::Cartesian<Number_type> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::X_monotone_curve_2 Segment_2;

typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Vertex Vertex_2;
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> Trap_pl;

void print_neighboring_vertices(Arrangement_2::Vertex_const_handle v) {
    if (v->is_isolated()) {
        std::cout << "The vertex (" << v->point() << ") is isolated" << std::endl;
        return;
    }
    Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
    first = curr = v->incident_halfedges();
    std::cout << "The neighbors of the vertex (" << v->point() << ") are:";
    do {
        // Note that the current halfedge is directed from u to v:
        Arrangement_2::Vertex_const_handle u = curr->source();
        std::cout << " (" << u->point() << ")";
    } while (++curr != first);
    std::cout << std::endl;
}

void print_ccb(Arrangement_2::Ccb_halfedge_const_circulator circ) {
    Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    std::cout << "(" << curr->source()->point() << ")";
    do {
        Arrangement_2::Halfedge_const_handle he = curr;
        std::cout << " [" << he->curve() << "] "
                  << "(" << he->target()->point() << ")";
    } while (++curr != circ);
    std::cout << std::endl;
}

void print_face(Arrangement_2::Face_const_handle f) {
    // Print the outer boundary.
    if (f->is_unbounded())
        std::cout << "Unbounded face. " << std::endl;
    else {
        std::cout << "Outer boundary: ";
        print_ccb(f->outer_ccb());
    }
    // Print the boundary of each of the holes.
    Arrangement_2::Hole_const_iterator hi;
    int index = 1;
    for (hi = f->holes_begin(); hi != f->holes_end(); ++hi, ++index) {
        std::cout << " Hole #" << index << ": ";
        print_ccb(*hi);
    }
    // Print the isolated vertices.
    Arrangement_2::Isolated_vertex_const_iterator iv;
    for (iv = f->isolated_vertices_begin(), index = 1;
         iv != f->isolated_vertices_end(); ++iv, ++index) {
        std::cout << " Isolated vertex #" << index << ": "
                  << "(" << iv->point() << ")" << std::endl;
    }
}

void add_polygon(Arrangement_2 *arr, const std::vector<Point_2> &points) {
    std::cout << "Adding new polygon with points { ";
    Segment_2 segments[points.size()];
    for (int i = 0; i < points.size(); i++) {
        segments[i] = Segment_2(points.at(i), points.at((i + 1) % points.size()));
        std::cout << "(" << points.at(i) << ") ";
    }
    std::cout << "}" << std::endl;
    CGAL::insert(*arr, &segments[0], &segments[points.size()]);
}

int main() {
    // Defines a polygon
//    Point points[] = {Point(0, 0), Point(5.1, 0), Point(1, 1), Point(0.5, 6)};
//    Polygon_2 pgn(points, points + 4);
    // check if the polygon is simple.
//    std::cout << "The polygon is " <<
//         (pgn.is_simple() ? "" : "not ") << "simple." << std::endl;
//    // check if the polygon is convex
//    std::cout << "The polygon is " <<
//         (pgn.is_convex() ? "" : "not ") << "convex." << std::endl;

    // An arrangement
    Arrangement_2 arr;
    // Add small square
    add_polygon(&arr, std::vector<Point_2>(
            {Point_2(0, 0), Point_2(5, 0), Point_2(5, 5), Point_2(0, 5)}));
    // Add larger trapezoid above square
    add_polygon(&arr, std::vector<Point_2>(
            {Point_2(-10, 10), Point_2(10, 10), Point_2(5, 15), Point_2(-1, 15)}));
    // Add a bounding box
    add_polygon(&arr, std::vector<Point_2>(
            {Point_2(-20, 20), Point_2(20, 20), Point_2(20, -20), Point_2(-20, -20)}));

    std::cout << std::endl;

//    print_neighboring_vertices(arr.vertices_begin());

    // Print the faces of the arrangement
//    for (Arrangement_2::Face_const_iterator fit = arr.faces_begin(); fit != arr.faces_end(); fit++) {
//        std::cout << "Theres an face! Bounded: " << !fit->is_unbounded() << std::endl;
//        print_face(fit);
//    }

    if (false) {
        typedef Arrangement_2::Vertex_const_handle Vertex_const_handle;
        typedef std::pair<Vertex_const_handle, std::pair<CGAL::Object, CGAL::Object>> Vert_decomp_entry;
        typedef std::list<Vert_decomp_entry> Vert_decomp_list;
        Vert_decomp_list vd_list;
        CGAL::decompose(arr, std::back_inserter(vd_list));

        typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
        typedef Arrangement_2::Face_const_handle Face_const_handle;
        // Print the results.
        Vert_decomp_list::const_iterator vd_iter;
        std::pair<CGAL::Object, CGAL::Object> curr;
        Vertex_const_handle vh;
        Halfedge_const_handle hh;
        Face_const_handle fh;

        for (vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
            curr = vd_iter->second;
            std::cout << "Vertex (" << vd_iter->first->point() << ") : ";

            std::cout << " feature below: ";
            if (CGAL::assign(vh, curr.first))
                std::cout << '(' << vh->point() << ')';
            else if (CGAL::assign(fh, curr.first))
                if (!fh->is_fictitious())
                    std::cout << '[' << fh->outer_ccb().ptr() << ']';
                else
                    std::cout << "NONE";
            else
                std::cout << "EMPTY";

            std::cout << "   feature above: ";
            if (CGAL::assign(vh, curr.second))
                std::cout << '(' << vh->point() << ')' << std::endl;
            else if (CGAL::assign(hh, curr.second))
                if (!fh->is_fictitious())
                    std::cout << '[' << fh->outer_ccb().ptr() << ']' << std::endl;
                else
                    std::cout << "NONE" << std::endl;
            else
                std::cout << "EMPTY" << std::endl;
        }
    }

    std::cout << "Added all segments to arrangement." << std::endl;
    print_arrangement_size(arr);

    // Add vertical edges that induce the vertical decomposition.
    Traits_2 traits;
    Kernel *kernel = &traits;
    vertical_decomposition(arr, *kernel);
    std::cout << std::endl << "Added vertical bullet paths." << std::endl;
    print_arrangement_size(arr);

    return 0;
}
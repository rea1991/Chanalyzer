//---------------------------------------------------------
/**    @file	Auxiliary.h
 *     @brief	Auxiliary.h is the header for CLASS
 *               Auxiliary.cpp								*/
//---------------------------------------------------------

#ifndef Auxiliary_h
#define Auxiliary_h

#include "globals.h"
#include <map>

#ifdef DBGMEM_CRT
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#endif

#ifdef ENABLE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <fstream>
#include <unistd.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#endif


namespace PMP = CGAL::Polygon_mesh_processing;

class PointCell; //Line added to avoid "unspecified identifier" error
class FacetCell; //Line added to avoid "unspecified identifier" error


// Regular triangulation and alpha-shapes:
typedef CGAL::Exact_predicates_inexact_constructions_kernel                 K;
typedef CGAL::Triangulation_vertex_base_with_info_3<PointCell*,K>           Vertex_with_info;
typedef CGAL::Regular_triangulation_vertex_base_3<K,Vertex_with_info>       Vb;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<K,Vb>                         FAS_Vb;
typedef CGAL::Triangulation_cell_base_with_info_3< FacetCell*[4],K>         Cell_with_info;
typedef CGAL::Regular_triangulation_cell_base_3<K,Cell_with_info>           RT_Cell_with_info;
typedef CGAL::Fixed_alpha_shape_cell_base_3<K,RT_Cell_with_info>            FAS_Cell_with_info;
typedef CGAL::Triangulation_data_structure_3<FAS_Vb,FAS_Cell_with_info>     Tds;
typedef Tds::Cell_circulator                                                Cell_circulator;
typedef CGAL::Regular_triangulation_3<K,Tds>                                Rt;
typedef CGAL::Fixed_alpha_shape_3<Rt>                                       Fixed_alpha_shape_3;
typedef Fixed_alpha_shape_3::Cell_handle                                    Alpha_Cell_handle;
typedef Fixed_alpha_shape_3::Vertex_handle                                  Alpha_Vertex_handle;
typedef Fixed_alpha_shape_3::Facet                                          Alpha_Facet;
typedef Fixed_alpha_shape_3::Edge                                           Alpha_Edge;
typedef K::Point_3                                                          Point3;
typedef K::Weighted_point_3                                                 Weighted_point;
typedef Rt::Vertex_iterator                                                 Vertex_iterator;
typedef Rt::Cell_iterator                                                   Cell_iterator;
typedef Rt::Finite_vertices_iterator                                        Finite_Vertex_Iterator;
typedef Rt::Finite_cells_iterator                                           Finite_Cells_Iterator;
typedef Rt::Finite_edges_iterator                                           Finite_Edges_Iterator;
typedef Rt::Finite_facets_iterator                                          Finite_Facets_Iterator;
typedef Rt::Vertex_handle                                                   Vertex_handle;
typedef Rt::Cell_handle                                                     Cell_handle;
typedef Rt::Facet                                                           Facet;
typedef Rt::Edge                                                            Edge;
typedef Rt::Cell                                                            Cell;


bool orthogonal_sphere(Weighted_point& v, Weighted_point& v1, Weighted_point& v2, Weighted_point& v3, Weighted_point& v4);
bool flow_relation(Point3& vb,  Weighted_point& oc, std::vector<Weighted_point>& common_vertices);
double vol_tet(Weighted_point& v1, Weighted_point& v2, Weighted_point& v3, Weighted_point& v4);
double surf_tet(Weighted_point& v1, Weighted_point& v2, Weighted_point& v3, Weighted_point& v4);
double sur_tets(std::map<int, Cell_handle>& int2tetrahedron, std::vector<int>& flagged_tetrahedra, std::map<std::vector<int>, std::vector<Weighted_point>>& tetrahedra2common_vertices, int flag);
double sur_tets(std::vector <int> indices_setOfTetrahedra, std::map<int, Cell_handle>& int2tetrahedron, std::map<std::vector<int>, std::vector<Weighted_point>>& tetrahedra2common_vertices);
double vol_tets(std::map<int, Cell_handle>& int2tetrahedron, std::vector<int>& flagged_tetrahedra, int flag);
bool write_triangulation2off(Rt& reg_triang, std::string root_file_name, const int d);
bool write_alphaFilt(Fixed_alpha_shape_3& alpha, const int flag = 0);
bool write_flagged_tetrahedra2off(Fixed_alpha_shape_3& alpha_shape, std::vector<int> flagged_tetrahedra, std::string file_name, int K);

#endif

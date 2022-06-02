//---------------------------------------------------------
/**    @file	Chanalyzer.h
 *     @brief	Chanalyzer.h is the header for CLASS
 *               Chanalyzer.cpp								*/
//---------------------------------------------------------

#ifndef Chanalyzer_h
#define Chanalyzer_h

#include "Auxiliary.h"
#include "ConnollySurface.h"
#include "globals.h"
#include "Graph.h"

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
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#endif


namespace PMP = CGAL::Polygon_mesh_processing;



class Chanalyzer
{
    
#ifdef ENABLE_CGAL
private:
    // Mesh processing -- skeletonization
    typedef CGAL::Simple_cartesian<double>                        Kernel;
    typedef Kernel::Point_3                                       Point;
    typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
    typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
    typedef Skeletonization::Skeleton                             Skeleton;
    typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
    typedef Skeleton::edge_descriptor                             Skeleton_edge;
    typedef Polyhedron::Vertex_handle                             Vertex_handle_P;
    typedef Polyhedron::Halfedge_handle                           Halfedge_handle_P;
    typedef Polyhedron::Facet_handle                              Facet_handle_P;
    
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
    
    //the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
    #define RESET   "\033[0m"
    #define BLACK   "\033[30m"      /* Black */
    #define RED     "\033[31m"      /* Red */
    #define GREEN   "\033[32m"      /* Green */
    #define YELLOW  "\033[33m"      /* Yellow */
    #define BLUE    "\033[34m"      /* Blue */
    #define MAGENTA "\033[35m"      /* Magenta */
    #define CYAN    "\033[36m"      /* Cyan */
    #define WHITE   "\033[37m"      /* White */
    
#endif
    
private:
    
    double probe_radius;
    
public:
    
    Chanalyzer()
    {
        probe_radius=0;
    }
    
    void set_probeRadius(double pr);
    double get_probeRadius();
    
    void Edelsbrunner(Rt& reg_triang, Fixed_alpha_shape_3& alpha_shape) const;
    
    void find_mouths(int seed, Fixed_alpha_shape_3& alpha_shape, std::vector<int>& flagged_tetrahedra, std::map<int, Cell_handle>& int2tetrahedron, std::map<Cell_handle, int>& tetrahedron2int,  std::map<int, Vertex_handle>& int2vertex, std::map<Vertex_handle, int>& vertex2int, std::map<std::vector<int>, std::vector<Weighted_point>>& tetrahedra2commonVertices, std::map<std::vector<int>, std::vector<Vertex_handle>>& tetrahedra2commonVertexHandles, std::vector<std::vector<int>>& mouths_atoms_idx, std::string output_folder) const;
    
};
#endif

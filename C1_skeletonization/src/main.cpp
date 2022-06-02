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


namespace PMP = CGAL::Polygon_mesh_processing;

class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
                
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

int main(int argc,char* argv[])
{
    
    
    // INPUT PATH:
    InputParser input_parser(argc, argv);
    const std::string &str_input_path = input_parser.getCmdOption("-i");
    std::string input_path = str_input_path.c_str();
    std::cout << input_path << std::endl;
    if (input_path.empty())
    {
        std::cout << "Please provide path to input file!";
        return 0;
    }
    
    
    // LOAD MESH:
    std::ifstream input(input_path + "channel.off");
    Polyhedron tmesh;
    input >> tmesh;
    if (!CGAL::is_triangle_mesh(tmesh))
    {
        std::cout << "Input geometry is not triangulated." << std::endl;
    }
    
    
    // MESH PREPROCESSING
    // Incrementally fill the holes
    unsigned int nb_holes = 0;
    for(Halfedge_handle_P h : halfedges(tmesh))
    {
        if(h->is_border())
        {
            std::vector<Facet_handle_P>  patch_facets;
            std::vector<Vertex_handle_P> patch_vertices;
            bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(tmesh,
                                                                             h,
                                                                             std::back_inserter(patch_facets),
                                                                             std::back_inserter(patch_vertices),
                                                                             PMP::parameters::vertex_point_map(get(CGAL::vertex_point, tmesh))
                                                                             .geom_traits(Kernel())));
            std::cout << " Number of facets in constructed patch: " << patch_facets.size() << std::endl;
            std::cout << " Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
            std::cout << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
            ++nb_holes;
        }
    }
    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;
    std::ofstream out("./output/CGAL_filled.off");
    out.precision(17);
    out << tmesh << std::endl;
    out.close();
    
    
    // SKELETONIZATION:
    Skeleton skeleton;
    Skeletonization mcs(tmesh);
    mcs.set_medially_centered_speed_tradeoff(1);
    mcs.contract_until_convergence();
    mcs.convert_to_skeleton(skeleton);
    std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
    std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
    // Output all the edges of the skeleton:
    std::ofstream output("./output/CGAL_skel-poly.polylines.txt");
    for(Skeleton_edge e : CGAL::make_range(edges(skeleton)))
    {
        const Point& s = skeleton[source(e, skeleton)].point;
        const Point& t = skeleton[target(e, skeleton)].point;
        output << "2 "<< s << " " << t << "\n";
    }
    
    output.close();
    // Output skeleton points:
    output.open("./output/CGAL_skeleton.txt");
    for(Skeleton_vertex v : CGAL::make_range(vertices(skeleton)))
        for(vertex_descriptor vd : skeleton[v].vertices)
            output << skeleton[v].point  << "\n";
    output.close();
    
}

    
    



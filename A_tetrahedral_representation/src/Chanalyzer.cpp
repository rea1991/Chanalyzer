
#include "Chanalyzer.h"

extern double cc_vol_min;
extern double cc_vol_max;
extern double cc_sur_min;
extern double cc_sur_max;
extern double cc_com_min;
extern double cc_com_max;
extern int num_candidates2return;

void Chanalyzer :: set_probeRadius(double pr)
{
    probe_radius = pr;
}

//

double Chanalyzer :: get_probeRadius()
{
    return probe_radius;
}

//

void Chanalyzer :: Edelsbrunner(Rt& reg_triang, Fixed_alpha_shape_3& alpha_shape) const{
    
    std::cout << INFO << GREEN << "\n\n\n START OF EDELSBRUNNER ALGORITHM" << RESET << std::endl;
    
    
    //1) Pre-processing...
    std::cout << INFO << GREEN << "1) Triangulation pre-processing..." << RESET;
    
    //1.0) Defining some variables of interest...
    
    //Some counters:
    int I1=0, I2=0, I3=0, I4=0, I5=0;
    
    //Undirected graph encoding CGAL tetrahedra (finite AND infinite):
    Graph regular_tetrahedra(alpha_shape.number_of_cells());
    
    //Map that associates, for each integer, the respective tetrahedron:
    std::map<int, Cell_handle> int2tetrahedron;
    std::map<Cell_handle, int> tetrahedron2int;
    
    //Map that associates, for each integer, the respective vertex:
    std::map<int, Vertex_handle> int2vertex;
    std::map<Vertex_handle, int> vertex2int;
    
    //Flags to identify a tetrahedron:
    //-)If it is infinite, the respective entry is -1;
    //-)If it is finite, the respective entry represents the number of appeareance (according to the CGAL ordering).
    std::vector<int> is_finite(alpha_shape.number_of_cells(),-1);
    
    //Flags to identify connected components in Edelsbrunners output...clarified later!
    std::vector<int> flagged_tetrahedra(alpha_shape.number_of_cells(),0);
    
    //Map that stores a couple of integers (identifying two adjacent cells) and
    //their common vertex handles:
    std::map<std::vector<int>, std::vector<Weighted_point>> tetrahedra2commonVertices;
    std::map<std::vector<int>, std::vector<Vertex_handle>> tetrahedra2commonVertexHandles;
    
    //Map each facet to its type w.r.t. the alpha-shape
    //(we use the fact that each facets is uniquely determined by its barycenter):
    std::map<std::vector<double>, int> facet2type;
    
    //Map each edge to its type w.r.t. the alpha-shape
    //(we use the fact that each edge is uniquely determined by its barycenter):
    std::map<std::vector<double>, int> edge2type;
    
    //List of orthogonal spheres:
    std::vector<Weighted_point> orthogonal_spheres={};
    
    //1.1.1) Saving facet type...
    for(Finite_Facets_Iterator ffit = alpha_shape.finite_facets_begin();ffit!=alpha_shape.finite_facets_end();ffit++)
    {
        const Facet& facetp = (*ffit);
        //Vertices of facetp:
        const Weighted_point& p1 = facetp.first->vertex((facetp.second+1)&3)->point();
        const Weighted_point& p2 = facetp.first->vertex((facetp.second+2)&3)->point();
        const Weighted_point& p3 = facetp.first->vertex((facetp.second+3)&3)->point();
        
        int ctype = alpha_shape.classify(*ffit);
        if(ctype == Fixed_alpha_shape_3::INTERIOR)
            I1++;
        if(ctype == Fixed_alpha_shape_3::EXTERIOR)
            I2++;
        if(ctype == Fixed_alpha_shape_3::REGULAR)
            I3++;
        if(ctype == Fixed_alpha_shape_3::SINGULAR)
            I4++;
        
        facet2type.insert(std::pair<std::vector<double>, int>({
            (p1.x() + p2.x() + p3.x())/3,
            (p1.y() + p2.y() + p3.y())/3,
            (p1.z() + p2.z() + p3.z())/3}, alpha_shape.classify((*ffit))));
    }
    std::cout << std::endl;
    std::cout << INFO << "Some information on facets..." << std::endl;
    std::cout << INFO << "*) There are " << I1 << " INTERIOR facets" << std::endl;
    std::cout << INFO << "*) There are " << I2 << " EXTERIOR facets" << std::endl;
    std::cout << INFO << "*) There are " << I3 << " REGULAR facets" << std::endl;
    std::cout << INFO << "*) There are " << I4 << " SINGULAR facets" << std::endl;
    std::cout << INFO << "*) Total number of finite facets: " << I1+I2+I3+I4 << std::endl;
    
    //1.1.2) Saving edge type...
    I1=0, I2=0, I3=0, I4=0, I5=0;
    
    for(Finite_Edges_Iterator feit = alpha_shape.finite_edges_begin();feit!=alpha_shape.finite_edges_end();feit++)
    {
        const Weighted_point& p1 = feit->first->vertex(feit->second)->point();
        const Weighted_point& p2 = feit->first->vertex(feit->third)->point();
        
        int ctype = alpha_shape.classify((*feit));
        
        if(ctype == Fixed_alpha_shape_3::INTERIOR)
            I1++;
        if(ctype == Fixed_alpha_shape_3::EXTERIOR)
            I2++;
        if(ctype == Fixed_alpha_shape_3::REGULAR)
            I3++;
        if(ctype == Fixed_alpha_shape_3::SINGULAR)
            I4++;
        
        edge2type.insert(std::pair<std::vector<double>, int>({
            (p1.x() + p2.x())/2,
            (p1.y() + p2.y())/2,
            (p1.z() + p2.z())/2}, alpha_shape.classify((*feit))));
    }
    std::cout << INFO << "Some information on edges..." << std::endl;
    std::cout << INFO << "*) There are " << I1 << " INTERIOR edges" << std::endl;
    std::cout << INFO << "*) There are " << I2 << " EXTERIOR edges" << std::endl;
    std::cout << INFO << "*) There are " << I3 << " REGULAR edges" << std::endl;
    std::cout << INFO << "*) There are " << I4 << " SINGULAR edges" << std::endl;
    std::cout << INFO << "*) Total number of finite edges: " << I1+I2+I3+I4 << std::endl;
    
    //1.2) Defining int2tetrahedron, is_finite, flagged_tetrahedra, facet2type, and tetrahedra2commonVertices...
    Cell_iterator fcit = alpha_shape.cells_begin();
    I1=0, I2=0;
    //For any tetrahedron I1>=0:
    for(;fcit!=alpha_shape.cells_end();fcit++)
    {
        
        //Bijectively map an integer to a tetrahedron:
        int2tetrahedron.insert(std::pair<int, Cell_handle>(I1,fcit));
        tetrahedron2int.insert(std::pair<Cell_handle, int>(fcit,I1));
        
        //Checking if tetrahedron I1 is finite
        Weighted_point v;
        if(alpha_shape.is_infinite(fcit)==0)
        {
            orthogonal_sphere(v, fcit->vertex(0)->point(), fcit->vertex(1)->point(), fcit->vertex(2)->point(), fcit->vertex(3)->point());
            is_finite[I1]=I2;
            orthogonal_spheres.push_back(v);
            I2++;
        }
        else
            flagged_tetrahedra[I1] = -1;
        I1++;
    }
    
    
    //Re-setting some auxiliary variables and defining a new one:
    I1=0, I2=0, I3=0, I4=0;
    
    //For any tetrahedron:
    fcit = alpha_shape.cells_begin();
    for(;fcit!=alpha_shape.cells_end();fcit++)
    {
        
//        std::cout << I5 << "/" << alpha_shape.number_of_cells() << std::endl;
//        I5++;
        
        //Index of current cell:
        int current_cell_idx = tetrahedron2int[fcit];
        
        // Finding all (unique) cells incident to each vertex of current cell (iterating through i)
        std::vector<int> incident_cells_idxs = {};
        for (int i=0; i<4; i++)
        {
            std::vector<Cell_handle> ll;
            alpha_shape.incident_cells(fcit->vertex(i),std::back_inserter(ll));
            
            // Checking each incident cell (iterating through j), pt. 1
            for (int j=0; j<ll.size(); j++)
                incident_cells_idxs.push_back(tetrahedron2int[ll[j]]);
        }
        sort(incident_cells_idxs.begin(), incident_cells_idxs.end());
        incident_cells_idxs.erase(unique(incident_cells_idxs.begin(), incident_cells_idxs.end()), incident_cells_idxs.end() );
        
        // Checking each incident cell (iterating through j), pt. 2
        for (int j=0; j<incident_cells_idxs.size(); j++)
        {
            int current_incident_cell_idx = incident_cells_idxs[j];
            Cell_handle current_incident_cell = int2tetrahedron[incident_cells_idxs[j]];
            
            // Checking how many vertices of the incident cell are in the original one (iterating through k)
            // Check performed only once per pair!
            if(current_cell_idx<current_incident_cell_idx)
            {
                // Store the common vertices:
                std::vector<Weighted_point> common_vertices={};
                std::vector<Vertex_handle> common_vertexHandles={};
                for (int k=0; k<4; k++)
                    if (fcit->has_vertex(current_incident_cell->vertex(k)))
                    {
                        common_vertices.push_back(current_incident_cell->vertex(k)->point());
                        common_vertexHandles.push_back(current_incident_cell->vertex(k));
                    }
                
                //If such a facet is EXTERIOR, then add edges connecting the tetrahedra:
                if (common_vertexHandles.size()==3)
                {
                    if(facet2type[{
                        (common_vertexHandles[0]->point().x() + common_vertexHandles[1]->point().x() + common_vertexHandles[2]->point().x())/3,
                        (common_vertexHandles[0]->point().y() + common_vertexHandles[1]->point().y() + common_vertexHandles[2]->point().y())/3,
                        (common_vertexHandles[0]->point().z() + common_vertexHandles[1]->point().z() + common_vertexHandles[2]->point().z())/3
                        }]==Fixed_alpha_shape_3::EXTERIOR)
                    {
                        regular_tetrahedra.addEdge(current_cell_idx,current_incident_cell_idx);
                        regular_tetrahedra.addEdge(current_incident_cell_idx,current_cell_idx);
                        I1++;
                    }
                    else
                        I2++;
                    tetrahedra2commonVertices.insert(std::pair<std::vector<int>, std::vector<Weighted_point>>({current_cell_idx,current_incident_cell_idx},common_vertices));
                    tetrahedra2commonVertices.insert(std::pair<std::vector<int>, std::vector<Weighted_point>>({current_incident_cell_idx,current_cell_idx},common_vertices));
                    tetrahedra2commonVertexHandles.insert(std::pair<std::vector<int>, std::vector<Vertex_handle>>({current_cell_idx,current_incident_cell_idx},common_vertexHandles));
                    tetrahedra2commonVertexHandles.insert(std::pair<std::vector<int>, std::vector<Vertex_handle>>({current_incident_cell_idx,current_cell_idx},common_vertexHandles));
                }

                //If the tetrahedra have an edge in common....
                if(common_vertexHandles.size()==2)
                {
                    //... and such an edge is EXTERIOR, then add edges connecting the tetrahedra:
                    if(
                       edge2type[{
                           (common_vertexHandles[0]->point().x() + common_vertexHandles[1]->point().x())/2,
                           (common_vertexHandles[0]->point().y() + common_vertexHandles[1]->point().y())/2,
                           (common_vertexHandles[0]->point().z() + common_vertexHandles[1]->point().z())/2
                           }]==Fixed_alpha_shape_3::EXTERIOR)
                    {
                        regular_tetrahedra.addEdge(current_cell_idx,current_incident_cell_idx);
                        regular_tetrahedra.addEdge(current_incident_cell_idx,current_cell_idx);
                        I3++;
                    }
                    else
                        I4++;
                }
            }
            
        }
        
    }
    
    
    std::cout << INFO << "Some information regarding connection of adjacent tetrahedra..." << std::endl;
    std::cout << INFO << "*) Connection by EXTERIOR CGAL facets: " << I1 << "/" << I1+I2 << "!" << std::endl;
    std::cout << INFO << "*) Connection by EXTERIOR CGAL edges: " << I3 << "/" << I3+I4 << "!" << std::endl;
    std::cout << INFO << "Done!" << std::endl;
    
    
    
    
    //1.3) Defining int2vertex
    Finite_Vertex_Iterator fvit = alpha_shape.finite_vertices_begin();
    I1=0;
    for(;fvit!=alpha_shape.finite_vertices_end();fvit++)
    {
        int2vertex.insert(std::pair<int, Vertex_handle>(I1,fvit));
        vertex2int.insert(std::pair<Vertex_handle, int>(fvit,I1));
        I1++;
    }
    
    


    //2) Compute or load the flow relation:
    auto step2_start = std::chrono::high_resolution_clock::now();

    //Directed graph for flow relation:
    Graph Edelsbrunner(alpha_shape.number_of_cells());


    std::cout << INFO << GREEN << "2) Computing or loading the flow relation..." << RESET << std::endl;
    for(I1=0; I1<alpha_shape.number_of_cells(); I1++)
    {
        //For any finite cell I1:
        if(is_finite[I1]>-1)
        {

            Cell_handle current_cell_handle = int2tetrahedron[I1];

            Weighted_point v11= (*current_cell_handle).vertex(0)->point();
            Weighted_point v21= (*current_cell_handle).vertex(1)->point();
            Weighted_point v31= (*current_cell_handle).vertex(2)->point();
            Weighted_point v41= (*current_cell_handle).vertex(3)->point();

            //Barycenter of tetrahedron I1:
            double vb_x = (v11.x()+v21.x()+v31.x()+v41.x())/4;
            double vb_y = (v11.y()+v21.y()+v31.y()+v41.y())/4;
            double vb_z = (v11.z()+v21.z()+v31.z()+v41.z())/4;
            Point3 vb = Point3(vb_x, vb_y, vb_z);


            for (int I2=0; I2<4; I2++)
            {
                Cell_handle neighbor_I2 = current_cell_handle->neighbor(I2);
                if(tetrahedra2commonVertices.find({I1,tetrahedron2int[neighbor_I2]})!=tetrahedra2commonVertices.end())
                    if(flow_relation(vb, orthogonal_spheres[is_finite[I1]], tetrahedra2commonVertices[{I1,tetrahedron2int[neighbor_I2]}]))
                        Edelsbrunner.addEdge(tetrahedron2int[neighbor_I2],I1);
            }

        }

        //Connect any pair of infinite tetrahedra with edges having weight 0...
        //NB. these edges are only intended to apply Dijkstra algorithm only once!
        if(is_finite[I1]==-1)
        {

            Cell_handle current_cell_handle = int2tetrahedron[I1];
            for (int I2=0; I2<4; I2++)
            {
                Cell_handle neighbor_I2 = current_cell_handle->neighbor(I2);
                if (is_finite[tetrahedron2int[neighbor_I2]]==-1 && I1!=tetrahedron2int[neighbor_I2])
                {
                    Edelsbrunner.addEdge(I1,tetrahedron2int[neighbor_I2],0);
                    Edelsbrunner.addEdge(tetrahedron2int[neighbor_I2],I1,0);
                }
            }

        }

    }
    Edelsbrunner.save_graph("flow_relation.txt");
    auto step2_stop = std::chrono::high_resolution_clock::now();
    auto step2_duration = std::chrono::duration_cast<std::chrono::seconds>(step2_stop - step2_start);
    std::cout << INFO << "Time taken to compute the flow relation... "
    << step2_duration.count() << " seconds!" << std::endl;
    std::cout << INFO << "Done!" << std::endl;


    //3) Removing ancestors of infinity and alpha-shape...
    std::cout << INFO << GREEN << "3) Cleaning up regular triangulation..." << RESET << std::endl;
    std::cout << INFO << "Edelsbrunner's graph has " << Edelsbrunner.getNumberOfVertices() << " vertices." << std::endl;

    //Store the index of one of the infinity tetrahedra (in alpha_shape):
    //(useful to get all ancestors of infinity by Dijkstra's algorithm).
    int s;
    for(I1=0; I1<alpha_shape.number_of_cells(); I1++)
        if(is_finite[I1]==-1)
        {
            s = I1;
            break;
        }

    //Dijkstra algorithm to find the ancestors of infinity:
    vector<int> path(Edelsbrunner.getNumberOfVertices());
    vector<int> dist = Edelsbrunner.dijkstraDist(s, path);

    //Removing infinite ancestors...
    //*)flagged_tetrahedra[i]=-1 identifies infinite tetrahedra;
    //*)flagged_tetrahedra[i]=-2 identifies their parents;
    //*)flagged_tetrahedra[i]=-3 identifies the parents of the parents, i.e., the grandparents of the infinity tetrahedra;
    //*)...
    std::vector<bool> visited(flagged_tetrahedra.size(), 0);
    for (int i = 0; i < dist.size(); i++)
        if(dist[i] != infi)
        {
            flagged_tetrahedra[i] = -dist[i]-1;
            visited[i] = 1;
        }

    //Removing tetrahedra from the alpha shape:
    fcit = alpha_shape.cells_begin();
    I1=0;
    for(;fcit!=alpha_shape.cells_end();fcit++){
        //Check if a tetrahedron belongs to the alpha-shape and if it is infinite:
        int ctype = alpha_shape.classify(fcit);
        if(ctype != Fixed_alpha_shape_3::EXTERIOR)
        {
            flagged_tetrahedra[I1] = 1;
            visited[I1] = 1;
        }
        I1++;
    }
    std::cout << INFO << "Done!" << std::endl;


    //4) Edelsbrunner surviving tetrahedra:
    std::cout << INFO << GREEN << "4) Saving Edelsbrunner's surviving tetrahedra..." << RESET << std::endl;
    std::string str = "./ED/ED.OFF";
    write_flagged_tetrahedra2off(alpha_shape,flagged_tetrahedra, str, 0);
    std::cout << INFO << "Done!" << std::endl;


    //5) Computing of connected components:
    std::cout << INFO << GREEN << "5) Computing the connected components..." << RESET << std::endl;
    int K=1;
    for (int i=0; i<flagged_tetrahedra.size(); i++)
    {
        if(flagged_tetrahedra[i]==0)
        {
            K++;
            regular_tetrahedra.DFS(i, visited, flagged_tetrahedra, K);
        }
    }
    K--;

    //Re-flagging the connected components...
    for (int i=0; i<flagged_tetrahedra.size(); i++)
        if(flagged_tetrahedra[i]>0)
            flagged_tetrahedra[i]--;
    std::cout << INFO << "The number of total (Edelsbrunnian) connected components is " << K << std::endl;
    std::cout << INFO << "Done!" << std::endl;


    std::cout << GREEN << "END OF EDELSBRUNNER ALGORITHM\n" << RESET << std::endl << std::endl;


    /// POST-PROCESSING
    std::cout << INFO << GREEN << "START OF POST-PROCESSING...\n" << RESET;


    //1) Computing surface area, volume, and compactness...
    std::vector<double> volume_CCs(K,0);
    std::vector<double> surface_CCs(K,0);
    std::vector<double> compactness_CCs(K,0);

    for (int k=1; k<K+1; k++)
    {
        volume_CCs[k-1]  = vol_tets(int2tetrahedron, flagged_tetrahedra, k);
        std::vector<int> indices_pocket = {};
        for(int i=0; i<flagged_tetrahedra.size(); i++)
            if(flagged_tetrahedra[i]==k)
                indices_pocket.push_back(i);
        surface_CCs[k-1] = sur_tets(indices_pocket, int2tetrahedron, tetrahedra2commonVertices);
        compactness_CCs[k-1] = pow(surface_CCs[k-1],3)/(pow(volume_CCs[k-1],2)*36*M_PI);
    }


    //2) Saving the above-computed descriptors...
    std::ofstream fout;

    fout.open("./ED/shape_descriptors/ED_volumes.txt");
    for (int l = 0; l < volume_CCs.size(); l++)
        fout << volume_CCs.at(l) << std::endl;
    fout.close();

    fout.open("./ED/shape_descriptors/ED_surface_areas.txt");
    for (int l = 0; l < surface_CCs.size(); l++)
        fout << surface_CCs.at(l) << std::endl;
    fout.close();

    fout.open("./ED/shape_descriptors/ED_compactness.txt");
    for (int l = 0; l < compactness_CCs.size(); l++)
        fout << compactness_CCs.at(l) << std::endl;
    fout.close();


    //3) Saving only those connected components within certain ranges...
    std::cout << INFO << "Saving those connected components such that:" << std::endl;
    std::cout << INFO << "*) Volume is between "<< cc_vol_min << " and " << cc_vol_max << ";" <<  std::endl;
    std::cout << INFO << "*) Surface area is between "<< cc_sur_min << " and " << cc_sur_max << ";" << std::endl;
    std::cout << INFO << "*) Compactness between "<< cc_com_min << " and " << cc_com_max << "." << std::endl;
    
    // Find the n CCs with highest volume:
    std::multimap<double, std::size_t> multiVolume_CCs;
    for (std::size_t i = 0; i != volume_CCs.size(); ++i)
        multiVolume_CCs.insert({volume_CCs[i], i});

    std::vector<double> sorted_Volume_CCs;
    std::vector<std::size_t> indices_sort;
    for (const auto & kv : multiVolume_CCs) {
        sorted_Volume_CCs.push_back(kv.first);
        indices_sort.push_back(kv.second);
    }
    std::cout << RED << std::endl << "Listing the first XX connected components w.r.t. volume:" << RESET <<std::endl;
    for (int i=indices_sort.size()-1; i>indices_sort.size()-(num_candidates2return+1); i--)
    {
        int k = indices_sort[i]+1;

        if(
           volume_CCs[k-1]>=cc_vol_min && volume_CCs[k-1]<=cc_vol_max &&
           surface_CCs[k-1]>=cc_sur_min && surface_CCs[k-1]<=cc_sur_max &&
           compactness_CCs[k-1]>=cc_com_min && compactness_CCs[k-1]<=cc_com_max)
        {
            str = "./ED/ED_CCs/CC_" + to_string((indices_sort.size()-1)-i) + ".OFF";
            write_flagged_tetrahedra2off(alpha_shape,flagged_tetrahedra, str, k);

            //Summary...
            std::cout << INFO << RED << to_string((indices_sort.size()-1)-i) << ") Connected component k="<< k << RESET << std::endl;
            std::cout << INFO << "Volume is \t\t" << volume_CCs[k-1] << std::endl;
            std::cout << INFO << "Surface area is \t" << surface_CCs[k-1] << std::endl;
            std::cout << INFO << "Compactness is \t" << compactness_CCs[k-1] << std::endl;

            //I want to save a txt which contains the atom centers and radii...
            std::vector<std::vector<double>> pocket_atoms {};
//            std::ofstream output_file;
//            output_file.open("./ED/ED_atoms/CC_" + to_string((indices_sort.size()-1)-i) + "_" +  to_string(k) + ".xyzr");
            for (int ii=0; ii<flagged_tetrahedra.size(); ii++)
            {
                if(flagged_tetrahedra[ii]==k)
                {
                    for (int jj=0;jj<4;jj++)
                    {
                        Cell& cellp = *int2tetrahedron[ii];
                        pocket_atoms.push_back({
                            cellp.vertex(jj)->point().x(),
                            cellp.vertex(jj)->point().y(),
                            cellp.vertex(jj)->point().z(),
                            sqrt(cellp.vertex(jj)->point().weight())-probe_radius});
                    }
                }
            }
            std::sort(pocket_atoms.begin(), pocket_atoms.end());
            pocket_atoms.erase(std::unique(pocket_atoms.begin(), pocket_atoms.end()), pocket_atoms.end());
//            for (int ii=0; ii<pocket_atoms.size(); ii++)
//                output_file << pocket_atoms[ii][0] << " " << pocket_atoms[ii][1] << " " << pocket_atoms[ii][2] << " " << pocket_atoms[ii][3] << std::endl;
//            output_file.close();

//            //I also want to save a txt which contains the indices of tetrahedra belonging to each connected component...
//            std::vector<int> pocket_tetrahedra {};
//            output_file.open("./ED/ED_tet_index/CC_" + to_string((indices_sort.size()-1)-i) + "_" + to_string(k) + ".txt");
//            for (int ii=0; ii<flagged_tetrahedra.size(); ii++)
//                if(flagged_tetrahedra[ii]==k)
//                    output_file << ii << std::endl;
//            output_file.close();
            
            // Save the tetrahedra involved in each CC:
            std::ofstream output_file_tet;
            output_file_tet.open("ED/ED_tet/CC_" + to_string((indices_sort.size()-1)-i) + ".txt");
            
            for (int ii=0; ii<flagged_tetrahedra.size(); ii++)
            {
                if(flagged_tetrahedra[ii]==k)
                {
                    for (int jj=0;jj<4;jj++)
                    {
                        Cell& cellp = *int2tetrahedron[ii];
                        pocket_atoms.push_back({
                            cellp.vertex(jj)->point().x(),
                            cellp.vertex(jj)->point().y(),
                            cellp.vertex(jj)->point().z(),
                            sqrt(cellp.vertex(jj)->point().weight())-probe_radius});
                        output_file_tet << cellp.vertex(jj)->point().x() << " " << cellp.vertex(jj)->point().y() << " " << cellp.vertex(jj)->point().z() << " ";
                    }
                    output_file_tet << std::endl;
                }
            }
            output_file_tet.close();

//            // Save the tetrahedra involved in each CC:
//            std::ofstream output_file_facets;
//            output_file_facets.open("./ED/ED_tet_vert/CC_" + to_string((indices_sort.size()-1)-i) + "_" + to_string(k) + ".txt");
//
//            for (int i=0; i<flagged_tetrahedra.size(); i++)
//            {
//                if(flagged_tetrahedra[i]==k)
//                {
//                    for (int j=0;j<4;j++)
//                    {
//                        Cell& cellp = *int2tetrahedron[i];
//                        pocket_atoms.push_back({cellp.vertex(j)->point().x(),
//                            cellp.vertex(j)->point().y(),
//                            cellp.vertex(j)->point().z(),
//                            sqrt(cellp.vertex(j)->point().weight())-probe_radius});
//                        output_file_facets << cellp.vertex(j)->point().x() << " " << cellp.vertex(j)->point().y() << " " << cellp.vertex(j)->point().z() << " ";
//                    }
//                    output_file_facets << std::endl;
//                }
//            }
//            output_file_facets.close();
            
            
            // Mouths of a specific connected component
            int seed;
            std::cout << INFO << "Finding mouths..." << RESET << std::endl;
            std::vector<std::vector<int>> mouths_atoms_idx = {};
            for (int ii=0; ii<flagged_tetrahedra.size(); ii++)
                //if (flagged_tetrahedra[ii]==indices_sort[indices_sort.size()-1]+1)
                if (flagged_tetrahedra[ii]==k)
                {
                    seed = ii;
                    break;
                }
            std::cout << INFO << "Current seed is " << seed << std::endl;
            std::cout << INFO << "Current flag is " << flagged_tetrahedra[seed] << std::endl;
            write_flagged_tetrahedra2off(alpha_shape,flagged_tetrahedra, "ED/mouths/" + to_string((indices_sort.size()-1)-i) + "/CC.off", flagged_tetrahedra[seed]);
            find_mouths(seed, alpha_shape, flagged_tetrahedra, int2tetrahedron, tetrahedron2int, int2vertex, vertex2int, tetrahedra2commonVertices, tetrahedra2commonVertexHandles, mouths_atoms_idx, "ED/mouths/" + to_string((indices_sort.size()-1)-i) + "/");
            std::cout << std::endl;

        }
    }
    
    std::cout << std::endl << INFO << GREEN << "END OF POST-PROCESSING...\n" << RESET;
    
}

//

void Chanalyzer :: find_mouths(int seed, Fixed_alpha_shape_3& alpha_shape, std::vector<int>& flagged_tetrahedra, std::map<int, Cell_handle>& int2tetrahedron, std::map<Cell_handle, int>& tetrahedron2int, std::map<int, Vertex_handle>& int2vertex, std::map<Vertex_handle, int>& vertex2int, std::map<std::vector<int>, std::vector<Weighted_point>>& tetrahedra2commonVertices, std::map<std::vector<int>, std::vector<Vertex_handle>>& tetrahedra2commonVertexHandles, std::vector<std::vector<int>>& mouths_atoms_idx, std::string output_folder) const{
    
    
    // (1) Find tetrahedra in the connected component of interest
    int FoI = flagged_tetrahedra[seed];
    //std::cout << "Flag of interest (FoI) is " << FoI << std::endl;
    std::vector<Cell_handle> tetrahedra = {};
    for(int i=0; i<flagged_tetrahedra.size(); i++)
        if(flagged_tetrahedra[i]==FoI)
            tetrahedra.push_back(int2tetrahedron[i]);
    
    
    // (2) Find all facets forming the mouths
    std::vector <std::vector <Weighted_point>> mouth_facets = {};
    std::vector <std::vector <Vertex_handle>>  mouth_facetHandles = {};
    std::vector<int> tetrahedra_outerPeel_idx(alpha_shape.number_of_cells(),0); // Identifies tetrahedra for mouth!
    for (int i=0; i<tetrahedra.size(); i++)
    {
        Cell_handle current_cell_handle = tetrahedra[i];
        for (int k=0; k<4; k++)
        {
            Cell_handle neighbor_k = current_cell_handle->neighbor(k);
            if (flagged_tetrahedra[tetrahedron2int[neighbor_k]]<0)      // be careful, <0 is advised!
            {
                tetrahedra_outerPeel_idx[tetrahedron2int[neighbor_k]]=1;
                mouth_facets.push_back(tetrahedra2commonVertices[{tetrahedron2int[current_cell_handle],tetrahedron2int[neighbor_k]}]);
                mouth_facetHandles.push_back(tetrahedra2commonVertexHandles[{tetrahedron2int[current_cell_handle],tetrahedron2int[neighbor_k]}]);
            }
        }
    }
    write_flagged_tetrahedra2off(alpha_shape,tetrahedra_outerPeel_idx, output_folder + "tetrahedra_outer_peel.off", 1);
    
    
    // (3) Find mouths (i.e., connected components in the set of facets we have just computed)
    Graph mouths(mouth_facets.size());
    for (int i=0; i<mouth_facets.size(); i++)
    {
        std::vector <Weighted_point> T_i = mouth_facets[i];
        for (int j=i+1; j<mouth_facets.size(); j++)
        {
            std::vector <Weighted_point> T_j = mouth_facets[j];
            int number_commonVertices = 0;
            
            for (int k=0; k<3; k++)
                if (count(T_j.begin(), T_j.end(), T_i[k])>0)
                    number_commonVertices++;
            
            if (number_commonVertices==2){
                mouths.addEdge(i,j);
                mouths.addEdge(j,i);
            }
            
        }
    }
    
    std::vector<int>  flagged_facets(mouth_facets.size(), 0);
    std::vector<bool> visited_facets(mouth_facets.size(), 0);
    int number_of_mouths=0;
    for (int i=0; i<flagged_facets.size(); i++)
    {
        if(flagged_facets[i]==0)
        {
            number_of_mouths++;
            mouths.DFS(i, visited_facets, flagged_facets, number_of_mouths);
        }
    }
    std::cout << INFO << "Total number of facets determining the mouths is " << mouth_facets.size() << std::endl;
    std::cout << INFO << "The number of mouths is " << number_of_mouths << std::endl;
    
    
    // (4) Saving mouth atoms for the connected components
    // N.B. each triplet of rows identify a facet!
    std::ofstream file_mouth_atoms;
    file_mouth_atoms.open(output_folder + "mouth_atoms_" + to_string(seed) + ".txt");
    std::vector<int> mouth_atoms_idx = {};
    
    int max_flag = -1;  // auxiliary variable, to construct mouths_atoms_idx
    for (int i=0; i<mouth_facets.size(); i++)
    {
        
        //Recover the three points corresponding to facet i:
        std::vector<Vertex_handle>& common_vertexHandles = mouth_facetHandles[i];
        Weighted_point v1 = common_vertexHandles[0]->point();
        Weighted_point v2 = common_vertexHandles[1]->point();
        Weighted_point v3 = common_vertexHandles[2]->point();
        mouth_atoms_idx.push_back(vertex2int[common_vertexHandles[0]]);
        mouth_atoms_idx.push_back(vertex2int[common_vertexHandles[1]]);
        mouth_atoms_idx.push_back(vertex2int[common_vertexHandles[2]]);
        
        // Save atoms to file:
        file_mouth_atoms << v1.x() << " " << v1.y() << " " << v1.z() << " " << sqrt(v1.weight())-probe_radius << std::endl;
        file_mouth_atoms << v2.x() << " " << v2.y() << " " << v2.z() << " " << sqrt(v2.weight())-probe_radius << std::endl;
        file_mouth_atoms << v3.x() << " " << v3.y() << " " << v3.z() << " " << sqrt(v3.weight())-probe_radius << std::endl;
        
        if (max_flag<flagged_facets[i])
        {
            max_flag = flagged_facets[i];
            mouths_atoms_idx.push_back({
                vertex2int[common_vertexHandles[0]],
                vertex2int[common_vertexHandles[1]],
                vertex2int[common_vertexHandles[2]]});
        }
        else
        {
            mouths_atoms_idx[flagged_facets[i]-1].push_back(vertex2int[common_vertexHandles[0]]);
            mouths_atoms_idx[flagged_facets[i]-1].push_back(vertex2int[common_vertexHandles[1]]);
            mouths_atoms_idx[flagged_facets[i]-1].push_back(vertex2int[common_vertexHandles[2]]);
        }
        
    }
    file_mouth_atoms.close();
    std::sort(mouth_atoms_idx.begin(), mouth_atoms_idx.end());
    mouth_atoms_idx.resize(std::distance(mouth_atoms_idx.begin(), std::unique(mouth_atoms_idx.begin(), mouth_atoms_idx.end())));
    
    //for (int i=0; i<mouth_atoms_idx.size(); i++)
    //    std::cout << "Atom " << mouth_atoms_idx[i] << std::endl;
    
    std::vector<Vertex_handle> mouth_atoms = {};
    for (int i=0; i<mouth_atoms_idx.size(); i++)
        mouth_atoms.push_back(int2vertex[mouth_atoms_idx[i]]);
    
    // (5) Save mouths to OFF files
    std::ofstream file_mouths;
    for (int i=1; i<=number_of_mouths; i++){
        
        std::cout << INFO << "Mouth number " << i-1 << " has " << mouths_atoms_idx[i-1].size() << " facets!" << std::endl;
        
        file_mouths.open(output_folder + "_mouth_" + to_string(i) + ".off");
        
        std::vector<int> current_mouth_atoms_idx = mouths_atoms_idx[i-1];
        std::sort(current_mouth_atoms_idx.begin(), current_mouth_atoms_idx.end());
        current_mouth_atoms_idx.resize(std::distance(current_mouth_atoms_idx.begin(), std::unique(current_mouth_atoms_idx.begin(), current_mouth_atoms_idx.end())));
        
        int number_of_mouth_facets = count(flagged_facets.begin(), flagged_facets.end(), i);
        file_mouths << "OFF" << std::endl;
        file_mouths << current_mouth_atoms_idx.size() << " " << number_of_mouth_facets << " 0" << std::endl;
        
        for(int j=0; j<current_mouth_atoms_idx.size(); j++)
        {
            Vertex_handle current_vertex = int2vertex[current_mouth_atoms_idx[j]];
            file_mouths << current_vertex->point().x() << " " << current_vertex->point().y() << " " << current_vertex->point().z() << std::endl;
        }
        
        for(int j=0; j<mouths_atoms_idx[i-1].size(); j++)
        {
            
            int current_vertex_idx = mouths_atoms_idx[i-1][j];
            //current_mouth_atoms_idx
            std::vector<int>::iterator it = std::find(current_mouth_atoms_idx.begin(), current_mouth_atoms_idx.end(), current_vertex_idx);
            int current_atom_idx = std::distance(current_mouth_atoms_idx.begin(), it);
            
            if (j%3==0)
                file_mouths << "3 " << current_atom_idx << " ";
            else if (j%3==1)
                file_mouths << current_atom_idx << " ";
            else
                file_mouths << current_atom_idx << std::endl;
        }
        
        file_mouths.close();
    }
}


#include "Auxiliary.h"

//

bool orthogonal_sphere(Weighted_point& v, Weighted_point& v1, Weighted_point& v2, Weighted_point& v3, Weighted_point& v4){
    /*info:
     */
    
    double x1 = v1.x(), x2 = v2.x(), x3 = v3.x(), x4 = v4.x();
    double y1 = v1.y(), y2 = v2.y(), y3 = v3.y(), y4 = v4.y();
    double z1 = v1.z(), z2 = v2.z(), z3 = v3.z(), z4 = v4.z();
    
    double t1 = -(std::pow(x1,2.0) + std::pow(y1,2.0) + std::pow(z1,2.0));
    double t2 = -(std::pow(x2,2.0) + std::pow(y2,2.0) + std::pow(z2,2.0));
    double t3 = -(std::pow(x3,2.0) + std::pow(y3,2.0) + std::pow(z3,2.0));
    double t4 = -(std::pow(x4,2.0) + std::pow(y4,2.0) + std::pow(z4,2.0));
    
    double T=0;
    T=T-x2*(y3*z4-y4*z3)+x3*(y2*z4-y4*z2)-x4*(y2*z3-y3*z2);
    T=T+x1*(y3*z4-y4*z3)-x3*(y1*z4-y4*z1)+x4*(y1*z3-y3*z1);
    T=T-x1*(y2*z4-y4*z2)+x2*(y1*z4-y4*z1)-x4*(y1*z2-y2*z1);
    T=T+x1*(y2*z3-y3*z2)-x2*(y1*z3-y3*z1)+x3*(y1*z2-y2*z1);
    
    double D=0;
    D=D-t2*(y3*z4-y4*z3)+t3*(y2*z4-y4*z2)-t4*(y2*z3-y3*z2);
    D=D+t1*(y3*z4-y4*z3)-t3*(y1*z4-y4*z1)+t4*(y1*z3-y3*z1);
    D=D-t1*(y2*z4-y4*z2)+t2*(y1*z4-y4*z1)-t4*(y1*z2-y2*z1);
    D=D+t1*(y2*z3-y3*z2)-t2*(y1*z3-y3*z1)+t3*(y1*z2-y2*z1);
    D=D/T;
    
    double E=0;
    E=E-x2*(t3*z4-t4*z3)+x3*(t2*z4-t4*z2)-x4*(t2*z3-t3*z2);
    E=E+x1*(t3*z4-t4*z3)-x3*(t1*z4-t4*z1)+x4*(t1*z3-t3*z1);
    E=E-x1*(t2*z4-t4*z2)+x2*(t1*z4-t4*z1)-x4*(t1*z2-t2*z1);
    E=E+x1*(t2*z3-t3*z2)-x2*(t1*z3-t3*z1)+x3*(t1*z2-t2*z1);
    E=E/T;
    
    double F=0;
    F=F-x2*(y3*t4-y4*t3)+x3*(y2*t4-y4*t2)-x4*(y2*t3-y3*t2);
    F=F+x1*(y3*t4-y4*t3)-x3*(y1*t4-y4*t1)+x4*(y1*t3-y3*t1);
    F=F-x1*(y2*t4-y4*t2)+x2*(y1*t4-y4*t1)-x4*(y1*t2-y2*t1);
    F=F+x1*(y2*t3-y3*t2)-x2*(y1*t3-y3*t1)+x3*(y1*t2-y2*t1);
    F=F/T;
    
    double G=0;
    G=G-t1*(z2*(x3*y4-x4*y3)-z3*(x2*y4-x4*y2)+z4*(x2*y3-x3*y2));
    G=G+t2*(z1*(x3*y4-x4*y3)-z3*(x1*y4-x4*y1)+z4*(x1*y3-x3*y1));
    G=G-t3*(z1*(x2*y4-x4*y2)-z2*(x1*y4-x4*y1)+z4*(x1*y2-x2*y1));
    G=G+t4*(z1*(x2*y3-x3*y2)-z2*(x1*y3-x3*y1)+z3*(x1*y2-x2*y1));
    G=G/T;
    
    double x = -D/2;
    double y = -E/2;
    double z = -F/2;
    double r = 0.5*pow( D*D + E*E + F*F - 4*G, 0.5);
    v=Weighted_point(Point3(x,y,z),(r*r));
    
    const bool FAIL =0;
    const bool SUCCESS =1;
    
    
    
    return SUCCESS;
}

//

bool flow_relation(Point3& vb,  Weighted_point& oc, std::vector<Weighted_point>& common_vertices) {


    //We must check whether the interior of tetrahedron I1 and its orthogonal center lies on
    //opposite sides with respect to the plane defined by the three common vertices.

    double cv1_x = common_vertices[0].x();
    double cv1_y = common_vertices[0].y();
    double cv1_z = common_vertices[0].z();
    double cv2_x = common_vertices[1].x();
    double cv2_y = common_vertices[1].y();
    double cv2_z = common_vertices[1].z();
    double cv3_x = common_vertices[2].x();
    double cv3_y = common_vertices[2].y();
    double cv3_z = common_vertices[2].z();

    //Orthogonal center of first tetrahedron:
    double voc_x = oc.x();
    double voc_y = oc.y();
    double voc_z = oc.z();

    //Barycenter of first tetrahedron:
    double vb_x = vb.x();
    double vb_y = vb.y();
    double vb_z = vb.z();

    //Plane passing through the three common vertices:
    double a = (cv2_y-cv1_y)*(cv3_z-cv1_z)-(cv2_z-cv1_z)*(cv3_y-cv1_y);
    double b = -((cv2_x-cv1_x)*(cv3_z-cv1_z)-(cv2_z-cv1_z)*(cv3_x-cv1_x));
    double c = (cv2_x-cv1_x)*(cv3_y-cv1_y)-(cv2_y-cv1_y)*(cv3_x-cv1_x);
    double d = -(a*cv1_x+b*cv1_y+c*cv1_z);

    //Check if flow relation holds:
    if((a*vb_x + b*vb_y + c*vb_z + d)*(a*voc_x + b*voc_y + c*voc_z + d)<0)
        return 1;
    else
        return 0;
}

//

double vol_tet(Weighted_point& v1, Weighted_point& v2, Weighted_point& v3, Weighted_point& v4) {

    double V= abs(
                  (v1.z()-v2.z())*(v3.x()*v4.y()-v4.x()*v3.y())+
                  (v3.z()-v1.z())*(v2.x()*v4.y()-v4.x()*v2.y())+
                  (v1.z()-v4.z())*(v2.x()*v3.y()-v3.x()*v2.y())+
                  (v2.z()-v3.z())*(v1.x()*v4.y()-v4.x()*v1.y())+
                  (v4.z()-v2.z())*(v1.x()*v3.y()-v3.x()*v1.y())+
                  (v3.z()-v4.z())*(v1.x()*v2.y()-v2.x()*v1.y())
                  );
    return V/6;
}

//

double surf_tet(Weighted_point& v1, Weighted_point& v2, Weighted_point& v3, Weighted_point& v4) {

    double S = 0;
    double a, b, c, s;

    //Heron's formula...

    //Facet {v1,v2,v3}
    a = sqrt(pow(v2.x()-v1.x(),2) + pow(v2.y()-v1.y(),2) + pow(v2.z()-v1.z(),2));
    b = sqrt(pow(v3.x()-v1.x(),2) + pow(v3.y()-v1.y(),2) + pow(v3.z()-v1.z(),2));
    c = sqrt(pow(v3.x()-v2.x(),2) + pow(v3.y()-v2.y(),2) + pow(v3.z()-v2.z(),2));
    s = 0.5*(a + b + c);
    S += sqrt(s*(s-a)*(s-b)*(s-c));

    //Facet {v1,v2,v4}
    a = sqrt(pow(v2.x()-v1.x(),2) + pow(v2.y()-v1.y(),2) + pow(v2.z()-v1.z(),2));
    b = sqrt(pow(v4.x()-v1.x(),2) + pow(v4.y()-v1.y(),2) + pow(v4.z()-v1.z(),2));
    c = sqrt(pow(v4.x()-v2.x(),2) + pow(v4.y()-v2.y(),2) + pow(v4.z()-v2.z(),2));
    s = 0.5*(a + b + c);
    S += sqrt(s*(s-a)*(s-b)*(s-c));

    //Facet {v1,v3,v4}
    a = sqrt(pow(v3.x()-v1.x(),2) + pow(v3.y()-v1.y(),2) + pow(v3.z()-v1.z(),2));
    b = sqrt(pow(v4.x()-v1.x(),2) + pow(v4.y()-v1.y(),2) + pow(v4.z()-v1.z(),2));
    c = sqrt(pow(v4.x()-v3.x(),2) + pow(v4.y()-v3.y(),2) + pow(v4.z()-v3.z(),2));
    s = 0.5*(a + b + c);
    S += sqrt(s*(s-a)*(s-b)*(s-c));

    //Facet {v2,v3,v4}
    a = sqrt(pow(v3.x()-v2.x(),2) + pow(v3.y()-v2.y(),2) + pow(v3.z()-v2.z(),2));
    b = sqrt(pow(v4.x()-v2.x(),2) + pow(v4.y()-v2.y(),2) + pow(v4.z()-v2.z(),2));
    c = sqrt(pow(v4.x()-v3.x(),2) + pow(v4.y()-v3.y(),2) + pow(v4.z()-v3.z(),2));
    s = 0.5*(a + b + c);
    S += sqrt(s*(s-a)*(s-b)*(s-c));

    return S;

}

//

double sur_tets(std::map<int, Cell_handle>& int2tetrahedron, std::vector<int>& flagged_tetrahedra, std::map<std::vector<int>, std::vector<Weighted_point>>& tetrahedra2commonVertices, int flag) {

    double surface_tetrahedra = 0;
    for(int i=0; i<flagged_tetrahedra.size(); i++)
        if(flagged_tetrahedra[i]==flag)
        {

            Cell& cellp = (*int2tetrahedron[i]);
            Weighted_point v1 = cellp.vertex(0)->point();
            Weighted_point v2 = cellp.vertex(1)->point();
            Weighted_point v3 = cellp.vertex(2)->point();
            Weighted_point v4 = cellp.vertex(3)->point();

            //Total surface -- upper bound...
            surface_tetrahedra += surf_tet(v1, v2, v3, v4);

        }


    //Remove internal facets from surface area computation...
    for(int i=0; i<flagged_tetrahedra.size(); i++)
    {
        for(int j=i+1; j<flagged_tetrahedra.size(); j++)
        {

            if(flagged_tetrahedra[i]==flag && flagged_tetrahedra[j]==flag)
            {
                //If I1 and I2 share a common facet:
                if(tetrahedra2commonVertices.find({i,j})!=tetrahedra2commonVertices.end())
                {

                    std::vector<Weighted_point>& common_vertices = tetrahedra2commonVertices[{i,j}];
                    if (common_vertices.size()>2){
                        Weighted_point v1 = common_vertices[0];
                        Weighted_point v2 = common_vertices[1];
                        Weighted_point v3 = common_vertices[2];
                        //Heron's formula...
                        double a = sqrt(pow(v2.x()-v1.x(),2) + pow(v2.y()-v1.y(),2) + pow(v2.z()-v1.z(),2));
                        double b = sqrt(pow(v3.x()-v1.x(),2) + pow(v3.y()-v1.y(),2) + pow(v3.z()-v1.z(),2));
                        double c = sqrt(pow(v3.x()-v2.x(),2) + pow(v3.y()-v2.y(),2) + pow(v3.z()-v2.z(),2));
                        double s = 0.5*(a + b + c);
                        surface_tetrahedra -= 2*sqrt(s*(s-a)*(s-b)*(s-c));
                    }

                }
            }
        }
    }

    return surface_tetrahedra;

}

//

double sur_tets(std::vector <int> indices_setOfTetrahedra, std::map<int, Cell_handle>& int2tetrahedron, std::map<std::vector<int>, std::vector<Weighted_point>>& tetrahedra2commonVertices) {

    double surface_tetrahedra = 0;
    for(int i=0; i<indices_setOfTetrahedra.size(); i++)
    {
        Cell& cellp = (*int2tetrahedron[indices_setOfTetrahedra[i]]);
        Weighted_point v1 = cellp.vertex(0)->point();
        Weighted_point v2 = cellp.vertex(1)->point();
        Weighted_point v3 = cellp.vertex(2)->point();
        Weighted_point v4 = cellp.vertex(3)->point();

        //Total surface -- upper bound...
        surface_tetrahedra += surf_tet(v1, v2, v3, v4);
    }

    //Remove internal facets from surface area computation...
    for(int i=0; i<indices_setOfTetrahedra.size(); i++)
    {
        for(int j=i+1; j<indices_setOfTetrahedra.size(); j++)
        {

            //If I1 and I2 share a common facet:
            if(tetrahedra2commonVertices.find({indices_setOfTetrahedra[i],indices_setOfTetrahedra[j]})!=tetrahedra2commonVertices.end())
            {
                std::vector<Weighted_point>& common_vertices = tetrahedra2commonVertices[{indices_setOfTetrahedra[i],indices_setOfTetrahedra[j]}];
                Weighted_point v1 = common_vertices[0];
                Weighted_point v2 = common_vertices[1];
                Weighted_point v3 = common_vertices[2];
                //Heron's formula...
                double a = sqrt(pow(v2.x()-v1.x(),2) + pow(v2.y()-v1.y(),2) + pow(v2.z()-v1.z(),2));
                double b = sqrt(pow(v3.x()-v1.x(),2) + pow(v3.y()-v1.y(),2) + pow(v3.z()-v1.z(),2));
                double c = sqrt(pow(v3.x()-v2.x(),2) + pow(v3.y()-v2.y(),2) + pow(v3.z()-v2.z(),2));
                double s = 0.5*(a + b + c);
                surface_tetrahedra -= 2*sqrt(s*(s-a)*(s-b)*(s-c));

            }
        }
    }


    return surface_tetrahedra;
}

//

double vol_tets(std::map<int, Cell_handle>& int2tetrahedron, std::vector<int>& flagged_tetrahedra, int flag) {

    double volume_tetrahedra = 0;
    for(int i=0; i<flagged_tetrahedra.size(); i++)
        if(flagged_tetrahedra[i]==flag)
        {
            Cell& cellp = (*int2tetrahedron[i]);
            Weighted_point v1 = cellp.vertex(0)->point();
            Weighted_point v2 = cellp.vertex(1)->point();
            Weighted_point v3 = cellp.vertex(2)->point();
            Weighted_point v4 = cellp.vertex(3)->point();
            //Total volume
            volume_tetrahedra += vol_tet(v1, v2, v3, v4);
        }
    return volume_tetrahedra;
}

//

bool write_triangulation2off(Rt& reg_triang, std::string root_file_name, const int d) {
    /*info:
     Given a regular triangulation, the simplices of dimension d=0,1,2,3 are saved in a file (.off for d=0,1,3, .ply for d=2).
     For d=-1, we store all simplices (one file per dimension).
     */

    const bool FAIL =0;
    const bool SUCCESS =1;

    if(d==-1 || d==0){

        std::string str0 =  root_file_name + "simplices_0.OFF";
        FILE* fp = fopen(str0.c_str() ,"w");
        fprintf(fp,"COFF\n%d 0 0",((int)reg_triang.number_of_vertices()));

        //PRINT VERTICES:
        std :: map <Weighted_point,int> vert_map;
        Weighted_point v1,v2,v3;
        Finite_Vertex_Iterator fvit = reg_triang.finite_vertices_begin();
        int i = 0;
        for (;reg_triang.finite_vertices_end()!=fvit;fvit++){
            vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
            i++;
            fprintf(fp,"\n%f %f %f 0 255 255 255",fvit->point().x(),fvit->point().y(),fvit->point().z());
        }
        std :: cout << endl << INFO << "Done: Number of points: "<< i;
        rewind(fp);
        fclose(fp);

    }

    if(d==-1 || d==1){

        std::string str1 =  root_file_name + "simplices_1.ply";
        FILE* fp = fopen(str1.c_str() ,"w");
        fprintf(fp,"ply\n");
        fprintf(fp,"format ascii 1.0\n");
        fprintf(fp,"comment object: %d edges with color per vertex\n", ((int)reg_triang.number_of_finite_edges()));
        fprintf(fp,"element vertex %d\n", ((int)reg_triang.number_of_vertices()));
        fprintf(fp,"property float x\n");
        fprintf(fp,"property float y\n");
        fprintf(fp,"property float z\n");
        fprintf(fp,"property uchar red\n");
        fprintf(fp,"property uchar green\n");
        fprintf(fp,"property uchar blue\n");
        fprintf(fp,"element edge %d\n", ((int)reg_triang.number_of_finite_edges()));
        fprintf(fp,"property int vertex1\n");
        fprintf(fp,"property int vertex2\n");
        fprintf(fp,"end_header\n\n");

        //PRINT VERTICES:
        std :: map <Weighted_point,int> vert_map;
        Weighted_point v1,v2,v3;
        Finite_Vertex_Iterator fvit = reg_triang.finite_vertices_begin();
        int i = 0;
        for (;reg_triang.finite_vertices_end()!=fvit;fvit++){
            vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
            i++;
            fprintf(fp,"%f %f %f 255 0 0\n",fvit->point().x(),fvit->point().y(),fvit->point().z());
        }

        //PRINT EDGES:
        Finite_Edges_Iterator feit = reg_triang.finite_edges_begin();
        int j=0;
        for(;feit!=reg_triang.finite_edges_end();feit++){
            v1 = feit->first->vertex(feit->second)->point();
            v2 = feit->first->vertex(feit->third)->point();
            fprintf(fp,"\n%d %d",vert_map[v1],vert_map[v2]);
            if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end()) )
            {
                std :: cout << endl<< WARN << "ERROR, vertex not mapped\n";
                return FAIL;
            }
            j++;

        }
        std :: cout << endl << INFO << "Done: Number of edges: "<< j;

        rewind(fp);
        fclose(fp);
    }


    if(d==-1 || d==2){

        std::string str2 =  root_file_name + "simplices_2.OFF";
        FILE* fp = fopen(str2.c_str() ,"w");
        fprintf(fp,"COFF\n%d %d 0",((int)reg_triang.number_of_vertices()), ((int)reg_triang.number_of_finite_facets()) );

        //PRINT VERTICES:
        std :: map <Weighted_point,int> vert_map;
        Weighted_point v1,v2,v3;
        Finite_Vertex_Iterator fvit = reg_triang.finite_vertices_begin();
        int i = 0;
        for (;reg_triang.finite_vertices_end()!=fvit;fvit++){
            vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
            i++;
            fprintf(fp,"\n%f %f %f 0 255 0 255",fvit->point().x(),fvit->point().y(),fvit->point().z());
        }

        //PRINT FACETS:
        Finite_Facets_Iterator ffit = reg_triang.finite_facets_begin();
        int j=0;
        for(;ffit!=reg_triang.finite_facets_end();ffit++){
            const Facet& facetp = (*ffit);
            v1= facetp.first->vertex((facetp.second+1)&3)->point();
            v2= facetp.first->vertex((facetp.second+2)&3)->point();
            v3= facetp.first->vertex((facetp.second+3)&3)->point();
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v2], vert_map[v3] );
            if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
            {
                std :: cout << endl<< WARN << "ERROR, vertex not mapped\n";
                return FAIL;
            }
            j++;
        }
        std :: cout << endl << INFO << "Done: Number of facets: "<< j;

        rewind(fp);
        fclose(fp);

    }

    if(d==-1 || d==3){

        std::string str3 =  root_file_name + "simplices_3.OFF";
        FILE* fp = fopen(str3.c_str() ,"w");
        fprintf(fp,"COFF\n%d %d 0",((int)reg_triang.number_of_vertices()), 4*((int)reg_triang.number_of_finite_cells()) );

        //PRINT VERTICES:
        std :: map <Weighted_point,int> vert_map;
        Weighted_point v1,v2,v3,v4;
        Finite_Vertex_Iterator fvit = reg_triang.finite_vertices_begin();
        int i = 0;
        for (;reg_triang.finite_vertices_end()!=fvit;fvit++){
            vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
            i++;
            fprintf(fp,"\n%f %f %f 0 0 255 255",fvit->point().x(),fvit->point().y(),fvit->point().z());
        }

        //PRINT FACETS OF TETRAHEDRA:
        Finite_Cells_Iterator fcit = reg_triang.finite_cells_begin();
        int k=0;
        for(;fcit!=reg_triang.finite_cells_end();fcit++){
            const Cell& cellp = (*fcit);
            v1= cellp.vertex(0)->point();
            v2= cellp.vertex(1)->point();
            v3= cellp.vertex(2)->point();
            v4= cellp.vertex(3)->point();
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v2], vert_map[v3]);
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v2], vert_map[v4]);
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v3], vert_map[v4]);
            fprintf(fp,"\n3 %d %d %d",vert_map[v2],vert_map[v3], vert_map[v4]);
            if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end())||(vert_map.find(v4) == vert_map.end()) )
            {
                std :: cout << endl<< WARN << "ERROR, vertex not mapped\n";
                return FAIL;
            }
            k++;
        }
        std :: cout << endl << INFO << "Done: Number of cells: "<< k;

        rewind(fp);
        fclose(fp);

    }

    return SUCCESS;
}

//

bool write_alphaFilt(Fixed_alpha_shape_3& alpha, const int flag){
    /*info:
     CGAL classifiers:
     REGULAR = not singular AND on the boundary: that is alpha-exposed
     EXTERIOR = I think, not alpha-exposed and not interior.
     SINGULAR = Exterior and ONLY connected to other exterior elements
     INTERIOR = Within the boundary (REGULAR), not alpha exposed
     OTHER OPTIONS THAN REGULAR AND INTERIOR NOT SUPPORTED YET!
     Rationale:
     The function creates an OFF file. All possible vertex coordinates (independently from the CGAL classification)
     are printed and a map associating a unique label to each vertex is created.
     Then the facets corresponding to the chosen classification are printed on the file
     by listing the associated vertices labels (3 per facet, one facet per line)
     */
    const int EXTE = 0;
    const int SINGU = 1;
    const int REGU = 2;
    const int INTE = 3;
    const bool FAIL =0;
    const bool SUCCESS =1;
    char temp[BUFLEN];
    
    std :: cout << endl << INFO<< "Fetching alpha-simplex" ;
    std::string chosen_class;
    if(flag==REGU){
        std :: cout << endl << INFO<< "tag = REGULAR ";
        sprintf(temp,"# Boundary of the 0-alpha shape (alpha-exposed) with atomic radii augmented with the probe radius, if any.");
        chosen_class = "REGULAR";
    }
    else if(flag==INTE) {
        std :: cout << endl << INFO<< "tag = INTERNAL";
        sprintf(temp,"# Interior of the 0-alpha shape (alpha-exposed) with atomic radii augmented with the probe radius, if any.");
        chosen_class = "INTERIOR";
    }
    else if(flag==EXTE) {
        std :: cout << endl << INFO<< "tag = EXTERIOR";
        sprintf(temp,"# Exterior of the 0-alpha shape (alpha-exposed) with atomic radii augmented with the probe radius, if any.");
        chosen_class = "EXTERIOR";
    }
    else if(flag==SINGU) {
        std :: cout << endl << INFO<< "tag = SINGULAR";
        sprintf(temp,"# Singular of the 0-alpha shape (alpha-exposed) with atomic radii augmented with the probe radius, if any.");
        chosen_class = "SINGULAR";
    }
    std::string output_file_name = "./alpha_shape/" + chosen_class + ".OFF";
    FILE* fp = fopen(output_file_name.c_str(), "w");


    std::string title = temp;
    std :: string buff(title.size() + 15, ' ');

    fprintf(fp,"%s",buff.c_str());

    // here store unique map for vertices
    std :: map <Weighted_point,int> vert_map;
    // 3 vertices per triangle
    Weighted_point v1,v2,v3;

    Finite_Vertex_Iterator fvit = alpha.finite_vertices_begin();
    int i = 0;

    std :: cout << endl << INFO << "Building vertex list..";
    //creating "verbose vertex list" some won't be used in the OFF file
    // with some changes, should not be too hard to just print vertices that will be used..
    for (;alpha.finite_vertices_end()!=fvit;fvit++){
        int vtype = alpha.classify(fvit);
        vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
        i++;
        //save top half of OFF file, WHAT ABOUAT WEIGHTS?
        fprintf(fp,"\n %f %f %f",fvit->point().x(),fvit->point().y(),fvit->point().z());
    }
    std :: cout << endl << INFO << "Total Number of vertices: "<< i;

    // assign vertex to correct facet(s)
    std :: cout << endl << INFO << "Iterating through facets..";
    Finite_Facets_Iterator ffit = alpha.finite_facets_begin();
    int k=0;
    for(;ffit!=alpha.finite_facets_end();ffit++){
        int ttype = alpha.classify(*ffit);
        if(flag == REGU){
            if(ttype == Fixed_alpha_shape_3:: REGULAR){
                const Alpha_Facet& facetp = (*ffit);
                v1= facetp.first->vertex((facetp.second+1)&3)->point();
                v2= facetp.first->vertex((facetp.second+2)&3)->point();
                v3= facetp.first->vertex((facetp.second+3)&3)->point();

                if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
                {
                    std :: cout << endl<< WARN << "ERROR, vertex not mapped while retrieving alpha shape\n";
                    return FAIL;
                }
                fprintf(fp,"\n3 %d %d %d",vert_map.at(v1),vert_map.at(v2), vert_map.at(v3) );
                k++;
            }
        }
        else if (flag==INTE){
            if(ttype == Fixed_alpha_shape_3:: INTERIOR)
            {
                const Alpha_Facet& facetp = alpha.mirror_facet((*ffit)); //for consistency to what done when building patches
                //const Alpha_Facet& facetp = (*ffit);
                v1= facetp.first->vertex((facetp.second+1)&3)->point();
                v2= facetp.first->vertex((facetp.second+2)&3)->point();
                v3= facetp.first->vertex((facetp.second+3)&3)->point();

                if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
                {
                    std :: cout << endl << WARN << "ERROR, vertex not mapped while retrieving alpha shape\n";
                    return FAIL;
                }
                fprintf(fp,"\n3 %d %d %d",vert_map.at(v1),vert_map.at(v2), vert_map.at(v3) );
                k++;
            }
        }
        else if (flag==EXTE){
            if(ttype == Fixed_alpha_shape_3:: EXTERIOR)
            {
                const Alpha_Facet& facetp = alpha.mirror_facet((*ffit)); //for consistency to what done when building patches
                //const Alpha_Facet& facetp = (*ffit);
                v1= facetp.first->vertex((facetp.second+1)&3)->point();
                v2= facetp.first->vertex((facetp.second+2)&3)->point();
                v3= facetp.first->vertex((facetp.second+3)&3)->point();

                if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
                {
                    std :: cout << endl << WARN << "ERROR, vertex not mapped while retrieving alpha shape\n";
                    return FAIL;
                }
                fprintf(fp,"\n3 %d %d %d",vert_map.at(v1),vert_map.at(v2), vert_map.at(v3) );
                k++;
            }
        }
        else if (flag==SINGU){
            if(ttype == Fixed_alpha_shape_3:: SINGULAR)
            {
                const Alpha_Facet& facetp = alpha.mirror_facet((*ffit)); //for consistency to what done when building patches
                //const Alpha_Facet& facetp = (*ffit);
                v1= facetp.first->vertex((facetp.second+1)&3)->point();
                v2= facetp.first->vertex((facetp.second+2)&3)->point();
                v3= facetp.first->vertex((facetp.second+3)&3)->point();

                if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
                {
                    std :: cout << endl << WARN << "ERROR, vertex not mapped while retrieving alpha shape\n";
                    return FAIL;
                }
                fprintf(fp,"\n3 %d %d %d",vert_map.at(v1),vert_map.at(v2), vert_map.at(v3) );
                k++;
            }
        }

        else{
            std :: cout << endl << WARN<<"Invalid option for printing alpha shape.\n"<<std::endl;
            return FAIL;
        }
    }
    //move to the top and add extra line with number of vertices and facets, as required by OFF format
    rewind(fp);
    //fprintf(fp,"OFF %d %d 0\n",i,k);
    fprintf(fp,"%s\nOFF %d %d 0\n",title.c_str(),i,k);
    //fp << "OFF" << i<< k << std::endl;
    fclose(fp);
    //print off file
    std :: cout << endl << INFO << "Done: Number of facets with indicated attributes: "<< k << "\n "<< std ::endl;
    return SUCCESS;
}

//

bool write_flagged_tetrahedra2off(Fixed_alpha_shape_3& alpha_shape, std::vector<int> flagged_tetrahedra, std::string file_name, int K){

    const bool FAIL =0;
    const bool SUCCESS =1;


    //Number of remaining tetrahedra:
    int num_tetrahedra=0;
    for(int I=0;I<flagged_tetrahedra.size();I++){
        if(flagged_tetrahedra[I]==K)
            num_tetrahedra++;
    }
    //std::cout << INFO << "There are " << num_tetrahedra << " tetrahedra with flag " << K << std::endl;

    //Save into an .OFF file
    FILE* fp = fopen(file_name.c_str() ,"w");
    fprintf(fp,"COFF\n%d %d 0",((int)alpha_shape.number_of_vertices()), 4*num_tetrahedra );

    //PRINT VERTICES:
    std :: map <Weighted_point,int> vert_map;
    Weighted_point v1,v2,v3,v4;
    Finite_Vertex_Iterator fvit = alpha_shape.finite_vertices_begin();
    int i = 0;
    for (;alpha_shape.finite_vertices_end()!=fvit;fvit++){
        vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
        i++;
        fprintf(fp,"\n%f %f %f 0 0 255 255",fvit->point().x(),fvit->point().y(),fvit->point().z());
    }

    //PRINT FACETS OF TETRAHEDRA:
    Cell_iterator fcit = alpha_shape.cells_begin();
    int k=0;
    for(;fcit!=alpha_shape.cells_end();fcit++){
        if(flagged_tetrahedra[k]==K){

            const Cell& cellp = (*fcit);
            v1= cellp.vertex(0)->point();
            v2= cellp.vertex(1)->point();
            v3= cellp.vertex(2)->point();
            v4= cellp.vertex(3)->point();
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v2], vert_map[v3]);
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v2], vert_map[v4]);
            fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v3], vert_map[v4]);
            fprintf(fp,"\n3 %d %d %d",vert_map[v2],vert_map[v3], vert_map[v4]);
            if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) ==vert_map.end())||(vert_map.find(v4) == vert_map.end()) )
            {
                std :: cout << endl<< WARN << "ERROR, vertex not mapped\n";
                return FAIL;
            }
        }
        k++;
    }

    rewind(fp);
    fclose(fp);

    return SUCCESS;
}

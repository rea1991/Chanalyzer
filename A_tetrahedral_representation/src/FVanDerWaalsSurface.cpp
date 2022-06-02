	
#include "FVanDerWaalsSurface.h"

void FVanDerWaalsSurface::init()
{
	savePovRay = false;
	surfType = MOLECULAR_SURFACE;
	isRCbased = false;	
}

FVanDerWaalsSurface::FVanDerWaalsSurface():Surface()
{
	init();
}

FVanDerWaalsSurface::FVanDerWaalsSurface(DelPhiShared* ds)
{
	init();
	// set environment
	delphi = ds;	
}

void FVanDerWaalsSurface::init(ConfigFile* cf)
{
	bool buildStatus = cf->read<bool>( "Build_status_map", false );
	bool tri = cf->read<bool>( "Triangulation", false );

	if (!buildStatus && tri)
	{
		// status map is needed to deduced in/out vertices				
		cout << endl << ERR << "FVDW Surfce needs status map to perform triangulation";
		cout << endl << REMARK << "Please set Build_status_map = true";
		cout << endl;
		exit(-1);
	}
	setTriangulationFlag(false);	
	inside = 5;	
}

FVanDerWaalsSurface::FVanDerWaalsSurface(ConfigFile* cf,DelPhiShared* ds):Surface(cf)
{
	init();
	init(cf);
	// set environment
	delphi = ds;	
}


/** Colour the grid by using an obvious space fill model*/
bool FVanDerWaalsSurface::build()
{
	if (delphi == NULL)
	{
		cout << endl << WARN << "Cannot get surface without DelPhi environment!";
		return false;
	}

	int na = delphi->numAtoms;
	
	// force not accurate triangulation
	accurateTriangulation = false;

	#ifdef ENABLE_BOOST_THREADS
	int num_thd;		
	if (num_cores==-1)
		num_thd = na;		
	else
	{
		cout << endl << INFO << "Forcing user selected num threads " << num_cores;
		num_thd = MIN(num_cores,na);
	}
	#else
	int num_thd = 1;
	#endif

	int cycle = na/num_thd;
	int rem = na%num_thd;
			
	#ifdef ENABLE_BOOST_THREADS
	boost::thread_group thdGroup; 
	#endif
			
	// setup split
	int j=0;
					
	for (j=0;j<num_thd;j++)
	{
		int start = cycle*j;
		int stop = cycle*(j+1);
		#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&FVanDerWaalsSurface::painter,this,start,stop));
		#else
			painter(start,stop);
		#endif
	}
			
	// remaining part
	if (rem>0)
	{
		int start = cycle*j;
		int stop = start+rem;
		#ifdef ENABLE_BOOST_THREADS
			thdGroup.create_thread(boost::bind(&FVanDerWaalsSurface::painter,this,start,stop));
		#else
			painter(start,stop);
		#endif
	}
	// end setup
			
	#ifdef ENABLE_BOOST_THREADS
	// join
	thdGroup.join_all();
	#endif
	return true;
}	

void FVanDerWaalsSurface::painter(int start,int stop)
{
	double downx,downy,downz,upx,upy,upz;
	Atom** at = delphi->atoms;
	double xmin = delphi->xmin;
	double ymin = delphi->ymin;
	double zmin = delphi->zmin;
	double side = delphi->side;
	double hside = delphi->hside;
	//short*** status = delphi->status;
	short* status = delphi->status;
	double* x = delphi->x;
	double* y = delphi->y;
	double* z = delphi->z;
	int NX = delphi->nx, NY = delphi->ny, NZ = delphi->nz;
	double p[3];
	double* center, dist2;
	int* epsmap = delphi->epsmap;
	bool bstat = delphi->buildStatus;
	bool beps = delphi->buildEpsMap;

	// colour the grid
	for (int i=start;i<stop;i++)
	{
		double radius = at[i]->radius;
		double radius2 = at[i]->radius2;
		center = at[i]->pos;

		downx = center[0]-radius;
		downy = center[1]+radius;
		downz = center[2]-radius;

		upx = center[0]+radius;
		upy = center[1]-radius;
		upz = center[2]+radius;
	          
		// resolve which are the grid cubes that
		// are cut by the bounding box object. These
		// cubes see the cell bounding box    
		int ix_start = (int)rintp((downx-xmin)/side);
		int iy_start = (int)rintp((downy-ymin)/side);
		int iz_start = (int)rintp((downz-zmin)/side);
	     	     
		int ix_end = (int)rintp((upx-xmin)/side);
		int iy_end = (int)rintp((upy-ymin)/side);
		int iz_end = (int)rintp((upz-zmin)/side);

		if (ix_start<0) ix_start=0;
		if (iy_start<0) iy_start=0;
		if (iz_start<0) iz_start=0;

		if (ix_end<0) ix_end=0;
		if (iy_end<0) iy_end=0;
		if (iz_end<0) iz_end=0;

		if (ix_end>=(int)NX) ix_end=NX-1;
		if (iy_end>=(int)NY) iy_end=NY-1;	          
		if (iz_end>=(int)NZ) iz_end=NZ-1;	          

		if (ix_start>=(int)NX) ix_start=NX-1;
		if (iy_start>=(int)NY) iy_start=NY-1;	          
		if (iz_start>=(int)NZ) iz_start=NZ-1;	          

		for (int iz=iz_start;iz<=iz_end;iz++)
		 for (int iy=iy_start;iy>=iy_end;iy--)
		  for (int ix=ix_start;ix<=ix_end;ix++)
		  {
			p[0]=x[ix];
			p[1]=y[iy];
			p[2]=z[iz];
					
			// set inside for status map if not already inside
			//if (bstat && status[ix][iy][iz]!=STATUS_POINT_INSIDE)
			//if (bstat && (STATUSMAP(ix,iy,iz,NX,NY))!=STATUS_POINT_INSIDE)
			if (bstat && read3DVector<short>(status,ix,iy,iz,NX,NY,NZ)!=STATUS_POINT_INSIDE)
			{
				DIST2(dist2,p,center)
				//if (dist2<radius2)	status[ix][iy][iz]=STATUS_POINT_INSIDE;
				if (dist2<radius2)	
					//STATUSMAP(ix,iy,iz,NX,NY)=STATUS_POINT_INSIDE;
					write3DVector<short>(status,STATUS_POINT_INSIDE,ix,iy,iz,NX,NY,NZ);
			}
					
			// set inside for epsmap
			//if (beps && (EPSMAP(ix,iy,iz,0,NX,NY,NZ))!=inside)
			if (beps && read4DVector<int>(epsmap,ix,iy,iz,0,NX,NY,NZ,3)!=inside)
			{
				p[0]+=hside;
				DIST2(dist2,p,center)
				if (dist2<radius2) 
					//EPSMAP(ix,iy,iz,0,NX,NY,NZ)=inside;
					write4DVector<int>(epsmap,inside,ix,iy,iz,0,NX,NY,NZ,3);
			}

			//if (beps && (EPSMAP(ix,iy,iz,1,NX,NY,NZ))!=inside)
			if (beps && read4DVector<int>(epsmap,ix,iy,iz,1,NX,NY,NZ,3)!=inside)
			{
				p[0]-=hside;
				p[1]+=hside;
				DIST2(dist2,p,center)
				if (dist2<radius2) 
					//EPSMAP(ix,iy,iz,1,NX,NY,NZ)=inside;
					write4DVector<int>(epsmap,inside,ix,iy,iz,1,NX,NY,NZ,3);
			}

			//if (beps && (EPSMAP(ix,iy,iz,2,NX,NY,NZ))!=inside)
			if (beps && read4DVector<int>(epsmap,ix,iy,iz,2,NX,NY,NZ,3)!=inside)
			{						
				p[1]-=hside;
				p[2]+=hside;
				DIST2(dist2,p,center)
				if (dist2<radius2) 
					//EPSMAP(ix,iy,iz,2,NX,NY,NZ)=inside;
					write4DVector<int>(epsmap,inside,ix,iy,iz,2,NX,NY,NZ,3);
			}
		  }						
	}
}

/**TODO*/
bool FVanDerWaalsSurface::preBoundaryProjection()
{
	/*
	bool flag;
	// useful to project bgp
	if (projBGP)
	{
		flag = buildAuxiliaryGrid();
		if (!flag)
		{
			cout << endl << ERR << "Cannot build 3D auxiliary grid";
			return flag;
		}
	}	
	*/
	return true;
	
}


bool FVanDerWaalsSurface::save(char* fileName)
{
	return false;
}

bool FVanDerWaalsSurface::load(char* fileName)
{
	return false;
}

void FVanDerWaalsSurface::printSummary()
{
}

/** TODO*/
bool FVanDerWaalsSurface::getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)
{	
	/*
	// get the cells that are associated to this grid point
	// by querying the auxiliary grid
	double dist;
	set<int>cells;
	double hside = delphi->hside;

	double aposx = p[0];
	double aposy = p[1];
	double aposz = p[2];

	// move from delphi grid to auxiliary grid
	int irefx = (int)rintp((aposx-xmin)/side);
	int irefy = (int)rintp((aposy-ymin)/side);
	int irefz = (int)rintp((aposz-zmin)/side);

	int ixaux = irefx;
	int iyaux = irefy;
	int izaux = irefz;

	for (int i=0;i<ind[ixaux][iyaux][izaux];i++)
		cells.insert((GRID_CONNOLLY_CELL_MAP(ixaux,iyaux,izaux,i,nx,ny,nz)));

	// keep the nearest patch
	set<int>::iterator it;
	double locProj[3], minDist=INFINITY, locNorm[3];
	bool ff=false;

	dist = 0;
	locProj[0]=0;
	locProj[1]=0;
	locProj[2]=0;

	locNorm[0]=0;
	locNorm[1]=0;
	locNorm[2]=0;

	bool isFeasy = false;
	int ptype;

	for (it=cells.begin();it!=cells.end();it++)
	{	
		ConnollyCell* cc = sesComplex[(*it)];

		if (cc->patch_type==POINT_CELL)
		{
			PointCell* pc = (PointCell*)cc;
			projectToSphere(p,delphi->atoms[pc->id]->pos,delphi->atoms[pc->id]->radius,locProj,(double*)NULL,dist);
		}
		else if (cc->patch_type==REGULAR_FACE_CELL || cc->patch_type==SINGULAR_FACE_CELL)
		{
			FacetCell* fc = (FacetCell*)cc;
			projectToSphere(p,fc->center,probe_radius,locProj,(double*)NULL,dist);
		}
		else
		{
			EdgeCell* ec = (EdgeCell*)cc;
			projectToTorus(p,ec,locProj,(double*)NULL,dist);
		}
		
		bool flag = isFeasible(cc,locProj);

		if (!flag)
			continue;

		if (dist<minDist)
		{
			minDist = dist;
			(*proj1)=locProj[0];
			(*proj2)=locProj[1];
			(*proj3)=locProj[2];

			(*normal1)=locNorm[0];
			(*normal2)=locNorm[1];
			(*normal3)=locNorm[2];
			isFeasy = flag;
			ptype = cc->patch_type;
		}
	}

	if (cells.size()==0)
	{
		cout << endl << WARN << "Empty cell in getProjection!";
		return false;
	}

	// Connolly surface is not that nice. Projections can't be analytical everywhere.
	// We are trying to project a point into a singularity. Old DelPhi style.
	if (minDist == INFINITY)
	{
		bool fixed = false;

		double** sampledPoints;									
		int coiNum = 500;
		sampledPoints = allocateMatrix2D<double>(coiNum,3);

		for (it=cells.begin();it!=cells.end();it++)			
		{	
			ConnollyCell* cc = sesComplex[(*it)];					
			if (cc->patch_type==SINGULAR_EDGE_CELL ||  cc->patch_type==REGULAR_EDGE_CELL)
			{
				EdgeCell* ec = (EdgeCell*)cc;
				
				// try to manage singularity. Explicitly project to probes and keep the nearest one
				if (ec->isSelfIntersecting)
				{
					// identify the nearest the common probes and check the feasibility of the projection
					// for at least one of them.
					
					//  set of incident probes to a given edge
					FacetCell* fcv[MAX_INCIDENT_PROBES];					
					int index = 0;

					// manage SINGULAR EDGE?
					if (ec->patch_type == SINGULAR_EDGE_CELL)
					{
						{
						#ifdef ENABLE_BOOST_THREADS
						boost::mutex::scoped_lock scopedLock(mutex);
						#endif
						(*errorStream) << endl << WARN << "Singular edge projection";
						}
					}

					// get all the sourrounding probe stations
					if (ec->patch_type == REGULAR_EDGE_CELL)
					{
						vector<FacetCell*>& f1 = atomPatches[ec->id[0]]->incidentProbes;
										
						int ref1 = ec->id[0];
						int ref2 = ec->id[1];				
						unsigned int np = f1.size();
																			
						// get all the incident probes. Most of the cases will be 2.
						for (unsigned int indf=0;indf<np;indf++)
						{				
							int i1 = (f1)[indf]->id[0];
							int i2 = (f1)[indf]->id[1];
							int i3 = (f1)[indf]->id[2];

							if	(i1==ref1 && i2==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i2==ref1 && i1==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i1==ref1 && i3==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i3==ref1 && i1==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i2==ref1 && i3==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i3==ref1 && i2==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}					
						}
					}	
				
					getCoi(ec->center,ec->rcoi,sampledPoints,coiNum,ec->u,ec->v);
					
					// for all the points get projection. Keep the nearest feasible
					for (int kk=0;kk<coiNum;kk++)
					{
						projectToSphere(p,sampledPoints[kk],probe_radius,locProj,(double*)NULL,dist);
						for (int ll=0;ll<index;ll++)
						{						
							bool flag = isFeasible(fcv[ll],locProj);
							
							if (!flag)
								continue;

							if (dist<minDist)
							{
								minDist = dist;
								(*proj1)=locProj[0];
								(*proj2)=locProj[1];
								(*proj3)=locProj[2];
								(*normal1)=0;
								(*normal2)=0;
								(*normal3)=0;
								fixed = true;
							}
						}
					}					
				}				
			}			
		} // end singular projection
		deleteMatrix2D<double>(coiNum,sampledPoints);
		
		if (!fixed)
		{
			(*proj1)=p[0];
			(*proj2)=p[1];
			(*proj3)=p[2];

			{  
			#ifdef ENABLE_BOOST_THREADS
				boost::mutex::scoped_lock scopedLock(mutex);
			#endif			
			(*errorStream) << endl << WARN << "Approximating bgp with grid point "; 					
			}
		}
		else
		{
			// non tangent continuity has been properly managed. 
			// the obtained projection is the best thing allowed by Connolly definition.
		}	
	}
	*/
	return true;	
}

void FVanDerWaalsSurface::preProcessPanel()
{
}

void FVanDerWaalsSurface::postRayCasting()
{
}

void FVanDerWaalsSurface::getRayIntersection(double pa[3],double pb[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)
{
}
	
FVanDerWaalsSurface::~FVanDerWaalsSurface()
{
	clear();
}

void FVanDerWaalsSurface::clear()
{

}

void FVanDerWaalsSurface::saveSphere(ostream& of,double* center,double radius)
{
	char buff2[BUFLEN];
	sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf  pigment{color Green}}",center[0],center[1],center[2],radius);	  				
	of << buff2;
}


/**TODO: use it to reposition vertices! (avoid ray casting)*/
/*
bool ConnollySurface::raySphere(double* orig,double* dir,double* sphere_center,double sphere_radius,double* t1,double* t2)
{
	// perform sphere intersection test
	double A,B,C,temp[3],tt;
	A = DOT(dir,dir);
	SUB(temp,orig,sphere_center)
		B = DOT(temp,dir);
	B*=2;
	C= DOT(temp,temp)-sphere_radius*sphere_radius;
	double det = B*B-4*A*C;

	// no intersection
	if (det<0)
		return false;

	det = sqrt(det);

	(*t1)= (-B-det)/(2*A);
	(*t2)= (-B+det)/(2*A);

	if ((*t2)<(*t1))
	{
		tt = (*t1);
		(*t1)=(*t2);
		(*t2)=tt;
	}

	return true;
}
*/
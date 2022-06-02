#include "CoulombicSurface.h"

void CoulombicSurface::init()
{
	_cutoff = DEFAULT_CUTOFF_COULOMBIC;
	surfType = MOLECULAR_SURFACE;
	_iso = DEFAULT_ISO;
}

void CoulombicSurface::init(ConfigFile* cf)
{
}

CoulombicSurface::CoulombicSurface()
{
	init();
}

CoulombicSurface::CoulombicSurface(DelPhiShared* ds):MeshSurface(ds)
{	
	init();
}

CoulombicSurface::CoulombicSurface(ConfigFile* cf,DelPhiShared* ds):MeshSurface(cf,ds)
{
	init();
	init(cf);
}

CoulombicSurface::~CoulombicSurface()
{
	clear();
}

void CoulombicSurface::clear()
{
}

void CoulombicSurface::setIsoValue(double iso)
{		
	_iso=iso;
}

double CoulombicSurface::getIsoValue()
{
	return _iso;
}


bool CoulombicSurface::build()
{	
	// numgrid points that are neighbour of an atom 
	int numgrid = (int)(_cutoff/delphi->side+0.5);
	cout << endl << INFO << "Using cut-off " << _cutoff << " num neighbour grid points " << numgrid;

	if (scalarField!=NULL)
		deleteMatrix3D<double>(last_nx,last_ny,scalarField);
	
	scalarField = allocateMatrix3D<double>(delphi->nx,delphi->ny,delphi->nz);

	if (scalarField==NULL)
	{
		cout << endl << ERR << "Cannot allocate scalar field!";
		exit(-1);
	}

	for (int i=0;i<delphi->nx;i++)
		for (int j=0;j<delphi->nx;j++)
			for (int k=0;k<delphi->nx;k++)
				scalarField[i][j][k]=0.0;

	Atom** atoms = delphi->atoms;	
	printf("\n");
	int na = delphi->numAtoms;

	buildAtomsMap();
	// fill surface class structures
	for (int i=0;i<na;i++)
	{
		printf("\r%sCoulombic %.2f%%        ",INFO,((float)i+1)/na*100.0);
		double* pos = atoms[i]->pos;
		double r = atoms[i]->radius;
		double r2 = r*r;
		// get ref grid point
		int ix = (int)rintp((pos[0]-delphi->xmin)/delphi->side);
		int iy = (int)rintp((pos[1]-delphi->ymin)/delphi->side);
		int iz = (int)rintp((pos[2]-delphi->zmin)/delphi->side);

		int start_x = MAX(0,(ix-numgrid));
		int start_y = MAX(0,(iy-numgrid));
		int start_z = MAX(0,(iz-numgrid));

		int end_x = MIN((delphi->nx),(ix+numgrid));
		int end_y = MIN((delphi->ny),(iy+numgrid));
		int end_z = MIN((delphi->nz),(iz+numgrid));

		for (int ii=start_x;ii<end_x;ii++)
			for (int jj=start_y;jj<end_y;jj++)
				for (int kk=start_z;kk<end_z;kk++)
				{
					double dist2 = 0;
					double p[3];
					p[0] = delphi->x[ii]-delphi->hside;
					p[1] = delphi->y[jj]-delphi->hside;
					p[2] = delphi->z[kk]-delphi->hside;
					
					DIST2(dist2,p,pos);

					//dist2 = sqrt(dist2);						
					// TODO add charges weight
					bool inside = false;
					double signed_dist = 0;
					DIST2(signed_dist,pos,p)
					double rad = r;
					signed_dist-=(rad*rad);
					if (signed_dist<0)
						inside = true;
					// check if not visited and if inside vdw and protect it
					if (inside)
						scalarField[ii][jj][kk] = INSIDE;
					else
						scalarField[ii][jj][kk] += 1./dist2;
				}
	}	

	disposeAtomsMap();

	bool old = accurateTriangulation;
	accurateTriangulation = true;	
	isAvailableScalarField = true;
	triangulateSurface(_iso,"coulombic.off");

	//cout << endl << INFO << "Smoothing coulombic surface...";
	//smoothSurface("coulombic.off");
	
	isAvailableScalarField = false;
	accurateTriangulation = old;
	cout << "ok!";

	deleteMatrix3D<double>(delphi->nx,delphi->ny,scalarField);
	// no more necessary
	//scalarField = NULL;

	char buff[BUFLEN];
	//strcpy(buff,"triangulatedSurf_fixed.off");
	strcpy(buff,"coulombic.off");

	// the triangulation tools of surface where used. Now they are cleaned

	////////////// reset surface status ////////////////////////////
	if (triList.size()>0)
	{
		vector<int*>::iterator it;
		for (it = triList.begin();it!=triList.end();it++)
			deleteVector<int>((*it));
	}

	if (vertList.size()>0) 
	{
		vector<double*>::iterator it;
		for (it = vertList.begin();it!=vertList.end();it++)
			deleteVector<double>((*it));
	}

	triList.clear();
	vertList.clear();
	totalSurfaceArea=0;
	totalVolume=0;	
	//////////////////////////////////////////////////////////

	// load its mesh
	load(buff);
	// the rest continue as a mesh surface
	//bool flag = MeshSurface::build();
	MeshSurface::printSummary();
	return true;
}

void  CoulombicSurface::printSummary()
{	
	cout << endl << INFO << "Cut-off distance " << _cutoff << " [A]";
}
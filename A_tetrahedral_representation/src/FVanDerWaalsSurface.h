//---------------------------------------------------------
/**    @file	FVanDerWaalsSurface.h
*     @brief	FVanDerWaalsSurface.h is the header for CLASS
*               FVanDerWaalsSurface.cpp								*/
//---------------------------------------------------------

#ifndef FVanDerWaalsSurface_h
#define FVanDerWaalsSurface_h

#include "Surface.h"
#include "SurfaceFactory.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

/** @brief This class builds the Van Der Waals surface as the union of the set of spheres
in a fast way: the surface is not reduced in analytical patches (Voronoi diagram) so
rays cannot be cast, but instead the grid is only coloured based on in/out information. 
The surface is not analytically sampled thus the triangulation can be only non accurate (bisecting marching cubes). 
The advantage is the the phase of surface construction is not necessary and that, still, 
projection of boundary grid points is fully analytical. Only triangulation is affected by this approximation.
Regarding FD PDE solution no approximation is performed.

@author Sergio Decherchi 
@date 22/03/2013
*/
class FVanDerWaalsSurface: public Surface
{


private:

	bool savePovRay;

public:
	/** Default constructor*/
	FVanDerWaalsSurface();
	/** set DelPhi environment*/
	FVanDerWaalsSurface(DelPhiShared* ds);
	/** set configuration and DelPhi environment*/
	FVanDerWaalsSurface(ConfigFile* cf,DelPhiShared* ds);

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	/** Compute VDW surface. Call it after load*/
	virtual bool build();
	/** Save it in a simple ASCII format (.ses)*/
	virtual bool save(char* fileName);
	/**Load the surface from a file in .ses format*/
	virtual bool load(char* fileName);
	/** Print summary of VDWS*/
	virtual void printSummary();
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3);
	/** Not implemented */
	virtual void getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals);
	/** function for the constructor without arguments*/
	virtual void init();
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);
	/**function for the denstructor*/
	virtual void clear();
	/** Not implmented*/
	virtual void preProcessPanel();
	/** Not implmented*/
	virtual void postRayCasting();
	virtual bool preBoundaryProjection();
	/////////////////////////////////////////////////////////////

	void setSavePovRay(bool ss) 
	{ 
		savePovRay = ss;
	}
	
	bool getSavePovRay() 
	{
		return savePovRay;
	}

	
	virtual ~FVanDerWaalsSurface();

private:

	/** save path (sphere) of that atom*/
	void saveSphere(ostream& of,double* center,double radius);
	/** the painter thread of the grid*/
	void painter(int start,int end);
	
};

//REGISTER_SURFACE(FVanDerWaalsSurface,"fvdw")
/*
static class FVanDerWaalsSurfaceRegister{ 
	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
		return new FVanDerWaalsSurface(conf,ds); 
	} 
	public: 
		FVanDerWaalsSurfaceRegister() 
		{ 
			surfaceFactory().add("fvdw",createSurface); 
		} 
} FVanDerWaalsSurfaceRegisterObject;
*/

static SurfaceRecorder<FVanDerWaalsSurface> fvdwRecorder("fvdw");

#endif
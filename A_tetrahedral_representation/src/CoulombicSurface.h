
//---------------------------------------------------------
/**    @file	CoulombicSurface.h
*     @brief	CoulombicSurface.h is the header for CLASS
                CoulombicSurface.cpp								*/
//---------------------------------------------------------

#ifndef CoulombicSurface_h
#define CoulombicSurface_h

#include "Surface.h"
#include "MeshSurface.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif

#define DEFAULT_CUTOFF_COULOMBIC 12

#define INSIDE (6)
#define DEFAULT_ISO (6-0.1)
/** @brief This class builds the Coulombic surface, triangulate it and the use mesh surface routines
for the rest of the work.

@author Sergio Decherchi 
@date 30/04/2013
*/
class CoulombicSurface: public MeshSurface
{
private:
	/** cut-off distance in Angstrom to speed-up computations. */
	double _cutoff;
	double _iso;
	
public:
	/** Default constructor*/
	CoulombicSurface();
	/** set DelPhi environment*/
	CoulombicSurface(DelPhiShared* ds);		
	/** set configuration and DelPhi environment*/
	CoulombicSurface(ConfigFile* cf,DelPhiShared* ds);	

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	/** Compute blobby mesh and run MeshSurface::build()*/
	virtual bool build();
	/** Print number a summary of the blobby*/
	virtual void printSummary();		
	/** function for the constructor without arguments*/
	virtual void init();
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);
	/**function for the denstructor*/
	virtual void clear();
	/////////////////////////////////////////////////////////////

	//////////////////////////// Coulombic SURFACE METHODS ///////////////////////////////
	void setIsoValue(double b);
	double getIsoValue();

	virtual ~CoulombicSurface();
};


static class CoulombicSurfaceRegister{ 
	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
		return new CoulombicSurface(conf,ds); 
	} 
	public: 
		CoulombicSurfaceRegister() 
		{ 
			surfaceFactory().add("coulombic",createSurface); 
		} 
} CoulombicSurfaceRegisterObject;

//static SurfaceRecorder<CoulombicSurface> coulombRecorder("coulombic");

#endif
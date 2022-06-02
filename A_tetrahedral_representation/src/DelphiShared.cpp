
#include "DelphiShared.h"

void DelPhiShared::init()
{
	epsmap = NULL;
	buildEpsMap = false;
	buildStatus = false;
	x = NULL;
	y = NULL;
	z = NULL;
	idebmap = NULL;
	atoms = NULL;
	scspos = NULL;
	scsnor = NULL;
	scsarea = NULL;
	ibgp = NULL;
	status = NULL;
	cavitiesVec = NULL;
	delphiBinding = false;
	multi_diel = false;
	tempStatus=NULL;
	isAvailableAtomInfo = false;
	atsurf = NULL;
}

DelPhiShared::DelPhiShared()
{
	init();
}

void DelPhiShared::init(double scale,double perfill, string fn,bool eps_flag,bool stat_flag,bool multi,bool atinfo)
{
	buildEpsMap = eps_flag;
	buildStatus = stat_flag;
	multi_diel = multi;
	isAvailableAtomInfo = atinfo;
	
	cout << endl << INFO << "Loading atoms....";

	bool flag = loadAtoms(fn);

	if (!flag)
	{
		cout << endl << ERR << "Missing or corrupt atoms file. Initialization failed";
		exit(-1);
	}

	flag = buildGrid(scale,perfill);
	if (!flag)
	{
		cout << endl << ERR << "Initialization failed";
		exit(-1);
	}
	cout << endl << INFO << "Initialization completed";
}

DelPhiShared::DelPhiShared(double scale,double perfill, string fn,bool eps_flag,bool stat_flag,bool multi,bool atinfo)
{
	init();
	init(scale,perfill,fn,eps_flag,stat_flag,multi,atinfo);
}

void DelPhiShared::DelPhiBinding(	double xmin,double ymin,double zmin,
									double xmax,double ymax,double zmax,
									double c1,double c2,double c3,double rmax,double perf,
									int* local_i_epsmap,int igrid,double scale,
									double* local_i_scspos,double* local_i_scsnor,
									bool* local_i_idebmap,int* local_i_ibgp,int maxbgp,
									bool bstatus,int* _atsurf)
{
	cout << endl << INFO << "DelPhi Binding...";	
	delphiBinding = true;
	buildStatus = bstatus;
	atsurf = _atsurf;
	
	if (x!=NULL)
		deleteVector<double>(x);
	if (y!=NULL)
		deleteVector<double>(y);
	if (z!=NULL)
		deleteVector<double>(z);

	if (buildStatus)
		if (status!=NULL) 
			deleteVector<short>(status);

	this->xmin = xmin;
	this->ymin = ymin;
	this->zmin = zmin;

	this->rmaxdim=rmax;
	this->perfill = perf;
	this->xmax = xmax;
	this->ymax = ymax;
	this->zmax = zmax;

	this->baricenter[0]=c1;
	this->baricenter[1]=c2;
	this->baricenter[2]=c3;

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	this->scale = scale;
	side = 1/scale;
	hside = side/2;
	A = side*side;
	nx = igrid;
	ny = igrid;
	nz = igrid;

	if (atoms!=NULL)
	{
		for (int i=0;i<numAtoms;i++)
			delete atoms[i];
		delete[] atoms;
	}

	atoms = NULL;
	
	// assuming all maps are already allocated before binding 
	this->epsmap  = local_i_epsmap;
	this->idebmap = local_i_idebmap;
	this->scsnor = local_i_scsnor;
	this->scsarea = NULL;
	this->scspos = local_i_scspos;
	this->ibgp = local_i_ibgp;

	int tot = nx*ny*nz;
	for (int i=0;i<tot;i++)
		idebmap[i]=true;

	for (int i=0;i<nx;i++)
		x[i]=xmin+i*side;

	for (int i=0;i<ny;i++)
		y[i]=ymin+i*side;

	for (int i=0;i<nz;i++)
		z[i]=zmin+i*side;

	// allocate grid point status memory
	// and initialize all to temporary outside status
	if (buildStatus)
	{
		/*
		status = allocateMatrix3D<short>(nx,ny,nz);
		for (int i=0;i<nx;i++)
			for (int j=0;j<ny;j++)
				for (int k=0;k<nz;k++)
					status[i][j][k]=2;
		*/
		int tot = nx*ny*nz;
		status = allocateVector<short>(tot);
		for (int i=0;i<tot;i++)
			status[i]=STATUS_POINT_TEMPORARY_OUT;
	}
	// set max bgp
	this->maxbgp = maxbgp;
	
	cout << "ok!";
	cout << endl << INFO << "Max number of bgps was set from DelPhi to " << maxbgp;
}


void DelPhiShared::finalizeBinding(int* ibnum)
{
	// update indexes to be conformant to Fortran
	for (int i=0;i<nbgp;i++)
	{
	  ibgp[3*i]++;
	  ibgp[3*i+1]++;
	  ibgp[3*i+2]++;
	}

	// update indexes to be conformant to Fortran
	if (atsurf!=NULL)
		for (int i=0;i<nbgp;i++)
			atsurf[i]++;

	epsmap = NULL; //avoid deletion
	scspos = NULL; //avoid deletion
	scsnor = NULL; //avoid deletion
	idebmap = NULL; //avoid deletion
	//scsarea = NULL; //for now not used, so it must be distroied
	ibgp = NULL;   //avoid deletion	
	// save ibnum
	(*ibnum) = nbgp;
}


bool DelPhiShared::loadAtoms(string fn)
{
	// load atoms file
    ifstream fin;
    fin.open(fn.c_str(),ios::in);
	char buffer[BUFLEN];
	double max_rad = 0;

	if (multi_diel)
		cout << endl << INFO << "Atom dielectric info is available..";

	if (isAvailableAtomInfo)
		cout << endl << INFO << "Atom Info is available..";

    if (fin.fail())
	{
		cout << endl << WARN << "Cannot read file " << fn;
		return false;
	}

	vector<Atom*> list;
	Atom* atm;

    while (!fin.eof( ))  
    {
        fin.getline(buffer,BUFLEN);
		if (strlen(buffer)>1)
		{
			atm = new Atom();

			if (!multi_diel && !isAvailableAtomInfo)
			{
				sscanf(buffer,"%lf %lf %lf %lf",&(atm->pos[0]),&(atm->pos[1]),&(atm->pos[2]),&(atm->radius));               
				atm->radius2 = (atm->radius)*(atm->radius);
			}
			else if (multi_diel)
			{
				int d = -1;
				sscanf(buffer,"%lf %lf %lf %lf %d",&(atm->pos[0]),&(atm->pos[1]),&(atm->pos[2]),&(atm->radius),&d);
				atm->radius2 = (atm->radius)*(atm->radius);
				if (d==-1)
				{
					cout << endl << ERR << "Cannot get the dielectric value of the atom number " << list.size() << "; please add it after radius entry to the xyzr input file";
					cout << endl;
					exit(-1);
				}
				else
					atm->dielectric = d;
			}
			else if (isAvailableAtomInfo)
			{
				int d = -1;
				char name[10],resName[10],chain[2];
				int resid;
				sscanf(buffer,"%lf %lf %lf %lf %d %s %s %d %s",&(atm->pos[0]),&(atm->pos[1]),&(atm->pos[2]),&(atm->radius),&d,name,resName,&resid,chain);
				atm->radius2 = (atm->radius)*(atm->radius);
				if (d==-1)
				{
					cout << endl << ERR << "Cannot get the dielectric value of the atom number " << list.size() << "; please add it after radius entry to the xyzr input file";
					cout << endl;
					exit(-1);
				}
				else
					atm->dielectric = d;			

				atm->ai.setName(string(name));
				atm->ai.setResName(string(resName));
				atm->ai.setResNum(resid);
				atm->ai.setChain(string(chain));
			}
			max_rad = MAX(atm->radius,max_rad);
			// debug
			//printf("\n %lf %lf %lf %lf",atm->pos[0],atm->pos[1],atm->pos[2],atm->radius);               
			list.push_back(atm);
		}
		else
		{
			// skip empty line
		}
    }
    fin.close();    
	
	// freeze list
	vector<Atom*>::iterator it = list.begin();
	if (atoms!=NULL)
	{
		for (int i=0;i<numAtoms;i++)
			delete atoms[i];
		delete[] atoms;
	}

	atoms = new Atom*[list.size()];
	
	for (int i=0;it!=list.end();it++,i++)
		atoms[i] = (*it);

	numAtoms = (int)list.size();
	cout << endl << INFO << "Read " << numAtoms << " atoms";
	list.clear();

	if (numAtoms<4)
	{
		cout << endl << ERR << "NanoShaper needs at least 4 atoms to work.";
		cout << endl << REMARK << "To emulate 4 atoms you can place dummy atoms with null radius";
		cout << endl << REMARK << "at the same centers of the real atoms";
		exit(-1);
	}

	if (max_rad==0)
	{
		cout << endl << ERR << "All null radii? If you are using place holders atoms please set at least one radius > 0";
		exit(-1);
	}

	return true;

}

bool DelPhiShared::loadAtoms(int na,double* x,double* r,double* q,int* d,char* atinf)
{
	cout << endl << INFO << "Read " << na << " atoms";
	if (atoms!=NULL)
	{
		for (int i=0;i<numAtoms;i++)
			delete atoms[i];
		delete[] atoms;
	}

	if (r==NULL || x==NULL)
	{
		cout << endl << WARN << "Atoms info not available!";
		return false;
	}

	atoms = new Atom*[na];
	numAtoms = na;
	int lastChar = 0;
	for (int i=0;i<na;i++)
	{
	//	printf("\nRadius is %f",r[i]);
		double qq;
		int dd;
		if (q==NULL)
			qq = 0;
		else
			qq=q[i];

		if (d==NULL)
			dd = 0;
		else
			dd=d[i];

		// parse atinf variable and store if available		
		if (atinf!=NULL)
		{
			isAvailableAtomInfo = true;
			char buff[100];
			for (int kk=0;kk<15;kk++)
				buff[kk] = atinf[lastChar++];
			buff[15]='\0';
			//cout << endl << buff;
			string str(buff);
			string name = str.substr(0,5);
			// remove white spaces
			name.erase( remove_if( name.begin(), name.end(), static_cast<int(*)(int)>(isspace) ), name.end() );		
			string resName = str.substr(6,3);
			// remove white spaces
			resName.erase( remove_if( resName.begin(), resName.end(), static_cast<int(*)(int)>(isspace) ), resName.end() );
			atoms[i] = new Atom(x[i*3],x[i*3+1],x[i*3+2],r[i],qq,dd,name,resName,atoi(str.substr(11,4).c_str()),str.substr(10,1));
			//atoms[i]->print();
		}
		else
			atoms[i] = new Atom(x[i*3],x[i*3+1],x[i*3+2],r[i],qq,dd);		
	}

	return true;
}
	
void DelPhiShared::getBounds(double* cmin,double* cmax)
{
	cmin[0] = INFINITY;
	cmin[1] = INFINITY;
	cmin[2] = INFINITY;
	
	cmax[0] = -INFINITY;
	cmax[1] = -INFINITY;
	cmax[2] = -INFINITY;

	for (int ia=0;ia<numAtoms;ia++)
	{
		cmin[0] = min(cmin[0],atoms[ia]->pos[0]-atoms[ia]->radius);		
		cmin[1] = min(cmin[1],atoms[ia]->pos[1]-atoms[ia]->radius);
		cmin[2] = min(cmin[2],atoms[ia]->pos[2]-atoms[ia]->radius);

		cmax[0] = max(cmax[0],atoms[ia]->pos[0]+atoms[ia]->radius);		
		cmax[1] = max(cmax[1],atoms[ia]->pos[1]+atoms[ia]->radius);
		cmax[2] = max(cmax[2],atoms[ia]->pos[2]+atoms[ia]->radius);
	}
}
bool DelPhiShared::buildGrid(double scale,double perfill)
{
	if (atoms==NULL)
	{
		cout << endl << ERR << "Cannot build grid with no atoms!";
		return false;
	}

	if (scale<0)
	{
		cout << endl << WARN << "Cannot use a <0 scale: setting " << DEFAULT_SCALE;
		scale = DEFAULT_SCALE;
	}

	if (perfill < 0 || perfill > 100)
	{
		cout << endl << WARN << "Perfil is in (0,100). Setting " << DEFAULT_PERFIL;
		scale = DEFAULT_PERFIL;	
	}

	this->scale = scale;
	this->perfill = perfill;

	double cmin[3],cmax[3],oldmid[3];

	getBounds(cmin,cmax);

	oldmid[0]=(cmax[0]+cmin[0])/2.;
	oldmid[1]=(cmax[1]+cmin[1])/2.;
	oldmid[2]=(cmax[2]+cmin[2])/2.;

	cout << endl << INFO << "Geometric baricenter ->  " << oldmid[0] << 
		" " << oldmid[1] << " " << oldmid[2];

	double v[6];
	v[0] = fabs(cmax[0]-oldmid[0]);
	v[1] = fabs(cmin[0]-oldmid[0]);
	v[2] = fabs(cmax[1]-oldmid[1]);
	v[3] = fabs(cmin[1]-oldmid[1]);
	v[4] = fabs(cmax[2]-oldmid[2]);
	v[5] = fabs(cmin[2]-oldmid[2]);

	rmaxdim = -INFINITY;

	for (int i=0;i<6;i++)
		if (rmaxdim<v[i])
			rmaxdim = v[i];

	rmaxdim = 2*rmaxdim;

	       
	unsigned int igrid = (unsigned int)floor(scale*100./perfill*rmaxdim);

	if ((igrid%2)==0)
		igrid++;
	
	cout << endl << INFO << "Grid is " << igrid;

	baricenter[0] = oldmid[0];
	baricenter[1] = oldmid[1];
	baricenter[2] = oldmid[2];

	xmin = oldmid[0]-(igrid-1)/(2*scale);
	ymin = oldmid[1]-(igrid-1)/(2*scale);
	zmin = oldmid[2]-(igrid-1)/(2*scale);

	xmax = oldmid[0]+(igrid-1)/(2*scale);
	ymax = oldmid[1]+(igrid-1)/(2*scale);
	zmax = oldmid[2]+(igrid-1)/(2*scale);

	cout << endl << INFO << "MAX " << xmax << " " << ymax << " " << zmax;
	cout << endl << INFO << "MIN " << xmin << " " << ymin << " "  << zmin;
	cout << endl << INFO << "Perfil " << perfill << " %";
	cout << endl << INFO << "Rmaxdim " << rmaxdim;

	cout << endl << INFO << "Allocating memory...";
	cout.flush();

	if (x!=NULL)
		deleteVector<double>(x);
	if (y!=NULL)
		deleteVector<double>(y);
	if (z!=NULL)
		deleteVector<double>(z);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	side = 1/scale;
	hside = side/2;
	A = side*side;

	nx = igrid;
	ny = igrid;
	nz = igrid;

	for (int i=0;i<nx;i++)
		x[i]=xmin+i*side;

	for (int i=0;i<ny;i++)
		y[i]=ymin+i*side;

	for (int i=0;i<nz;i++)
		z[i]=zmin+i*side;

	// allocate epsmap memory	
	bool flag = true;
	
	if (buildEpsMap)
	{
		flag = clearAndAllocEpsMaps();
		
		if (!flag)
			return false;

		if (idebmap!=NULL)
			deleteVector<bool>(idebmap);

		int tot = nx*ny*nz;
		idebmap = allocateVector<bool>(tot);
		for (int i=0;i<tot;i++)
			idebmap[i]=true;
	}

	// allocate grid point status memory
	// and initialize all to temporary outside status
	if (buildStatus)
	{		
		if (status!=NULL)
			deleteVector<short>(status);

		int tot = nx*ny*nz;
		status = allocateVector<short>(tot);
		//int alignementBits = 64;
		//status = allocateAlignedVector<short>(tot,alignementBits);

		if (status==NULL)
		{
			cout << endl << ERR << "Not enough memory to allocate status map";
			cout << endl;
			exit(-1);
		}
	
		for (int i=0;i<tot;i++)
			status[i]=STATUS_POINT_TEMPORARY_OUT;
	}

	cout << "ok!";
	return true;

}



bool DelPhiShared::clearAndAllocEpsMaps()
{
	if (epsmap!=NULL)
		deleteVector<int>(epsmap);

	size_t tot = size_t(nx)*size_t(ny)*size_t(nz)*3;	
	epsmap = allocateVector<int>(tot);
	
	if (epsmap == NULL)
	{
		cout << endl << ERR << "Error allocating epsmap!";
		return false;
	}
	
	// init epsmap
	for (size_t i=0;i<tot;i++)
		epsmap[i]=0;
	
	return true;

}

void DelPhiShared::clearEpsMaps()
{

	if (epsmap == NULL)
	{
		cout << endl << ERR << "Cannot clear null map!";
		return;
	}
	
	size_t n = nx*ny*nz*3;
	for (size_t i= 0;i<n;i++)
		epsmap[i] = 0;
	/*
	for (int i=0;i<nx;i++)
			for (int j=0;j<ny;j++)
				for (int k=0;k<nz;k++)
				{
					EPSMAP(i,j,k,0,nx,ny,nz)=0;
					EPSMAP(i,j,k,1,nx,ny,nz)=0;
					EPSMAP(i,j,k,2,nx,ny,nz)=0;
				}
	*/
	return;
}

void DelPhiShared::saveIdebMap(char* fname)
{
	if  (idebmap == NULL)
	{
		cout << endl << WARN << "Cannot save null idebmap!";
		return;
	}
	
	char f1[BUFLEN];	
	//if (debug_stern == -1)
		sprintf(f1,"%s.idebmap.txt",fname);
	//else
	//	sprintf(f1,"%s.idebmap.stern.txt",fname);
	FILE* fp = fopen(f1,"w");
	for (int i=0;i<nx;i++)
	{		
		for (int j=0;j<ny;j++)
		{
			for (int k=0;k<nz;k++)
				//fprintf(fp,"%d ",(int)IDEBMAP(i,j,k,nx,ny));
				fprintf(fp,"%d ",read3DVector<bool>(idebmap,i,j,k,nx,ny,nz));

			fprintf(fp,"\n");
		}		
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void DelPhiShared::saveEpsMaps(char* fname)
{
	if  (epsmap == NULL)
	{
		cout << endl << WARN << "Cannot save null epsmap!";
		return;
	}
	
	char f1[BUFLEN],f2[BUFLEN],f3[BUFLEN];
	sprintf(f1,"%s.epsmapx.txt",fname);
	sprintf(f2,"%s.epsmapy.txt",fname);
	sprintf(f3,"%s.epsmapz.txt",fname);

	FILE* fp = fopen(f1,"w");
	for (int i=0;i<nx;i++)
	{		
		for (int j=0;j<ny;j++)
		{
			for (int k=0;k<nz;k++)
				//fprintf(fp,"%d ",EPSMAP(i,j,k,0,nx,ny,nz));
				fprintf(fp,"%d ",read4DVector<int>(epsmap,i,j,k,0,nx,ny,nz,3));

			fprintf(fp,"\n");
		}		
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen(f2,"w");
	for (int i=0;i<nx;i++)
	{		
		for (int j=0;j<ny;j++)
		{
			for (int k=0;k<nz;k++)
				//fprintf(fp,"%d ",EPSMAP(i,j,k,1,nx,ny,nz));
				fprintf(fp,"%d ",read4DVector<int>(epsmap,i,j,k,1,nx,ny,nz,3));

			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}		
	fclose(fp);

	fp = fopen(f3,"w");
	for (int i=0;i<nx;i++)
	{		
		for (int j=0;j<ny;j++)
		{
			for (int k=0;k<nz;k++)
				//fprintf(fp,"%d ",EPSMAP(i,j,k,2,nx,ny,nz));
				fprintf(fp,"%d ",read4DVector<int>(epsmap,i,j,k,2,nx,ny,nz,3));

			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}		
	fclose(fp);
}

DelPhiShared::~DelPhiShared()
{
	clear();
}

void DelPhiShared::clear()
{
	if (epsmap!=NULL)
		deleteVector<int>(epsmap);
	if  (ibgp!=NULL)
		deleteVector<int>(ibgp);
	if  (scspos!=NULL)
		deleteVector<double>(scspos);
	if  (scsnor!=NULL)
		deleteVector<double>(scsnor);
	if  (scsarea!=NULL)
		deleteVector<double>(scsarea);
	if (x!=NULL)
		deleteVector<double>(x);
	if (y!=NULL)
		deleteVector<double>(y);
	if (z!=NULL)
		deleteVector<double>(z);	
	if (status!=NULL)
		//deleteMatrix3D<short>(nx,ny,status);
		deleteVector<short>(status);
		//deleteAlignedVector<short>(status);
	if (idebmap!=NULL)
		deleteVector<bool>(idebmap);

	if (atoms!=NULL)
	{
		for (int i=0;i<numAtoms;i++)
			delete atoms[i];
		delete[] atoms;
	}

	// free memory 
	if (cavitiesVec !=NULL)
	{
		vector< vector < int* >* >::iterator it;
		for (it=cavitiesVec->begin();it!=cavitiesVec->end();it++)
		{
			vector < int* >* inner = (*it);
			vector < int* >::iterator it2;
			for (it2=inner->begin();it2!=inner->end();it2++)
				free((*it2));
			delete inner;
		}		
	}
	delete cavitiesVec;
}

void DelPhiShared::saveBGP(char* fname)
{
	if  (scspos == NULL || scsnor==NULL)
	{
		cout << endl << WARN << "Cannot save null scspos or scsnor!";
		return;
	}

	char ff[BUFLEN];
	sprintf(ff,"%s.projections.txt",fname);		
	FILE* fp = fopen(ff,"w");
	FILE* fp2 = fopen("surf.txt","w");
	fprintf(fp,"%d\n",nbgp);
	for (int iibgp=0;iibgp<nbgp;iibgp++)
	{
		fprintf(fp,"%d %d %d %lf %lf %lf %lf %lf %lf",ibgp[iibgp*3]+1,ibgp[iibgp*3+1]+1,ibgp[iibgp*3+2]+1,
			scspos[iibgp*3],scspos[iibgp*3+1],scspos[iibgp*3+2],
			scsnor[iibgp*3],scsnor[iibgp*3+1],scsnor[iibgp*3+2]);

		fprintf(fp2,"%lf %lf %lf %lf %lf %lf",scspos[iibgp*3],scspos[iibgp*3+1],scspos[iibgp*3+2],
			scsnor[iibgp*3],scsnor[iibgp*3+1],scsnor[iibgp*3+2]);

		if (iibgp<nbgp-1)
		{
			fprintf(fp,"\n");
			fprintf(fp2,"\n");
		}
	}		
	fclose(fp);
	fclose(fp2);
}

void DelPhiShared::saveStatus(char* fname)
{
	if  (status == NULL)
	{
		cout << endl << WARN << "Cannot save null status!";
		return;
	}
	char ff[BUFLEN];
	sprintf(ff,"%s.status.txt",fname);
	FILE* fp = fopen(ff,"w");
	
	for (int i=0;i<nx;i++)
	{		
		for (int j=0;j<ny;j++)
		{
			for (int k=0;k<nz;k++)
			{
				//fprintf(fp,"%d ",status[i][j][k]);
				//fprintf(fp,"%d ",STATUSMAP(i,j,k,nx,ny));
				fprintf(fp,"%d ",read3DVector<short>(status,i,j,k,nx,ny,nz));
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	saveCavities(false);
}

void DelPhiShared::initCav2Atoms()
{
	int nc = (int)cavitiesVec->size();

	if (cav2atoms.size()!=0)
		clearCav2Atoms();

	cav2atoms.reserve(nc);
	for (int ic=0;ic<nc;ic++)
		cav2atoms.push_back(new set<int>());
}

void DelPhiShared::clearCav2Atoms()
{
	if (cav2atoms.size()!=0)
	{
		for (unsigned int ic=0;ic<cav2atoms.size();ic++)
			delete cav2atoms[ic];
	}
	cav2atoms.clear();
}

void DelPhiShared::markPockets(short* stat,vector<bool>& isPocket)
{
	vector < vector <int* >* >::iterator it;
	int num = (int)cavitiesVec->size();
	
	if (stat==NULL)
	{
		cout << endl << WARN << "Cannot mark pockets with a reference status map!";
		return;
	}

	if (stat==status)
	{
		cout << endl << WARN << "In marking pockets the status map passed must be from another reference object; here the same was passed";
		return;
	}

	short* status = stat;

	if (num==0)
		return;

	int i=0,sync_i=0;
	
	for (i=0,sync_i=0,it=cavitiesVec->begin();it!=cavitiesVec->end();it++,i++,sync_i++)
	{
		if (cavitiesFlag[sync_i]==true)
		{
			i--;
			continue;
		}
		
		vector <int* >::iterator it2;
		vector<int*>* vec = (*it);

		bool isPocketFlag = false;
		
		for (it2=vec->begin();it2!=vec->end();it2++)
		{
			int* vec = *it2;
			short value = read3DVector<short>(status,vec[0],vec[1],vec[2],nx,ny,nz);
			// if at least one was out, then that's a pocket
			//if (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny)==STATUS_POINT_TEMPORARY_OUT || STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny)==STATUS_POINT_OUT)
			if (value==STATUS_POINT_TEMPORARY_OUT || value==STATUS_POINT_OUT)
			{
				isPocketFlag = true;
				break;
			}
		}
		isPocket.push_back(isPocketFlag);
	}	
	return;
}

void DelPhiShared::parseAtomInfo(char* atinf,string mol,int natom,double* xn1,double* rad)
{
	// get atoms and info from delphi
	FILE* ff = fopen(mol.c_str(),"w");
	int lastChar = 0;
		
	// parse atinfo into a file
	for (int i=0;i<natom;i++)
	{
		char buff[100];
		for (int kk=0;kk<15;kk++)
			buff[kk] = atinf[lastChar++];

		buff[15]='\0';
		string str(buff);
		string name = str.substr(0,5);
		// avoid disambiguation
		name.erase( remove_if( name.begin(), name.end(), static_cast<int(*)(int)>(isspace) ), name.end() );			
		string resName = str.substr(6,3);			
		resName.erase( remove_if( resName.begin(), resName.end(), static_cast<int(*)(int)>(isspace) ), resName.end() );	
		fprintf(ff,"%f %f %f %f 0 %s %s %d %s\n",xn1[3*i],xn1[3*i+1],xn1[3*i+2],rad[i],name.c_str(),resName.c_str(),atoi(str.substr(11,4).c_str()),str.substr(10,1).c_str());
	}
	fclose(ff);
}

int DelPhiShared::cavitiesToAtoms(double rad)
{
	vector < vector <int* >* >::iterator it;
	int num = (int)cavitiesVec->size();

	if (num==0)
		return 0;

	char buff[100];
	int i=0,sync_i=0;
	int nonActive = 0;
	FILE* fp;

	for (i=0,sync_i=0,it=cavitiesVec->begin();it!=cavitiesVec->end();it++,i++,sync_i++)
	{
		if (cavitiesFlag[sync_i]==true)
		{
			nonActive++;
			i--;
			continue;
		}
		sprintf(buff,"cav%d.txt",i);
		fp = fopen(buff,"w");
		sprintf(buff,"all_cav%d.txt",i);
		FILE* fp2 = fopen(buff,"w");
		
		vector <int* >::iterator it2;
		vector<int*>* vec = (*it);
		int cavityId = sync_i+STATUS_FIRST_CAV;

		int count = 0;
		double lastx,lasty,lastz;
		for (it2=vec->begin();it2!=vec->end();it2++)
		{
			int* vec = *it2;
			fprintf(fp2,"%f %f %f %f\n",x[vec[0]],y[vec[1]],z[vec[2]],rad);
			
			// save only support cavity points
			//if (STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny)==-cavityId)
			if (read3DVector<short>(status,vec[0],vec[1],vec[2],nx,ny,nz)==-cavityId)
			{
				//printf("\n %d",STATUSMAP((vec[0]),(vec[1]),(vec[2]),nx,ny));
				fprintf(fp,"%f %f %f %f\n",x[vec[0]],y[vec[1]],z[vec[2]],rad);
				//printf("\n%d %d %d %f %f %f",vec[0],vec[1],vec[2],x[vec[0]],y[vec[1]],z[vec[2]]);
				lastx = x[vec[0]];
				lasty = y[vec[1]];
				lastz = z[vec[2]];
				count++;
			}
		}

		if (count == 0)
		{
			cout << endl << ERR << "Zero support atoms to save";
			exit(-1);
		}
		
		// clone atoms to let NanoShaper work on this set of points
		// NanoShaper needs at least 4 atoms for triangulation
		if (count<4)
		{
			count = 4 - count;
			for (int i=0;i<count;i++)
				fprintf(fp,"%f %f %f %f\n",lastx,lasty,lastz,rad);
		}
		fclose(fp);
		fclose(fp2);
	}	
	fp = fopen("numcav.txt","w");
	fprintf(fp,"%d",num-nonActive);
	fclose(fp);

	return num-nonActive;
}

void DelPhiShared::saveCavities2(bool onlySaveNonFilled,string sysName)
{
	if (isAvailableAtomInfo)
	{
		FILE* fp;
		char buff[100];
		sprintf(buff,"%s.pocket",sysName.c_str());
		fp = fopen(buff,"w");
		int savedIndex = 1;
		for (unsigned int i=0;i<cav2atoms.size();i++)
		{		
			set<int>::iterator setIt;		
			set<int>* setPt = cav2atoms[i];

			// save only non filled if required, or save all of them
			if ((onlySaveNonFilled && !cavitiesFlag[i]) || (!onlySaveNonFilled))
			{	
				for (setIt = setPt->begin(); setIt!=setPt->end();setIt++)
				{
					//printf("\nIndex %d",i);
					fflush(stdout);
					int ii = (*setIt);
					Atom* at = atoms[ii];
					AtomInfo& ai = at->ai;
					if ((ai.getName()).size()==4)
						fprintf(fp,"ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f                 %6d POC\n",ii+1,ai.getName().c_str(),ai.getResName().c_str(),ai.getChain().c_str(),ai.getResNum(),at->pos[0],at->pos[1],at->pos[2],savedIndex);
					else
						fprintf(fp,"ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f                 %6d POC\n",ii+1,ai.getName().c_str(),ai.getResName().c_str(),ai.getChain().c_str(),ai.getResNum(),at->pos[0],at->pos[1],at->pos[2],savedIndex);
				}	
				savedIndex++;
			}
		}
		fclose(fp);		
	}
}
void DelPhiShared::saveCavities(bool onlySaveNonFilled)
{
	if  (status == NULL)
	{
		cout << endl << WARN << "Cannot save cavities with null status!";
		return;
	}

	FILE* fp2 = fopen("cavities.txt","w");
	for (int k=0;k<nz;k++)	
	{		
		for (int j=0;j<ny;j++)
		{
			for (int i=0;i<nx;i++)
			{
				//if (status[i][j][k]>3)
				//if ((STATUSMAP(i,j,k,nx,ny))>STATUS_POINT_OUT)
				if (read3DVector<short>(status,i,j,k,nx,ny,nz)>STATUS_POINT_OUT)
					fprintf(fp2,"\n %f %f %f",x[i],y[j],z[k]);
			}
		}		
	}	
	fclose(fp2);

	fp2 = fopen("cavitiesSize.txt","w");	
	vector <double>::iterator itCav; 	
	vector <int*>::iterator itCavVec;
	int i=0;

	fprintf(fp2,"%f\n",((float)cavitiesSize.size()));

	// for each cavity save the baricenter and the estimated volume
	// save all cavities; add a flag for the filled one
	for (itCav =  cavitiesSize.begin(),i=0; itCav!=cavitiesSize.end();itCav++,i++)
	{
		double barx,bary,barz;
		barx = 0;
		bary = 0;
		barz = 0;
		for (itCavVec =  cavitiesVec->at(i)->begin(); itCavVec!=cavitiesVec->at(i)->end();itCavVec++)
		{
			barx += x[(*itCavVec)[0]];
			bary += y[(*itCavVec)[1]];
			barz += z[(*itCavVec)[2]];
		}
		double volSize = cavitiesSize.at(i);
		int size = (int)cavitiesVec->at(i)->size();
		barx/=size;
		bary/=size;
		barz/=size;
		int fillFlag;

		if (cavitiesFlag[i])
			fillFlag = 1;
		else
			fillFlag = 0;
		// save baricenter,volume, and filling flag
		fprintf(fp2,"%lf %lf %lf %lf %d\n",barx,bary,barz,volSize,fillFlag);	
	}
	fclose(fp2);		 

	fp2 = fopen("cavAtomsSerials.txt","w");
	for (unsigned int i=0;i<cav2atoms.size();i++)
	{		
		set<int>::iterator setIt;		
		set<int>* setPt = cav2atoms[i];
		
		// save only non filled if required, or save all of them
		if ((onlySaveNonFilled && !cavitiesFlag[i]) || (!onlySaveNonFilled))
		{			
			for (setIt = setPt->begin(); setIt!=setPt->end();setIt++)
				fprintf(fp2,"%d ",(*setIt)+1);
			fprintf(fp2,"\n");
		}
	}
	fclose(fp2);

	if (isAvailableAtomInfo)
	{
		FILE* fp3;
		//cout << endl << INFO << "Atoms info is available I am saving also residues";
		fp2 = fopen("residues.txt","w");
		fp3 = fopen("residues_vmd.txt","w");

		for (unsigned int i=0;i<cav2atoms.size();i++)
		{		
			set<int>::iterator setIt;		
			set<int>* setPt = cav2atoms[i];

			// only save the non filled cavities
			if (onlySaveNonFilled && !cavitiesFlag[i])
			{	
				// detected the involved residues in the cavity
				set<int> involvedResidues;
				for (setIt = setPt->begin(); setIt!=setPt->end();setIt++)
				{
					int atomIndex = (*setIt);
					involvedResidues.insert(atoms[atomIndex]->ai.getResNum());
				}
				
				int ll=0,tot=(int)involvedResidues.size();

				for (setIt = involvedResidues.begin(); setIt!=involvedResidues.end();setIt++,ll++)
				{
					int resIndex = (*setIt);
					fprintf(fp2,"%d ",resIndex);
					fprintf(fp3,"(resid %d and not water) ",resIndex);

					if (ll!=(tot-1))
						fprintf(fp3," or ");
				}
				fprintf(fp2,"\n");
				fprintf(fp3,"\n");
			}
		}
		fclose(fp2);
		fclose(fp3);
	}
}

#pragma once
#include"ImageF.h"
#include "vol_math_Raw3D_Independt.h"
#include<string>
using  std::string;
class LevelSet
{
public:
	ImageF regFunction(ImageF &s,int m,int n);
	ImageF& initialg(ImageF &g);
	//LevelSet(void);
	~LevelSet(void);
	ImageF& drlse_edge(ImageF &phi_0,ImageF &g,float lambda,float mu,float alfa,float epsilon,int timestep, int iter,const char * sada);
	void testout(Raw2D *ret);
};


#pragma once
#include"ImageF.h"
#include<string>
using  std::string;
class LevelSet
{
public:
	ImageF regFunction(ImageF &s,int m,int n);
	LevelSet(void);
	~LevelSet(void);
	ImageF drlse_edge(ImageF phi_0,ImageF g,float lambda,float mu,float alfa,float epsilon,int timestep, int iter,string sada);



};


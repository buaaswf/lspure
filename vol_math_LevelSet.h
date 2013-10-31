#ifndef _VOL_MATH_LEVELSET_H_
#define _VOL_MATH_LEVELSET_H_

#include "vol_math_Raw3D_Independt.h"
#include "test.h"

#include<string>
using  std::string;
class LevelSet
{
public:
	void initialg(Raw2D *g);
	//LevelSet(void);
	~LevelSet(void);
	void drlse_edge(Raw2D *phi,Raw2D *g,float lambda,float mu,float alfa,float epsilon,int timestep, int iter, const char * sada);
	Raw2D& uchar2double(Raw2D &img);
	void testout(Raw2D *ret);
};

Raw2D regFunction(Raw2D &s,double m,double n);

#endif  //_VOL_MATH_LEVELSET_H_
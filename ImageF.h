#pragma once
#include"Raw3D_Independt.h"
class ImageF :public Raw2D
{	
public:
	int length;
	 int width;
public:
	ImageF(void);
	ImageF(int,int,PIXTYPE *);
	~ImageF(void);
	friend ImageF operator/(ImageF x,float s);
};


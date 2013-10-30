#ifndef _TEST_H_
#define _TEST_H_

#include<iostream>
#include <crtdbg.h> 
#include "vol_math_Raw3D_Independt.h"
#include "CImg.h" 
using namespace cimg_library;

static void IShowImg(Raw2D& img)
{
	CImg <double> sourceimage(img.getXsize(),img.getYsize(),1,1,0);
	cimg_for_insideXY(sourceimage,x,y,0)
	{
		PIXTYPE val=img.get(x,y);
		/*if (val>0)
		{*/
			sourceimage(x,y,0)=(double)(val);
		//}
		//
		
	}
	sourceimage.display("hello");
}
#endif
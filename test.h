#ifndef _TEST_H_
#define _TEST_H_

#include "cv.h"
#include "highgui.h"
#include<iostream>
#include <crtdbg.h> 
#include "vol_math_Raw3D_Independt.h"
#include "CImg.h" 
using namespace cimg_library;
static void showImg(Raw2D& img)
{
	IplImage *source = cvCreateImage(cvSize(img.getYsize(), img.getXsize()), 8, 1);
	for (int i=0;i<source->height;i++)//
	{
		for (int j=0;j<source->width;j++)//
		{
				CV_IMAGE_ELEM(source,uchar,i,j) = img.get(i,j);
		}
	}
	if (source!=NULL)
	{
		cvNamedWindow("source",0);
		cvShowImage("source",source);
		cvWaitKey(0);
		cvDestroyWindow("source");
	}
	cvReleaseImage(&source);
}
static void IShowImg(Raw2D& img)
{
	CImg <double> sourceimage(img.getYsize(),img.getXsize(),1,1,0);
	cimg_for_insideXY(sourceimage,x,y,0)
	{
		PIXTYPE val=img.get(y,x);
		/*if (val>0)
		{*/
			sourceimage(x,y,0)=(double)(val);
	/*	}*/
		
	}
	sourceimage.display("hello");
}
#endif
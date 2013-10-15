#include"ImageF.h"
#include "vol_math_LevelSet.h"
#include "cv.h"
#include "highgui.h"
#include<iostream>
#include <crtdbg.h> 
using namespace std;
int main()
{
	int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
	tmpFlag |= _CRTDBG_LEAK_CHECK_DF; 
	_CrtSetDbgFlag( tmpFlag ); 
	IplImage *source=NULL;
	source=cvLoadImage("E:\\geo\\levelset\\DRLSE_v0\\gourd.bmp",0);
	//source=cvLoadImage("E:\\geo\\testbmps\\ls1.bmp",0);
	if (source!=NULL)
	{
		cvNamedWindow("source",0);
		cvShowImage("source",source);
		cvWaitKey(0);
		cvDestroyWindow("source");
	}
	int x=source->height;
	int y=source->width;
	PIXTYPE *p=new PIXTYPE[x*y];
	for (int i=0;i<x;i++)
	{
		for (int j=0;j<y;j++)
		{
			p[i*y+j]=CV_IMAGE_ELEM(source,uchar,i,j);
		}
	}

	Raw2D *initial=new Raw2D(x,y);
	Raw2D *raw2d=new Raw2D(x,y,p);
	//cvReleaseImage(&source);

	for (int i=0;i<x;i++)
	{
		for (int j=0;j<y;j++)
		{
			if(i<25&&i>20&&j>20&&j<25)
				initial->putXY(i*y+j,255);
			else initial->putXY(i*y+j,0);
			//cout<<initial->get(i,j)<<endl;
			
		}
		cout<<endl;
	}

	//raw2d->guassConv(raw2d,2);
	LevelSet *ls=new LevelSet();
	char const *pt="single-well";
	Raw2D *ret=new Raw2D(ls->drlse_edge(*initial,*raw2d,1.0,1.0,1.0,3.0,1,1,pt));
	PIXTYPE *result;
	//result=ret->gety();
	
	for (int i=0;i<source->height;i++)
	{
		for (int j=0;j<source->width;j++)
		{
			CV_IMAGE_ELEM(source,uchar,i,j)=ret->get(i,j);
		}
	}
	if (source!=NULL)
	{
		cvNamedWindow("source",0);
		cvShowImage("source",source);
		cvWaitKey(0);
		cvDestroyWindow("source");
	}

	 system("pause");
}
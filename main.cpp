#include"ImageF.h"
#include "vol_math_LevelSet.h"
#include"levelset.h"
#include "test.h"
//#include "cv.h"
//#include "highgui.h"
#include<iostream>
#include <crtdbg.h> 
#include "CImg.h" 
using namespace cimg_library;
using namespace std;

int main()
{
	int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
	tmpFlag |= _CRTDBG_LEAK_CHECK_DF; 
	_CrtSetDbgFlag( tmpFlag ); 
	//IplImage *source=NULL;
	//source=cvLoadImage("G:\\geo hackthon\\geo\\levelset\\DRLSE_v0\\gourd.bmp",0);
	//source=cvLoadImage("E:\\geo\\levelset\\DRLSE_v0\\gourddouble.bmp",0);
	CImg<unsigned char> source("E:\\geo\\levelset\\DRLSE_v0\\gourd.bmp");
	//int x=source->height;
	//int y=source->width;
	int x=source.height();
	int y=source.width();
	PIXTYPE *p=new PIXTYPE[x*y];
	PIXTYPE *q=new PIXTYPE[x*y];
	for (int i=0;i<x;i++)
	{
		for (int j=0;j<y;j++)
		{
			//p[i*y+j]=CV_IMAGE_ELEM(source,uchar,i,j);
			p[i*y+j]=source.atXY(i,j)/255.0;
		}
	}

	Raw2D *initial=new Raw2D(x,y,q);
	Raw2D *raw2d=new Raw2D(x,y,p);
	//cvReleaseImage(&source);

	for (int i=0;i<x;i++)
	{
		for (int j=0;j<y;j++)
		{
			if(i>25&&i<35&&j>40&&j<50)
				initial->putXY(i*y+j,1.0);
			else if (i<25||i>35||j<40||j>50)
			{
				initial->putXY(i*y+j,-1.0);
			}

			else	initial->putXY(i*y+j,0);
			//cout<<initial->get(i,j)<<endl;
			
		}
		//cout<<endl;
	}
	//showImg(*initial);
	IShowImg(*initial);
	//showImg(*raw2d);
	IShowImg(*raw2d);
	//raw2d->guassConv(raw2d,2);
	//Raw2D *data=new Raw2D(raw2d);
	LevelSet *ls=new LevelSet();
	ls->initialg(*raw2d);
	
	char const *pt="single-well";
	int iter_outer=1;
	
	ls->drlse_edge(*initial,*raw2d,5.0,0.2,-3,255,1,iter_outer,pt);

	 system("pause");
}
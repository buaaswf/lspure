#include "vol_math_LevelSet.h"
#include "test.h"
#include <iostream>
#include <crtdbg.h> 
#include "CImg.h" 
#include "Filter.h"
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
	CImg<unsigned char> source;//("E:\\geo\\levelset\\DRLSE_v0\\gourd.bmp");
	source.load("gourd.bmp");
	/*source.load("E:\\geo\\levelset\\DRLSE_v0\\gourd.bmp");*/
	//source.display();
	//int x=source->height;
	//int y=source->width;
	int x=source.width();
	int y=source.height();
	PIXTYPE *p=new PIXTYPE[x*y];
	PIXTYPE *q=new PIXTYPE[x*y];
	for (int i=0;i<x;i++)
	{
		for (int j=0;j<y;j++)
		{
			//p[i*y+j]=CV_IMAGE_ELEM(source,uchar,i,j);
			p[i*y+j]=(double)source.atXY(i,j);
		}
	}

	Raw2D *initial=new Raw2D(x,y,q);
	Raw2D *raw2d=new Raw2D(x,y,p);
	//int inx=40;
	//int inxl=50;
	//int iny=25;
	//int inyl=35;
	int inx=10;//ÉÏÏÂ
	int inxl=70;
	int iny=10;//×óÓÒ
	int inyl=50;
	//cvReleaseImage(&source);
/************************************************************************/
/*  initial contour out of  region                                                                    */
/************************************************************************/

	//for (int i=0;i<x;i++)
	//{
	//	for (int j=0;j<y;j++)
	//	{
	//		if(i>10&&i<60&&j>9&&j<50)
	//			initial->putXY(i*y+j,2.0);
	//		else if (i<10||i>60||j<9||j>50)
	//		{
	//			initial->putXY(i*y+j,-2.0);
	//		}

	//		else	initial->putXY(i*y+j,-2.0);
	//		//cout<<initial->get(i,j)<<endl;
	//		
	//	}
	//	//cout<<endl;
	//}
	/************************************************************************/
	/* initial inside the region                                                                     */
	/************************************************************************/
	for (int i=0;i<x;i++)
	{
		for (int j=0;j<y;j++)
		{
			if(i>inx&&i<inxl&&j>iny&&j<inyl)
				initial->putXY(i*y+j,-2.0);
			else if (i<inx||i>inxl||j<iny||j>inyl)
			{
				initial->putXY(i*y+j,2.0);
			}

			else	initial->putXY(i*y+j,-2.0);
			//cout<<initial->get(i,j)<<endl;

		}
		//cout<<endl;
	}
	//showImg(*initial);
	IShowImg(*initial);
	//showImg(*raw2d);
	IShowImg(*raw2d);
	Filter filter;
	filter.guassFilter(raw2d,2);
	IShowImg(*raw2d);
	//Raw2D *data=new Raw2D(raw2d);
	LevelSet *ls=new LevelSet();
	ls->initialg(raw2d);
	//enum PotentialFunction{single_well = 1,double_well};
	char const *pt="single_well";
	int iter_outer=120;
	
	ls->drlse_edge(initial,raw2d,5.0,0.2,3,1.5,1,iter_outer,pt);

	 system("pause");
}
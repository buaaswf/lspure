//#include"Canny.h"
#include"Levelset.h"
#include"ImageF.h"
#include"LevelSet.h"
#include "Raw3D_Independt.h"
#include<iostream>;
using namespace std;
main()
{
 	PIXTYPE *p=new PIXTYPE[100];
	Raw2D *raw2d=new Raw2D(10,10,p);
	//Raw2d *raw2d =new Raw2d(10,10,p);	
	PIXTYPE *y=NULL;
	int i=0,j=0;

	for (i=0;i<10;i++)
	{
		for (j=0;j<10;j++)
		{
			raw2d->putXY(i*10+j,i+1);
			y=raw2d->gety();
			PIXTYPE test=p[i*10+j];

		}
		
	}
	//Raw2D *raw2dderive=static_cast<ImageF> *raw2d;
	//raw2d->guassConv(raw2d,2);
	LevelSet *ls=new LevelSet();
	Raw2D *raw2dderive=raw2d->guassConv(raw2d,2);
	ls->drlse_edge(*raw2d,*raw2dderive,1,1,1,1,1,2,"single-well");
	p=raw2d->gety();
	 i=0;
/*
	while (p[i])
	{
		PIXTYPE test = p[i];
		i++;
		cout<<"no:"<<i<<"p[i]="<<p[i]<<endl;
	}*/

	 system("pause");
}
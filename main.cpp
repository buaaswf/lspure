//#include"Canny.h"
#include"Levelset.h"
#include<iostream>;
using namespace std;
main()
{
 	PIXTYPE *p=new PIXTYPE[100];
// 	float g[5][5]={ {0.0000 ,   0.0000 ,   0.0002  ,  0.0000  ,  0.0000},
// 		{0.0000 ,   0.0113 ,   0.0837   , 0.0113  ,0.0000},
// 	{0.0002  ,  0.0837  ,  0.6187   , 0.0837 ,   0.0002},
// 	{0.0000 ,   0.0113 ,   0.0837 ,   0.0113  ,  0.0000},
// 	{0.0000 ,   0.0000  ,  0.0002 ,   0.0000 ,   0.0000}};

	ImageF *raw2d=new ImageF(10,10,p);	
	PIXTYPE *temp=new PIXTYPE[100];
	PIXTYPE *y=NULL;
	int i=0,j=0;

	for (i=0;i<10;i++)
	{
		for (j=0;j<10;j++)
		{
			raw2d->putXY(i*10+j,i+1);
			PIXTYPE test = p[i*10+j];
			y=raw2d->gety();
			test=y[i*10+j];
		/*	*(p+i*10+j)=j;*/
			
/*
			temp=raw2d->gety();
			int t=temp[i*10+j];
			cout<<t<<endl;*/

		}
		
	}
	raw2d->guassConv(raw2d,2);
	/*
	for ()
		{
		}*/
	memcpy(temp,p,100);
	//p=raw2d->gety();
/*	cout<<*p<<endl;*/
/*
	int num=sizeof(*p);
	for (i=0;i<num;i++)
	{
		cout<<p[i]<<endl;

	}*/
	 i=0;
	while (p[i])
	{
		PIXTYPE test = p[i];
		i++;
		cout<<"no:"<<i<<"p[i]="<<temp[i]<<endl;
	}
	 system("pause");
}
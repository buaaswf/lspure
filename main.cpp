//#include"Canny.h"
#include"Levelset.h"
#include<iostream>;
using namespace std;
int main()
{
 	PIXTYPE *p=new PIXTYPE[100];
	ImageF *raw2d=new ImageF(10,10,p);	
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
	raw2d->guassConv(raw2d,2);
	LevelSet *ls=new LevelSet();
	ls->drlse_edge(*raw2d,*raw2d,1,1,1,1,1,2,"single-well");
	p=raw2d->gety();
	 i=0;
	while (p[i])
	{
		PIXTYPE test = p[i];
		i++;
		cout<<"no:"<<i<<"p[i]="<<p[i]<<endl;
	}
	 system("pause");
}
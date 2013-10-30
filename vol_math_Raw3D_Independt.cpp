#include "vol_math_Raw3D_Independt.h"
#include<math.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include "Filter.h"
using namespace std;
Raw3D_Independt::Raw3D_Independt(void)
{
}


Raw3D_Independt::~Raw3D_Independt(void)
{
}

static float lgtt=log10(2.0f);
#define M_EXP 2.7182818284590452353602874713527

Raw2D::Raw2D()
{
	xsize=0;
	ysize=0;
	data=NULL;
}

Raw2D::Raw2D(int xsize,int ysize,PIXTYPE *y)
{
	this->xsize=xsize;
	this->ysize=ysize;
	this->data=y;
}

Raw2D::Raw2D( int xsize,int ysize )
{
	this->xsize=xsize;
	this->ysize=ysize;
	this->data=new PIXTYPE[xsize*ysize];
}

Raw2D::Raw2D( Raw2D *r)
{
	this->xsize=r->xsize;
	this->ysize=r->ysize;
	this->data=new PIXTYPE[xsize*ysize];
	//for (int i=0;i<xsize;i++)
	//{
	// for (int j=0;j<ysize;j++)
	// {
	//	 data[i*ysize+j]=r->data[i*(r->ysize)+j];
	// }
	//}
	memcpy(this->data,r->data,sizeof(PIXTYPE)*xsize*ysize);
}

Raw2D::Raw2D(const Raw2D& r)
{
	//if (this->xsize != r.getXsize() && this->ysize != r.getYsize())
	//{
		this->xsize=r.xsize;
		this->ysize=r.ysize;
	//	//if (this->data != NULL)
	//	//	delete[] this->data;
		this->data=new PIXTYPE[xsize*ysize];
	//}
	memcpy(this->data, r.data, sizeof(PIXTYPE)*xsize*ysize);
}

Raw2D::~Raw2D(void)
{
	if(this->data!=NULL)
		delete [] this->data;
	data=NULL;
}

/* Raw3D::Raw3D(void)
{
z=new Raw2D();
zsize=0;
}*/

Raw3D::Raw3D(int zsize,Raw2D *src)
{
	this->zsize=zsize;
	this->z=src;
}
Raw3D::Raw3D(void)
{
	z=0;
}

Raw3D::~Raw3D(void){
	if(this->z!=NULL)
		delete [] this->z;
	z=NULL;
	//cout<<"Raw3D is deconstruct"<<endl;
}

//added FOR TEST DATA
void rawarray(int xsize,int ysize,int const zsize,PIXTYPE *yy)
{
	int i=0,j=0;
	//Raw2D *raw2D[10000];//=new Raw2D [zsize];
	//Raw2D Raw2D[zsize];
	//PIXTYPE *p=new PIXTYPE[135161];
	try
	{

		for(i=180;i<zsize;i++)
		{	//p+=zsize;
			PIXTYPE *p=new PIXTYPE[xsize*ysize];
			for(j=0;j<xsize*ysize;j++)
			{
				//for(j=0;j<ysize;j++)
				//{*(p+j*xsize+j)=*(y+j*xsize+j);}
				//if (j==3991)
				//	printf("3991");
				PIXTYPE pp=*(yy+i*(xsize*ysize)+j);
				p[j]=pp;
			}//for j=...

			try
			{
				Raw2D *raw2D=new Raw2D(xsize,ysize,p);
				//Filter *filter=new Filter();
				//filter->TrilateralFilter(raw2D,1);
				delete raw2D;
			}
			catch(std::bad_alloc)
			{
				delete [] p;
				p=NULL;
				i--;
			}

			//raw2D[i]->TrilateralFilter(raw2D[i],3);



			//delete raw2D;
			printf("no  %d\n" ,i);
			//raw2D++;
			//delete [] p;
			//delete raw2D;
			//delete p;
			//p=NULL;
		}//for i=..


	}//try ...
	catch (const std::bad_alloc& e)
	{
		printf("p alloc failed");
	}



	//return raw2D[zsize];
};


//================================================================================
//Some helper functions for Raw2D and Raw3D class.
//================================================================================

/**
Memory allocator for a rectangular 2D object. Takes an 'empty' Raw2D
object and reserves a rectangular field of ixsize by iysize pixels.
**/
void Raw2D::sizer(int ixsize, int iysize) {
	if(data!=NULL)
		delete [] this->data;

	data=NULL;
	this->data = new PIXTYPE[size()];	// & allocate memory.
	xsize = ixsize;				// set new image size,
	ysize = iysize;
}

/**
Memory allocator for a 2D object.  Allocates space for an object of
same size as 'src'.
**/
void Raw2D::sizer(Raw2D* src) 
{
	int ix, iy;
	ix = src->getXsize();
	iy = src->getYsize();
	sizer(ix,iy);
}

/**
Copy all the pixels from 'src'.  If size of 'src' is different from
us, 'wipe' out our current pixels, and change our size to match 'src'.
If a 'wipe' was necessary, return 'FALSE', else return 'TRUE'.
-**/

bool Raw2D::wipecopy(Raw2D* src) {
	bool out = true;
	if(getYsize() != src->getYsize() || getXsize()!=src->getXsize()) { // resize to fit 'src'.
		out = false;
	}
	*this = *src;
	return(out);
}

/**
Sets the size and allocates memory. Discards any existing contents if
the 3D object is not already 'empty'
**/

void Raw3D::sizer(int ixsize, int iysize, int izsize) {
	int i;
	if(z!=NULL)
		delete[]this->z;
	z = new Raw2D[izsize];			// make room for the 2D objects,
	zsize = izsize;
	for(i=0; i< zsize; i++) 
		z[i].sizer(ixsize,iysize);	// use the Raw2D sizer.	
}

void Raw3D::sizer(Raw3D* src)
{
	z=src->z;
	zsize=src->zsize;

}
/**
Copy all the pixels from 'src', resizing as needed. Do not change
name of this object.
**/ 
void Raw3D::wipecopy(Raw3D& src) {
	int k,kmax;

	if(&src==NULL)return;
	if(src.zsize==0) return;		// ignore empty inputs.
	if(src.getZsize() != zsize || src.getYsize() != z[0].getYsize() ||
		src.getXsize() != getXsize()) {
			sizer(&src);
	}
	kmax = getZsize();
	for(k=0; k< kmax; k++) {		// copy each field;
		z[k].wipecopy(&(src.z[k]));
	}
}

Raw2D operator/(const PIXTYPE val, const Raw2D& img)
{
	Raw2D res(img);
	for (int i = 0; i < img.size(); ++i)
		res.data[i] = val/res.data[i];
	return res;
}

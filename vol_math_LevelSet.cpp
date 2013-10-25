#include "vol_math_LevelSet.h"

#include <iostream>
#define pi 3.14
#include<math.h>

using namespace std;
PIXTYPE *p;
PIXTYPE	*q;
int i,j;

LevelSet::~LevelSet(void)
{
	cout<<"Destructor called"<<endl;


}
//*ImageF Dirac(ImageF x,float sigma);*/
Raw2D& regFunction(ImageF &s,int m,int n);
ImageF& operator+(ImageF &x,ImageF &y);
ImageF& operator+(ImageF &x,float s);
ImageF& cos(ImageF &x);
ImageF& operator *(int p1,ImageF &x);
ImageF& operator * (ImageF &x,ImageF &y);
ImageF* operator * (ImageF &x,ImageF *y);
ImageF& operator / (ImageF &x,float);
ImageF& sin(ImageF &s);
ImageF* Dirac(ImageF *x,float sigma)
{
	//showImg(*x);
	ImageF *ret=new ImageF(x);

	//Raw2D *f=new Raw2D((255*255*(1.0/2.0)/sigma)*(cos(3.14*(*ret/sigma))+1));
	Raw2D *f=new Raw2D(((1.0/2.0)/sigma)*(cos(3.14*(*ret/sigma))+1));
	//Raw2D *f=new Raw2D(((-1)**ret**ret)/4.0+255*2/3);
	ImageF *b=new Raw2D(regFunction(*ret,-sigma,sigma));
	*f=*f**b;
	ret=f;
	//showImg(*f);
	//IShowImg(*ret);
	return ret ;

}

ImageF* gradientx(ImageF *g)
{
	int n=g->getXsize();
	int m=g->getYsize();
	Raw2D* ret=new Raw2D(n, m);
	int i,j;

	for(int i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			if (i+1<n && j+1<m)
			{
				if (g->get(i+1,j)-g->get(i,j)>0)
				{
					ret->put(i,j,(g->get(i+1,j)-g->get(i,j)));
				}
				else if (g->get(i+1,j)-g->get(i,j)<0)
				{
					ret->put(i,j,(g->get(i+1,j)-g->get(i,j)));
				}
				else 
				{
					ret->put(i,j,0);
				}
			}
			else 
				ret->put(i,j,0);
		}
	}
	return ret;
}
ImageF* gradientxgorignal(ImageF *g)
{
	int n=g->getXsize();
	int m=g->getYsize();
	Raw2D* ret=new Raw2D(g);
	int i,j;
	int temp1,temp2;

	for(int i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{

			if(j>0)
				temp1=j-1;
			else 
				temp1=0;
			/*if(j<m-1)
				temp2=j+1;
			else 
				temp2=m-1;*/
			ret->put(i,j,(g->get(i,temp1+1)-g->get(i,temp1)));




		}
	}
	return ret;
}
ImageF* gradientxgc(ImageF *g)
{
	int n=g->getXsize();
	int m=g->getYsize();
	Raw2D* ret=new Raw2D(g);
	int i,j;
	int temp1,temp2;

	for(int i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{

			if(j>0)
				temp1=j-1;
			else 
				temp1=0;
			if(j<m-1)
				temp2=j+1;
			else 
				temp2=m-1;
			ret->put(i,j,(g->get(i,temp2)-g->get(i,temp1))/2.0);




		}
	}
	return ret;
}
ImageF* gradientxg(ImageF *g)
{
	int n=g->getXsize();
	int m=g->getYsize();
	Raw2D* ret=new Raw2D(g);
	int i,j;

	for(int i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{

			if (i+1<n&&j+1<m)
			{

				if (g->get(i+1,j)-g->get(i,j)>0)
				{
					ret->put(i,j,(g->get(i+1,j)-g->get(i,j)));
				}
				else if (g->get(i+1,j)-g->get(i,j)<0)
				{
					ret->put(i+1,j,(-1*g->get(i+1,j)+g->get(i,j)));
				}
				//else ret->put(i,j,0);
			}
		}
	}
		return ret;
}
ImageF* gradienty(ImageF* g )
{
	int n=g->getXsize();
	int m=g->getYsize();
	ImageF* ret=new ImageF(g->getXsize(),g->getYsize());
	int i,j;
	//PIXTYPE *gg=(*g).gety();
	//PIXTYPE *p=(*ret).gety();
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			if (i+1<n&&j+1<m)
			{
				if (g->get(i,j+1)-g->get(i,j)>0)
				{
					ret->put(i,j,(g->get(i,j+1)-g->get(i,j)));
				}
				else if (g->get(i,j+1)-g->get(i,j)<0)
				{
					ret->put(i,j+1,(g->get(i,j)-g->get(i,j+1)));
				}
				else ret->put(i,j,0);
			}
			else ret->put(i,j,0);
		}
	}
	return ret;

}
ImageF* gradientyg(ImageF* g )
{
	int n=g->getXsize();
	int m=g->getYsize();
	ImageF* ret=new ImageF(g);
	int i,j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			if (i+1<n&&j+1<m)
			{
				if (g->get(i,j+1)-g->get(i,j)>0)
				{
					ret->put(i,j,(g->get(i,j+1)-g->get(i,j)));
				}
				else if (g->get(i,j+1)-g->get(i,j)<0)
				{
					ret->put(i,j+1,(-1*g->get(i,j+1)+g->get(i,j)));
				}
				//else ret->put(i,j,0);
			}

		}
	}
	return ret;

}

ImageF* gradientygc(ImageF *g)
{
	int n=g->getXsize();
	int m=g->getYsize();
	Raw2D* ret=new Raw2D(g);
	int i,j;
	int temp1,temp2;

	for(int i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{

			if(i>0)
				temp1=i-1;
			else 
				temp1=0;
			if(i<n-1)
				temp2=i+1;
			else 
				temp2=n-1;
			ret->put(i,j,(g->get(temp2,j)-g->get(temp1,j))/2.0);




		}
	}
	return ret;
}
ImageF* gradientygcorignal(ImageF *g)
{
	int n=g->getXsize();
	int m=g->getYsize();
	Raw2D* ret=new Raw2D(g);
	int i,j;
	int temp1,temp2;

	for(int i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{

			if(i>0)
				temp1=i-1;
			else 
				temp1=0;
			//if(i<n-1)
			//	temp2=i+1;
			//else 
			//	temp2=n-1;
			int val=(g->get(temp1+1,j)-g->get(temp1,j));
			if (val<0)
			{
				cout<<val<<endl;
			}
			ret->put(i,j,val);




		}
	}
	return ret;
}
ImageF& cos(ImageF &xdata)
{
	ImageF *x=new ImageF(xdata);
	int m=x->getXsize();
	int n=x->getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			x->put(i,j,cos(float(x->get(i,j))));
		}

	}
	return *x;
}

ImageF* NeumannBoundCond(ImageF *img)
{
	int nrow=img->getXsize();
	int ncol=img->getYsize();
	int i,j;
	//the four point SDF
	img->putXY(0,img->getXY(2*ncol+2));
	img->putXY(ncol-1,img->get(2,ncol-3));
	img->putXY((nrow-1)*ncol,img->getXY((nrow-3)*ncol+2));
	img->putXY(nrow*ncol-1,img->get(nrow-3,ncol-3));
	//first and the last column SDF
	for(i=2;i<nrow-2;i++)
	{
		img->putXY((i-1)*ncol,img->getXY((i-1)*ncol+2));
		img->putXY(ncol*(i-1)+ncol-1,img->getXY(ncol*(i-1)-3));

	}
	//first and last row SDF
	for(j=2;j<ncol-2;j++)
	{
		img->putXY(j,img->getXY(2*ncol+j));
		img->putXY(ncol*(nrow-2)+j,img->getXY(ncol*(nrow-4)+j));
	}
	for (i=0;i<nrow;i++)
	{
		for(j=0;j<ncol;j++)
		{
			img->put(i,j,img->get(i,j));
		}
	}
	return img;

}

ImageF& ImageFSqrt(ImageF &x,ImageF &y)
{
	// 	ImageF x=xdata;
	// 	ImageF y=ydata;
	//p=x.gety();
	//q=y.gety();
	ImageF  *ret=new ImageF (x);
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			//*(p+i*n+j)=sqrt();
			ret->put(i,j,sqrt(float(x.get(i,j)*x.get(i,j)+y.get(i,j)*y.get(i,j))));
		}

	}
	return *ret;


}
ImageF& operator *(int p1,ImageF &x)
{

	//p=x.gety();
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			//*(p+i*n+j)=*(p+i*n+j)*p1;
			ret->put(i,j,x.get(i,j)*p1);
		}

	}
	return *ret;

}

ImageF& operator *(float p1,ImageF &x)
{

	//p=x.gety();
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{

			ret->put(i,j,x.get(i,j)*p1);
		}

	}
	return *ret;

}
ImageF& operator * (ImageF &x,ImageF &y)
{
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			ret->put(i,j,ret->get(i,j)*y.get(i,j));
		}

	}
	return *ret;

}
ImageF* operator * (ImageF &x,ImageF *y)
{
	ImageF *ret=new ImageF(x);
	int n=x.getXsize();
	int m=x.getYsize();
	for (int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			
			ret->put(i,j,x.get(i,j)*y->get(i,j));
		}

	}
	return ret;

}
//ImageF& operator/(ImageF &xdata,ImageF& ydata)
//{
//	/*array the corresponding position*/
//	//p=xdata.gety();
//	//q=ydata.gety();
//	ImageF *ret=new ImageF(xdata);
//	int m=xdata.getXsize();
//	int n=xdata.getYsize();
//	for (i=0;i<m;i++)
//	{
//		for(j=0;j<n;j++)
//		{
//			/*if(*(q+i*n+j)!=0)
//				*(p+i*n+j)=*(p+i*n+j)/(*(q+i*n+j));*/
//			if (ydata.get(i,j)!=0)
//			{
//				ret->put(i,j,xdata.get(i,j)/ydata.get(i,j));
//
//			}
//		}
//
//	}
//	return *ret;
//
//}

ImageF& operator/(ImageF &x,float t)

{
	//p=x.gety();
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	PIXTYPE temp;
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			temp=x.get(i,j)/t;
			ret->put(i,j,temp);
		}

	}
	return *ret;

}
ImageF& operator/(float t,ImageF &x)

{
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			ret->put(i,j,t/(x.get(i,j)+0.001));
		}

	}
	return x;

}
ImageF&  operator+( ImageF  &x,float s)
{
	//p=x.gety();
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	double val=0;
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			val=x.get(i,j)+s;
			if (0<val&&val<255)
			{
				ret->put(i,j,val);
			}
			else
			{
				ret->put(i,j,val);
			}
		}

	}
	return *ret;

}

//curvature
ImageF& div(ImageF &x, ImageF &y)
{
	ImageF *ret=new ImageF(x.getXsize(),x.getYsize());
	ImageF *gradx = gradientxg(&x),
		   *grady = gradientyg(&y);
	*ret = (*gradx + *grady);

	return *ret;
}

ImageF& operator+(ImageF &x,ImageF &y)
{
	int m=x.getXsize();
	int n=x.getYsize();
	ImageF *ret=new ImageF(m, n);
	double  val=0;
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			val=x.get(i,j)+y.get(i,j);
			if (val<256)
			{
				ret->put(i, j, val);
			}
			else 
			{
				ret->put(i,j,val);
			}
			
		}
	}
	return *ret;
}

ImageF& operator -(ImageF &x,ImageF &y)
{
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	int val=0;
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			val=x.get(i,j)-y.get(i,j);
			if (val>0)
			{
				ret->put(i,j,x.get(i,j)-y.get(i,j));
			}
			else 
			{
				ret->put(i,j,val);
			}
		
		}

	}
	return *ret;
}
ImageF& operator -(ImageF *x,ImageF &y)
{
	ImageF *ret=new ImageF(x);
	int m=x->getXsize();
	int n=x->getYsize();
	int val=0;
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			val=x->get(i,j)-y.get(i,j);
			if (val>0)
			{
				ret->put(i,j,x->get(i,j)-y.get(i,j));
			}
			else 
			{
				ret->put(i,j,val);
			}

		}

	}
	return *ret;
}
ImageF& operator -(ImageF &x,int s)
{
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	int val=0;
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			val=x.get(i,j)-s;
			if(val>0)
			ret->put(i,j,x.get(i,j)-s);
			else 
				ret->put(i,j,val);
		}

	}
	return *ret;
}
ImageF & operator/(ImageF &x,ImageF &y)
{
	ImageF *ret=new ImageF(x);
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			ret->put(i,j,(255*(float)x.get(i,j))/((float)y.get(i,j)));
		}

	}
	return *ret;


}
ImageF* del2(ImageF *phi)
{
	int m=phi->getXsize();
	int n=phi->getYsize();
	Raw2D *ret2=new Raw2D(m, n);

	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			
			if (i+1<m && j+1<n && i-1>=0 && j-1>=0)
			{
				int value = phi->get(i+1, j) + phi->get(i-1, j) + phi->get(i, j-1) + phi->get(i, j+1) - 4*(phi->get(i,j));
				if (value >= 0)
					ret2->put(i, j, (PIXTYPE)value);
				else
					ret2->put(i, j, value);
			}
			else 
			{
				ret2->put(i, j, 0);
			}
		}
	}
	return ret2;
}

ImageF& regFunction(ImageF &s,int m,int n)
{
	ImageF *ss=new ImageF(s);
	//p=s.gety();
	int l=ss->getXsize();
	int k=ss->getYsize();
	PIXTYPE val=0;
	for (i=0;i<l;i++)
	{
		for(j=0;j<k;j++)
		{
			val=ss->get(i,j);
			if(val>=m && val<=n)
			{
				//ss->put(i,j,255);//unsigned char version
				ss->put(i,j,1);
			}
			else if(ss->get(i,j)==m||ss->get(i,j)==n)
				/*	ss->put(i,j,val);
				else ss->put(i,j,val);*/
				ss->put(i,j,0);
			else ss->put(i,j,0);
		}

	}

	return *ss;

}
ImageF* distReg_p2(ImageF *phi)
{
	Raw2D *ret=new Raw2D(phi);
	int m=phi->getXsize();
	int n=phi->getYsize();
	Raw2D *phi_x=new Raw2D(gradientx(phi));
	ImageF *phi_y=new ImageF(gradienty(phi));
	//ImageF s=ImageFSqrt(((*phi_x)*(*phi_x)) + ((*phi_y)*(*phi_y)));
	ImageF  s=ImageFSqrt(*phi_x,*phi_y);
	ImageF a=regFunction(s,0,1);
	ImageF b=regFunction(s,0,1);//need to be changed.
	ImageF ps=a*sin(2*s)/(2*pi)+b*(s-1);
	ImageF dps=regFunction(ps,0,0)*ps+regFunction(ps,0,0)/(regFunction(s,0,0)+regFunction(s,0,0));
	ImageF f=div(dps**phi_x-*phi_x,dps**phi_y-*phi_y)+4**del2(phi);
	*ret=f;
	return ret;
}

ImageF& sin(ImageF &s)
{
	//overload sin
	//p=s.gety();
	ImageF *ss=new ImageF(s);
	int m=s.getXsize();
	int n=s.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			ss->put(i,j,255*(float)sin(float(s.get(i,j))));
		}

	}
	return *ss;
}

void write(Raw2D &destImg)
{
	//FILE *p=fopen("F:\\ls1.raw","ab+");
	//	fwrite(destImg.gety(),sizeof(PIXTYPE),destImg.getXsize()*destImg.getYsize(),p);
	//fclose(p);
}

ImageF& LevelSet::initialg(ImageF &g)
{
	ImageF *ret=new ImageF(&g);
	ImageF *gx=gradientxg(ret);
	ImageF *gy=gradientyg(ret);
	//g=*ret;
	g=*gx**gx+*gy**gy;
	//g=(255/(g+1));//or  unsigned char version
	//g=-255/(g+1)+255;unsigned char version
	//g=-1*g+255; unsigned char version
	g=1.0/(g+1.0);
	delete ret;
	return g;
}
//
//ImageF & pow(PIXTYPE e,ImageF &x)
//{
//	int m=x.getXsize();
//	int n=x.getYsize();
//	ImageF * ret=new ImageF(x);
//	for (int i=0;i<m;i++)
//	{
//		for (int j=0;j<n;j++)
//		{
//			ret->put(i,j,pow(e,x.get(i,j)));
//
//		}
//	}
//	return *ret;
//}

/// \brief
ImageF& LevelSet::drlse_edge(ImageF &phi_0,ImageF &g,float lambda,float mu,float alfa,float epsilon,int timestep, int iter,const char * potentialFunction)
{
	ImageF *phi=new ImageF(phi_0);
	int m=g.getXsize();
	int n=g.getYsize();
	ImageF *vx=gradientx(&g);
	ImageF *vy=gradienty(&g);
	Raw2D *diracPhi;

	for(int i=0;i<iter;i++)
	{
	
		//NeumannBoundCond(phi);
		IShowImg(*phi);
		ImageF &phi_x=*gradientx(phi);
		ImageF &phi_y=*gradienty(phi);
		Raw2D *s=new Raw2D(m,n);
		*s=ImageFSqrt(phi_x,phi_y);
		float smallNumber=1e-10;
		ImageF *Nx=new ImageF(phi_x);
		ImageF *Ny=new ImageF(phi_y);
		//*Nx=pow(255,phi_x/(*s+smallNumber));
		//*Ny=pow(255,phi_y/(*s+smallNumber));
		*Nx=(phi_x/(*s+smallNumber));
		*Ny=(phi_y/(*s+smallNumber));
		ImageF * curvature=new ImageF(m,n);
		*curvature=div(*Nx,*Ny);
		ImageF *distRegTerm=new ImageF(m,n);
		char *p1="single-well";
		if (strcmp(potentialFunction,p1)+1)
		{
		/*
		 compute distance regularization term in equation (13) 
		 with the single-well potential p1.
		 */
			*distRegTerm= (4*(*del2(phi)) - (*curvature));
		}
		
		else if (strcmp(potentialFunction,"double-well")+1)
		{
			*distRegTerm=distReg_p2(phi);  // compute the distance regularization term in eqaution (13) with the double-well potential p2.
		}
		else printf("EEROR");
		
		diracPhi=new Raw2D(Dirac(phi,epsilon));
		Raw2D* areaTerm=new Raw2D((g*(*diracPhi))); 
		ImageF  *edgeTerm=new ImageF(m,n);
		*edgeTerm = (*diracPhi) * ((*vx) * (*Nx))+((*vy) * (*Ny)) + (*diracPhi) * ( g * (*curvature));
		//*phi=*phi + timestep*(mu*(*distRegTerm) +lambda*(*edgeTerm) + alfa*(*areaTerm));
		//IShowImg(*distRegTerm);
		//IShowImg(*edgeTerm);
		//IShowImg(*areaTerm);
		//IShowImg(*phi);
		//IShowImg(4*(*del2(phi)));
		//IShowImg(*curvature);
		//IShowImg(phi_x);
		//IShowImg(phi_y);
		//IShowImg(*Nx);
		//IShowImg(*Ny);
		//IShowImg(g);
		//IShowImg(*s);
		phi_0=*phi;
	
	
	}	
	//showImg(phi_0);
	IShowImg(phi_0);
	return phi_0; 
}
//ImageF& LevelSet::uchar2double(ImageF &img)
//{
//	int m=img.getXsize();
//	int n=img.getYsize();
//	for (int i=0;i<m;i++)
//	{
//		for (int j=0;j<n;j++)
//		{
//			img.put(i,j,double(img.get(i,j)));
//		}
//	}
//}
void LevelSet::testout(Raw2D *ret){

	int i=0;

	while (i<100)
	{
		//PIXTYPE test = p[i];
		i++;
		cout<<"no====>>:"<<i<<"p[i]="<<ret->getXY(i)<<endl;
	}
}
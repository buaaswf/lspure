#include "LevelSet.h"
#include"ImageF.h"
#define pi 3.14
#include<math.h> ;
PIXTYPE  *p,*q;
int i,j;
LevelSet::LevelSet(void)
{
}


LevelSet::~LevelSet(void)
{
}
/*ImageF Dirac(ImageF x,float sigma);*/
ImageF regFunction(ImageF &s,int m,int n);
ImageF operator+(ImageF x,ImageF y);
ImageF operator+(ImageF x,float s);
ImageF cos(ImageF &x);
ImageF operator *(int p1,ImageF &x);
ImageF operator * (ImageF x,ImageF y);
ImageF sin(ImageF &s);
ImageF Dirac(ImageF &x,float sigma)
{
	ImageF f=(1/2/sigma)*(cos(pi*x/sigma)+1);
	ImageF b=regFunction(x,sigma,-sigma);
	f=f*b;
	x=f;
	return x;

}

ImageF gradientx(ImageF *g )
{
	int n=(*g).getXsize();
	int m=(*g).getYsize();
	int i,j;
	PIXTYPE *p=(*g).gety();
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			//*(p+i)=*(p+i-1)-*(p+i);
			//*(g.gety()+i)=*(g.gety()+i-1)-*(g.gety()+i);
			*(p+j+i*n)=*(p+j-1+i*n)-*(p+j+i*n);
			(*g).put(i,j,*(p+i*n+j));
		}
	}
	return *g;


}

ImageF gradienty(ImageF* g )
{
	int n=(*g).getXsize();
	int m=(*g).getYsize();
	int i,j;
	PIXTYPE *p=(*g).gety();
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			*(p+j+i*n)=*(p+j-1+i*n)-*(p+j+(i-1)*n);
			//*(g.gety()+i)=*(g.gety()+i-1)-*(g.gety()+i);
			(*g).put(i,j,*(p+i*n+j));
		}
	}
	return *g;



}
ImageF operator &(ImageF &x,ImageF &y)
{
	int i,j,n,m;
	n=x.width;
	m=x.length;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(x.gety()+i*m+j)=*(x.gety()+i*m+j)&(*(y.gety()+i*m+j));
		}

	}
	return x;

}
ImageF cos(ImageF &x)
{
	p=x.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=cos(float(*(p+i*n+j)));
		}
		
	}
	return x;
}

ImageF NeumannBoundCond(ImageF *img)
{
	int nrow=(*img).length;
	int ncol=(*img).width;
	int i,j;
	PIXTYPE *p=(*img).gety();
	//g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
	p[0]=p[3*nrow+2];
	p[ncol-1]=p[3*nrow+ncol-3];
	p[nrow*(ncol-1)]=p[(ncol-2)*nrow+2];
	p[nrow*ncol-1]=p[(nrow-2)*(ncol-2)-1];
	//g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);  
	for(i=1;i<nrow-2;i++)
	{
		p[i]=p[3*ncol+i];
		p[(ncol-1)*nrow+i-1]=p[(nrow-2)*ncol+i];
	}
	//g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  	
	for(j=2;j<ncol-2;j++)
	{
		p[ncol*(j-1)+1]=p[nrow*(j-1)+1];
		p[ncol*(j-1)+nrow-1]=p[nrow*(j-1)+nrow-1];
	
	}
	for (i=0;i<ncol;i++)
	{
		for(j=0;j<nrow;j++)
		{
			(*img).put(i,j,p[i*ncol+j]);
		}
	}
	return *img;

}

ImageF ImageFSqrt(ImageF &x,ImageF &y)
{
	p=x.gety();
	q=y.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=sqrt(float(*(p+i*n+j)^2+*(q+i*n+j)^2));
		}

	}
	return x;


}
ImageF operator *(int p1,ImageF &x)
{

	p=x.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)*p1;
		}

	}
	return x;

}
ImageF operator * (ImageF x,ImageF y)
{
	p=x.gety();
	q=y.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)*(*(q+i*n+j));
		}

	}
	return x;

}
ImageF operator/(ImageF &x,ImageF& y)
{
/*array the corresponding position*/
	
	p=x.gety();
	q=y.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)/(*(q+i*n+j));
		}

	}
	return x;

}

ImageF operator/(ImageF x,float t)

{
	p=x.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)/(t+0.001);
		}

	}
	return x;

}
ImageF operator+(ImageF x,float s)
{
	p=x.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)+s;
		}

	}
	return x;

}

ImageF div(ImageF x,ImageF y)
{
	return (x+y);
}
ImageF operator+(ImageF x,ImageF y)
{
	p=x.gety();
	q=y.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)+(*(q+i*n+j));
		}

	}
	return x;
}

ImageF operator -(ImageF x,ImageF y)
{
	p=x.gety();
	q=y.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)-(*(q+i*n+j));
		}

	}
	return x;
}
	
ImageF operator -(ImageF x,int s)
{
	p=x.gety();
	int m=x.getXsize();
	int n=x.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(p+i*n+j)=*(p+i*n+j)-s;
		}

	}
	return x;
}
ImageF del2(ImageF *phi)
{
	PIXTYPE *pt;
	p=phi->gety();
	
	int m=phi->getXsize();
	int n=phi->getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(pt+i*n+j)=*(p+i*n+j-1)+*(p+(i-1)*n+j)+*(p+i*n+j+1)+*(p+(i+1)*n+j)-4**(p+i*n+j);
		}

	}
	//ImageF phi_x=gradientx(phi);
	//ImageF phi_y=gradienty(phi);
/*
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			*(pt+i*n+j)=*(p+i*n+j-1)+*(p+(i-1)*n+j)+*(p+i*n+j+1)+*(p+(i+1)*n+j)-4**(p+i*n+j);
		}

	}*/
	

	return *phi;
}
ImageF regFunction(ImageF &s,int m,int n)
{
	p=s.gety();
	int l=s.getXsize();
	int k=s.getYsize();
	for (i=0;i<l;i++)
	{
		for(j=0;j<k;j++)
		{
			if(*(p+i*n+j)>m&&*(p+i*n+j)<n)
				*(p+i*n+j)=1;
			else if(*(p+i*n+j)==m||*(p+i*n+j)==n)
				*(p+i*n+j)=1;
			else *(p+i*n+j)=0;
		}

	}

		return s;

}
ImageF distReg_p2(ImageF *phi)
{
	ImageF *phi_x=new ImageF();
	ImageF *phi_y=new ImageF();
	*phi_x=gradientx(phi);
	*phi_y=gradienty(phi);
	//ImageF s=ImageFSqrt(((*phi_x)*(*phi_x)) + ((*phi_y)*(*phi_y)));
	ImageF  s=ImageFSqrt(*phi_x,*phi_y);
	ImageF a=regFunction(s,0,1);
	ImageF b=regFunction(s,0,1);//need to be changed.
	ImageF ps=a*sin(2*pi*s)/(2*pi)+b*(s-1);
	ImageF dps=regFunction(ps,0,0)*ps+regFunction(ps,0,0)/(regFunction(s,0,0)+regFunction(s,0,0));
	ImageF f=div(dps**phi_x-*phi_x,dps**phi_y-*phi_y)+4*del2(phi);
	*phi=f;
//ImageF b=(s>1);
	return *phi;
}
ImageF sin(ImageF &s)
{
	//overload sin
	p=s.gety();
	int m=s.getXsize();
	int n=s.getYsize();
	for (i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			 (*(p+i*n+j))=(float)sin(float(*(p+i*n+j)));
		}

	}
	return s;
}
ImageF drlse_edge(ImageF phi_0,ImageF g,float lambda,float mu,float alfa,float epsilon,int timestep, int iter,string sada,char * potentialFunction)
{
	ImageF *phi=new ImageF();
		*phi=phi_0;

	for(int i=0;i<iter;i++)
	{
		*phi=NeumannBoundCond(phi);
		ImageF *vx=new ImageF();
		ImageF *vy=new ImageF();
		(*vy)*(*vx);
		ImageF *phi_x=new ImageF();
		ImageF *phi_y=new ImageF();
		*vx=gradientx(&g);
		*vy=gradienty(&g);
		*phi_x=gradientx(&g);
		*phi_y=gradienty(&g);
		ImageF *s=new ImageF();
		*s=ImageFSqrt( *vx, *vy);
		float smallNumber=1e-10;
		ImageF *Nx,*Ny;
		*Nx=*phi_x/(*s+smallNumber);
		*Ny=*phi_y/(*s+smallNumber);
		ImageF * curvature=new ImageF();
		*curvature=div(*Nx,*Ny);
		ImageF distRegTerm;
		if (strcmp(potentialFunction,"single-well"))
			/*
			compute distance regularization term in equation (13) 
			with the single-well potential p1.
			*/
			distRegTerm= 4*del2(phi)-*curvature;
		//printf("");

		else if (strcmp(potentialFunction,"double-well"))
		{
			distRegTerm=distReg_p2(phi);  // compute the distance regularization term in eqaution (13) with the double-well potential p2.

			printf("asda");

		}
		else printf("EEROR");

		//float eplsion=0.1;
		ImageF diracPhi=Dirac(*phi,epsilon);
		ImageF areaTerm=diracPhi*g; 
		ImageF  *edgeTerm=new ImageF();
		*edgeTerm=diracPhi*(*vx**Nx+*vy**Ny) + diracPhi*g*(*curvature);
		//vx*vy;
		 *phi=*phi + timestep*(mu*distRegTerm + lambda**edgeTerm + alfa*areaTerm);
		return phi_0; 
	}
}
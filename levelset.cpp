#include "levelset.h"
#define _AFXDLL
#include <math.h>

#define PI       3.14159265358979323846
#define M_E      2.71828182845904523536

levelset::levelset(void)
:lambda(5.0)
,epsilon(1.5)
,alf(1.5)
,delt(5.0)
,numIter(30)
{
	mu = 0.2/delt;
	u0 = NULL;
	g = NULL;
	vx = NULL;
	vy = NULL;
	Nx = NULL;
	Ny = NULL;
	nxx = NULL;
	nyy = NULL;
	diracU = NULL;
	curvature = NULL;
	weightedLengthTerm = NULL;
	penalizingTerm = NULL;
	weightedAreaTerm = NULL;
}

levelset::levelset(CvMat *reference_area, CvMat *init_contour)
:lambda(5.0)
,epsilon(1.5)
,alf(1.5)
,delt(5.0)
,numIter(30)
{
	mu = 0.2/delt;
	u0 = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	g = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	vx =cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	vy = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	Nx = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	Ny = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	nxx = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	nyy = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	diracU = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	curvature =cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	weightedLengthTerm =cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	penalizingTerm = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);
	weightedAreaTerm = cvCreateMat(reference_area->rows,reference_area->cols,CV_64FC1);


	//cvSmooth(reference_area,g,CV_GAUSSIAN,3,0,0);
	u0 = cvCloneMat(init_contour);
	g = cvCloneMat(reference_area);

}


void levelset::initialize()
{

	// 	cvSobel(g,Nx,1,0);
	// 	cvSobel(g,Ny,0,1);

	//	sobel_x(g,Nx);
	//	sobel_y(g,Ny);

	gradient_x(g,Nx);
	gradient_y(g,Ny);

	double x,y;
	for (int i=0;i<g->rows;i++)
	{
		for (int j=0;j<g->cols;j++)
		{

			x = CV_MAT_ELEM(*Nx,double,i,j);
			y = CV_MAT_ELEM(*Ny,double,i,j);
			CV_MAT_ELEM(*g,double,i,j) = 1.0/(x*x+y*y+1);
		}
	}


	//cvSobel(g,vx,1,0);  // here the value of g is already changed!!!
	//cvSobel(g,vy,0,1);

	//sobel_x(g,vx);
	//sobel_y(g,vy);

	gradient_x(g,vx);
	gradient_y(g,vy);


}

levelset::~levelset(void)
{
	if (u0)
		cvReleaseMat(&u0);

	if (g)
		cvReleaseMat(&g);

	if (vx)
	{
		cvReleaseMat(&vx);
	}

	if (vy)
	{
		cvReleaseMat(&vy);
	}

	if (Nx)
	{
		cvReleaseMat(&Nx);
	}

	if (Ny)
	{
		cvReleaseMat(&Ny);
	}

	if (nxx)
	{
		cvReleaseMat(&nxx);
	}

	if (nyy)
	{
		cvReleaseMat(&nyy);
	}

	if (diracU)
	{
		cvReleaseMat(&diracU);
	}

	if (curvature)
	{
		cvReleaseMat(&curvature);
	}

	if (weightedAreaTerm)
	{
		cvReleaseMat(&weightedAreaTerm);
	}

	if (penalizingTerm)
	{
		cvReleaseMat(&penalizingTerm);
	}

	if (weightedLengthTerm)
	{
		cvReleaseMat(&weightedLengthTerm);
	}

}

void levelset::sobel_y(CvMat *source, CvMat *dest)
{
	int width, height;
	CvMat *temp_zero_padding;
	width = source->width;
	height= source->height;
	temp_zero_padding = cvCreateMat(height+2,width+2,CV_64FC1);
	cvSetZero(temp_zero_padding);

	for (int i=1;i<height+1;i++)
	{
		for (int j=1;j<width+1;j++)
		{
			CV_MAT_ELEM(*temp_zero_padding,double,i,j) = CV_MAT_ELEM(*source,double,i-1,j-1);
		}
	}


	int i_prime;
	int j_prime;

	for (int i=0;i<height;i++)
	{
		for (int j=0;j<width;j++)
		{
			i_prime = i+1;
			j_prime = j+1;

			double zuoxia = CV_MAT_ELEM(*temp_zero_padding,double,i_prime+1,j_prime-1);
			double xia = CV_MAT_ELEM(*temp_zero_padding,double,i_prime+1,j_prime);
			double youxia = CV_MAT_ELEM(*temp_zero_padding,double,i_prime+1,j_prime+1);

			double zuoshang = CV_MAT_ELEM(*temp_zero_padding,double,i_prime-1,j_prime-1);
			double shang = CV_MAT_ELEM(*temp_zero_padding,double,i_prime-1,j_prime);
			double youshang =CV_MAT_ELEM(*temp_zero_padding,double,i_prime-1,j_prime+1);

			double result = zuoshang+2*shang+youshang-zuoxia-2*xia-youxia;

			CV_MAT_ELEM(*dest,double,i,j)= result;



		}
	}

	cvReleaseMat(&temp_zero_padding);

}


void levelset::sobel_x(CvMat *source, CvMat *dest)
{

	int width, height;
	CvMat *temp_zero_padding;
	width = source->width;
	height= source->height;
	temp_zero_padding = cvCreateMat(height+2,width+2,CV_64FC1);
	cvSetZero(temp_zero_padding);

	for (int i=1;i<height+1;i++)
	{
		for (int j=1;j<width+1;j++)
		{
			CV_MAT_ELEM(*temp_zero_padding,double,i,j) = CV_MAT_ELEM(*source,double,i-1,j-1);
		}
	}


	int i_prime;
	int j_prime;

	for (int i=0;i<height;i++)
	{
		for (int j=0;j<width;j++)
		{
			i_prime = i+1;
			j_prime = j+1;

			double zuoshang = CV_MAT_ELEM(*temp_zero_padding,double,i_prime-1,j_prime-1);
			double zuo = CV_MAT_ELEM(*temp_zero_padding,double,i_prime,j_prime-1);
			double zuoxia = CV_MAT_ELEM(*temp_zero_padding,double,i_prime+1,j_prime-1);

			double youshang = CV_MAT_ELEM(*temp_zero_padding,double,i_prime-1,j_prime+1);
			double you = CV_MAT_ELEM(*temp_zero_padding,double,i_prime,j_prime+1);
			double youxia =CV_MAT_ELEM(*temp_zero_padding,double,i_prime+1,j_prime+1);

			double result = zuoshang+2*zuo+zuoxia-youshang-2*you-youxia;

			CV_MAT_ELEM(*dest,double,i,j)= result;



		}
	}

	cvReleaseMat(&temp_zero_padding);



}


void levelset::laplace(CvMat *source, CvMat *dest)
{

	int width, height;
	CvMat *temp_zero_padding;
	width = source->width;
	height= source->height;
	temp_zero_padding = cvCreateMat(height+2,width+2,CV_64FC1);
	cvSetZero(temp_zero_padding);

	for (int i=1;i<height+1;i++)
	{
		for (int j=1;j<width+1;j++)
		{
			CV_MAT_ELEM(*temp_zero_padding,double,i,j) = CV_MAT_ELEM(*source,double,i-1,j-1);
		}
	}


	int i_prime;
	int j_prime;

	for (int i=0;i<height;i++)
	{
		for (int j=0;j<width;j++)
		{
			i_prime = i+1;
			j_prime = j+1;

			double shang = CV_MAT_ELEM(*temp_zero_padding,double,i_prime-1,j_prime);
			double xia = CV_MAT_ELEM(*temp_zero_padding,double,i_prime+1,j_prime);
			double zuo = CV_MAT_ELEM(*temp_zero_padding,double,i_prime,j_prime-1);
			double you = CV_MAT_ELEM(*temp_zero_padding,double,i_prime,j_prime+1);
			double zhong = CV_MAT_ELEM(*temp_zero_padding,double,i_prime,j_prime);

			double result = shang + xia + zuo + you - 4*zhong;

			CV_MAT_ELEM(*dest,double,i,j)= result;

		}
	}

	cvReleaseMat(&temp_zero_padding);


}

void levelset::gradient_x(CvMat *source , CvMat *dest)
{
	CvMat stub1,stub2,stub3;

	int rows = source->rows;

	for (int i=0;i<rows;i++)
	{
		if (i==0)
		{
			cvGetRow(source,&stub1,0);
			cvGetRow(source,&stub2,1);

			cvGetRow(dest,&stub3,0);
			cvmSub(&stub2,&stub1,&stub3);

		}
		else if (i==rows-1)
		{
			cvGetRow(source,&stub1,rows-2);
			cvGetRow(source,&stub2,rows-1);

			cvGetRow(dest,&stub3,rows-1);
			cvmSub(&stub2,&stub1,&stub3);
		}
		else
		{
			cvGetRow(source,&stub1,i-1);
			cvGetRow(source,&stub2,i+1);

			cvGetRow(dest,&stub3,i);
			cvmSub(&stub2,&stub1,&stub3);
			cvConvertScale(&stub3,&stub3,0.5);
		}


	}

}

void levelset::gradient_y(CvMat *source , CvMat *dest)
{

	CvMat stub1,stub2,stub3;

	int cols = source->cols;

	for (int j=0;j<cols;j++)
	{
		if (j==0)
		{
			cvGetCol(source,&stub1,0);
			cvGetCol(source,&stub2,1);

			cvGetCol(dest,&stub3,0);
			cvmSub(&stub2,&stub1,&stub3);

		}
		else if (j==cols-1)
		{
			cvGetCol(source,&stub1,cols-2);
			cvGetCol(source,&stub2,cols-1);

			cvGetCol(dest,&stub3,cols-1);
			cvmSub(&stub2,&stub1,&stub3);
		}
		else
		{
			cvGetCol(source,&stub1,j-1);
			cvGetCol(source,&stub2,j+1);

			cvGetCol(dest,&stub3,j);
			cvmSub(&stub2,&stub1,&stub3);
			cvConvertScale(&stub3,&stub3,0.5);
		}


	}

}


void levelset::Dirac(CvMat *x, double sigma,CvMat *result)
{
	int nrows = x->rows;
	int ncols = x->cols;
	double temp;
	double temp2;

	for (int i=0;i<nrows;i++)
	{
		for (int j=0;j<ncols;j++)
		{
			temp = CV_MAT_ELEM(*x,double,i,j);
			temp2 = (1.0/2.0/sigma)*(1+cos(PI*temp/sigma));

			if (temp2<=sigma && temp2>=-sigma)
			{
				CV_MAT_ELEM(*result,double,i,j) = temp2;
			}
			else
				CV_MAT_ELEM(*result,double,i,j) = 0;

		}
	}


}

void levelset::curvature_central  (CvMat *nx , CvMat *ny , CvMat *result)
{
	// 	cvSobel(nx,nxx,1,0);
	// 	cvSobel(ny,nyy,0,1);
	//	sobel_x(nx,nxx);
	//	sobel_y(ny,nyy);
	gradient_x(nx,nxx);
	gradient_y(ny,nyy);
	cvmAdd(nxx,nyy,result);
}

void levelset::NeumannBoundCond(CvMat *f)
{
	int nrow = f->rows;
	int ncol = f->cols;

	CV_MAT_ELEM(*f,double,0,0) = CV_MAT_ELEM(*f,double,2,2);
	CV_MAT_ELEM(*f,double,0,ncol-1) = CV_MAT_ELEM(*f,double,2,ncol-3);
	CV_MAT_ELEM(*f,double,nrow-1,0) = CV_MAT_ELEM(*f,double,nrow-3,2);
	CV_MAT_ELEM(*f,double,nrow-1,ncol-1) = CV_MAT_ELEM(*f,double,nrow-3,ncol-3);

	for (int j=1;j<ncol-1;j++)
	{
		CV_MAT_ELEM(*f,double,0,j) = CV_MAT_ELEM(*f,double,2,j);
		CV_MAT_ELEM(*f,double,nrow-1,j) = CV_MAT_ELEM(*f,double,nrow-3,j);
	}

	for (int i=1;i<nrow-1;i++)
	{
		CV_MAT_ELEM(*f,double,i,0) = CV_MAT_ELEM(*f,double,i,2);
		CV_MAT_ELEM(*f,double,i,ncol-1) = CV_MAT_ELEM(*f,double,i,ncol-3);
	}



}

void  levelset::evolution_levelset(CvMat *u0, CvMat *g, double lambda, double mu, double alf, double epsilon, double delt, int numIter)
{
	CvScalar myscalar = cvScalar(pow(10,-10.0));

	for (int i=0;i<numIter;i++)
	{
		NeumannBoundCond(u0);
		cvSetZero(curvature);

		// 		cvSobel(u0,nxx,1,0);
		// 		cvSobel(u0,nyy,0,1);

		//		sobel_x(u0,nxx);
		//		sobel_y(u0,nyy);

		gradient_x(u0,nxx);
		gradient_y(u0,nyy);

		cvMul(nxx,nxx,Nx);
		cvMul(nyy,nyy,Ny);
		cvmAdd(Nx,Ny,curvature);

		cvAddS(curvature,myscalar,curvature);

		cvPow(curvature,curvature,0.5);
		cvDiv(nxx,curvature,Nx);
		cvDiv(nyy,curvature,Ny);

		Dirac(u0,epsilon,diracU);
		curvature_central(Nx,Ny,curvature);

		cvMul(vx,Nx,nxx);
		cvMul(vy,Ny,nyy);
		cvmAdd(nxx,nyy,nxx);
		cvMul(g,curvature,nyy);
		cvmAdd(nxx,nyy,nxx);

		cvMul(diracU,nxx,weightedLengthTerm,lambda);

		//	cvLaplace(u0,nxx);
		laplace(u0,nxx);
		cvmSub(nxx,curvature,nyy);

		cvConvertScale(nyy,penalizingTerm,mu);
		cvMul(diracU,g,weightedAreaTerm,alf);    //calculate the weightedAreaTerm

		cvmAdd(weightedLengthTerm,weightedAreaTerm,nxx);
		cvmAdd(nxx,penalizingTerm,nyy);

		cvConvertScale(nyy,nxx,delt);
		cvmAdd(u0,nxx,u0);

	}
}



void levelset::evolution_for_outsidecall()
{
	evolution_levelset(u0 , g, lambda, mu, alf , epsilon , delt,30);
}
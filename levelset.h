#pragma once
#ifndef LEVEL_H
#define LEVEL_H
#include "cv.h"
#include"highgui.h"

class levelset
{
public:
	levelset(void);
	levelset(CvMat *reference_area,CvMat *init_contour);
	~levelset(void);


	CvMat *u0; //level set function to be updated
	CvMat *g;  //edge indicator function

	CvMat *vx;
	CvMat *vy;
	CvMat *Nx;
	CvMat *Ny;
	CvMat *nxx;
	CvMat *nyy;

	CvMat *diracU;
	CvMat *curvature;
	CvMat *weightedLengthTerm;
	CvMat *penalizingTerm;
	CvMat *weightedAreaTerm;


	double lambda; //coefficient of the weighted length term L(\phi)
	double mu; //coefficient of the internal (penalizing) energy term P(\phi)
	double alf;  //coefficient of the weighted area term A(\phi), choose smaller alf 
	double epsilon;  //the parameter in the definition of smooth Dirac function, default value 1.5
	double delt;  //time step of iteration, see the paper for the selection of time step and mu 
	int numIter; //number of iterations. 

private:
	void sobel_x(CvMat *source, CvMat *dest);
	void sobel_y(CvMat *source, CvMat *dest);
	void laplace(CvMat *source, CvMat *dest);
	void gradient_x(CvMat *source , CvMat *dest);
	void gradient_y(CvMat *source , CvMat *dest);
	void Dirac(CvMat *x , double sigma, CvMat *result);
	void curvature_central(CvMat *nx , CvMat *ny , CvMat *result);
	void NeumannBoundCond(CvMat *f);
	void evolution_levelset(CvMat *u0 , CvMat *g, double lambda, double mu, double alf , double epsilon , double delt, int numIter);

public:
	void initialize(); //calculate the each value of the matrices
	void evolution_for_outsidecall();
};
#endif
#ifndef _FILTER_H_
#define _FILTER_H_

#include "vol_math_Raw3D_Independt.h"

class Filter
{
public:
	static Raw2D* guassFilter(Raw2D* src, int halfsize)
	{
		int i=0,j=0,m=0,n=0,width=0,length=0,index=0;
		float sum = 0;
		int delta=1;
		width=src->getXsize();
		length=src->getYsize();

		Raw2D *guass = new Raw2D(width, length);	///< The result of guass filter. 

		for (i=0;i<width;i++)
		{
			for (j=0;j<length;j++)
			{
				sum=0;
				float weight=0, total=0;
				for( m=i-halfsize; m<=i+halfsize; m++)
				{
					for(n=j-halfsize; n<=j+halfsize; n++)
					{
						if(m >= 0 && m < width && n>=0 && n < length) 
						{
							//weight=1.0f/((m-i)*(m-i)+(n-i)*(n-i)+1);
							weight=1.0f/exp((float)((m-i)*(m-i)+(n-i)*(n-i)));
							sum += weight*(src->get(m, n));
							total += weight;
						}
					}
				}

				if(total!=0)
				{	
					sum /= total;//total is 1,regulation
					guass->put(i, j, (PIXTYPE)sum);		
				}
				else  //should never come here
				{
					cout << "total==0" << endl;
				}
			}
		}
		
		return guass;
	}

	//Trilateral filter consisting of gradient filter, adaptive neighborhood
	//computation and detail filter
	void TrilateralFilter(Raw2D* srcImg, PIXTYPE sigmaC); 

	//Computes X and Y gradients of the input image
	void ComputeGradients(Raw2D* pX, Raw2D* pY); 

	//Bilaterally filters  gradients pX and pY 
	void BilateralGradientFilter(Raw2D* pX, Raw2D* pY, Raw2D* pSmoothX, 
		Raw2D* pSmoothY, PIXTYPE sigmaC, PIXTYPE sigmaS, int filterSize); 

	//Builds the stack of min-max image gradients; returns the range variance
	PIXTYPE buildMinMaxImageStack(Raw2D* pX, Raw2D* pY, Raw3D* pMinStack,
		Raw3D* pMaxStack , int levelMax, PIXTYPE beta); 

	//Finds the adaptive neighborhood size (stack level) 
	//from the min-max gradient stack
	void findAdaptiveRegion(Raw3D* pMinStack, Raw3D* pMaxStack, PIXTYPE R, int levelMax); 

	//Filters the detail signal and computes the final output image	
	void DetailBilateralFilter(Raw2D* srcImg, Raw2D* pSmoothX, Raw2D* pSmoothY, 
		Raw2D* fTheta, float sigmaCTheta, float sigmaRTheta); 
};

#endif
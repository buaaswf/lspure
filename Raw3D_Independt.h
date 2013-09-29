#pragma once
#define PIXTYPE unsigned char
class Raw3D;
class Raw2D  {
	private:   			//-----------------DATA----------------- 
	int xsize;		// # of pixels per scanline,
	int ysize;		// # of scanlines in this Raw2D.
	PIXTYPE *y;		// 1D array of PIXTYPE that are accessed as a 2D array.
	
public:				//---------------init fcns-------------
	Raw2D(int,int,PIXTYPE *);	
	Raw2D(void);// constructor for 'empty' Raw2Ds
	~Raw2D(void);		// destructor; releases memory
	void sizer(int ixsize, int iysize);	// get mem for rectangle of pixels
void guassConv(Raw2D *raw2d,int halfsize);	//gauss filter
	void sizer(Raw2D* src);					// get same amt. of mem as 'src'
	int getXsize(void) {return xsize;};		// get # pixels per scanline
	int getYsize(void) {return ysize;};		// get # of scanlines.
	int wipecopy(Raw2D* src);			// copy, even with size mismatch change from bool swf 2013 4 16
	void put(int ix, int iy, PIXTYPE val) {	// write 'val' at location ix,iy.
		y[ix + xsize*iy] = val; 
	};
	inline PIXTYPE get(int ix, int iy) {	// read the value at ix,iy.
		int index=ix + xsize*iy;
		return(y[index]); 
	};
	PIXTYPE getXY(int ixy){		// read value at 1D address ixy
		return y[ixy];
	};
	void putXY(int ixy,PIXTYPE val){// write value at 1D address ixy
		y[ixy] = val;
	};
	PIXTYPE *gety()
	{
		return y;

	};
	//---------------Trilateral Filter fcns-------------

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
class Raw3D {
public:
	Raw2D *z;	// dynam. allocated space for a set of Raw2D objects.
	int zsize;	// # of Raw2D objects stored.

public:							
	Raw3D(void);// 'empty' Raw3D constructor.
	Raw3D(int zsize,Raw2D *src);//swf add for read data 
	~Raw3D(void);	// destructor.
	void sizer(int ixsize, int iysize, int izsize); // reserve memory
	void sizer(Raw3D* src);			// get same amt. of mem as 'src
	int getZsize(void) {				// return # of Raw2D's we hold;
		return(zsize); 
	};
	int getYsize() {					// # of Raw1D's in zval-th Raw2D;
		return(z[0].getYsize()); 
	};
	int getXsize(){						// # of pixels on yval,zval-th line
		return(z[0].getXsize()); 
	};
	PIXTYPE get(int ix, int iy, int iz) {
		return(z[iz].get(ix,iy));	// write 'val' at location ix,iy,iz. 
	};
	void put(int ix, int iy, int iz, PIXTYPE val) { 
		z[iz].put(ix,iy,val);		//write 'val' at location ix,iy,iz.
	};
	void wipecopy(Raw3D& src);			// copy, resize as needed.
};
class Raw3D_Independt
{
public:
	Raw3D_Independt(void);
	~Raw3D_Independt(void);
};


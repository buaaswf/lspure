#ifndef _VOL_MATH_RAW3D_INDEPENDT_H_
#define _VOL_MATH_RAW3D_INDEPENDT_H_

#ifndef Raw2D_H
#define Raw2D_H

#include <iostream>
using namespace std;

typedef double PIXTYPE;
class Raw3D;
class Raw2D
{
private:   			//-----------------DATA----------------- 
	int xsize;		// # of pixels per scanline,
	int ysize;		// # of scanlines in this Raw2D.
	PIXTYPE *data;		// 1D array of PIXTYPE that are accessed as a 2D array.

	void putData(const PIXTYPE *res)
	{
		if (data != NULL)
			delete [] data;
		memcpy(data, res, xsize*ysize);
	}

public:				//---------------init fcns-------------
	const static int MAXPIXEL = 255;
	const static int  MINPIXEL = 0;

	Raw2D(int, int, PIXTYPE*);
	Raw2D(int, int);
	Raw2D(Raw2D* r);
	Raw2D(const Raw2D& r);
	Raw2D(void);		// constructor for 'empty' Raw2Ds
	~Raw2D(void);		// destructor; releases memory

	int size() const { return xsize*ysize; }

	/// \brief Swap data with img
	Raw2D& swap(Raw2D& img)
	{
		std::swap(this->xsize, img.xsize);
		std::swap(this->ysize, img.ysize);
		std::swap(this->data, img.data);
		return img;
	}

	Raw2D& operator=(Raw2D img)
	{
		img.swap(*this);
		return *this;
	}

	Raw2D& operator+=(const Raw2D &img)
	{
		for (int i = 0; i < size(); ++i)
			this->data[i] += img.data[i];
		return *this;
	}

	Raw2D & operator+=(const PIXTYPE val)
	{
		for (int i = 0; i < size(); ++i)
		{
			this->data[i] += val;
		}
		return *this;
	}

	Raw2D operator+(const Raw2D &img)
	{
		return Raw2D(*this) += img;
	}

	Raw2D operator+(const PIXTYPE val)
	{
		return Raw2D(*this) += val;
	}

	Raw2D & operator-=( const Raw2D &img)
	{
		for (int i = 0; i<size();++i)
			this->data[i] -= img.data[i];
		return *this;
	}

	Raw2D & operator -= ( const PIXTYPE val)
	{
		for (int i=0;i<size();i++)
			this->data[i]-=val;
		return *this;
	}

	Raw2D operator -(const Raw2D &img)
	{
		return Raw2D(*this)-=img;
	}

	Raw2D operator - (const PIXTYPE val)
	{
		return Raw2D(*this) -= val;
	}

	Raw2D& operator*=(const Raw2D& img)
	{
		for (int i = 0; i < size(); ++i)
			this->data[i] *= img.data[i];
		return *this;
	}

	Raw2D& operator*=(const PIXTYPE val)
	{
		for (int i = 0; i < size(); ++i)
			this->data[i] *= val;
		return *this;
	}

	Raw2D operator*(const Raw2D& img)
	{
		return  Raw2D(*this) *= img;
	}

	Raw2D operator*(const PIXTYPE val)
	{
		return Raw2D(*this) *= val;
	}

	Raw2D& operator/=(const Raw2D& img)
	{
		for (int i = 0; i < size(); ++i)
			this->data[i] /= img.data[i];
		return *this;
	}

	Raw2D& operator/=(const PIXTYPE val)
	{
		for (int i = 0; i < size(); ++i)
			this->data[i] /= val;
		return *this;
	}

	Raw2D operator/(const Raw2D& img)
	{
		return Raw2D(*this) /= img;
	}

	Raw2D operator/(const PIXTYPE val)
	{
		return Raw2D(*this) *= val;
	}

	friend Raw2D operator/(const PIXTYPE val, const Raw2D& img);

	void sizer(int ixsize, int iysize);	// get mem for rectangle of pixels
	void sizer(Raw2D* src);					// get same amt. of mem as 'src'
	int getXsize(void) const {return xsize;}		// get # pixels per scanline
	int getYsize(void) const {return ysize;}		// get # of scanlines.

	inline void put(int ix, int iy, PIXTYPE val)	// write 'val' at location ix,iy.
	{
#ifdef _DEBUG  //only check under debug mode
		if (iy + ysize*ix < size())
		{
#endif
			data[iy + ysize*ix] = val;
#ifdef _DEBUG
		}
		else 
			cout<<"out of size put"<<endl;
#endif
	}

	inline PIXTYPE get(int ix, int iy) {	// read the value at ix,iy.
#ifdef _DEBUG
		if(iy + ysize*ix<=size())
		{
#endif
			return data[iy + ysize*ix]; 
#ifdef _DEBUG
		}
		else 
		{
			cout<<"out of size get"<<endl;
		}
#endif
	}

	PIXTYPE getXY(int ixy)
	{		// read value at 1D address ixy
		if(ixy<xsize*ysize)
		{
			return data[ixy];
		}
		else 
		{
			cout<<"out of size get "<<endl;
		}
	}

	void putXY(int ixy,PIXTYPE val){// write value at 1D address ixy
		if (ixy<xsize*ysize)
		{
			data[ixy] = val;
		}
		else cout<<"out of size putxy"<<endl;

	}
	bool Raw2D::wipecopy(Raw2D* src);

};

#endif  //Raw2D_H

#ifndef Raw3D_H
#define Raw3D_H

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
#endif

#ifndef Raw3D_Independt_H
#define Raw3D_Independt_H
class Raw3D_Independt
{
public:
	Raw3D_Independt(void);
	~Raw3D_Independt(void);
};
#endif

#endif //_VOL_MATH_RAW3D_INDEPENDT_H_
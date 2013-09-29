#include "ImageF.h"
#include "Raw3D_independt.h"

ImageF::ImageF(int length,int width,PIXTYPE *p):Raw2D(length,width,p)
{


}
ImageF::ImageF()
{
	this->length=length;
	this->width=width;

}

ImageF::~ImageF(void)
{
}

#include "image.h"
//#include <mbstring.h>
#include <stdio.h>
#include <math.h>
//# include "vol_math_trilateralfilter.h"
#include<iostream>
using namespace std;
image::image(int length,int height,int width)
{
	this->height=height;
	this->length=length;
	this->width=width;
	int length1=this->getlength();
	this->buf=(u_char*)malloc(length1);
}

void image::readImage(u_char* buf,char const *file ,int size)
{
	//this->height=height;
	//this->length=length;
	//this->width=width;
	//long size=size;
	FILE *op;
		op=fopen(file,"rb");
	if(op==NULL)
	{
		printf("open fail");
	}
	fread(buf,sizeof(u_char),size,op);
	fclose(op);
	printf("read is ok");
}
float* image::buf2float(u_char *buf)
{
	u_char *p;
	p=buf;
	int i=0;
	long length=this->getlength();
	//p=p[length]
	//int *imageval=new int[length];
	float *imgf=new float[length];
	while(p)
	{
		
		//*(buf+i)=*p;
		imgf[i]= (float)p[i] ;
		imgf+=i;
		p++;i++;
	}
	return imgf;
}	

/** 
** method to remove sharp the raw image with unsharp mask 
* @param gray input grayscale binary array  
* @param smooth output data for smooth result, the memory need to be allocated outside of the function 
* @param width width of the input grayscale image 
* @param height height of the input grayscale image 
*/ 
void sharpenImage  (unsigned char* gray, unsigned char* smooth, int width, int height,int length)  
{  
      
    int templates[25] = { -1, -4, -7, -4, -1,   
        -4, -16, -26, -16, -4,   
        -7, -26, 505, -26, -7,  
        -4, -16, -26, -16, -4,   
        -1, -4, -7, -4, -1 };         
    memcpy ( smooth, gray, width*height*sizeof(unsigned char) );  
    for (int j=2;j<height-2;j++)  
    {  
        for (int i=2;i<width-2;i++)  
        {  
            int sum = 0;  
            int index = 0;  
            for ( int m=j-2; m<j+3; m++)  
            {  
                for (int n=i-2; n<i+3; n++)  
                {  
                    sum += gray [ m*width + n] * templates[index++] ;  
                }  
            }  
            sum /= 273;  
            if (sum > 255)  
                sum = 255;  
            if (sum <0)  
                sum = 0;  
            smooth [ j*width+i ] = sum;  
        }  
    }  
} 
//float image::checkresult()
image::~image(void)
{
	//free(buf);
	cout<<"free buf"<<endl;
}

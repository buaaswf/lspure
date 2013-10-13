#pragma once
#define u_char unsigned char
class image
{
private:
	int length;
	int width;
	int height;
	public:
	u_char * buf;
	image(int length,int height,int width);
	int getlength()
	{
		return (height*length*width);
	}
	~image(void);
	void readImage( u_char* buf,char const *file ,int size);
	void sharpenImage  (unsigned char* gray, unsigned char* smooth, int width, int height,int length) ;
	float * buf2float(u_char *buf);
	void save();

};


// levelset_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "levelset_test.h"
#include "levelset.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// The one and only application object

CWinApp theApp;

using namespace std;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;
	IplImage *levelset_test;
	CvMat *reference_image;
	CvMat *init_contour;

	// initialize MFC and print and error on failure
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: change error code to suit your needs
		_tprintf(_T("Fatal Error: MFC initialization failed\n"));
		nRetCode = 1;
	}
	else
	{
		levelset_test = cvLoadImage("twocells.bmp",0);
		cvSmooth(levelset_test,levelset_test,CV_GAUSSIAN,3,0,0);
		reference_image = cvCreateMat(levelset_test->height,levelset_test->width,CV_64FC1);
		init_contour = cvCreateMat(levelset_test->height,levelset_test->width,CV_64FC1);

		int row = reference_image->rows;
		int col = init_contour->cols;
		int widthstep = levelset_test->widthStep;

		for (int i=0;i<row;i++)
		{
			for (int j=0;j<col;j++)
			{
				int index = i*widthstep+j;
				int pixel_value = (unsigned char)levelset_test->imageData[index];
				CV_MAT_ELEM(*reference_image,double,i,j) = (double)pixel_value;
				if ( i < 8 || i > row-9 || j < 8 || j > col -9 )
				{
					CV_MAT_ELEM(*init_contour,double,i,j) = 4.0;
				}
				else
				{
					CV_MAT_ELEM(*init_contour,double,i,j) = -4.0;

				}


			}
		}

		for (int j=8;j<col-8;j++)
		{
			CV_MAT_ELEM(*init_contour,double,8,j) = 0;
			CV_MAT_ELEM(*init_contour,double,row-9,j) = 0;
		}

		for (int i=8;i<row-8;i++)
		{
			CV_MAT_ELEM(*init_contour,double,i,8) = 0;
			CV_MAT_ELEM(*init_contour,double,i,col-9) = 0;

		}
		levelset test_levelset(reference_image,init_contour);
		test_levelset.initialize();
		test_levelset.evolution_for_outsidecall();
		
		init_contour = cvCloneMat(test_levelset.u0);
		for (int i=0;i<row;i++)
		{
			for (int j=0;j<col;j++)
			{
				double aa = CV_MAT_ELEM(*init_contour,double,i,j);
				int index =  i*widthstep+j;
				if (aa < 0)
				{
					levelset_test->imageData[index]=255;
					
				}
				else
					levelset_test->imageData[index]=0;

			}
		}

		cvNamedWindow("result");
		for (;;)
		{

			cvShowImage("result",levelset_test);
			char c = cvWaitKey(10);
			if (c=='q')
				break;
		}
		
		
		cvDestroyWindow("result");
		cvReleaseImage(&levelset_test);
		cvReleaseMat(&reference_image);
		cvReleaseMat(&init_contour);
		
	}

	return nRetCode;
}

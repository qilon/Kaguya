/*******************************************************
\Author:	Qingxiong Yang
\Function:	Basic functions.
********************************************************/
#ifndef QX_CVPR09_CTBF_BASIC_H
#define QX_CVPR09_CTBF_BASIC_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <vector>
//#include <process.h>
// #include <direct.h>
// #include <io.h>
#include <time.h>
#include <string>
#include <memory.h>
#include <algorithm>
#include <functional>      // For greater<int>()
#include <iostream>
#if _MSC_VER > 1020   // if VC++ version is > 4.2
   using namespace std;  // std c++ libs implemented in std
#endif
#define QX_DEF_PADDING					10
#define QX_DEF_THRESHOLD_ZERO			1e-6
class   qx_timer		{public: void start();	float stop(); void time_display(char *disp=""); void fps_display(char *disp=""); private: clock_t m_begin; clock_t m_end;};

inline float qx_max_f3(float*a){return(max(max(a[0],a[1]),a[2]));}
inline float qx_min_f3(float*a){return(min(min(a[0],a[1]),a[2]));}
/*Box filter*/
void boxcar_sliding_window(double **out,double **in,double **temp,int h,int w,int radius);
void boxcar_sliding_window_x(double *out,double *in,int h,int w,int radius);
void boxcar_sliding_window_y(double *out,double *in,int h,int w,int radius);
/*Gaussian filter*/
int gaussian_recursive(double **image,double **temp,double sigma,int order,int h,int w);
void gaussian_recursive_x(double **od,double **id,int w,int h,double a0,double a1,double a2,double a3,double b1,double b2,double coefp,double coefn);
void gaussian_recursive_y(double **od,double **id,int w,int h,double a0,double a1,double a2,double a3,double b1,double b2,double coefp,double coefn);
/*basic functions*/

//inline float min(float a,float b){if(a<b) return(a); else return(b);}
//inline float max(float a,float b){if(a>b) return(a); else return(b);}
inline int qx_sum_u3(unsigned char *a) {return(a[0]+a[1]+a[2]);}
inline double qx_sum_d3(double*a){return(a[0]+a[1]+a[2]);}
inline unsigned char qx_min_u3(unsigned char *a){return(min(min(a[0],a[1]),a[2]));}
inline unsigned char qx_max_u3(unsigned char *a){return(max(max(a[0],a[1]),a[2]));}
inline void image_zero(float **in,int h,int w,float zero=0){memset(in[0],zero,sizeof(float)*h*w);}
inline void image_zero(double **in,int h,int w,double zero=0){memset(in[0],zero,sizeof(double)*h*w);}
inline void image_zero(unsigned char**in,int h,int w,unsigned char zero=0){memset(in[0],zero,sizeof(unsigned char)*h*w);}
inline void image_zero(double ***in,int h,int w,int d,double zero=0){memset(in[0][0],zero,sizeof(double)*h*w*d);}
inline unsigned char rgb_2_gray(unsigned char*in)
{
	return (unsigned char)(0.299*in[0]+0.587*in[1]+0.114*in[2]+0.5);
}
inline int qx_square_difference_u3(unsigned char *a,unsigned char *b){int d1,d2,d3; d1=(*a++)-(*b++); d2=(*a++)-(*b++);	d3=(*a++)-(*b++); return(int(d1*d1+d2*d2+d3*d3));}
void qx_specular_free_image(unsigned char ***image_specular_free,unsigned char ***image_normalized,float **diffuse_chromaticity_max,int h,int w);
inline double *get_color_weighted_table(double sigma_range,int len)
{
	double *table_color,*color_table_x; int y;
	table_color=new double [len];
	color_table_x=&table_color[0];
	for(y=0;y<len;y++) (*color_table_x++)=exp(-double(y*y)/(2*sigma_range*sigma_range));
	return(table_color);
}
inline void color_weighted_table_update(double *table_color,double dist_color,int len)
{
	double *color_table_x; int y;
	color_table_x=&table_color[0];
	for(y=0;y<len;y++) (*color_table_x++)=exp(-double(y*y)/(2*dist_color*dist_color));
}

inline void vec_min_val(float &min_val,float *in,int len)
{
	min_val=in[0];
	for(int i=1;i<len;i++) if(in[i]<min_val) min_val=in[i];	
}
inline void vec_min_val(unsigned char &min_val,unsigned char *in,int len)
{
	min_val=in[0];
	for(int i=1;i<len;i++) if(in[i]<min_val) min_val=in[i];	
}
inline void vec_max_val(float &max_val,float *in,int len)
{
	max_val=in[0];
	for(int i=1;i<len;i++) if(in[i]>max_val) max_val=in[i];	
}
inline void vec_max_val(unsigned char &max_val,unsigned char *in,int len)
{
	max_val=in[0];
	for(int i=1;i<len;i++) if(in[i]>max_val) max_val=in[i];	
}
inline void down_sample_1(unsigned char **out,unsigned char **in,int h,int w,int scale_exp)
{
	int y,x; int ho,wo; unsigned char *out_y,*in_x;
	ho=(h>>scale_exp); wo=(w>>scale_exp); 
	for(y=0;y<ho;y++)
	{
		out_y=&out[y][0]; in_x=in[y<<scale_exp];
		for(x=0;x<wo;x++) *out_y++=in_x[x<<scale_exp];
	}
}
inline void down_sample_1(float**out,float**in,int h,int w,int scale_exp)
{
	int y,x; int ho,wo; float *out_y,*in_x;
	ho=(h>>scale_exp); wo=(w>>scale_exp); 
	for(y=0;y<ho;y++)
	{
		out_y=&out[y][0]; in_x=in[y<<scale_exp];
		for(x=0;x<wo;x++) *out_y++=in_x[x<<scale_exp];
	}
}
inline double qx_linear_interpolate_xy(double **image,double x,double y,int h,int w)
{
	int x0,xt,y0,yt; double dx,dy,dx1,dy1,d00,d0t,dt0,dtt;
	x0=int(x); xt=min(x0+1,w-1); y0=int(y); yt=min(y0+1,h-1);
	dx=x-x0; dy=y-y0; dx1=1-dx; dy1=1-dy; d00=dx1*dy1; d0t=dx*dy1; dt0=dx1*dy; dtt=dx*dy;
	return(d00*image[y0][x0]+d0t*image[y0][xt]+dt0*image[yt][x0]+dtt*image[yt][xt]);
}
/*memory*/
inline double *** qx_allocd_3(int n,int r,int c,int padding=QX_DEF_PADDING)
{
	double *a,**p,***pp;
    int rc=r*c;
    int i,j;
	a=(double*) malloc(sizeof(double)*(n*rc+padding));
	if(a==NULL) {printf("qx_allocd_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
    p=(double**) malloc(sizeof(double*)*n*r);
    pp=(double***) malloc(sizeof(double**)*n);
    for(i=0;i<n;i++) 
        for(j=0;j<r;j++) 
            p[i*r+j]=&a[i*rc+j*c];
    for(i=0;i<n;i++) 
        pp[i]=&p[i*r];
    return(pp);
}
inline void qx_freed_3(double ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline unsigned char** qx_allocu(int r,int c,int padding=QX_DEF_PADDING)
{
	unsigned char *a,**p;
	a=(unsigned char*) malloc(sizeof(unsigned char)*(r*c+padding));
	if(a==NULL) {printf("qx_allocu() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(unsigned char**) malloc(sizeof(unsigned char*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freeu(unsigned char **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline unsigned char *** qx_allocu_3(int n,int r,int c,int padding=QX_DEF_PADDING)
{
	unsigned char *a,**p,***pp;
    int rc=r*c;
    int i,j;
	a=(unsigned char*) malloc(sizeof(unsigned char )*(n*rc+padding));
	if(a==NULL) {printf("qx_allocu_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
    p=(unsigned char**) malloc(sizeof(unsigned char*)*n*r);
    pp=(unsigned char***) malloc(sizeof(unsigned char**)*n);
    for(i=0;i<n;i++) 
        for(j=0;j<r;j++) 
            p[i*r+j]=&a[i*rc+j*c];
    for(i=0;i<n;i++) 
        pp[i]=&p[i*r];
    return(pp);
}
inline void qx_freeu_3(unsigned char ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline float** qx_allocf(int r,int c,int padding=QX_DEF_PADDING)
{
	float *a,**p;
	a=(float*) malloc(sizeof(float)*(r*c+padding));
	if(a==NULL) {printf("qx_allocf() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(float**) malloc(sizeof(float*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freef(float **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline float *** qx_allocf_3(int n,int r,int c,int padding=QX_DEF_PADDING)
{
	float *a,**p,***pp;
    int rc=r*c;
    int i,j;
	a=(float*) malloc(sizeof(float)*(n*rc+padding));
	if(a==NULL) {printf("qx_allocf_3() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
    p=(float**) malloc(sizeof(float*)*n*r);
    pp=(float***) malloc(sizeof(float**)*n);
    for(i=0;i<n;i++) 
        for(j=0;j<r;j++) 
            p[i*r+j]=&a[i*rc+j*c];
    for(i=0;i<n;i++) 
        pp[i]=&p[i*r];
    return(pp);
}
inline void qx_freef_3(float ***p)
{
	if(p!=NULL)
	{
		free(p[0][0]);
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline int** qx_alloci(int r,int c,int padding=QX_DEF_PADDING)
{
	int *a,**p;
	a=(int*) malloc(sizeof(int)*(r*c+padding));
	if(a==NULL) {printf("qx_alloci() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(int**) malloc(sizeof(int*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freei(int **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}
inline double** qx_allocd(int r,int c,int padding=QX_DEF_PADDING)
{
	double *a,**p;
	a=(double*) malloc(sizeof(double)*(r*c+padding));
	if(a==NULL) {printf("qx_allocd() fail, Memory is too huge, fail.\n"); getchar(); exit(0); }
	p=(double**) malloc(sizeof(double*)*r);
	for(int i=0;i<r;i++) p[i]= &a[i*c];
	return(p);
}
inline void qx_freed(double **p)
{
	if(p!=NULL)
	{
		free(p[0]);
		free(p);
		p=NULL;
	}
}

#include "qx_basic.h"
inline void qx_timer::start(){ m_begin = clock(); }
inline float qx_timer::stop(){ m_end = clock(); return (float(m_end - m_begin) / CLOCKS_PER_SEC); }
inline void qx_timer::time_display(char *disp){ printf("Running time (%s) is: %5.5f Seconds.\n", disp, stop()); }
inline void qx_timer::fps_display(char *disp){ printf("Running time (%s) is: %5.5f frame per second.\n", disp, 1.0f / stop()); }
inline void boxcar_sliding_window_x(double *out, double *in, int h, int w, int radius)
{
	double scale = 1.0f / (2 * radius + 1);
	for (int y = 0; y < h; y++) {
		double t;
		// do left edge
		t = in[y*w] * radius;
		for (int x = 0; x < radius + 1; x++) {
			t += in[y*w + x];
		}
		out[y*w] = t * scale;
		for (int x = 1; x < radius + 1; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[y*w];
			out[c] = t * scale;
		}
		// main loop
		for (int x = radius + 1; x < w - radius; x++) {
			int c = y*w + x;
			t += in[c + radius];
			t -= in[c - radius - 1];
			out[c] = t * scale;
		}
		// do right edge
		for (int x = w - radius; x < w; x++) {
			int c = y*w + x;
			t += in[(y*w) + w - 1];
			t -= in[c - radius - 1];
			out[c] = t * scale;
		}

	}
}
inline void boxcar_sliding_window_y(double *out, double *in, int h, int w, int radius)
{
	double scale = 1.0f / (2 * radius + 1);
	for (int x = 0; x < w; x++)
	{
		double t;
		// do left edge
		t = in[x] * radius;
		for (int y = 0; y < radius + 1; y++) {
			t += in[y*w + x];
		}
		out[x] = t * scale;
		for (int y = 1; y < radius + 1; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[x];
			out[c] = t * scale;
		}
		// main loop
		for (int y = radius + 1; y < h - radius; y++) {
			int c = y*w + x;
			t += in[c + radius*w];
			t -= in[c - (radius*w) - w];
			out[c] = t * scale;
		}
		// do right edge
		for (int y = h - radius; y < h; y++) {
			int c = y*w + x;
			t += in[(h - 1)*w + x];
			t -= in[c - (radius*w) - w];
			out[c] = t * scale;
		}
	}
}
inline void boxcar_sliding_window(double **out, double **in, double **temp, int h, int w, int radius)
{
	boxcar_sliding_window_x(temp[0], in[0], h, w, radius);
	boxcar_sliding_window_y(out[0], temp[0], h, w, radius);
}
inline void gaussian_recursive_x(double **od, double **id, int w, int h, double a0, double a1, double a2, double a3, double b1, double b2, double coefp, double coefn)
{
	double xp = 0.0f;  // previous input
	double yp = 0.0f;  // previous output
	double yb = 0.0f;  // previous output by 2
	for (int y = 0; y<h; y++)
	{
		xp = id[y][0]; yb = coefp*xp; yp = yb;
		for (int x = 0; x < w; x++)
		{
			double xc = id[y][x];
			double yc = a0*xc + a1*xp - b1*yp - b2*yb;
			od[y][x] = yc;
			xp = xc; yb = yp; yp = yc;
		}
	}
	// reverse pass
	// ensures response is symmetrical
	double xn = 0.f;
	double xa = 0.f;
	double yn = 0.f;
	double ya = 0.f;
	for (int y = 0; y<h; y++)
	{
		xn = xa = id[y][w - 1]; yn = coefn*xn; ya = yn;
		for (int x = w - 1; x >= 0; x--) {
			double xc = id[y][x];
			double yc = a2*xn + a3*xa - b1*yn - b2*ya;
			xa = xn; xn = xc; ya = yn; yn = yc;
			od[y][x] = od[y][x] + yc;
		}
	}
}
inline void gaussian_recursive_y(double **od, double **id, int w, int h, double a0, double a1, double a2, double a3, double b1, double b2, double coefp, double coefn)
{
	double xp = 0.0f;  // previous input
	double yp = 0.0f;  // previous output
	double yb = 0.0f;  // previous output by 2
	for (int x = 0; x < w; x++)
	{
		xp = id[0][x]; yb = coefp*xp; yp = yb;
		for (int y = 0; y<h; y++)
		{
			double xc = id[y][x];
			double yc = a0*xc + a1*xp - b1*yp - b2*yb;
			od[y][x] = yc;
			xp = xc; yb = yp; yp = yc;
		}
	}


	// reverse pass
	// ensures response is symmetrical
	double xn = 0.f;
	double xa = 0.f;
	double yn = 0.f;
	double ya = 0.f;
	for (int x = 0; x < w; x++)
	{
		xn = xa = id[h - 1][x]; yn = coefn*xn; ya = yn;
		for (int y = h - 1; y >= 0; y--)
		{
			double xc = id[y][x];
			double yc = a2*xn + a3*xa - b1*yn - b2*ya;
			xa = xn; xn = xc; ya = yn; yn = yc;
			od[y][x] = od[y][x] + yc;
		}
	}
}
inline int gaussian_recursive(double **image, double **temp, double sigma, int order, int h, int w)
{
	const double
		nsigma = sigma < 0.1f ? 0.1f : sigma,
		alpha = 1.695f / nsigma,
		ema = exp(-alpha),
		ema2 = exp(-2 * alpha),
		b1 = -2 * ema,
		b2 = ema2;
	double a0 = 0, a1 = 0, a2 = 0, a3 = 0, coefp = 0, coefn = 0;
	switch (order) {
	case 0: {
		const double k = (1 - ema)*(1 - ema) / (1 + 2 * alpha*ema - ema2);
		a0 = k;
		a1 = k*(alpha - 1)*ema;
		a2 = k*(alpha + 1)*ema;
		a3 = -k*ema2;
	} break;

	case 1: {
		const double k = (1 - ema)*(1 - ema) / ema;
		a0 = k*ema;
		a1 = a3 = 0;
		a2 = -a0;
	} break;

	case 2: {
		const double
			ea = exp(-alpha),
			k = -(ema2 - 1) / (2 * alpha*ema),
			kn = (-2 * (-1 + 3 * ea - 3 * ea*ea + ea*ea*ea) / (3 * ea + 1 + 3 * ea*ea + ea*ea*ea));
		a0 = kn;
		a1 = -kn*(1 + k*alpha)*ema;
		a2 = kn*(1 - k*alpha)*ema;
		a3 = -kn*ema2;
	} break;

	default:
		fprintf(stderr, "gaussianFilter: invalid order parameter!\n");
		return 0;
	}
	coefp = (a0 + a1) / (1 + b1 + b2);
	coefn = (a2 + a3) / (1 + b1 + b2);
	//timer.start();
	gaussian_recursive_x(temp, image, w, h, a0, a1, a2, a3, b1, b2, coefp, coefn);
	//image_display(temp,h,w);
	gaussian_recursive_y(image, temp, w, h, a0, a1, a2, a3, b1, b2, coefp, coefn);
	//timer.fps_display();
	return(0);
}
inline void qx_specular_free_image(unsigned char ***image_specular_free, unsigned char ***image_normalized, float **diffuse_chromaticity_max, int h, int w)
{
	int y, x;
	unsigned char *image_specular_free_x, *image_normalized_x; float *diffuse_chromaticity_max_x;
	unsigned char r, g, b; double imax, isum; float rf, gf, bf, c, t0, t1, t2, t3, diffuse, specular;
	//*image_sum_x,*image_max_x,*chromaticity_max_x,
	image_specular_free_x = image_specular_free[0][0];
	image_normalized_x = image_normalized[0][0];
	diffuse_chromaticity_max_x = diffuse_chromaticity_max[0];
	for (y = 0; y<h; y++)
	{
		for (x = 0; x<w; x++)
		{
			t1 = 3.f*(*diffuse_chromaticity_max_x++) - 1.f;
			t3 = 1.0f / 3.0f;
			r = (*image_normalized_x++);
			g = (*image_normalized_x++);
			b = (*image_normalized_x++);
			if (t1>0)
			{
				isum = r + g + b;
				if (isum == 0) c = 0;
				else
				{
					imax = max(max(r, g), b);
					c = (float)(imax / isum);
				}
				t0 = t1*c;
				if (fabs(t0)<QX_DEF_THRESHOLD_ZERO)
				{
					*image_specular_free_x++ = r;
					*image_specular_free_x++ = g;
					*image_specular_free_x++ = b;
				}
				else
				{
					t2 = (3.0f*c - 1.f);
					diffuse = float(imax*t2 / t0);
					specular = float(t3*(isum - diffuse));
					rf = r - specular;
					gf = g - specular;
					bf = b - specular;
					if (rf<0.f) rf = 0.f; else if (rf>255.f) rf = 255.f;
					if (gf<0.f) gf = 0.f; else if (gf>255.f) gf = 255.f;
					if (bf<0.f) bf = 0.f; else if (bf>255.f) bf = 255.f;
					*image_specular_free_x++ = (unsigned char) (rf + 0.5f);
					*image_specular_free_x++ = (unsigned char) (gf + 0.5f);
					*image_specular_free_x++ = (unsigned char) (bf + 0.5f);
				}
			}
			else
			{
				*image_specular_free_x++ = r;
				*image_specular_free_x++ = g;
				*image_specular_free_x++ = b;
			}
		}
	}
}
#endif
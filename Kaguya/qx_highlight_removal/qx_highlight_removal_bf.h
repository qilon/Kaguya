/***************************************************************************************************
\Author:	Qingxiong Yang
\Function:	Highlight removal using bilateral filter.
\reference: Qingxiong Yang, Shengnan Wang and Narendra Ahuja, Real-time Specular Highlight Removal 
			Using Bilateral Filtering, European Conference on Computer Vision (ECCV) 2010.
****************************************************************************************************/
#ifndef QX_HIGHLIGHT_REMOVAL_BF_H
#define QX_HIGHLIGHT_REMOVAL_BF_H
#include "qx_ctbf_ss.h"
#define QX_DEF_DARK_PIXEL								20
#define QX_DEF_THRESHOLD_SIGMA_CHANGE					0.03f
class qx_highlight_removal_bf
{
public:
	qx_highlight_removal_bf();
	~qx_highlight_removal_bf();
	void clean();
	int init(int h,int w,
		unsigned char threshold_dark_pixel=QX_DEF_DARK_PIXEL,
		float threshold_sigma_change=QX_DEF_THRESHOLD_SIGMA_CHANGE);
	int diffuse(unsigned char***image_diffuse,unsigned char***image,int nr_iter=0);
private:
	qx_ctbf_ss m_bf;
	int m_h,m_w; unsigned char m_threshold_dark_pixel; float m_threshold_sigma_change; int m_nr_iteration;
	unsigned char ***m_image_diffuse,**m_image_sf,***m_image_backup,**m_mask,**m_mask_dark_pixel,**m_temp;
	float ***m_max_chrom,**m_max_chrom_backup;
	void compute_approximated_maximum_diffuse_chromaticity(unsigned char **image_approximated_max_diffuse_chromaticity,
		unsigned char ***image_normalized,float**image_max_chrom,unsigned char**mask,unsigned char threshold_dark_pixel,int h,int w);
	void compute_diffuse_reflection_from_maximum_diffuse_chromaticity(unsigned char ***image_approximated_max_diffuse_chromaticity,unsigned char ***image_normalized,
		float **max_diffuse_chromaticity,unsigned char**mask,int h,int w);
};

inline qx_highlight_removal_bf::qx_highlight_removal_bf()
{
	m_image_sf = NULL;
	m_mask = NULL;
	m_mask_dark_pixel = NULL;
	m_image_backup = NULL;
	m_image_diffuse = NULL;
	m_max_chrom = NULL;
	m_max_chrom_backup = NULL;
	m_temp = NULL;
}
inline qx_highlight_removal_bf::~qx_highlight_removal_bf()
{
	clean();
}
inline void qx_highlight_removal_bf::clean()
{
	qx_freeu(m_image_sf); m_image_sf = NULL;
	qx_freeu(m_mask); m_mask = NULL;
	qx_freeu(m_temp); m_temp = NULL;
	qx_freeu(m_mask_dark_pixel); m_mask_dark_pixel = NULL;
	qx_freeu_3(m_image_diffuse); m_image_diffuse = NULL;
	qx_freeu_3(m_image_backup); m_image_backup = NULL;
	qx_freef_3(m_max_chrom); m_max_chrom = NULL;
	qx_freef(m_max_chrom_backup); m_max_chrom_backup = NULL;
}
inline int qx_highlight_removal_bf::init(int h, int w, unsigned char threshold_dark_pixel, float threshold_sigma_change)
{
	m_h = h; m_w = w; m_threshold_dark_pixel = threshold_dark_pixel; m_threshold_sigma_change = threshold_sigma_change;
	m_image_sf = qx_allocu(m_h, m_w);
	m_mask = qx_allocu(m_h, m_w);
	m_temp = qx_allocu(m_h, m_w);
	m_mask_dark_pixel = qx_allocu(m_h, m_w);
	m_image_diffuse = qx_allocu_3(m_h, m_w, 3);
	m_image_backup = qx_allocu_3(m_h, m_w, 3);
	m_max_chrom = qx_allocf_3(2, m_h, m_w);
	m_max_chrom_backup = qx_allocf(m_h, m_w);
	m_bf.init(m_h, m_w);
	return(0);
}
inline int qx_highlight_removal_bf::diffuse(unsigned char***image_diffuse, unsigned char***image, int nr_iter)
{
	float **max_chrom_bf, **max_chrom;
	max_chrom = m_max_chrom[0];
	max_chrom_bf = m_max_chrom[1];
	m_nr_iteration = nr_iter;
	if (m_nr_iteration>0)/*if # of iteration is pre-define*/
	{
		for (int i = 0; i<nr_iter; i++)
		{
			if (i == 0)
			{
				compute_approximated_maximum_diffuse_chromaticity(m_image_sf, image, max_chrom, m_mask_dark_pixel, m_threshold_dark_pixel, m_h, m_w);
			}
			m_bf.joint_bilateral_filter(max_chrom_bf, max_chrom, m_image_sf, m_mask_dark_pixel, 10, 0.25, 10.f / 255);
			for (int y = 0; y<m_h; y++) for (int x = 0; x<m_w; x++) max_chrom[y][x] = max(max_chrom[y][x], max_chrom_bf[y][x]);
		}
	}
	else/*Iterative bilateral filtering until converge*/
	{
		int i = 0;
		while (i >= 0)
		{
			if (i == 0)/*compute our approximation ($\lambda_{max}$) of the maximum diffuse chromaticity ($\Lambda_{max}$)*/
			{
				compute_approximated_maximum_diffuse_chromaticity(m_image_sf, image, max_chrom, m_mask_dark_pixel, m_threshold_dark_pixel, m_h, m_w);
				//image_display(m_image_sf,m_h,m_w);
			}
			//image_display(max_chrom,m_h,m_w);
			m_bf.joint_bilateral_filter(max_chrom_bf, max_chrom, m_image_sf, m_mask_dark_pixel, 10, 0.25, 10.f / 255);/*O(1) bilateral filtering*/
			int counter = 0;
			for (int y = 0; y<m_h; y++) for (int x = 0; x<m_w; x++)
			{
				if (max_chrom_bf[y][x] - max_chrom[y][x]>m_threshold_sigma_change) counter++;
				max_chrom[y][x] = max(max_chrom[y][x], max_chrom_bf[y][x]);
			}
			//image_display(max_chrom,m_h,m_w);
			i++;
			//if(counter>0) printf("[%02d] unconverge: [%010d: %010d]\n",i,counter,m_h*m_w);
			if (counter == 0)
			{
				m_nr_iteration = i;
				i = -1;
			}
		}
	}
	//image_display(image,m_h,m_w);
	compute_diffuse_reflection_from_maximum_diffuse_chromaticity(image_diffuse, image, max_chrom, m_mask_dark_pixel, m_h, m_w);/*compute diffuse reflection*/
	//image_display(image_diffuse,m_h,m_w);
	return(m_nr_iteration);
}
inline void qx_highlight_removal_bf::compute_approximated_maximum_diffuse_chromaticity(unsigned char **image_approximated_max_diffuse_chromaticity,
	unsigned char ***image_normalized, float**image_max_chrom, unsigned char**mask,
	unsigned char threshold_dark_pixel, int h, int w)
{
	int y, x;
	unsigned char *image_approximated_x, *image_normalized_x, *mask_x;
	unsigned char r, g, b; int imax, isum;
	//double t; float rf,gf,bf,c,diffuse,specular;
	float c;//diffuse,specular;
	image_approximated_x = image_approximated_max_diffuse_chromaticity[0];
	image_normalized_x = image_normalized[0][0];
	mask_x = mask[0];
	for (y = 0; y<h; y++)
	{
		for (x = 0; x<w; x++)
		{
			r = (*image_normalized_x++);
			g = (*image_normalized_x++);
			b = (*image_normalized_x++);
			imax = max(max(r, g), b);
			if (imax>threshold_dark_pixel)
			{
				*mask_x++ = 255;
				isum = r + g + b;
				c = (float)imax / isum;
				if (isum>10) image_max_chrom[y][x] = c;
				else image_max_chrom[y][x] = 0;

				float chrom[3] = { (float)r / isum, (float)g / isum, (float)b / isum };
				float cmin = qx_min_f3(chrom);
				float cmax = qx_max_f3(chrom);
				float cappro = (cmax - cmin) / (1 - 3 * cmin);
				*image_approximated_x++ = (unsigned char)(255 * (cappro*1.5 - 0.5) + 0.5);
				//m_mask[y][x]=unsigned char(rf*0.3+gf*0.6+bf*0.1+0.5);

			}
			else
			{
				*mask_x++ = 0;
				*image_approximated_x++ = 0;
			}
		}
	}
}
inline void qx_highlight_removal_bf::compute_diffuse_reflection_from_maximum_diffuse_chromaticity(unsigned char ***image_approximated_max_diffuse_chromaticity, unsigned char ***image_normalized,
	float **max_diffuse_chromaticity, unsigned char**mask, int h, int w)
{
	int y, x;
	unsigned char *image_specular_free_x, *image_normalized_x; float *diffuse_chromaticity_max_x;
	unsigned char r, g, b; double imax, isum; float rf, gf, bf, c, t0, t1, t2, t3, diffuse, specular;
	//*image_sum_x,*image_max_x,*chromaticity_max_x,
	image_specular_free_x = image_approximated_max_diffuse_chromaticity[0][0];
	image_normalized_x = image_normalized[0][0];
	diffuse_chromaticity_max_x = max_diffuse_chromaticity[0];
	t3 = 1.0f / 3.0f;
	for (y = 0; y<h; y++)
	{
		for (x = 0; x<w; x++)
		{
			if (mask[y][x])
			{
				//t1=3.f*(*diffuse_chromaticity_max_x++)-1.f; 
				//r=(*image_normalized_x++);
				//g=(*image_normalized_x++);
				//b=(*image_normalized_x++);
				t1 = 3.f*max_diffuse_chromaticity[y][x] - 1.f;
				r = image_normalized[y][x][0];
				g = image_normalized[y][x][1];
				b = image_normalized[y][x][2];
				if (t1>0)
				{
					isum = r + g + b;
					if (isum == 0) c = 0;
					else
					{
						imax = max(max(r, g), b);
						c = float(imax / isum);
					}
					t0 = t1*c;
					if (fabs(t0)<QX_DEF_THRESHOLD_ZERO)
					{
						//*image_specular_free_x++=r;
						//*image_specular_free_x++=g;
						//*image_specular_free_x++=b;
						image_approximated_max_diffuse_chromaticity[y][x][0] = r;
						image_approximated_max_diffuse_chromaticity[y][x][1] = g;
						image_approximated_max_diffuse_chromaticity[y][x][2] = b;
					}
					else
					{
						t2 = (3.0f*c - 1.f);
						diffuse = float(imax*t2 / t0);
						//specular=max(0,t3*(isum-diffuse));
						specular = float(t3*(isum - diffuse));
						rf = r - specular;
						gf = g - specular;
						bf = b - specular;
						if (rf<0.f) rf = 0.f; else if (rf>255.f) rf = 255.f;
						if (gf<0.f) gf = 0.f; else if (gf>255.f) gf = 255.f;
						if (bf<0.f) bf = 0.f; else if (bf>255.f) bf = 255.f;
						//*image_specular_free_x++=unsigned char(rf+0.5f);
						//*image_specular_free_x++=unsigned char(gf+0.5f);
						//*image_specular_free_x++=unsigned char(bf+0.5f);
						image_approximated_max_diffuse_chromaticity[y][x][0] = (unsigned char)(rf + 0.5f);
						image_approximated_max_diffuse_chromaticity[y][x][1] = (unsigned char)(gf + 0.5f);
						image_approximated_max_diffuse_chromaticity[y][x][2] = (unsigned char)(bf + 0.5f);
					}
				}
			}
			else
			{
				//*image_specular_free_x++=(*image_normalized_x++);
				//*image_specular_free_x++=(*image_normalized_x++);
				//*image_specular_free_x++=(*image_normalized_x++);
				//*diffuse_chromaticity_max_x++;

				for (int c = 0; c<3; c++) image_approximated_max_diffuse_chromaticity[y][x][c] = image_normalized[y][x][c];
			}

		}
	}
}


#endif

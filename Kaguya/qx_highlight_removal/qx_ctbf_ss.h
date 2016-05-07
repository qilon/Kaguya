/*************************************************************************************************
\Author:	Qingxiong Yang
\Function:	This class implements the O(1) bilateral filtering method presented in the reference.
\reference: Qingxiong Yang, Kar-Han Tan and Narendra Ahuja, Real-time O(1) Bilateral Filtering, 
               IEEE Conference on Computer Vision and Pattern Recognition (CVPR) 2009.
**************************************************************************************************/
#ifndef QX_CTBF_SS_H
#define QX_CTBF_SS_H
#define QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER		0
#define QX_DEF_CTBF_BOX_BILATERAL_FILTER			1
#define QX_DEF_CTBF_INTENSITY_RANGE					256
#define QX_DEF_CTBF_SIGMA_SPATIAL_DEFAULT			0.03
#define QX_DEF_CTBF_SIGMA_RANGE_DEFAULT				0.08
#define QX_DEF_CTBF_H_MIN							24
#define QX_DEF_CTBF_W_MIN							32
class qx_ctbf_ss
{
public:
    qx_ctbf_ss();
    ~qx_ctbf_ss();
    void clean();
	int init(int h_original,int w_original,
		int spatial_filter=QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER,//default: Gaussian BF
		double sigma_spatial=QX_DEF_CTBF_SIGMA_SPATIAL_DEFAULT,//0.03~16/512
		double sigma_range=QX_DEF_CTBF_SIGMA_RANGE_DEFAULT);//0.08~20/255
	int joint_bilateral_filter(unsigned char**image_filtered,unsigned char**image,unsigned char**texture,unsigned char**mask,int nr_scale,double sigma_spatial,double sigma_range);
	int joint_bilateral_filter(float**image_filtered,float**image,unsigned char**texture,unsigned char**mask,int nr_scale=8,double sigma_spatial=0,double sigma_range=0);
private:
	char m_str[300];
	int m_h,m_w,m_h_original,m_w_original,m_nr_shift,m_radius; double m_sigma_range,m_sigma_spatial; int m_nr_scale; int m_spatial_filter;
	double ***m_jk,**m_wk; double **m_box; 
	double *m_grayscale;
	double *m_table;
	unsigned char **m_image_y,**m_image_y_downsampled,**m_image_y_downsampled_texture;
	float**m_image_y_downsampled_f;
	void get_down_sampled_y_component(unsigned char **out,unsigned char ***in,int h,int w,int scale_exp);
	void rgb_2_yuv(unsigned char **out_y,unsigned char **out_u,unsigned char **out_v,unsigned char ***in,int h,int w);
	//void yuv_2_rgb_upsample(unsigned char ***out,unsigned char **in_y,unsigned char **in_u,unsigned char **in_v,int h,int w);
	void yuv_2_rgb(unsigned char ***out,unsigned char **in_y,unsigned char ***in_rgb,int h,int w);
};

#include "qx_basic.h"
//#include "qx_ppm.h"
#include "qx_ctbf_ss.h"
inline qx_ctbf_ss::qx_ctbf_ss()
{
	m_box = NULL;
	m_jk = NULL;
	m_wk = NULL;
	m_grayscale = NULL;
	m_image_y_downsampled = NULL;
	m_image_y_downsampled_f = NULL;
	m_image_y_downsampled_texture = NULL;
}
inline qx_ctbf_ss::~qx_ctbf_ss()
{
	clean();
}
inline void qx_ctbf_ss::clean()
{
	qx_freed(m_box); m_box = NULL;
	qx_freed_3(m_jk); m_jk = NULL;
	qx_freed(m_wk); m_wk = NULL;
	if (m_grayscale != NULL) delete[] m_grayscale; m_grayscale = NULL;
	
	if (m_table != NULL) delete[] m_table; m_table = NULL;

	qx_freeu(m_image_y_downsampled); m_image_y_downsampled = NULL;
	qx_freef(m_image_y_downsampled_f); m_image_y_downsampled_f = NULL;
	qx_freeu(m_image_y_downsampled_texture); m_image_y_downsampled_texture = NULL;
}
inline int qx_ctbf_ss::init(int h_original, int w_original, int spatial_filter, double sigma_spatial, double sigma_range)
{
	m_h_original = h_original; m_w_original = w_original;
	int shift_h = int(log((double)m_h_original / QX_DEF_CTBF_H_MIN) / log(2.0) + 0.5);
	int shift_w = int(log((double)m_w_original / QX_DEF_CTBF_W_MIN) / log(2.0) + 0.5);
	m_nr_shift = max(1, min(shift_h, shift_w));
	printf("m_nr_shift: [%d]\n", m_nr_shift);
	m_h = (h_original >> m_nr_shift); m_w = (w_original >> m_nr_shift);
	m_spatial_filter = spatial_filter; m_sigma_spatial = sigma_spatial; m_sigma_range = sigma_range; m_nr_scale = QX_DEF_CTBF_INTENSITY_RANGE;
	if (m_spatial_filter == QX_DEF_CTBF_BOX_BILATERAL_FILTER) m_radius = int(m_sigma_spatial*min(m_h, m_w) + 0.5);
	else if (m_spatial_filter != QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER)
	{
		printf("Note: ONLY support box and Gaussian spatial filter!");
		printf("Switching to Gaussian spatial filter automatically.....");
		m_spatial_filter = QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER;
	}
	/*memory allocation*/
	m_box = qx_allocd(m_h, m_w);
	m_jk = qx_allocd_3(2, m_h, m_w);
	m_wk = qx_allocd(m_h, m_w);
	m_table = get_color_weighted_table(m_sigma_range*QX_DEF_CTBF_INTENSITY_RANGE, QX_DEF_CTBF_INTENSITY_RANGE);
	m_grayscale = new double[m_nr_scale];

	m_image_y_downsampled = qx_allocu(m_h, m_w);
	m_image_y_downsampled_f = qx_allocf(m_h, m_w);
	m_image_y_downsampled_texture = qx_allocu(m_h, m_w);
	return(0);
}

inline int qx_ctbf_ss::joint_bilateral_filter(unsigned char**image_filtered, unsigned char**image, unsigned char**texture, unsigned char**mask, int nr_scale,
	double sigma_spatial, double sigma_range)
{
	unsigned char image_min, image_max;
	int y, x, jk_0, jk_1;
	if (sigma_range>QX_DEF_THRESHOLD_ZERO)
	{
		m_sigma_range = sigma_range;
		color_weighted_table_update(m_table, m_sigma_range*QX_DEF_CTBF_INTENSITY_RANGE, QX_DEF_CTBF_INTENSITY_RANGE);
	}
	if (sigma_spatial>QX_DEF_THRESHOLD_ZERO) m_sigma_spatial = sigma_spatial;
	//image_display(image,m_h_original,m_w_original);
	down_sample_1(m_image_y_downsampled, image, m_h_original, m_w_original, m_nr_shift);
	down_sample_1(m_image_y_downsampled_texture, texture, m_h_original, m_w_original, m_nr_shift);

	//image_display(m_image_y_downsampled,m_h,m_w);
	vec_min_val(image_min, texture[0], m_h_original*m_w_original);
	vec_max_val(image_max, texture[0], m_h_original*m_w_original);
	m_nr_scale = nr_scale;
	//m_nr_scale=max(2,int(double(image_max-image_min)/(255*m_sigma_range)+0.5));
	//printf("[max,min]:[%5.5f,%5.5f]\n",(float)image_max,(float)image_min);
	//printf("[sigma_range: %1.3f]\n",m_sigma_range);
	//printf("[nr_scale: %d]\n",m_nr_scale);
	m_grayscale[0] = (double)image_min;
	m_grayscale[m_nr_scale - 1] = (double)image_max;
	double delta_scale = double(m_nr_scale - 1) / (image_max - image_min);
	double shift_inv = 1.0 / (1 << m_nr_shift);
	for (int i = 1; i<m_nr_scale - 1; i++) m_grayscale[i] = (double)image_min + i*double(image_max - image_min) / (m_nr_scale - 1);
	for (int i = 0; i<m_nr_scale; i++)
	{
		double **jk;
		if (i == 0)
		{
			jk_0 = 0;
			jk_1 = 1;
			jk = m_jk[jk_0];
		}
		else jk = m_jk[jk_1];
		for (y = 0; y<m_h; y++)
		{
			for (x = 0; x<m_w; x++)
			{
				int index = int(abs(m_grayscale[i] - m_image_y_downsampled_texture[y][x]) + 0.5f);
				jk[y][x] = m_table[index] * m_image_y_downsampled[y][x];
				m_wk[y][x] = m_table[index];
			}
		}
		if (m_spatial_filter == QX_DEF_CTBF_BOX_BILATERAL_FILTER)
		{
			boxcar_sliding_window(jk, jk, m_box, m_h, m_w, m_radius);
			boxcar_sliding_window(m_wk, m_wk, m_box, m_h, m_w, m_radius);
		}
		else if (m_spatial_filter == QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER)
		{
			gaussian_recursive(jk, m_box, m_sigma_spatial*min(m_h, m_w), 0, m_h, m_w);
			gaussian_recursive(m_wk, m_box, m_sigma_spatial*min(m_h, m_w), 0, m_h, m_w);
		}
		for (y = 0; y<m_h; y++)
		{
			for (x = 0; x<m_w; x++)
			{
				jk[y][x] /= m_wk[y][x];
			}
		}
		//image_display(jk,m_h,m_w);
		unsigned char *image_y = texture[0];
		if (i>0)
		{
			for (y = 0; y<m_h_original; y++)
			{
				for (x = 0; x<m_w_original; x++)
				{
					if (mask[y][x])
					{
						double kf;
						kf = double((*image_y++) - image_min)*delta_scale;
						int k = int(kf);
						double ys = min(y*shift_inv, double(m_h - 1));
						double xs = min(x*shift_inv, double(m_w - 1));
						if (k == (i - 1))
						{
							double alpha = (k + 1) - kf;
							image_filtered[y][x] = (unsigned char) (min(255.0, alpha*qx_linear_interpolate_xy(m_jk[jk_0], xs, ys, m_h, m_w) + (1.f - alpha)*qx_linear_interpolate_xy(m_jk[jk_1], xs, ys, m_h, m_w)));
						}
						else if (k == i&&i == (m_nr_scale - 1))
						{
							image_filtered[y][x] = (unsigned char) (min(255.0, qx_linear_interpolate_xy(m_jk[jk_1], xs, ys, m_h, m_w) + 0.5f));
						}
					}
					else
					{
						*image_y++;
					}
				}
			}
			jk_1 = jk_0;
			jk_0 = (jk_0 + 1) % 2;
		}
	}
	return(0);
}

inline int qx_ctbf_ss::joint_bilateral_filter(float**image_filtered, float**image, unsigned char **texture, unsigned char**mask, int nr_scale,
	double sigma_spatial, double sigma_range)
{
	unsigned char image_min, image_max;
	int y, x, jk_0, jk_1;
	if (sigma_range>QX_DEF_THRESHOLD_ZERO)
	{
		m_sigma_range = sigma_range;
		color_weighted_table_update(m_table, m_sigma_range*QX_DEF_CTBF_INTENSITY_RANGE, QX_DEF_CTBF_INTENSITY_RANGE);
	}
	if (sigma_spatial>QX_DEF_THRESHOLD_ZERO) m_sigma_spatial = sigma_spatial;
	//image_display(image,m_h_original,m_w_original);
	down_sample_1(m_image_y_downsampled_f, image, m_h_original, m_w_original, m_nr_shift);
	down_sample_1(m_image_y_downsampled_texture, texture, m_h_original, m_w_original, m_nr_shift);

	//image_display(m_image_y_downsampled,m_h,m_w);
	vec_min_val(image_min, texture[0], m_h_original*m_w_original);
	vec_max_val(image_max, texture[0], m_h_original*m_w_original);
	m_nr_scale = nr_scale;
	//m_nr_scale=max(2,int(double(image_max-image_min)/(255*m_sigma_range)+0.5));
	//printf("[max,min]:[%5.5f,%5.5f]\n",(float)image_max,(float)image_min);
	//printf("[sigma_range: %1.3f]\n",m_sigma_range);
	//printf("[nr_scale: %d]\n",m_nr_scale);
	m_grayscale[0] = (double)image_min;
	m_grayscale[m_nr_scale - 1] = min(255, (int)image_max + 1);
	double delta_scale = double(m_nr_scale - 1) / (image_max - image_min);
	double shift_inv = 1.0 / (1 << m_nr_shift);
	for (int i = 1; i<m_nr_scale - 1; i++) m_grayscale[i] = (double)image_min + i*double(image_max - image_min) / (m_nr_scale - 1);
	for (int i = 0; i<m_nr_scale; i++)
	{
		double **jk;
		if (i == 0)
		{
			jk_0 = 0;
			jk_1 = 1;
			jk = m_jk[jk_0];
		}
		else jk = m_jk[jk_1];
		for (y = 0; y<m_h; y++)
		{
			for (x = 0; x<m_w; x++)
			{
				int index = int(abs(m_grayscale[i] - m_image_y_downsampled_texture[y][x]) + 0.5f);
				jk[y][x] = m_table[index] * m_image_y_downsampled_f[y][x];
				m_wk[y][x] = m_table[index];
			}
		}
		if (m_spatial_filter == QX_DEF_CTBF_BOX_BILATERAL_FILTER)
		{
			boxcar_sliding_window(jk, jk, m_box, m_h, m_w, m_radius);
			boxcar_sliding_window(m_wk, m_wk, m_box, m_h, m_w, m_radius);
		}
		else if (m_spatial_filter == QX_DEF_CTBF_GAUSSIAN_BILATERAL_FILTER)
		{
			gaussian_recursive(jk, m_box, m_sigma_spatial*min(m_h, m_w), 0, m_h, m_w);
			gaussian_recursive(m_wk, m_box, m_sigma_spatial*min(m_h, m_w), 0, m_h, m_w);
		}
		for (y = 0; y<m_h; y++)
		{
			for (x = 0; x<m_w; x++)
			{
				jk[y][x] /= m_wk[y][x];
			}
		}
		//image_display(jk,m_h,m_w);
		unsigned char *image_y = texture[0];
		if (i>0)
		{
			for (y = 0; y<m_h_original; y++)
			{
				for (x = 0; x<m_w_original; x++)
				{

					if (mask[y][x])
					{
						double kf;
						kf = double((*image_y++) - image_min)*delta_scale;
						int k = int(kf);
						double ys = min(y*shift_inv, (double)m_h - 1);
						double xs = min(x*shift_inv, (double)m_w - 1);
						if (k == (i - 1))
						{
							double alpha = (k + 1) - kf;
							image_filtered[y][x] = float(alpha*qx_linear_interpolate_xy(m_jk[jk_0], xs, ys, m_h, m_w) + (1.f - alpha)*qx_linear_interpolate_xy(m_jk[jk_1], xs, ys, m_h, m_w));
						}
						else if (k == i&&i == (m_nr_scale - 1))
						{
							image_filtered[y][x] = (float)qx_linear_interpolate_xy(m_jk[jk_1], xs, ys, m_h, m_w);
						}
					}
					else *image_y++;
				}
			}
			jk_1 = jk_0;
			jk_0 = (jk_0 + 1) % 2;
		}
	}
	return(0);
}

#endif
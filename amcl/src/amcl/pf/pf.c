/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
/**************************************************************************
 * Desc: Simple particle filter for localization.
 * Author: Andrew Howard
 * Date: 10 Dec 2002
 * CVS: $Id: pf.c 6345 2008-04-17 01:36:39Z gerkey $
 *************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "amcl/pf/pf.h"
#include "amcl/pf/pf_pdf.h"
#include "amcl/pf/pf_kdtree.h"
#include "portable_utils.hpp"

// Compute the required number of samples, given that there are k bins
// with samples in them.
static int pf_resample_limit(pf_t *pf, int k);

// Create a new filter
pf_t *pf_alloc(int min_samples, int max_samples,
			   double alpha_slow, double alpha_fast,					  // 前四个是对滤波器的配置，约束了最少和最多的样本数量，权重过滤值衰减速率
			   pf_init_model_fn_t random_pose_fn, void *random_pose_data) // 参数random_pose_fn是一个函数指针，提供了生成随机样本的函数
																		  //  对于参数random_pose_data，我的理解是粒子的样本空间， amcl实际使用的是地图对象
{
	int i, j;
	pf_t *pf;
	pf_sample_set_t *set; // set: 指向某个粒子集合          set是用来遍历2个粒子集的
	pf_sample_t *sample;  // sample: 指向某个粒子           sample是用来遍历粒子集中的每个粒子的

	srand48(time(NULL));

	pf = calloc(1, sizeof(pf_t)); 
	// calloc 初始化为0   ； malloc 只申请内存，不初始化  所以 calloc 的执行会比 malloc 稍微费时，因为它多了初始化的步骤
								  // pf 分配内存
	pf->random_pose_fn = random_pose_fn;
	pf->random_pose_data = random_pose_data;

	pf->min_samples = min_samples;
	pf->max_samples = max_samples; // 设置最大粒子数量

	// Control parameters for the population size calculation.  [err] is
	// the max error between the true distribution and the estimated
	// distribution.  [z] is the upper standard normal quantile for (1 -
	// p), where p is the probability that the error on the estimated
	// distrubition will be less than [err].
	pf->pop_err = 0.01;
	pf->pop_z = 3;
	pf->dist_threshold = 0.5;

	// Number of leaf nodes is never higher than the max number of samples
	pf->limit_cache = calloc(max_samples, sizeof(int));

	pf->current_set = 0;	// 表示当前激活的粒子集；共有两个粒子集，=0，表示当前激活第一个
	for (j = 0; j < 2; j++) // 遍历2个粒子集
	{
		set = pf->sets + j; // pf->sets 是pf_samples_set_t的数组，set是一个pf_samples_set_t的指针  ，此处set依次获得指向两个数组元素的指针
							// set 是遍历的2个粒子集 pf_sample_set_t,  set指针，分别逐个指向操作
		set->sample_count = max_samples;
		set->samples = calloc(max_samples, sizeof(pf_sample_t)); // 根据参数传入的最大粒子数, 为每个set中的samples指针，分配内存
																 // 内存中包含【max_samples】个粒子对象
																 // 一个 set->samples 的内存占用： 32*max_samples 字节
																 // 设置每个粒子集的：粒子数，所有粒子的初始化；
		for (i = 0; i < set->sample_count; i++)		// 每个set中的例子数=max_samples
		{
			sample = set->samples + i;
			sample->pose.v[0] = 0.0;
			sample->pose.v[1] = 0.0;
			sample->pose.v[2] = 0.0;
			sample->weight = 1.0 / max_samples; 		// 为粒子集中的每个粒子设置位姿(0,0,0) x=0,y=0,theta=0，  权重= 1/sample_count
		}

		// HACK: is 3 times max_samples enough?
		set->kdtree = pf_kdtree_alloc(3 * max_samples); 
				// 构建了三倍于样本集合尺寸的直方图对象kdtree

		set->cluster_count = 0;
		set->cluster_max_count = max_samples; // 最大粒子簇的数量，同粒子数相等；
		set->clusters = calloc(set->cluster_max_count, sizeof(pf_cluster_t));
		// 一个 set->clusters 的内存占用：176*max_samples 字节
		//一个pf滤波器有2个set，粗略计算内存占用，以5000个粒子为例：5000*(32+176)*2/1024/1024 = 1.98 m,         2M的内存，对于嵌入式来说可能很大吧
		set->mean = pf_vector_zero();		
		set->cov = pf_matrix_zero();			// 当前粒子集合的 位姿均值与协方差 都设置为 零向量和零矩阵
	}

	pf->w_slow = 0.0;
	pf->w_fast = 0.0;

	pf->alpha_slow = alpha_slow;
	pf->alpha_fast = alpha_fast;

	// set converged to 0
	pf_init_converged(pf); // 设置pf-> converged = 0;  pf->sets[current_set].converged = 0；设置为不收敛

	return pf;
}

// Free an existing filter
void pf_free(pf_t *pf)
{
	int i;

	free(pf->limit_cache); // free pf对象中申请的内存

	for (i = 0; i < 2; i++)
	{
		free(pf->sets[i].clusters);
		pf_kdtree_free(pf->sets[i].kdtree);
		free(pf->sets[i].samples);
	}
	free(pf); // 释放pf滤波器

	return;
}

// Initialize the filter using a gaussian
void pf_init(pf_t *pf, pf_vector_t mean, pf_matrix_t cov)
// 该函数有三个参数，其中pf是将要初始化的滤波器对象，mean和cov分别是机器人初始位姿和协方差描述
{
	int i;
	pf_sample_set_t *set;
	pf_sample_t *sample;
	pf_pdf_gaussian_t *pdf;

	set = pf->sets + pf->current_set; // 得到激活的 粒子set

	// Create the kd tree for adaptive sampling
	pf_kdtree_clear(set->kdtree);

	set->sample_count = pf->max_samples; // pf_alloc中已经操作过，这里有对当前set重新设置了一下

	pdf = pf_pdf_gaussian_alloc(mean, cov);		// 根据参数mean,cov 创建一个pdf对象；

	// Compute the new sample poses
	for (i = 0; i < set->sample_count; i++)
	{
		sample = set->samples + i;
		sample->weight = 1.0 / pf->max_samples;		// 设置current set中每个粒子的权重，pf_alloc中也操作过
		sample->pose = pf_pdf_gaussian_sample(pdf); // 具体细节没有看懂，功能是：位姿是服从给定均值mean的，按照高斯分布生成的位姿，方差是通过协方差矩阵计算出来的

		// Add sample to histogram
		pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
	}

	pf->w_slow = pf->w_fast = 0.0;

	pf_pdf_gaussian_free(pdf);

	// Re-compute cluster statistics
	pf_cluster_stats(pf, set);

	// set converged to 0
	pf_init_converged(pf);

	return;
}

// Initialize the filter using some model
void pf_init_model(pf_t *pf, pf_init_model_fn_t init_fn, void *init_data)
{
	int i;
	pf_sample_set_t *set;
	pf_sample_t *sample;

	set = pf->sets + pf->current_set;

	// Create the kd tree for adaptive sampling
	pf_kdtree_clear(set->kdtree);

	set->sample_count = pf->max_samples;

	// Compute the new sample poses
	for (i = 0; i < set->sample_count; i++)
	{
		sample = set->samples + i;
		sample->weight = 1.0 / pf->max_samples;
		sample->pose = (*init_fn)(init_data);

		// Add sample to histogram
		pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
	}

	pf->w_slow = pf->w_fast = 0.0;

	// Re-compute cluster statistics
	pf_cluster_stats(pf, set);

	// set converged to 0
	pf_init_converged(pf);

	return;
}

void pf_init_converged(pf_t *pf)
{
	pf_sample_set_t *set;
	set = pf->sets + pf->current_set; // pf->current_set 初始化后默认=0，表示在第一个粒子集； pf->sets表示的是粒子集数组sets[2]的首地址；
									  // 执行pf->sets + i, 就表示获得第i个粒子集的首地址；
	set->converged = 0;				  // 将当前激活的粒子集的   converged = 0
	pf->converged = 0;				  // 将粒子滤波器的 converged 设置为0；—— 不收敛
}

int pf_update_converged(pf_t *pf)
{
	int i;
	pf_sample_set_t *set;
	pf_sample_t *sample;
	double total;

	set = pf->sets + pf->current_set;
	double mean_x = 0, mean_y = 0;

	for (i = 0; i < set->sample_count; i++)
	{
		sample = set->samples + i;

		mean_x += sample->pose.v[0];
		mean_y += sample->pose.v[1];
	}
	mean_x /= set->sample_count;
	mean_y /= set->sample_count;

	for (i = 0; i < set->sample_count; i++)
	{
		sample = set->samples + i;
		if (fabs(sample->pose.v[0] - mean_x) > pf->dist_threshold ||
			fabs(sample->pose.v[1] - mean_y) > pf->dist_threshold)
		{
			set->converged = 0;
			pf->converged = 0;
			return 0;
		}
	}
	set->converged = 1;
	pf->converged = 1;
	return 1;
}

// Update the filter with some new action
void pf_update_action(pf_t *pf, pf_action_model_fn_t action_fn, void *action_data)
{
	pf_sample_set_t *set;

	set = pf->sets + pf->current_set;

	(*action_fn)(action_data, set);

	return;
}

#include <float.h>
// Update the filter with some new sensor observation
void pf_update_sensor(pf_t *pf, pf_sensor_model_fn_t sensor_fn, void *sensor_data)
{
	int i;
	pf_sample_set_t *set;
	pf_sample_t *sample;
	double total;

	set = pf->sets + pf->current_set;

	// Compute the sample weights
	total = (*sensor_fn)(sensor_data, set);

	set->n_effective = 0;

	if (total > 0.0)
	{
		// Normalize weights
		double w_avg = 0.0;
		for (i = 0; i < set->sample_count; i++) // 遍历所有粒子
		{
			sample = set->samples + i;							 // 操作的当前粒子保存在sample指针中
			w_avg += sample->weight;							 // 权重求和
			sample->weight /= total;							 // ??? 当前粒子的权重修改为：权重/total ???
			set->n_effective += sample->weight * sample->weight; // 权重的平方和
		}
		// Update running averages of likelihood of samples (Prob Rob p258)
		w_avg /= set->sample_count; // 所有粒子的平均权重
		if (pf->w_slow == 0.0)
			pf->w_slow = w_avg;
		else
			pf->w_slow += pf->alpha_slow * (w_avg - pf->w_slow);
		if (pf->w_fast == 0.0)
			pf->w_fast = w_avg;
		else
			pf->w_fast += pf->alpha_fast * (w_avg - pf->w_fast);
		// printf("w_avg: %e slow: %e fast: %e\n",
		// w_avg, pf->w_slow, pf->w_fast);
	}
	else
	{
		// Handle zero total
		for (i = 0; i < set->sample_count; i++)
		{
			sample = set->samples + i;
			sample->weight = 1.0 / set->sample_count;
		}
	}

	set->n_effective = 1.0 / set->n_effective;
	return;
}

// copy set a to set b
void copy_set(pf_sample_set_t *set_a, pf_sample_set_t *set_b)
{
	int i;
	double total;
	pf_sample_t *sample_a, *sample_b;

	// Clean set b's kdtree
	pf_kdtree_clear(set_b->kdtree);

	// Copy samples from set a to create set b
	total = 0;
	set_b->sample_count = 0;

	for (i = 0; i < set_a->sample_count; i++)
	{
		sample_b = set_b->samples + set_b->sample_count++;

		sample_a = set_a->samples + i;

		assert(sample_a->weight > 0);

		// Copy sample a to sample b
		sample_b->pose = sample_a->pose;
		sample_b->weight = sample_a->weight;

		total += sample_b->weight;

		// Add sample to histogram          // 把该粒子插入到直方图中
		pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);
	}

	// Normalize weights
	for (i = 0; i < set_b->sample_count; i++)
	{
		sample_b = set_b->samples + i;
		sample_b->weight /= total;
	}

	set_b->converged = set_a->converged;
}

// Resample the distribution
void pf_update_resample(pf_t *pf)
{
	int i;
	double total;
	pf_sample_set_t *set_a, *set_b;
	pf_sample_t *sample_a, *sample_b;

	// double r,c,U;
	// int m;
	// double count_inv;
	double *c;

	double w_diff;

	set_a = pf->sets + pf->current_set;
	set_b = pf->sets + (pf->current_set + 1) % 2;

	if (pf->selective_resampling != 0)
	{
		if (set_a->n_effective > 0.5 * (set_a->sample_count))
		{
			// copy set a to b
			copy_set(set_a, set_b);

			// Re-compute cluster statistics
			pf_cluster_stats(pf, set_b);

			// Use the newly created sample set
			pf->current_set = (pf->current_set + 1) % 2;
			return;
		}
	}

	// Build up cumulative probability table for resampling.
	// TODO: Replace this with a more efficient procedure
	// (e.g., http://www.network-theory.co.uk/docs/gslref/GeneralDiscreteDistributions.html)
	c = (double *)malloc(sizeof(double) * (set_a->sample_count + 1));
	c[0] = 0.0;
	for (i = 0; i < set_a->sample_count; i++)
		c[i + 1] = c[i] + set_a->samples[i].weight;

	// Create the kd tree for adaptive sampling
	pf_kdtree_clear(set_b->kdtree);

	// Draw samples from set a to create set b.
	total = 0;
	set_b->sample_count = 0;

	w_diff = 1.0 - pf->w_fast / pf->w_slow;
	if (w_diff < 0.0)
		w_diff = 0.0;
	// printf("w_diff: %9.6f\n", w_diff);

	// Can't (easily) combine low-variance sampler with KLD adaptive
	// sampling, so we'll take the more traditional route.
	/*
	// Low-variance resampler, taken from Probabilistic Robotics, p110
	count_inv = 1.0/set_a->sample_count;
	r = drand48() * count_inv;
	c = set_a->samples[0].weight;
	i = 0;
	m = 0;
	*/
	while (set_b->sample_count < pf->max_samples)
	{
		sample_b = set_b->samples + set_b->sample_count++;

		if (drand48() < w_diff)
			sample_b->pose = (pf->random_pose_fn)(pf->random_pose_data);
		else
		{
			// Can't (easily) combine low-variance sampler with KLD adaptive
			// sampling, so we'll take the more traditional route.
			/*
			// Low-variance resampler, taken from Probabilistic Robotics, p110
			U = r + m * count_inv;
			while(U>c)
			{
			  i++;
			  // Handle wrap-around by resetting counters and picking a new random
			  // number
			  if(i >= set_a->sample_count)
			  {
				r = drand48() * count_inv;
				c = set_a->samples[0].weight;
				i = 0;
				m = 0;
				U = r + m * count_inv;
				continue;
			  }
			  c += set_a->samples[i].weight;
			}
			m++;
			*/

			// Naive discrete event sampler
			double r;
			r = drand48();
			for (i = 0; i < set_a->sample_count; i++)
			{
				if ((c[i] <= r) && (r < c[i + 1]))
					break;
			}
			assert(i < set_a->sample_count);

			sample_a = set_a->samples + i;

			assert(sample_a->weight > 0);

			// Add sample to list
			sample_b->pose = sample_a->pose;
		}

		sample_b->weight = 1.0;
		total += sample_b->weight;

		// Add sample to histogram
		pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);

		// See if we have enough samples yet
		if (set_b->sample_count > pf_resample_limit(pf, set_b->kdtree->leaf_count))
			break;
	}

	// Reset averages, to avoid spiraling off into complete randomness.
	if (w_diff > 0.0)
		pf->w_slow = pf->w_fast = 0.0;

	// fprintf(stderr, "\n\n");

	// Normalize weights
	for (i = 0; i < set_b->sample_count; i++)
	{
		sample_b = set_b->samples + i;
		sample_b->weight /= total;
	}

	// Re-compute cluster statistics
	pf_cluster_stats(pf, set_b);

	// Use the newly created sample set
	pf->current_set = (pf->current_set + 1) % 2;

	pf_update_converged(pf);

	free(c);
	return;
}

// Compute the required number of samples, given that there are k bins
// with samples in them.  This is taken directly from Fox et al.
int pf_resample_limit(pf_t *pf, int k)
{
	double a, b, c, x;
	int n;

	// Return max_samples in case k is outside expected range, this shouldn't
	// happen, but is added to prevent any runtime errors
	if (k < 1 || k > pf->max_samples)
		return pf->max_samples;

	// Return value if cache is valid, which means value is non-zero positive
	if (pf->limit_cache[k - 1] > 0)
		return pf->limit_cache[k - 1];

	if (k == 1)
	{
		pf->limit_cache[k - 1] = pf->max_samples;
		return pf->max_samples;
	}

	a = 1;
	b = 2 / (9 * ((double)k - 1));
	c = sqrt(2 / (9 * ((double)k - 1))) * pf->pop_z;
	x = a - b + c;

	n = (int)ceil((k - 1) / (2 * pf->pop_err) * x * x * x);

	if (n < pf->min_samples)
	{
		pf->limit_cache[k - 1] = pf->min_samples;
		return pf->min_samples;
	}
	if (n > pf->max_samples)
	{
		pf->limit_cache[k - 1] = pf->max_samples;
		return pf->max_samples;
	}

	pf->limit_cache[k - 1] = n;
	return n;
}

// Re-compute the cluster statistics for a sample set
void pf_cluster_stats(pf_t *pf, pf_sample_set_t *set)
{
	int i, j, k, cidx;
	pf_sample_t *sample;
	pf_cluster_t *cluster;

	// Workspace
	double m[4], c[2][2];		// 这是聚类的统计结果存储位置？ cluster statistics???
	size_t count;
	double weight;

	// Cluster the samples
	pf_kdtree_cluster(set->kdtree);

	// Initialize cluster stats
	set->cluster_count = 0;

	for (i = 0; i < set->cluster_max_count; i++)
	{
		cluster = set->clusters + i;
		cluster->count = 0;
		cluster->weight = 0;
		cluster->mean = pf_vector_zero();
		cluster->cov = pf_matrix_zero();

		for (j = 0; j < 4; j++)
			cluster->m[j] = 0.0; // cluster中的 m[4] ,c[2][2] 一个4维向量，和2x2矩阵的目的是什么？
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
				cluster->c[j][k] = 0.0;
	}

	// Initialize overall filter stats
	count = 0;
	weight = 0.0;
	set->mean = pf_vector_zero();
	set->cov = pf_matrix_zero();
	for (j = 0; j < 4; j++)
		m[j] = 0.0;
	for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
			c[j][k] = 0.0;

	// Compute cluster stats
	for (i = 0; i < set->sample_count; i++)
	{
		sample = set->samples + i;

		// printf("%d %f %f %f\n", i, sample->pose.v[0], sample->pose.v[1], sample->pose.v[2]);

		// Get the cluster label for this sample
		cidx = pf_kdtree_get_cluster(set->kdtree, sample->pose);
		assert(cidx >= 0);
		if (cidx >= set->cluster_max_count)
			continue;
		if (cidx + 1 > set->cluster_count)
			set->cluster_count = cidx + 1;

		cluster = set->clusters + cidx;

		cluster->count += 1;
		cluster->weight += sample->weight;

		count += 1;
		weight += sample->weight;

		// Compute mean
		cluster->m[0] += sample->weight * sample->pose.v[0];
		cluster->m[1] += sample->weight * sample->pose.v[1];
		cluster->m[2] += sample->weight * cos(sample->pose.v[2]);
		cluster->m[3] += sample->weight * sin(sample->pose.v[2]);

		m[0] += sample->weight * sample->pose.v[0];
		m[1] += sample->weight * sample->pose.v[1];
		m[2] += sample->weight * cos(sample->pose.v[2]);
		m[3] += sample->weight * sin(sample->pose.v[2]);

		// Compute covariance in linear components
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
			{
				cluster->c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
				c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
			}
	}

	// Normalize
	for (i = 0; i < set->cluster_count; i++)
	{
		cluster = set->clusters + i;

		cluster->mean.v[0] = cluster->m[0] / cluster->weight;
		cluster->mean.v[1] = cluster->m[1] / cluster->weight;
		cluster->mean.v[2] = atan2(cluster->m[3], cluster->m[2]);

		cluster->cov = pf_matrix_zero();

		// Covariance in linear components
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
				cluster->cov.m[j][k] = cluster->c[j][k] / cluster->weight -
									   cluster->mean.v[j] * cluster->mean.v[k];

		// Covariance in angular components; I think this is the correct
		// formula for circular statistics.
		cluster->cov.m[2][2] = -2 * log(sqrt(cluster->m[2] * cluster->m[2] +
											 cluster->m[3] * cluster->m[3]) /
										cluster->weight);

		// printf("cluster %d %d %f (%f %f %f)\n", i, cluster->count, cluster->weight,
		// cluster->mean.v[0], cluster->mean.v[1], cluster->mean.v[2]);
		// pf_matrix_fprintf(cluster->cov, stdout, "%e");
	}

	assert(fabs(weight) >= DBL_EPSILON);
	if (fabs(weight) < DBL_EPSILON)
	{
		printf("ERROR : divide-by-zero exception : weight is zero\n");
		return;
	}
	// Compute overall filter stats
	set->mean.v[0] = m[0] / weight;
	set->mean.v[1] = m[1] / weight;
	set->mean.v[2] = atan2(m[3], m[2]);

	// Covariance in linear components
	for (j = 0; j < 2; j++)
		for (k = 0; k < 2; k++)
			set->cov.m[j][k] = c[j][k] / weight - set->mean.v[j] * set->mean.v[k];

	// Covariance in angular components; I think this is the correct
	// formula for circular statistics.
	set->cov.m[2][2] = -2 * log(sqrt(m[2] * m[2] + m[3] * m[3]));

	return;
}

void pf_set_selective_resampling(pf_t *pf, int selective_resampling)
{
	pf->selective_resampling = selective_resampling;
}

// Compute the CEP statistics (mean and variance).
void pf_get_cep_stats(pf_t *pf, pf_vector_t *mean, double *var)
{
	int i;
	double mn, mx, my, mrr;
	pf_sample_set_t *set;
	pf_sample_t *sample;

	set = pf->sets + pf->current_set;

	mn = 0.0;
	mx = 0.0;
	my = 0.0;
	mrr = 0.0;

	for (i = 0; i < set->sample_count; i++)
	{
		sample = set->samples + i;

		mn += sample->weight;
		mx += sample->weight * sample->pose.v[0];
		my += sample->weight * sample->pose.v[1];
		mrr += sample->weight * sample->pose.v[0] * sample->pose.v[0];
		mrr += sample->weight * sample->pose.v[1] * sample->pose.v[1];
	}

	assert(fabs(mn) >= DBL_EPSILON);
	if (fabs(mn) < DBL_EPSILON)
	{
		printf("ERROR : divide-by-zero exception : mn is zero\n");
		return;
	}

	mean->v[0] = mx / mn;
	mean->v[1] = my / mn;
	mean->v[2] = 0.0;

	*var = mrr / mn - (mx * mx / (mn * mn) + my * my / (mn * mn));

	return;
}

// Get the statistics for a particular cluster.
// 获取某个指定的【clabel】cluster的统计信息
int pf_get_cluster_stats(pf_t *pf, int clabel, double *weight,		
						 pf_vector_t *mean, pf_matrix_t *cov)
{
	pf_sample_set_t *set;
	pf_cluster_t *cluster;

	set = pf->sets + pf->current_set;

	if (clabel >= set->cluster_count)
		return 0;
	cluster = set->clusters + clabel;

	*weight = cluster->weight;
	*mean = cluster->mean;
	*cov = cluster->cov;

	return 1;
}

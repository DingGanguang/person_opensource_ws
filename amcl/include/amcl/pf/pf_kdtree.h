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
 * Desc: KD tree functions
 * Author: Andrew Howard
 * Date: 18 Dec 2002
 * CVS: $Id: pf_kdtree.h 6532 2008-06-11 02:45:56Z gbiggs $
 *************************************************************************/

#ifndef PF_KDTREE_H
#define PF_KDTREE_H

#ifdef INCLUDE_RTKGUI
#include "rtk.h"
#endif


// Info for a node in the tree
typedef struct pf_kdtree_node
{
  // Depth in the tree
  int leaf, depth;					//树的深度，叶子

  // Pivot dimension and value		 //分叉的维度和分叉值
  int pivot_dim;
  double pivot_value;

  // The key for this node			//这个节点的key，这里这个key有3个元素，在AMCL的情形下表示位姿
  int key[3];

  // The value for this node		//这个节点的value，在AMCL的情形下表示位姿的 权重
  double value;

  // The cluster label (leaf nodes)		// cluster 的标签，表示第几个cluster?
  int cluster;

  // Child nodes	
  struct pf_kdtree_node *children[2];		// 每个node 都有2个子节点

} pf_kdtree_node_t;


// A kd tree
typedef struct				// 定义一棵kd树：
{
  // Cell size
  double size[3];		//树的size

  // The root node of the tree     //树的根节点
  pf_kdtree_node_t *root;

  // The number of nodes in the tree		 //kd树的节点数
  int node_count, node_max_count;
  pf_kdtree_node_t *nodes;

  // The number of leaf nodes in the tree		  //这棵树的叶子节点数
  int leaf_count;

} pf_kdtree_t;


// Create a tree
extern pf_kdtree_t *pf_kdtree_alloc(int max_size);

// Destroy a tree
extern void pf_kdtree_free(pf_kdtree_t *self);

// Clear all entries from the tree
extern void pf_kdtree_clear(pf_kdtree_t *self);

// Insert a pose into the tree
extern void pf_kdtree_insert(pf_kdtree_t *self, pf_vector_t pose, double value);

// Cluster the leaves in the tree
extern void pf_kdtree_cluster(pf_kdtree_t *self);

// Determine the probability estimate for the given pose
extern double pf_kdtree_get_prob(pf_kdtree_t *self, pf_vector_t pose);

// Determine the cluster label for the given pose
extern int pf_kdtree_get_cluster(pf_kdtree_t *self, pf_vector_t pose);


#ifdef INCLUDE_RTKGUI

// Draw the tree
extern void pf_kdtree_draw(pf_kdtree_t *self, rtk_fig_t *fig);

#endif

#endif

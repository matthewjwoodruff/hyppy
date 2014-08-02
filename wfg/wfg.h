#ifndef _WFG_H_
#define _WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "avl.h"

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
      struct avl_node_t * tnode;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
} FRONT;

typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

FILECONTENTS *readFile(char[]);

void printContents(FILECONTENTS *);

int compare_tree_asc( const void *p1, const void *p2);
int greater(const void *v1, const void *v2);
void cleanup_filecontents(FILECONTENTS* filecontents);
double compute_hypervolume(FRONT* front);

extern int maxm;
extern int maxn;
extern FRONT* fs;
extern avl_tree_t *tree;
extern POINT ref;
extern int nobj;

#endif

#ifndef _WFG_H_
#define _WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
	int n_allocated_points;
	POINT *allocated_points;
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
void cleanup_point(POINT* point);
double compute_hypervolume(FRONT* front, POINT* referencepoint);
double hypervolume(int nobj, int npoints, double* values);

#endif

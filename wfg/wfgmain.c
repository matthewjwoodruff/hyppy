#include <sys/time.h>
#include <sys/resource.h>
#include "wfg.h"
#include "avl.h"

static void trimLine(char line[])
{
  int i = 0;

  while(line[i] != '\0')
  {
    if (line[i] == '\r' || line[i] == '\n')
    {
      line[i] = '\0';
      break;
    }
    i++;
  }
}

void printContents(FILECONTENTS *f)
{
  for (int i = 0; i < f->nFronts; i++)
  {
    printf("Front %d:\n", i+1);
    for (int j = 0; j < f->fronts[i].nPoints; j++)
    {
      printf("\t");
      for (int k = 0; k < f->fronts[i].n; k++)
      {
        printf("%f ", f->fronts[i].points[j].objectives[k]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

FILECONTENTS *readFile(char filename[])
{
  FILE *fp;
  char line[BUFSIZ];
  int front = 0, point = 0, objective = 0;

  FILECONTENTS *fc = malloc(sizeof(FILECONTENTS));
  fc->nFronts = 0;
  fc->fronts = NULL;

  fp = fopen(filename, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "File %s could not be opened\n", filename);
    exit(EXIT_FAILURE);
  }

  while(fgets(line, sizeof line, fp) != NULL)
  {
    trimLine(line);
    if (strcmp(line, "#") == 0)
    {
      front = fc->nFronts;
      fc->nFronts++;
      fc->fronts = realloc(fc->fronts, sizeof(FRONT) * fc->nFronts);
      fc->fronts[front].nPoints = 0;
      fc->fronts[front].points = NULL;
    }
    else
    {
      FRONT *f = &fc->fronts[front];
      point = f->nPoints;
      f->nPoints++;
      f->points = realloc(f->points, sizeof(POINT) * f->nPoints);
      f->n = 0;
      f->points[point].objectives = NULL;
      f->points[point].tnode = malloc(sizeof(avl_node_t));
      char *tok = strtok(line, " \t\n");
      do
      {
        POINT *p = &f->points[point];
        objective = f->n;
        f->n++;
        p->objectives = realloc(p->objectives, sizeof(OBJECTIVE) * f->n);
        p->objectives[objective] = atof(tok);
      } while ((tok = strtok(NULL, " \t\n")) != NULL);
    }
  }

  fc->nFronts--;
  // for (int i = 0; i < fc->nFronts; i++) fc->fronts[i].n = fc->fronts[i].points[0].nObjectives;
        fclose(fp);
  /* printf("Read %d fronts\n", fc->nFronts);
     printContents(fc); */
  return fc;
}

int main(int argc, char *argv[]) 
// processes each front from the file 
{
  // read.c allocates the FRONTs
  FILECONTENTS *f = readFile(argv[1]); // do this in python

  // find the biggest fronts
  // don't do this at all.  allocation at this stage is not what needs to be
  // optimized.  We burn all of our time within, not between fronts.
  // We should just have a method to compute HV for one front.
  // So our wrapper needs to figure out how big the front is and how many
  // objectives it has.
  // That's fine.  From the C perspective we're getting a simple-minded array, 
  // and we're just going to follow read.c to build the front structure.
  //

  // Here, we're figuring out how much space to allocate for the "stack"
  // of fronts.  I think this is taking the allocate-once optimization
  // overboard.  The big expense is not allocating space for these fronts,
  // as long as you do it once per hypervolume computation.
  for (int i = 0; i < f->nFronts; i++)
    {if (f->fronts[i].nPoints > maxm) maxm = f->fronts[i].nPoints;
     if (f->fronts[i].n       > maxn) maxn = f->fronts[i].n;
    }

  // Allocate memory for the FRONTs: enough for as many FRONTS 
  // as there are points in the biggest front.
  // What I'm going to change is that we'll just allocate enough 
  // for the current front.
  #if opt == 0
    fs = malloc(sizeof(FRONT) * maxm);
  #else

    // slicing (opt > 1) saves a level of recursion
    int maxd = maxn - (opt / 2 + 1); 
    fs = malloc(sizeof(FRONT) * maxd);

    // 3D base (opt = 3) needs space for the sentinels
    int maxp = maxm + 2 * (opt / 3);
    //int maxp = 100000;
    for (int i = 0; i < maxd; i++) 
      {fs[i].points = malloc(sizeof(POINT) * maxp); 
       for (int j = 0; j < maxp; j++) 
       {
         fs[i].points[j].tnode = malloc(sizeof(avl_node_t));
         // slicing (opt > 1) saves one extra objective at each level
         fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (maxn - (i + 1) * (opt / 2)));
       }
      }
  #endif

  tree = avl_alloc_tree ((avl_compare_t) compare_tree_asc,
                         (avl_freeitem_t) free);

  // initialise the reference point
  // This should be sent in because it is computed by the wrapper if not a
  // parameter of the wrapper.
  ref.objectives = malloc(sizeof(OBJECTIVE) * maxn);
  ref.tnode = malloc(sizeof(avl_node_t));
  if (argc == 2)
    {printf("No reference point provided: using the origin\n");
     for (int i = 0; i < maxn; i++) ref.objectives[i] = 0;
    }
  else if (argc - 2 != maxn)
    {printf("Your reference point should have %d values\n", maxn);
     return 0;
    }
  else 
  for (int i = 2; i < argc; i++) ref.objectives[i - 2] = atof(argv[i]);

  double hypervolume;
  for (int i = 0; i < f->nFronts; i++) 
    {      
      struct timeval tv1, tv2;
      struct rusage ru_before, ru_after;
      getrusage (RUSAGE_SELF, &ru_before);
 
      hypervolume = compute_hypervolume(&(f->fronts[i]));
      printf("hv(%d) = %1.10f\n", i+1, hypervolume);

      getrusage (RUSAGE_SELF, &ru_after);
      tv1 = ru_before.ru_utime;
      tv2 = ru_after.ru_utime;
      printf("Time: %f (s)\n", tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);
    }
  cleanup_filecontents(f);
  free(f);

  return 0;
  // Bottom line: we need to write, in C, a function that accepts the array of
  // objectives as a flat array of doubles (for simplicity), its length and
  // width as ints, and its reference point.  This function will construct a
  // FRONT f and then return hv(f).  It really shouldn't leak memory, so I also
  // need to figure out how to clean up a FRONT.
  // Valgrind is my friend here.
  // So it turns out read.c is my big leaker.
  // I wrote cleanup_filecontents and thereby got rid of all of my 
  // indirect leaks.  what's left is fs, the avl tree, the reference point, and
  // the points in fs.
}

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
  char** given_reference;
  given_reference = &(argv[2]);

  // initialise the reference point
  // This should be sent in because it is computed by the wrapper if not a
  // parameter of the wrapper.
  int refobj = argc - 2;
  POINT given_refpoint;
  given_refpoint.objectives = malloc(sizeof(OBJECTIVE) * refobj);
  given_refpoint.tnode = malloc(sizeof(avl_node_t));
  for (int i = 0; i < refobj; i++) {
    given_refpoint.objectives[i] = atof(given_reference[i]);
  }

  double hypervolume;
  for (int i = 0; i < f->nFronts; i++) {      
    struct timeval tv1, tv2;
    struct rusage ru_before, ru_after;
    getrusage (RUSAGE_SELF, &ru_before);

    FRONT* currentfront = &(f->fronts[i]);
    POINT* ref;
    POINT zero_refpoint;
    zero_refpoint.objectives = malloc(sizeof(OBJECTIVE) * currentfront->n);
    zero_refpoint.tnode = malloc(sizeof(avl_node_t));

    for (int i = 0; i < refobj; i++) {
      zero_refpoint.objectives[i] = 0.0;
    }

    if(currentfront->n <= refobj){
      ref = &given_refpoint;
    } else {
      ref = &zero_refpoint;
    }

    hypervolume = compute_hypervolume(&(f->fronts[i]), ref);
    cleanup_point(&zero_refpoint);
    printf("hv(%d) = %1.10f\n", i+1, hypervolume);

    getrusage (RUSAGE_SELF, &ru_after);
    tv1 = ru_before.ru_utime;
    tv2 = ru_after.ru_utime;
    printf("Time: %f (s)\n", tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);
  }
  cleanup_filecontents(f);
  free(f);

  return 0;
}

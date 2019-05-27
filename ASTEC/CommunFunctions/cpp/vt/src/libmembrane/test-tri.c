#include <vtmalloc.h>

#include <vt_tree.h>
#include <vt_common.h>

int main( int argc  __attribute__ ((unused)), char *argv[] )
{

    fprintf(stdout, "result\n");
    vt_pointList p;
    int i, axis;
    vt_point *point;

    VT_AllocPointList(&p, 10);


    double x[] = {1., 5., 4., 8 , 2, 6, 7, 5, 0, 9};
    double y[] = {1., 5., 4., 3 , 2, 6, 7, 5, 0, 9};
    double z[] = {1., 5., 4., 8 , 2, 3, 7, 5, 0, 9};

    for (i=0 ; i<p.n; i++)
    {
        point = &(p.list[i]);
        point->weight = 1;
        point->x = x[i];
        point->y = y[i];
        point->z = z[i];
    }
    int *ind = vtmalloc(10*sizeof(int), "ind", argv[0] );
    for(i = 0; i<p.n ; i++) ind[i]=i;

    axis = 0;
    tri_rapide(p, axis, 0, p.n-1, ind);
    fprintf(stdout, "indX = [%d", ind[0]);
    for(i = 1; i<p.n ; i++)
    {
        fprintf(stdout, ", %d", ind[i]);
    }
    fprintf(stdout, "]\n");

    axis = 1;
    for(i = 0; i<p.n ; i++) ind[i]=i;
    tri_rapide(p, axis, 0, p.n-1, ind);
    fprintf(stdout, "indY = [%d", ind[0]);
    for(i = 1; i<p.n ; i++)
    {
        fprintf(stdout, ", %d", ind[i]);
    }
    fprintf(stdout, "]\n");

    axis = 2;
    for(i = 0; i<p.n ; i++) ind[i]=i;
    tri_rapide(p, axis, 0, p.n-1, ind);
    fprintf(stdout, "indZ = [%d", ind[0]);
    for(i = 1; i<p.n ; i++)
    {
        fprintf(stdout, ", %d", ind[i]);
    }
    fprintf(stdout, "]\n");

fprintf(stdout, "expected:\n");

    char *txt="indX = [8, 0, 4, 2, 1, 7, 5, 6, 3, 9]\n\
indY = [8, 0, 4, 3, 2, 7, 1, 5, 6, 9]\n\
indZ = [8, 0, 4, 5, 2, 7, 1, 6, 3, 9]";
    fprintf(stdout, "%s\n", txt);

    return(1);
}

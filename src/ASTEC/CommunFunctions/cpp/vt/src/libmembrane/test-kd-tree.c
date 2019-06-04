#include <vt_tree.h>
#include <vt_common.h>

int main()
{

    vt_pointList p;
    int i;
    vt_point *point;

    tree_node *node;

    VT_AllocPointList(&p, 10);

    fprintf(stdout, "result with list:\n");

    double x[] = {1., 5., 4., 8 , 2, 6, 7, 5, 0, 9};
    double y[] = {1., 4., 4., 3 , 2, 6, 7, 5, 0, 9};
    double z[] = {1., 5., 4., 8 , 2, 3, 7, 5, 0, 9};

    for (i=0 ; i<p.n; i++)
    {
        point = &(p.list[i]);
        point->weight = 1;
        point->attribute = i;
        point->x = x[i];
        point->y = y[i];
        point->z = z[i];
    }

    vt_kdtree(&node, p, 0);

    VT_FreePointList(&p);

    printTree(node);

    /*
    clearTree(&node);

    fprintf(stdout, "result with image:\n");

    vt_image *img;

    img = _VT_Inrimage("/home/gmicheli/tmp.hdr");

    VT_AllocPointListWithImage(&p, img);


    VT_ExtractPointListWithImage(*img, &p);

    VT_FreeImage(img);


    vt_kdtree(&node, p, 0);

    VT_FreePointList(&p);

    printTree(node);

    */

    fprintf(stdout, "\nexpected: \n");

    char *txt="median = { 0.0, 0.0, 0.0 }, weight = 1.0\n\
median = { 1.0, 1.0, 1.0 }, weight = 1.0\n\
median = { 2.0, 2.0, 2.0 }, weight = 1.0\n\
median = { 4.0, 4.0, 4.0 }, weight = 1.0\n\
median = { 5.0, 4.0, 5.0 }, weight = 1.0\n\
median = { 5.0, 5.0, 5.0 }, weight = 1.0\n\
median = { 6.0, 6.0, 3.0 }, weight = 1.0\n\
median = { 8.0, 3.0, 8.0 }, weight = 1.0\n\
median = { 7.0, 7.0, 7.0 }, weight = 1.0\n\
median = { 9.0, 9.0, 9.0 }, weight = 1.0";
    fprintf(stdout, "%s\n", txt);



    /*  Nearest Neighbor Search */
    vt_point q;
    q.attribute = 0;
    q.weight = 1;
    q.x = 8;
    q.y = 0;
    q.z = 7;



    vt_point nearest;

    double d;

    fprintf(stdout, "\nTest VT_NNS (q = { %1.1lf, %1.1lf, %1.1lf }):\n\n", q.x, q.y, q.z);


    VT_NNS(q, node, &nearest, &d);

    fprintf(stdout, "nearest = { %1.1lf, %1.1lf, %1.1lf }, attribute = %d, d = %1.1f\n", nearest.x, nearest.y, nearest.z, nearest.attribute, d);

    fprintf(stdout, "\nexpected: \n");
    char *t="nearest = { 8.0, 3.0, 8.0 }";
    fprintf(stdout, "%s\n\n", t);

    clearTree(&node);





    return(1);
}

/* Poisson.c
 * Author: Jonathan Bao
 * contact: jjb8021@gmail.com
 * 2017.07.20
 * Description: A coding exercise from Burden and Faires - Numerical Analysis 9E.
 * To solve the Poisson equation with a finite difference algorithm.
 * Chapter 12.1 Example 2, p722. ISBN-13: 978-0-538-73351-9
 */

/* Poisson Equation Finite-Difference Solver
 *
 * to approximate the solution u(x,y) to the Poisson equation
 * u_xx (x,y) + u_yy (x,y) = f(x,y), a<=x<=b, c<=y<=d,
 * subject to the boundary conditions
 * u(x,y) = g(x,y) if x=a or x=b and c<=y<=d
 * and
 * u(x,y) = g(x,y) if y=c or y=d and a<=x<=b
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// grid
double a,b,c,d; //  endpoints (a,c) bottom left, (a,d) top left, (b,c) bottom right, (b,d) top right
int m,n; // resolution. m intervals in vertical (y) by n in horizontal (x) and m,n >= 3
double TOL; // tolerance
int N; // max number of iterations

double f(double x, double y)
{ // function on RHS of Poisson equation
    return x*exp(y);
}

int g(double x, double y, double* output)
{   // function describing the boundary conditions:
    /*
     *  linearly increasing and exponentially increasing
     */
    double e = 0.00005; // test for equality within a margin of error
    if (x < a+e && x > a-e)  // left edge
    {
        *output= 0.0;
        return 0;
    }
    else if (x < b+e && x> b-e)  // right edge
    {
        *output=  2.0*exp(y);
        return 0;
    }
    else if (y < c+e && y > c-e)  // bottom edge
    {
        *output=  x;
        return 0;
    }
    else if (y < d+e && y> d-e)  // top edge
    {
        *output=  x*M_E;
        return 0;
    }
    else
    {
        //printf("error: point(%3.2f,%3.2f) is not on the boundary.\n",x,y);
        return -1; // error: point(x,y) is not on the boundary
    }
}

// output solution to screen
int print_onearray(double* array, int size, char* label);
int print_twoarray(double** array, int n_x, int m_y, char* label);

int poisson()
{ // main solver
    /* grid specification
     * a rectangular region 0<=x<=2, 0<=y<=1
     * with dimensions 6x5 width by height in resolution
     */
    a=0, b=2;
    c=0, d=1;
    n=6, m=5;
    // iteration parameters
    TOL = 1e-10;
    N = 1e2;

    int i,j; // dummy index
    double h,k; // step size in x and y respectively
    double *x,*y; // x and y position coordinates
    double **w; // approximate solution of the grid

    h=(b-a)/(double)n;
    k=(d-c)/(double)m;
    x = (double*)calloc(sizeof(double),n+1);
    y = (double*)calloc(sizeof(double),m+1);

    // create mesh points
    for (i=0; i<n+1; i++)
    {
        x[i]=a+i*h;
    }
    for (i=0; i<m+1; i++)
    {
        y[i]=c+i*k;
    }
    print_onearray(x,n+1, "x coord");
    print_onearray(y,m+1, "y coord");

    w=(double**)calloc(sizeof(double*),n+1);
    for (i=0; i<n+1; i++)
    {
        w[i]=calloc(sizeof(double),m+1);
    }

    // calculate the boundary conditions
    double *bc = calloc(sizeof(double),1);
    for (i=0; i<n+1; i++)
    {
        for (j=0; j<m+1; j++)
        {
            *bc=0;
            g(x[i],y[j],bc);
            w[i][j]=*bc;
        }
    }
    print_twoarray(w,n+1,m+1, "boundary conditions");

    double lambda = h*h/(k*k);
    double mu = 2.0*(1.0+lambda);
    int L = 1; // iteration counter
    double z; // value of the approximate solution at a point
    double *g1= calloc(sizeof(double),1); // stores value of boundary points used in calculation
    double *g2 = calloc(sizeof(double),1);
    double NORM = 0;

    printf("Verbose mode to show intermediate solutions? 1 for yes.\n");
    int *verbose_mode=calloc(sizeof(int),1);
    scanf ("%i", verbose_mode);
    printf("\n");

    while (L<=N) // step 6
    {
        // Gauss-Seidel iterations
        // step 7
        // top left corner of interior mesh
        g(a,y[m-1],g1);
        g(x[1],d,g2);
        z = (-h*h*f(x[1],y[m-1])+(*g1)+lambda*(*g2)+lambda*w[1][m-2]+w[2][m-1])/mu;
        NORM = fabs(z-w[1][m-1]);
        w[1][m-1]=z;
        // step 8
        // top edge of interior mesh
        for (i=2; i < n-1; i++)
        {
            g(x[i],d,g1);
            z=(-h*h*f(x[i],y[m-1])+lambda*(*g1)+w[i-1][m-1]+w[i+1][m-1]+lambda*w[i][m-2])/mu;
            if( fabs(w[i][m-1]-z)> NORM)
            {
                NORM=fabs(w[i][m-1]-z);
            }
            w[i][m-1]=z;
        }
        // step 9
        // top right corner of interior mesh
        g(b,y[m-1],g1);
        g(x[n-1],d,g2);
        z = (-h*h*f(x[n-1],y[m-1])+(*g1)+lambda*(*g2)+w[n-2][m-1]+lambda*w[n-1][m-2])/mu;
        if(fabs(w[n-1][m-1]-z)> NORM)
        {
            NORM=fabs(w[n-1][m-1]-z);
        }
        w[n-1][m-1]=z;

        // step 10
        for (j=m-2; j>1; j--)
        {
            // step 11
            // leftmost endpoint of a row
            g(a,y[j],g1);
            z=(-h*h*f(x[1],y[j])+(*g1)+lambda*w[1][j+1]+lambda*w[1][j-1]+w[2][j])/mu;
            if (fabs(w[1][j]-z)>NORM)
            {
                NORM=fabs(w[1][j]-z);
            }
            w[1][j]=z;
            // step 12
            // interior points in a row
            for (i=2; i<n-1; i++)
            {
                z=(-h*h*f(x[i],y[j])+w[i-1][j]+lambda*w[i][j+1]+w[i+1][j]+lambda*w[i][j-1])/mu;
                if (fabs(w[i][j])>NORM)
                {
                    NORM=fabs(w[i][j]-z);
                }
                w[i][j]=z;
            }
            // step 13
            // rightmost endpoint of a row
            g(b,y[j],g1);
            z=(-h*h*f(x[n-1],y[j])+(*g1)+w[n-2][j]+lambda*w[n-1][j+1]+lambda*w[n-1][j-1])/mu;
            if(fabs(w[n-1][j]-z)>NORM)
            {
                NORM=fabs(w[n-1][j]-z);
            }
            w[n-1][j]=z;
        } // end steps 11 through 13
        // step 14
        // bottom left corner of interior mesh
        g(a,y[1],g1);
        g(x[1],c,g2);
        z=(-h*h*f(x[1],y[1])+(*g1)+lambda*(*g2)+lambda*w[1][2]+w[2][1])/mu;
        if(fabs(w[1][1]-z)>NORM)
        {
            NORM=fabs(w[1][1]-z);
        }
        w[1][1]=z;
        // step 15
        // bottom edge of interior mesh
        for(i=2; i<n-1; i++)
        {
            g(x[i],c,g1);
            z=(-h*h*f(x[i],y[1])+lambda*(*g1)+w[i-1][1]+lambda*w[i][2]+w[i+1][1])/mu;
            if(fabs(w[i][1]-z)>NORM)
            {
                NORM=fabs(w[i][1]-z);
            }
            w[i][1]=z;
        }
        // step 16
        // bottom right corner of interior mesh
        g(b,y[1],g1);
        g(x[n-1],c,g2);
        z=(-h*h*f(x[n-1],y[1])+(*g1)+lambda*(*g2)+w[n-2][1]+lambda*w[n-1][2])/mu;
        if(fabs(w[n-1][1]-z)>NORM)
        {
            NORM=fabs(w[n-1][1]-z);
        }
        w[n-1][1]=z;
        // step 17
        if(NORM<TOL)
        {
            // step 19
            printf("Congratulations. Solver converged to a solution!\nResults after %i iterations within a tolerance of %.10e\n\n",L,TOL);
            printf("x, y, approximate solution\n");
            // step 18
            // print solution
            for(i=1; i<n; i++)
            {
                for(j=1; j<m; j++)
                {
                    printf("%3.2f, %3.2f, %3.8f\n", x[i],y[j],w[i][j]);
                }
            }
            print_twoarray(w,n+1,m+1, "\n");
            printf("The procedure was successful.\n");
            return 0;
        }

        if (*verbose_mode==1)
        {
            // print intermediate solution
            printf("iteration %i\n",L);
            print_twoarray(w,n+1,m+1,"intermediate solution");
        }
        // step 20
        L++;
    } // end Gauss-Seidel iterations
    // step 21
    printf("Maximum number of iterations exceeded.\n");
    return 1;
}

int print_onearray(double* array, int size, char* label)
{ // prints one-dimensional numerical array (of double type) to screen
    if (label != NULL)
    {
        printf("%s\n", label);
    }
    for (int i=0; i<size; i++)
    {
        printf("%3.2f ", array[i]);
    }
    printf("\n\n");
    return 0;
}

int print_twoarray(double** array, int n_x, int m_y, char* label)
{ // prints two-dimensional numerical array (of double type) to screen
    if (label != NULL)
    {
        printf("%s\n", label);
    }
    for (int j=m_y-1; j>-1; j--)
    {
        for (int i=0; i<n_x; i++)
        {
            printf("%3.2f ", array[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return 0;
}

int main()
{
    poisson();
    return 0;
}

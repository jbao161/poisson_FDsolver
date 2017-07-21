# poisson_FDsolver
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

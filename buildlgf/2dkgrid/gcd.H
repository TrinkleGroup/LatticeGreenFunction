#ifndef __GCD_H
#define __GCD_H

// **** Some introductory number theory algorithms:
// iterative gcd, using Euclid's algorithm:
int gcd (int a, int b) 
{ int r, x=abs(a), y=abs(b);
 while (y!=0) {r = (x%y);  x=y;  y=r;}
 return x;}

// gcd(a,b,c) is the gcd of gcd(a,b) and c...
int gcd (int a, int b, int c) 
{ return gcd(gcd(a,b),c);}

#endif

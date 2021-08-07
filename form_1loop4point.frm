*Problem: Decompose the 4-leg loop integral into tensor integrals (Passarino-Veltman loop integrals) or A,B,C,D,... functions.
* The mass of all mediating particles inside the loop = m ;
*
**********<r1**********a***************p>**************d********r4=-(r1+r2+r3)>*********
*			*				*				
*			*				*
*			*				*				
*			^				p+r1+rr2+r3	
*			p+r1				v				
*			*				*
*			*				*				
*			*				*
*			*				*
**********<r2***********b**********<p+r1+r2*************c**************r3>***************
*1. Declaration:
symbols n,m,i1,i2,i3,i4,j,t,x;
dimension n;
**n is the number of dimension in dimensional regularization

indices a,b,c,d,o,s,u,v;

vectors r1,r2,r3,r4,p,k,q,r,g,k1,k2,k3;

cfunctions f,fiv,f1,f2,f3,f4,I;

cfunctions G,A,B,C,D,M;
functions A0,B0,C0,D0;
cfunction Dot(symmetric), Dot4;
ctensor ts;

set nt:i1,i2,i3,i4;
set ind:a,b,c,d;

******************************************************************************************
*2. Definitions:

local S = G(p+r1,p+r1+r2,p+r1+r2+r3,p,a,b,c,d);
** G is the 4-point function considered; the momenta and the gamma matrices' indices {a,b,c,d} are denoted following the figure above; 

id G(p?,q?,k?,r?,a?,b?,c?,d?)=( g_(1,a)*(g_(1,p)+m*gi_(1))*g_(1,b)*(g_(1,q)+m*gi_(1))*g_(1,c)*
(g_(1,k)+m*gi_(1))*g_(1,d)*(g_(1,r)+m*gi_(1)) )*f(p,q,k,r);
** f is the denominator part of the function ; 1 is chosen as the spin index
 
Tracen,1;

print +s S;
.sort


******************************************************************************************
*3. Simplifying the expression:
** Transform to tensor integrals:

*** Transform terms~p^2 
repeat;
id f(?a,q?)*p.p=f(?a,q)*(fiv(q)-2*Dot(p,q-p)-Dot(q-p,q-p)+m^2);
id Dot(p?,q?)=p.q;
id Dot(p?,0)=0;
id Dot(0,0)=0;
endrepeat;
**fiv(q)=q^2-m^2;


***Transform dot product p.ri --> ri.rj
repeat;
id f(?a)*p.r1 = f(?a)*(fiv(p+r1)-fiv(p)-r1.r1)/2;
id f(?a)*p.r2 = f(?a)*(fiv(p+r1+r2)- fiv(p+r1)-2*r1.r2-r2.r2)/2;
id f(?a)*p.r3 = f(?a)*(fiv(p+r1+r2+r3)-fiv(p+r1+r2)-2*r1.r3-2*r2.r3-r3.r3)/2;
endrepeat;

repeat;
id f(?a,p?,?b)*fiv(p?)=f(?a,?b);
endrepeat;
print +s S;
.sort

symmetrize f;


*** Change the momentum variable:
id f(p,?a)=f1(p,?a);
id f(?a, r1 + p)=f2(?a, r1 + p);

argument f2;
id p=p-r1;
endargument;

repeat;
id p(a?)*f2(?a)=(r(a)-r1(a))*f2(?a);
endrepeat;

id r=p;

print +s S;
.sort

id f1(?a) = f(?a);
id f2(?a) = f(?a);
symmetrize f;
print +s S;
.sort

argument f;
id p=0;
*f(p+r1, p+r1+r2, ....) --> f(r1, r1+r2, ...)
endargument;

id f(0,?a)=f(?a);


print +s S;
.sort

*** Transform to terms of tensor integral functions: A,B,C,D
totensor p,ts;

id f(?a)=f(?a)*ts();
id ts(?a)*ts()=ts(?a);


id f(p?)*ts(?a)=B(?a,M(p));
id f(p?,q?)*ts(?a)=C(?a,M(p,q));
id f(p?,q?,k?)*ts(?a)=D(?a,M(p,q,k));
id ts(?a)=A(?a);

id A?(?a,M(?b))=A(?a,M(?b))*M(?b);
print +s S;
.sort


*** Tensor integral coefficients Ai,Bi,...Aii,Bii,...:

**** Scalar:
id B(M(?a))=B0;
id C(M(?a))=C0;
id D(M(?a))=D0;

print +s S;
.sort


**** Tensor:

***Idea for the algorithm: .ex B(a,b,M(k1,k2)-->B(a,M(k1,k2))*B(b,M(k1,k2))-->[B(a,k1)+B(a,k2)]*[B(b,k1)+B(b,k2)]
*** --> [B(1)*k1(a)+B(2)*k2(a)]*[B(1)*k1(b)+B(2)*k2(b)] --> B(1,1)*k1(a)*k1(b)+B(1,2)*k1(a)*k2(b)+B(2,1)*k2(a)*k1(b)+B(2,2)*k2(a)*k2(b);

*** Transform B(a,b,c,...) --> B(a)*B(b)*B(c)...
repeat;
id A?(a?,?b,b?,M(?a))=A(a,M(?a))*A(?b,b,M(?a));
endrepeat;

*** Transform B(a,M(k1,k2,...))=B(1)*k1(a)+B(2)*k2(2)+....

**** g:additional element for tensor metric, then g(a)*g(b)=d_(a,b) 
argument A,B,C,D;
id M(?a)=M(g,?a,k);
endargument;
*k is just an intermediate variable


**** nt[t] gives information about the order of momentum arguments for each tensor component {a,b,c,d}: M(k1,k2,k3)--> k1-1; k2-2; k3-3;...:
repeat;
id A?(a?ind[t],M(p?,q?,?b))=nt[t]*(A(a,M(p))+A(a,M(q,?b)));
endrepeat;

id A?(a?,M(k))=0;

print +s S;
.sort

****...then k1-1 --> ex. B(1)*k1
id x?nt[t?]^j?*A?(a?ind[t?],M(p?))=A(j-1)*p(a);


**** Giving back the right form of tensor metric:
id once g(a?)*g(b?)*g(c?)*g(d?)=d_(a,b)*d_(c,d)+d_(a,c)*d_(b,d)+d_(a,d)*d_(b,c);
id once g(a?)*g(b?)=d_(a,b);

id g(a?)=0;*Because the tensor metric just replace for 2 indices

print +s S;
.sort

*** COLLIER's form:
repeat;
id A?(3)=A(0,0,0,1);
id A?(2)=A(0,0,1,0);
id A?(1)=A(0,1,0,0);
id A?(0)^2=A(1,0,0,0);
id A?(0)^4=A(2,0,0,0);
endrepeat;

repeat;
id A?(a?,b?,c?,d?)*A?(o?,s?,u?,v?)=A(a+o,b+s,c+u,d+v);
endrepeat;


id A?{B,C}(?a,a?,0)=A(?a,a);

*M brings the information of the momenta for the functions.
print +s S;
.end


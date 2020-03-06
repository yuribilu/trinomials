/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PARI macros and scripts used for the article 

"Trinomials, singular moduli and Riffaut's conjecture"

Yuri Bilu, Florian Luca, Amalia Pizarro-Madariaga


Table of contents


A. General scripts

maxabs(v) biggest  absolute value in a vector
minabs(v) smallest  absolute value in a vector
kpower(n,k) biggest kth power dividing an integer
prodprimes(x) product of primes \le x
primeprod(X) biggest prime P such that \Prod_{p\le P}<X 


B. Scripts on discriminants and singular moduli

taus(Delta) finding the tau_x of given discriminant Delta
sigmods(Delta) finding the singular moduli of given discriminant Delta
deltas(X) discriminants up to X
deltas_ge_h(X,h) list of discriminants up to X with class number at least h


C. Scripts for Subsection 6.1

quadnonreszero(p) quadratic non-residues mod p including 0
negzerkron(x) odd residue classes with negative or zero Kronecker at all primes \le x
negzerokrondiscrs(X,x) odd discriminants up to X which have negative or zero Kronecker at all primes \le x 
suit(Delta) finding the suitable integers of given discriminant Delta
antitrin(x,y,z) max(0, 1-|z/y|-(2*|y/x|+2*|y/x|^3))  
antitrindis(Delta) a lower bound for max{antitrin(x,y,z) : x,y,z singular moduli of discriminant Delta}
antitrindislist(v) finding the smallest antitrindis\ge 0.1 for the discriminants in a list, and printing the list of discriminants with antitrindis<0.1


D. Scripts for Subsection 6.2

deltas_three() generating the list of discriminants of class number 3
fk(F,k) calculating polynomials F_k
prnu(Delta) finding entries of Table 1 for a given discriminant
wrlatex(v,filename) write in a file vector in latex array row format
allintable(filename) generate latex code for Table 1 and write it in a file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/ 


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A. General scripts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


/*--------------------------------------------------------------------------------
biggest  absolute value in a vector
input: vector with real or complex entries
output: biggest absolute value of its entries
-----------------------------------------------------------------------------------*/

maxabs(v)=vecmax(abs(v));


/*--------------------------------------------------------------------------------
smallest  absolute value in a vector
input: vector with real or complex entries
output: smallest absolute value of its entries
-----------------------------------------------------------------------------------*/

minabs(v)=vecmin(abs(v)); 


/*-------------------------------------------------------
biggest kth power dividing an integer
input: an integer n, a positive integer k
output: the biggest positive integers m such that m^k|n
------------------------------------------------------------*/

kpower(n,k) =

{
my(m,f,p,nu, s);
f=factor(n);
p=f[,1];  /* prime factors of n */
nu=f[,2]; /* exponents of the prime factors */ 

/* euclidean division of the exponents by k and suppressing the remainders */ 

s=length(p);
for (i=1,s,
nu[i]=(nu[i]-(nu[i]%k))/k;
);

/* calculating m */ 

m=1;
for(i=1,s,
m=m*p[i]^nu[i];
);

return(m);
}


/*--------------------------------------------------------
product of primes \le x
input: x\ge 2
output: \prod_{p\le x}p
------------------------------------------------------*/

prodprimes(x)=
{
my(p,P);
if (x<2, print("x too small"); return(0););
P=1;
p=2;
for (i=1,primepi(x),
P=P*p;
p=nextprime(p+1););
return(P);
}


/*------------------------------------------------------
biggest prime p such that prodprimes(p)<X 
input X \ge 3
output p
---------------------------------------------------------*/

primeprod(X)=
{
my(p,q,P,n);
if(X<3,print("X too small");return(0););
n=10*log(X);
p=2;
P=2;
for(i=1,n,
q=nextprime(p+1);
P=P*q;
if (P>=X,return(p););
p=q;);
print("internal error");
return(0);
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
B. Scripts on discriminants and singular moduli
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*--------------------------------------------------------------------------
finding the tau_x of given discriminant Delta
input: Delta (must be <0, =0,1mod4, otherwise error message)
output list of all \tau_x
----------------------------------------------------------------------------- */

taus(Delta) =
{
    my(a,b,c,clno,m,v);

     if (Delta>0, 
    print("Bad Delta");
    return(0););


    if (Delta%4==2, 
    print("Bad Delta");
    return(0););
	
	    if (Delta%4==3, 
    print("Bad Delta");
    return(0););
    
   
    clno  =qfbclassno(Delta);

    X = abs(Delta);
	
	v=vector(clno,i,0); 
	
    maxa = floor(sqrt(X/3));



    for(a=1,maxa,
	for(b=-a+1,a,
	    if((b^2-Delta)%(4*a) != 0, next;);
	    c = (b^2-Delta)/(4*a);
	    
	    /*
	      Note that b>-a, next two ifs check that
	      we are in the fundamental domain.
	    */	 

	    if(a>c,next;);
	    if( (a==c) && (b<0),next;);

	    /* (a,b,c) must be coprime. */
	    if(gcd(a,gcd(b,c))!=1,next;);
	    
	    
	    m++;
		
			
v[m]=(b+I*sqrt(X))/(2*a);		 
	);
    );

/* Sanity check. Make sure our n matches with Pari's computation */
/* Used only for debugging purposes */
    if(m != clno,
        return("FAIL (internal error)");
    	return(0);
    ); 
    return(v);
}




/*------------------------------------------------------------
finding the singular moduli of given discriminant Delta
--------------------------------------------------------------*/

sigmods(Delta) = vector(qfbclassno(Delta),i, ellj(taus(Delta)[i])); 





/*--------------------------------------------------------------------
list of discriminants up to X
input: X\ge 3
output: discriminants Delta  in [-X,0]
----------------------------------------------------------------------*/

deltas(X)=
{
my(m,n,i);
if (X<3, print("X too small");return(0););
m=floor((X+1)/4);
n=floor(X/4);
v=vector(m+n,i,0);
for(i=1,m,
v[2*i-1]=1-4*i;);
for(i=1,n,
v[2*i]=-4*i;);
return(v);
}

/*--------------------------------------------------------------------
list of discriminants up to X with class number at least h
input: X\ge 3, h\ge 1
output: discriminants Delta  in [-X,0] with h(Delta) \ge h
----------------------------------------------------------------------*/

deltas_ge_h(X,h)=
{
my(u,v,w,n,i,j);
if (X<3, print("X too small");return(0););
if (h<1, print("h too small");return(0););
u=deltas(X);
n=length(u);
v=vector(n,i,0);
j=1;
for(i=1,n,



if(qfbclassno(u[i])>=h, 
v[j]=u[i];
j++;););




if(j==1, print("X too small");return(0););
w=vector(j-1,i,0);
w=v[1..j-1];
return(w);
}








/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C. Scripts for Subsection 6.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/



/*------------------------------------------------------------
list of quadratic non-residues mod p including 0
input: odd prime p
output: the list
-------------------------------------------------------------*/

quadnonreszero(p)=
{
my(v,m,n);
n=2;
v=vector((p+1)/2,i,0);
v[1]=Mod(0,p);
for(i=1,p-1,
m=Mod(i,p);
if(issquare(m)==0, v[n]=m; n++););
return(v);
}





/*--------------------------------------------------
odd residue classes with negative or zero Kronecker at all primes \le x
input: x\ge 2
output: classes a mod 4\prod_{p\le x}p which are 0,4,5 mod 8 and (a/p)=0,-1 for all odd p\le x
------------------------------------------------*/
negzerkron(x) = 
{
my(u,v,w,p,pi,n,m);

if(x<2,print("x too small"); return(0););

v=[Mod(5,8)];
p=2;
pi=primepi(x);
for(i=1,pi-1,
p=nextprime(p+1);
u=quadnonreszero(p);
n=1;
w=vector(length(u)*length(v),i,0);
for(j=1,length(u),
for(k=1,length(v),
w[n]=chinese(u[j],v[k]);
n++;
););
v=w;
);

return(v);
}


/*-------------------------------------------------------
Odd discriminants up to X which have negative or zero Kronecker at all primes \le x 
input: X\ge 12, x\ge nextprime(primeprod(X/4)+1)
output: list of discriminants
----------------------------------------------------------*/


negzerokrondiscrs(X,x) =

{
my(p,q,r,c,d,Q,u,v,V,w,lu, lw,n);
if(X<12,print("x too small"); return(0););
p=primeprod(X/4);
q=nextprime(p+1);
if(x<q,print("x too small"); return(0););

/* forming the list of residues */

u=negzerkron(p);

print("List of "length(u)" residues formed");
print("p="p);

/* sieving mod q */

w=quadnonreszero(q);
lu=length(u);
lw=length(w);
v=vector(lu*lw,i,0);
Q=4*prodprimes(q);
n=0;

for(i=1,lu,

for(j=1,lw,

/* writing a discriminant into the list */ 

c=chinese(u[i],w[j]);

d=Q-lift(c);

if(d<=X,n++;v[n]=-d;);  

);

);



if(n==0,
print("No discriminants");
return(););

u=v[1..n];


print("List of "n" discriminants formed");
print("q="q);



r=nextprime(q+1);
if (r>x, return(u);
);



/* siving mod bigger primes */


for(k=1,primepi(x)-primepi(q),
print("r="r);
n=0;
lu=length(u);
lv=lu;
v=vector(lv,i,0);

for(i=1,lu,
if(issquare(Mod(u[i],r))==0, n++; v[n]=u[i];);
if(u[i]%r==0, n++; v[n]=u[i];);
);

if(n==0,
print("No discriminants");
return();

);

print(n" discriminants left");

u=v[1..n];
r=nextprime(r+1);
);

return(u);
}




/*--------------------------------------------------------------------------
finding the suitable integers of given discriminant Delta
input: Delta (must be <0, =0,1mod4, otherwise error message)
output list of [a,b], with each suitable integers occuring exactly once as a
----------------------------------------------------------------------------- */

suit(Delta) =
{
    my(a,b,c,clno,m,v);

     if (Delta>0, 
    print("Bad Delta");
    return(0););


    if (Delta%4==2, 
    print("Bad Delta");
    return(0););
	
	    if (Delta%4==3, 
    print("Bad Delta");
    return(0););
    
   
    clno  =qfbclassno(Delta);

    X = abs(Delta);
	
	v=vector(clno,i,0); 
	
    maxa = floor(sqrt(X/3));

m=0;

    for(a=1,maxa,
	for(b=0,a,
	    if((b^2-Delta)%(4*a) != 0, next;);
	    c = (b^2-Delta)/(4*a);
	    
 

	    if(a>c,next;);
	    if( (a==c) && (b<0),next;);

	    /* (a,b,c) must be coprime. */
	    if(gcd(a,gcd(b,c))!=1,next;);
	    
	    
	    m++;
		
			
v[m]=[a,b];
next(2);		 
	);
    );
v=v[1..m];
    return(v);
}







/*-------------------------------------------------------------------------------
max(0, 1-|z/y|-(2*|y/x|+2*|y/x|^3))  
input: complex numbers x,y,z with x\ne 0
output: the value
------------------------------------------------------------------------------------*/

antitrin(x,y,z) = max(0, 1-abs(z/y)-(2*abs(y/x)+2*(abs(y/x))^3))


/*------------------------------------------------------------------------------
a lower bound for max{antitrin(x,y,z) : x,y,z singular moduli of discriminant Delta}
input: an imaginary quadratic discriminant Delta
output: some d>0 satisfying d\le antitrin(x,y,z) for some choice of singular moduli x,y,z of discriminant Delta; 
0 if antitrin(x,y,z)=0 for all choices of x,y,z that we tried
remark: we do not claim to calculate the actual maximum
-----------------------------------------------------------------------------------------*/

antitrindis(Delta) =

{my(v,w,X,m,a,b,x,y,z,t,d);

if (qfbclassno(Delta) < 4,

return(0);

);

v=suit(Delta);

X=abs(Delta);

m=length(v);


/* The case when there are more than 2 suitable integers */

if (m>2,


a=v[1][1];
b=v[1][2];
x=ellj((b+sqrt(X)*I)/(2*a));


a=v[2][1];
b=v[2][2];
y=ellj((b+sqrt(X)*I)/(2*a));


a=v[m][1];
b=v[m][2];
z=ellj((b+sqrt(X)*I)/(2*a));

/* for some discriminants (like -667) accidentally |y|<|z| */

if (abs(y)<abs(z),

t=y; 
y=z;
z=t;

);

d=antitrin(x,y,z);
return(d); 

);

/* The case of exactly 2 suitable integers 
(there are very few discriminants like this) */ 

w=sigmods(Delta); 

x=maxabs(w);
y=maxabs(w[2..length(w)]);
z=minabs(w);

d=antitrin(x,y,z);

return(d); 

}



/*-------------------------------------------------------------------------------
finding the smallest antitrindis\ge 0.1 for the discriminants in a list, and printing the list of discriminants with antitrindis<0.1
input: a list of imaginary quadratic discriminants
output: the smallest value of antitrindis \ge 0.1
------------------------------------------------------------------------------------*/

antitrindislist(v) = 

{

my(m,d);

m=length(v);

d=1;

for (i=1,m,

if (antitrindis(v[i]) <0.1, 
print(v[i]" is a bad discriminant");
next;
);

d=min(d, antitrindis(v[i]));
);

return(d); 

}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C. Scripts for Subsection 6.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


/*---------------------------------------------------------------------------
generating the list of discriminants of class number 3
input: none
output: the list
------------------------------------------------------------------------*/

deltas_three()=
{
my(u,v,n,i,j);


u=deltas(1000);
n=length(u);
v=vector(25,i,0);
j=0;

for(i=1,n,

if(qfbclassno(u[i])==3, 
j++;
v[j]=u[i];

);

);

if(j<25, print("Something went wrong");return(0););

return(v);
}


/*-----------------------------------------
calculating polynomials F_k
input: F a monic polynomial of degree 3 with main variable x, k a positive integer
output: F_k
---------------------------------------------*/

fk(F,k) =

{my(a,b,c,a2,b2,A,B,ak,bk);

if (poldegree(F)!=3,
print("F must be of degree 3");
return(0);
);

if (polcoef(F,3)!=1,
print("F must be monic");
return(0);
);

if (k==0,
return((x-1)^3);
);

if (k==1,
return(F);
);

a=polcoef(F,2); 

b=polcoef(F,1); 

c=polcoef(F,0); 

a2=-a^2+2*b;

b2=b^2-2*a*c;

if (k==2,
return
(x^3+a2*x^2+b2*x-c^2);
);

A=[-3,a,a2];
B=[3,b,b2];

for(i=3,k,

ak=-a*A[3]-b*A[2]-c*A[1];

bk=b*B[3]-a*c*B[2]+c^2*B[1];

A=[A[2],A[3],ak];

B=[B[2],B[3],bk];

);

return(x^3+ak*x^2+bk*x-(-c)^k);

}




/*-----------------------------------------------------
finding entries of Table 1 for a given discriminant
input: a discriminant Delta with h=3
output: 9-component vector [Delta,d,p,f,nu_p(c),r_0,nu_0,lambda,mu]
Remark: for all discriminants we have f=2 and nu_p(c)=3; hence they do not appear in Table 1
-----------------------------------------------------*/

prnu(Delta) = 

{

my(F,a,b,c,d,v,m,p,f,Fk,D,r,nu,lambda,mu);

if(qfbclassno(Delta)!=3,
print("Delta must be of class number 3"); 

return(0);
);

F=polclass(Delta);
c=polcoef(F,0);
b=polcoef(F,1);
a=polcoef(F,2);

/* generating polynomial F */ 

d=gcd(a,gcd(kpower(b,2),kpower(c,3)));
F=(d^(-3))*subst(F,x, d*x);

c=polcoef(F,0);
b=polcoef(F,1);
a=polcoef(F,2);

/* prime divisors of c */
v=factor(c)[,1];

/* finding biggest p dividing c but not b */ 

m=length(v);

for (i=0,m-1,
p=v[m-i];
if((Delta*b)%p!=0,
break;
);

if (i=m-1,
print("no p found");
return(0);
);
);

/* finding f */ 

f=1;

if(kronecker(Delta,p)==-1,
f=2;
);


/* finding r_0  */

fordiv (p^f-1, k, 

Fk=fk(F,k);
D=poldisc(Fk);

if(D%p==0,
r=k;

break;

);

);

/* finding nu_0 */ 

nu=valuation(D,p)/f;

/* the case p=2 is special */ 

if (p==2,

D=poldisc(fk(F,2*r));

nu=valuation(D,p)/f-1;

); 

/* finding lambda & mu */

roots= d^(-1)*sigmods(Delta);
X=roots[1];
Y=roots[2];

lambda=(3/2)*((log(abs(X)))/(log(abs(X/Y))))*(p^nu)/r; 

mu=(3*log(2)+log(2.01))/(log(abs(X/Y)))*(p^nu)/r; 


return([Delta,d,p,valuation(c,p),f,r,nu,lambda,mu]);

}


/*----------------------------------------------------------------
write vector in a file in latex array row format
input: v vector of real numbers, filename a file
output: same vector printed in the file in the formal of one row of latex array
------------------------------------------------------------------*/

wrlatex(v,filename) = 

{my(m);

m=length(v);

for(i=1,m,
write(filename,v[i]"&");
);

return();

}


/*-------------------------------------------------------------------------
generate latex code for Table 1 and write it in a file
input: filename
output: code written in filename
------------------------------------------------------------------------*/

allintable(filename) =

{

my(deltas,ds,ps,nupcs,fs,rs,nus,lambdas,mus,v,Delta,roots,p,r,nu,X,Y,d,lambda,mu);

deltas=deltas_three();

ds=vector(25,i,0); 

ps=vector(25,i,0); 

nupcs=vector(25,i,0); 

fs=vector(25,i,0); 

rs=vector(25,i,0); 

nus=vector(25,i,0); 

lambdas=vector(25,i,0); 

mus=vector(25,i,0); 

for(i=1,25,

Delta=deltas[i];

v=prnu(Delta); 

ds[i]=v[2];

ps[i]=v[3];

nupcs[i]=v[4];

fs[i]=v[5];

rs[i]=v[6]; 

nus[i]=v[7];

lambdas[i]=v[8];

mus[i]=v[9];


);

wrlatex(deltas,filename);

wrlatex(ds,filename);

wrlatex(ps,filename);

wrlatex(nupcs,filename);

wrlatex(fs,filename);

wrlatex(rs,filename);

wrlatex(nus,filename);

wrlatex(lambdas,filename);

wrlatex(mus,filename);

return();

}


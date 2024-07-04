#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define BASE 10

//-----------Prototipos------------
void inverso(int *v, int *r, int n);
void valor_inicial(int x, int *r, int n);
int * karatsuba(int *a, int *b, int n);
int *copia(int *a,int i,int j,int k);
int compara(int *a,int *b,int n);
int restal(int *a,int n,int *b,int m,int k);
void sumal(int *a,int n,int *b,int m,int k);
void fixmul(int *a, int *b, int *r, int n);
void recortar_karatsuba(int *c, int n);
void raiz(int *a, int *r, int n);
void divide(int * a, int *d, int *r, int n);
void inverso_raiz(int *v, int *r, int n);
void valor_inicial_raiz(int x, int *r, int n);
int * arquimedes(int n_lados, int n);

void main(int argc, char *argv[]){
	if(argc!=3){
		printf("\nUso: Numero de decimales, Numero de iteraciones");
		exit(1);
	}
    int n = atoi(argv[1]);
	int n_iter = atoi(argv[2]);
	int i,j;
	n++;
    int *r = arquimedes(n_iter, n);
    printf("\nAproximacion de Pi con %d decimales y %d iteraciones: ", n-1, n_iter);
    for(i=n-1,j=1; i>=0; i--, j++){
		if(j==2)printf(".");
		printf("%d", r[i]);
	}
}

int * arquimedes(int n_iter, int n){
	
	int i, j;
	int *tres = (int *)calloc(n, sizeof(int));
	int *dos = (int *)calloc(n, sizeof(int));
	int *p = (int *)calloc(n, sizeof(int));
	int *p1 = (int *)calloc(n, sizeof(int));
	int *q = (int *)calloc(n, sizeof(int));
	int *q1 = (int *)calloc(n, sizeof(int));
    int *r = (int *)calloc(n, sizeof(int));

	for(i=0; i<n-1; i++){tres[i]=0;dos[i]=0;}
	tres[n-1]=3;
	dos[n-1]=2;
	
	raiz(tres, p, n);
	fixmul(p, tres, q, n);

	free(tres);

	divide(q, dos, p, n);

	for(i=1;i<=n_iter;i++)
	{
		printf("\n%da iteracion ------", i);
		printf("\nP: ");for(j=n-1;j>=0;j--){printf("%d", p[j]);}
		printf("\nQ: ");for(j=n-1;j>=0;j--){printf("%d", q[j]);}
		inverso(p, p1, n);
		inverso(q, q1, n);
		sumal(p1, n, q1, n, 0);
		divide(p1, dos, r, n);
		inverso(r, q, n);
		fixmul(q, p, r, n);
		raiz(r, p, n);
	}
	sumal(p,n,q,n,0);
	divide(p, dos, r, n);
	free(dos);free(p);free(p1);free(q);free(q1);
	return r;

}

//DEVUELVE EL INVERSO (1/X) SIN REALIZAR UNA DIVISION
//1/d xk+1= x(2-dx) 
void inverso(int *v, int *r, int n){
	int i, j;
	int *x = (int *)calloc(n, sizeof(int));
	int *dos = (int*)calloc(n, (sizeof(int)));
	
	if(v[n-1] == 0){ valor_inicial(2, x, n); }
	else{ valor_inicial(v[n-1], x, n); }	
	//2^k>=n

	int iter = ((int)log2(n))+1;
	for(i=0; i<iter; i++){
		for(j=0; j<n-1; j++){dos[j]=0;}dos[n-1]=2;
		fixmul(x, v, r, n);
		restal(dos, n, r, n, 0);
		fixmul(x, dos, x, n);
    }
	free(dos);
	for(i=0; i<n; i++){
		r[i]=x[i];
	}

	free(x);
}

void valor_inicial(int x, int *r, int n){
	for(int i=0; i<n; i++){
		r[i]=0;
	}
	if(n<=5){
		printf("\nEl numero de decimales debe ser > 4.");
		exit(1);
	}
	else if(x>=8){
		r[n-2]=1;
		r[n-3]=2;
		r[n-4]=5;	
	} 
  	else if(x>=4){
		r[n-2]=2;
		r[n-3]=5;
	} 
    else {
		r[n-2]=5;
	}
}

void recortar_karatsuba(int *c, int n){
	int i, j;
	for(i=n-1, j=(2*n)-2; i>=0; j--, i--){
		c[i] = c[j];
	}
}

//MULTIPLICACION CON COMPROBACION DE OVERFLOW
void fixmul(int *a, int *b, int *r, int n){
    if(a[n-1]*b[n-1]>=10){
        printf("\nOverflow ");
		exit(1);
    }
	else if((a[n-1]&&b[n-1]>=3)&&(a[n-2]&&b[n-2]>=2)&&(a[n-3]&&b[n-3]>=7)){
		printf("\nOverflow ");
		exit(1);
	}
	else{
		int *x = karatsuba(a, b, n);
		recortar_karatsuba(x, n);
    	for(int i=0; i<n; i++){r[i]=x[i];}
	}
}


// k es el desplazamiento a la izquierda, n el tamaño de a y m el de b
// deja el resultado en a
void sumal(int *a,int n,int *b,int m,int k){
  int l,carry=0;
  for(l=0;l<m;++l) 
  {
    a[l+k]+=b[l]+carry;
    carry=a[l+k]/BASE;
    a[l+k]%=BASE;
  }
  if(carry>0) a[m+k]+=carry;	
}

// k es el desplazamiento a la izquierda, n el tamaño de a y m el de b
// deja el resultado en a
// devuelve el acarreo
int restal(int *a,int n,int *b,int m,int k){
  int l,carry=0;
  for(l=0;l<m;++l) 
  {
    a[l+k]-=b[l]+carry;
    if(a[l+k]<0) { a[l+k]+=BASE;carry=1; }
    else         carry=0;
  }
  while((l+k<n)&&carry)
  {
    a[l+k]-=carry;
    if(a[l+k]<0) { a[l+k]+=BASE;carry=1; }
    else         carry=0;
    ++l;
  }
  return carry;	
}

// devuelve +1, 0 o -1 si a>b, a=b o a<b 
int compara(int *a,int *b,int n){
  int i;
  for(i=n-1;i>=0;--i)
  {
    if(a[i]>b[i]) return 1;
    if(a[i]<b[i]) return -1;
  }
  return 0;
}

// el ultimo parametro es para la memoria a pedir
// y debe ser mayor o igual que la necesaria: (j-i+1)
int *copia(int *a,int i,int j,int k){
	int l,*b=(int *)calloc(k,sizeof(int));
	for(l=i;l<=j;++l) b[l-i]=a[l];
	return b;
}

int * karatsuba(int *a, int *b, int n){
	int k=(n+1)/2, restar=1;
	int *c, *a1,*a2, *b1, *b2, *t1, *t2, *p1, *p2, *p3;
	if(n>1){
		a1=copia(a,0,k-1,k);
		a2=copia(a,k,n-1,k);
		b1=copia(b,0,k-1,k);
		b2=copia(b,k,n-1,k);
		p1=karatsuba(a1,b1,k);
		p2=karatsuba(a2,b2,n-k);
		
		if(compara(a1,a2,k)>=0){
			restal(a1,k,a2,k,0);
			t1=a1;
		}
		else{
			restal(a2,k,a1,k,0);
			t1=a2;
			restar=1-restar;
		}
		
		if(compara(b1,b2,k)>=0){
			restal(b1,k,b2,k,0);
			t2=b1;
		}
		else{
			restal(b2,k,b1,k,0);
			t2=b2;
			restar=1-restar;
		}
		
		p3=karatsuba(t1,t2,k);
		
		free(a1);free(a2);free(b1);free(b2);
		c = copia(p1, 0, 2*k-1, 2*n);
		sumal(c,2*n,p2,2*(n-k),k);
		sumal(c,2*n,p1,2*k,k);
		if(restar){
			restal(c,2*n,p3,2*k,k);
		}
		else{
			sumal(c,2*n,p3,2*k,k);
		}
		sumal(c,2*n,p2,2*(n-k),2*k);
		free(p1);free(p2);free(p3);
	}
	else{
		c=(int *)calloc(2, sizeof(int));
		c[0]=a[0]*b[0];
		c[1]=c[0]/BASE;
		c[0]%=BASE;
	}
	return c;
}

// division {xk+1 = xk + xk(1-dxk)}
// recibe dos vectores 'a' y 'd' y guarda el resultado en 'r'
void divide(int * a, int *d, int *r, int n){
	// 1/'d'
    inverso(d, r, n);
	// (1/d) * a
	fixmul(r, a, r, n);
}

// raiz cuadrada de un vector 'a'
// guarda el resultado en 'r'
void raiz(int *a, int *r, int n){
    // 1/raiz(a)
    inverso_raiz(a, r, n);
    // 1/raiz(a) * a
	fixmul(a, r, r, n);
}

// funcion para calcular 1/raiz(v)
// guarda el resultado en 'r'
void inverso_raiz(int *v, int *r, int n){
	int i, j;
	int *x = (int *)calloc(n, sizeof(int));
	int *uno = (int *)calloc(n, sizeof(int));
	int *dos = (int *)calloc(n, sizeof(int));

	// creamos el vector dos
	for(i=0; i<n; i++){dos[i]=0;}dos[n-1]=2;
	valor_inicial_raiz(v[n-1], x, n);

	int iter = ((int)log2(n))+1;
	for(i=0; i<iter; i++){
		//Aqui termina la comprobacion cutre
		// x = x^2
		fixmul(x, x, r, n);
		// r = a*(x^2)
		fixmul(v, r, r, n);
		// 1-a*(x^2)
		for(j=0; j<n; j++){uno[j]=0;}uno[n-1]=1;
		restal(uno, n, r, n, 0);
		// el resultado de la resta se guarda en el vector "uno"
		// (1-a*(x^2))/2
		divide(uno, dos, r, n);
		// x * (1-a*(x^2))/2
		fixmul(x, r, r, n);
		// x + x * (1-a*(x^2))/2
		sumal(x, n, r, n, 0);
	}
	free(uno);free(dos);
	for(i=0; i<n; i++){
		r[i]=x[i];
	}
	free(x);
}  

void valor_inicial_raiz(int x, int *r, int n){
	for(int i=0; i<n; i++){
		r[i]=0;
	}
	r[n-2]=2;
}



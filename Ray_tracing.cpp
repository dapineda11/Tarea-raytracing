/*
Pontificia Universidad Javeriana
Daniel Andrés Pineda Suárez
Proyecto Primer corte DisPro: Ray_Tracing
Ingeniería Electrónica
2018
                                                               
                 ,(%(.     (   *(%(.   *#(.,(#((((#. *##((#*   .*%#,     *(     *##    .*#*     ./                     
                  .%,     #%#    #%    ,(   *#.       #/   %%    %/     .#%,     ,#%(    (     ,#%/                    
                  .%,    #..%*    #(   (    *%,..(,   #/ ./%,    %/    .( ,%,    ,( (%,  (     ( *%.                   
                  .%,   ////(%*    %( (     *#.  .,   #/ (%      %/    #///#%    ,(  ,#( (    #///(%,                  
                . .%.  (,    /%.   ,%#.     *#.    /  #/  .%(    %/   (     (#   ,(    /%(   #.    (%                  
                /#/  .//*    ///*   ..     *///////  *//*   *//,*//,,//,   .///,,//,     / ,//*   .///*                
                                                 //,/*,              *.                                                  */
                                                 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SIZE 512
#define SKY 1.0e+8

typedef struct VEC3_T VEC3_T;
struct VEC3_T{
	double x,y,z;
};

typedef struct RAY_T RAY_T;
struct RAY_T{
	VEC3_T direction, origin;
};

typedef struct SPH_T SPH_T;
struct SPH_T{
	VEC3_T center;
	double radius;
};

VEC3_T Vec_Create (double, double, double);
int Vec_Print (VEC3_T);
VEC3_T Vec_Addition (VEC3_T, VEC3_T);
VEC3_T Vec_Subtraction (VEC3_T, VEC3_T);
VEC3_T Vec_Scale (VEC3_T, double);
double Vec_Dot_Product (VEC3_T, VEC3_T);
VEC3_T Vec_Normalize (VEC3_T);
int Sphere_Ray_Intersection (SPH_T, RAY_T, double*);


 static SPH_T my_sphere = { { 0.0, 0.0, 0.0 }, 4.2 };
 static RAY_T ray = { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, -10.0 } };
 static VEC3_T light = { 0.0, 1.0, 0.0};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 int main ()
{
  double distance = 0.0;

  fprintf (stdout, "P2\n%d %d\n255\n", SIZE, SIZE ); //Conversión a P2 (Escala de Grises)

  VEC3_T position;

  for (int i=SIZE/2; i>-SIZE/2 ;i--)
  {
  	for (int j=-SIZE/2; j<SIZE/2 ; j++)
  	{
  		position = Vec_Create(i*(10.0/SIZE), j*(10.0/SIZE), 0.0 );
  		ray.direction = Vec_Subtraction (position, ray.origin);
  		ray.direction = Vec_Normalize (ray.direction);

  		Sphere_Ray_Intersection (my_sphere, ray, &distance);
  	

	  }
	  fprintf(stdout,"\n");
  }

  return 0;
} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
VEC3_T Vec_Create (double x, double y, double z)
{
	VEC3_T tmp;
	tmp.x=x;
	tmp.y=y;
	tmp.z=z;

	return tmp;
} 

int Vec_Print (VEC3_T vector)
{
	fprintf(stderr, "x: %f\n", vector.x);
	fprintf(stderr, "y: %f\n", vector.y);
	fprintf(stderr, "z: %f\n", vector.z);

	return 1;
} //Vec_Print/

VEC3_T Vec_Addition (VEC3_T summand_1, VEC3_T summand_2)
{
	VEC3_T tmp;
	tmp.x=summand_1.x+summand_2.x;
	tmp.y=summand_1.y+summand_2.y;
	tmp.z=summand_1.z+summand_2.z;

	return tmp;
}//Vec_Addition/

VEC3_T Vec_Subtraction (VEC3_T subtrahend, VEC3_T minuend)
{
	VEC3_T tmp;
	tmp.x=subtrahend.x-minuend.x;
	tmp.y=subtrahend.y-minuend.y;
	tmp.z=subtrahend.z-minuend.z;

	return tmp;
}//Vec_Subtraction/

VEC3_T Vec_Scale (VEC3_T vector, double factor)
{
	VEC3_T tmp;
	tmp.x= vector.x * factor;
	tmp.y= vector.y * factor;
	tmp.z= vector.z * factor;

	return tmp;
}//Vec_Scale/

double Vec_Dot_Product (VEC3_T vector_1, VEC3_T vector_2)
{
	double dot_product = 0.0;
	dot_product  += vector_1.x * vector_2.x;
	dot_product  += vector_1.y * vector_2.y;
	dot_product  += vector_1.z * vector_2.z;

	return dot_product;
} //Vec_Dot_Product/

VEC3_T Vec_Normalize (VEC3_T vector)
{
	double tmp= Vec_Dot_Product(vector, vector);
	if (tmp==0.0)
	{ exit(1);
	}
	return 	Vec_Scale(vector, 1.0/sqrt(tmp));

} //Vec_Normalize/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Sphere_Ray_Intersection (SPH_T sphere, RAY_T ray, double* distance) //Gráfico de Intersección
{
	VEC3_T voc= Vec_Subtraction (ray.origin, sphere.center);
	double vd = Vec_Dot_Product(voc, ray.direction);
	double voc_mag_sqr = Vec_Dot_Product (voc, voc);
	double r_sqr = sphere.radius*sphere.radius;

	double discriminant = vd*vd-voc_mag_sqr+r_sqr;
		unsigned int ans1 = (255*(-vd + sqrt(discriminant)))/10;
		unsigned int ans2 = (255*(-vd - sqrt(discriminant)))/10;

        if(ans1<=ans2){
                fprintf(stdout, "%d ", ans1);
        }
		        else if(ans2==0){
           fprintf(stdout, "%d ", ans1);
        }
		         else if(ans1>ans2){
                 fprintf(stdout, "%d ", ans2);
                }
				        else if(ans1==0 & ans2==0){
                  fprintf(stdout, "0 ");
        }                      
		                      else if(ans1==0){
                  fprintf(stdout, "%d ", ans2);
        }


}



//Integra SPICE con el integrador Bulirsch-Stoer para solucionar la EOM de partículas de prueba

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<SpiceUsr.h>
//#include<malloc.h>

#define X 0
#define Y 1
#define Z 2

#define VX 3
#define VY 4
#define VZ 5

//=============================================================== 
//                   Función de graficación.
//=============================================================== 
int  grafica(double tiempo, int t);
//=============================================================== 
//           _-_ Importación de Kernels de SPICE -_-
//=============================================================== 
void load_kernel(){
  //furnsh_c("de430.bsp"); 
  //furnsh_c("de421.bsp"); 
  //furnsh_c("de432s.bsp"); furnsh_c("de414.bsp"); furnsh_c("de418.bsp");  
  furnsh_c("de431_part-1.bsp");  
  furnsh_c("de431_part-2.bsp");  
  furnsh_c("naif0011.tls");
}
//=============================================================== 
//-----Definición de estructuras-------
//=============================================================== 
typedef struct
{
  int id;
  SpiceDouble pos[3];
  SpiceDouble vel[3];
  SpiceDouble mass;
}Bodies;    //Nombre de la variable

Bodies *plansys,*plansys_aux, *minimoon;

int n_bodies = 2;
int n_mmoons = 1;
int Npart;
int N_bulirsch_stoer = 12; // this numbers should be multiple of 2
int Ndivitions = 6;
int Nfiles;


#include "bs_integ.c"



//=============================================================== 
//---------Funcion para contar líneas ------------
//=============================================================== 

int counterLines(char *infile) 
{
  int NDAT,c;
  static FILE *pf;
  if((pf=fopen(infile,"r"))==NULL)
    printf("no puedo abrir archivo %s\n",infile);
  NDAT=0;
  while((c=fgetc(pf))!= EOF)
    {
      if(c == '\n')
        {
          ++NDAT;
        }
    }
  fclose(pf);
  return NDAT;
}

//=============================================================== 
//           -----------Función Principal ------------------
//=============================================================== 

int main(int argc, char *argv[]){

  int i,j; //
  int n, m = 7; //
  //int m = 7; //id,x,y,z,vx,vy,vz.
  Npart=n_mmoons;



  //=============================================================== 
  //------- Carga condiciones iniciales de la simulacion ----------
  //Cond iniciales respecto a la Tierra X,Y,Z en AU, Vx,Vy,Vz en m/s
  //=============================================================== 

  char namefiles[1000]; 
  FILE *pf; 
  Nfiles = atoi(argv[1]); 
  
  int c;
  for(i=0; i<Nfiles; i++) 
    { 
      sprintf(namefiles,"../inf_sphere/inifiles/file%.3d.dat",i); 
      pf = fopen(namefiles,"r"); 
      n_mmoons = counterLines(namefiles);
      n = n_bodies + n_mmoons; 
      minimoon=malloc(n_mmoons*sizeof(Bodies));     //Alocación de memoria 
     
      //while((c=fgetc(pf))!= '\n')

      for(j=0; j<Nfiles; j++) 
	{
	  fscanf(pf,"%lf %lf %lf %lf %lf %lf",
		 &minimoon[j].pos[X],&minimoon[j].pos[Y],&minimoon[j].pos[Z], 
		 &minimoon[j].vel[X],&minimoon[j].vel[Y],&minimoon[j].vel[Z]);
	  minimoon[j].id = j;
	  minimoon[j].mass = 1.0;
	}
      fclose(pf);
    }


  plansys=malloc(n_bodies*sizeof(Bodies));      //Alocación de memoria 
  plansys_aux=malloc(n_bodies*sizeof(Bodies));  //Alocación de memoria  
     

  //  printf("File loaded. Example minimoon 4:\n");
  int k=4;
  // printf("%.17e %.17e %.17e %.17e %.17e %.17e \n",
  //	 minimoon[k].pos[X],minimoon[k].pos[Y],minimoon[k].pos[Z],
  //	 minimoon[k].vel[X],minimoon[k].vel[X],minimoon[k].vel[X]);

  //sizeof()
  

  //-----------------------------------------------
  // System dynamical properties
  //-----------------------------------------------
  SpiceDouble UM=1.989e30;  //kg
  SpiceDouble UL=149597870700 ;  //UA in m
  SpiceDouble UTs=sqrt(UL*UL*UL /(6.67384e-11 * UM) ); //UT en segundos
  SpiceDouble UT=UTs/3.154e+7; //  //UT en años
  SpiceDouble UV=UL/UTs; 
  int t;
  
  printf("UV: %e %e\n", UV/1000,   UL*UL*UL / (UTs*UTs *UM) );

  //UT=5.03e6 s / (365*24*3600)

  //-----------------------------------------------
  // ---------   Spice Initialization   ----------
  //-----------------------------------------------  

  load_kernel();
  SpiceDouble tiempoluz, time, tini, tend, tiempo, step;
  ConstSpiceChar* fecha_ini= "01/01/2015 00:00:00.000";
  ConstSpiceChar* fecha_end= "01/01/2016 00:00:00.000";
  SpiceDouble sun_state[6], mercury_state[6],venus_state[6], earth_state[6], moon_state[6], mars_state[6];   
  
  //-----------------------------------------------
  //------  Gravitational Scenery SPICE  --------
  //-----------------------------------------------    
  // Reference System
  char  RefSys[20]="SUN";
  SpiceDouble ti, te;
  str2et_c(fecha_ini,&ti);
  str2et_c(fecha_end,&te);

  tini=0.0;
  tend=(te-ti)/UTs;
  step=(tend-tini)/100.0;

  printf("Tini= %lf, \t Tend= %lf, step %f\n",tini,tend, step);
  printf("Tini= %lf, \t Tend= %lf, step %f\n",tini*UT,tend*UT, step*UT);

//  exit(0);
  SpiceDouble spice_time;
  
  //-----------------------------------------------
  // ---------  Frame  Transformation    ----------
  //-----------------------------------------------
  
  spice_time=tini;
  spkezr_c("EARTH"  ,spice_time, "J2000", "NONE", RefSys, earth_state  , &tiempoluz);

  earth_state[X]*=1000.0;
  earth_state[Y]*=1000.0;
  earth_state[Z]*=1000.0;
  earth_state[VX]*=1000.0;
  earth_state[VY]*=1000.0;
  earth_state[VZ]*=1000.0;

 /*
  
  for (k=0; k<= n_mmoons;k++)    
    {
      minimoon[k].pos[X]= (minimoon[k].pos[X] + earth_state[X]);
      minimoon[k].pos[Y]= (minimoon[k].pos[Y] + earth_state[Y]);
      minimoon[k].pos[Z]= (minimoon[k].pos[Z] + earth_state[Z]);
      minimoon[k].vel[X]= (minimoon[k].pos[X] + earth_state[VX]);
      minimoon[k].vel[Y]= (minimoon[k].pos[Y] + earth_state[VY]);
      minimoon[k].vel[Z]= (minimoon[k].pos[Z] + earth_state[VZ]);
    }
*/

  for (k=0; k<= n_mmoons;k++)    
    {
      minimoon[k].pos[X]= -earth_state[X];
      minimoon[k].pos[Y]= -earth_state[Y];
      minimoon[k].pos[Z]= -earth_state[Z];
      minimoon[k].vel[X]= -earth_state[VX];
      minimoon[k].vel[Y]= -earth_state[VY];
      minimoon[k].vel[Z]= -earth_state[VZ];
    }


  i=2;

  printf("R_Tierra=%e \n",
	 sqrt(earth_state[X]*earth_state[X] + earth_state[Y]*earth_state[Y] +earth_state[Z]*earth_state[Z]));
  
  printf("R_Mmoon=%e \n",
	 sqrt(minimoon[i].pos[X]*minimoon[i].pos[X] + minimoon[i].pos[Y]*minimoon[i].pos[Y] + minimoon[i].pos[Z]*minimoon[i].pos[Z]));

  //-----------------------------------------------  
  // Masses assignament
  //-----------------------------------------------
  plansys[0].mass=1.0;
  plansys[1].mass=3.30e23/UM;
  plansys[2].mass=4.87e24/UM;
  plansys[3].mass=5.97e24/UM;
  plansys[4].mass=7.35e22/UM;
  plansys[5].mass=6.42e23/UM;
  
  FILE *state, *excluded, *accepted;
  FILE *Forbit;
  
  Forbit = fopen("./output/minimoon.dat","a");  
  excluded=fopen("./output/excluded.txt", "a");
  accepted=fopen("./output/accepted.txt", "a");
  
  char nombre[1000];

  SpiceDouble dtl, rcol;


  rcol=6400000; //        m  -distancia de colision
  //-----------------------------------------------  
  // Temporal Evolution Each Particle
  //-----------------------------------------------

  
  for(i=1; i<2; i++) {
    //n_mmoons; i++) {
    t=0;
    
    for(tiempo=tini;tiempo<=tend;tiempo+=step){
      // printf("iteracion: %f\n",tiempo);
      
      t+=1;
      sprintf(nombre,"./output/system_state_%.4d.txt",t);
      state=fopen(nombre, "w");
      
      spice_time=tiempo;
      
      //        Target,  time  ,Frame   ,Correction, RefSystem, get_state , light-correct
      spkezr_c("0"      ,spice_time, "J2000", "NONE", RefSys, sun_state    , &tiempoluz);
      spkezr_c("MERCURY",spice_time, "J2000", "NONE", RefSys, mercury_state, &tiempoluz);
      spkezr_c("VENUS"  ,spice_time, "J2000", "NONE", RefSys, venus_state  , &tiempoluz);
      spkezr_c("EARTH"  ,spice_time, "J2000", "NONE", RefSys, earth_state  , &tiempoluz);
      spkezr_c("MOON"   ,spice_time, "J2000", "NONE", RefSys, moon_state   , &tiempoluz);
      spkezr_c("4"      ,spice_time, "J2000", "NONE", RefSys, mars_state   , &tiempoluz);

      /*     
      //for (k=0; k<=n_mmoons;k++){
      printf("%.17e %.17e %.17e %.17e %.17e %.17e \n",
	     minimoon[i].pos[X],minimoon[i].pos[Y],minimoon[i].pos[Z],
	     minimoon[i].vel[X],minimoon[i].vel[Y],minimoon[i].vel[Z]);
	     //}
	     
	     */
      
      
      //      printf("Time: %f %f\n", spice_time, UT);
      
      //      if (Spice_time>=30
     

  for (k=0; k<= n_mmoons;k++)    
    {
      minimoon[k].pos[X]= -earth_state[X];
      minimoon[k].pos[Y]= -earth_state[Y];
      minimoon[k].pos[Z]= -earth_state[Z];
      minimoon[k].vel[X]= -earth_state[VX];
      minimoon[k].vel[Y]= -earth_state[VY];
      minimoon[k].vel[Z]= -earth_state[VZ];
    }


      minimoon[i].pos[X] = minimoon[i].pos[X]/UL *1000;
      minimoon[i].pos[Y] = minimoon[i].pos[Y]/UL *1000;
      minimoon[i].pos[Z] = minimoon[i].pos[Z]/UL *1000;
      minimoon[i].vel[X] = minimoon[i].vel[X]/UV *1000;
      minimoon[i].vel[Y] = minimoon[i].vel[Y]/UV *1000; 
      minimoon[i].vel[Z] = minimoon[i].vel[Z]/UV *1000;
      
      
      // SUN
      plansys[0].pos[0]=sun_state[X]/UL  *1000;
      plansys[0].pos[1]=sun_state[Y]/UL  *1000;
      plansys[0].pos[2]=sun_state[Z]/UL  *1000;
      plansys[0].vel[0]=sun_state[VX]/UV *1000;
      plansys[0].vel[1]=sun_state[VY]/UV *1000;
      plansys[0].vel[2]=sun_state[VZ]/UV *1000;
      
      // Mercury
      plansys[1].pos[0]=mercury_state[X]/UL  *1000;
      plansys[1].pos[1]=mercury_state[Y]/UL  *1000;
      plansys[1].pos[2]=mercury_state[Z]/UL  *1000;
      plansys[1].vel[0]=mercury_state[VX]/UV *1000;
      plansys[1].vel[1]=mercury_state[VY]/UV *1000;
      plansys[1].vel[2]=mercury_state[VZ]/UV *1000;
      
      // Venus
      plansys[2].pos[0]=venus_state[X]/UL  *1000;
      plansys[2].pos[1]=venus_state[Y]/UL  *1000;
      plansys[2].pos[2]=venus_state[Z]/UL  *1000;
      plansys[2].vel[0]=venus_state[VX]/UV *1000;
      plansys[2].vel[1]=venus_state[VY]/UV *1000;
      plansys[2].vel[2]=venus_state[VZ]/UV *1000;
      
      // Earth
      plansys[3].pos[0]=earth_state[X]/UL  *1000;
      plansys[3].pos[1]=earth_state[Y]/UL  *1000;
      plansys[3].pos[2]=earth_state[Z]/UL  *1000;
      plansys[3].vel[0]=earth_state[VX]/UV *1000;
      plansys[3].vel[1]=earth_state[VY]/UV *1000;
      plansys[3].vel[2]=earth_state[VZ]/UV *1000;
      
      // Moon
      plansys[4].pos[0]=moon_state[X]/UL  *1000;
      plansys[4].pos[1]=moon_state[Y]/UL  *1000;
      plansys[4].pos[2]=moon_state[Z]/UL  *1000;
      plansys[4].vel[0]=moon_state[VX]/UV *1000;
      plansys[4].vel[1]=moon_state[VY]/UV *1000;
      plansys[4].vel[2]=moon_state[VZ]/UV *1000;
      
      // Mars
      plansys[5].pos[0]=mars_state[X]/UL  *1000;
      plansys[5].pos[1]=mars_state[Y]/UL  *1000;
      plansys[5].pos[2]=mars_state[Z]/UL  *1000;
      plansys[5].vel[0]=mars_state[VX]/UV *1000;
      plansys[5].vel[1]=mars_state[VY]/UV *1000;
      plansys[5].vel[2]=mars_state[VZ]/UV *1000;
      


  printf("R_Tierra=%e \n",
	 sqrt(earth_state[X]*earth_state[X] + earth_state[Y]*earth_state[Y] +earth_state[Z]*earth_state[Z]));
  
  printf("R_Mmoon=%e \n",
	 sqrt(minimoon[i].pos[X]*minimoon[i].pos[X] + minimoon[i].pos[Y]*minimoon[i].pos[Y] + minimoon[i].pos[Z]*minimoon[i].pos[Z]));

  //exit(0);      

      fprintf(state,"%.17e %.17e %.17e %.17e %.17e %.17e \n",
	      sun_state[0], sun_state[1], sun_state[2], 
	      sun_state[3], sun_state[4], sun_state[5]);
      
      fprintf(state,"%.17e %.17e %.17e %.17e %.17e %.17e \n",
	      mercury_state[0], mercury_state[1], mercury_state[2], 
	      mercury_state[3], mercury_state[4], mercury_state[5]);
      
      fprintf(state,"%.17e %.17e %.17e %.17e %.17e %.17e \n",
	      venus_state[0], venus_state[1], venus_state[2], 
	      venus_state[3], venus_state[4], venus_state[5]);

      fprintf(state,"%.17e %.17e %.17e %.17e %.17e %.17e \n",
	      earth_state[0], earth_state[1], earth_state[2], 
	      earth_state[3], earth_state[4], earth_state[5]);

      fprintf(state,"%.17e %.17e %.17e %.17e %.17e %.17e \n",
	      moon_state[0], moon_state[1], moon_state[2], 
	      moon_state[3], moon_state[4], moon_state[5]);

      fprintf(state,"%.17e %.17e %.17e %.17e %.17e %.17e \n",
	      mars_state[0], mars_state[1], mars_state[2], 
	      mars_state[3], mars_state[4], mars_state[5]);

      fprintf(Forbit,"%.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e \n",
	      minimoon[i].pos[X],minimoon[i].pos[Y],minimoon[i].pos[Z],
	      minimoon[i].vel[X],minimoon[i].vel[Y],minimoon[i].vel[Z],
	      sqrt(minimoon[i].pos[X]*minimoon[i].pos[X] +
		   minimoon[i].pos[Y]*minimoon[i].pos[Y] +
		   minimoon[i].pos[Z]*minimoon[i].pos[Z]));
      
      //Condiciones de parada de las integraciones individuales:
      //Distancia asteroide 
      
      dtl= sqrt((earth_state[X] - minimoon[i].pos[X]) * (earth_state[X] - minimoon[i].pos[X]) +
		(earth_state[Y] - minimoon[i].pos[Y]) * (earth_state[Y] - minimoon[i].pos[Y]) +
		(earth_state[Z] - minimoon[i].pos[Z]) * (earth_state[Z] - minimoon[i].pos[Z]));


      if (dtl<=rcol) {
	fprintf(excluded," %d \n", minimoon[i].id );
	break;
      }
      

      bulirsch_stoer(step);
      fclose(state);

      grafica(tiempo,t);   
      
    }
    
    fclose(excluded);
    fclose(accepted);
    fclose(Forbit);    
  }    
  free(plansys);
  free(plansys_aux);
  free(minimoon);
  printf("Done\n");
  
    
  return 0;
  
}

  
int  grafica(double tiempo, int t){
  FILE *gplt;
  gplt=fopen("./output/animation.gpl", "w");
  fprintf(gplt,"reset \n");
  fprintf(gplt,"UA=149597871.0; sc=2 \n");
  fprintf(gplt,"set terminal png size 800,600 \n");
  fprintf(gplt,"set size ratio -1 \n");
  fprintf(gplt,"set  xlabel 'x [AU] '\n");
  fprintf(gplt,"set  xlabel 'y [AU] '\n");
  fprintf(gplt,"set output './plots/st%.4d.png' \n ",t);
  fprintf(gplt,"set title 'Time=%lf s '\n ",tiempo);
  //fprintf(gplt,"splot[-2e9:2e9][-2e9:2e9][-2e9:2e9] './state_%lf.txt' u 1:2:3 lt 7 ps 1 not \n ",tiempo);
  fprintf(gplt,"plot[-sc:sc][-sc:sc] './output/system_state_%.4d.txt' u ($1/UA):($2/UA) lt 7 ps 1 t 'Solar System', './output/minimoon.dat' u ($1/UA):($2/UA) w p lt 7 lc 11  t 'minimoons' \n",t);
  //fprintf(gplt,"plot './output/system_state_%.4d.txt' u 1:2 lt 7 ps 1 t 'Solar System', './output/minimoon.dat' u 1:2 w p lt 7 lc 11  t 'minimoons' \n",t);
  fprintf(gplt,"reset \n");
  fclose(gplt);
  system("gnuplot ./output/animation.gpl");
  return 0;
}

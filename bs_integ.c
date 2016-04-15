//-------------------------------------------------------------------------//
//Definición de estructuras 
//-------------------------------------------------------------------------//

typedef struct
{
  double pos[3];
  double vel[3];
  double accel[3];
  double mass;
}be_particles;

be_particles *be_part;

double H;

//-------------------------------------------------------------------------//
//Llamado a funciones del programa: aceleracion, integrador e interpolador.
//-------------------------------------------------------------------------//
int midpoint_modified(int N, 
		      double be_poly_pos[Ndivitions][Npart][3],
		      double be_poly_vel[Ndivitions][Npart][3],
		      double be_poly_hk[Ndivitions]);
//-------------------------------------------------------------------------//
int Compute_acceleration(int i);
//-------------------------------------------------------------------------//
//Rutina de integración Bulirsch-Stoer: Ingresa paso del tiempo dt.
//-------------------------------------------------------------------------//

int bulirsch_stoer(double dt){

  int i,j,k,N;
  double t,h;

  double pos[Npart][3],vel[Npart][3],a[Ndivitions];
  double be_poly_pos[Ndivitions][Npart][3];
  double be_poly_vel[Ndivitions][Npart][3];
  double be_poly_hk[Ndivitions];
  double Poscm[3],Velcm[3];

  //FILE *Forbit;
  //Forbit = fopen("./output/minimoon.dat","a");

  be_part = (be_particles * )malloc( (Npart)*sizeof(be_particles) );
  if(be_part == NULL){
    printf("Allocation of be_part failed\n");
    exit(0);
  }

  H = dt;
    
  for(N=2; N<=N_bulirsch_stoer; N = N+2)
    {
      midpoint_modified(N, be_poly_pos, be_poly_vel, be_poly_hk);
    }
  
  h = 0.0;
  for(k=0; k<Npart; k++)
    {
      pos[k][0] = pos[k][1] = pos[k][2] = 0.0;
      vel[k][0] = vel[k][1] = vel[k][2] = 0.0;
    }  
  
  for(i=0; i<Ndivitions; i++)
    {
      a[i] = 1.0;
      
      for(j=0; j<Ndivitions; j++)
	{
	  if(i!=j)
	    a[i] = a[i] * ( h - be_poly_hk[j] )/( be_poly_hk[i] - be_poly_hk[j] ); 
	}
      
      for(k=0; k<Npart; k++)
	{
	  vel[k][0] = vel[k][0] + a[i] * be_poly_vel[i][k][0];
	  vel[k][1] = vel[k][1] + a[i] * be_poly_vel[i][k][1];
	  vel[k][2] = vel[k][2] + a[i] * be_poly_vel[i][k][2];
	  
	  pos[k][0] = pos[k][0] + a[i] * be_poly_pos[i][k][0];
	  pos[k][1] = pos[k][1] + a[i] * be_poly_pos[i][k][1];
	  pos[k][2] = pos[k][2] + a[i] * be_poly_pos[i][k][2];
	}
    }
  
  //printf("control 3 t = %lf N = %d\n",t,N);    
  
  for(i=0; i<Npart; i++)
      
  {
    minimoon[i].vel[0] = vel[i][0];
    minimoon[i].vel[1] = vel[i][1];
    minimoon[i].vel[2] = vel[i][2];
    
    minimoon[i].pos[0] =pos[i][0];
    minimoon[i].pos[1] =pos[i][1];
    minimoon[i].pos[2] =pos[i][2];
    }
  

  //Almacena posiciones y velocidades 
  /*

  fprintf(Forbit,"%.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e \n", 
	  minimoon[i].pos[X],minimoon[i].pos[Y],minimoon[i].pos[Z],minimoon[i].vel[X],minimoon[i].vel[Y],minimoon[i].vel[Z], 
	  sqrt(minimoon[i].pos[X]*minimoon[i].pos[X] + minimoon[i].pos[Y]*minimoon[i].pos[Y]+minimoon[i].pos[Z]*minimoon[i].pos[Z]));

  
    fprintf(Forbit,"%.17e  %.17e  %.17e  %.17e  %.17e  %.17e  %.17e \n",
    minimoon[i].pos[X], minimoon[i].pos[Y], minimoon[i].pos[Z],
    minimoon[i].vel[X], minimoon[i].vel[Y], minimoon[i].vel[Z],
    sqrt(minimoon[i].pos[X]*minimoon[i].pos[X] + 
    minimoon[i].pos[Y]*minimoon[i].pos[Y] +  
    minimoon[i].pos[Z]*minimoon[i].pos[Z]));  
  */
  //fclose(Forbit);  

  free(be_part); 
  return 0;
}
int midpoint_modified(int N, 
		      double be_poly_pos[Ndivitions][Npart][3],
		      double be_poly_vel[Ndivitions][Npart][3],
		      double be_poly_hk[Ndivitions])
{
  int i,k,n;
  double h,rminus2[N+1][Npart][3],vminus2[N+1][Npart][3];
  
  k = N/2 - 1;
  //printf("k = %d\n",k);
  h = be_poly_hk[k] = H/N;
  
  //printf("calculo hk = %lf\n",be_poly_hk[k]);
  // y_0 = y(t)
  for(i=0; i<Npart; i++)
    { 
      rminus2[0][i][0] = be_part[i].pos[0] = minimoon[i].pos[0];
      rminus2[0][i][1] = be_part[i].pos[1] = minimoon[i].pos[1];
      rminus2[0][i][2] = be_part[i].pos[2] = minimoon[i].pos[2];
      
      vminus2[0][i][0] = be_part[i].vel[0] = minimoon[i].vel[0];
      vminus2[0][i][1] = be_part[i].vel[1] = minimoon[i].vel[1];
      vminus2[0][i][2] = be_part[i].vel[2] = minimoon[i].vel[2]; 
    }
  
  //printf("calculo y_0\n");
  // y_1 = y_0 + hF_0
  for(i=0; i<Npart; i++)
    {

      Compute_acceleration(i);

      vminus2[1][i][0] = be_part[i].vel[0] = be_part[i].vel[0] + h*be_part[i].accel[0];  
      vminus2[1][i][1] = be_part[i].vel[1] = be_part[i].vel[1] + h*be_part[i].accel[1];  
      vminus2[1][i][2] = be_part[i].vel[2] = be_part[i].vel[2] + h*be_part[i].accel[2];

      rminus2[1][i][0] = be_part[i].pos[0] = be_part[i].pos[0] + h*be_part[i].vel[0];  
      rminus2[1][i][1] = be_part[i].pos[1] = be_part[i].pos[1] + h*be_part[i].vel[1];  
      rminus2[1][i][2] = be_part[i].pos[2] = be_part[i].pos[2] + h*be_part[i].vel[2];
    }
  //printf("calculo y_1\n");
  // y_n = y_{n-2} + 2hF_{n-1} 
  for(n=2; n<=N; n++)
    {
      for(i=0; i<Npart; i++)
	{	  
	  Compute_acceleration(i);
	  
	  vminus2[n][i][0] = be_part[i].vel[0] = vminus2[n-2][i][0] + 2.0*h*be_part[i].accel[0];  
	  vminus2[n][i][1] = be_part[i].vel[1] = vminus2[n-2][i][1] + 2.0*h*be_part[i].accel[1];  
	  vminus2[n][i][2] = be_part[i].vel[2] = vminus2[n-2][i][2] + 2.0*h*be_part[i].accel[2];
	  
	  rminus2[n][i][0] = be_part[i].pos[0] = rminus2[n-2][i][0] + 2.0*h*be_part[i].vel[0];  
	  rminus2[n][i][1] = be_part[i].pos[1] = rminus2[n-2][i][1] + 2.0*h*be_part[i].vel[1];  
	  rminus2[n][i][2] = be_part[i].pos[2] = rminus2[n-2][i][2] + 2.0*h*be_part[i].vel[2];
	}
      
    }
  
  //printf("calculo y_n\n");
  for(i=0; i<Npart; i++)
    {
      Compute_acceleration(i);
      
      be_part[i].vel[0] = 0.5*( be_part[i].vel[0] + ( vminus2[N-1][i][0] + h*be_part[i].accel[0]) );
      be_part[i].vel[1] = 0.5*( be_part[i].vel[1] + ( vminus2[N-1][i][1] + h*be_part[i].accel[1]) );
      be_part[i].vel[2] = 0.5*( be_part[i].vel[2] + ( vminus2[N-1][i][2] + h*be_part[i].accel[2]) );
      
      be_part[i].pos[0] = 0.5*( be_part[i].pos[0] + ( rminus2[N-1][i][0] + h*be_part[i].vel[0]) );
      be_part[i].pos[1] = 0.5*( be_part[i].pos[1] + ( rminus2[N-1][i][1] + h*be_part[i].vel[1]) );
      be_part[i].pos[2] = 0.5*( be_part[i].pos[2] + ( rminus2[N-1][i][2] + h*be_part[i].vel[2]) );
      
      be_poly_vel[k][i][0] = be_part[i].vel[0]; 
      be_poly_vel[k][i][1] = be_part[i].vel[1]; 
      be_poly_vel[k][i][2] = be_part[i].vel[2]; 
      
      be_poly_pos[k][i][0] = be_part[i].pos[0]; 
      be_poly_pos[k][i][1] = be_part[i].pos[1]; 
      be_poly_pos[k][i][2] = be_part[i].pos[2]; 
    }

  return 0;
}

//-------------------------------------------------------------------------//
//                   Calcula aceleraciones debidas a los planetas
//-------------------------------------------------------------------------//

int Compute_acceleration(int i)
{
  int j;
  double r;
  
  //printf("i =%d\n",i);
  be_part[i].accel[0] = 0.0;
  be_part[i].accel[1] = 0.0;
  be_part[i].accel[2] = 0.0;
  //printf("Npart = %d N_part_ext = %d\n",Npart,Npart_ext);

  for(j=0; j<n_bodies; j++)
    {
      r = sqrt ( pow( (be_part[i].pos[0] - plansys[j].pos[0]) ,2) + 
		 pow( (be_part[i].pos[1] - plansys[j].pos[1]) ,2) + 
		 pow( (be_part[i].pos[2] - plansys[j].pos[2]) ,2) );
      
      be_part[i].accel[0] = be_part[i].accel[0] - plansys[j].mass * ( ( be_part[i].pos[0] - plansys[j].pos[0] ) / pow(r,3) );
      be_part[i].accel[1] = be_part[i].accel[1] - plansys[j].mass * ( ( be_part[i].pos[1] - plansys[j].pos[1] ) / pow(r,3) );
      be_part[i].accel[2] = be_part[i].accel[2] - plansys[j].mass * ( ( be_part[i].pos[2] - plansys[j].pos[2] ) / pow(r,3) );
    }

  return 0;
}

#include <cmath>
#include <iostream>
#include <ctime>
#include <fstream>
#include <iostream>

#define R_T 1.
#define dt 0.01
#define boxdim 8.

const int grid=200;
const int num_pts_mesh=grid*grid*grid;

double LT = 18; //! Local Time (=Orbit Position)
double SSL = M_PI/2.; //! SubSolarLatitude [rad]: the theta coordinate of subsolar point in TIIS, depends on Saturns orbit status
double SLT = -M_PI/12.*LT+M_PI*0.5; //! SubsolarLongiTude, the phi coordinate of the subsolar point in TIIS, depends on Titan local time
int flyby = 99;

double origin[3]={0.501*boxdim,0.501*boxdim,0.501*boxdim};
double x[grid],y[grid],z[grid];
double pos_array[3*num_pts_mesh];
float suv[grid][grid][grid];
float suv_array[num_pts_mesh];
// double n0=3.21e+22; //! 1/m^3 1e+22 backes N2
// double H=66.3e+3; //! 0.014 37km normalized to R_T for N2
// double n0=1.06e+19; //! 1/m^3 CH4 
// double H=90.6e+3;  //! CH4 0.0291
// double sigma=19.4035*1.e-22; //! 1e-21 schunk&nagy N2 abs
double sigma=25.7923*1.e-22; //! 2e-21 schunk&nagy CH4 abs
// double sigma=4.00475*1.e-22; //! 1e-22 schunk&nagy H2 abs
double I_inf=1.42059*1.e+13; //! 1e+13
double r;
double vec_x[3];
double sun[3]={sin(SSL)*cos(SLT),sin(SSL)*sin(SLT) ,cos(SSL) };
//! sigma_i * 10^-22 for m^2, cross-sections for N2
double sigma_N2_abs[37]={0.72,2.261,4.958,8.392,10.210,10.900,10.493,11.670,11.700,13.857,16.910,16.395,21.675,23.160,23.471,24.501,24.130,22.400,22.787,22.790,23.370,23.339,31.755,26.540,24.662,120.490,14.180,16.487,33.578,16.992,20.249,9.680,2.240,50.988,0.,0.,0.};
double sigma_N2_ion[37]={0.443,1.479,3.153,5.226,6.781,8.100,7.347,9.180,9.210,11.600,15.350,14.669,20.692,22.100,22.772,24.468,24.130,22.400,22.787,22.790,23.370,23.339,29.235,25.480,15.060,65.800,8.500,8.860,14.274,0.,0.,0.,0.,0.,0.,0.,0.};
//! sigma_i *10^-22 for m^2, cross-sections for CH4
double sigma_CH4_abs[37]={0.204,0.593,1.496,2.794,3.857,5.053,4.360,6.033,6.059,7.829,10.165,9.776,14.701,18.770,21.449,24.644,27.924,31.052,30.697,33.178,35.276,34.990,39.280,41.069,42.927,45.458,45.716,46.472,45.921,48.327,48.968,48.001,41.154,38.192,32.700,30.121,29.108};
double sigma_CH4_ion[37]={0.051,0.147,0.387,0.839,1.192,1.681,1.398,2.095,2.103,2.957,3.972,3.820,6.255,8.442,9.837,11.432,13.398,14.801,14.640,15.734,17.102,16.883,19.261,20.222,21.314,22.599,22.763,23.198,22.886,25.607,24.233,13.863,0.136,0.475,0.,0.,0.};
//! sigma_i *10^-22 for m^2, cross-sections for H2
double sigma_H2_abs[37]={0.0108,0.0798,0.2085,0.4333,0.6037,0.8388,0.7296,1.0180,1.0220,1.4170,1.9420,1.9010,3.0250,3.8700,4.5020,5.3560,6.1680,7.0210,6.8640,7.8110,8.4640,8.4450,9.9000,10.7310,11.3720,10.7550,8.6400,7.3390,8.7480,8.2530,0.4763,0.1853,0.,0.0456,0.,0.,0.};
double sigma_H2_ion[37]={0.0097,0.0758,0.2009,0.4028,0.5509,0.7454,0.6538,0.8999,0.9041,1.2960,1.7840,1.7420,2.8900,3.7780,4.0470,5.2540,6.0500,6.9000,6.7410,7.6680,8.2990,8.2880,9.7020,10.7310,9.7610,8.6240,7.0710,5.0720,6.6290,0.0899,0.,0.,0.,0.,0.,0.,0.};
double dlambda[37]={50, 50, 50, 50, 1, 1, 50, 1, 1, 50, 1, 50, 50, 1, 50, 50, 1, 1, 50,1, 1, 50, 50, 1, 50, 1, 1, 1, 50, 50, 50, 50, 1, 50, 1, 1, 50};


//! I_Inf_i * 10^13 for SI
double I_inf_i[37]={1.2,0.45,4.8,3.1,0.46,0.21,1.679,0.8,6.9,0.965,0.650,0.314,0.383,0.29,0.285,0.452,0.72,1.27,0.357,0.53,1.59,0.342,0.32,0.36,0.141,0.17,0.26,0.702,0.758,1.625,3.537,3,4.4,1.475,3.5,2.1,2.467};
double A[37]={1.0017e-2,7.1250e-3,1.3375e-2,1.9450e-2,2.7750e-3,1.3768e-1,2.6467e-2,2.5000e-2,3.3333e-3,2.2450e-2,6.5917e-3,3.6542e-2,7.4082e-3,7.4917e-3,2.0225e-2,8.7583e-3,3.2667e-3,5.1583e-3,3.6583e-3,1.6175e-2,3.3250e-3,1.18e-2,4.2667e-3,3.0417e-3,4.7500e-3,3.8500e-3,1.2808e-2,3.2750e-3,4.7667e-3,4.8167e-3,5.6750e-3,4.9833e-3,3.9417e-3,4.4168e-3,5.1833e-3,5.2833e-3,4.3750e-3};
double F107, F107A, flyby_dist_sun_saturn;


double R_cut=8.; //! R_t is not included here
double max_rate=0;
double pos_max[3]={0.,0.,0.};


int count_neg_value=0;
int count_beneath_surface=0;
int count_shadow=0;
int count_point_ionosphere=0;
int count_progress=0;

using namespace std;



void calc_solar_flux_flyby()
{
  switch(flyby)
    {//! Quick standard secnario
      case 0: 
	F107=70.7; //! daily flux
	F107A=75.29; //! averaged flux +-1 month
	flyby_dist_sun_saturn=9.469 ;
	break;
      //! T9 on 26.12.2005, Distant Wakeflyby
      case 9: 
	F107=89.5; //! daily flux
	F107A=85.39; //! averaged flux +-1 month
	flyby_dist_sun_saturn=9.106;  //!AE
	LT=21.;
	SSL=-(-19.06*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);
	break;
      //! T63 on 12.12.2009, Distand Wakeflyby CA 01:03
      case 63:
	F107=70.7; //! daily flux
	F107A=75.29; //! averaged flux +-1 month
	flyby_dist_sun_saturn=9.469;  //!AE
	LT=17.;
	SSL=-(2.04*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);
	break;
	//! T70 on 21.06.2010
      case 70: 
	F107=72.2; //! daily flux
	F107A=78.00; //! averaged flux +-1 month  NOTE this is approximated
	flyby_dist_sun_saturn=9.527;  //!AE
	LT=16.1;
	SSL=-(4.913*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);	
	//! T75 on 19.4.2011, Distant Wakeflyby CA 05:00
      case 75:
	F107=113.9; //! daily flux
	F107A=106.68; //! averaged flux +-1 month
	flyby_dist_sun_saturn=9.617;  //!AE
	LT=14.2;
	SSL=-(9.273*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);
	break;
      case 99: //New-Horizons Pluto encounter on 14.7.2015
	F107=113.9; //! daily flux
	F107A=106.68; //! averaged flux +-1 month
	flyby_dist_sun_saturn=32.9;  //!AE
	LT=18;
	SSL=-(0*M_PI)/180.+M_PI/2.;
	SLT=-M_PI/12.*LT+M_PI*0.5;
	sun[0]=sin(SSL)*cos(SLT);
	sun[1]=sin(SSL)*sin(SLT);
	sun[2]=cos(SSL);
	break;
	
    }
    
  double P=(F107+F107A)/2.;
  for(int i=0;i<37;i++)
  {
    I_inf_i[i]*=(1+A[i]*(P-80))/(flyby_dist_sun_saturn*flyby_dist_sun_saturn);   
  }
  
  
  
}


//! Neutral N2 profile
double N2(double z, bool &in_atmosphere, double *z_vec)
{
  if((abs(z_vec[0])<=R_cut && abs(z_vec[1])<=R_cut && abs(z_vec[2])<=R_cut))
//   {
//    if(z<=R_cut)
//    return 3.21e+22*exp(-2575*(z-R_T)/66.3); //!fixed scale height N2, and CH4, use respective parameters H,n0
//    return 1.5e+20*exp(-1184*(z-R_T)/71.5);
      return 1.e+20*exp(-1184*(z-R_T)/82.4)+1.e+16*exp(-1184*(z-1.5)/54.7)+1.e+13*exp(-1184*(z-2)/318.5);
  else
  {
      in_atmosphere=false;
    //return 1.5e+20*exp(-1184*(z-R_T)/71.5);
      return 1.e+20*exp(-1184*(z-R_T)/82.4)+1.e+16*exp(-1184*(z-1.5)/54.7)+1.e+13*exp(-1184*(z-2)/318.5);
  }
}

//! Neutral CH4 profile
double CH4(double z, bool &in_atmosphere, double *z_vec)
{
  if((abs(z_vec[0])<=R_cut && abs(z_vec[1])<=R_cut && abs(z_vec[2])<=R_cut))
//   {
//    if(z<=R_cut)
//    return 1.06e+19*exp(-2575*(z-R_T)/90.6); //!fixed scale height N2, and CH4, use respective parameters H,n0
//    return 1.5e+17*exp(-1184*(z-R_T)/150.6); 
      return 1.e+18*exp(-1184*(z-R_T)/88.3)+7.64e+13*exp(-1184*(z-1.5)/161.5)+9.47e+12*exp(-1184*(z-2)/509.4);
 else
  {
    in_atmosphere=false;
//    return 1.5e+17*exp(-1184*(z-R_T)/150.6); 
      return 1.e+18*exp(-1184*(z-R_T)/88.3)+7.64e+13*exp(-1184*(z-1.5)/161.5)+9.47e+12*exp(-1184*(z-2)/509.4);
  }
}

//! Neutral H2 profile
double H2(double z, bool &in_atmosphere, double *z_vec)
{
  if((abs(z_vec[0])<=R_cut && abs(z_vec[1])<=R_cut && abs(z_vec[2])<=R_cut))
//   {
//   if(z<=R_cut)
  {
    if(z<=R_T+0.4757)
// 	if(z<=R_T+0.6)
      return 5.82e+19*exp(-2575*(z-R_T)/69.8); //! H1=69.8 km
// 	  return 0.;
    else
      return 2.52439e+13*exp(-2575*(z-R_T)/297.6);  //!H2=297,6km
//       return 2.564e+12*exp(-2575*(z-R_T)/2000.6);  
//       return 3.1443e+12*exp(-2575*(z-R_T)/1500.6);  
//       return 4.72807e+12*exp(-2575*(z-R_T)/1000.6);  
//       return 1.60598e+13*exp(-2575*(z-R_T)/500.6);  

  }
  else
  {
    in_atmosphere=false;
    return 2.52439e+13*exp(-2575*(z-R_T)/297.6);
//     return 3.1443e+12*exp(-2575*(z-R_T)/1500.6);
//     return 4.72807e+12*exp(-2575*(z-R_T)/1000.6); 
//        return 1.60598e+13*exp(-2575*(z-R_T)/500.6);  

  }
}

double calc_photo_rate_trapezrule(double *z_vec, double z0, int i, int j, int k)
{
  double not_in_shadow=1.;
  long double Prod=0.;
  bool do_integrate=true;
  double N2x0=N2(z0,do_integrate,z_vec);
  double CH4x0=CH4(z0,do_integrate,z_vec);
    //! We get the starting height and the grid coordinates
  //! We need to approximate the integral by many trapezes
  double temp=0.;
  double temp1=0.;
  double temp2=0.;
  
  
  double a=z0;
  z_vec[0]+=dt*sun[0];
  z_vec[1]+=dt*sun[1];
  z_vec[2]+=dt*sun[2];
  double b=sqrt(z_vec[0]*z_vec[0]+z_vec[1]*z_vec[1]+z_vec[2]*z_vec[2]);
  do
  {
    if(a>=R_T && b>=R_T)
    {
//       temp+=dt/2*(N(a,z_vec,in_atmosphere)+N(b,z_vec,in_atmosphere));
      temp1+=dt/2*(N2(a,do_integrate,z_vec)+N2(b,do_integrate,z_vec));
      temp2+=dt/2*(CH4(a,do_integrate,z_vec)+CH4(b,do_integrate,z_vec));
      
      z_vec[0]+=dt*sun[0];
      z_vec[1]+=dt*sun[1];
      z_vec[2]+=dt*sun[2];
      a=b;
      b=sqrt(z_vec[0]*z_vec[0]+z_vec[1]*z_vec[1]+z_vec[2]*z_vec[2]);
    }
    else
    {
      count_shadow+=1;
      not_in_shadow=0.;
      break;
    }
  }
  while(do_integrate);
//   tau[i][j][k]=2575000*temp*sigma;
//   cout<<"Integral was "<<temp<<endl;
//   temp=pow(M_E,-n0*temp);
  //! Complete absorption + EUVAC wavelength model
  for(int i=1;i<=37;i++)
   // Prod+=sigma_N2_ion[i]*1.e-22*I_inf_i[i]*1.e+13*exp(-1.e-22*1184000*(sigma_N2_abs[i]*temp1+sigma_CH4_abs[i]*temp2));
     Prod+=sigma_CH4_ion[i]*1.e-22*I_inf_i[i]*1.e+13*exp(-1.e-22*1184000*(sigma_N2_abs[i]*temp1+sigma_CH4_abs[i]*temp2));
//     Prod+=dlambda[i]*sigma_H2_ion[i]*1.e-22*I_inf_i[i]*1.e+13*exp(-1.e-22*2575000*(sigma_N2_abs[i]*temp1+sigma_CH4_abs[i]*temp2+sigma_H2_abs[i]*temp3));
  
  if(z0<1.5)
    return not_in_shadow*CH4x0*Prod+100;
  else
    return not_in_shadow*CH4x0*Prod;
    
  //!For Baumjohann comparison without 37 bin model#
  
//   P=exp(-2575000*temp2*sigma); //! => e^(-tau)
//   return not_in_shadow*sigma*I_inf/(9.582*9.582)*CH4x0*P+(1.-not_in_shadow);
  //! Combined two species
//   tau[i][j][k]=2575000*temp;
//   P=exp(-2575000*temp);
//   return not_in_shadow*I_inf*P/(9.582*9.582)*(1.e-21*N2x0+2.e-21*CH4x0+1e-22*H2x0)+(1.-not_in_shadow);
//   return not_in_shadow*I_inf/(9.582*9.582)*(1.e-21*N2x0*P)+(1.-not_in_shadow);
//   return not_in_shadow*I_inf/(9.582*9.582)*(2.e-21*CH4x0*P)+(1.-not_in_shadow);
  //!without bin model, complete absorption
//   P=exp(-2575000*1.e-22*(19.4035*temp1+25.7923*temp2+4.00475*temp3));
//   return not_in_shadow*I_inf/(9.582*9.582)*1.e-22*(19.4035*N2x0+25.7923*CH4x0+4.00475*H2x0)*P+(1.-not_in_shadow);
//     P=exp(-2575000*1.e-22*(19.4035*temp1+25.7923*temp2+4.00475*temp3));
//   return not_in_shadow*I_inf/(9.582*9.582)*1.e-22*(/*19.4035*N2x0+*//*25.7923*CH4x0+*/4.00475*H2x0)*P+(1.-not_in_shadow);
}


int main()
{
  clock_t start,finish;
  double time;
  start = clock();
  double zvec[3];
  //!Fill grid with coordinates, Box with 20^3 R_T, Origin in center, This is the global mesh
  for(int i=0;i<grid;i++)
  {
    x[i]=-origin[0]+i*boxdim/(grid-1);
    y[i]=-origin[1]+i*boxdim/(grid-1);
    z[i]=-origin[2]+i*boxdim/(grid-1);
  }
  int i_j_k;
  
  calc_solar_flux_flyby();
  
  //! Calculate SUV rate at every grid point
  for(int i=0;i<grid;i++)
    for(int j=0;j<grid;j++)
      for(int k=0;k<grid;k++)
      {
	count_progress++;

	
	zvec[0]=x[i];zvec[1]=y[j];zvec[2]=z[k]; //! this vector has global cartesic coordinates
	r=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);  //! R_T not included meaning this is NOT altitude above ground
	if(r>=R_T)//! This guarantees that our starting point is at least on the surface
	{
	  i_j_k=i*grid*grid+j*grid+k;	  
	  //!SUV-rate function call
	  suv[i][j][k]=calc_photo_rate_trapezrule(zvec,r,i,j,k); 
// 	  pos_array[0*num_pts_mesh+i_j_k]=x[i];
// 	  pos_array[1*num_pts_mesh+i_j_k]=y[j];
// 	  pos_array[2*num_pts_mesh+i_j_k]=z[k];
	  if(suv[i][j][k]<0)count_neg_value++; //!negative values are wrong, check the code
	  suv_array[i_j_k]=suv[i][j][k];
	  //!Global Maxium determination
	  if(suv[i][j][k]>max_rate)
	  {
	    max_rate=suv[i][j][k];
	    pos_max[0]=x[i];
	    pos_max[1]=y[j];
	    pos_max[2]=z[k];
	  }

	}
	else //! obviously we are in the planet with our point here, no ionization rate here
	{
	  count_beneath_surface++;
	  suv[i][j][k]=0.; //! Value of 1. for better Visualization
	}
		
	if(count_progress%(num_pts_mesh/100)==0)
	  cout<<"\rDone with "<<100*count_progress/num_pts_mesh<<" % of total mesh points"<<flush;

      }
  count_point_ionosphere=grid*grid*grid-count_beneath_surface;
  cout<<"\nNumber points beneath Titans surface is: "<<count_beneath_surface<<endl;
  cout<<"Number points in Titans shadow is: "<<count_shadow<<endl;
  cout<<"Resolution: "<<boxdim*1184/grid<<"km"<<endl;
  cout<<"Max value at grid: "<<max_rate<<endl;
  cout<<"Position of max value: R: "<<sqrt(pos_max[0]*pos_max[0]+pos_max[1]*pos_max[1]+pos_max[2]*pos_max[2])<<", theta: "<<acos(pos_max[2]/sqrt(pos_max[0]*pos_max[0]+pos_max[1]*pos_max[1]+pos_max[2]*pos_max[2]))<<", phi: "<<atan(pos_max[1]/pos_max[0])<<endl;
  cout<<"neg values : "<<count_neg_value<<endl;
 
//!Output stuff, if grid is changed, output must also be changed for same data
  //!Output z=0 ebene 
//   ofstream out;
//   out.open("N2-3d.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//       for(int k=0;k<grid;k++)
// 	out<<x[i]<<" "<<y[j]<<" "<<z[k]<<" "<<suv[i][j][k]<<endl;
//   out.close();
  
  //! 3D output in txt, length renormed to meters
//   ofstream outSI;
//   outSI.open("H2-3dSI.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//       for(int k=0;k<grid;k++)
// 	outSI<<x[i]*2575000<<" "<<y[j]*2575000<<" "<<z[k]*2575000<<" "<<suv[i][j][k]<<endl;
//   outSI.close();
//   
//   
// //! 2d output
//   ofstream out_2d;
//   out_2d.open("H2-2dSI.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
// 	out_2d<<x[i]*2575000<<" "<<y[j]*2575000<<" "<<suv[i][j][grid/2]<<endl;
//   out_2d.close();
  
//!Binary output, for ptracer usage, The writing procedure must exactly match the one to be read in CHybrid_IonProfiles read_extern_IonProfile_pTracer()
//! Write header for 
  char Tracer_Info[42]="test";
  short num_extern_Nodes_short[3]={grid,grid,grid};
  int TL=0;
	float extern_Origin_FILEREAL[3]={static_cast<float>(origin[0]),static_cast<float>(origin[1]),static_cast<float>(origin[2])};
  float extern_Length_FILEREAL[3]={boxdim,boxdim,boxdim};
  float Radius=1;
  float SI_Quantities[4]={1,1,1,1};//!
  
  short num_comps=1;
  int ftype=0;
  char Field_Name[50]="ion_prod_rate";
//   float field_array_box=

  
  ofstream out_bin;
  out_bin.open("CH4_Pluto_woffset_fine",ofstream::binary);
  out_bin.write( reinterpret_cast<char*> (Tracer_Info),42*sizeof(char));
  out_bin.write(reinterpret_cast<char*> (num_extern_Nodes_short), 3*sizeof(short));
  out_bin.write(reinterpret_cast<char*> (&TL), sizeof(int));
  out_bin.write(reinterpret_cast<char*> (extern_Origin_FILEREAL), 3*sizeof(float));
  out_bin.write(reinterpret_cast<char*> (extern_Length_FILEREAL), 3*sizeof(float));
  out_bin.write(reinterpret_cast<char*> (&Radius), sizeof(float));
  out_bin.write(reinterpret_cast<char*> (SI_Quantities), 4*sizeof(float));
  
//   out_bin.write(reinterpret_cast<char*> (pos_array),num_pts_mesh*3*sizeof(float));

  out_bin.write( reinterpret_cast<char*> (&num_comps),sizeof(short));  //!Number of components of the field (1 or 3)
  out_bin.write( reinterpret_cast<char*> (&ftype),sizeof(int));  	//!field type 9,1,2,3 isnt even used after read?
  out_bin.write( reinterpret_cast<char*> (Field_Name), 50*sizeof(char));  //! also not used after read
  out_bin.write( reinterpret_cast<char*> (suv_array), num_comps*num_pts_mesh*sizeof(float));
  
  
  out_bin.close();

//! Output along single line (parallel to grid)
  ofstream lineoutx;
  lineoutx.open("CH4prod-x.txt");
  for(int i=0;i<grid;i++)
//       if(sqrt(x[i]*x[i]+y[j]*y[j])>=1.9 && sqrt(x[i]*x[i]+y[j]*y[j])<=2.1)
      lineoutx<<x[i]<<" "<<suv[i][grid/2][grid/2]<<endl;
  lineoutx.close();
//   
//   ofstream lineouty;
//   lineoutx.open("H2-lineouty.txt");
//   for(int i=0;i<grid;i++)
// //       if(sqrt(x[i]*x[i]+y[j]*y[j])>=1.9 && sqrt(x[i]*x[i]+y[j]*y[j])<=2.1)
//       lineoutx<<x[i]<<" "<<suv[grid/2][i][grid/2]<<endl;
//   lineoutx.close();
  
//! Output cos_chi for checking
//   ofstream out_chi;
//   out_chi.open("cos_chi.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//     {
//       vec_x[0]=x[i];vec_x[1]=y[j];vec_x[2]=0.;
//       out_chi<<x[i]<<" "<<y[j]<<" "<<cos_chi(vec_x,sun)<<endl;
//     }
//   out_chi.close(); 
  
//   ofstream tau_out;
//   tau_out.open("tau.txt");
//   for(int i=0;i<grid;i++)
//     for(int j=0;j<grid;j++)
//     {
//       vec_x[0]=x[i];vec_x[1]=y[j];vec_x[2]=0.;
//       tau_out<<x[i]<<" "<<y[j]<<" "<<tau[i][j][150]<<endl;
//     }
//   tau_out.close(); 
//   (x[i]*sun[1]+y[j]*sun[2]+a[3]*b[3])/sqrt((a[1]*a[1]+a[2]*a[2]+a[3]*a[3])*(b[1]*b[1]+b[2]*b[2]+b[3]*b[3]))
//   ofstream test;
//   test.open("dat4.txt");
// //       if(sqrt(x[i]*x[i]+y[j]*y[j])>=1.9 && sqrt(x[i]*x[i]+y[j]*y[j])<=2.1)
//       test<<x[40+16]<<" "<<z[40+16]<<" "<<suv2[40+16][0][40+16]<<endl;
//       test<<x[40+18]<<" "<<z[40+14]<<" "<<suv2[40+18][0][40+14]<<endl;
//   test.close();
  
  cout<<"done"<<endl;   
  finish = clock();
  time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  cout<< " Time: " << time << "s." << endl << endl;   
  
  
  
}

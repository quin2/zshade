#include <FPT.h>
#include "M3d_matrix_tools.c"
#include <math.h>
#include <stdlib.h>
#include <float.h>

//time the program
#include <time.h>

//define window dimensions
#define X_DEVICE 800
#define Y_DEVICE 800

//step sizes for renderer
double du = 0.006;
double dv = 0.006;

//step size for normals
double small = 0.001;

//hither line for clipping
double hither = 0.01;

//define z buffer
double zbuf[X_DEVICE][Y_DEVICE];

//clear buffer
int clearBuffer(){
	for(int i = 0; i < X_DEVICE; i++){
		for(int j = 0; j < Y_DEVICE; j++){
			zbuf[i][j] = DBL_MAX; //set to max value
		}
	}
	return 0;
}

//sphere drawing function
int sphere(double u, double v, double xyz[3]){
	xyz[0] = cos(u) *cos(v);
	xyz[1] = sin(v);
	xyz[2] = sin(u) * cos(v);

	return 0;
}

int cyl(double u, double v, double xyz[3]){
	xyz[0] = 0.25 * cos(u);
	xyz[1] = 0.25 * sin(u);
	xyz[2] = v;

	return 0;
}

int taurus(double u, double v, double xyz[3]){
	double inner = 0.75;
	double outer = 1.0;

	xyz[0] = cos(u) * (outer + (inner * cos(v)));
	xyz[1] = sin(u) * (outer + (inner * cos(v)));
	xyz[2] = sin(v);
	return 0;
}

// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  light_in_eye_space[0] = 100; light_in_eye_space[1] = 200; light_in_eye_space[2] = 0;

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}

int zPlot(double ulo, double uhi, double vlo, double vhi, int (*func)(double u, double v, double xyz[3]), double mat[4][4], double color[3]){
	double H = X_DEVICE / 2;	//half of device
  	double half = 3.14159 / 4;	//half angle
  	double tanhalf = tan(half);
  	double c = H / tanhalf; //half of device over tan(h)

	for(double u = ulo; u <= uhi; u += du){
		for(double v = vlo; v <= vhi; v += dv){
			double p[3];
			func(u, v, p);
			double q[3];
			func(u + small, v, q);
			double r[3];
			func(u, v + small, r);

			M3d_mat_mult_pt(p, mat, p);
			M3d_mat_mult_pt(q, mat, q);
			M3d_mat_mult_pt(r, mat, r);

			//clip points
			//fabs(p[1]/p[2]) > hither && fabs(p[0]/p[2]) > hither && p[0] > hither
			//clip algo is a little screwy so I'm bypassing it now!
			if(p[2] > hither && fabs(p[1]/p[2]) <= tanhalf && fabs(p[0]/p[2]) <= tanhalf){

				//display for screen
				double ybar = p[1] / p[2];
				double xbar = p[0] / p[2];

				double xbarbar = (c * xbar) + H;
				double ybarbar = (c * ybar) + H;

			    int xbb = (int)xbarbar;
				int ybb = (int)ybarbar;
				
				if(p[2] < zbuf[xbb][ybb]){
					zbuf[xbb][ybb] = p[2];

					//calc normals
					q[0] = p[0] - q[0]; q[1] = p[1] - q[1]; q[2] = p[2] - q[2];
					r[0] = p[0] - r[0]; r[1] = p[1] - r[1]; r[2] = p[2] - r[2];

					double f[3];
					f[0] = (q[1] * r[2]) - (q[2] * r[1]);
					f[1] = (q[0] * r[2]) - (q[2] * r[0]);
					f[2] = (q[0] * r[1]) - (q[1] * r[0]);

					//turn into unit
					double flen = sqrt((f[0] * f[0]) + (f[1] * f[1]) + (f[2] * f[2]));
					f[0] /= flen; f[1] /= flen; f[2] /= flen;

					double eye[3];
					eye[0] = 0 - p[0]; eye[1] = 0 - p[1]; eye[2] = 0 - p[2];
					double output[3];

					//find color here!

					Light_Model(color, eye, p, f, output);

					G_rgb(output[0], output[1], output[2]);
					G_point(xbarbar, ybarbar);
				} 
				
			} 

		}
	}

	return 0;
}

//all scene graph specific code goes here
int drawScene(double V[4][4]){
	//init transformation matrix
   	double color[3];
   	double m[4][4], minv[4][4] ;
	int tlist[100] ;
	double plist[100] ;
	int num;

	//get rid of this shape for now
	//eye in center 
	/*
	num = 0;
	color[0] = 0.5; color[1] = 0.5; color[2] = 0.0;
	zPlot(0, 2 * M_PI, -M_PI/2, M_PI/2, sphere, V, color);
	*/

	//taurus
	num = 0;
	color[0] = 0.5; color[1] = 0.0; color[2] = 0.0;
	zPlot(0, 2 * M_PI, 0, 2 * M_PI, taurus, V, color);

	//boutta make a movie independant 
	//need new trucks independant
	return 1;
}

int main(){
	//init stuff here
	G_init_graphics(X_DEVICE, Y_DEVICE);
	G_rgb(0,0,0);
   	G_clear();

   	//clear buffer
   	clearBuffer();	

	double eye[3]; double coi[3]; double up[3];

    int numFrames = 100;
    int fnum = 0;
    float t = 0.0;
    double V[4][4], Vi[4][4];
    int q = 0;

    clock_t start, end;
    double cpu_time_used;

    start = clock();

	while(fnum < numFrames){ 
	    t = 1.0 * fnum/numFrames;

	    eye[0] = 15 * cos(2*M_PI*t);
	    eye[1] = 6 * t;
	    eye[2] = 7 * sin(2*M_PI*t);

	    coi[0] = 0; coi[1] = 0; coi[2] = 0;

	    up[0] = eye[0]; up[1] = eye[1] + 1; up[2] = eye[2];

	    M3d_view(V, Vi, eye, coi, up);

	    drawScene(V);

	    //movie related code
	    /*
	    char str[20];
	    sprintf(str, "pic%04d.xwd", fnum);
	    G_save_image_to_file(str);
	    */

	    //q = G_wait_key(); if (q == 'q') { break; }
	    G_display_image();

	    G_rgb(0,0,0);
	    G_clear();
	    clearBuffer();	

	    //drawScene(V);

	    fnum++;
  	}

  	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Code took %lf sec to run,\n", cpu_time_used);
    printf("Averaging %lf FPS\n", (double)numFrames/cpu_time_used);

    return 0;
}
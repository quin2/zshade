static float3 acc_lightmodel (float3 irgb,
                 float3 s,
                 float3 p,
                 float3 n,
                 float AMBIENT,
                 float MAX_DIFFUSE,
                 float SPECPOW)
{
  float3 argb = (float3)(0.0, 0.0, 0.0);

  float3 light_in_eye_space ;
  light_in_eye_space.x = 100; light_in_eye_space.y = 200; light_in_eye_space.z = 0;

  float len ;
  float3 N ; 
  len = length(n) ;
  if (len == 0) return argb ;
  N = normalize(n) ;

  float3 E ;
  E = s - p;
  len = length(E);
  if (len == 0) return argb ;
  E = normalize(E);
  float NdotE = dot(N,E) ;

  float3 L ;
  L = light_in_eye_space - p;
  len = length(L) ;
  if (len == 0) return argb ;
  L = normalize(L) ;
  float NdotL = dot(N, L) ;

  float max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;

  //is goto safe kernal syntax?
  float intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N *= (float3)(-1.0);
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }

  // ignore Blinn's variant
  float3 R ; // Reflection vector of incoming light
  R = (2 * NdotL * N) - L;

  float EdotR = dot(E, R) ;

  float diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  float specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  intensity = AMBIENT + diffuse + specular ;



  LLL : ;

  float f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb = f * irgb;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb = g * irgb + f;
  }

  return argb;
}

//possibly make an object type switch in here? but then the range would be fucked up...
static float3 taurus(float u, float v){
  float3 xyz;

  float inner = 0.75;
  float outer = 1.0;

  xyz.x = cos(u) * (outer + (inner * cos(v)));
  xyz.y = sin(u) * (outer + (inner * cos(v)));
  xyz.z = sin(v);
  return xyz;
}

static float3 sphere(double u, double v){
  float3 xyz;

  xyz.x = cos(u) *cos(v);
  xyz.y = sin(v);
  xyz.z = sin(u) * cos(v);

  return xyz;
}

static float3 matmult(float3 in, float16 mat){
	float3 out;

	out.x = mat.s0*in.x + mat.s1*in.y + mat.s2*in.z + mat.s3;
	out.y = mat.s4*in.x + mat.s5*in.y + mat.s6*in.z + mat.s7;
	out.z = mat.s8*in.x + mat.s9*in.y + mat.sA*in.z + mat.sB;

	return out; 
}

__kernel void zplot(float du,
	float dv,
	float H,
  int imW,
  int shape, 
  float halfA,
  float3 color,
	float16 mat,
	float AMBIENT,
	float MAX_DIFFUSE,
	int SPECPOW,
	__global float * buf,
  __write_only image2d_t out)
{
  float small = 0.001;

  float tanhalf = (float)tan(halfA);
  float c = H / tanhalf; 
    
  //delete du and dv if it causes too much trouble
  float u = (float)get_global_id(0) * du;
  float v = (float)get_global_id(1) * dv;

  //remove mult by du and dv to get original...
  //therotically better to use seperate kernals (we're cheating, fuck rewritng everything)
  float3 p; float3 q; float3 r;
  if(shape == 1){
    p = sphere(u, v);
    q = sphere(u + small, v);
    r = sphere(u, v + small);
  }
  else if(shape == 0){
    p = taurus(u, v);
    q = taurus(u + small, v);
    r = taurus(u, v + small);
  }
  

  p = matmult(p, mat);
  q = matmult(q, mat);
  r = matmult(r, mat);


  int st1 = 1; //(fabs(p.x/p.z), tanhalf)
  int st2 = islessequal(fabs(p.y/p.z), tanhalf);
  int st3 = islessequal(fabs(p.x/p.z), tanhalf);

  if(st1 && st2 && st3){
    //display for screen
    float ybar = p.y / p.z;
    float xbar = p.x / p.z;

    float xbarbar = (c * xbar) + H;
    float ybarbar = (c * ybar) + H;

    int xbb = (int)xbarbar;
    int ybb = (int)ybarbar;

    //insert fancy accessing math here
    int index = (ybb * imW) + xbb;
    float cp = buf[index];
    if(isless(p.z, cp)){
      //this is effecent lol...only look at ambient!!!
      buf[index] = (float)p.z;

      //calc normals 
      q = p - q;
      r = p - r;

      float3 f = cross(q, r);

      //turn into unit
      f = normalize(f);

      float3 eye = (float3)0.0 - p;

      //find color here!
      float3 output = acc_lightmodel(color, eye, p, f, AMBIENT, MAX_DIFFUSE, SPECPOW);

      //set output image stuff here!!!
      write_imagef(out, (int2)(xbb, ybb), (float4)(output, 0.0));
    } 
  }
  return;
}


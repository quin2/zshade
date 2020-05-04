#include <FPT.h>
#include "M3d_matrix_tools.c"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <unistd.h>

//openCL helper lib for C
#include <cf4ocl2.h>

//handle error assertions
#include <assert.h>

//error handeling macros, borrowed from cf4ocl example project
#define ERROR_MSG_AND_EXIT(msg) \
    do { fprintf(stderr, "\n%s\n", msg); exit(EXIT_FAILURE); } while(0)
#define HANDLE_ERROR(err) \
    if (err != NULL) { ERROR_MSG_AND_EXIT(err->message); }

//define window dimensions
#define X_DEVICE 800
#define Y_DEVICE 800

//shape primitives
#define SPHERE 1
#define TAURUS 0

//step sizes for renderer
float du = 0.006;
float dv = 0.006;

//step size for normals
float small = 0.001;

//hither line for clipping
float hither = 0.01;

//vision parameters
float H = X_DEVICE / 2;	//half of device
float half = 3.14159 / 4;	//half angle

// To support the light model :
float light_in_eye_space[3] ;
float AMBIENT      = 0.2 ;
float MAX_DIFFUSE  = 0.5 ;
float SPECPOW      = 50 ;

//output array
float output[X_DEVICE][Y_DEVICE][4];

//openCL variables
CCLContext* ctx;
CCLDevice* dev;
CCLQueue* queue;
CCLProgram* prg;
CCLKernel* krnl;

//image format
cl_image_format image_format = { CL_RGBA, CL_FLOAT };

//Z buffer
cl_float* zbuf = NULL;

//setup context and build program for doughnut shader
int setup_acc(){
    CCLErr* err = NULL;

    //device selected here is my internal intel GFX card, you probably need to change this
    int dev_idx = 1;

    //create compute context from selected compute device
    ctx = ccl_context_new_from_menu_full(&dev_idx, &err);
    HANDLE_ERROR(err);

    //get specified device in context 
    dev = ccl_context_get_device(ctx, 0, &err);
    HANDLE_ERROR(err);

    //create command queue
    queue = ccl_queue_new(ctx, dev, 0, &err);
    HANDLE_ERROR(err);

    prg = ccl_program_new_from_source_file(ctx, "zshade.cl", &err);
    HANDLE_ERROR(err);

    const char * bldlog;
    CCLErr* err_bld = NULL;
    ccl_program_build(prg, NULL, &err);
    if ((err) && (err->code == CL_BUILD_PROGRAM_FAILURE)) {
        bldlog = ccl_program_get_build_log(prg, &err_bld);
        fprintf(stderr, "Error building program: \n%s", bldlog);
    }

    //get kernal wrapper
    krnl = ccl_program_get_kernel(prg, "zplot", &err);
    HANDLE_ERROR(err);

    return 0;
}

int close_acc(){
    ccl_program_destroy(prg);
    ccl_queue_destroy(queue);
    ccl_context_destroy(ctx);

    //check that all wrappers destoryed
    assert(ccl_wrapper_memcheck());

    return 0;
}

int clearBuf(){
    zbuf = (cl_float*) malloc(sizeof(cl_float) * X_DEVICE * Y_DEVICE);
    for(int i = 0; i < (X_DEVICE * Y_DEVICE); i++){
        zbuf[i] = CL_INFINITY; 
    }
    return 0;
}


int zPlot_acc(float ulo, float uhi, float vlo, float vhi, double mat[4][4], float color[3], int shape){
	//use step size and u/v int. to describe step size....
	int wX = (int)fabs((uhi - ulo) / du);
	int wY = (int)fabs((vhi - vlo) / dv);

	//still pass in du, dv
	//should be able to rehydrate these by doing u = x * du where x is loop position...
	//only really works in this case, otherwise u = (x * du) + ulo haha 
    //will need to change for cyl!

	//image
    CCLImage* img_out;

    //queue waitlist for z buffer
    CCLEventWaitList ewl = NULL;

    //error handeling object
    CCLErr* err = NULL;
    
    //make an image output var for us :)
    img_out = ccl_image_new(ctx, CL_MEM_WRITE_ONLY,
        &image_format, NULL, &err,
        "image_type", (cl_mem_object_type) CL_MEM_OBJECT_IMAGE2D,
        "image_width", (size_t) X_DEVICE,
        "image_height", (size_t) Y_DEVICE,
        NULL);
    HANDLE_ERROR(err);
    
    //estimate local and global worksizes for device
    size_t gws[2]; //global work size, don't touch this
    size_t lws[2]; //local work size, don't touch this
    size_t real_ws[2]; //real work size

    real_ws[0] = wX; //set these to work sizes that we defined
    real_ws[1] = wY;

    ccl_kernel_suggest_worksizes(krnl, dev, 2, real_ws, gws, lws, &err);
    HANDLE_ERROR(err);

    //variable translation
    cl_float3 color_k = (cl_float3){color[0], color[1], color[2]};
    cl_float16 mat_k = (cl_float16){mat[0][0], mat[0][1], mat[0][2], mat[0][3], 
    							  mat[1][0], mat[1][1], mat[1][2], mat[1][3],
    							  mat[2][0], mat[2][1], mat[2][2], mat[2][3],
    							  mat[3][0], mat[3][1], mat[3][2], mat[3][3]};

    cl_float du_s = (cl_float)uhi/wX;
    cl_float dv_s = (cl_float)vhi/wY;
    cl_float H_s = (cl_float)H;
    cl_int xd_s = (cl_int)X_DEVICE;
    cl_int shape_s = (cl_int)shape; 
    cl_float half_s = (cl_float)half;
    cl_float ambient_s = (cl_float)AMBIENT;
    cl_float diffuse_s = (cl_float)MAX_DIFFUSE;
    cl_int specpow_s = (cl_int)SPECPOW;

    // buffer buffer (no for real)
	CCLBuffer* zbufbuf;
	zbufbuf = ccl_buffer_new(ctx, CL_MEM_READ_WRITE, X_DEVICE * Y_DEVICE * sizeof(cl_float), NULL, &err);
	HANDLE_ERROR(err);

    //queue buffer write
    //we're ending up with 'invalid' values here.... no good
	CCLEvent* evt_write1;
    size_t zSize = X_DEVICE * Y_DEVICE * sizeof(cl_float);
	evt_write1 = ccl_buffer_enqueue_write(zbufbuf, queue, CL_FALSE, 0,
        zSize, zbuf, NULL, &err);
    HANDLE_ERROR(err);

    ccl_event_wait_list_add(&ewl, evt_write1, NULL);

    //move all the stuff to the shader
    ccl_kernel_set_args_and_enqueue_ndrange(
        krnl, queue, 2, NULL, gws, lws, NULL, &err,
        ccl_arg_priv(du_s, cl_float), 
        ccl_arg_priv(dv_s, cl_float), 
        ccl_arg_priv(H_s, cl_float),
        ccl_arg_priv(xd_s, cl_int), 
        ccl_arg_priv(shape_s, cl_int),
        ccl_arg_priv(half_s, cl_float),
        ccl_arg_priv(color_k, cl_float3),
        ccl_arg_priv(mat_k, cl_float16),
        ccl_arg_priv(ambient_s, cl_float),
        ccl_arg_priv(diffuse_s, cl_float),
        ccl_arg_priv(specpow_s, cl_int),
        zbufbuf,
        img_out, 
        NULL);
    HANDLE_ERROR(err);

    //wait until we're done
    ccl_enqueue_barrier(queue, &ewl, &err);
    HANDLE_ERROR(err);

    //origin + region of complete image to read from
    size_t origin[3] = { 0, 0, 0 }; //start from the bottom and read the whole thing
    size_t region[3] = {X_DEVICE, Y_DEVICE, 1}; //read the while thing

    //get image
    unsigned char* output_image;
    output_image = (unsigned char*)
        malloc(X_DEVICE * Y_DEVICE * 4 * sizeof(float));

    ccl_image_enqueue_read(img_out, queue, CL_TRUE, origin, region,
       0, 0, output_image, NULL, &err);
    HANDLE_ERROR(err);

    //read Z buffer back. 
    ccl_buffer_enqueue_read(zbufbuf, queue, CL_TRUE, 0,
        zSize, zbuf, NULL, &err);
    HANDLE_ERROR(err);

    //get Z buffer here
    /*
    to draw 2 shapes:
    usually just draws to graphics buffer seperatly (so we're all good)
    so if we ARE computing a point, we need to tell the painter to fuck right off 
    ...if we think we're drawing straight black!!!!!
    otherwise draw the point

    the other challenge is downloading the Z buffer and making it read write.(done)
    */

    //wait for shader to finish here
    ccl_queue_finish(queue, &err);
    HANDLE_ERROR(err);

    //copy openCL image to a nice comfy float array
    memcpy(output, output_image, sizeof( float ) * 800 * 800 * 4 );

    //and show it!
    for(int i = 0; i < X_DEVICE; i++){
        for(int j = 0; j < Y_DEVICE; j++){
           if(!(output[i][j][0] == 0.0 && output[i][j][1] == 0.0 && output[i][j][2] == 0.0)){
                G_rgb(output[i][j][0], output[i][j][1], output[i][j][2]);
                G_point(i, j);
           }
        }
    }

    //release image and buffers
    free(output_image);
    ccl_buffer_destroy(zbufbuf);

    ccl_image_destroy(img_out);

    //done
    return 0;
}

//all scene graph specific code goes here
int drawScene(double V[4][4]){
	//init transformation matrix
   	float color[3];
   	float m[4][4], minv[4][4] ;
	int tlist[100] ;
	float plist[100] ;
	int num;

    //clear buffer
    clearBuf();

    //solid of revolution MAY be revolving around the wrong thing...just my best guess

	//taurus
	num = 0;
	color[0] = 0.5; color[1] = 0.0; color[2] = 0.0;
	zPlot_acc(0, 2 * M_PI, 0, 2 * M_PI, V, color, TAURUS);

    //sphere
    color[0] = 0.5; color[1] = 0.5; color[2] = 0.0;
    zPlot_acc(0, 2 * M_PI, -M_PI/2, M_PI/2, V, color, SPHERE);

	return 1;
}

int main(){
    G_init_graphics(X_DEVICE, Y_DEVICE);
    setup_acc();

	double eye[3]; double coi[3]; double up[3];

    int numFrames = 100;
    int fnum = 0;
    double t = 0.0;
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

	    //q = G_wait_key(); if (q == 'q') { break; }
        G_display_image();

	    fnum++;
        G_rgb(0,0,0);
        G_clear();
  	}

    close_acc();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Code took %lf sec to run,\n", cpu_time_used);
    printf("Averaging %lf FPS\n", (double)numFrames/cpu_time_used);

  	return 0;
}

/*
acom `pkg-config --cflags cf4ocl2` acc_zbuf.c `pkg-config --libs cf4ocl2`
biggest bottleneck seems to be interop...
*/

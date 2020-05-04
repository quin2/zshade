#include <FPT.h>

//openCL helper lib for C
#include <cf4ocl2.h>

//handle error assertions
#include <assert.h>

#include <unistd.h>

//error handeling macros, borrowed from cf4ocl example project
#define ERROR_MSG_AND_EXIT(msg) \
    do { fprintf(stderr, "\n%s\n", msg); exit(EXIT_FAILURE); } while(0)
#define HANDLE_ERROR(err) \
    if (err != NULL) { ERROR_MSG_AND_EXIT(err->message); }

//image dimensions
#define X_DEVICE 800
#define Y_DEVICE 800

float output[X_DEVICE][Y_DEVICE][4];

int main(){
	//wipe array
	for(int i = 0; i < X_DEVICE; i++){
		for(int j = 0; j < Y_DEVICE; j++){
			for(int k = 0; k < 4; k++){
				output[i][j][k] = 0.5f;
			}
		}
	}
	//enviroment vars
	CCLContext* ctx;
    CCLDevice* dev;
    CCLImage* img_out;
    CCLQueue* queue;
    CCLProgram* prg;
    CCLKernel* krnl;
    CCLEventWaitList ewl = NULL;

    //i'll use the GPU :)
    int dev_idx = 1;

    //error handeling object
    CCLErr* err = NULL;

    //fancee image
    cl_image_format image_format = { CL_RGBA, CL_FLOAT };

    //context
    ctx = ccl_context_new_from_menu_full(&dev_idx, &err);
    HANDLE_ERROR(err);

    //get specified device in context 
    dev = ccl_context_get_device(ctx, 0, &err);
    HANDLE_ERROR(err);

    //create command queue
    queue = ccl_queue_new(ctx, dev, 0, &err);
    HANDLE_ERROR(err);

    //make an image output var for us :)
    img_out = ccl_image_new(ctx, CL_MEM_WRITE_ONLY,
        &image_format, NULL, &err,
        "image_type", (cl_mem_object_type) CL_MEM_OBJECT_IMAGE2D,
        "image_width", (size_t) X_DEVICE,
        "image_height", (size_t) Y_DEVICE,
        NULL);
    HANDLE_ERROR(err);

    //load shader from file
    prg = ccl_program_new_from_source_file(ctx, "shard.cl", &err);
    HANDLE_ERROR(err);

    //build shader
    const char * bldlog;
    CCLErr* err_bld = NULL;
    ccl_program_build(prg, NULL, &err);
    if ((err) && (err->code == CL_BUILD_PROGRAM_FAILURE)) {
		bldlog = ccl_program_get_build_log(prg, &err_bld);
		fprintf(stderr, "Error building program: \n%s", bldlog);
	}

	//get kernal wrapper
    krnl = ccl_program_get_kernel(prg, "art", &err);
    HANDLE_ERROR(err);

    //estimate local and global worksizes for device
    size_t gws[2]; //global work size, don't touch this
    size_t lws[2]; //local work size, don't touch this
    size_t real_ws[2]; //real work size

    real_ws[0] = X_DEVICE; //set these to work sizes that we defined
    real_ws[1] = Y_DEVICE;

    ccl_kernel_suggest_worksizes(krnl, dev, 2, real_ws, gws, lws, &err);
    HANDLE_ERROR(err);

    //variable transmission
    float toSend = 0.5f;
    cl_float toSend_k = (cl_float)toSend;

    //let's make up another to send (this one will be harder)
    cl_float3 color_k = (cl_float3){0.5f, 0.0f, 0.0f};
    //we might need to transmit this one with the help of a buffer...
    //make one here

    /*
    CCLBuffer* colorBuf;
	colorBuf = ccl_buffer_new(ctx, CL_MEM_READ_ONLY, sizeof(cl_float3), NULL, &err);
	HANDLE_ERROR(err);
	//annnnnnnd we gotta copy to it
	CCLEvent* evt_write;
    evt_write = ccl_buffer_enqueue_write(colorBuf, queue, CL_FALSE, 0,
        sizeof(cl_float3), &color_k, NULL, &err);
    HANDLE_ERROR(err);
    //and wait for us to finish
    ccl_event_wait_list_add(&ewl, evt_write, NULL);
    */


    //run it
    //could be ccl_arg_local(toSend_k, cl_float)
    ccl_kernel_set_args_and_enqueue_ndrange(
        krnl, queue, 2, NULL, gws, lws, NULL, &err,
        ccl_arg_priv(toSend_k, cl_float), ccl_arg_priv(color_k, cl_float3), img_out, NULL);
    HANDLE_ERROR(err);

    //and waitttttttt
    ccl_enqueue_barrier(queue, &ewl, &err);
    HANDLE_ERROR(err);

    //get zee image
    size_t origin[3] = { 0, 0, 0 };
    size_t region[3] = {X_DEVICE, Y_DEVICE, 1};

    unsigned char* output_image;
    //allocate space
    output_image = (unsigned char*)
        malloc(X_DEVICE * Y_DEVICE * 4 * sizeof(float));

    ccl_image_enqueue_read(img_out, queue, CL_TRUE, origin, region,
        0, 0, output_image, NULL, &err);
    HANDLE_ERROR(err);

    //should b finishing here
   	ccl_queue_finish(queue, &err);
   	HANDLE_ERROR(err);

   	//get data buffer to on device memory
    //get the array one row at a time
    /*
    int sz = X_DEVICE * 4; //where 4 is number of floats in a float4
    for(int j = 0; j < Y_DEVICE; j++){
    	memcpy(&output[j], &output_image[j * sz], sizeof(float) * sz);
	}
    */

    //this version only copies a small amount of the stuff..
    memcpy(output, output_image, sizeof( float ) * 800 * 800 * 4 );

    G_init_graphics(X_DEVICE, Y_DEVICE);

    //either way first ever CL/FPT interop!
    for(int i = 0; i < X_DEVICE; i++){
    	for(int j = 0; j < Y_DEVICE; j++){
    		//the hard way
    		/*
    		float rgba[4];
    		int index = (j * X_DEVICE) + i;
    		memcpy(&rgba[0], &output_image[index *4], sizeof(float) * 4);
    		G_rgb(rgba[0], rgba[1], rgba[2]);
    		*/
    		
    		G_rgb(output[i][j][0], output[i][j][1], output[i][j][2]);
    		G_point(i, j);
    	}
	}

	int q;
	q = G_wait_key();

    //free image
    free(output_image);
    //release wrappers
    ccl_image_destroy(img_out);
    ccl_program_destroy(prg);
    ccl_queue_destroy(queue);
    ccl_context_destroy(ctx);
    //check for release
    assert(ccl_wrapper_memcheck());
    //exit
    return EXIT_SUCCESS;
}

/*
build status:
shader output image to FPT translation layer working
passing in floats working
passing in float vectors working!
passing in arrays is almost canon at this point, so no need to fuck with it
*/

/*
    ccl_kernel_set_arg(krnl, 0, ccl_arg_local(du_s, cl_float));

    ccl_kernel_set_arg(krnl, 1, ccl_arg_local(dv_s, cl_float));
    ccl_kernel_set_arg(krnl, 2, ccl_arg_local(H_s, cl_float));
    ccl_kernel_set_arg(krnl, 3, ccl_arg_local(xd_s, cl_int));
    ccl_kernel_set_arg(krnl, 4, ccl_arg_local(yd_s, cl_int));

    ccl_kernel_set_arg(krnl, 5, ccl_arg_local(half_s, cl_float));

    ccl_kernel_set_arg(krnl, 6, colorBuf);
    ccl_kernel_set_arg(krnl, 7, matBuf);

    ccl_kernel_set_arg(krnl, 8, ccl_arg_local(ambient_s, cl_float));
    ccl_kernel_set_arg(krnl, 9, ccl_arg_local(diffuse_s, cl_float));
    ccl_kernel_set_arg(krnl, 10, ccl_arg_local(specpow_s, cl_int));

    ccl_kernel_set_arg(krnl, 11, zbufbuf);
    ccl_kernel_set_arg(krnl, 12, img_out);

    ccl_kernel_enqueue_ndrange(krnl, queue, 2, 0, gws, lws, NULL, &err);
    HANDLE_ERROR(err);
    */
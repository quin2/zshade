//openCL helper lib for C
#include <cf4ocl2.h>

//handle error assertions
#include <assert.h>

#define X_DEVICE 800
#define Y_DEVICE 800

float shOUT[X_DEVICE][Y_DEVICE]

//can CL shaders be used for art? let's see!
int main(){
	CCLContext* ctx;
    CCLDevice* dev;
    CCLImage* img_out;
    CCLQueue* queue;
    CCLProgram* prg;
    CCLKernel* krnl;

    //we'll have the GPU, please :)
    int dev_idx = 1;

    //error handeling object
    CCLErr* err = NULL;

    //image format
    cl_image_format image_format = { CL_RGBA, CL_FLOAT };

    
}
__kernel void art(float inv, float3 inVec, __write_only image2d_t output_img){
    int x = get_global_id(0);
    int y = get_global_id(1);

    float th = float(x)/800.0f;
    float4 col = {th, 0.0f, 0.5f, 0.0f};

    
    float thing = (float)((x * x) + (y * y));
    if(x < 20){
    	col.y = 0.5f;
    }



    //put it out there for the world to see
    write_imagef(output_img, (int2)(x, y), col);
}

//"pink seas" quinno
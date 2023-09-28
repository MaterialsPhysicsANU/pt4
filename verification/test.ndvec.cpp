#include <iostream>
#include <omp.h>

#include "quick_rng.h"
#include "ndvec.h"

int getl(void){
    auto a = slice<float>({4,4});
    for(size_t i = 0; i < a.data.size(); ++i){
        a.data[i] = i;
    }

    float f1 = a.getl({1.3,3.7});
    float f2 = a.getl({2,3});
    float f3 = a.getl({0.4,2.1});
    float f4 = a.getl({-1.2,0.6});

    std::cout << a << std::endl;
    std::cout << f1 << std::endl;
    std::cout << f2 << std::endl;
    std::cout << f3 << std::endl;
    std::cout << f4 << std::endl;

    const float tol = 1E-5;
    if(std::abs(f1 - 13.3f) > tol){return 1;}
    if(std::abs(f2 - 14.0f) > tol){return 1;}
    if(std::abs(f3 - 8.8f) > tol){return 1;}
    if(std::abs(f4 - 2.4f) > tol){return 1;}

    return 0;
}

int addl(void){
    auto a = slice<float>({4,4},0);
    std::vector<float> reference = {
                             0.21f,0.09f,0.0f,0.0f,
                             0.49f,0.21f,0.0f,0.0f,
                             0.0f,0.7f,1.3f,0.0f,
                             0.0f,0.0f,0.0f,0.0f,
                             };

    a.addl({1.3,2},1);
    a.addl({2,2},1);
    a.addl({0.3,0.7},1);

    std::cout << a << std::endl;

    bool match = true;
    for(size_t i = 0; i < a.data.size(); ++i){
        match &= (std::abs(a.data[i] - reference[i]) < 1.0E-7);
    }
    if(!match){
        return 1;
    }

    return 0;
}

int atomic(void){
    const size_t sz = 16;
    auto a = volume<float>({sz,sz,sz},0);

    xorshift xorshift0;
    xorshift0.srand(0xE8F2'AC24);
    xorshift0.rand(37u);
    #pragma omp parallel for
    for(size_t i = 0; i < 1'000'000; ++i){
        xorshift xorshift_private = xorshift0;
        float x = (xorshift_private.rand()/UINT32_MAX);
        float y = (xorshift_private.rand()/UINT32_MAX);
        float z = (xorshift_private.rand()/UINT32_MAX);
        a.addl({x,y,z},1);

    }

    if(a.sum() != 1'000'000){
        std::cout << a.sum() << std::endl;
        return 1;
    }

    return 0;
}

int gradient(void){
    volume<float> texture = ncf2vol("../verification/ndvec/texture.nc");
    volume<float> reference = ncf2vol("../verification/ndvec/reference.nc");

    volume<float3> grad = texture.grad();

    auto gradx = volume<float>(texture.size);

    for(size_t i = 0; i < gradx.data.size(); ++i){
        gradx.data[i] = grad.data[i][0];
    }

    bool match = true;
    for(size_t i = 0; i < gradx.data.size(); ++i){
        match &= (std::abs(gradx.data[i] - reference.data[i]) < 1.0E-9);
    }
    if(!match){
        vol2ncf("../verification/ndvec/gradx.nc", gradx);
        return 1;
    }
    return 0;
}

int convolve(void){
    volume<float> texture = ncf2vol("../verification/ndvec/texture.nc");
    volume<float> reference = ncf2vol("../verification/ndvec/ref_convolved.nc");


    std::array<float,3> kernel = {0.3333f, 0.3333f, 0.3333f};//uniform kernel
    volume<float> convolved = texture.convolve(kernel);

    bool match = true;
    for(size_t i = 0; i < convolved.data.size(); ++i){
        match &= (std::abs(convolved.data[i] - reference.data[i]) < 1.0E-9);
    }
    if(!match){
        vol2ncf("../verification/ndvec/convolved.nc", convolved);
        return 1;
    }
    return 0;
}

int slice2tiff(void){
    auto vecslice = slice<float3>({64,64});

    for(size_t ix = 0; ix < vecslice.size[0]; ++ix){
    for(size_t iy = 0; iy < vecslice.size[1]; ++iy){
        vecslice[{ix,iy}] = {sinf((float)ix/8.0f), cosf((float)iy/8.0f), 1.0f};
    }
    }

    slice2tiff("../verification/ndvec/rgba8.tiff", vecslice);
    return 0;
}

int main(int argc, char* argv[]){

#if _OPENMP
    omp_set_max_active_levels(8);
#endif

    if (!(argc > 1)) {
    std::cerr << "fatal error: no test selected" << std::endl;
    return -2;
    }
    std::string arg = argv[1];

    if(arg == "getl"){ return getl();}
    if(arg == "addl"){ return addl();}
    if(arg == "atomic"){ return atomic();}
    if(arg == "gradient"){ return gradient();}
    if(arg == "convolve"){ return convolve();}
    if(arg == "slice2tiff"){ return slice2tiff();}

    return -1;
}
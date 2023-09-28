#ifndef NDVEC_H
#define NDVEC_H

#include <array>
#include <queue>
#include <type_traits>

#include "vec.h"
#include "util.h"

template <typename T>
struct is_std_array : std::false_type {};

template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {};

template <typename T>
inline constexpr bool is_std_array_v = is_std_array<T>::value;

template<typename T, std::size_t N>
struct ndvec{
    std::vector<T> data;
    std::array<size_t, N> size;//x,y,z,...

    ndvec(void) = default;

    ndvec(std::array<size_t, N> sz){
        size = sz;
        size_t len = product(sz);
        data.resize(len);
    }

    ndvec(std::array<size_t, N> sz, T val){
        size = sz;
        size_t len = product(sz);
        data = std::vector<T>(len, val);
    }

    const inline T& operator[](const std::array<size_t, N>& index) const{
        return data[flatten(index, size)];
    }

    inline T& operator[](const std::array<size_t, N>& index){
        return data[flatten(index, size)];
    }

    template <typename U>
    std::enable_if_t<is_std_array_v<U>>
    atomic_add(const std::array<size_t, N>& index, const U& v){
            size_t i = flatten(index, size);
            const size_t M = v.size();
            for(size_t j = 0; j < M; ++j){
                #pragma omp atomic update
                data[i][j] += v[j];
            }
    }

    template <typename U>
    std::enable_if_t<!is_std_array_v<U>>
    atomic_add(const std::array<size_t, N>& index, const U v){
        size_t i = flatten(index, size);
        #pragma omp atomic update
        data[i] += v;
    }

    inline ndvec<T,N>& operator+=(const ndvec<T,N>& f){
        data += f.data;
        return *this;
    }

    inline ndvec<T,N>& operator+=(const T v){
        data += v;
        return *this;
    }

    inline ndvec<T,N>& operator-=(const ndvec<T,N>& f){
	    data -= f.data;
        return *this;
    }

    inline ndvec<T,N>& operator*=(const T v){
		data *= v;
        return *this;
    }
	
	template <typename U = T, typename = std::enable_if_t<is_std_array_v<U>>>
    inline ndvec<T, N>& operator*=(const ndvec<typename U::value_type, N>& v) {
        data *= v.data;
        return *this;
    }

    inline ndvec<T,N>& operator*=(const ndvec<T,N> f){
		data *= f.data;
        return *this;
    }

    inline ndvec<T,N>& operator/=(const ndvec<T,N>& f){
        data /= f.data;
        return *this;
    }

    inline void fill(T v){
        ::fill(data,v);
    }
    
    inline void clamp(T min, T max){
        clamp(data, min, max);
    }

    inline T max() const{
        return std::max_element(data.begin(), data.end())[0];
    }

    inline T min() const{
        return std::min_element(data.begin(), data.end())[0];
    }

    int save(const char* fp){
        std::ofstream wf(fp, std::ios::out |std::ios::binary);
        if(!wf){
            return 1;
        }
        wf.write((char *) &size,sizeof(std::array<size_t,N>));
        wf.write((char *) data.data(),data.size()*sizeof(T));
        
        wf.close();
        if(!wf.good()){
            return 1;
        }
    return 0;
    }

    static ndvec<T,N> load(const char* fp){
        auto fld = ndvec<T,N>();

        std::ifstream rf(fp, std::ios::in | std::ios::binary);
        if(!rf) {
        }

        rf.read((char *) &(fld.size), sizeof(std::array<size_t,N>));
        fld.data.resize(product(fld.size));
        rf.read((char *) fld.data.data(), fld.data.size()*sizeof(T));

        rf.close();
        if(!rf.good()) {
        }
        return fld;
    }

    friend std::ostream& operator<< (std::ostream& os,const ndvec<T,N>& f);

    //determines if index pair is in the slice
    inline bool in(const std::array<size_t, N>& index) const{
        bool isin = true;
        for(int dim = 0; dim < N; dim++){
        isin &= (0 <= index[dim] && index[dim] < size[dim]);
        }

        return isin;
    }

    //set value or does nothing
    inline void setz(const std::array<size_t,N>& index, T z){
        if(in(index)){
            (*this)[index] = z;
        }
    }

    //adds value n-linearly, ie with [ ,bi,tri]-linear interpolation
    inline void addl(const std::array<float,N>& index,const T z){
        std::array<float,  N> clamped_index;
        for(uint_fast8_t i = 0; i < N; ++i){
            clamped_index[i] = std::clamp(index[i], 0.0f, size[i]-1.0f);
        } 
        std::array<float,  N> frac;
        std::array<size_t, N> ncube = floor<size_t,float,N>(clamped_index,frac);
        
        //this will only work up till N = 64 //we keep the current index in the bits of vertex
        for(uint64_t vertex = 0; vertex < (1u << N); ++vertex){
            std::array<size_t,N> offset;//holds the offset for which vertex we are at
            for(uint_fast8_t i = 0; i < N; ++i){
                offset[i] = 0b1 & (vertex >> i);
            }

            std::array<size_t,N> vertex_index = ncube + offset;
            T vertex_v = z;//vertex value
            for(uint_fast8_t i = 0; i < N; ++i){//add weight
                vertex_v *=  offset[i] ? frac[i] : (1.0f - frac[i]);
            }
            if(in(vertex_index)){
                atomic_add(vertex_index, vertex_v);
            }
        }
    }

    //gets an element from the slice clamping to be in range
    inline T get(const std::array<size_t, N>& index) const{
        std::array<size_t,N> clamped;
        for(int dim = 0; dim < N; dim++){
            clamped[dim] = std::clamp(index[dim], size_t(0), size[dim]-size_t(1));
        }
        return (*this)[clamped];
    }

    //gets value or returns z
    inline T getz(const std::array<size_t, N>& index, T z) const{
        if(in(index)){
            return (*this)[index];
        }else{
            return z;
        }
    }

    //gets value n-linearly, ie with [ ,bi,tri]-linear interpolation
    inline T getl(const std::array<float,N>& index) const{
        std::array<float,  N> clamped_index;
        for(uint_fast8_t i = 0; i < N; ++i){
                clamped_index[i] = std::clamp(index[i], 0.0f, size[i]-1.0f);
        } 
        std::array<float,  N> frac;
        const std::array<size_t, N> ncube = floor<size_t,float,N>(clamped_index,frac);
        
        T v = T();//makes a zero scalar/vector etc
        
        //this will only work up till N = 64 //we keep the current index in the bits of vertex
        for(uint64_t vertex = 0; vertex < (1u << N); ++vertex){
            std::array<size_t,N> offset;//holds the offset for which vertex we are at
            for(uint_fast8_t i = 0; i < N; ++i){
                offset[i] = 0b1 & (vertex >> i);
            }

            std::array<size_t,N> vertex_index = ncube + offset;
            T vertex_v = get(vertex_index);//vertex value
            for(uint_fast8_t i = 0; i < N; ++i){//add weight
                vertex_v *=  offset[i] ? frac[i] : (1.0f - frac[i]);
            }
            v += vertex_v;
        }

        return v;
    }
    
    //fold with addition
    inline double sum(void) const {
        double acc = 0;
        for(T val : data){
            acc += val;
        }
        return acc;
    }
 
    inline void grad(ndvec<std::array<T,N>,N>& grd) const{
        #pragma omp parallel for
        for(size_t i = 0; i < data.size(); ++i){
            std::array<size_t,N> index = unflatten(i,size);
            for(size_t d = 0; d < N; d++){
                T dVi = T(0); 
                T dVidi = T(0);
                auto pos_step = std::array<size_t,N>();//initialises to zero
                auto neg_step = std::array<size_t,N>();//initialises to zero
                pos_step[d] += 1;
                neg_step[d] -= 1;

                dVi = get(pos_step + index) - get(neg_step + index);
                dVidi = dVi*size[d]/2;//normalise to [-1,1]^N space.
                grd.data[i][d] = dVidi;
            }
        }
    }

    inline ndvec<std::array<T,N>,N> grad(void) const{
        auto grd = ndvec<std::array<T,N>,N>(size);
        grad(grd);        
        return grd; 
    }

    template <typename U = T, std::size_t S>
    ndvec<T,N> convolve(const std::array<U,S>& kernel, std::array<bool,N> axes = {1,1,1}) const{
        const size_t kernal_size = kernel.size();
        const size_t data_size = data.size();
        auto convolved = ndvec<T,N>(*this);
        auto convolved_buffer = ndvec<T,N>(size);

        for(size_t d = 0; d < N; d++){
        if(axes[d]){
            #pragma omp parallel for if (data_size > (1UL << 20))
            for(size_t i = 0; i < data_size; ++i){
                std::array<size_t,N> index = unflatten(i,size);
                T sum = T();
                for(size_t j = 0; j < kernal_size; ++j){
                    std::array<size_t,N> sumindex = index;//adjust kernal to fit;
                    int sumid = static_cast<int>(sumindex[d]);
                    sumid = std::clamp(sumid + (int)j - (int)kernal_size/2, 0, (int)size[d]); 
                    sumindex[d] = static_cast<size_t>(sumid);

                    sum += kernel[j]*convolved.get(sumindex);
                }
                convolved_buffer.data[i] = sum;


            }
            convolved = convolved_buffer;
        }
        }
        
        return convolved; 
    }

    //supsamples a slice to a higher resolution
    ndvec<T,N> supsample(int level) const{
        std::array<size_t, N> sup_size = (this->size);
        if(level == 1){return ndvec<T,N>(*this);}
        sup_size *= (size_t)level;

        auto sup_vec = ndvec<T,N>(sup_size, T());

        const size_t ssz = sup_vec.data.size();
        for(size_t i = 0; i < ssz; i++){
            float3 src_index = array_static_cast<float>(unflatten(i, sup_size))/(float)level;
            sup_vec.data[i] = (*this).getl(src_index);
        }
        return sup_vec;
    }

    //subsamples a slice to a lower resolution
    ndvec<T,N> subsample(int level) const{
        std::array<size_t, N> sub_size = (this->size);
        if(level == 1){return ndvec<T,N>(*this);}
        sub_size /= (size_t)level;

        auto sub_vec = ndvec<T,N>(sub_size, T());

        const size_t ssz = sub_vec.data.size();
        for(size_t i = 0; i < ssz; i++){
            size_t3 src_index = unflatten(i, sub_size)*(size_t)level;
            sub_vec.data[i] += (*this)[src_index];
        }
        return sub_vec;
    }
};

template<typename T, std::size_t N>
inline ndvec<T,N>& operator*(ndvec<T,N>& f, const T v){
    f *= v;
    return f;
}

template<typename T, std::size_t N>
inline ndvec<T,N>& operator*( const T v, ndvec<T,N>& f){
    return f*v;
}

template<typename T, std::size_t N>
inline ndvec<T,N>& operator*(ndvec<T,N>& f, const ndvec<typename T::value_type,N> v){
    f *= v;
    return f;
}

template<typename T, std::size_t N>
inline ndvec<T,N>& operator*( const ndvec<typename T::value_type,N> v, ndvec<T,N>& f){
    return f*v;
}

template<typename T, std::size_t N>
inline ndvec<T,N>& operator+(ndvec<T,N>& f, const T v){
    f += v;
    return f;
}

template<typename T, std::size_t N>
inline ndvec<T,N>& operator+(const T v, ndvec<T,N>& f){
    return f+v;
}

template<typename T, std::size_t N>
inline  ndvec<T,N>& operator+(ndvec<T,N>& f, const ndvec<T,N>& h){
    f += h;
    return f;
}

template<typename T, std::size_t N>
inline  ndvec<T,N>& operator-(ndvec<T,N>& f, const ndvec<T,N>& h){
    f -= h;
    return f;
}

template<typename T, std::size_t N>
inline  ndvec<T,N>& operator/(ndvec<T,N>& f, const ndvec<T,N>& h){
    f /= h;
    return f;
}

template<typename T, std::size_t N>
    ndvec<typename T::value_type,N> dot(const ndvec<T,N>& f, const ndvec<T,N>& h){
    auto inner_product  = ndvec<typename T::value_type,N>(f.size);
    const size_t sz = f.data.size();
    for(size_t i = 0; i < sz; ++i){
        inner_product.data[i] = dot(f.data[i],h.data[i]);
    }
    return inner_product;
}

template<typename T, std::size_t N>
inline std::array<T,N> screen2voxel(const std::array<size_t,N>& size, const std::array<T,N>& screen_point){
    std::array<T, N> size_T;
    for(size_t i = 0; i < N; i++){
        size_T[i] = static_cast<T>(size[i]);
    }

    std::array<T,N> voxel_point = (screen_point+T(1))*(size_T/T(2));
    return voxel_point;
}

template<typename T, std::size_t N>
inline std::array<T,N> voxel2screen(const std::array<size_t,N>& size, const std::array<T,N>& voxel_point){
    std::array<T, N> size_T;
    for(size_t i = 0; i < N; i++){
        size_T[i] = static_cast<T>(size[i]);
    }
    return (voxel_point)*(T(2)/size_T)+(-T(1));
}

template <typename T>
struct volume : ndvec<T, 3>{
    using ndvec<T, 3>::ndvec;
    volume(ndvec<T, 3> v) : ndvec<T, 3>(std::move(v)) {}

    friend std::ostream& operator<< (std::ostream& os,const ndvec<T,3>& f){
        char delimiter = ',';

        for(size_t kpos = 0; kpos < f.size[2]; kpos++){
            for(size_t jpos = 0; jpos < f.size[1]; jpos++){
                for(size_t ipos = 0; ipos < f.size[0]; ipos++){
                    T val = f[{ipos,jpos,kpos}];
                    os << val << delimiter;

                }
                os.put('\n');
            }
            os.put('\n');
            os << "= = = = = = = = =\n";
            os.put('\n');
        }
        return os;
    }
};

int vol2ncf(const char* fp, const volume<float> vol);
volume<float> ncf2vol(const char* fp);

template <typename T>
struct slice : ndvec<T, 2> {
    using ndvec<T, 2>::ndvec;
    slice(ndvec<T, 2> v) : ndvec<T, 2>(std::move(v)) {}

    friend std::ostream& operator<< (std::ostream& os,const ndvec<T,2>& f){
        char delimiter = ',';

        for(size_t jpos = 0; jpos < f.size[1]; jpos++){
            for(size_t ipos = 0; ipos < f.size[0]; ipos++){
                T val = f[{ipos,jpos}];
                os << val << delimiter;

            }
            os.put('\n');
        }

        return os;
    }

    int y_stack(std::vector<T>& row){
        if (this->size[1] == 0){
            this->data = row;
            this->size[1] = 1;
            this->size[0] = this->data.size();
        }else{
            if (row.size() != this->size[0]) {return 0;}
            auto it = this->data.end();
            this->data.insert(it,row.begin(),row.end());
            this->size[1]++;
        }
        return 1;
    }

    //gets value bilinearly, clamping to be in range
    T get_bilinear(const std::array<size_t,2>& index) const{
        size_t x = index[0];
        size_t y = index[1];
        if(0 <= y && y < this->size[1]-1 && 0 <= x && x < this->size[0]-1){
            double yf,xf;//fractional parts
            double yi,xi;//integral parts
            size_t i,j;
            
            yf = modf(y, &yi);
            xf = modf(x, &xi);
            i = (int)yi;
            j = (int)xi;
            
            T v = 0.0;
            
            v += (*this)[{j,i}]     * (1.0 - xf) * (1.0 - yf);
            v += (*this)[{j,i+1}]   * (1.0 - xf) * (      yf);
            v += (*this)[{j+1,i}]   * (      xf) * (1.0 - yf);
            v += (*this)[{i+1,j+1}] * (      xf) * (      yf);
            
            return v;
        }else{
            size_t j = std::clamp<size_t>(x, 0.0, this->size[0]-1.0);
            size_t i = std::clamp<size_t>(y, 0.0, this->size[1]-1.0);
            return (*this)[{j,i}];
        }
    }

    //subsamples a slice to a lower resolution
    slice<T> subsample(int level) const{
        std::array<size_t, 2> sub_size = (this->size);
        if(level == 1){return slice<T>(*this);}
        sub_size /= (size_t)level;

        auto subf = slice<T>(sub_size, T(0));

        for(size_t i = 0; i < this->size[0]; i++){
        for(size_t j = 0; j < this->size[1]; j++){
            subf[{i/level,j/level}] += (*this)[{i,j}];
        }
        }
        return subf;
    }

    //interpolates a slice from a slice of a different size
    void interpolate(slice<T>& f){
        int w = f.size[0];
        int h = f.size[1];

        double xscale = (double)w/(double)this->size[0];
        double yscale = (double)h/(double)this->size[1];

        for(size_t j = 0; j< this->size[0]; j++){
        for(size_t i = 0; i< this->size[1]; i++){

            //exact position
            double jpos = j*xscale;
            double ipos = i*yscale;

            (*this)[{j,i}] = get_bilinear({static_cast<size_t>(jpos),static_cast<size_t>(ipos)});
        }
        }

    }

    //rotates the slice //radians
    void rot(double angle){
        slice<T> tmp = slice<T>(this->size);
        
        double c = cos(angle);
        double s = sin(angle);
        
        for(size_t j = 0; j< this->size[0]; j++){
        for(size_t i = 0; i< this->size[1]; i++){
            //rotate (j,i) on to coordinate on old image (x,y)
            size_t jc  = j-this->size[0]/2;
            size_t ic  = i-this->size[1]/2;
            double xc = c * jc - s * ic;  
            double yc = s * jc + c * ic;  
            
            tmp[{j,i}] = get_bilinear({
                    static_cast<size_t>(xc+this->size[0]/2.0),
                    static_cast<size_t>(yc+this->size[1]/2.0)
                    });


        }
        }            
        this->data = tmp.data;
    }

};

template<typename T>
int pol2cart(const slice<T>& p, slice<T>& c){
    double r_sf = p.size[0]/(0.5*c.size[0]);//radial scale factor
    for(int i = 0; i< c.size[1]; i++){
    for(int j = 0; j< c.size[0]; j++){
        //from rect (i,j) to polar (r,theta)
        int jc  = (j-c.size[0]/2) ;
        int ic  = (i-c.size[1]/2) ;
        
        size_t r = r_sf*sqrt(jc*jc + ic*ic);
        if(r >= p.size[0]){continue;}

        size_t theta = int(p.size[1]*(atan2(ic,jc)/(2*pi) + 1.0)) % p.size[1];
        
        size_t jw = (j+c.size[0]/2) % c.size[0];
        size_t iw = (i+c.size[1]/2) % c.size[1];
       
        c[{jw,iw}] = p.get_bilinear({r,theta});
        
    }
    } 
    return 1;
}

template<typename T> slice<T> tiff2slice(const char* fp);
template<typename T> int  slice2tiff(const char* fp, slice<T>& slice0);
template<typename T> int  slice2tiff(const char* fp, slice<std::array<T,3>>& slice0);

#endif

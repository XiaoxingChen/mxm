#if !defined(_CV_PIXEL_H_)
#define _CV_PIXEL_H_

#include <array>
#include <cmath>
#include <string>
#include "common.h"
#include "accessor.h"
#include "linalg_mat.h"


namespace mxm
{

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, uint8_t>::type
quantizeToU8(DType val)
{
    return val < 0 ? 0 : val > 1 ? 255 : val * 255.99;
}

template<typename DType>
typename std::enable_if<std::is_integral<DType>::value, uint8_t>::type
quantizeToU8(DType val) { return static_cast<uint8_t>(val); }

template<typename DeriveType, typename DType>
// class RGBA_F32
class RGBA
{
public:
    // using DType = float;
    // using DType = typename DeriveType::dtype;
    DeriveType* deriveThis() {return reinterpret_cast<DeriveType*>(this);}
    const DeriveType* deriveThis() const {return reinterpret_cast<const DeriveType*>(this);}
    const DType& at(size_t i) const {return (*deriveThis())(i);}

    const DType & r() const {return (*deriveThis())(0);}
    const DType & g() const {return (*deriveThis())(1);}
    const DType & b() const {return (*deriveThis())(2);}
    const DType & a() const {return deriveThis()->size() > 2 ? (*deriveThis())(3) : 0.;}

    int rU8() const { return quantizeToU8<DType>(r()); }
    int gU8() const { return quantizeToU8<DType>(g()); }
    int bU8() const { return quantizeToU8<DType>(b()); }
    int aU8() const { return quantizeToU8<DType>(a()); }
};

template<typename DType, size_t NChannel>
class PixelType: public LinearData<PixelType<DType, NChannel>, DType>, public RGBA<PixelType<DType,NChannel>, DType>
{
public:
    // static const size_t NChannel = 3;
    using InitialType = std::array<DType, NChannel>;
    using ThisType = PixelType<DType,NChannel>;
    // using DType = float;
    using BaseType0 = LinearData<PixelType<DType,NChannel>, DType>;
    PixelType(const InitialType& v): data_(v){}
    PixelType(){ }
    PixelType(DType val){ BaseType0::traverse([&](size_t i){(*this)(i) = val;}); }

    const DType& operator () (size_t i) const { return data_.at(i); }
    DType& operator () (size_t i) { return data_.at(i); }

    static ThisType black() {return ThisType(0);}
    static ThisType white() {return ThisType(1);}
    static constexpr size_t size() {return NChannel;}

    enum{
        N_CHANNEL = NChannel
    };

    // constexpr DType dtype() const;
    // using dtype = DType;

private:
    std::array<DType, NChannel> data_;
};

template<typename DType, size_t NChannel>
inline typename std::enable_if<std::is_same<DType, uint8_t>::value, std::vector<uint8_t>>::type
serialize(const Matrix<PixelType<DType, NChannel>>& img)
{
    std::vector<uint8_t> src_mem(NChannel * img.shape(0) * img.shape(1));
    img.traverse([&](auto i, auto j){
        for(size_t k = 0; k < NChannel; k++)
        {
            src_mem.at(i * img.shape(1) * NChannel + j * NChannel + k) = img(i,j)(k);
        }
    });
    return src_mem;
}

template<typename DType, size_t NChannel>
inline typename std::enable_if<std::is_floating_point<DType>::value, std::vector<uint8_t>>::type
serialize(const Matrix<PixelType<DType, NChannel>>& img)
{
    std::vector<uint8_t> src_mem(NChannel * img.shape(0) * img.shape(1));
    img.traverse([&](auto i, auto j){
        for(size_t k = 0; k < NChannel; k++)
        {
            src_mem.at(i * img.shape(1) * NChannel + j * NChannel + k) = quantizeToU8(img(i,j)(k));
        }
    });
    return src_mem;
}

using Pixel = PixelType<float, 3>;

template<typename DType, size_t N>
std::string to_string(const PixelType<DType, N>& px, size_t prec=6)
{
    std::string ret("(");
    px.traverse([&](size_t i){
        ret += (px(i) >= 0 && i > 0 ? "," : "");
        ret += mxm::to_string(px(i), prec); });

    ret += ")";
    return ret;
}

namespace random
{

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type
uniform();

template<typename DType>
typename std::enable_if<
    std::is_same<
        DType, PixelType< typename DType::EntryType , DType::size()>
    >::value, DType
>::type
uniform()
{
    DType ret(0);
    for(size_t i = 0; i < DType::size(); i++) ret(i) = uniform<typename DType::EntryType>();
    return ret;
}

} // namespace random


} // namespace mxm


#endif // _CV_PIXEL_H_

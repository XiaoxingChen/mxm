#if !defined(_CV_PIXEL_H_)
#define _CV_PIXEL_H_

#include <array>
#include <cmath>
#include <string>
#include <limits>
#include <type_traits>
#include "common.h"
#include "accessor.h"
#include "linalg_mat.h"


namespace mxm
{

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value , std::string>::type
to_string(const DType& v, size_t prec);

template<typename DType>
typename std::enable_if<std::is_integral<DType>::value , std::string>::type
to_string(const DType& v, size_t prec);

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

    template<typename DTypeRHS>
    PixelType(const PixelType<DTypeRHS, NChannel>& rhs){ BaseType0::traverse([&](size_t i){(*this)(i) = DType(rhs(i));}); }

    const DType& operator () (size_t i) const { return data_.at(i); }
    DType& operator () (size_t i) { return data_.at(i); }

    static ThisType black() {return ThisType(0);}
    static ThisType white() {return ThisType(1);}
    static ThisType red() {auto ret = black(); ret(0) = 1; return ret;}
    static ThisType green() {auto ret = black(); ret(1) = 1; return ret;}
    static ThisType blue() {auto ret = black(); ret(2) = 1; return ret;}
    static constexpr size_t size() {return NChannel;}

#if 1
    template<typename Tin=DType, typename Tout=float>
    operator typename std::enable_if_t<std::is_integral<Tin>::value && NChannel == 1, Tout>() const
    {
        return Tout(data_[0]);
    }

    template<typename Tin=DType, typename Tout=uint8_t>
    operator typename std::enable_if_t<std::is_floating_point<Tin>::value && NChannel == 1, Tout>() const
    {
        return Tout(data_[0]);
    }

    template<typename Tin=DType, typename Tout=float>
    operator typename std::enable_if_t<std::is_same<Tin, Tout>::value && NChannel == 1, Tout>() const
    {
        return data_[0];
    }
#endif

    enum{
        N_CHANNEL = NChannel
    };

    // constexpr DType dtype() const;
    // using dtype = DType;

private:
    std::array<DType, NChannel> data_;
};

template<typename T>
constexpr std::enable_if_t<std::is_arithmetic<T>::value, size_t> channelNum() { return 1; }

template<typename PType>
constexpr std::enable_if_t<std::is_same<PType, PixelType< typename PType::EntryType, PType::size()>>::value, size_t> channelNum() { return PType::size(); }

template<typename DType, size_t N>
std::string to_string(const PixelType<DType, N>& px, size_t prec)
{
    std::string ret("(");
    px.traverse([&](size_t i){
        ret += (px(i) >= 0 && i > 0 ? "," : "");
        ret += mxm::to_string(px(i), prec); });

    ret += ")";
    return ret;
}

template<typename DType, size_t N>
std::string to_string(const PixelType<DType, N>& px)
{
    return to_string(px, 6);
}

template<typename DType, size_t N>
PixelType<DType, N> inv(const PixelType<DType, N>& val)
{
    PixelType<DType, N> ret;
    ret.traverse([&](size_t i){ret(i) = (DType(1)/val(i));});
    return ret;
}

template<typename PType>
Matrix<PixelType<uint8_t, channelNum<PType>()>> quantize(const Matrix<PType>& img)
{
    Matrix<PixelType<uint8_t, channelNum<PType>()>> tmp(img.shape(), {}, ROW);
    tmp = img;
    return tmp;
}

template<typename DType>
std::array<DType, 2> elementwiseBounds(const Matrix<DType>& mat);
inline void normalize(Matrix<float>& img, float range=255.)
{
    auto min_max = mxm::elementwiseBounds(img);
    img -= min_max[0];
    img *= (range/ (min_max[1] - min_max[0]));
}

template<typename PType>
typename std::enable_if_t<
    std::is_same<
        PType, typename mxm::PixelType< typename PType::EntryType, PType::size()>
    >::value,
void>
normalize(Matrix<PType>& img, float range=255.)
{
    using DType = typename PType::EntryType;
    std::array<DType, 2> min_max{
        std::numeric_limits<DType>::max(),
        std::numeric_limits<DType>::min()};

    for(size_t ch = 0; ch < PType::size(); ch++)
    {
        img.traverse([&](auto i, auto j){
            min_max[0] = std::min(min_max[0], img(i,j)(ch));
            min_max[1] = std::max(min_max[1], img(i,j)(ch));
        });
    }
    // auto min_max = elementwiseBounds(img);
    // std::cout << mxm::to_string(min_max[0]) << "," << mxm::to_string(min_max[1]) << std::endl;
    img -= min_max[0];
    img *= (range /(min_max[1] - min_max[0]));
    // std::cout << mxm::to_string(mxm::inv(min_max[1] - min_max[0]))  << std::endl;
}

template<typename PType>
Matrix<PType> normalized(const Matrix<PType>& img, float range=255.)
{
    Matrix<PType> ret(img);
    normalize(ret, range);
    return ret;
}

using Pixel = PixelType<float, 3>;

template<typename LType, typename DType>
typename std::enable_if_t<
    std::is_same<
        typename PixelType< typename LType::EntryType, LType::size()>::value, LType
    >::value , Matrix<DType>
> operator * (LType lhs, const Matrix<DType>& rhs) { return rhs * lhs;}

// template<typename DType, size_t N>
// struct NormTraits<PixelType<DType, N>>{using type = typename PixelType<DType, N>::EntryType;};
template<typename DType, size_t N>
struct Traits<PixelType<DType, N>>{
    using ArithType = DType;
    using EntryType = DType;
    static PixelType<DType, N> identity() { return PixelType<DType, N>(1); }
    static PixelType<DType, N> zero() { return PixelType<DType, N>(0); }
};

namespace random
{

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type
uniform();

template<typename DType>
typename std::enable_if<
    std::is_same<
        DType, ::mxm::PixelType< typename DType::EntryType , DType::size()>
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

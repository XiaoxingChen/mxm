#if !defined(_LINALG_MAT_BLOCK_H_)
#define _LINALG_MAT_BLOCK_H_

#include <initializer_list>
#include "linalg_mat.h"

namespace mxm
{

class Indexer
{
public:
    Indexer():index_(0), from_right_(false), as_end_(false){}

    bool empty() const
    {
        return (0 == index_) && (false == from_right_);
    }

    Indexer(size_t idx)
    :index_(idx), from_right_(false){}

    Indexer(size_t idx, bool from_right, bool as_end)
    :index_(idx),from_right_(from_right), as_end_(as_end){}

    void setEnd(bool v) { as_end_ = v; }

    size_t eval(size_t size) const
    {
        if(empty()) return as_end_ ? size : 0;

        if(from_right_) return size - index_;
        return index_;
    }

    Indexer operator - (size_t i) const
    {
        if(from_right_)
            return Indexer(index_ + i, true, as_end_);
        return Indexer(index_ - i, false, as_end_);
    }

    Indexer operator + (size_t i) const
    {
        if(from_right_)
            return Indexer(index_ - i, true, as_end_);
        return Indexer(index_ + i, false, as_end_);
    }

private:
    size_t index_;
    bool from_right_;
    bool as_end_;

};

inline Indexer end()
{
    return Indexer(0, true, false);
}

struct Range
{

    Range(std::initializer_list<Indexer> indices)
    {
        if(indices.size() > 0) begin = *indices.begin();
        if(indices.size() > 1) end = *(indices.begin() + 1);
        begin.setEnd(false);
        end.setEnd(true);
    }

    Range()
    {
        begin.setEnd(false);
        end.setEnd(true);
    }
    size_t size(size_t full_size) const { return end.eval(full_size) - begin.eval(full_size); }
    Indexer begin;
    Indexer end;
};

class Block
{
public:
    Block(const Range& row, const Range& col)
        :row_(row), col_(col) {}

    template<typename DType>
    std::array<size_t, 4> getBlock(const Matrix<DType>& mat) const
    {
        std::array<size_t, 4> ret;
        ret[0] = row_.begin.eval(mat.shape(0));
        ret[1] = row_.end.eval(mat.shape(0));
        ret[2] = col_.begin.eval(mat.shape(1));
        ret[3] = col_.end.eval(mat.shape(1));
        return ret;
    }

    std::array<std::array<size_t, 2>, 2> deduct(const Shape& shape) const
    {
        std::array<std::array<size_t, 2>, 2> offset_shape;
        offset_shape[0] = {row_.begin.eval(shape[0]), col_.begin.eval(shape[1])};
        offset_shape[1] = {
            row_.end.eval(shape[0]) - row_.begin.eval(shape[0]),
            col_.end.eval(shape[1]) - col_.begin.eval(shape[1])};
        return offset_shape;
    }
protected:
    Range row_;
    Range col_;
};

inline std::array<std::array<size_t, 2>, 2> deduct(const Block& b, const Shape& mat)
{
    return b.deduct(mat);
}

// inline Mat& Mat::set(const Block& s, const Mat& rhs)
// {
//     auto row01_col01 = s.getBlock(*this);
//     return Mat::setBlock(row01_col01[0], row01_col01[2], rhs);
// }

inline Block Row(const Indexer& idx)
{
    return Block({idx, idx + 1}, {0,end()});
}


inline Block Col(const Indexer& idx)
{
    return Block({0,end()}, {idx, idx + 1});
}

} // namespace mxm
#endif // _LINALG_MAT_BLOCK_H_

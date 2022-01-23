#if !defined(_TEST_METRIC_TREE_H_)
#define _TEST_METRIC_TREE_H_

#include "mxm/spatial_metric_tree.h"
#include "mxm/linalg.h"

using namespace mxm;

class EuclideanPointMetricTree
:public MetricTreeBase<Matrix<float>, Matrix<float>, float >
{
public:
    virtual float distance(const Matrix<float>& p0, const Matrix<float>& p1) const override
    {
        return mxm::norm(p0 - p1);
    }

    virtual Matrix<float> point(size_t i) const override
    {
        return point_set_(Col(i));
    }

    virtual size_t size() const { return point_set_.shape(1); };

    EuclideanPointMetricTree(const Matrix<float>& points): point_set_(points) {}

    std::multimap<float, size_t> naiveRadiusSearch(const Matrix<float>& query_point, float radius) const
    {
        std::multimap<float, size_t> ret;
        for(size_t i = 0; i < size(); i++)
        {
            auto dist = distance(point(i), query_point);
            if(dist < radius)
                ret.insert({dist, i});
        }
        return ret;
    }

    std::multimap<float, size_t> naiveNearestNeighbourSearch(const Matrix<float>& query_point, size_t k) const
    {
        std::multimap<float, size_t> ret;
        for(size_t i = 0; i < size(); i++)
        {
            auto dist = distance(point(i), query_point);
            ret.insert({dist, i});
        }

        while(ret.size() > k)
        {
            ret.erase(std::prev(ret.end()));
        }
        return ret;
    }

private:
    Matrix<float> point_set_;
};

void testMetricTree()
{
    Matrix<float> points(fixRow(2),
    {
        0.749941, 0.910648, 0.132996, 0.181847, 0.982361, 0.263803, 0.095355, 0.145539, 0.282673, 0.136069, 0.802111, 0.869292, 0.077557, 0.579705, 0.627384, 0.549860, 0.008094, 0.144955, 0.680287, 0.853031, 0.533933, 0.622055, 0.438667, 0.350952, 0.199551, 0.513250, 0.138001, 0.401808, 0.382333, 0.075967, 0.762421, 0.239916, 0.040471, 0.123319, 0.252956, 0.183908, 0.504771, 0.239953, 0.822605, 0.417267, 0.981723, 0.049654, 0.823455, 0.902716, 0.301827, 0.944787, 0.047944, 0.490864, 0.247848, 0.489253, 0.544056, 0.337719, 0.887726, 0.900054, 0.343840, 0.369247, 0.033269, 0.111203, 0.162872, 0.780252, 0.877364, 0.389739, 0.210302, 0.241691, 0.272753, 0.403912, 0.492442, 0.096455, 0.223420, 0.131973, 0.490301, 0.942051, 0.954943, 0.956135, 0.651969, 0.575209, 0.757504, 0.059780, 0.435246, 0.234780, 0.552881, 0.353159, 0.053153, 0.821194, 0.322511, 0.015403, 0.404981, 0.043024, 0.908434, 0.168990, 0.809137, 0.649116, 0.258296, 0.731722, 0.129847, 0.647746, 0.493327, 0.450924, 0.378500, 0.547009,
        0.718470, 0.296321, 0.287805, 0.744693, 0.623436, 0.188955, 0.811874, 0.686775, 0.312508, 0.183511, 0.385711, 0.368485, 0.343930, 0.625619, 0.815769, 0.780227, 0.669285, 0.081126, 0.379419, 0.929386, 0.458497, 0.775713, 0.303614, 0.486792, 0.910565, 0.435859, 0.488618, 0.446784, 0.755790, 0.306349, 0.108062, 0.508509, 0.393590, 0.510772, 0.870187, 0.817628, 0.360861, 0.794831, 0.917118, 0.644318, 0.543805, 0.378609, 0.140144, 0.811580, 0.199873, 0.532826, 0.948925, 0.350727, 0.990110, 0.939002, 0.240076, 0.875943, 0.016521, 0.550156, 0.388615, 0.622475, 0.779689, 0.587045, 0.476638, 0.207742, 0.540138, 0.301246, 0.018451, 0.470923, 0.900183, 0.230488, 0.194495, 0.844309, 0.885998, 0.194764, 0.441223, 0.225922, 0.147829, 0.170708, 0.239502, 0.227664, 0.799653, 0.435699, 0.473015, 0.311102, 0.089823, 0.923380, 0.644551, 0.430207, 0.633064, 0.184816, 0.584382, 0.904881, 0.726654, 0.979748, 0.354638, 0.438870, 0.680407, 0.111119, 0.707322, 0.258065, 0.162329, 0.408720, 0.133736, 0.594896
    }, ROW);
    EuclideanPointMetricTree tree(points);
    tree.build(true);
    bool valid_tree = tree.verifyTree(0);
    if(!valid_tree)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));


    {
        auto expected = tree.naiveRadiusSearch(Vector<float>{0.3, 0.3}, 0.3);
        auto result = tree.radiusSearch(Vector<float>{0.3, 0.3}, 0.3);

        if(expected != result)
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        Vec target_pt({.5, .5});
        // size_t k = 5;
        for(size_t k = 5; k < 20; k+= 3)
        {
            auto result = tree.nearestNeighborSearch(target_pt, k);
            std::multimap<float, size_t> expected = tree.naiveNearestNeighbourSearch(target_pt, k);

            if(expected != result)
            {
                throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
            }
        }


    }


}



#endif // _TEST_METRIC_TREE_H_

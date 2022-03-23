#include "test_config.h"

#if TEST_AVAILABLE_CV_BASIC
#include <iostream>
#include "mxm/cv_basic.h"
#include "mxm/linalg.h"
#include "mxm/cv_kernel.h"
#include "mxm/random.h"
#include "mxm/cv_corner.h"
#include "mxm/cv_3d.h"
#include "mxm/model_camera.h"
#include "mxm/lie_special_orthogonal.h"
#include "mxm/geometry_cube.h"
#include "mxm/cv_calibration.h"


using namespace mxm;

// Data source:
// EuroC V1_01_easy cam0 1403715287062142976.png (200,303)
const Matrix<float>& dataImagePatch()
{
    static Matrix<float> patch({41,41}, {
    0.533333, 0.545098, 0.549020, 0.552941, 0.549020, 0.560784, 0.560784, 0.560784, 0.545098, 0.552941, 0.588235, 0.615686, 0.639216, 0.654902, 0.654902, 0.658824, 0.650980, 0.666667, 0.666667, 0.737255, 0.905882, 0.921569, 0.831373, 0.815686, 0.792157, 0.847059, 0.972549, 0.960784, 0.796078, 0.631373, 0.619608, 0.658824, 0.686275, 0.705882, 0.705882, 0.690196, 0.690196, 0.674510, 0.662745, 0.686275, 0.701961,
    0.447059, 0.466667, 0.494118, 0.509804, 0.521569, 0.533333, 0.552941, 0.556863, 0.552941, 0.560784, 0.529412, 0.529412, 0.564706, 0.580392, 0.588235, 0.607843, 0.603922, 0.607843, 0.607843, 0.666667, 0.776471, 0.827451, 0.827451, 0.819608, 0.807843, 0.831373, 0.905882, 1.000000, 0.917647, 0.698039, 0.631373, 0.631373, 0.690196, 0.686275, 0.690196, 0.670588, 0.690196, 0.666667, 0.682353, 0.674510, 0.709804,
    0.262745, 0.270588, 0.294118, 0.325490, 0.349020, 0.388235, 0.423529, 0.462745, 0.482353, 0.490196, 0.486275, 0.486275, 0.505882, 0.560784, 0.600000, 0.615686, 0.647059, 0.647059, 0.647059, 0.694118, 0.756863, 0.874510, 0.878431, 0.839216, 0.803922, 0.807843, 0.843137, 0.976471, 1.000000, 0.796078, 0.647059, 0.619608, 0.674510, 0.666667, 0.662745, 0.682353, 0.690196, 0.678431, 0.670588, 0.666667, 0.678431,
    0.254902, 0.227451, 0.215686, 0.211765, 0.211765, 0.231373, 0.294118, 0.352941, 0.364706, 0.333333, 0.305882, 0.282353, 0.333333, 0.478431, 0.603922, 0.678431, 0.678431, 0.682353, 0.670588, 0.690196, 0.760784, 0.886275, 0.996078, 0.850980, 0.796078, 0.796078, 0.811765, 0.886275, 1.000000, 0.933333, 0.733333, 0.631373, 0.635294, 0.647059, 0.674510, 0.670588, 0.678431, 0.686275, 0.674510, 0.654902, 0.666667,
    0.400000, 0.313726, 0.227451, 0.211765, 0.207843, 0.262745, 0.345098, 0.325490, 0.337255, 0.309804, 0.235294, 0.215686, 0.239216, 0.419608, 0.741176, 0.850980, 0.878431, 0.854902, 0.823529, 0.733333, 0.764706, 0.800000, 0.870588, 0.823529, 0.796078, 0.796078, 0.811765, 0.835294, 0.972549, 1.000000, 0.850980, 0.670588, 0.615686, 0.643137, 0.686275, 0.682353, 0.682353, 0.686275, 0.678431, 0.670588, 0.682353,
    0.580392, 0.486275, 0.298039, 0.215686, 0.215686, 0.239216, 0.254902, 0.262745, 0.266667, 0.247059, 0.207843, 0.192157, 0.215686, 0.392157, 0.768627, 0.917647, 0.964706, 0.980392, 0.960784, 0.827451, 0.776471, 0.764706, 0.800000, 0.803922, 0.811765, 0.807843, 0.803922, 0.796078, 0.894118, 1.000000, 0.972549, 0.745098, 0.627451, 0.635294, 0.658824, 0.674510, 0.678431, 0.694118, 0.686275, 0.674510, 0.678431,
    0.662745, 0.564706, 0.329412, 0.247059, 0.219608, 0.207843, 0.207843, 0.211765, 0.200000, 0.196078, 0.192157, 0.192157, 0.200000, 0.313726, 0.627451, 0.807843, 0.870588, 0.894118, 0.925490, 0.847059, 0.768627, 0.756863, 0.780392, 0.827451, 0.882353, 0.815686, 0.796078, 0.815686, 0.815686, 0.929412, 1.000000, 0.898039, 0.678431, 0.607843, 0.635294, 0.666667, 0.650980, 0.674510, 0.674510, 0.666667, 0.662745,
    0.686275, 0.552941, 0.392157, 0.305882, 0.223529, 0.203922, 0.196078, 0.188235, 0.184314, 0.192157, 0.192157, 0.203922, 0.223529, 0.274510, 0.372549, 0.427451, 0.447059, 0.494118, 0.639216, 0.780392, 0.752941, 0.733333, 0.764706, 0.870588, 0.972549, 0.874510, 0.815686, 0.796078, 0.792157, 0.843137, 0.972549, 0.980392, 0.737255, 0.615686, 0.607843, 0.654902, 0.647059, 0.666667, 0.690196, 0.690196, 0.670588,
    0.647059, 0.866667, 0.898039, 0.545098, 0.266667, 0.215686, 0.192157, 0.184314, 0.188235, 0.192157, 0.200000, 0.231373, 0.309804, 0.478431, 0.701961, 0.584314, 0.325490, 0.282353, 0.325490, 0.478431, 0.643137, 0.701961, 0.725490, 0.772549, 0.870588, 0.827451, 0.796078, 0.811765, 0.807843, 0.823529, 0.901961, 0.992157, 0.854902, 0.647059, 0.600000, 0.635294, 0.647059, 0.658824, 0.678431, 0.670588, 0.674510,
    0.478431, 0.976471, 1.000000, 0.819608, 0.341176, 0.227451, 0.200000, 0.184314, 0.192157, 0.196078, 0.219608, 0.364706, 0.588235, 0.662745, 0.611765, 0.513726, 0.345098, 0.254902, 0.247059, 0.278431, 0.388235, 0.529412, 0.635294, 0.733333, 0.756863, 0.772549, 0.796078, 0.815686, 0.803922, 0.819608, 0.835294, 0.933333, 0.964706, 0.745098, 0.615686, 0.592157, 0.635294, 0.662745, 0.682353, 0.690196, 0.694118,
    0.368627, 0.721569, 0.984314, 0.803922, 0.341176, 0.215686, 0.200000, 0.196078, 0.188235, 0.192157, 0.227451, 0.333333, 0.364706, 0.376471, 0.349020, 0.305882, 0.262745, 0.239216, 0.235294, 0.243137, 0.270588, 0.317647, 0.458824, 0.611765, 0.701961, 0.760784, 0.862745, 0.847059, 0.800000, 0.807843, 0.803922, 0.882353, 0.976471, 0.882353, 0.670588, 0.603922, 0.615686, 0.650980, 0.662745, 0.662745, 0.682353,
    0.396078, 0.392157, 0.419608, 0.349020, 0.243137, 0.203922, 0.200000, 0.192157, 0.192157, 0.192157, 0.203922, 0.219608, 0.254902, 0.278431, 0.250980, 0.235294, 0.227451, 0.231373, 0.250980, 0.278431, 0.278431, 0.286275, 0.294118, 0.380392, 0.564706, 0.733333, 0.905882, 0.909804, 0.803922, 0.784314, 0.780392, 0.815686, 0.933333, 0.960784, 0.772549, 0.603922, 0.588235, 0.615686, 0.654902, 0.654902, 0.658824,
    0.682353, 0.435294, 0.321569, 0.262745, 0.231373, 0.219608, 0.203922, 0.203922, 0.196078, 0.200000, 0.207843, 0.227451, 0.239216, 0.235294, 0.231373, 0.215686, 0.215686, 0.243137, 0.384314, 0.564706, 0.490196, 0.380392, 0.290196, 0.282353, 0.372549, 0.576471, 0.733333, 0.819608, 0.788235, 0.776471, 0.780392, 0.788235, 0.839216, 0.949020, 0.898039, 0.690196, 0.600000, 0.596078, 0.635294, 0.647059, 0.639216,
    0.905882, 0.749020, 0.541176, 0.337255, 0.262745, 0.215686, 0.203922, 0.200000, 0.196078, 0.196078, 0.227451, 0.317647, 0.337255, 0.313726, 0.274510, 0.250980, 0.239216, 0.286275, 0.658824, 0.956863, 0.968627, 0.807843, 0.380392, 0.309804, 0.482353, 0.603922, 0.572549, 0.611765, 0.709804, 0.756863, 0.760784, 0.764706, 0.788235, 0.870588, 0.968627, 0.815686, 0.627451, 0.592157, 0.603922, 0.623529, 0.631373,
    0.964706, 0.929412, 0.835294, 0.615686, 0.376471, 0.274510, 0.231373, 0.207843, 0.203922, 0.196078, 0.274510, 0.588235, 0.749020, 0.674510, 0.541176, 0.392157, 0.266667, 0.278431, 0.450980, 0.615686, 0.760784, 0.737255, 0.372549, 0.298039, 0.403922, 0.400000, 0.329412, 0.380392, 0.596078, 0.694118, 0.721569, 0.756863, 0.780392, 0.811765, 0.937255, 0.921569, 0.717647, 0.592157, 0.611765, 0.631373, 0.631373,
    0.972549, 0.933333, 0.941176, 0.854902, 0.690196, 0.443137, 0.298039, 0.239216, 0.215686, 0.207843, 0.231373, 0.435294, 0.831373, 0.909804, 0.905882, 0.709804, 0.317647, 0.239216, 0.250980, 0.278431, 0.298039, 0.301961, 0.247059, 0.219608, 0.239216, 0.243137, 0.250980, 0.266667, 0.333333, 0.458824, 0.600000, 0.701961, 0.737255, 0.776471, 0.854902, 0.956863, 0.823529, 0.643137, 0.596078, 0.607843, 0.639216,
    0.949020, 0.972549, 0.945098, 0.901961, 0.874510, 0.752941, 0.537255, 0.341176, 0.235294, 0.219608, 0.215686, 0.270588, 0.384314, 0.466667, 0.509804, 0.407843, 0.250980, 0.223529, 0.231373, 0.250980, 0.262745, 0.231373, 0.215686, 0.207843, 0.219608, 0.266667, 0.290196, 0.274510, 0.278431, 0.286275, 0.364706, 0.517647, 0.654902, 0.745098, 0.788235, 0.929412, 0.921569, 0.717647, 0.600000, 0.580392, 0.611765,
    0.960784, 0.945098, 0.976471, 0.952941, 0.937255, 0.917647, 0.803922, 0.600000, 0.392157, 0.250980, 0.211765, 0.219608, 0.235294, 0.247059, 0.247059, 0.239216, 0.223529, 0.223529, 0.298039, 0.411765, 0.396078, 0.305882, 0.243137, 0.207843, 0.235294, 0.407843, 0.556863, 0.419608, 0.349020, 0.309804, 0.313726, 0.364706, 0.482353, 0.639216, 0.745098, 0.843137, 0.945098, 0.839216, 0.619608, 0.560784, 0.588235,
    0.972549, 0.952941, 0.980392, 0.972549, 0.956863, 0.949020, 0.917647, 0.819608, 0.662745, 0.462745, 0.298039, 0.215686, 0.188235, 0.219608, 0.219608, 0.211765, 0.211765, 0.239216, 0.494118, 0.945098, 0.937255, 0.647059, 0.301961, 0.219608, 0.250980, 0.572549, 0.886275, 0.898039, 0.725490, 0.584314, 0.505882, 0.478431, 0.447059, 0.494118, 0.611765, 0.733333, 0.866667, 0.921569, 0.725490, 0.580392, 0.537255,
    0.968627, 0.933333, 0.964706, 0.952941, 0.964706, 0.964706, 0.960784, 0.933333, 0.866667, 0.717647, 0.541176, 0.376471, 0.243137, 0.203922, 0.207843, 0.215686, 0.215686, 0.250980, 0.596078, 1.000000, 1.000000, 0.807843, 0.325490, 0.223529, 0.219608, 0.305882, 0.545098, 0.854902, 0.988235, 1.000000, 1.000000, 0.972549, 0.886275, 0.709804, 0.443137, 0.505882, 0.709804, 0.909804, 0.858824, 0.615686, 0.517647,
    0.941176, 0.917647, 0.913725, 0.925490, 0.937255, 0.945098, 0.941176, 0.952941, 0.937255, 0.878431, 0.780392, 0.643137, 0.458824, 0.301961, 0.223529, 0.215686, 0.219608, 0.258824, 0.541176, 1.000000, 1.000000, 0.596078, 0.266667, 0.219608, 0.219608, 0.247059, 0.298039, 0.392157, 0.509804, 0.615686, 0.713726, 0.831373, 0.870588, 0.650980, 0.341176, 0.341176, 0.505882, 0.803922, 0.921569, 0.741176, 0.525490,
    0.937255, 0.901961, 0.913725, 0.913725, 0.905882, 0.909804, 0.937255, 0.960784, 0.941176, 0.929412, 0.890196, 0.850980, 0.709804, 0.533333, 0.364706, 0.254902, 0.231373, 0.262745, 0.549020, 0.996078, 0.909804, 0.392157, 0.235294, 0.219608, 0.239216, 0.329412, 0.388235, 0.372549, 0.301961, 0.298039, 0.309804, 0.349020, 0.384314, 0.352941, 0.356863, 0.447059, 0.556863, 0.756863, 0.890196, 0.792157, 0.521569,
    0.890196, 0.890196, 0.901961, 0.913725, 0.921569, 0.941176, 0.952941, 0.952941, 0.964706, 0.933333, 0.921569, 0.921569, 0.878431, 0.784314, 0.647059, 0.447059, 0.301961, 0.270588, 0.376471, 0.596078, 0.525490, 0.270588, 0.223529, 0.215686, 0.278431, 0.643137, 0.890196, 0.717647, 0.423529, 0.301961, 0.298039, 0.333333, 0.403922, 0.505882, 0.631373, 0.709804, 0.737255, 0.729412, 0.701961, 0.556863, 0.478431,
    0.780392, 0.776471, 0.780392, 0.803922, 0.831373, 0.862745, 0.878431, 0.913725, 0.929412, 0.937255, 0.952941, 0.933333, 0.925490, 0.898039, 0.854902, 0.745098, 0.552941, 0.356863, 0.286275, 0.270588, 0.262745, 0.235294, 0.207843, 0.211765, 0.254902, 0.498039, 0.921569, 1.000000, 0.827451, 0.494118, 0.513726, 0.619608, 0.705882, 0.737255, 0.713726, 0.658824, 0.588235, 0.470588, 0.360784, 0.356863, 0.431373,
    0.776471, 0.768627, 0.752941, 0.745098, 0.756863, 0.780392, 0.784314, 0.803922, 0.803922, 0.835294, 0.850980, 0.905882, 0.894118, 0.917647, 0.898039, 0.886275, 0.815686, 0.627451, 0.396078, 0.286275, 0.247059, 0.223529, 0.215686, 0.231373, 0.258824, 0.349020, 0.658824, 1.000000, 0.984314, 0.831373, 0.756863, 0.756863, 0.705882, 0.647059, 0.513726, 0.380392, 0.294118, 0.262745, 0.274510, 0.345098, 0.427451,
    0.847059, 0.780392, 0.752941, 0.776471, 0.784314, 0.784314, 0.760784, 0.749020, 0.784314, 0.772549, 0.760784, 0.776471, 0.788235, 0.811765, 0.858824, 0.882353, 0.905882, 0.847059, 0.666667, 0.439216, 0.294118, 0.274510, 0.294118, 0.337255, 0.427451, 0.560784, 0.733333, 0.815686, 0.803922, 0.756863, 0.690196, 0.568627, 0.415686, 0.305882, 0.258824, 0.231373, 0.227451, 0.247059, 0.290196, 0.376471, 0.439216,
    0.768627, 0.760784, 0.760784, 0.772549, 0.807843, 0.792157, 0.760784, 0.749020, 0.784314, 0.807843, 0.752941, 0.760784, 0.749020, 0.764706, 0.772549, 0.803922, 0.823529, 0.870588, 0.850980, 0.733333, 0.490196, 0.427451, 0.552941, 0.650980, 0.701961, 0.745098, 0.709804, 0.674510, 0.588235, 0.454902, 0.345098, 0.270588, 0.243137, 0.219608, 0.211765, 0.203922, 0.223529, 0.337255, 0.486275, 0.427451, 0.450980,
    0.780392, 0.768627, 0.772549, 0.741176, 0.752941, 0.733333, 0.725490, 0.741176, 0.792157, 0.807843, 0.800000, 0.733333, 0.772549, 0.788235, 0.776471, 0.749020, 0.741176, 0.760784, 0.803922, 0.811765, 0.764706, 0.678431, 0.721569, 0.741176, 0.682353, 0.607843, 0.478431, 0.349020, 0.282353, 0.250980, 0.243137, 0.235294, 0.215686, 0.203922, 0.211765, 0.207843, 0.250980, 0.549020, 0.694118, 0.466667, 0.447059,
    0.854902, 0.823529, 0.835294, 0.788235, 0.768627, 0.733333, 0.737255, 0.717647, 0.733333, 0.745098, 0.737255, 0.745098, 0.800000, 0.866667, 0.854902, 0.772549, 0.760784, 0.752941, 0.760784, 0.756863, 0.800000, 0.800000, 0.701961, 0.588235, 0.431373, 0.333333, 0.274510, 0.235294, 0.215686, 0.231373, 0.286275, 0.313726, 0.278431, 0.243137, 0.235294, 0.235294, 0.290196, 0.576471, 0.600000, 0.466667, 0.447059,
    0.858824, 0.866667, 0.878431, 0.874510, 0.831373, 0.756863, 0.737255, 0.749020, 0.760784, 0.760784, 0.745098, 0.760784, 0.760784, 0.776471, 0.784314, 0.780392, 0.819608, 0.835294, 0.760784, 0.768627, 0.772549, 0.776471, 0.764706, 0.635294, 0.419608, 0.298039, 0.250980, 0.219608, 0.219608, 0.266667, 0.592157, 0.678431, 0.525490, 0.400000, 0.364706, 0.341176, 0.329412, 0.427451, 0.482353, 0.450980, 0.454902,
    0.800000, 0.800000, 0.803922, 0.831373, 0.827451, 0.776471, 0.729412, 0.772549, 0.839216, 0.819608, 0.756863, 0.756863, 0.756863, 0.764706, 0.760784, 0.784314, 0.803922, 0.850980, 0.823529, 0.784314, 0.772549, 0.756863, 0.756863, 0.733333, 0.650980, 0.454902, 0.301961, 0.231373, 0.219608, 0.278431, 0.623529, 0.890196, 0.933333, 0.866667, 0.749020, 0.623529, 0.631373, 0.611765, 0.478431, 0.458824, 0.458824,
    0.913725, 0.890196, 0.862745, 0.854902, 0.784314, 0.760784, 0.752941, 0.764706, 0.866667, 0.874510, 0.760784, 0.745098, 0.749020, 0.807843, 0.823529, 0.811765, 0.796078, 0.788235, 0.772549, 0.768627, 0.803922, 0.839216, 0.784314, 0.741176, 0.713726, 0.662745, 0.498039, 0.321569, 0.235294, 0.223529, 0.305882, 0.498039, 0.776471, 0.901961, 0.866667, 0.694118, 0.721569, 0.560784, 0.447059, 0.454902, 0.458824,
    0.850980, 0.850980, 0.866667, 0.839216, 0.772549, 0.741176, 0.737255, 0.745098, 0.784314, 0.780392, 0.760784, 0.745098, 0.756863, 0.866667, 0.929412, 0.917647, 0.866667, 0.815686, 0.792157, 0.776471, 0.764706, 0.823529, 0.839216, 0.776471, 0.717647, 0.694118, 0.670588, 0.552941, 0.376471, 0.250980, 0.239216, 0.286275, 0.494118, 0.662745, 0.623529, 0.517647, 0.388235, 0.364706, 0.419608, 0.458824, 0.458824,
    0.823529, 0.811765, 0.811765, 0.780392, 0.760784, 0.760784, 0.749020, 0.768627, 0.776471, 0.792157, 0.764706, 0.752941, 0.741176, 0.874510, 0.952941, 0.968627, 0.949020, 0.882353, 0.847059, 0.788235, 0.764706, 0.760784, 0.752941, 0.733333, 0.733333, 0.713726, 0.713726, 0.682353, 0.596078, 0.415686, 0.282353, 0.262745, 0.360784, 0.419608, 0.313726, 0.231373, 0.235294, 0.349020, 0.447059, 0.454902, 0.458824,
    0.898039, 0.925490, 0.886275, 0.882353, 0.815686, 0.764706, 0.741176, 0.745098, 0.819608, 0.933333, 0.815686, 0.725490, 0.752941, 0.866667, 0.937255, 0.964706, 0.949020, 0.929412, 0.898039, 0.870588, 0.796078, 0.760784, 0.752941, 0.741176, 0.737255, 0.780392, 0.784314, 0.737255, 0.686275, 0.603922, 0.462745, 0.305882, 0.250980, 0.219608, 0.196078, 0.188235, 0.223529, 0.439216, 0.541176, 0.509804, 0.474510,
    0.827451, 0.831373, 0.831373, 0.839216, 0.788235, 0.760784, 0.764706, 0.745098, 0.772549, 0.870588, 0.815686, 0.752941, 0.741176, 0.807843, 0.905882, 0.952941, 0.945098, 0.921569, 0.925490, 0.905882, 0.898039, 0.843137, 0.784314, 0.752941, 0.745098, 0.776471, 0.807843, 0.749020, 0.709804, 0.670588, 0.611765, 0.482353, 0.321569, 0.223529, 0.184314, 0.184314, 0.258824, 0.490196, 0.607843, 0.607843, 0.541176,
    0.858824, 0.870588, 0.870588, 0.858824, 0.815686, 0.752941, 0.741176, 0.713726, 0.741176, 0.800000, 0.807843, 0.749020, 0.725490, 0.800000, 0.894118, 0.921569, 0.941176, 0.929412, 0.925490, 0.921569, 0.937255, 0.917647, 0.854902, 0.784314, 0.717647, 0.741176, 0.733333, 0.729412, 0.698039, 0.686275, 0.662745, 0.619608, 0.505882, 0.352941, 0.235294, 0.200000, 0.321569, 0.537255, 0.631373, 0.650980, 0.627451,
    0.772549, 0.772549, 0.792157, 0.811765, 0.784314, 0.749020, 0.733333, 0.737255, 0.752941, 0.850980, 0.874510, 0.764706, 0.713726, 0.745098, 0.854902, 0.905882, 0.937255, 0.909804, 0.909804, 0.909804, 0.921569, 0.921569, 0.933333, 0.870588, 0.819608, 0.760784, 0.721569, 0.713726, 0.698039, 0.737255, 0.749020, 0.694118, 0.619608, 0.533333, 0.392157, 0.274510, 0.368627, 0.513726, 0.619608, 0.647059, 0.619608,
    0.733333, 0.760784, 0.737255, 0.729412, 0.741176, 0.733333, 0.721569, 0.713726, 0.733333, 0.764706, 0.784314, 0.733333, 0.713726, 0.713726, 0.823529, 0.890196, 0.874510, 0.878431, 0.890196, 0.898039, 0.929412, 0.941176, 0.933333, 0.925490, 0.890196, 0.819608, 0.752941, 0.729412, 0.729412, 0.729412, 0.745098, 0.701961, 0.670588, 0.631373, 0.572549, 0.466667, 0.431373, 0.447059, 0.541176, 0.607843, 0.600000,
    0.729412, 0.745098, 0.733333, 0.733333, 0.741176, 0.725490, 0.717647, 0.713726, 0.725490, 0.733333, 0.811765, 0.764706, 0.698039, 0.733333, 0.815686, 0.890196, 0.862745, 0.831373, 0.854902, 0.913725, 0.929412, 0.929412, 0.937255, 0.921569, 0.909804, 0.898039, 0.831373, 0.745098, 0.733333, 0.717647, 0.698039, 0.678431, 0.686275, 0.682353, 0.658824, 0.584314, 0.474510, 0.396078, 0.447059, 0.501961, 0.552941,
    0.752941, 0.749020, 0.749020, 0.737255, 0.741176, 0.737255, 0.709804, 0.713726, 0.725490, 0.768627, 0.909804, 0.831373, 0.717647, 0.725490, 0.784314, 0.866667, 0.878431, 0.878431, 0.890196, 0.901961, 0.941176, 0.921569, 0.937255, 0.937255, 0.933333, 0.917647, 0.878431, 0.835294, 0.776471, 0.737255, 0.698039, 0.690196, 0.678431, 0.701961, 0.729412, 0.611765, 0.533333, 0.392157, 0.368627, 0.462745, 0.505882
    });
    return patch;
}

void testPixel()
{
    Pixel color({1,1,1});
    if(color.rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Pixel px2(color);
    if(px2.rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::vector<Pixel> img(10, px2);
    if(img.back().rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::vector<Pixel> img2(img);
    img.insert(img.end(), img2.begin(), img2.begin() + 5);
    if(img.back().rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));


    {
        auto mat_a = random::uniform<Pixel>({3,3});
        if(mat_a(0,0).r() > 1)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testImageResize()
{
    {// case 1
        Matrix<Pixel> img({2,2});
        img(0,0) = Pixel::black();
        img(0,1) = Pixel::white();
        img(1,0) = Pixel::black();
        img(1,1) = Pixel::white();

        auto img3x3 = resize(img, Shape({3,3}));
        if(img3x3(1,1) != Pixel({.5, .5, .5}))
        {
            std::cout << to_string(img3x3,2) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {//case 2
        Matrix<float> img({2,2},{0,1,2,3});
        Matrix<float> expected = Matrix<float>::ones({8,8});
        expected(Block({0,4},{0,4}))*=0;
        expected(Block({0,4},{4,8}))*=1;
        expected(Block({4,8},{0,4}))*=2;
        expected(Block({4,8},{4,8}))*=3;

        if(mxm::resize(img, 4, "nearest") != expected)
        {
            std::cout << to_string(mxm::resize(img, 4, "nearest"),2) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

void testPixelMemory()
{
    std::vector<Pixel> pxs;
    pxs.push_back(Pixel({0.5,0.25,0.125}));
    pxs.push_back(Pixel({0.5,0.25,0.125}));
    pxs.push_back(Pixel({0.5,0.25,0.125}));

    for(size_t i = 0;i < 9; i+= Pixel::size())
    {
        if(
            ((float*)pxs.data())[i+0] != 0.5 ||
            ((float*)pxs.data())[i+1] != 0.25 ||
            ((float*)pxs.data())[i+2] != 0.125)
        {
            std::cout << ((float*)pxs.data())[i+0] << ","
                << ((float*)pxs.data())[i+1] << ","
                << ((float*)pxs.data())[i+2] << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

void testKernels()
{
    if(0){//case 01
        std::cout << mxm::to_string(kernel::gauss<float>(3)) << std::endl;
        std::cout << mxm::sum(kernel::gauss<float>(3)) << std::endl;
    }
}

void testBresenhamCircle()
{
    {
        // auto ret = halfQuarterBresenhamCircle(8.3);
        // expendFullBresenhamCircle(ret);
        // auto ret = bresenhamCircle(8.3);
        auto ret = bresenhamCircle(6);
        // std::cout << "ret size: " << ret.shape(1) << std::endl;
        // std::cout << mxm::to_string(ret) << std::endl;
    }

    {
        // auto cnt = 0;
        // traverseBresenhamCircleArea(3, {5,5}, [&](auto i, auto j){
        //     std::cout << cnt++ << " " ;//<< std::endl;
        //     std::cout << "(" << i << "," <<  j << ")" << std::endl;
        //     });
    }
}

void testBriefDescriptor()
{
    if(1){

        Matrix<float> raw(convolute(dataImagePatch(), kernel::gauss(3)));
        auto coord_buff = fastCorners(raw, 0.02);
        Matrix<size_t> fp_0 = gridPartitionNonMaximalSuppression(coord_buff, 1, 100);
        auto desc = calculateBriefDescriptor<8>(raw, fp_0);

        Matrix<float> rotated = flip(raw).T();
        coord_buff = fastCorners(rotated, 0.02);
        Matrix<size_t> fp = gridPartitionNonMaximalSuppression(coord_buff, 1, 100);
        auto desc_r = calculateBriefDescriptor<8>(rotated, fp);

        float delta_orientation = SO2ToAngle(
            mxm::rodrigues2D(desc(0,0).orientation()).T().matmul(
                mxm::rodrigues2D(desc_r(0,0).orientation())));
        if( abs(delta_orientation - M_PI_2) > 1 * eps() )
        {
            std::cout << delta_orientation << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(desc_r(0,0).distance(desc(0,0)) > 0.2)
        {
            std::cout << "raw: " << mxm::to_string(fp_0.T()) << mxm::to_string(desc(0,0)) << std::endl;
            std::cout << "rotated: " << mxm::to_string(fp.T()) << mxm::to_string(desc_r(0,0)) << std::endl;
            std::cout << "distance: " << desc_r(0,0).distance(desc(0,0)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    // std::cout << mxm::to_string(flip (Matrix<float>::identity(3))) << std::endl;
}

void testHomography()
{
    {
        Matrix<float> pts_src({2,4},{2,3, 4,1, 5,5, 3,6}, COL);
        Matrix<float> pts_dst({2,4},{1,1, 3,1, 3,3, 1,3}, COL);
        Matrix<float> homo({1,4},{1,1,1,1});
        Matrix<float> expected = vstack(pts_dst, homo);
        auto mat_h = findHomographyMatrix(pts_src, pts_dst);
        // std::cout << mxm::to_string(mat_h) << std::endl;
        // return ;
        auto result = mat_h.matmul( vstack(pts_src, homo) );
        for(size_t i = 0; i < 4; i++) result(Col(i)) *= (1./result(2, i));

        if(norm(result - expected) > std::numeric_limits<float>::epsilon() * 50)
        {
            std::cout << "error norm:" << norm(result - expected) << std::endl;
            std::cout << "WARNING: todo fix.\n" << std::string(__FILE__) + ":" + std::to_string(__LINE__) << std::endl;
            if(norm(result - expected) > std::numeric_limits<float>::epsilon() * 150)
            {
                std::cout << mxm::to_string(mat_h) << std::endl;
                std::cout << mxm::to_string(result) << std::endl;
                std::cout << norm(result - expected) << std::endl;
                throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
            }
        }
    }

    {
        Matrix<float> pts_src(fixRow(2),{2,3, 3,2, 4,1, 5,5, 4,5.5, 3,6}, COL);
        Matrix<float> pts_dst(fixRow(2),{1,1, 2,1, 3,1, 3,3, 2,3.0, 1,3}, COL);
        float error = 0.02;
        pts_dst += (random::uniform<float>(pts_dst.shape()) * error);
        Matrix<float> homo = Matrix<float>::ones({1, pts_dst.shape(1)});
        Matrix<float> expected = vstack(pts_dst, homo);

        auto mat_h = findHomographyMatrix(pts_src, pts_dst);
        auto result = mat_h.matmul( vstack(pts_src, homo) );
        for(size_t i = 0; i < pts_dst.shape(1); i++) result(Col(i)) *= (1./result(2, i));

        if(norm(result - expected) > error * 10 )
        {
            std::cout << (norm(result - expected)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

void testICP()
{
    {
        size_t pt_num = 30;
        size_t dim = 3;
        Rotation<float> rot = Rotation<float>::fromAxisAngle({1,0,0}, M_PI * 0.25);
        auto pts = random::uniform<float>({dim, pt_num});
        for(size_t i = 0; i < dim; i++)
        {
            pts(Row(i)) -= mxm::sum(pts(Row(i))) / pt_num;
        }
        auto pts2 = rot.apply(pts);
        auto result = icpFindRotation(pts, pts2);
        if(!isIdentity(rot.asMatrix().matmul(result.asMatrix().T()), nullptr, 5 * std::numeric_limits<float>::epsilon()))
        {
            std::cout << mxm::to_string(rot.asMatrix()) << std::endl;
            std::cout << mxm::to_string(result.asMatrix()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

void testEpipolarGeometry()
{
#if 0
    { // test eight point
        const size_t PT_NUM = 8;
        Camera<float> cam;
        cam.setFocalLength({500,500}).setPrincipalOffset({240, 320});
        Matrix<float> pts({3,PT_NUM});
        for(uint8_t i = 0; i < PT_NUM; i++)
        {
            pts(Col(i)) = binaryToVector<float>(3, i);
        }
        pts = Rotation<float>::fromAxisAngle({1.,1,1}, M_PI / 4.).apply(pts);
        pts(Row(2)) += 5;
        cam.setPosition({-0.5, 0, 0});
        auto pose1 = cam.pose();
        auto pts1 = cam.project(pts);
        auto homo1 = cam.pose().inv().apply(pts);
        homo1 /= homo1(Row(2));
        auto mat_norm1 = findNormalizeMatrix(pts1);
        auto homo_px1 = vstack(pts1, decltype(pts1)::ones({1, pts1.shape(1)}));
        auto norm_pts1 = mat_norm1.matmul(homo_px1);

        cam.setPosition({0.5, 0, 0});
        auto pose2 = cam.pose();
        auto pts2 = cam.project(pts);
        auto homo2 = cam.pose().inv().apply(pts);
        homo2 /= homo2(Row(2));
        auto mat_norm2 = findNormalizeMatrix(pts2);
        auto homo_px2 = vstack(pts2, decltype(pts2)::ones({1, pts2.shape(1)}));
        auto norm_pts2 = mat_norm2.matmul(homo_px2);

        auto delta_pose = pose1.inv() * pose2;

        auto actual_essential = so::wedge(delta_pose.translation()).matmul(delta_pose.rotation().asMatrix());
        auto actual_fundamental = cam.invMatrix().T().matmul(actual_essential).matmul(cam.invMatrix());

        std::cout << "pts1: \n" << mxm::to_string(pts1) << std::endl;
        std::cout << "pts2: \n" << mxm::to_string(pts2) << std::endl;
        std::cout << "All data ready." << std::endl;

        std::cout << "norm_pts1: \n" << mxm::to_string(norm_pts1) << std::endl;
        std::cout << "norm_pts2: \n" << mxm::to_string(norm_pts2) << std::endl;


        Matrix<float> mat_f_normalized = epipolarEightPoint<double>(norm_pts1, norm_pts2);

        auto expect_zero = checkEpipolarConstraints(mat_f_normalized, norm_pts1, norm_pts2);
        std::cout << "constraits:\n" << mxm::to_string(expect_zero, 10) << std::endl;

        auto u_d_vh = svd(mat_f_normalized);
        std::cout << "mat_f_normalized svd: \n" << mxm::to_string(u_d_vh) << std::endl;
        u_d_vh[1](2,0) = 0;
        mat_f_normalized = u_d_vh[0].matmul(diagonalMatrix(u_d_vh[1])).matmul(u_d_vh[2]);

        auto mat_f = mat_norm2.T().matmul(mat_f_normalized).matmul(mat_norm1);
        expect_zero = checkEpipolarConstraints(mat_f, homo_px1, homo_px2);
        std::cout << "constraits:\n" << mxm::to_string(expect_zero) << std::endl;


        std::cout << "mat_f:\n" << mxm::to_string(mat_f) << std::endl;
        std::cout << "actual mat_f:\n" << mxm::to_string(actual_fundamental) << std::endl;

        // auto u_d_vh = svd(mat_f);
        // std::cout << "fundamental svd: \n" << mxm::to_string(u_d_vh) << std::endl;
        // u_d_vh[1](2,0) = 0;
        // mat_f = u_d_vh[0].matmul(diagonalMatrix(u_d_vh[1])).matmul(u_d_vh[2]);
        // std::cout << mxm::to_string(u_d_vh[0].matmul(diagonalMatrix(u_d_vh[1])).matmul(u_d_vh[2]) - mat_f) << std::endl;

        // mat_f *= (1.f / mat_f(2,2));




        std::cout << "mat_f projected:\n" << mxm::to_string(mat_f) << std::endl;
        expect_zero = checkEpipolarConstraints(mat_f, homo_px1, homo_px2);
        std::cout << "final error:\n" << mxm::to_string(expect_zero, 10) << std::endl;
    }
#endif


    {
        // test eight point
        // two questions remained to be answered:
        // 1. After replacing the singular values of the solution of Direct Linear Method with [1,1,0],
        //    does the essential matrix that projected to the manifold still preserves the epipolar constraints?
        // 2. While extracting the rotation and translation from the essential matrix, the rotation was
        //    composed from R = (U R_z V^h). R_z was gauranteed that det(R_z) = 1, but the orthogonal matrices U and V
        //    have determinant 1 or -1. Therefore, it's possible that the final rotation matrix R is not in SO3.
        //    It happens when det(U) = -1, det(V) = 1 and det(R) = -1.
#if 0
        const size_t PT_NUM = 8;
        Camera<float> cam;
        cam.setFocalLength({500,500}).setPrincipalOffset({240, 320});
        Matrix<float> pts({3,PT_NUM});
        for(uint8_t i = 0; i < PT_NUM; i++)
        {
            pts(Col(i)) = binaryToVector<float>(3, i);
        }
        // std::cout << "original points:\n" << mxm::to_string(pts) << std::endl;
        pts = Rotation<float>::fromAxisAngle({1.,1,1}, M_PI / 4.).apply(pts);

        pts(Row(2)) += 5;

        cam.setPosition({-0.5, 0, 0});
        auto pose1 = cam.pose();
        auto pts1 = cam.project(pts);
        auto homo1 = cam.pose().inv().apply(pts);
        homo1 /= homo1(Row(2));
        auto mat_norm1 = findNormalizeMatrix(homo1(Block({0, end()-1}, {})));
        auto norm_pts1 = mat_norm1.matmul(homo1);

        cam.setPosition({0.5, 0, 0});
        auto pose2 = cam.pose();
        auto pts2 = cam.project(pts);
        auto homo2 = cam.pose().inv().apply(pts);
        homo2 /= homo2(Row(2));
        auto mat_norm2 = findNormalizeMatrix(homo2(Block({0, end()-1}, {})));
        auto norm_pts2 = mat_norm2.matmul(homo2);

        auto actual_essential = so::wedge(Vector<float>{1.f, 0, 0});

        std::cout << "pts1: \n" << mxm::to_string(pts1) << std::endl;
        std::cout << "pts2: \n" << mxm::to_string(pts2) << std::endl;
        std::cout << "All data ready." << std::endl;

        auto mat_e_normalized = epipolarEightPoint(norm_pts1, norm_pts2);

        auto expect_zero = checkEpipolarConstraints(mat_e_normalized, norm_pts1, norm_pts2);
        std::cout << "constraits: " << mxm::to_string(expect_zero) << std::endl;
        auto mat_e = mat_norm2.T().matmul(mat_e_normalized).matmul(mat_norm1);
        expect_zero = checkEpipolarConstraints(mat_e, homo1, homo2);
        std::cout << "constraits: " << mxm::to_string(expect_zero) << std::endl;

        std::cout << "mat_e:\n" << mxm::to_string(mat_e) << std::endl;

        expect_zero = checkEpipolarConstraints(mat_e, homo1, homo2);
        std::cout << "constraits: " << mxm::to_string(expect_zero) << std::endl;

        std::cout << "mat_e:\n" << mxm::to_string(mat_e) << std::endl;
        std::cout << "actual_essential:\n" << mxm::to_string(actual_essential) << std::endl;

        auto u_d_vh = svd(mat_e);
        std::cout << "essential svd: \n" << mxm::to_string(u_d_vh) << std::endl;

        Matrix<float> rot_zp({3,3},{0,1,0, -1,0,0, 0,0,1}, ROW);
        auto rot_zn = rot_zp.T();

        auto r1 = u_d_vh[0].matmul(rot_zp).matmul(u_d_vh[2]);
        auto r2 = u_d_vh[0].matmul(rot_zn).matmul(u_d_vh[2]);
        auto t1 = so::vee(u_d_vh[0].matmul(rot_zp).matmul(u_d_vh[0].T()));
        auto t2 = so::vee(u_d_vh[0].matmul(rot_zn).matmul(u_d_vh[0].T()));
        std::cout << "r1:\n" << mxm::to_string(r1) << std::endl;
        std::cout << "r2:\n" << mxm::to_string(r2) << std::endl;
        std::cout << "t1:\n" << mxm::to_string(t1) << std::endl;
        std::cout << "t2:\n" << mxm::to_string(t2) << std::endl;
        std::cout << "det:\n" << mxm::det(u_d_vh[0].matmul(u_d_vh[2])) << std::endl;
        std::cout << "angle: " << SO::findAngle<3>(r2) * 180 / M_PI << std::endl;
        // auto sigma = (u_d_vh[1](0,0) + u_d_vh[1](1,0)) * 0.5f;
        // u_d_vh[1](2,0) = 0;
        // mat_e = u_d_vh[0].matmul(diagonalMatrix(u_d_vh[1])).matmul(u_d_vh[2]);
#endif
    }

#if 0
    // test normalize points
    {
        Matrix<float> pts1(fixRow(2), {
            190., 272.49231, 166.38622, 239.43414, 240.50644, 313.77134, 215.85516, 281.66667,
            320., 373.93882, 393.07985, 446.13581, 293.2445 , 337.77011, 359.14993, 403.33333});
        auto homo_pts1 = vstack(pts1, decltype(pts1)::ones({1, pts1.shape(1)}));
        // Matrix<float> pts2(fixRow(2), {
        //     290., 379.11615, 257.19823, 335.67568, 326.64298, 404.77772, 295.08671, 365.,
        //     320., 373.93882, 393.07985, 446.13581, 293.2445 , 337.77011, 359.14993, 403.33333});

        auto norm_mat = findNormalizeMatrix(pts1);
        float error(0);
        if(!isZero(invDiagHomogeneous(norm_mat).matmul(norm_mat.matmul(homo_pts1)) - homo_pts1, &error))
        {
            std::cout << "===error: " << error << std::endl;
        }

    }
#endif

#if 1
    {
        using DType = double;
        const size_t PT_NUM = 8;
        Camera<DType> cam;
        cam.setFocalLength({500,500}).setPrincipalOffset({0, 0});
        Matrix<DType> pts({3,PT_NUM});
        for(uint8_t i = 0; i < PT_NUM; i++)
        {
            pts(Col(i)) = binaryToVector<DType>(3, i);
        }
        pts *= 2.f;
        pts = Rotation<DType>::fromAxisAngle({1.,1,1}, M_PI / 4.).apply(pts);

        pts(Row(2)) += 5;

        cam.setPosition({0, 0, 0});
        auto pose1 = cam.pose();
        auto pts1 = cam.project(pts);
        auto homo1 = cam.pose().inv().apply(pts);
        homo1 /= homo1(Row(2));

        cam.setPosition({1, 0.2, -0.4});
        cam.setOrientation(Rotation<DType>::fromAxisAngle({1,0,1}, -0.3f));
        auto pose2 = cam.pose();
        auto pts2 = cam.project(pts);
        auto homo2 = cam.pose().inv().apply(pts);
        homo2 /= homo2(Row(2));

        // std::cout << "pts1: \n" << mxm::to_string(pts1) << std::endl;
        // std::cout << "pts2: \n" << mxm::to_string(pts2) << std::endl;

        auto good_guess = Vector<DType>{-1.f,0,0};
        auto local_minimum_guess = Vector<DType>{-0.f,0,1};
        auto t_r = epipolarLeastSquare(homo1, homo2, good_guess);
        auto rotation_expected = (pose1.rotation() * pose2.rotation().inv()).asMatrix();
        auto translation_expected = (pose1 * pose2.inv()).translation().normalized();

        DType t_error(0);
        DType r_error(0);
        DType tol(1e-8);

        if(!isZero(so::wedge(t_r[0]).matmul(translation_expected), &t_error, tol)
        || !isZero(rotation_expected - t_r[1], &r_error, tol))
        {
            std::cout << "t_error: " << t_error << std::endl;
            std::cout << "expected:\n" << mxm::to_string(translation_expected) << std::endl;
            std::cout << "get:\n" << mxm::to_string(t_r[0]) << std::endl;

            std::cout << "r_error: " << r_error << std::endl;
            std::cout << "expected:\n" << mxm::to_string(rotation_expected) << std::endl;
            std::cout << "get:\n" << mxm::to_string(t_r[1]) << std::endl;

            auto cost = checkEpipolarConstraints(so::wedge(t_r[0]).matmul(t_r[1]) , homo1, homo2);
            std::cout << "residual:\n" << mxm::to_string (cost) << std::endl;
            auto cost_true = checkEpipolarConstraints(so::wedge(translation_expected).matmul(rotation_expected) , homo1, homo2);
            std::cout << "residual ref:\n" << mxm::to_string (cost) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
#endif
}


void testPnP()
{
#if 1
    std::vector<size_t> cali_board_reso{8,8};
    Matrix<float> pts3d({3, cali_board_reso.at(0) * cali_board_reso.at(1)});
    pts3d.setBlock(0,0, generateNCubeVertices({5e-1, 5e-1}, cali_board_reso) );
    Camera<float, 3> cam;
    cam.setResolution({500, 500}).setFov({M_PI /3,M_PI /3}).setPosition({0, 0, -2});

    Matrix<float> pts2d = cam.project(pts3d);
    // std::cout << mxm::to_string(pts3d) << std::endl;
    // std::cout << mxm::to_string(pts2d) << std::endl;

    PerspectiveNPointOptimizer<float> problem(pts3d, pts2d, cam);
    problem.initialGuess(Vector<float>{0.2,0.2, -2.1}, Vector<float>{0.1, -0.02, 0.05});
    problem.solve(6, 0, "gn");

    if(! isZero(problem.tState() - cam.pose().translation(), nullptr, eps<float>() * 5))
    {
        std::cout << mxm::to_string(problem.tState()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(! isZero(problem.rState() - so::vee(SO::log<3>(cam.pose().rotation().asMatrix())), nullptr, eps<float>() * 5))
    {
        std::cout << mxm::to_string(problem.rState()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

#endif
}

void testPinholeCalibration()
{
#if 1
    using DType = double;
    std::vector<size_t> cali_board_reso{8,8};
    size_t frame_num = 4;
    Vector<RigidTransform<DType, 3>> poses(frame_num);
    poses(0).translation() = Vector<DType>{0.0,0.0,-2};
    poses(1).translation() = Vector<DType>{0.0,0.2,-2};
    poses(2).translation() = Vector<DType>{0.2,0.0,-2};
    poses(3).translation() = Vector<DType>{0.1,0.3,-2};

    poses(0).rotation() = Rotation<DType, 3>::fromAxisAngle({0,0,1}, 0.0);
    poses(1).rotation() = Rotation<DType, 3>::fromAxisAngle({-1,0,0}, 0.1);
    poses(2).rotation() = Rotation<DType, 3>::fromAxisAngle({0,1,0}, 0.1);
    poses(3).rotation() = Rotation<DType, 3>::fromAxisAngle({0,1,0}, 0.1);

    Matrix<DType> pts3d({3, cali_board_reso.at(0) * cali_board_reso.at(1)});
    pts3d.setBlock(0,0, generateNCubeVertices({5e-1, 5e-1}, cali_board_reso) );
    Camera<DType, 3> cam;
    auto p_distor = Distortion<DType>::radialTangential({-0.28340811, 0.07395907, 0.00019359, 1.76187114e-05});
    cam.setDistortion(p_distor);
    cam.setFocalLength({457.296, 458.654}).setResolution({480, 752}).setPrincipalOffset({248.375, 367.215});

    Vector<Matrix<DType>> pts2d(frame_num);
    for(size_t i = 0; i < frame_num; i++)
    {
        cam.setPose(poses(i));
        pts2d(i) = cam.project(pts3d);
    }
    // std::cout << mxm::to_string(pts3d) << std::endl;
    // std::cout << mxm::to_string(pts2d) << std::endl;
    Vector<RigidTransform<DType, 3>> pose_guess(poses);
    for(size_t i = 0; i < frame_num; i++)
    {
        pose_guess(i) = RigidTransform<DType, 3>::identity();
        pose_guess(i).translation()(2) = -2;
    }

    Camera<DType, 3> cam_guess(cam);
    cam_guess.setFocalLength(Vector<DType>{500, 500});
    cam_guess.setPrincipalOffset(Vector<DType>{250, 0.5 * 752});
    cam_guess.setDistortion(Distortion<DType>::radialTangential({0,0,0,0}));

    PinholeCameraIntrinsicEstimator<DType> problem(pts3d, pts2d);
    problem.initialGuess(pose_guess, cam_guess);
    problem.solve(5, 0, "gn");

    if(!isZero(cam.focalLength() - problem.camera().focalLength(), nullptr, 1e-4))
    {
        std::cout << "f: " << mxm::to_string(problem.camera().focalLength().T()) << std::endl;
        std::cout << "expected: " << mxm::to_string(cam.focalLength().T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(!isZero(cam.principalOffset() - problem.camera().principalOffset(), nullptr, 1e-4))
    {
        std::cout << "c: " << mxm::to_string(problem.camera().principalOffset().T()) << std::endl;
        std::cout << "expected: " << mxm::to_string(cam.principalOffset().T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    auto distort_result = *static_cast< RadialTangentialDistortion<DType>*>(problem.camera().distortion().get());
    auto distort_expect = *static_cast< RadialTangentialDistortion<DType>*>(cam.distortion().get());

    if(!isZero(distort_result.k() - distort_expect.k(), nullptr, 1e-5))
    {
        std::cout << "k: " << mxm::to_string(distort_result.k().T()) << std::endl;
        std::cout << "expected: " << mxm::to_string(distort_expect.k().T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(!isZero(distort_result.p() - distort_expect.p(), nullptr, 1e-8))
    {
        std::cout << "k: " << mxm::to_string(distort_result.p().T()) << std::endl;
        std::cout << "expected: " << mxm::to_string(distort_expect.p().T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

#endif
}


void testCvBasic()
{
    testPixel();
    testImageResize();
    // testPixelMemory();
    testKernels();
    testBresenhamCircle();
    testBriefDescriptor();
    testHomography();
    testICP();
    testEpipolarGeometry();
    testPnP();
    testPinholeCalibration();
}
#else
void testCvBasic(){}
#endif


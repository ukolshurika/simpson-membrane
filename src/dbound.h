#ifndef DBOUND_H_
#define DBOUND_H_

#include <vector>
#include "membrane.h"

using namespace std;
class DBound{
public:

vector<double>  yc = {0.0231837, 0.0256211, 0.0280518, 0.030476, 0.0328936, 0.0353047, 
0.0377091, 0.0401071, 0.0424985, 0.0448834, 0.0472618, 0.0496337, 
0.0519991, 0.054358, 0.0567105, 0.0590565, 0.0613961, 0.0637293, 
0.0660561, 0.0683764, 0.0706904, 0.072998, 0.0752993, 0.0775941, 
0.0798827, 0.0821649, 0.0844408, 0.0867104, 0.0889737, 0.0912307, 
0.0934815, 0.095726, 0.0979642, 0.100196, 0.102422, 0.104642, 
0.106855, 0.109062, 0.111263, 0.113458, 0.115647, 0.117829, 0.120006, 
0.122176, 0.12434, 0.126498, 0.12865, 0.130796, 0.132936, 0.13507, 
0.137197, 0.139319, 0.141435, 0.143544, 0.145648, 0.147746, 0.149837, 
0.151923, 0.154003, 0.156077, 0.158145, 0.160207, 0.162263, 0.164313, 
0.166357, 0.168396, 0.170428, 0.172455, 0.174476, 0.176491, 0.1785, 
0.180503, 0.182501, 0.184493, 0.186479, 0.188459, 0.190433, 0.192402, 
0.194365, 0.196322, 0.198274, 0.20022, 0.20216, 0.204094, 0.206023, 
0.207946, 0.209864, 0.211776, 0.213682, 0.215582, 0.217477, 0.219367, 
0.22125, 0.223129, 0.225001, 0.226868, 0.22873, 0.230586, 0.232436, 
0.234281, 0.236121, 0.237954, 0.239783, 0.241606, 0.243423, 0.245235, 
0.247042, 0.248843, 0.250639, 0.252429, 0.254214, 0.255994, 0.257768, 
0.259537, 0.2613, 0.263058, 0.264811, 0.266558, 0.2683, 0.270037, 
0.271768, 0.273495, 0.275215, 0.276931, 0.278641, 0.280346, 0.282046, 
0.283741, 0.28543, 0.287114, 0.288793, 0.290467, 0.292136, 0.293799, 
0.295457, 0.297111, 0.298759, 0.300401, 0.302039, 0.303672, 0.305299, 
0.306922, 0.308539, 0.310151, 0.311758, 0.31336, 0.314957, 0.31655, 
0.318137, 0.319719, 0.321296, 0.322867, 0.324435, 0.325997, 0.327554, 
0.329106, 0.330653, 0.332195, 0.333732, 0.335265, 0.336792, 0.338315, 
0.339832, 0.341345, 0.342853, 0.344356, 0.345854, 0.347348, 0.348836, 
0.35032, 0.351799, 0.353273, 0.354742, 0.356207, 0.357666, 0.359121, 
0.360571, 0.362017, 0.363458, 0.364894, 0.366325, 0.367751, 0.369173, 
0.37059, 0.372003, 0.37341, 0.374814, 0.376212, 0.377606, 0.378995, 
0.38038, 0.381759, 0.383135, 0.384505, 0.385872, 0.387233, 0.38859, 
0.389942, 0.39129, 0.392633, 0.393972, 0.395306, 0.396636, 0.397961, 
0.399282, 0.400598, 0.40191, 0.403217, 0.40452, 0.405818, 0.407112, 
0.408401, 0.409686, 0.410967, 0.412243, 0.413515, 0.414782, 0.416045, 
0.417303, 0.418558, 0.419807, 0.421053, 0.422294, 0.423531, 0.424763, 
0.425992, 0.427215, 0.428435, 0.42965, 0.430861, 0.432068, 0.43327, 
0.434468, 0.435662, 0.436852, 0.438038, 0.439219, 0.440396, 0.441569, 
0.442737, 0.443902, 0.445062, 0.446218, 0.44737, 0.448518, 0.449661, 
0.450801, 0.451936, 0.453067, 0.454194, 0.455317, 0.456436, 0.457551, 
0.458662, 0.459769, 0.460871, 0.46197, 0.463064, 0.464155, 0.465241, 
0.466324, 0.467402, 0.468476, 0.469547, 0.470613, 0.471676, 0.472734, 
0.473788, 0.474839, 0.475886, 0.476928, 0.477967, 0.479002, 0.480033, 
0.48106, 0.482083, 0.483102, 0.484117, 0.485129, 0.486136, 0.48714, 
0.48814, 0.489136, 0.490128, 0.491116, 0.492101, 0.493081, 0.494058, 
0.495031, 0.496001, 0.496966, 0.497928, 0.498886, 0.49984, 0.500791, 
0.501737, 0.50268, 0.50362, 0.504555, 0.505487, 0.506415, 0.507339, 
0.50826, 0.509177, 0.510091, 0.511, 0.511906, 0.512809, 0.513708, 
0.514603, 0.515494, 0.516382, 0.517266, 0.518147, 0.519024, 0.519898, 
0.520768, 0.521634, 0.522497, 0.523356, 0.524212, 0.525064, 0.525912, 
0.526757, 0.527599, 0.528437, 0.529271, 0.530102, 0.53093, 0.531754, 
0.532574, 0.533391, 0.534205, 0.535015, 0.535822, 0.536625, 0.537425, 
0.538221, 0.539014, 0.539803, 0.54059, 0.541372, 0.542152, 0.542928, 
0.5437, 0.544469, 0.545235, 0.545998, 0.546757, 0.547512, 0.548265, 
0.549014, 0.54976, 0.550502, 0.551241, 0.551977, 0.55271, 0.553439, 
0.554165, 0.554888, 0.555607, 0.556323, 0.557036, 0.557746, 0.558452, 
0.559156, 0.559856, 0.560552, 0.561246, 0.561936, 0.562623, 0.563307, 
0.563988, 0.564666, 0.56534, 0.566011, 0.566679, 0.567344, 0.568006, 
0.568665, 0.56932, 0.569973, 0.570622, 0.571268, 0.571911, 0.572551, 
0.573188, 0.573822, 0.574452, 0.57508, 0.575704, 0.576326, 0.576944, 
0.57756, 0.578172, 0.578781, 0.579388, 0.579991, 0.580591, 0.581188, 
0.581782, 0.582374, 0.582962, 0.583547, 0.584129, 0.584709, 0.585285, 
0.585858, 0.586429, 0.586996, 0.587561, 0.588122, 0.588681, 0.589237, 
0.58979, 0.59034, 0.590887, 0.591431, 0.591972, 0.59251, 0.593046, 
0.593578, 0.594108, 0.594635, 0.595159, 0.59568, 0.596199, 0.596714, 
0.597227, 0.597736, 0.598243, 0.598748, 0.599249, 0.599748, 0.600243, 
0.600736, 0.601227, 0.601714, 0.602199, 0.602681, 0.60316, 0.603636, 
0.60411, 0.604581, 0.605049, 0.605514, 0.605977, 0.606437, 0.606894, 
0.607349, 0.607801, 0.60825, 0.608697, 0.60914, 0.609581, 0.61002, 
0.610456, 0.610889, 0.611319, 0.611747, 0.612172, 0.612595, 0.613015, 
0.613432, 0.613847, 0.614259, 0.614668, 0.615075, 0.615479, 0.615881, 
0.61628, 0.616676, 0.61707, 0.617461, 0.61785, 0.618236, 0.61862, 
0.619001, 0.619379, 0.619755, 0.620129, 0.6205, 0.620868, 0.621234, 
0.621597, 0.621958, 0.622316, 0.622672, 0.623026, 0.623376, 0.623725, 
0.624071, 0.624414, 0.624755, 0.625094, 0.62543, 0.625763, 0.626094, 
0.626423, 0.626749, 0.627073, 0.627395, 0.627714, 0.62803, 0.628344, 
0.628656, 0.628965, 0.629272, 0.629577, 0.629879, 0.630179, 0.630476, 
0.630771, 0.631064, 0.631355, 0.631643, 0.631928, 0.632211, 0.632492, 
0.632771, 0.633047, 0.633321, 0.633593, 0.633862, 0.634129, 0.634394, 
0.634656, 0.634917, 0.635174, 0.63543, 0.635683, 0.635934, 0.636183, 
0.636429, 0.636674, 0.636916, 0.637155, 0.637393, 0.637628, 0.637861, 
0.638092, 0.63832, 0.638546, 0.63877, 0.638992, 0.639212, 0.639429, 
0.639644, 0.639857, 0.640068, 0.640277, 0.640483, 0.640688, 0.64089, 
0.64109, 0.641287, 0.641483, 0.641677, 0.641868, 0.642057, 0.642244, 
0.642429, 0.642612, 0.642792, 0.642971, 0.643147, 0.643321, 0.643494, 
0.643664, 0.643832, 0.643997, 0.644161, 0.644323, 0.644482, 0.64464, 
0.644795, 0.644949, 0.6451, 0.645249, 0.645396, 0.645542, 0.645685, 
0.645826, 0.645965, 0.646102, 0.646237, 0.646369, 0.6465, 0.646629, 
0.646756, 0.646881, 0.647004, 0.647124, 0.647243, 0.64736, 0.647475, 
0.647588, 0.647699, 0.647807, 0.647914, 0.648019, 0.648122, 0.648223, 
0.648322, 0.648419, 0.648515, 0.648608, 0.648699, 0.648788, 0.648876, 
0.648961, 0.649045, 0.649126, 0.649206, 0.649284, 0.64936, 0.649434, 
0.649506, 0.649576, 0.649644, 0.649711, 0.649775, 0.649838, 0.649898, 
0.649957, 0.650014, 0.650069, 0.650123, 0.650174, 0.650224, 0.650271, 
0.650317, 0.650361, 0.650403, 0.650444, 0.650482, 0.650519, 0.650554, 
0.650587, 0.650618, 0.650647, 0.650675, 0.650701, 0.650725, 0.650747, 
0.650767, 0.650786, 0.650803, 0.650818, 0.650831, 0.650843, 0.650852, 
0.65086, 0.650867, 0.650871, 0.650874, 0.650874, 0.650874};

vector<double> xc = {0.00265231, 0.00464052, 0.00662231, 0.00859768, 0.0105667, 
0.0125292, 0.0144855, 0.0164353, 0.0183788, 0.020316, 0.0222468, 
0.0241714, 0.0260896, 0.0280015, 0.0299072, 0.0318067, 0.0336998, 
0.0355868, 0.0374675, 0.039342, 0.0412104, 0.0430725, 0.0449285, 
0.0467783, 0.048622, 0.0504595, 0.0522909, 0.0541162, 0.0559355, 
0.0577486, 0.0595557, 0.0613567, 0.0631517, 0.0649406, 0.0667235, 
0.0685005, 0.0702714, 0.0720363, 0.0737953, 0.0755483, 0.0772954, 
0.0790366, 0.0807718, 0.0825012, 0.0842246, 0.0859422, 0.0876538, 
0.0893597, 0.0910597, 0.0927539, 0.0944422, 0.0961247, 0.0978015, 
0.0994725, 0.101138, 0.102797, 0.104451, 0.106099, 0.107741, 
0.109378, 0.111009, 0.112634, 0.114253, 0.115867, 0.117475, 0.119078, 
0.120675, 0.122266, 0.123852, 0.125432, 0.127007, 0.128576, 0.13014, 
0.131698, 0.13325, 0.134797, 0.136339, 0.137875, 0.139406, 0.140931, 
0.14245, 0.143965, 0.145474, 0.146977, 0.148475, 0.149968, 0.151455, 
0.152937, 0.154414, 0.155885, 0.157351, 0.158812, 0.160267, 0.161718, 
0.163163, 0.164602, 0.166037, 0.167466, 0.16889, 0.170309, 0.171722, 
0.173131, 0.174534, 0.175932, 0.177325, 0.178713, 0.180096, 0.181474, 
0.182846, 0.184214, 0.185576, 0.186934, 0.188286, 0.189634, 0.190976, 
0.192314, 0.193646, 0.194974, 0.196296, 0.197614, 0.198926, 0.200234, 
0.201537, 0.202835, 0.204128, 0.205416, 0.2067, 0.207978, 0.209252, 
0.210521, 0.211785, 0.213045, 0.214299, 0.215549, 0.216794, 0.218034, 
0.21927, 0.220501, 0.221727, 0.222949, 0.224166, 0.225378, 0.226585, 
0.227788, 0.228987, 0.23018, 0.23137, 0.232554, 0.233734, 0.234909, 
0.23608, 0.237247, 0.238409, 0.239566, 0.240719, 0.241867, 0.243011, 
0.24415, 0.245285, 0.246416, 0.247542, 0.248664, 0.249781, 0.250894, 
0.252003, 0.253107, 0.254207, 0.255303, 0.256394, 0.257481, 0.258564, 
0.259642, 0.260716, 0.261786, 0.262852, 0.263913, 0.26497, 0.266023, 
0.267072, 0.268117, 0.269157, 0.270194, 0.271226, 0.272254, 0.273278, 
0.274298, 0.275313, 0.276325, 0.277333, 0.278336, 0.279336, 0.280331, 
0.281322, 0.28231, 0.283293, 0.284273, 0.285248, 0.28622, 0.287187, 
0.288151, 0.289111, 0.290066, 0.291018, 0.291966, 0.29291, 0.29385, 
0.294787, 0.295719, 0.296648, 0.297573, 0.298494, 0.299411, 0.300325, 
0.301235, 0.302141, 0.303043, 0.303941, 0.304836, 0.305727, 0.306614, 
0.307498, 0.308378, 0.309254, 0.310127, 0.310996, 0.311861, 0.312723, 
0.313581, 0.314436, 0.315287, 0.316134, 0.316978, 0.317818, 0.318655, 
0.319488, 0.320317, 0.321144, 0.321966, 0.322785, 0.323601, 0.324413, 
0.325222, 0.326027, 0.326829, 0.327627, 0.328422, 0.329214, 0.330002, 
0.330787, 0.331568, 0.332346, 0.333121, 0.333893, 0.334661, 0.335425, 
0.336187, 0.336945, 0.3377, 0.338451, 0.3392, 0.339945, 0.340687, 
0.341425, 0.342161, 0.342893, 0.343622, 0.344347, 0.34507, 0.345789, 
0.346506, 0.347219, 0.347929, 0.348635, 0.349339, 0.35004, 0.350737, 
0.351431, 0.352123, 0.352811, 0.353496, 0.354178, 0.354857, 0.355533, 
0.356206, 0.356876, 0.357543, 0.358207, 0.358868, 0.359526, 0.360181, 
0.360833, 0.361482, 0.362128, 0.362772, 0.363412, 0.364049, 0.364684, 
0.365316, 0.365944, 0.36657, 0.367193, 0.367814, 0.368431, 0.369045, 
0.369657, 0.370266, 0.370872, 0.371475, 0.372076, 0.372673, 0.373268, 
0.37386, 0.37445, 0.375036, 0.37562, 0.376201, 0.37678, 0.377356, 
0.377929, 0.378499, 0.379067, 0.379632, 0.380194, 0.380753, 0.381311, 
0.381865, 0.382417, 0.382966, 0.383512, 0.384056, 0.384597, 0.385136, 
0.385672, 0.386206, 0.386737, 0.387265, 0.387791, 0.388314, 0.388835, 
0.389353, 0.389869, 0.390382, 0.390893, 0.391401, 0.391907, 0.39241, 
0.392911, 0.393409, 0.393905, 0.394399, 0.39489, 0.395378, 0.395864, 
0.396348, 0.396829, 0.397308, 0.397785, 0.398259, 0.398731, 0.3992, 
0.399667, 0.400132, 0.400594, 0.401054, 0.401512, 0.401967, 0.40242, 
0.402871, 0.40332, 0.403766, 0.40421, 0.404651, 0.40509, 0.405528, 
0.405962, 0.406395, 0.406825, 0.407253, 0.407679, 0.408103, 0.408524, 
0.408943, 0.40936, 0.409775, 0.410187, 0.410598, 0.411006, 0.411412, 
0.411816, 0.412218, 0.412617, 0.413015, 0.41341, 0.413803, 0.414194, 
0.414583, 0.41497, 0.415355, 0.415737, 0.416118, 0.416496, 0.416873, 
0.417247, 0.417619, 0.41799, 0.418358, 0.418724, 0.419088, 0.41945, 
0.41981, 0.420168, 0.420524, 0.420878, 0.42123, 0.42158, 0.421928, 
0.422274, 0.422618, 0.42296, 0.423301, 0.423639, 0.423975, 0.424309, 
0.424642, 0.424972, 0.425301, 0.425627, 0.425952, 0.426275, 0.426595, 
0.426914, 0.427231, 0.427547, 0.42786, 0.428171, 0.428481, 0.428788, 
0.429094, 0.429398, 0.4297, 0.430001, 0.430299, 0.430596, 0.43089, 
0.431183, 0.431475, 0.431764, 0.432051, 0.432337, 0.432621, 0.432903, 
0.433184, 0.433462, 0.433739, 0.434014, 0.434288, 0.434559, 0.434829, 
0.435097, 0.435363, 0.435628, 0.435891, 0.436152, 0.436411, 0.436669, 
0.436925, 0.437179, 0.437432, 0.437683, 0.437932, 0.438179, 0.438425, 
0.438669, 0.438912, 0.439153, 0.439392, 0.439629, 0.439865, 0.4401, 
0.440332, 0.440563, 0.440792, 0.44102, 0.441246, 0.44147, 0.441693, 
0.441914, 0.442134, 0.442352, 0.442568, 0.442783, 0.442997, 0.443208, 
0.443418, 0.443627, 0.443834, 0.444039, 0.444243, 0.444445, 0.444646, 
0.444845, 0.445043, 0.445239, 0.445433, 0.445626, 0.445818, 0.446008, 
0.446196, 0.446383, 0.446568, 0.446752, 0.446935, 0.447116, 0.447295, 
0.447473, 0.447649, 0.447824, 0.447998, 0.44817, 0.448341, 0.44851, 
0.448677, 0.448843, 0.449008, 0.449171, 0.449333, 0.449494, 0.449653, 
0.44981, 0.449966, 0.450121, 0.450274, 0.450426, 0.450576, 0.450725, 
0.450873, 0.451019, 0.451164, 0.451307, 0.451449, 0.45159, 0.451729, 
0.451867, 0.452004, 0.452139, 0.452272, 0.452405, 0.452536, 0.452665, 
0.452794, 0.452921, 0.453046, 0.453171, 0.453293, 0.453415, 0.453535, 
0.453654, 0.453772, 0.453888, 0.454003, 0.454116, 0.454229, 0.45434, 
0.454449, 0.454558, 0.454665, 0.45477, 0.454875, 0.454978, 0.45508, 
0.455181, 0.45528, 0.455378, 0.455475, 0.45557, 0.455664, 0.455757, 
0.455849, 0.455939, 0.456029, 0.456116, 0.456203, 0.456289, 0.456373, 
0.456456, 0.456537, 0.456618, 0.456697, 0.456775, 0.456852, 0.456927, 
0.457001, 0.457075, 0.457146, 0.457217, 0.457287, 0.457355, 0.457422, 
0.457488, 0.457552, 0.457616, 0.457678, 0.457739, 0.457799, 0.457857, 
0.457915, 0.457971, 0.458026, 0.45808, 0.458133, 0.458185, 0.458235, 
0.458285, 0.458333, 0.45838, 0.458425, 0.45847, 0.458514, 0.458556, 
0.458597, 0.458637, 0.458676, 0.458714, 0.458751, 0.458786, 0.458821, 
0.458854, 0.458886, 0.458917, 0.458947, 0.458976, 0.459003, 0.45903, 
0.459055, 0.45908, 0.459103, 0.459125, 0.459146, 0.459166, 0.459185, 
0.459203, 0.459219, 0.459235, 0.459249, 0.459263, 0.459275, 0.459286, 
0.459296, 0.459305, 0.459313, 0.45932, 0.459326, 0.459331, 0.459335, 
0.459337, 0.459339, 0.459339, 0.459339};


vector<double> x1 ={0.00265231, 0.00464052, 0.00662231, 0.00859768, 0.0105667, 
0.0125292, 0.0144855, 0.0164353, 0.0183788, 0.020316, 0.0222469, 
0.0241714, 0.0260897, 0.0280016, 0.0299074, 0.0318068, 0.0337001, 
0.0355871, 0.0374679, 0.0393426, 0.041211, 0.0430733, 0.0449295, 
0.0467795, 0.0486235, 0.0504613, 0.0522931, 0.0541188, 0.0559384, 
0.0577521, 0.0595597, 0.0613614, 0.0631571, 0.0649468, 0.0667307, 
0.0685086, 0.0702806, 0.0720467, 0.073807, 0.0755614, 0.0773101, 
0.0790529, 0.08079, 0.0825213, 0.0842469, 0.0859668, 0.0876809, 
0.0893894, 0.0910923, 0.0927895, 0.0944812, 0.0961672, 0.0978477, 
0.0995227, 0.101192, 0.102856, 0.104514, 0.106167, 0.107815, 
0.109457, 0.111094, 0.112725, 0.114352, 0.115972, 0.117588, 0.119198, 
0.120803, 0.122403, 0.123997, 0.125587, 0.127171, 0.12875, 0.130324, 
0.131893, 0.133457, 0.135016, 0.13657, 0.138119, 0.139662, 0.141201, 
0.142735, 0.144265, 0.145789, 0.147308, 0.148823, 0.150333, 0.151838, 
0.153338, 0.154834, 0.156325, 0.157811, 0.159292, 0.160769, 0.162242, 
0.163709, 0.165172, 0.166631, 0.168085, 0.169535, 0.17098, 0.172421, 
0.173857, 0.175289, 0.176716, 0.17814, 0.179558, 0.180973, 0.182383, 
0.183789, 0.185191, 0.186588, 0.187982, 0.189371, 0.190756, 0.192137, 
0.193513, 0.194886, 0.196255, 0.197619, 0.19898, 0.200336, 0.201689, 
0.203037, 0.204382, 0.205723, 0.20706, 0.208393, 0.209722, 0.211047, 
0.212369, 0.213686, 0.215, 0.21631, 0.217617, 0.218919, 0.220218, 
0.221514, 0.222805, 0.224093, 0.225378, 0.226659, 0.227936, 0.22921, 
0.23048, 0.231746, 0.233009, 0.234269, 0.235525, 0.236778, 0.238027, 
0.239273, 0.240516, 0.241755, 0.242991, 0.244223, 0.245452, 0.246678, 
0.2479, 0.24912, 0.250336, 0.251548, 0.252758, 0.253964, 0.255167, 
0.256367, 0.257564, 0.258758, 0.259949, 0.261136, 0.26232, 0.263502, 
0.26468, 0.265855, 0.267027, 0.268197, 0.269363, 0.270526, 0.271686, 
0.272843, 0.273998, 0.275149, 0.276298, 0.277443, 0.278586, 0.279726, 
0.280863, 0.281997, 0.283129, 0.284257, 0.285383, 0.286506, 0.287627, 
0.288744, 0.289859, 0.290971, 0.29208, 0.293187, 0.294291, 0.295392, 
0.296491, 0.297587, 0.29868, 0.299771, 0.300859, 0.301945, 0.303028, 
0.304108, 0.305186, 0.306262, 0.307334, 0.308405, 0.309472, 0.310538, 
0.311601, 0.312661, 0.313719, 0.314774, 0.315827, 0.316878, 0.317926, 
0.318972, 0.320015, 0.321056, 0.322094, 0.323131, 0.324165, 0.325196, 
0.326225, 0.327252, 0.328277, 0.329299, 0.330319, 0.331337, 0.332352, 
0.333366, 0.334377, 0.335385, 0.336392, 0.337396, 0.338398, 0.339398, 
0.340396, 0.341391, 0.342385, 0.343376, 0.344365, 0.345352, 0.346337, 
0.347319, 0.3483, 0.349278, 0.350255, 0.351229, 0.352201, 0.353171, 
0.354139, 0.355105, 0.356069, 0.357031, 0.357991, 0.358949, 0.359905, 
0.360859, 0.361811, 0.362761, 0.363708, 0.364654, 0.365598, 0.366541, 
0.367481, 0.368419, 0.369355, 0.370289, 0.371222, 0.372152, 0.373081, 
0.374008, 0.374932, 0.375855, 0.376777, 0.377696, 0.378613, 0.379529, 
0.380442, 0.381354, 0.382264, 0.383173, 0.384079, 0.384984, 0.385886, 
0.386787, 0.387687, 0.388584, 0.38948, 0.390374, 0.391266, 0.392156, 
0.393045, 0.393932, 0.394817, 0.395701, 0.396582, 0.397462, 0.398341, 
0.399217, 0.400092, 0.400966, 0.401837, 0.402707, 0.403575, 0.404442, 
0.405307, 0.40617, 0.407032, 0.407892, 0.40875, 0.409607, 0.410462, 
0.411315, 0.412167, 0.413018, 0.413866, 0.414713, 0.415559, 0.416403, 
0.417245, 0.418086, 0.418925, 0.419763, 0.420599, 0.421433, 0.422266, 
0.423098, 0.423928, 0.424756, 0.425583, 0.426408, 0.427232, 0.428054, 
0.428875, 0.429694, 0.430512, 0.431328, 0.432143, 0.432956, 0.433768, 
0.434579, 0.435388, 0.436195, 0.437001, 0.437806, 0.438609, 0.43941, 
0.440211, 0.441009, 0.441807, 0.442603, 0.443397, 0.44419, 0.444982, 
0.445772, 0.446561, 0.447349, 0.448135, 0.448919, 0.449703, 0.450484, 
0.451265, 0.452044, 0.452822, 0.453598, 0.454374, 0.455147, 0.45592, 
0.456691, 0.45746, 0.458229, 0.458996, 0.459761, 0.460526, 0.461289, 
0.46205, 0.462811, 0.46357, 0.464327, 0.465084, 0.465839, 0.466593, 
0.467345, 0.468097, 0.468847, 0.469595, 0.470343, 0.471089, 0.471834, 
0.472577, 0.47332, 0.474061, 0.474801, 0.475539, 0.476276, 0.477012, 
0.477747, 0.478481, 0.479213, 0.479944, 0.480674, 0.481403, 0.48213, 
0.482856, 0.483581, 0.484305, 0.485028, 0.485749, 0.486469, 0.487188, 
0.487906, 0.488623, 0.489338, 0.490052, 0.490765, 0.491477, 0.492188, 
0.492897, 0.493606, 0.494313, 0.495019, 0.495724, 0.496427, 0.49713, 
0.497831, 0.498531, 0.49923, 0.499928, 0.500625, 0.501321, 0.502015, 
0.502709, 0.503401, 0.504092, 0.504782, 0.505471, 0.506159, 0.506845, 
0.507531, 0.508215, 0.508899, 0.509581, 0.510262, 0.510942, 0.511621, 
0.512299, 0.512975, 0.513651, 0.514326, 0.514999, 0.515672, 0.516343, 
0.517013, 0.517683, 0.518351, 0.519018, 0.519684, 0.520349, 0.521013, 
0.521676, 0.522337, 0.522998, 0.523658, 0.524317, 0.524974, 0.525631, 
0.526286, 0.526941, 0.527594, 0.528247, 0.528898, 0.529549, 0.530198, 
0.530847, 0.531494, 0.53214, 0.532786, 0.53343, 0.534073, 0.534716, 
0.535357, 0.535997, 0.536637, 0.537275, 0.537912, 0.538549, 0.539184, 
0.539818, 0.540452, 0.541084, 0.541716, 0.542346, 0.542976, 0.543604, 
0.544232, 0.544859, 0.545484, 0.546109, 0.546733, 0.547355, 0.547977, 
0.548598, 0.549218, 0.549837, 0.550455, 0.551072, 0.551688, 0.552303, 
0.552918, 0.553531, 0.554143, 0.554755, 0.555365, 0.555975, 0.556584, 
0.557192, 0.557798, 0.558404, 0.559009, 0.559613, 0.560217, 0.560819, 
0.56142, 0.562021, 0.56262, 0.563219, 0.563817, 0.564414, 0.56501, 
0.565605, 0.566199, 0.566792, 0.567385, 0.567976, 0.568567, 0.569157, 
0.569746, 0.570334, 0.570921, 0.571507, 0.572093, 0.572677, 0.573261, 
0.573844, 0.574426, 0.575007, 0.575587, 0.576166, 0.576745, 0.577322, 
0.577899, 0.578475, 0.57905, 0.579625, 0.580198, 0.580771, 0.581342, 
0.581913, 0.582483, 0.583053, 0.583621, 0.584189, 0.584755, 0.585321, 
0.585886, 0.586451, 0.587014, 0.587577, 0.588138, 0.588699, 0.58926, 
0.589819, 0.590378, 0.590935, 0.591492, 0.592048, 0.592604, 0.593158, 
0.593712, 0.594265, 0.594817, 0.595368, 0.595919, 0.596469, 0.597018, 
0.597566, 0.598113, 0.59866, 0.599206, 0.599751, 0.600295, 0.600838, 
0.601381, 0.601923, 0.602464, 0.603005, 0.603544, 0.604083, 0.604621, 
0.605158, 0.605695, 0.606231, 0.606766, 0.6073, 0.607834, 0.608366, 
0.608898, 0.60943, 0.60996, 0.61049, 0.611019, 0.611547, 0.612075, 
0.612601, 0.613127, 0.613653, 0.614177, 0.614701, 0.615224, 0.615746, 
0.616268, 0.616789, 0.617309, 0.617828, 0.618347, 0.618865, 0.619382, 
0.619899, 0.620415, 0.62093, 0.621444, 0.621958, 0.622471, 0.622983, 
0.623495, 0.624005, 0.624516, 0.625025, 0.625534, 0.626042, 0.626549, 
0.627056, 0.627561, 0.628067, 0.628571, 0.629075, 0.629578, 0.63008, 
0.630582, 0.631083, 0.631584, 0.632083};

vector<double> x0 = {0.9596, 0.9591, 0.9586, 0.9581, 0.9576, 0.9571, 0.9566, 0.9561, 
0.9556, 0.9551, 0.9546, 0.9541, 0.9536, 0.9531, 0.9526, 0.9521, 
0.9516, 0.9511, 0.9506, 0.9501, 0.9496, 0.9491, 0.9486, 0.9481, 
0.9476, 0.9471, 0.9466, 0.9461, 0.9456, 0.9451, 0.9446, 0.9441, 
0.9436, 0.9431, 0.9426, 0.9421, 0.9416, 0.9411, 0.9406, 0.9401, 
0.9396, 0.9391, 0.9386, 0.9381, 0.9376, 0.9371, 0.9366, 0.9361, 
0.9356, 0.9351, 0.9346, 0.9341, 0.9336, 0.9331, 0.9326, 0.9321, 
0.9316, 0.9311, 0.9306, 0.9301, 0.9296, 0.9291, 0.9286, 0.9281, 
0.9276, 0.9271, 0.9266, 0.9261, 0.9256, 0.9251, 0.9246, 0.9241, 
0.9236, 0.9231, 0.9226, 0.9221, 0.9216, 0.9211, 0.9206, 0.9201, 
0.9196, 0.9191, 0.9186, 0.9181, 0.9176, 0.9171, 0.9166, 0.9161, 
0.9156, 0.9151, 0.9146, 0.9141, 0.9136, 0.9131, 0.9126, 0.9121, 
0.9116, 0.9111, 0.9106, 0.9101, 0.9096, 0.9091, 0.9086, 0.9081, 
0.9076, 0.9071, 0.9066, 0.9061, 0.9056, 0.9051, 0.9046, 0.9041, 
0.9036, 0.9031, 0.9026, 0.9021, 0.9016, 0.9011, 0.9006, 0.9001, 
0.8996, 0.8991, 0.8986, 0.8981, 0.8976, 0.8971, 0.8966, 0.8961, 
0.8956, 0.8951, 0.8946, 0.8941, 0.8936, 0.8931, 0.8926, 0.8921, 
0.8916, 0.8911, 0.8906, 0.8901, 0.8896, 0.8891, 0.8886, 0.8881, 
0.8876, 0.8871, 0.8866, 0.8861, 0.8856, 0.8851, 0.8846, 0.8841, 
0.8836, 0.8831, 0.8826, 0.8821, 0.8816, 0.8811, 0.8806, 0.8801, 
0.8796, 0.8791, 0.8786, 0.8781, 0.8776, 0.8771, 0.8766, 0.8761, 
0.8756, 0.8751, 0.8746, 0.8741, 0.8736, 0.8731, 0.8726, 0.8721, 
0.8716, 0.8711, 0.8706, 0.8701, 0.8696, 0.8691, 0.8686, 0.8681, 
0.8676, 0.8671, 0.8666, 0.8661, 0.8656, 0.8651, 0.8646, 0.8641, 
0.8636, 0.8631, 0.8626, 0.8621, 0.8616, 0.8611, 0.8606, 0.8601, 
0.8596, 0.8591, 0.8586, 0.8581, 0.8576, 0.8571, 0.8566, 0.8561, 
0.8556, 0.8551, 0.8546, 0.8541, 0.8536, 0.8531, 0.8526, 0.8521, 
0.8516, 0.8511, 0.8506, 0.8501, 0.8496, 0.8491, 0.8486, 0.8481, 
0.8476, 0.8471, 0.8466, 0.8461, 0.8456, 0.8451, 0.8446, 0.8441, 
0.8436, 0.8431, 0.8426, 0.8421, 0.8416, 0.8411, 0.8406, 0.8401, 
0.8396, 0.8391, 0.8386, 0.8381, 0.8376, 0.8371, 0.8366, 0.8361, 
0.8356, 0.8351, 0.8346, 0.8341, 0.8336, 0.8331, 0.8326, 0.8321, 
0.8316, 0.8311, 0.8306, 0.8301, 0.8296, 0.8291, 0.8286, 0.8281, 
0.8276, 0.8271, 0.8266, 0.8261, 0.8256, 0.8251, 0.8246, 0.8241, 
0.8236, 0.8231, 0.8226, 0.8221, 0.8216, 0.8211, 0.8206, 0.8201, 
0.8196, 0.8191, 0.8186, 0.8181, 0.8176, 0.8171, 0.8166, 0.8161, 
0.8156, 0.8151, 0.8146, 0.8141, 0.8136, 0.8131, 0.8126, 0.8121, 
0.8116, 0.8111, 0.8106, 0.8101, 0.8096, 0.8091, 0.8086, 0.8081, 
0.8076, 0.8071, 0.8066, 0.8061, 0.8056, 0.8051, 0.8046, 0.8041, 
0.8036, 0.8031, 0.8026, 0.8021, 0.8016, 0.8011, 0.8006, 0.8001, 
0.7996, 0.7991, 0.7986, 0.7981, 0.7976, 0.7971, 0.7966, 0.7961, 
0.7956, 0.7951, 0.7946, 0.7941, 0.7936, 0.7931, 0.7926, 0.7921, 
0.7916, 0.7911, 0.7906, 0.7901, 0.7896, 0.7891, 0.7886, 0.7881, 
0.7876, 0.7871, 0.7866, 0.7861, 0.7856, 0.7851, 0.7846, 0.7841, 
0.7836, 0.7831, 0.7826, 0.7821, 0.7816, 0.7811, 0.7806, 0.7801, 
0.7796, 0.7791, 0.7786, 0.7781, 0.7776, 0.7771, 0.7766, 0.7761, 
0.7756, 0.7751, 0.7746, 0.7741, 0.7736, 0.7731, 0.7726, 0.7721, 
0.7716, 0.7711, 0.7706, 0.7701, 0.7696, 0.7691, 0.7686, 0.7681, 
0.7676, 0.7671, 0.7666, 0.7661, 0.7656, 0.7651, 0.7646, 0.7641, 
0.7636, 0.7631, 0.7626, 0.7621, 0.7616, 0.7611, 0.7606, 0.7601, 
0.7596, 0.7591, 0.7586, 0.7581, 0.7576, 0.7571, 0.7566, 0.7561, 
0.7556, 0.7551, 0.7546, 0.7541, 0.7536, 0.7531, 0.7526, 0.7521, 
0.7516, 0.7511, 0.7506, 0.7501, 0.7496, 0.7491, 0.7486, 0.7481, 
0.7476, 0.7471, 0.7466, 0.7461, 0.7456, 0.7451, 0.7446, 0.7441, 
0.7436, 0.7431, 0.7426, 0.7421, 0.7416, 0.7411, 0.7406, 0.7401, 
0.7396, 0.7391, 0.7386, 0.7381, 0.7376, 0.7371, 0.7366, 0.7361, 
0.7356, 0.7351, 0.7346, 0.7341, 0.7336, 0.7331, 0.7326, 0.7321, 
0.7316, 0.7311, 0.7306, 0.7301, 0.7296, 0.7291, 0.7286, 0.7281, 
0.7276, 0.7271, 0.7266, 0.7261, 0.7256, 0.7251, 0.7246, 0.7241, 
0.7236, 0.7231, 0.7226, 0.7221, 0.7216, 0.7211, 0.7206, 0.7201, 
0.7196, 0.7191, 0.7186, 0.7181, 0.7176, 0.7171, 0.7166, 0.7161, 
0.7156, 0.7151, 0.7146, 0.7141, 0.7136, 0.7131, 0.7126, 0.7121, 
0.7116, 0.7111, 0.7106, 0.7101, 0.7096, 0.7091, 0.7086, 0.7081, 
0.7076, 0.7071, 0.7066, 0.7061, 0.7056, 0.7051, 0.7046, 0.7041, 
0.7036, 0.7031, 0.7026, 0.7021, 0.7016, 0.7011, 0.7006, 0.7001, 
0.6996, 0.6991, 0.6986, 0.6981, 0.6976, 0.6971, 0.6966, 0.6961, 
0.6956, 0.6951, 0.6946, 0.6941, 0.6936, 0.6931, 0.6926, 0.6921, 
0.6916, 0.6911, 0.6906, 0.6901, 0.6896, 0.6891, 0.6886, 0.6881, 
0.6876, 0.6871, 0.6866, 0.6861, 0.6856, 0.6851, 0.6846, 0.6841, 
0.6836, 0.6831, 0.6826, 0.6821, 0.6816, 0.6811, 0.6806, 0.6801, 
0.6796, 0.6791, 0.6786, 0.6781, 0.6776, 0.6771, 0.6766, 0.6761, 
0.6756, 0.6751, 0.6746, 0.6741, 0.6736, 0.6731, 0.6726, 0.6721, 
0.6716, 0.6711, 0.6706, 0.6701, 0.6696, 0.6691, 0.6686, 0.6681, 
0.6676, 0.6671, 0.6666, 0.6661, 0.6656, 0.6651, 0.6646, 0.6641, 
0.6636, 0.6631, 0.6626, 0.6621, 0.6616, 0.6611, 0.6606, 0.6601, 
0.6596, 0.6591, 0.6586, 0.6581, 0.6576, 0.6571, 0.6566, 0.6561, 
0.6556, 0.6551, 0.6546, 0.6541, 0.6536, 0.6531, 0.6526, 0.6521, 
0.6516, 0.6511, 0.6506, 0.6501, 0.6496, 0.6491, 0.6486, 0.6481, 
0.6476, 0.6471, 0.6466, 0.6461, 0.6456, 0.6451, 0.6446, 0.6441, 
0.6436, 0.6431, 0.6426, 0.6421, 0.6416, 0.6411, 0.6406, 0.6401, 
0.6396, 0.6391, 0.6386, 0.6381, 0.6376, 0.6371, 0.6366, 0.6361, 
0.6356, 0.6351, 0.6346, 0.6341, 0.6336, 0.6331, 0.6326, 0.6321, 
0.6316, 0.6311};

const static int kSimpsonStep = 99;

  DBound(const Membrane& m);

  double operator()(int i) const;
  double B1(int i) const;
  double B2(int i) const;

  double Rho(int i) const;
  double RhodRho(int i) const;
  double S(int i) const;
  double SdS(int i) const;
  double Alpha(int i) const;
  double AlphadAlpha(int i) const;
  void PrintX0X1(int q) const;

  double SigmaE(int i) const;
  double H(int i) const;

  Membrane m_;
};

#endif //BOUND_H_

//
// Traffic Sign Recognizer
// 算法
// 作者：SunnyCase
// 创建日期：2016-02-26
#include "common.h"
#include <amp.h>
#include <amp_graphics.h>

DEFINE_NS_TSR_PRCSR

concurrency::graphics::unorm Grayscale(const concurrency::graphics::unorm_4& color) restrict(cpu, amp);
concurrency::graphics::unorm Threshold(concurrency::graphics::unorm color, float threshold) restrict(cpu, amp);
concurrency::graphics::unorm SusanTest(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index, unsigned int radius, float threshold) restrict(amp);
float CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index) restrict(amp);
bool FitEllipse(concurrency::index<2> p1, concurrency::index<2> p2, float p1Tan, float p2Tan, const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float, 2>& tangentView, float& v) restrict(amp);
concurrency::graphics::float_2 coord(concurrency::graphics::float_2 point, const concurrency::extent<2>& extent) restrict(cpu, amp);
bool Solve(float(&mat)[5][6]) restrict(cpu, amp);

template<typename T, uint32_t N>
void swap(T(&v1)[N], T(&v2)[N]) restrict(cpu, amp)
{
	T tmp[N];
	for (uint32_t i = 0; i < N; i++)
		tmp[i] = v1[i];
	for (uint32_t i = 0; i < N; i++)
		v1[i] = v2[i];
	for (uint32_t i = 0; i < N; i++)
		v2[i] = tmp[i];
}

END_NS_TSR_PRCSR
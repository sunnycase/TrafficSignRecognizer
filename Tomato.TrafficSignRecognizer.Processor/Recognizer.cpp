//
// Traffic Sign Recognizer
// Recognizer
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "pch.h"
#include "Recognizer.h"
#include <amp_math.h>
#include "algorithms.h"
#include "amprng/amp_sobol_rng.h"

using namespace NS_TSR_PRCSR;
using namespace Platform;
using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Foundation;
using namespace Windows::Graphics::Imaging;
using namespace WRL;

static const uint32_t zero = 0;

Recognizer::Recognizer(unsigned int width, unsigned int height)
	:_targetImageExtent(height, width),
	_acc_view(concurrency::direct3d::create_accelerator_view(accelerator(), true)),
	_targetImage(_targetImageExtent, 16U, _acc_view), _outputTex(_targetImageExtent, _acc_view),
	_edgeTex(_targetImageExtent, 8, _acc_view), _tangentTex(_targetImageExtent, 32, _acc_view),
	_redTex(_targetImageExtent, 8, _acc_view)
{

}

Windows::Foundation::IAsyncAction ^ Recognizer::SetTarget(BitmapFrame ^ frame)
{
	return create_async([=]
	{
		return PrepareTargetImage(frame);
	});
}

IAsyncOperation<bool>^ Recognizer::Recognize(Windows::Storage::Streams::IRandomAccessStream^ outputStream)
{
	return create_async([=]() -> task<bool>
	{
		return FindContours().then([=]
		{
			return create_task(BitmapEncoder::CreateAsync(BitmapEncoder::JpegEncoderId, outputStream))
				.then([=](BitmapEncoder^ encoder)
			{
				auto data = std::make_shared<std::vector<uint>>(_targetImageExtent.size());
				return copy_async(_outputTex, data->begin()).to_task().then([=]
				{
					encoder->SetPixelData(BitmapPixelFormat::Bgra8, BitmapAlphaMode::Premultiplied, _targetImageExtent[1], _targetImageExtent[0],
						96., 96., ArrayReference<byte>((byte*)data->data(), data->size() * 4));
					return create_task(encoder->FlushAsync());
				});
			});
		}).then([]
		{
			return false;
		});
	});
}

concurrency::task<void> Recognizer::PrepareTargetImage(BitmapFrame ^ frame)
{
	auto transform = ref new BitmapTransform();
	return create_task(frame->GetPixelDataAsync(BitmapPixelFormat::Bgra8, BitmapAlphaMode::Premultiplied, transform,
		ExifOrientationMode::RespectExifOrientation, ColorManagementMode::ColorManageToSRgb))
		.then([=](PixelDataProvider^ dataProvider)
	{
		auto data = dataProvider->DetachPixelData();
		array<uint, 2> src(_targetImageExtent, (uint*)data->begin(), (uint*)data->end(), _acc_view);
		texture_view<unorm_4, 2> targetImage(_targetImage);
		parallel_for_each(_acc_view, src.extent, [&src, targetImage](index<2> index) restrict(amp)
		{
			auto pixel = src[index];

			auto r = ((pixel & 0xFF0000) >> 16) / 255.f;
			auto g = ((pixel & 0xFF00) >> 8) / 255.f;
			auto b = (pixel & 0xFF) / 255.f;
			targetImage.set(index, unorm_4(r, g, b, 1.f));
		});
	});
}

concurrency::task<void> Recognizer::FindContours()
{
	//FindEdgesAndTangent();
	//texture_view<const unorm, 2> edges(_edgeTex);
	//texture_view<const unorm_4, 2> input(_targetImage);

	AbsorbRedTexels();
	FindEdgesAndTangent();

	texture_view<const unorm, 2> redTex(_redTex);
	texture_view<const unorm, 2> edges(_edgeTex);
	texture_view<const float, 2> tangents(_tangentTex);
	array_view<uint, 2> outputTex(_outputTex);

	return FindEllipses().then([=]
	{
		parallel_for_each(_acc_view, _targetImageExtent, [outputTex, tangents, edges](index<2> index) restrict(amp)
		{
			//const auto gray = fast_math::sin(fast_math::atan(tangents[index]));
			//if (gray < 0)
			//{
			//	auto lastGray = uint(-gray * 255);
			//	outputTex[index] = uint(0xFF0000FF | (lastGray << 16));
			//}
			//else if (gray > 0)
			//{
			//	auto lastGray = uint(gray * 255);
			//	outputTex[index] = uint(0xFF0000FF | (lastGray << 8));
			//}
			//const auto r = input[index].r * 255.f;
			//const auto g = input[index].g * 255.f;
			//const auto b = input[index].b * 255.f;
			//if (((g - r) / r >= 1.f) && (g / b >= 1.f) && (g - b <= 40))
			//{
			//	outputTex[index] = uint(0xFF000000);
			//}
			//else
			//{
			//	if (r - g > 40 && r - b > 50)
			//	{
			//		outputTex[index] = uint(0xFFFF0000);
			//	}
			//	else
			//		outputTex[index] = uint(0xFF000000);
			//}
		});
	});
}

void Recognizer::FindEdgesAndTangent()
{
	uint radius = 3;
	texture_view<unorm, 2> edgeView(_edgeTex);
	texture_view<float, 2> tangentView(_tangentTex);
	texture_view<const unorm, 2> redTex(_redTex);
	texture_view<const unorm_4, 2> inputTex(_targetImage);
	parallel_for_each(_acc_view, redTex.extent, [redTex, edgeView, inputTex, tangentView, radius](index<2> index) restrict(amp)
	{
		auto edge = SusanTest(redTex, index, radius, 0.1f);
		edgeView.set(index, edge);

		if (edge > 0.f)
			tangentView.set(index, CalculateTangent(redTex, index));
	});
}

void Recognizer::AbsorbRedTexels()
{
	texture_view<const unorm_4, 2> input(_targetImage);
	texture_view<unorm, 2> redTex(_redTex);

	parallel_for_each(_acc_view, _targetImageExtent, [input, redTex](index<2> index) restrict(amp)
	{
		const float_2 p(index[1] + 0.5f, index[0] + 0.5f);
		const float_2 points[] = {
			p + float_2(-1, -1),
			p + float_2(+0, -1),
			p + float_2(+1, -1),

			p + float_2(-1, +0),
			p + float_2(+0, +0),
			p + float_2(+1, +0),

			p + float_2(-1, +1),
			p + float_2(+0, +1),
			p + float_2(+1, +1)
		};
		uint32_t hasCount = 0;
		for (auto&& point : points)
		{
			if (IsRed(input.sample<filter_point>(coord(point, input.extent))))
				if (++hasCount >= 5)
				{
					redTex.set(index, unorm(1.f));
					break;
				}
		}
	});

	texture_view<const unorm, 2> redReader(_redTex);
	texture<unorm, 2> tmpTex(_targetImageExtent, 8, _acc_view);
	texture_view<unorm, 2> tmpWriter(tmpTex);
	parallel_for_each(_acc_view, _targetImageExtent, [tmpWriter, redReader](index<2> index) restrict(amp)
	{
		const float_2 p(index[1] + 0.5f, index[0] + 0.5f);
		if (redReader.sample<filter_point>(coord(p, redReader.extent)) == 0.f)
		{
			if (redReader.sample<filter_point>(coord(p + float_2(-1, +0), redReader.extent)) != 0.f ||
				redReader.sample<filter_point>(coord(p + float_2(+1, +0), redReader.extent)) != 0.f)
				tmpWriter.set(index, unorm(1.f));
			else if (redReader.sample<filter_point>(coord(p + float_2(+0, -1), redReader.extent)) != 0.f ||
				redReader.sample<filter_point>(coord(p + float_2(+0, +1), redReader.extent)) != 0.f)
				tmpWriter.set(index, unorm(1.f));
		}
		else
			tmpWriter.set(index, unorm(1.f));
	});

	copy(tmpTex, _redTex);
}

task<void> Recognizer::FindEllipses()
{
	float height = _targetImageExtent[0];

	array<uint, 1> edgePointsCount(1, &zero, _acc_view, access_type_read);
	array<index<2>, 2> edgePositions(_targetImageExtent, _acc_view);
	texture_view<const unorm, 2> edgeView(_edgeTex);

	parallel_for_each(_acc_view, edgeView.extent, [&edgePointsCount, &edgePositions, edgeView](index<2> index) restrict(amp)
	{
		if (edgeView[index] != 0.f)
		{
			auto idx = atomic_fetch_add(&edgePointsCount(0), 1);
			auto x = idx % edgeView.extent[1];
			auto y = idx / edgeView.extent[1];
			edgePositions(y, x) = index;
		}
	});

	const auto edgeSize = extent<1>(edgePointsCount(0));
	const auto maxEllipses = edgeSize * 100;

	typedef sobol_rng<5>::sobol_number<float> sobol_float_number;
	sobol_rng_collection<sobol_rng<5>, 1> sc_rng(_acc_view, maxEllipses, 1);
	concurrency::graphics::texture_view<const float, 2> tangentView(_tangentTex);
	array<uint32_t, 1> fitsCount(1, &zero, _acc_view, access_type_read);
	array<EllipseParam, 1> ellipses(maxEllipses, _acc_view);

	parallel_for_each(_acc_view, maxEllipses, [=, &edgePositions, &fitsCount, &ellipses](index<1> idx) restrict(amp)
	{
		const auto id1 = idx[0];
		const auto x1 = id1 % edgeView.extent[1];
		const auto y1 = id1 / edgeView.extent[1];
		const auto p1Index = edgePositions(y1, x1);

		// Get the sobol RNG 
		auto rng = sc_rng[idx];
		// Skip ahead to the right position
		rng.skip(sc_rng.direction_numbers(), idx[0]);

		float_2 points[5];
		for (uint32_t i = 0; i < 5; i++)
		{
			const auto id2 = int(rng.get_single(i + 1) * edgeSize[0]);
			const auto x2 = id2 % edgeView.extent[1];
			const auto y2 = id2 / edgeView.extent[1];
			const auto p2Index = edgePositions(y2, x2);
			points[i] = float_2(p2Index[1], edgeView.extent[0] - p2Index[0]);
		}

		FitEllipse(points, edgeView.extent[1], edgeView.extent[0], fitsCount, ellipses);
	});
	const auto ellipseCount = fitsCount(0);

	parallel_for_each(_acc_view, edgeSize, [=, &fitsCount, &ellipses, &edgePositions](index<1> index)restrict(amp)
	{
		const auto id1 = index[0];
		const auto x1 = id1 % edgeView.extent[1];
		const auto y1 = id1 / edgeView.extent[1];
		const auto p1 = edgePositions(y1, x1);

		const auto ellipseCount = fitsCount(0);
		for (uint32_t i = 0; i < ellipseCount; i++)
		{
			auto& ellipse = ellipses(i);
			const float_2 point(p1[1], height - p1[0]);
			if (OnEllipse(ellipse, point, 0.8f))
				atomic_fetch_add(&ellipse.rank, 1);
		}
	});

	if (ellipseCount)
	{
		std::vector<EllipseParam> ellipsesSort(ellipseCount);
		copy(ellipses.section(extent<1>(ellipseCount)), ellipsesSort.begin());
		std::sort(ellipsesSort.begin(), ellipsesSort.end(), [](const EllipseParam& left, const EllipseParam& right)
		{
			return left.rank > right.rank;
		});

		_ellipsesFit.clear();
		const float lambda = 0.6f;
		for (auto&& it : ellipsesSort)
		{
			auto min = lambda * 3.14f * (1.5f * (it.a + it.b) - fast_math::sqrt(it.a * it.b));
			if (it.a > 10.f && it.b > 10.f && it.rank >= min)
				_ellipsesFit.emplace_back(it);
		}

		FillCircleSignTargetsSource();
		CalculateZernike();

		//for (auto&& el : _ellipsesFit)
		//{
		//	texture_view<const unorm_4, 2> inputTex(_targetImage);
		//	array_view<uint, 2> outputTex(_outputTex);
		//	array<uint32_t, 1> graySum(1, &zero, _acc_view);
		//	parallel_for_each(_acc_view, _targetImageExtent, [=, &graySum](index<2> index) restrict(amp)
		//	{
		//		const float_2 pt(index[1], height - index[0]);
		//		if (InEllipse(el, pt))
		//		{
		//			atomic_fetch_add(&graySum[0], Grayscale(inputTex[index]) * 255);
		//		}
		//	});
		//	parallel_for_each(_acc_view, _targetImageExtent, [=, &graySum](index<2> index) restrict(amp)
		//	{
		//		const float threhold = graySum[0] / float(el.area) / 255.f;
		//		const float_2 pt(index[1], height - index[0]);
		//		if (InEllipse(el, pt))
		//		{
		//			const uint r = inputTex[index].r * 255;
		//			const uint g = inputTex[index].g * 255;
		//			const uint b = inputTex[index].b * 255;

		//			const float x = (pt.x - el.x);
		//			const float y = (pt.y - el.y);
		//			float x2 = x * fast_math::cos(-el.theta) - y * fast_math::sin(-el.theta);
		//			float y2 = x * fast_math::sin(-el.theta) + y * fast_math::cos(-el.theta);
		//			x2 = x2 / el.a * 50.f;
		//			y2 = y2 / el.b * 50.f;
		//			if (Grayscale(inputTex[index]) < threhold)
		//				outputTex(index) = uint(0xFFFFFFFF);
		//			//outputTex(index) = uint(0xFF000000 | (r << 16) | (g << 8) | b);
		//		}
		//	});
		//}
	}

	return task_from_result();
}

void Recognizer::FillCircleSignTargetsSource()
{
	float height = _targetImageExtent[0];
	texture_view<const unorm_4, 2> inputTex(_targetImage);

	_circleSignTargetsSource.clear();
	for (auto&& el : _ellipsesFit)
	{
		array<UnitCirclePoint, 1> points(el.area, _acc_view);
		array<uint32_t, 1> pointsCount(1, &zero, _acc_view, access_type_read);
		array<uint32_t, 1> graySum(1, &zero, _acc_view);
		parallel_for_each(_acc_view, _targetImageExtent, [=, &graySum](index<2> index) restrict(amp)
		{
			const float_2 pt(index[1], height - index[0]);
			if (InEllipse(el, pt))
				atomic_fetch_add(&graySum[0], Grayscale(inputTex[index]) * 255);
		});

		parallel_for_each(_acc_view, _targetImageExtent, [=, &points, &pointsCount, &graySum](index<2> index) restrict(amp)
		{
			const float threhold = graySum[0] / float(el.area) / 255.f / 1.5f;
			const float_2 pt(index[1], height - index[0]);
			if (InEllipse(el, pt))
			{
				const auto sum = float(atomic_fetch_add(&graySum[0], Grayscale(inputTex[index]) * 255));
				if (Grayscale(inputTex[index]) < threhold)
				{
					const float_2 coord1(pt.x - el.x, pt.y - el.y);
					float_2 coord2(coord1.x * fast_math::cos(-el.theta) - coord1.y * fast_math::sin(-el.theta),
						coord1.x * fast_math::sin(-el.theta) + coord1.y * fast_math::cos(-el.theta));
					coord2 /= float_2(el.a, el.b);

					const auto id = atomic_fetch_add(&pointsCount[0], 1);
					points[id] = { unorm(fast_math::sqrt(coord2.x * coord2.x + coord2.y * coord2.y)),
								   fast_math::atan2(coord2.y, coord2.x) };
				}
			}
		});
		_circleSignTargetsSource.emplace_back(pointsCount[0], std::move(points), el.area);
	}
}

static std::array<uint_2, 8> mnPairs = {
	uint_2(2, 2),uint_2(3, 1),uint_2(3, 3),uint_2(4, 2),
	uint_2(4, 4),uint_2(5, 1),uint_2(5, 3),uint_2(5, 5),
};

void Recognizer::CalculateZernike()
{
	circleSignZernikes.clear();
	circleSignZernikes.resize(_circleSignTargetsSource.size());
	size_t cntId = 0;
	for (auto&& target : _circleSignTargetsSource)
	{
		const float area = target.area;
		const auto pointsCount = extent<1>(target.pointsCount);
		array_view<UnitCirclePoint, 1> pointsView(target.points);
		array_view<uint, 2> outputTex(_outputTex);
		array<float, 1> V(pointsCount, _acc_view), Vj(pointsCount, _acc_view);
		for (size_t i = 0; i < mnPairs.size(); i++)
		{
			const auto p = mnPairs[i].x;
			const int q = mnPairs[i].y;
			parallel_for_each(_acc_view, pointsCount, [=, &V, &Vj](index<1> index) restrict(amp)
			{
				const auto point = pointsView[index];
				auto r = ZernikeR(p, q, point.rho);
				auto v = ZernikeV(r, q, point.theta);
				V[index] = v.x;
				Vj[index] = v.y;
				const auto x = point.rho * fast_math::cos(point.theta) * 20.f + outputTex.extent[1] / 2.f;
				const auto y = point.rho * fast_math::sin(point.theta) * 20.f + outputTex.extent[0] / 2.f;
				outputTex(outputTex.extent[0] - y, x) = uint(0xFFFFFFFF);
			});

			const std::vector<float> cpuV(V), cpuVj(Vj);
			concurrency::combinable<double> cbV, cbVj;
			parallel_for_each(cpuV.begin(), cpuV.end(), [&](float value) {cbV.local() += value;});
			parallel_for_each(cpuVj.begin(), cpuVj.end(), [&](float value) {cbVj.local() += value;});
			const auto sV = cbV.combine(std::plus<double>()) * (2 * p + 2) / area;
			const auto sVj = cbVj.combine(std::plus<double>()) * (2 * p + 2) / area;
			const auto Z = sqrt(sV * sV + sVj * sVj);
			circleSignZernikes[cntId][i] = { p, q, (float)Z };
		}
		break;
		cntId++;
	}
}

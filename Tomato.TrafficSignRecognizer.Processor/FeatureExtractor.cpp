//
// Traffic Sign Recognizer
// 特征提取器
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "pch.h"
#include "FeatureExtractor.h"
#include <amp_math.h>
#include "algorithms.h"
#include "amprng/amp_sobol_rng.h"

using namespace NS_TSR_PRCSR;
using namespace Platform;
using namespace Platform::Collections;
using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Foundation;
using namespace Windows::Foundation::Collections;
using namespace Windows::Graphics::Imaging;
using namespace WRL;

static const uint32_t zero = 0;

FeatureExtractor::FeatureExtractor(unsigned int width, unsigned int height)
	:_targetImageExtent(height, width),
	_acc_view(concurrency::direct3d::create_accelerator_view(accelerator(), true)),
	_targetImage(_targetImageExtent, 16U, _acc_view), _outputTex(_targetImageExtent, _acc_view),
	_edgeTex(_targetImageExtent, 8, _acc_view), _tangentTex(_targetImageExtent, 32, _acc_view),
	_redTex(_targetImageExtent, 8, _acc_view), _edgePositions(_targetImageExtent.size(), _acc_view)
{

}

Windows::Foundation::IAsyncAction ^ FeatureExtractor::SetTarget(BitmapFrame ^ frame)
{
	return create_async([=]
	{
		return PrepareTargetImage(frame);
	});
}

IAsyncOperation<bool>^ FeatureExtractor::Recognize(Windows::Storage::Streams::IRandomAccessStream^ outputStream)
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

Windows::Foundation::IAsyncOperation<Windows::Foundation::Collections::IIterable<Windows::Foundation::Collections::IVectorView<ZernikeResult>^>^>^ FeatureExtractor::CaculateZernikes()
{
	return create_async([=]()
	{
		return task_from_result().then([=]() -> IIterable<IVectorView<ZernikeResult>^>^
		{
			auto ellipses = ref new Vector<IVectorView<ZernikeResult>^>();
			for (auto&& el : _circleSignZernikes)
				ellipses->Append(ref new VectorView<ZernikeResult>(el.begin(), el.end()));
			return ellipses;
		});
	});
}

concurrency::task<void> FeatureExtractor::PrepareTargetImage(BitmapFrame ^ frame)
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

concurrency::task<void> FeatureExtractor::FindContours()
{
	AbsorbRedTexels();
	FindEdgesAndTangent();

	texture_view<const unorm_4, 2> input(_targetImage);
	texture_view<const unorm, 2> redTex(_redTex);
	texture_view<const unorm, 2> edges(_edgeTex);
	texture_view<const float_2, 2> tangents(_tangentTex);
	array_view<uint, 2> outputTex(_outputTex);

	return FindEllipses().then([=]
	{
		parallel_for_each(_acc_view, _targetImageExtent, [redTex, outputTex, tangents, edges](index<2> index) restrict(amp)
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

void FeatureExtractor::FindEdgesAndTangent()
{
	uint radius = 3;
	texture_view<unorm, 2> edgeView(_edgeTex);
	texture_view<float_2, 2> tangentView(_tangentTex);
	texture_view<const unorm, 2> redTex(_redTex);
	texture_view<const unorm_4, 2> inputTex(_targetImage);
	parallel_for_each(_acc_view, redTex.extent, [redTex, edgeView, inputTex, tangentView, radius](index<2> index) restrict(amp)
	{
		auto edge = SusanTest(redTex, index, radius, 0.1f);
		edgeView.set(index, edge);

		if (edge > 0.f)
			tangentView.set(index, CalculateTangent(inputTex, index));
	});
}

void FeatureExtractor::AbsorbRedTexels()
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
			if (redReader.sample<filter_point>(coord(p + float_2(-1, +0), redReader.extent)) != 0.f &&
				redReader.sample<filter_point>(coord(p + float_2(+1, +0), redReader.extent)) != 0.f)
				tmpWriter.set(index, unorm(1.f));
			else if (redReader.sample<filter_point>(coord(p + float_2(+0, -1), redReader.extent)) != 0.f &&
				redReader.sample<filter_point>(coord(p + float_2(+0, +1), redReader.extent)) != 0.f)
				tmpWriter.set(index, unorm(1.f));
		}
		else
			tmpWriter.set(index, unorm(1.f));
	});

	copy(tmpTex, _redTex);
}

concurrency::array_view<const EllipsePoints, 1> FeatureExtractor::AbsorbEllipsePointPairs(uint32_t& count)
{
	const float height = _targetImageExtent[0];
	array_view<uint, 2> outputTex(_outputTex);

	texture_view<const unorm, 2> edgeView(_edgeTex);
	array<uint, 1> edgePointsCount(1, &zero, _acc_view, access_type_read);
	array_view<index<2>, 1> edgePositionsView(_edgePositions);

	parallel_for_each(_acc_view, edgeView.extent, [=, &edgePointsCount](index<2> index) restrict(amp)
	{
		if (edgeView[index] != 0.f)
		{
			auto idx = atomic_fetch_add(&edgePointsCount(0), 1);
			edgePositionsView(idx) = index;
		}
	});

	_edgePointsCount = edgePointsCount(0);
	const auto edgeSize = extent<1>(_edgePointsCount);
	const auto maxPointPairs = extent<1>(edgeSize * 10);
	array<uint32_t, 1> fitsCount(1, &zero, _acc_view, access_type_read);
	array<EllipsePoints, 1> ellipses(maxPointPairs, _acc_view);

	sobol_rng_collection<sobol_rng<2>, 1> sc_rng(_acc_view, maxPointPairs, 1);
	concurrency::graphics::texture_view<const float_2, 2> tangentView(_tangentTex);
	parallel_for_each(_acc_view, maxPointPairs, [=, &fitsCount, &ellipses](index<1> idx) restrict(amp)
	{
		// Get the sobol RNG 
		auto rng = sc_rng[idx];
		// Skip ahead to the right position
		rng.skip(sc_rng.direction_numbers(), idx[0]);

		float_2 points[2];
		float_2 tangent[2];
		for (uint32_t i = 0; i < 2; i++)
		{
			const auto id = int(fast_math::round(rng.get_single(i + 1) * edgeSize[0]));
			const auto ptIndex = edgePositionsView(id);
			points[i] = float_2(ptIndex[1], height - ptIndex[0]);
			tangent[i] = tangentView(ptIndex);
		}
		EllipsePoints ePoint;
		if (FindEllipsePoints(points, tangent, edgeView, tangentView, ePoint))
		{
			ellipses(atomic_fetch_add(&fitsCount(0), 1)) = ePoint;
			outputTex(index<2>(height - ePoint.p1.y, ePoint.p1.x)) = 0xFF00FF00;
			outputTex(index<2>(height - ePoint.p2.y, ePoint.p2.x)) = 0xFF00FF00;
			outputTex(index<2>(height - ePoint.p3.y, ePoint.p3.x)) = 0xFF00FF00;
		}
	});
	count = fitsCount[0];

	return ellipses;
}

task<void> FeatureExtractor::FindEllipses()
{
	const float height = _targetImageExtent[0];
	const float width = _targetImageExtent[1];
	array_view<uint, 2> outputTex(_outputTex);
	uint32_t pointsPairCount;
	auto pointsPairView = AbsorbEllipsePointPairs(pointsPairCount);
	const auto maxFitTimes = extent<1>(pointsPairCount * 100);

	array<uint32_t, 1> fitsCount(1, &zero, _acc_view, access_type_read);
	array<EllipseParam, 1> ellipses(maxFitTimes, _acc_view);

	static const uint32_t N = 5;
	sobol_rng_collection<sobol_rng<N>, 1> sc_rng(_acc_view, maxFitTimes, 1);

	parallel_for_each(_acc_view, maxFitTimes, [=, &fitsCount, &ellipses](index<1> idx) restrict(amp)
	{
		// Get the sobol RNG 
		auto rng = sc_rng[idx];
		// Skip ahead to the right position
		rng.skip(sc_rng.direction_numbers(), idx[0]);
		float_2 points[N];
		for (uint32_t i = 0; i < N; i++)
		{
			const auto id = int(fast_math::round(rng.get_single(i + 1) * pointsPairCount));
			const auto pair = pointsPairView(id);
			points[i] = pair.p1;
		}
		FitEllipse(points, width, height, fitsCount, ellipses);
	});
	const auto ellipseCount = fitsCount(0);
	array_view<index<2>, 1> edgePositionsView(_edgePositions);
	const auto edgeSize = extent<1>(_edgePointsCount);

	parallel_for_each(_acc_view, edgeSize, [=, &fitsCount, &ellipses](index<1> index)restrict(amp)
	{
		const auto p1 = edgePositionsView(index[0]);
		const float_2 point(p1[1], height - p1[0]);

		const auto ellipseCount = fitsCount(0);
		for (uint32_t i = 0; i < ellipseCount; i++)
		{
			auto& ellipse = ellipses(i);
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
		_ellipsesFit.reserve(ellipseCount);
		const float lambda = 0.5f;
		for (auto&& it : ellipsesSort)
		{
			const auto min = lambda * it.length;
			if (it.a > 10.f && it.b > 10.f && it.rank >= min)
				if (!std::any_of(_ellipsesFit.begin(), _ellipsesFit.end(), [&](const EllipseParam& el)
				{
					if (el.id != it.id &&
						fabs(el.x - it.x) < 3.f && fabs(el.y - it.y) < 3.f)
						return (el.a < it.a && el.b < it.b) || (fabs(el.a - it.a) < 3.f && fabs(el.b - it.b) < 3.f);
					return false;
				}))
					_ellipsesFit.emplace_back(it);
		}


		FillCircleSignTargetsSource();
		CalculateZernike();
	}

	return task_from_result();
}

void FeatureExtractor::FillCircleSignTargetsSource()
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

		const auto scale = std::round(std::min(el.a, el.b) * 2.f);
		texture<unorm, 2> ellipseTex(extent<2>(scale, scale), 8, _acc_view);
		texture_view<unorm, 2> ellipseTexWriter(ellipseTex);
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

					const auto rho = unorm(fast_math::sqrt(coord2.x * coord2.x + coord2.y * coord2.y));
					if (rho < 0.9f)
					{
						coord2 = (coord2 / 2 + 0.5f) * scale;
						ellipseTexWriter.set(concurrency::index<2>(fast_math::round(scale - coord2.y), fast_math::round(coord2.x)), unorm(1.0));
					}
				}
			}
		});
		_circleSignTargetsSource.emplace_back(ellipseTex);
		break;
	}
}

//static std::array<int_2, 9> mnPairs = {
//	int_2(1, 1),int_2(2, 0),int_2(3, 1),
//	int_2(4, 0),int_2(4, 2),int_2(5, 3),
//	int_2(6, 2),int_2(7, 3),int_2(7, 5),
//};
static std::array<int_2, 11> mnPairs = {
	int_2(1, 1),int_2(2, 0),int_2(2, 2),
	int_2(3, 1),int_2(3, 3),int_2(4, 0),
	int_2(4, 2),int_2(4, 4),int_2(5, 1),
	int_2(5, 3),int_2(5, 5)
};

void FeatureExtractor::CalculateZernike()
{
	static const float RhoScale = 50;
	static const float ThetaScale = 100;

	_circleSignZernikes.clear();
	_circleSignZernikes.resize(_circleSignTargetsSource.size());
	size_t cntId = 0;
	for (auto&& target : _circleSignTargetsSource)
	{
		extent<2> sampleExtent(1.f * RhoScale, 6.28f * ThetaScale);
		array_view<uint, 2> outputTex(_outputTex);
		array<double, 1> V(sampleExtent.size(), _acc_view), Vj(sampleExtent.size(), _acc_view);
		for (size_t i = 0; i < mnPairs.size(); i++)
		{
			const auto p = mnPairs[i].x;
			const int q = mnPairs[i].y;
			parallel_for_each(_acc_view, sampleExtent, [=, &V, &Vj](index<2> index) restrict(amp)
			{
				const auto rho = index[0] / RhoScale;
				const auto theta = index[1] / ThetaScale - 3.14f;
				const auto x = rho * fast_math::cos(theta) / 2.f + 0.5f;
				const auto y = 1.f - (rho * fast_math::sin(theta) / 2.f + 0.5f);
				concurrency::index<1> id(sampleExtent[1] * index[0] + index[1]);

				const auto mx = rho * fast_math::cos(theta) * 20.f + outputTex.extent[1] / 2.f;
				const auto my = rho * fast_math::sin(theta) * 20.f + outputTex.extent[0] / 2.f;
				if (target.sample<filter_point>(float_2(x, y)) == 1.f)
				{
					outputTex(outputTex.extent[0] - my, mx) = uint(0xFFFFFFFF);
					auto r = ZernikeR(p, q, rho);
					auto v = ZernikeV(r, q, theta);
					V[id] = v.x;
					Vj[id] = v.y;
				}
				else
				{
					V[id] = Vj[id] = 0;
				}
			});

			const std::vector<double> cpuV(V), cpuVj(Vj);
			concurrency::combinable<double> cbV, cbVj;
			parallel_for_each(cpuV.begin(), cpuV.end(), [&](double value) {cbV.local() += value;});
			parallel_for_each(cpuVj.begin(), cpuVj.end(), [&](double value) {cbVj.local() += value;});
			const auto sV = cbV.combine(std::plus<double>()) * (p + 1) / 3.1415 / (RhoScale * ThetaScale);
			const auto sVj = cbVj.combine(std::plus<double>()) * (p + 1) / 3.1415 / (RhoScale * ThetaScale);
			const auto Z = sqrt(sV * sV + sVj * sVj);
			_circleSignZernikes[cntId][i] = { uint32_t(p), q, (float)Z };
		}
		break;
		cntId++;
	}
}

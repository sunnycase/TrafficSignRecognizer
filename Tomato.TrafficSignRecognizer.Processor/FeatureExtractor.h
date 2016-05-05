//
// Traffic Sign Recognizer
// 特征提取器
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "common.h"
#include <amp.h>
#include <amp_graphics.h>
#include <ppltasks.h>
#include "models.h"
#pragma once

DEFINE_NS_TSR_PRCSR

public ref class FeatureExtractor sealed
{
public:
	FeatureExtractor(unsigned int width, unsigned int height);

	Windows::Foundation::IAsyncAction^ SetTarget(Windows::Graphics::Imaging::BitmapFrame ^ frame);
	Windows::Foundation::IAsyncOperation<bool>^ Recognize(Windows::Storage::Streams::IRandomAccessStream^ outputStream);
	Windows::Foundation::IAsyncOperation<bool>^ Recognize();
	Windows::Foundation::IAsyncOperation<Windows::Foundation::Collections::IIterable<Windows::Foundation::Collections::IVectorView<ZernikeResult>^>^>^ CaculateZernikes();

private:
	concurrency::task<void> PrepareTargetImage(Windows::Graphics::Imaging::BitmapFrame ^ frame);
	concurrency::task<void> FindContours();
	// 查找边缘并计算切线方向
	void FindEdgesAndTangent();
	void AbsorbRedTexels();
	concurrency::array_view<const EllipsePoints, 1> AbsorbEllipsePointPairs(uint32_t& count);
	concurrency::task<void> FindEllipses();
	void FillCircleSignTargetsSource();
	void CalculateZernike();

private:
	concurrency::extent<2> _targetImageExtent;
	concurrency::accelerator_view _acc_view;
	concurrency::graphics::texture<concurrency::graphics::unorm_4, 2> _targetImage;
	concurrency::array<concurrency::graphics::uint, 2> _outputTex;
	concurrency::graphics::texture<concurrency::graphics::unorm, 2> _edgeTex;
	concurrency::graphics::texture<concurrency::graphics::unorm, 2> _redTex;
	concurrency::graphics::texture<concurrency::graphics::float_2, 2> _tangentTex;
	std::vector<EllipseParam> _ellipsesFit;
	std::vector<concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>> _circleSignTargetsSource;
	std::vector<std::array<ZernikeResult, 11>> _circleSignZernikes;
	uint32_t _edgePointsCount;
	concurrency::array<concurrency::index<2>, 1> _edgePositions;
};

END_NS_TSR_PRCSR
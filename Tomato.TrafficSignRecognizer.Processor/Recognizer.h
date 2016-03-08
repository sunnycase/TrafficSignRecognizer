//
// Traffic Sign Recognizer
// Recognizer
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "common.h"
#include <amp.h>
#include <amp_graphics.h>
#include <ppltasks.h>
#pragma once

DEFINE_NS_TSR_PRCSR

public ref class Recognizer sealed
{
public:
	Recognizer(unsigned int width, unsigned int height);

	Windows::Foundation::IAsyncAction^ SetTarget(Windows::Graphics::Imaging::BitmapFrame ^ frame);
	Windows::Foundation::IAsyncOperation<bool>^ Recognize(Windows::Storage::Streams::IRandomAccessStream^ outputStream);
private:
	concurrency::task<void> PrepareTargetImage(Windows::Graphics::Imaging::BitmapFrame ^ frame);
	concurrency::task<void> FindContours();
	// 查找边缘并计算切线方向
	void FindEdgesAndTangent();
	concurrency::task<void> FindEllipses();
private:
	concurrency::extent<2> _targetImageExtent;
	concurrency::accelerator_view _acc_view;
	concurrency::graphics::texture<concurrency::graphics::unorm_4, 2> _targetImage;
	concurrency::array<concurrency::graphics::uint, 2> _outputTex;
	concurrency::graphics::texture<concurrency::graphics::unorm, 2> _edgeTex;
	concurrency::graphics::texture<float, 2> _tangentTex;
};

END_NS_TSR_PRCSR
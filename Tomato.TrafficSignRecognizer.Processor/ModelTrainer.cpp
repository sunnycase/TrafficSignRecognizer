//
// Traffic Sign Recognizer
// 模型训练器
// 作者：SunnyCase
// 创建日期：2016-05-5
#include "pch.h"
#include "ModelTrainer.h"
#include "algorithms.h"
#include "SignModel.h"

using namespace NS_TSR_PRCSR;
using namespace Platform;
using namespace Platform::Collections;
using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Foundation;
using namespace Windows::Foundation::Collections;
using namespace Windows::Graphics::Imaging;
using namespace WRL;

ModelTrainer::ModelTrainer()
{

}

void ModelTrainer::AddSamples(int className, Windows::Foundation::Collections::IVectorView<Windows::Storage::StorageFile^>^ samples)
{
	_samples.emplace_back(className, samples);
}

Windows::Foundation::IAsyncOperation<SignModel^>^ ModelTrainer::Train()
{
	return create_async([=] {return TrainIntern(); });
}

concurrency::task<SignModel^> ModelTrainer::TrainIntern()
{
	std::vector<double> classNames;
	classNames.reserve(_samples.size());
	for (auto&& sampleGroup : _samples)
	{
		classNames.emplace_back(sampleGroup.className);
		std::vector<std::array<svm_node, 11>> nodes;
		nodes.reserve(sampleGroup.samples->Size);
		for (auto&& sample : sampleGroup.samples)
		{
			auto rawFeatures = co_await GetFeature(sample);
			if (rawFeatures[0])
				nodes.emplace_back(std::array<svm_node, 11>{ { 0, rawFeatures[0]} });
		}
	}
	return ref new SignModel(nullptr);
}

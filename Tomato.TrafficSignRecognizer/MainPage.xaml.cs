using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using System.Threading.Tasks;
using Tomato.TrafficSignRecognizer.Processor;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.Graphics.Imaging;
using Windows.Storage;
using Windows.Storage.Streams;
using Windows.UI.Popups;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Controls.Primitives;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Navigation;

//“空白页”项模板在 http://go.microsoft.com/fwlink/?LinkId=402352&clcid=0x409 上有介绍

namespace Tomato.TrafficSignRecognizer
{
    /// <summary>
    /// 可用于自身或导航至 Frame 内部的空白页。
    /// </summary>
    public sealed partial class MainPage : Page
    {
        readonly Dictionary<int, List<float[]>> features = new Dictionary<int, List<float[]>>();
        int count = 0;
        public MainPage()
        {
            this.InitializeComponent();
            Loaded += MainPage_Loaded;
        }

        private async void MainPage_Loaded(object sender, RoutedEventArgs e)
        {
            pb_Feature.Maximum = 2 * 6;
            var watch = new Stopwatch();
            watch.Start();
            await LoadFeatures();
            WriteTrainData();
            watch.Stop();
            tb_Status.Text = $"共提取 {count} 张图片，耗时 {watch.Elapsed}";
            //var zernikes = (await ext.CaculateZernikes()).Select(o => o.ToList()).ToList();
            //zernikes.ToString();
            //var outputFile = await ApplicationData.Current.LocalCacheFolder.CreateFileAsync("output.jpg", CreationCollisionOption.ReplaceExisting);
            //using (var stream = await outputFile.OpenAsync(FileAccessMode.ReadWrite))
            //{
            //    await ext.Recognize(stream);
            //}
        }

        private async void WriteTrainData()
        {
            var outputFile = await ApplicationData.Current.LocalCacheFolder.CreateFileAsync("train.txt", CreationCollisionOption.ReplaceExisting);
            using (var dataWriter = new DataWriter(await outputFile.OpenAsync(FileAccessMode.ReadWrite)))
            {
                foreach (var label in features)
                {
                    foreach (var item in label.Value)
                    {
                        dataWriter.WriteString(string.Format($"{label.Key} {string.Concat(item.Select((v, i) => $"{i + 1}:{v} "))}\r\n"));
                    }
                }
                await dataWriter.StoreAsync();
                await dataWriter.FlushAsync();
            }
            var msg = new MessageDialog("Feature extraction completed.");
            await msg.ShowAsync();
        }

        private async Task LoadFeatures()
        {
            for (int i = 0; i < 2; i++)
            {
                await LoadFeatures(i);
            }
        }

        private async Task LoadFeatures(int label)
        {
            List<float[]> features = new List<float[]>();
            for (int i = 21; i < 27; i++)
            {
                var feature = await ExtractFeature(label, i);
                if (feature != null)
                    features.Add(feature);
                Debug.WriteLine($"Label: {label} Id: {i}");
                pb_Feature.Value++;
                count++;
            }
            this.features.Add(label, features);
        }

        async Task<float[]> ExtractFeature(int label, int id)
        {
            var file = await StorageFile.GetFileFromApplicationUriAsync(new Uri($"ms-appx:///Assets/{label}/{id}.jpg"));
            var decoder = await BitmapDecoder.CreateAsync(await file.OpenReadAsync());
            var frame = await decoder.GetFrameAsync(0);
            var ext = new FeatureExtractor(frame.OrientedPixelWidth, frame.OrientedPixelHeight);
            await ext.SetTarget(frame);
            var outputFile = await ApplicationData.Current.LocalCacheFolder.CreateFileAsync($"{label}-{id}.jpg", CreationCollisionOption.ReplaceExisting);
            using (var stream = await outputFile.OpenAsync(FileAccessMode.ReadWrite))
            {
                await ext.Recognize(stream);
            }
            return (await ext.CaculateZernikes()).FirstOrDefault()?.Select(o => o.z)?.ToArray();
        }
    }
}

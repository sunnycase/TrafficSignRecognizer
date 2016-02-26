using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using Tomato.TrafficSignRecognizer.Processor;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.Graphics.Imaging;
using Windows.Storage;
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
        public MainPage()
        {
            this.InitializeComponent();
            Loaded += MainPage_Loaded;
        }

        private async void MainPage_Loaded(object sender, RoutedEventArgs e)
        {
            var file = await StorageFile.GetFileFromApplicationUriAsync(new Uri("ms-appx:///Assets/stop.jpg"));
            var decoder = await BitmapDecoder.CreateAsync(await file.OpenReadAsync());
            var frame = await decoder.GetFrameAsync(0);
            var rec = new Recognizer(frame.OrientedPixelWidth, frame.OrientedPixelHeight);
            await rec.SetTarget(frame);
            var outputFile = await ApplicationData.Current.LocalCacheFolder.CreateFileAsync("output.jpg", CreationCollisionOption.ReplaceExisting);
            using (var stream = await outputFile.OpenAsync(FileAccessMode.ReadWrite))
            {
                await rec.Recognize(stream);
            }
        }
    }
}

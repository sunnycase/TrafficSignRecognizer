using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Tomato.Mvvm;
using Tomato.TrafficSignRecognizer.Views;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;

namespace Tomato.TrafficSignRecognizer.ViewModels
{
    class MainViewModel : BindableBase
    {
        private Action _menuSwitch;
        private Frame _navigationService;

        public MainViewModel()
        {

        }

        public void NavigateToTrain()
        {
            _navigationService.Navigate(typeof(TrainPage));
        }

        public void NavigateToRecognize()
        {
            _navigationService.Navigate(typeof(RecognizePage));
        }

        public void SwitchMenu()
        {
            _menuSwitch();
        }

        public void SetupNavigationService(object sender, RoutedEventArgs e)
        {
            _navigationService = (Frame)sender;
            NavigateToTrain();
        }

        public void SetupMenu(object sender, RoutedEventArgs e)
        {
            var menu = (SplitView)sender;
            _menuSwitch = () => menu.IsPaneOpen = !menu.IsPaneOpen;
        }
    }
}

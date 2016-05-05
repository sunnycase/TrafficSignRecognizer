using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Tomato.Mvvm;
using Windows.Storage;
using Windows.Storage.Pickers;

namespace Tomato.TrafficSignRecognizer.ViewModels
{
    class TrainViewModel : BindableBase
    {
        private string _sampleFolderPath;
        public string SampleFolderPath
        {
            get { return _sampleFolderPath; }
            set { SetProperty(ref _sampleFolderPath, value); }
        }

        private bool _canOperate = true;
        public bool CanOperate
        {
            get { return _canOperate; }
            set { SetProperty(ref _canOperate, value); }
        }

        private bool _analyzeCompleted;
        public bool AnalyzeCompleted
        {
            get { return _analyzeCompleted; }
            set { SetProperty(ref _analyzeCompleted, value); }
        }

        private bool _isBusy = false;
        public bool IsBusy
        {
            get { return _isBusy; }
            set { SetProperty(ref _isBusy, value); }
        }

        private int _classCount;
        public int ClassCount
        {
            get { return _classCount; }
            set { SetProperty(ref _classCount, value); }
        }

        private int _samplesCount;
        public int SamplesCount
        {
            get { return _samplesCount; }
            set { SetProperty(ref _samplesCount, value); }
        }

        private double _trainPercent = 0.7;
        public double TrainPercent
        {
            get { return _trainPercent * 100; }
            set { SetProperty(ref _trainPercent, value / 100); }
        }

        private bool _trainCompleted;
        public bool TrainCompleted
        {
            get { return _trainCompleted; }
            set { SetProperty(ref _trainCompleted, value); }
        }

        private double _accuracy;
        public double Accuracy
        {
            get { return _accuracy * 100; }
            private set { SetProperty(ref _accuracy, value); }
        }

        private TimeSpan _elapsedTrainTime;
        public TimeSpan ElapsedTrainTime
        {
            get { return _elapsedTrainTime; }
            set { SetProperty(ref _elapsedTrainTime, value); }
        }

        private StorageFolder _sampleFolder;
        private readonly Dictionary<string, IReadOnlyCollection<StorageFile>> _sampleFiles = new Dictionary<string, IReadOnlyCollection<StorageFile>>();

        public TrainViewModel()
        {

        }

        public async void Browse()
        {
            try
            {
                CanOperate = false;
                var picker = new FolderPicker() { ViewMode = PickerViewMode.Thumbnail };
                picker.FileTypeFilter.Add(".");
                var folder = await picker.PickSingleFolderAsync();
                if (folder != null && !(_sampleFolder?.IsEqual(folder) ?? false))
                {
                    _sampleFolder = folder;
                    SampleFolderPath = folder.Path;
                    await AnalyzeSampleFolder();
                }
            }
            finally
            {
                CanOperate = true;
            }
        }

        public async void Train()
        {
            try
            {
                IsBusy = true;
                AnalyzeCompleted = false;



                TrainCompleted = true;
            }
            finally
            {
                IsBusy = false;
            }
        }

        private async Task AnalyzeSampleFolder()
        {
            try
            {
                IsBusy = true;
                AnalyzeCompleted = false;
                TrainCompleted = false;
                 var subFolders = await _sampleFolder.GetFoldersAsync();
                _sampleFiles.Clear();
                int totalFiles = 0;
                foreach (var subFolder in subFolders)
                {
                    var files = await subFolder.GetFilesAsync();
                    _sampleFiles.Add(subFolder.Name, files);
                    totalFiles += files.Count;
                }
                SamplesCount = totalFiles;
                ClassCount = _sampleFiles.Count;
                AnalyzeCompleted = true;
            }
            finally
            {
                IsBusy = false;
            }
        }
    }
}

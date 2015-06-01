using System;
using System.Linq;
using System.Collections.ObjectModel;
using System.Threading.Tasks;
using System.Windows.Input;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.MobileApp.Infrastructure;
using RiskAnalysisTool.Requests;
using RiskAnalysisTool.Tasks;
using XLabs.Forms.Mvvm;
using System.Net;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class TaskDetailViewModel : FunctionViewModel<ComputationTask>, IModelViewModel<ComputationTask>
    {
        private readonly AsyncCommand _computeCommand;
        private PortfolioViewModel _portfolioViewModel;
        private MarketDataListViewModel _marketDataListViewModel;
        private TaskSettingsViewModel _taskSettingsViewModel;

        public TaskDetailViewModel()
        {
            _computeCommand = new AsyncCommand(this.Compute);
            _portfolioViewModel = new PortfolioViewModel();
            _marketDataListViewModel = new MarketDataListViewModel();
            _taskSettingsViewModel = new TaskSettingsViewModel();
        }

        public PortfolioViewModel PortfolioViewModel
        {
            get { return _portfolioViewModel; }
            set { SetProperty(ref _portfolioViewModel, value); }
        }

        public MarketDataListViewModel MarketDataListViewModel
        {
            get { return _marketDataListViewModel; }
            set { SetProperty(ref _marketDataListViewModel, value); }
        }

        public TaskSettingsViewModel TaskSettingsViewModel
        {
            get { return _taskSettingsViewModel; }
            set { SetProperty(ref _taskSettingsViewModel, value); }
        }

        public ICommand ComputeCommand
        {
            get { return _computeCommand; }
        }

        public async Task Compute()
        {
            Task exTask = null;
            try
            {
                if (this.TaskSettingsViewModel.ComputeCreditVaR)
                {
                    ComputationRequest request = new CreditVaRRequest();
                    PushRequest(request);
                    await ServiceHelper.CallService<ComputationRequest, object>("api/Compute/Request");
                }

                if (this.TaskSettingsViewModel.ComputeUCva)
                {
                    ComputationRequest request = new UCvaRequest();
                    PushRequest(request);
                    await ServiceHelper.CallService<ComputationRequest, object>("api/Compute/Request");
                }

                if (this.TaskSettingsViewModel.ComputeBva)
                {
                    ComputationRequest request = new BvaRequest();
                    PushRequest(request);
                    await ServiceHelper.CallService<ComputationRequest, object>("api/Compute/Request");
                }

                if (this.TaskSettingsViewModel.ComputeFva)
                {
                    ComputationRequest request = new TvaRequest();
                    PushRequest(request);
                    await ServiceHelper.CallService<ComputationRequest, object>("api/Compute/Request");
                }

                exTask = App.Current.MainPage.DisplayAlert("Tasks Added", "Please monitor your email address", "Ok");
            }
            catch (WebException ex)
            {
                exTask = App.Current.MainPage.DisplayAlert("Error", ex.Message, "Ok");
            }

            if (exTask != null)
            {
                await exTask;
            }
        }

        private void PushRequest(ComputationRequest request)
        {
            request.Portfolio = this.PortfolioViewModel.Instruments.ToArray();
            request.MarketData = this.MarketDataListViewModel.Instruments.ToArray();
            request.DomesticCurrency = this.TaskSettingsViewModel.DomesticCurrency;
            request.NotificationEmail = this.TaskSettingsViewModel.EmailAddress;

            foreach (var inst in request.Portfolio)
            {
                ProcessSymbols(inst, this.TaskSettingsViewModel.DomesticCurrency);
            }

            foreach (var inst in request.MarketData)
            {
                ProcessSymbols(inst, this.TaskSettingsViewModel.DomesticCurrency);
            }
        }

        private static void ProcessSymbols(Instrument inst, string domesticCurrency)
        {
            if (inst is Swap)
            {
                Swap swap = inst as Swap;
                if (string.IsNullOrEmpty(swap.PaySymbol))
                {
                    swap.PaySymbol = "$" + domesticCurrency;
                }
                if (string.IsNullOrEmpty(swap.ReceiveSymbol))
                {
                    swap.ReceiveSymbol = "$" + domesticCurrency;
                }
            }

            if (inst is Correlation)
            {
                Correlation corr = inst as Correlation;
                if (string.IsNullOrEmpty(corr.Symbol1))
                {
                    corr.Symbol1 = "$" + domesticCurrency;
                }
                if (string.IsNullOrEmpty(corr.Symbol2))
                {
                    corr.Symbol2 = "$" + domesticCurrency;
                }
            }

            if (string.IsNullOrEmpty(inst.Symbol))
            {
                inst.Symbol = "$" + domesticCurrency;
            }
        }

        public void LoadModel(ComputationTask task)
        {

            this.PortfolioViewModel.Instruments = new ObservableCollection<Instrument>(task.Request.Portfolio);
            this.MarketDataListViewModel.Instruments = new ObservableCollection<Instrument>(task.Request.MarketData);
            this.TaskSettingsViewModel.DomesticCurrency = task.Request.DomesticCurrency;

            if (task.Request.GetType() == typeof(CreditVaRRequest))
            {
                this.TaskSettingsViewModel.ComputeCreditVaR = true;
            }
            else if (task.Request.GetType() == typeof(UCvaRequest))
            {
                this.TaskSettingsViewModel.ComputeUCva = true;
            }
            else if (task.Request.GetType() == typeof(BvaRequest))
            {
                this.TaskSettingsViewModel.ComputeBva = true;
            }
            else if (task.Request.GetType() == typeof(TvaRequest))
            {
                this.TaskSettingsViewModel.ComputeFva = true;
            }

            foreach (var inst in task.Request.Portfolio)
            {
                if (inst is Swap)
                {
                    Swap swap = inst as Swap;
                    if (string.IsNullOrEmpty(swap.PaySymbol))
                    {
                        swap.PaySymbol = "$" + task.Request.DomesticCurrency;
                    }
                    if (string.IsNullOrEmpty(swap.ReceiveSymbol))
                    {
                        swap.ReceiveSymbol = "$" + task.Request.DomesticCurrency;
                    }
                }

                if (string.IsNullOrEmpty(inst.Symbol))
                {
                    inst.Symbol = "$" + task.Request.DomesticCurrency;
                }
            }
        }

        public void UpdateModel(ComputationTask task)
        {
            throw new NotImplementedException();
        }
    }
}
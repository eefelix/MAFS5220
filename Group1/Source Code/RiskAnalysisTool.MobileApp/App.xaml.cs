using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;
using RiskAnalysisTool.MobileApp.ViewModels;
using RiskAnalysisTool.MobileApp.Views;
using Xamarin.Forms;
using XLabs.Forms.Mvvm;
using XLabs.Forms.Services;
using XLabs.Ioc;
using XLabs.Platform.Device;
using XLabs.Platform.Mvvm;
using XLabs.Platform.Services;
using XLabs.Serialization;

namespace RiskAnalysisTool.MobileApp
{
    public partial class App
    {
        private static readonly Lazy<App> _instance = new Lazy<App>();

        public App()
        {
            InitializeComponent();
        }

        protected override void OnStart()
        {
            base.OnStart();
            ViewFactory.Register<BondDetailView, BondDetailViewModel>();
            ViewFactory.Register<CorrelationDetailView, CorrelationDetailViewModel>();
            ViewFactory.Register<CreditDefaultSwapDetailView, CreditDefaultSwapDetailViewModel>();
            ViewFactory.Register<CrossCurrencySwapDetailView, CrossCurrencySwapDetailViewModel>();
            ViewFactory.Register<EquitySwapDetailView, EquitySwapDetailViewModel>();
            ViewFactory.Register<InterestRateSwapDetailView, InterestRateSwapDetailViewModel>();
            ViewFactory.Register<LoginView, LoginViewModel>();
            ViewFactory.Register<MarketDataListView, MarketDataListViewModel>();
            ViewFactory.Register<PortfolioView, PortfolioViewModel>();
            ViewFactory.Register<TaskDetailView, TaskDetailViewModel>();
            ViewFactory.Register<TaskSettingsView, TaskSettingsViewModel>();
            ViewFactory.Register<VolatilityDetailView, VolatilityDetailViewModel>();

            this.MainPage = new NavigationPage(
                ViewFactory.CreatePage<TaskDetailViewModel, TaskDetailView>() as Page);

            Resolver.Resolve<IDependencyContainer>()
                .Register<INavigationService>(t => new NavigationService(this.MainPage.Navigation))
                .Register<INavigation>(t => this.MainPage.Navigation);        }

        public static App Instance
        {
            get { return _instance.Value; }
        }
    }
}

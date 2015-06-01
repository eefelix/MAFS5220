using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RiskAnalysisTool.MobileApp.ViewModels;
using Xamarin.Forms;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.Views
{
	public partial class TaskDetailView : TabbedPage
    {
        public TaskDetailView()
        {
            InitializeComponent();
            this.Children.Add(ViewFactory.CreatePage<TaskSettingsViewModel, TaskSettingsView>(
                (vm, v) =>
                {
                    v.SetBinding<TaskDetailViewModel>(
                        TaskSettingsView.BindingContextProperty,
                        e => e.TaskSettingsViewModel
                        );
                }) as TaskSettingsView);
            this.Children.Add(ViewFactory.CreatePage<PortfolioViewModel, PortfolioView>(
                (vm, v) =>
                {
                    v.SetBinding<TaskDetailViewModel>(
                        TaskSettingsView.BindingContextProperty,
                        e => e.PortfolioViewModel
                        );
                }) as PortfolioView);
            this.Children.Add(ViewFactory.CreatePage<MarketDataListViewModel, MarketDataListView>(
                (vm, v) =>
                {
                    v.SetBinding<TaskDetailViewModel>(
                        TaskSettingsView.BindingContextProperty,
                        e => e.MarketDataListViewModel
                        );
                }) as MarketDataListView);
        }
    }
}

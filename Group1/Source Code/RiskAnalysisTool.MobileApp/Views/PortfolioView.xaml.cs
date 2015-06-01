using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Xamarin.Forms;
using RiskAnalysisTool.MobileApp.ViewModels;

namespace RiskAnalysisTool.MobileApp.Views
{
	public partial class PortfolioView
    {
        public PortfolioView()
        {
            InitializeComponent();
        }

        protected virtual async void OnAddButtonClicked(object sender, EventArgs e)
        {
            PortfolioViewModel viewModel = (PortfolioViewModel)this.BindingContext;
            int? index = await this.DisplaySelectionBox(
                "Add Instrument",
                "Interest Rate Swap (IRS)", "Cross Currency Swap", "Equity Swap"
                );
            if (index == null)
            {
                return;
            }

            string instrumentType;
            switch (index)
            {
                case 0:
                    instrumentType = "IRS";
                    break;
                case 1:
                    instrumentType = "CCS";
                    break;
                case 2:
                    instrumentType = "Equity Swap";
                    break;
                default:
                    throw new NotSupportedException();
            }

            await viewModel.AddInstrument(instrumentType);
        }

        protected virtual void OnEditButtonClicked(object sender, EventArgs e)
        {
            MenuItem mi = (MenuItem)sender;
            PortfolioViewModel viewModel = (PortfolioViewModel)this.BindingContext;
            viewModel.EditInstrumentCommand.Execute(mi.CommandParameter);
        }

        protected virtual void OnDeleteButtonClicked(object sender, EventArgs e)
        {
            MenuItem mi = (MenuItem)sender;
            PortfolioViewModel viewModel = (PortfolioViewModel)this.BindingContext;
            viewModel.RemoveInstrumentCommand.Execute(mi.CommandParameter);
        }
        protected virtual void OnListItemTapped(object sender, ItemTappedEventArgs e)
        {
            PortfolioViewModel viewModel = (PortfolioViewModel)this.BindingContext;
            viewModel.EditInstrumentCommand.Execute(e.Item);
	        ListView list = sender as ListView;
            list.SelectedItem = null;
     }

	    protected virtual void OnListItemSelected(object sender, SelectedItemChangedEventArgs e)
	    {
	        ListView list = sender as ListView;
            list.SelectedItem = null;
	    }
    }
}

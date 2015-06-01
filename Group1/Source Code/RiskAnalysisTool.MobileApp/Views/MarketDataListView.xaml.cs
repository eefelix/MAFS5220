using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.MobileApp.ViewModels;
using Xamarin.Forms;

namespace RiskAnalysisTool.MobileApp.Views
{
    public partial class MarketDataListView : ContentPage
    {
        public MarketDataListView()
        {
            InitializeComponent();
        }

        protected virtual async void OnAddButtonClicked(object sender, EventArgs e)
        {
            MarketDataListViewModel viewModel = (MarketDataListViewModel)this.BindingContext;
            int? index = await this.DisplaySelectionBox(
                "Add Instrument",
                "Zero-coupon bond", "Credit default swap (CDS)", "Correlation"
                );
            if (index == null)
            {
                return;
            }

            string instrumentType;
            switch (index)
            {
                case 0:
                    instrumentType = "Bond";
                    break;
                case 1:
                    instrumentType = "CDS";
                    break;
                case 2:
                    instrumentType = "Correlation";
                    break;
                default:
                    throw new NotSupportedException();
            }

            await viewModel.AddInstrument(instrumentType);
        }

        protected virtual void OnEditButtonClicked(object sender, EventArgs e)
        {
            MenuItem mi = (MenuItem)sender;
            MarketDataListViewModel viewModel = (MarketDataListViewModel)this.BindingContext;
            viewModel.EditInstrumentCommand.Execute(mi.CommandParameter);
        }

        protected virtual void OnDeleteButtonClicked(object sender, EventArgs e)
        {
            MenuItem mi = (MenuItem)sender;
            MarketDataListViewModel viewModel = (MarketDataListViewModel)this.BindingContext;
            viewModel.RemoveInstrumentCommand.Execute(mi.CommandParameter);
        }
        protected virtual void OnListItemTapped(object sender, ItemTappedEventArgs e)
        {
            MarketDataListViewModel viewModel = (MarketDataListViewModel)this.BindingContext;
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

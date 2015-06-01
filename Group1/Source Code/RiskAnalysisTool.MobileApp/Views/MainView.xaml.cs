using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Xamarin.Forms;

namespace RiskAnalysisTool.MobileApp.Views
{
	public partial class MainView
	{
		public MainView ()
		{
			InitializeComponent ();
            this.Detail = new NavigationPage(new InterestRateSwapDetailView());
		    this.IsPresented = false;
		}
	}
}

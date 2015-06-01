using System.Threading.Tasks;
using RiskAnalysisTool.MobileApp.Infrastructure;
using Xamarin.Forms;
using XLabs.Forms.Controls;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class LoginViewModel : ViewModel
    {
        public LoginViewModel()
        {
            LoginWithWeiboCommand = new AsyncCommand(LoginWithWeibo);
        }

        public AsyncCommand LoginWithWeiboCommand { get; private set; }

        private async Task LoginWithWeibo()
        {
            await Application.Current.MainPage.DisplayAlert("TODO", "TODO: Login using Weibo OAuth 2.0", "Got it");
            var popup = new PopupLayout
            {
                Content = new WebView
                {
                    Source = "http://open.weibo.com/wiki/Oauth2/authorize"
                }
            };
            popup.ShowPopup(popup);
        }
    }
}
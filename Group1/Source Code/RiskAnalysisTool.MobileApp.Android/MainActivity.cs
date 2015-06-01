using System;
using Android.App;
using Android.Content.PM;
using Android.Runtime;
using Android.Views;
using Android.Widget;
using Android.OS;
using Newtonsoft.Json;
using Xamarin.Forms.Platform.Android;
using XLabs.Ioc;
using XLabs.Platform.Device;
using XLabs.Platform.Mvvm;
using XLabs.Serialization;

namespace RiskAnalysisTool.MobileApp.Android
{
    [Activity(Label = "RiskAnalysisTool", ConfigurationChanges = ConfigChanges.ScreenSize | ConfigChanges.Orientation)]
    public class MainActivity : FormsApplicationActivity
    {
        protected override void OnCreate(Bundle bundle)
        {
            base.OnCreate(bundle);
            Xamarin.Forms.Forms.Init(this, bundle);

            var resolverContainer = new SimpleContainer();

            resolverContainer
                .Register<IDevice>(t => AndroidDevice.CurrentDevice)
                .Register<IDisplay>(t => t.Resolve<IDevice>().Display)
                .Register<IDependencyContainer>(resolverContainer);

            Resolver.SetResolver(resolverContainer.GetResolver());
            
            var app = new XFormsApp<RiskAnalysisTool.MobileApp.App>(RiskAnalysisTool.MobileApp.App.Instance);
            resolverContainer.Register<IXFormsApp>(app);

            this.LoadApplication(RiskAnalysisTool.MobileApp.App.Instance);
        }
    }
}


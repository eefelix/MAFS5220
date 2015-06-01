//using System;
//using System.Collections.Generic;
//using System.ComponentModel;
//using System.Reflection;
//using System.Runtime.CompilerServices;
//using System.Text;
//using System.Text.RegularExpressions;
//using System.Threading.Tasks;
//using Xamarin.Forms;
//using XLabs;
//using XLabs.Forms.Mvvm;

//namespace RiskAnalysisTool.MobileApp.Infrastructure
//{
//    /// <summary>
//    /// Base class for all view models
//    /// </summary>
//    public abstract class ViewModelEx : XLabs.Forms.Mvvm.NavigationAwareViewModel
//    {
//        private ViewModel _pendingViewModel;
//        private Page _pendingView;
//        private TaskCompletionSource<ViewModel> _tcs;

//        public Task<ViewModel> InvokeAsync<TViewModel>(bool animation = true)
//            where TViewModel : ViewModel
//        {
//            return InvokeAsync<TViewModel>(null, animation);
//        }

//        public Task<ViewModel> InvokeAsync<TViewModel>(Action<TViewModel, Page> activateAction, bool animation = true)
//            where TViewModel : ViewModel
//        {
//            var tcs = new TaskCompletionSource<ViewModel>();
//            _tcs = tcs;
//            this.Navigation.PushAsync<TViewModel>(
//                (vm, v) =>
//                {
//                    if (activateAction != null)
//                    {
//                        activateAction(vm, v);
//                    }
//                    _pendingViewModel = vm;
//                    _pendingView = v;
//                    v.Disappearing += ViewDisapearing;
//                });
//            return tcs.Task;
//        }

//        public Task<ViewModel> InvokeModalAsync<TViewModel>(bool animation = true)
//            where TViewModel : ViewModel
//        {
//            return InvokeModalAsync<TViewModel>(null, animation);
//        }

//        public Task<ViewModel> InvokeModalAsync<TViewModel>(Action<TViewModel, Page> activateAction, bool animation = true)
//            where TViewModel : ViewModel
//        {
//            var tcs = new TaskCompletionSource<ViewModel>();
//            _tcs = tcs;
//            this.Navigation.PushModalAsync<TViewModel>(
//                (vm, v) =>
//                {
//                    if (activateAction != null)
//                    {
//                        activateAction(vm, v);
//                    }
//                    _pendingViewModel = vm;
//                    _pendingView = v;
//                    v.Disappearing += ViewDisapearing;
//                });
//            return tcs.Task;
//        }

//        private void ViewDisapearing(object sender, EventArgs e)
//        {
//            _tcs.SetResult(_pendingViewModel);
//            _pendingView.Disappearing -= ViewDisapearing;
//            _pendingView = null;
//            _pendingViewModel = null;
//            _tcs = null;
//        }
//    }
//}

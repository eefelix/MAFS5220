using System;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;
using RiskAnalysisTool.MobileApp.Infrastructure;
using Xamarin.Forms;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public static class ViewModelExtensions
    {
        private class InvokeHelper<TViewModel, TResult> : IDisposable
            where TViewModel : FunctionViewModel<TResult>
        {
            private TViewModel _pendingViewModel;
            private Page _pendingView;
            private TaskCompletionSource<TResult> _tcs;

            public void Attach(TViewModel viewModel, Page view)
            {
                this.Detatch();
                _pendingView = view;
                _pendingViewModel = viewModel;
                _tcs = new TaskCompletionSource<TResult>();
                _pendingView.Disappearing += ViewDisapearing;
            }

            private void Detatch()
            {
                if (_pendingViewModel != null)
                {
                    _pendingViewModel = null;
                    _pendingView = null;
                }

                if (_pendingView != null)
                {
                    _pendingView.Disappearing -= ViewDisapearing;
                    _pendingView = null;
                }


                _tcs = null;

            }

            public TViewModel ViewModel
            {
                get { return _pendingViewModel; }
            }

            public Page View
            {
                get { return _pendingView; }
            }

            public TaskCompletionSource<TResult> TaskCompletionSource
            {
                get { return _tcs; }
            }
            private void ViewDisapearing(object sender, EventArgs e)
            {
                _tcs.SetResult(_pendingViewModel.Result);
                this.Detatch();
            }
            public void Dispose()
            {
                this.Detatch();
            }
        }

        public static Task<TResult> InvokeAsync<TViewModel, TResult>(this ViewModel viewModel, bool animation = true)
            where TViewModel : FunctionViewModel<TResult>
        {
            return viewModel.InvokeAsync<TViewModel, TResult>(null, animation);
        }

        public static Task<TResult> InvokeAsync<TViewModel, TResult>(this ViewModel viewModel, Action<TViewModel, Page> activateAction, bool animation = true)
            where TViewModel : FunctionViewModel<TResult>
        {
            var helper = new InvokeHelper<TViewModel, TResult>();
            viewModel.Navigation.PushAsync<TViewModel>(
                (vm, v) =>
                {
                    vm.Result = default(TResult);
                    if (activateAction != null)
                    {
                        activateAction(vm, v);
                    }
                    helper.Attach(vm, v);
                });
            return helper.TaskCompletionSource.Task;
        }

        public static Task<TResult> InvokeModalAsync<TViewModel, TResult>(this ViewModel viewModel, bool animation = true)
            where TViewModel : FunctionViewModel<TResult>
        {
            return viewModel.InvokeModalAsync<TViewModel, TResult>(null, animation);
        }

        public static Task<TResult> InvokeModalAsync<TViewModel, TResult>(this ViewModel viewModel, Action<TViewModel, Page> activateAction, bool animation = true)
            where TViewModel : FunctionViewModel<TResult>
        {
            var helper = new InvokeHelper<TViewModel, TResult>();
            viewModel.Navigation.PushModalAsync<TViewModel>(
                (vm, v) =>
                {
                    vm.Result = default(TResult);
                    if (activateAction != null)
                    {
                        activateAction(vm, v);
                    }
                    helper.Attach(vm, v);
                });
            return helper.TaskCompletionSource.Task;
        }


    }
}

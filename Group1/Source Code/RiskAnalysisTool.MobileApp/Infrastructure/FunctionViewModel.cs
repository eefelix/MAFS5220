using System;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.Infrastructure
{
    public class FunctionViewModel<T>: ViewModel
    {
        private T _result;

        public T Result
        {
            get { return _result; }
            set { this.SetProperty(ref _result, value); }
        }

        public async Task CloseAsync(T result)
        {
            this.Result = result;
            await this.Navigation.PopAsync();
        } 

    }
}

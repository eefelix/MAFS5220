using System;
using System.Linq;
using System.Collections.Generic;
using RiskAnalysisTool.MobileApp.Infrastructure;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class TaskSettingsViewModel : ViewModel
    {
        public TaskSettingsViewModel()
        {
            _computeBva = true;
            _computeCreditVaR = true;
            _computeFva = true;
            _computeUcva = true;
        }

        private bool _computeBva;
        private bool _computeCreditVaR;
        private bool _computeFva;
        private bool _computeUcva;
        private string _domesticCurrency;
        private string _emailAddress;

        public bool ComputeBva
        {
            get { return _computeBva; }
            set { this.SetProperty(ref _computeBva, value); }
        }

        public string EmailAddress
        {
            get { return _emailAddress; }
            set { SetProperty(ref _emailAddress, value); }
        }

        public bool ComputeCreditVaR
        {
            get { return _computeCreditVaR; }
            set { this.SetProperty(ref _computeCreditVaR, value); }
        }

        public bool ComputeFva
        {
            get { return _computeFva; }
            set { this.SetProperty(ref _computeFva, value); }
        }

        public bool ComputeUCva
        {
            get { return _computeUcva; }
            set { this.SetProperty(ref _computeUcva, value); }
        }

        public string DomesticCurrency
        {
            get { return _domesticCurrency; }
            set { this.SetProperty(ref _domesticCurrency, value); }
        }
    }
}

using System;
using System.Collections.Generic;
using System.Text;
using RiskAnalysisTool.Instruments;
using XLabs.Forms.Mvvm;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class VolatilityDetailViewModel : InstrumentDetailViewModel<ImpliedVolatility>
    {
        private double _volatility;

        public double Volatility
        {
            get { return _volatility; }
            set { this.SetProperty(ref _volatility, value); }
        }

        public override void LoadModel(ImpliedVolatility instrument)
        {
            base.LoadModel(instrument);
            this.Volatility = instrument.Value;
        }

        public override void UpdateModel(ImpliedVolatility instrument)
        {
            base.UpdateModel(instrument);
            instrument.Value = this.Volatility;
        }
    }
}

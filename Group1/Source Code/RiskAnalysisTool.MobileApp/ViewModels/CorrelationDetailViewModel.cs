using System;
using System.Collections.Generic;
using RiskAnalysisTool.Instruments;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class CorrelationDetailViewModel : InstrumentDetailViewModel<Correlation>
    {
        private double _correlation;
        private string _symbol1;
        private string _symbol2;

        public double Correlation
        {
            get { return _correlation; }
            set { this.SetProperty(ref _correlation, value); }
        }

        public string Symbol1
        {
            get { return _symbol1; }
            set { this.SetProperty(ref _symbol1, value); }
        }

        public string Symbol2
        {
            get { return _symbol2; }
            set { this.SetProperty(ref _symbol2, value); }
        }

        public override void LoadModel(Correlation instrument)
        {
            base.LoadModel(instrument);

            this.Symbol1 = instrument.Symbol1;
            this.Symbol2 = instrument.Symbol2;
            this.Correlation = instrument.Value;
        }

        public override void UpdateModel(Correlation instrument)
        {
            base.UpdateModel(instrument);

            instrument.Symbol1 = this.Symbol1;
            instrument.Symbol2 = this.Symbol2;
            instrument.Value = this.Correlation;
        }
    }
}

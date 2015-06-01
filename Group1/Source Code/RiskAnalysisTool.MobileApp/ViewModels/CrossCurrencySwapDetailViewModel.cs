using System;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class CrossCurrencySwapDetailViewModel : InstrumentDetailViewModel<CrossCurrencySwap>
    {
        private int _domesticFrequency;
        private double _domesticRate;
        private int _foreignFrequency;
        private double _foreignNominal;
        private double _foriegnRate;
        private bool _isPayForeign;
        private DateTime _maturityDate;
        private DateTime _startDate;
        private double _swapFxRate;
        private string _foriegnCurrencySymbol;
        private double _fxVolatility;

        public DateTime StartDate
        {
            get { return _startDate; }
            set { SetProperty(ref _startDate, value); }
        }

        public DateTime MaturityDate
        {
            get { return _maturityDate; }
            set { SetProperty(ref _maturityDate, value); }
        }

        public bool IsPayForeign
        {
            get { return _isPayForeign; }
            set { SetProperty(ref _isPayForeign, value); }
        }

        public double DomesticRate
        {
            get { return _domesticRate; }
            set { SetProperty(ref _domesticRate, value); }
        }

        public int DomesticFrequency
        {
            get { return _domesticFrequency; }
            set { SetProperty(ref _domesticFrequency, value); }
        }

        public double ForeignNominal
        {
            get { return _foreignNominal; }
            set { SetProperty(ref _foreignNominal, value); }
        }

        public double SwapFXRate
        {
            get { return _swapFxRate; }
            set { SetProperty(ref _swapFxRate, value); }
        }

        public double ForiegnRate
        {
            get { return _foriegnRate; }
            set { SetProperty(ref _foriegnRate, value); }
        }

        public string ForeignCurrencySymbol
        {
            get { return _foriegnCurrencySymbol; }
            set { SetProperty(ref _foriegnCurrencySymbol, value); }
        }

        public int ForeignFrequency
        {
            get { return _foreignFrequency; }
            set { SetProperty(ref _foreignFrequency, value); }
        }

        public double FxVolatility
        {
            get { return _fxVolatility; }
            set { SetProperty(ref _fxVolatility, value); }
        }

        public override void LoadModel(CrossCurrencySwap instrument)
        {
            base.LoadModel(instrument);

            this.StartDate = instrument.StartDate.Date;
            this.MaturityDate = instrument.MaturityDate.Date;

            this.IsPayForeign = (instrument.PaySymbol != null);

            if (this.IsPayForeign)
            {
                this.DomesticRate = instrument.ReceiveRate;
                this.DomesticFrequency = instrument.FloatingLegFrequency.ToMonths();
                this.ForeignNominal = instrument.PayNominal;
                this.ForiegnRate = instrument.PayRate;
                this.ForeignFrequency = instrument.FixedLegFrequency.ToMonths();
                this.ForeignCurrencySymbol = instrument.PaySymbol.Substring(1);
                this.SwapFXRate = instrument.ReceiveNominal / instrument.PayNominal;
            }
            else
            {
                this.DomesticRate = instrument.PayRate;
                this.DomesticFrequency = instrument.FixedLegFrequency.ToMonths();
                this.ForeignNominal = instrument.ReceiveNominal;
                this.ForiegnRate = instrument.ReceiveRate;
                this.ForeignFrequency = instrument.FloatingLegFrequency.ToMonths();
                this.ForeignCurrencySymbol = instrument.ReceiveSymbol == null ? "" : instrument.ReceiveSymbol.Substring(1);
                this.SwapFXRate = instrument.PayNominal / instrument.ReceiveNominal;
            }
        }

        public override void UpdateModel(CrossCurrencySwap instrument)
        {
            base.UpdateModel(instrument);
            instrument.StartDate = this.StartDate.Date;
            instrument.MaturityDate = this.MaturityDate.Date;

            if (this.IsPayForeign)
            {
                instrument.ReceiveRate = this.DomesticRate;
                instrument.FloatingLegFrequency = new Period(PeriodUnit.Month, this.DomesticFrequency);
                instrument.ReceiveSymbol = "$" + this.ForeignCurrencySymbol;
                instrument.ReceiveNominal = this.ForeignNominal * this.SwapFXRate;
                instrument.PayNominal = this.ForeignNominal;
                instrument.PayRate = this.ForiegnRate;
                instrument.FixedLegFrequency = new Period(PeriodUnit.Month, this.ForeignFrequency);
                instrument.PaySymbol = null;
            }
            else
            {
                instrument.PayRate = this.DomesticRate;
                instrument.FixedLegFrequency = new Period(PeriodUnit.Month, this.DomesticFrequency);
                instrument.PaySymbol = "$" + this.ForeignCurrencySymbol;
                instrument.PayNominal = this.ForeignNominal * this.SwapFXRate;
                instrument.ReceiveNominal = this.ForeignNominal;
                instrument.ReceiveRate = this.ForiegnRate;
                instrument.FloatingLegFrequency = new Period(PeriodUnit.Month, this.ForeignFrequency);
                instrument.ReceiveSymbol = null;
            }
        }
    }
}
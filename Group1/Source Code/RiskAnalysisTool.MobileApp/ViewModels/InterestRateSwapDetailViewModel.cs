using System;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class InterestRateSwapDetailViewModel : InstrumentDetailViewModel<InterestRateSwap>
    {
        private int _fixedFrequency;
        private double _fixedRate;
        private int _floatingFrequency;
        private double _floatingSpread;
        private bool _isPayFloating;
        private DateTime _maturityDate;
        private double _nominal;
        private DateTime _startDate;

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

        public double Nominal
        {
            get { return _nominal; }
            set { SetProperty(ref _nominal, value); }
        }

        public bool IsPayFloating
        {
            get { return _isPayFloating; }
            set { SetProperty(ref _isPayFloating, value); }
        }

        public double FixedRate
        {
            get { return _fixedRate; }
            set { SetProperty(ref _fixedRate, value); }
        }

        public int FixedFrequency
        {
            get { return _fixedFrequency; }
            set { SetProperty(ref _fixedFrequency, value); }
        }

        public double FloatingSpread
        {
            get { return _floatingSpread; }
            set { SetProperty(ref _floatingSpread, value); }
        }

        public int FloatingFrequency
        {
            get { return _floatingFrequency; }
            set { SetProperty(ref _floatingFrequency, value); }
        }

        public override void LoadModel(InterestRateSwap instrument)
        {
            base.LoadModel(instrument);

            this.StartDate = instrument.StartDate.Date;
            this.MaturityDate = instrument.MaturityDate;
            this.FixedRate = instrument.FixedRate;
            if (instrument.IsPayFloating)
            {
                this.FixedFrequency = instrument.FloatingLegFrequency.ToMonths();
                this.FloatingFrequency = instrument.FixedLegFrequency.ToMonths();
            }
            else
            {
                this.FixedFrequency = instrument.FixedLegFrequency.ToMonths();
                this.FloatingFrequency = instrument.FloatingLegFrequency.ToMonths();
            }

            this.FixedRate = instrument.FixedRate;
            this.FloatingSpread = instrument.FloatingSpread;
            this.IsPayFloating = instrument.IsPayFloating;
            this.Nominal = instrument.Nominal;
        }

        public override void UpdateModel(InterestRateSwap instrument)
        {
            base.UpdateModel(instrument);

            instrument.StartDate = this.StartDate.Date;
            instrument.MaturityDate = this.MaturityDate.Date;
            instrument.FixedRate = this.FixedRate;
            if (instrument.IsPayFloating)
            {
                instrument.FloatingLegFrequency = new Period(PeriodUnit.Month, this.FixedFrequency);
                instrument.FixedLegFrequency = new Period(PeriodUnit.Month, this.FloatingFrequency);
            }
            else
            {
                instrument.FixedLegFrequency = new Period(PeriodUnit.Month, this.FixedFrequency);
                instrument.FloatingLegFrequency = new Period(PeriodUnit.Month, this.FloatingFrequency);
            }

            instrument.PaySymbol = null;
            instrument.ReceiveSymbol = null;
            
            instrument.FloatingSpread = this.FloatingSpread;
            instrument.IsPayFloating = this.IsPayFloating;
            instrument.Nominal = this.Nominal;
        }
    }
}
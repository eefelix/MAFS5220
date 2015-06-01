using System;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    public class EquitySwapDetailViewModel : InstrumentDetailViewModel<EquitySwap>
    {
        private double _dividendYield;
        private int _equityFrequency;
        private string _equitySymbol;
        private int _fixedFrequency;
        private double _fixedRate;
        private bool _isPayEquity;
        private DateTime _maturityDate;
        private int _numberOfShares;
        private DateTime _startDate;
        private double _swapPrice;
        private double _volatility;

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

        public bool IsPayEquity
        {
            get { return _isPayEquity; }
            set { SetProperty(ref _isPayEquity, value); }
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

        public string EquitySymbol
        {
            get { return _equitySymbol; }
            set { SetProperty(ref _equitySymbol, value); }
        }

        public double SwapPrice
        {
            get { return _swapPrice; }
            set { SetProperty(ref _swapPrice, value); }
        }

        public int NumberOfShares
        {
            get { return _numberOfShares; }
            set { SetProperty(ref _numberOfShares, value); }
        }

        public double DividendYield
        {
            get { return _dividendYield; }
            set { SetProperty(ref _dividendYield, value); }
        }

        public int EquityFrequency
        {
            get { return _equityFrequency; }
            set { SetProperty(ref _equityFrequency, value); }
        }

        public double Volatility
        {
            get { return _volatility; }
            set { SetProperty(ref _volatility, value); }
        }

        public override void LoadModel(EquitySwap instrument)
        {
            base.LoadModel(instrument);
            
            this.StartDate = instrument.StartDate.Date;
            this.MaturityDate = instrument.MaturityDate.Date;


            if (instrument.PaySymbol == null)
            {
                this.IsPayEquity = false;
                // Receive equity
                this.EquitySymbol = instrument.ReceiveSymbol;
                this.EquityFrequency = instrument.FloatingLegFrequency.ToMonths();
                this.FixedFrequency = instrument.FixedLegFrequency.ToMonths();
            }
            else
            {
                this.IsPayEquity = true;
                this.EquitySymbol = instrument.PaySymbol;
                this.FixedFrequency = instrument.FloatingLegFrequency.ToMonths();
                this.EquityFrequency = instrument.FixedLegFrequency.ToMonths();
            }

            this.Volatility = instrument.Volatility;
            this.NumberOfShares = instrument.Amount;
            this.SwapPrice = instrument.StartPrice;
            this.DividendYield = instrument.DividendYield;
        }

        public override void UpdateModel(EquitySwap instrument)
        {
            base.UpdateModel(instrument);

            instrument.StartDate = this.StartDate.Date;
            instrument.MaturityDate = this.MaturityDate.Date;

            if (this.IsPayEquity)
            {
                instrument.PaySymbol = this.EquitySymbol;
                instrument.ReceiveSymbol = null;
                instrument.FixedLegFrequency = new Period(PeriodUnit.Month, this.EquityFrequency);
                instrument.FloatingLegFrequency = new Period(PeriodUnit.Month, this.FixedFrequency);
            }
            else
            {
                instrument.ReceiveSymbol = this.EquitySymbol;
                instrument.PaySymbol = null;
                instrument.FloatingLegFrequency = new Period(PeriodUnit.Month, this.EquityFrequency);
                instrument.FixedLegFrequency = new Period(PeriodUnit.Month, this.FixedFrequency);
            }

            instrument.Volatility = this.Volatility;

            instrument.Amount = this.NumberOfShares;
            instrument.StartPrice = this.SwapPrice;
            instrument.DividendYield = this.DividendYield;
        }
    }
}
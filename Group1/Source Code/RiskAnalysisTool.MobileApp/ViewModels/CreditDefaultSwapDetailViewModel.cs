using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.MobileApp.ViewModels
{
    internal class CreditDefaultSwapDetailViewModel : InstrumentDetailViewModel<CreditDefaultSwap>
    {
        private int _frequency;
        private double _spread;
        private int _tenor;
        private string _symbol = "ISSUER";

        public int Tenor
        {
            get { return _tenor; }
            set { SetProperty(ref _tenor, value); }
        }

        public int Frequency
        {
            get { return _frequency; }
            set { SetProperty(ref _frequency, value); }
        }

        public double Spread
        {
            get { return _spread; }
            set { SetProperty(ref _spread, value); }
        }

        public string[] AvaialbeSymbols
        {
            get
            {
                return new string[]{
                    "ISSUER", "ISSUER"
                };
            }
        }

        public string Symbol
        {
            get { return _symbol; }
            set { SetProperty(ref _symbol, value); }
        }

        public override void LoadModel(CreditDefaultSwap instrument)
        {
            switch (instrument.Symbol)
            {
                case "INVESTOR":
                    instrument.Symbol = "INVESTOR";
                    break;
                case "ISSUER":
                    instrument.Symbol = "ISSUER";
                    break;
                default:
                    break;
            }

            Frequency = instrument.Tenor.ToMonths();
            Tenor = instrument.Tenor.ToMonths();
            Spread = instrument.Spread;
        }

        public override void UpdateModel(CreditDefaultSwap instrument)
        {
            switch (this.Symbol)
            {
                case "INVESTOR":
                    instrument.Symbol = "INVESTOR";
                    break;
                case "ISSUER":
                    instrument.Symbol = "ISSUER";
                    break;
                default:
                    break;
            }
            instrument.Tenor = new Period(PeriodUnit.Month, this.Frequency);
            instrument.Frequency = new Period(PeriodUnit.Month, this.Frequency);
            instrument.Spread = Spread;
        }
    }
}
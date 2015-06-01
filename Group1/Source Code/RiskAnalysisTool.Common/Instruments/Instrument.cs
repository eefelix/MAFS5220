using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    [KnownType(typeof(Bond))]
    [KnownType(typeof(Swap))]
    [KnownType(typeof(InterestRateSwap))]
    [KnownType(typeof(CrossCurrencySwap))]
    [KnownType(typeof(EquitySwap))]
    [KnownType(typeof(CreditDefaultSwap))]
    public abstract class Instrument
    {
        public Instrument()
        {
            this.StartDate = DateTime.Today;
            this.EvaluationData = DateTime.Today;
            this.MaturityDate = DateTime.Today.AddYears(5);
            this.Symbol = " ";
        }
        [DataMember]
        public DateTime EvaluationData { get; set; }

        [DataMember]
        public double Price { get; set; }

        [DataMember]
        public DateTime MaturityDate { get; set; }

        [DataMember]
        public DateTime StartDate { get; set; }

        [DataMember]
        public string Symbol { get; set; }
    }
}

using System;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class CreditDefaultSwap : Instrument
    {
        public CreditDefaultSwap()
        {
            this.Spread = 0.01;
            this.Tenor = new Period(PeriodUnit.Year, 1);
            this.Frequency = new Period(PeriodUnit.Month, 3);
        }
        [DataMember]
        public double Spread { get; set; }

        [DataMember]
        public Period Tenor { get; set; }

        [DataMember]
        public Period Frequency { get; set; }
    }
}

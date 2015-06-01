using System;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    [KnownType(typeof(InterestRateSwap))]
    [KnownType(typeof(CrossCurrencySwap))]
    [KnownType(typeof(EquitySwap))]
    public class Swap : Instrument
    {

        [DataMember]
        public Period PayLegFrequency { get; set; }

        [DataMember]
        public Period ReceiveLegFrequency { get; set; }

        [DataMember]
        public String PaySymbol { get; set; }

        [DataMember]
        public String ReceiveSymbol { get; set; }


    }
}

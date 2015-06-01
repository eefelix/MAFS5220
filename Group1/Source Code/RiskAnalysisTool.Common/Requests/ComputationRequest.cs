using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Instruments;

namespace RiskAnalysisTool.Requests
{
    [DataContract]
    [KnownType(typeof(CreditVaRRequest))]
    [KnownType(typeof(UCvaRequest))]
    [KnownType(typeof(BvaRequest))]
    [KnownType(typeof(TvaRequest))]
    public class ComputationRequest
    {
        [DataMember]
        public string NotificationEmail { get; set; }
        [DataMember]
        public double[,] Correlation { get; set; }

        [DataMember]
        public ICollection<Instrument> Portfolio { get; set; }

        [DataMember]
        public ICollection<Instrument> MarketData { get; set; }

        [DataMember]
        public double RiskFreeRate { get; set; }

        [DataMember]
        public String DomesticCurrency { get; set; }

        [DataMember]
        public DateTime EvaluationData { get; set; }

    }
}

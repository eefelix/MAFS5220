using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Requests
{
    [DataContract]
    public class TvaRequest : ComputationRequest
    {
        public TvaRequest()
        {
            this.CounterpartyRecoveryRate = 0.4;
            this.InvestorRecoveryRate = 0.4;
        }
        [DataMember]
        public double CounterpartyRecoveryRate { get; set; }

        [DataMember]
        public double InvestorRecoveryRate { get; set; }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Requests
{
    [DataContract]
    public class CreditVaRRequest : ComputationRequest
    {
        public CreditVaRRequest()
        {
            this.ConfidenceLevel = 0.99;
            this.RecoveryRate = 0.4;
        }

        [DataMember]
        public double ConfidenceLevel { get; set; }

        [DataMember]
        public double RecoveryRate { get; set; }
    }
}

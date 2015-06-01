using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Requests
{
    [DataContract]
    public class UCvaRequest : ComputationRequest
    {
        public UCvaRequest()
        {
            this.RecoveryRate = 0.4;
        }

        [DataMember]
        public double RecoveryRate { get; set; }
    }
}

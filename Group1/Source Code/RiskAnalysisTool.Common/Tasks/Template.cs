using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Requests;

namespace RiskAnalysisTool.Tasks
{
    [DataContract]
    public class Template
    {
        [DataMember]
        public ICollection<Instrument> MarketData { get; set; }

        [DataMember]
        public string Name { get; set; }

        [DataMember]
        public IDictionary<string, ComputationRequest> Request { get; set; }
    }
}

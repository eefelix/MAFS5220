using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Requests;
using RiskAnalysisTool.Results;

namespace RiskAnalysisTool.Tasks
{
    [DataContract]
    public class ComputationTask
    {
        [DataMember]
        public int ID { get; set; }

        [DataMember]
        public ComputationRequest Request { get; set; }

        [DataMember(IsRequired = false, EmitDefaultValue = false)]
        public ComputationResult Result { get; set; }

        [DataMember(IsRequired = false, EmitDefaultValue = false)]
        public string Error { get; set; }

    }
}

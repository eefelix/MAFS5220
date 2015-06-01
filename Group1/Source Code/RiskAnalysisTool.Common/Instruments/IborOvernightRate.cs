using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class IborOvernightRate : Instrument
    {
        [DataMember]
        public DateTime Date { get; set; }

        [DataMember]
        public double Rate { get; set; }
    }
}

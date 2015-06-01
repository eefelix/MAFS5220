using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using RiskAnalysisTool.Time;

namespace RiskAnalysisTool.Instruments
{
    [DataContract]
    public class Bond : Instrument
    {
        public Bond()
        {
            this.MaturityDate = DateTime.Today.AddYears(1);
            this.Principle = 100;

        }
        [DataMember]
        public Period CouponFrequency { get; set; }

        [DataMember]
        public double CouponRate { get; set; }

        [DataMember]
        public double Principle { get; set; }

        [DataMember]
        public String Currency { get; set; }

    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RiskAnalysisTool.Instruments;

namespace RiskAnalysisTool.Requests
{
    class BilateralCvaRequest
    {
        public Instrument Instrument { get; set; }

        public ICollection<Instrument> MarketData { get; set; }

        public double RecoveryRate { get; set; }

        public double RiskFreeRate { get; set; }

        public double[,] WwrCorrelation { get; set; }
    }
}

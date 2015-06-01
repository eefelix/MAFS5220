using System;
using System.Linq;

namespace RiskAnalysisTool.Time
{
    public class Period
    {
        public Period(PeriodUnit u, double v)
        {
            Unit = u;
            Value = (int)v;
        }

        public PeriodUnit Unit { get; set; }

        public int Value { get; set; }

        public int ToMonths()
        {
            switch (this.Unit)
            {
                case PeriodUnit.Month:
                    return this.Value;
                case PeriodUnit.Year:
                    return this.Value * 12;
                default:
                    throw new NotSupportedException();

            }
        }

        public override string ToString()
        {
            string unitString;
            switch (this.Unit)
            {
                case  PeriodUnit.Day:
                    unitString = "D";
                    break;
                case  PeriodUnit.Month:
                    unitString = "M";
                    break;
                case  PeriodUnit.Year:
                    unitString = "Y";
                    break;
                default:
                    unitString = "?";
                    break;

            }
            return string.Format("{0}{1}", this.Value, unitString);

        }
    }
}

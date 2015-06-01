using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RiskAnalysisTool.MobileApp.Playground
{
    class Program
    {
        static void Main(string[] args)
        {
            double x = -0.0356;

            string s = string.Format("{0:+#0.00000%;-#0.00000%;+#0.00000%}", x);

            Console.WriteLine(s);
        }
    }
}

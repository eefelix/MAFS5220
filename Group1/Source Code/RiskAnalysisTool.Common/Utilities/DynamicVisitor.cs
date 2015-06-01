using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;

namespace RiskAnalysisTool.Utilities
{
    public abstract class DynamicVisitor<T>
    {
        public void Visit(T target)
        {
            Type visitorType = this.GetType();
            var q = from m in visitorType.GetRuntimeMethods()
                    where m.Name == "OnVisit" && !m.IsStatic && m.ReturnType == typeof(void)
                    let parameters = m.GetParameters()
                    where 
                    parameters.Length == 1
                    let p = parameters[0]
                    where 
                    !p.IsOut && p.ParameterType.GetTypeInfo().IsAssignableFrom(target.GetType().GetTypeInfo())
                    select m;

            var method = q.FirstOrDefault();

            if (method == null)
            {
                throw new ArgumentException("Invalid target type.", "target");
            }

            try
            {
                method.Invoke(this, new object[] {target});
            }
            catch (TargetInvocationException ex)
            {
                throw ex.InnerException;
            }
        }
    }
}

using System;
using System.Collections.Generic;

namespace RiskAnalysisTool.MobileApp.Infrastructure
{
    public interface IModelViewModel<in TModel>
    {
        void LoadModel(TModel instrument);

        void UpdateModel(TModel instrument);
    }
}

using RiskAnalysisTool.Instruments;
using RiskAnalysisTool.Requests;
using RiskAnalysisTool.Results;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Web.Http;
using Microsoft.WindowsAzure.Storage.Blob;
using Microsoft.WindowsAzure.Storage;
using Microsoft.Azure;
using Newtonsoft.Json;


namespace RiskAnalysisTool.WebApp.Controllers
{
    public class ComputeController : ApiController
    {
        [HttpPost]
        public void Post([FromBody] ComputationRequest request)
        {
            CloudStorageAccount account = CloudStorageAccount.Parse(
                CloudConfigurationManager.GetSetting("TaskStorage")
                );

            var blobClient = account.CreateCloudBlobClient();
            var container = blobClient.GetContainerReference("Requests");
            container.CreateIfNotExists();
            Random rand = new Random(Environment.TickCount);
            var blob =
                container.GetBlockBlobReference(string.Format("{0}_{1}", DateTime.Now.ToString("yyyyMMdd"), rand.Next()));

            blob.UploadText(JsonConvert.SerializeObject(request));
        }

        //[HttpGet]
        //public CreditVaRResult Test()
        //{
        //    CreditVaRRequest request = new CreditVaRRequest();
        //    request.AsOfDate = new DateTime(2014, 1, 2);

        //    request.ProductType = ProductType.Equity;
        //    request.RecoveryRate = 0.4;
        //    request.RiskFreeRate = 0.03;
        //    request.UnderlyingExpectedReturn = 0.1;
        //    request.UnderlyingPrice = 100;
        //    request.UnderlyingVolatility = 0.3;
        //    request.WwrCorrelation = 0.0;
        //    request.VaRDate = new DateTime(2014, 6, 30);
        //    request.Maturity = new DateTime(2015, 1, 2);
        //    request.StrikePrice = 110;
        //    request.DefaultModel = new RiskAnalysisTool.Models.AT1PModelSetting() { DefaultBarrier = 0.4, LeverageRatio = 0.6 };

        //    using (ComputeEngine engine = new ComputeEngine())
        //    {
        //        return engine.Compute(request);
        //    }
        //}
    }
}

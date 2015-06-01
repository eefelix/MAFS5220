using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Mail;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using Microsoft.Azure;
using Microsoft.WindowsAzure.ServiceRuntime;
using Microsoft.WindowsAzure.Storage;
using Microsoft.WindowsAzure.Storage.Blob;
using Microsoft.WindowsAzure.Storage.Table;
using Newtonsoft.Json;
using RiskAnalysisTool.Calculation;
using RiskAnalysisTool.Requests;
using RiskAnalysisTool.Results;
using RiskAnalysisTool.Tasks;
using RiskAnalysisTool.Worker.DAO;
using SendGrid;

namespace RiskAnalysisTool.Worker
{
    public class WorkerRole : RoleEntryPoint
    {
        private bool _stopping = false;
        public override bool OnStart()
        {
            ServicePointManager.DefaultConnectionLimit = 12;
            _stopping = false;
            return base.OnStart();
        }

        public override void OnStop()
        {
            _stopping = true;
            base.OnStop();
        }

        public override void Run()
        {
            while (!_stopping)
            {
                CloudStorageAccount account = CloudStorageAccount.Parse(
                    CloudConfigurationManager.GetSetting("TaskStorage")
                    );

                var blobClient = account.CreateCloudBlobClient();
                var container = blobClient.GetContainerReference("Requests");
                container.CreateIfNotExists();

                var item = container.ListBlobs().FirstOrDefault();
                if (item == null)
                {
                    Thread.Sleep(5);
                }

                var blob = new CloudBlockBlob(item.Uri);
                var request = JsonConvert.DeserializeObject<ComputationRequest>(blob.DownloadText());

                string title = null;
                string text = null;

                title = "Task computing started";
                text = JsonConvert.SerializeObject(request, Formatting.Indented);

                SendEmail(request, title, text);

                ComputationResult result = null;
                try
                {

                    using (ComputeEngine engine = new ComputeEngine())
                    {
                        result = engine.calculate(request);
                    }
                    title = "Task computing complete";
                    text = result.ToString();
                }
                catch (Exception ex)
                {
                    title = "Task computing failed";
                    text = ex.Message;
                }

                SendEmail(request, title, text);

            }
        }

        private static void SendEmail(ComputationRequest request, string title, string text)
        {
            // Create network credentials to access your SendGrid account
            var username = "azure_57a2560e78ffeca3a24b33fb140a9eaa@azure.com";
            var pswd = "6Ujk6BWvIUfT213";

            var credentials = new NetworkCredential(username, pswd);

            // Create the email object first, then add the properties.
            SendGridMessage myMessage = new SendGridMessage();
            myMessage.AddTo(request.NotificationEmail);
            myMessage.From = new MailAddress("azure_57a2560e78ffeca3a24b33fb140a9eaa@azure.com", "Risk Analysis Tool");
            myMessage.Subject = title;
            myMessage.Text = text;

            // Create an Web transport for sending email.
            var transportWeb = new Web(credentials);

            // Send the email.
            // You can also use the **DeliverAsync** method, which returns an awaitable task.
            transportWeb.DeliverAsync(myMessage).Wait();
        }
    }
}

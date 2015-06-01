using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json.Linq;

namespace RiskAnalysisTool.MobileApp.Infrastructure
{
    internal class ServiceHelper
    {
        private static readonly Uri BaseUri = new Uri("http://riskanalysistool.cloudapp.net/");

        public static async Task<TResult> CallService<TData, TResult>(string relativeUri, string method = "GET", TData data = null)
            where TData : class
        {
            string requestData;
            switch (method.ToUpperInvariant())
            {
                case "GET":
                case "DELETE":
                    if (object.ReferenceEquals(data, null))
                    {
                        requestData = "";
                    }
                    else
                    {
                        if (typeof(TData) == typeof(string))
                        {
                            requestData = data.ToString();
                            ;
                        }
                        else
                        {
                            var properties = from p in typeof(TData).GetRuntimeProperties()
                                             where p.CanRead
                                             let getter = p.GetMethod
                                             where getter != null && !getter.IsStatic && getter.IsPublic
                                             select p;
                            requestData = string.Join("&", properties.Select(p => string.Format("{0}={1}", p.Name, p.GetValue(data))));
                        }
                    }
                    break;

                case "POST":
                case "PUT":
                    if (object.ReferenceEquals(data, null))
                    {
                        requestData = "";
                    }
                    else
                    {
                        if (typeof(TData) == typeof(string))
                        {
                            requestData = data.ToString();
                            ;
                        }
                        else if (typeof(TData).IsArray)
                        {
                            JArray jobject = new JArray(data);
                            requestData = jobject.ToString();
                        }
                        else
                        {
                            JObject jobject = new JObject(data);
                            requestData = jobject.ToString();
                        }
                    }
                    break;

                default:
                    throw new InvalidOperationException(string.Format("Unsupported HTTP method: {0}", method));
            }

            string response = await CallService(relativeUri, method, requestData, "application/json");

            JObject obj = JObject.Parse(response);
            return obj.Value<TResult>();
        }

        public static async Task<string> CallService(string relativeUri, string method = "GET", string data = "", string contentType = "application/json")
        {
            Uri uri = new Uri(BaseUri, relativeUri);
            using (HttpClient client = new HttpClient())
            {
                HttpResponseMessage response;
                switch (method.ToUpperInvariant())
                {
                    case "GET":
                        response = await client.GetAsync(uri + "?" + data);
                        break;
                    case "POST":
                        response = await client.PostAsync(uri, new ByteArrayContent(Encoding.UTF8.GetBytes(data)));
                        break;
                    case "PUT":
                        response = await client.PutAsync(uri, new ByteArrayContent(Encoding.UTF8.GetBytes(data)));
                        break;
                    case "DELETE":
                        response = await client.DeleteAsync(uri);
                        break;
                    default:
                        throw new ArgumentException("Unsupported HTTP method.", "method");
                }

                response.EnsureSuccessStatusCode();
                return await response.Content.ReadAsStringAsync();
            }
        }
    }
}

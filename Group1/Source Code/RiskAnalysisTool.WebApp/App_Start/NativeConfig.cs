using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Web;

namespace RiskAnalysisTool.WebApp
{
    public static class NativeConfig
    {
        //[DllImport("kernel32", SetLastError = true, CharSet = CharSet.Ansi)]
        //static extern IntPtr LoadLibrary([MarshalAs(UnmanagedType.LPStr)]string lpFileName);

        public static void ConfigureNativePath()
        {
            string path = Environment.GetEnvironmentVariable("PATH");
            if (!path.EndsWith(";"))
            {
                path += ";";
            }
            path += HttpRuntime.BinDirectory;
            Environment.SetEnvironmentVariable("PATH", path, EnvironmentVariableTarget.Process);

            //IntPtr hDll = LoadLibrary(Path.Combine(HttpRuntime.BinDirectory,"RiskAnalysisTool.CalculationLibrary.dll"));
            //int err = Marshal.GetLastWin32Error();
        }
    }
}